/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Yury Lysogorskiy (ICAMS)
   GRACE-3L KOKKOS implementation
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/3l/kk,             PairGRACE3LKokkos_Mixed_Device);
PairStyle(grace/3l/kk/device,      PairGRACE3LKokkos_Mixed_Device);
PairStyle(grace/3l/kk/fp32,        PairGRACE3LKokkos_FP32_Device);
PairStyle(grace/3l/kk/fp32/device, PairGRACE3LKokkos_FP32_Device);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRACE_3L_KOKKOS_H
#define LMP_PAIR_GRACE_3L_KOKKOS_H

#include "pair.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "pair_kokkos.h"

#include <string>
#include <vector>
#include <cstdint>

// Opt 15: cuBLAS true-fp32 SGEMM for the radial MLP (CUDA backend only).
#ifdef KOKKOS_ENABLE_CUDA
#include <cublas_v2.h>
#endif

namespace LAMMPS_NS {

// ============================================================================
// GRACE3LModel: Host-side weight container loaded from .npz file.
//
// Holds EVERY weight tensor grace_3l_weights.npz carries (exporter
// kokkos_export_3l.py), stored RAW (float64 on host) exactly as written.
// copy_weights_to_device() casts these into the NNScalar/GeomScalar device
// views. The forward/backward kernels (Tasks 2.3-3.4) consume the device
// layout; this struct is the CPU staging area + the name/index directory.
//
// NORM CONVENTION (PINNED, applied exactly ONCE downstream): every Linear-like
// weight is stored RAW together with a SEPARATE scalar norm = 1/sqrt(n_in).
// The consumer forms w*norm (and b*norm) once. The ONLY place this file folds
// the norm early is the per-element species-modulation precompute
// (SPBF::z_tr / z_proj built in copy_weights_to_device): there chem_embed is
// contracted with the (lin_transform / chem_linear) weight and the norm folded
// in once — mirroring GRACE-2L's d_z_tr_A. Everywhere else weights go to the
// device raw and the kernel multiplies the stored norm scalar once.
// ============================================================================
struct GRACE3LModel {
  // ---- Architecture metadata ----
  int n_elements = 0;
  int embedding_size = 0;
  int lmax = 0;                 // angular lmax (== SPBF _lmax; no bare npz key -> from A1_lmax)
  int nradbase = 0;             // Chebyshev basis count (== SPBF n_rad_basis)
  int radial_basis_p = 0;       // envelope polynomial order used by geometry (== SPBF p)
  int radial_basis_p_global = 0; // the (unused-by-geometry) global radial_basis_p key
  double rcut = 0.0;
  bool has_bond_specific_cutoff = false;
  std::vector<double> bond_cutoff_map;   // [n_elements * n_elements]

  // Per-layer radial ceilings (3 layers). Set from spbf[k].n_rad_max.
  int L1_nradmax = 0;
  int L2_nradmax = 0;
  int L3_nradmax = 0;

  // Chemical embedding [n_elements, embedding_size]
  std::vector<double> chem_embedding;

  // Per-element energy shifts and global output scale
  std::vector<double> shift_values;      // [n_elements]
  double output_scale = 1.0;             // ConstantScaleShiftTarget; 1.0 if absent

  // Element names (model index -> chemical symbol)
  std::vector<std::string> element_names;

  // ---- UQ / extrapolation-grade artifacts (optional; present iff has_uq) ----
  // Schema v6 (basis-RP, normalize+density): the GMM feature phi is built from the
  // invariant (l=0) energy-path B-basis of the three per-layer scalar reduces
  // (rho1, rho2, rho3) as proj = (B/||B||) * R (L2-normalized) with n_density
  // log-norm density channels appended: phi = concat(proj, dens). R is stored
  // verbatim; gamma is the only UQ signal (forward only, no force-error model).
  bool has_uq = false;
  int uq_schema_version = 0;
  int uq_n_elements = 0;
  int uq_max_clusters = 0;            // Kmax
  int uq_feature_dim = 0;             // D = rp_dim + n_density (full GMM feature dim)
  int uq_rp_dim = 0;                  // random-projection output dim (R cols; e.g. 128)
  int uq_n_density = 0;               // density channels appended after proj (D - rp_dim)
  int uq_normalize = 0;              // 1 -> L2-normalize basis before projection (v6)
  double uq_density_scale = 1.0;     // density_scale (=1.0 for v6)
  int uq_d_basis = 0;                 // concatenated invariant basis width (rows of R)
  std::vector<double> uq_centroids;          // [E*Kmax*D]
  std::vector<double> uq_inv_cov;            // [E*Kmax*D*D] (pre-inverted precision)
  std::vector<int>    uq_n_clusters;         // [E]
  std::vector<double> uq_interp_thresholds;  // [E*Kmax]
  std::vector<double> uq_rp_matrix;          // [D_basis*rp_dim] row-major projection R

  // ---- Generic MLP (radial MLP per SPBF; energy readout MLP) ----
  // Each layer stored raw: W flat [n_in*n_out] row-major, optional bias b
  // flat [n_out] (hidden layers only), scalar norm. silu on hidden layers,
  // linear (no bias/activation) on the output layer.
  struct MLP {
    int n_layers = 0;
    bool has_bias = false;
    std::vector<std::vector<double>> W;   // [n_layers][n_in*n_out]
    std::vector<std::vector<double>> b;   // [n_layers][n_out] (empty entry = no bias)
    std::vector<double> norm;             // [n_layers]
    std::vector<int> n_in;                // [n_layers]
    std::vector<int> n_out;               // [n_layers]
  };

  // ---- SPBF (Single Particle Basis Function) per layer: A1 (scalar), A2/A3 (equivariant) ----
  struct SPBF {
    std::string name;
    bool equivariant = false;
    int n_rad_max = 0, n_rad_basis = 0, lmax = 0, Lmax = 0, p = 0;
    double rcut = 0.0, inv_avg_n_neigh = 1.0;
    std::vector<int> l_tile;              // [nlm] radial-channel -> lm tile
    MLP mlp;                              // radial MLP (bias on hidden layers)
    // scalar mode (A1): species modulation z_tr = lin_transform(Z)
    std::vector<double> lin_transform_W;  // [embedding_size, n_rad_max]
    double lin_transform_norm = 1.0;
    // equivariant mode (A2/A3): L=0 chemistry injection + dense CG couple
    std::vector<double> chem_linear_W;    // [embedding_size, n_rad_max]
    double chem_linear_norm = 1.0;
    std::vector<double> chem_l0_mask;     // [n_lm_ind]
    std::vector<double> cg_W;             // [(lmax+1)^2 * n_lm_ind, nfunc]
    int nfunc = 0, n_lm_ind = 0;
  };

  // ---- GeneralProductFunction, cp_l mode (use_S=False => rank == n_out) ----
  struct CPL {
    std::string name;
    int rank = 0, nfunc = 0, n_left = 0, n_right = 0;
    int n_groups_left = 0, n_groups_right = 0, n_cg = 0;
    double norm_u = 1.0, norm_v = 1.0;
    std::vector<double> U;                // [n_groups_left, rank, n_left]
    std::vector<double> V;                // [n_groups_right, rank, n_right]
    std::vector<double> cg;               // [n_cg]
    std::vector<int> group_left, group_right;   // [n_lm_left], [n_lm_right]
    std::vector<int> left_ind, right_ind, m_sum_ind;  // [n_cg]
    std::vector<int> out_l, out_parity;   // [nfunc]
  };

  // ---- FunctionReduceN (metadata driven; per-instruction collectors) ----
  struct ReduceInstr {
    std::string name;
    std::vector<double> W;                // flat: elem_dep [n_types,n_out,n_in,w_shape] else [1,n_out,n_in,w_shape]
    double norm = 1.0;
    int n_in = 0, w_shape = 0, n_conn = 0;
    std::vector<int> collect_ind;         // [n_conn]
    std::vector<int> w_l_tile;            // [n_conn]
    std::vector<int> total_sum_ind;       // [n_conn]
  };
  struct Reduce {
    std::string name;
    int n_out = 0, n_funcs = 0;
    bool only_invar = false, elem_dep = false;
    std::vector<ReduceInstr> instr;       // 1 (product reduces) or 4 (eq*/rho*)
    std::vector<double> norm_map;         // [n_funcs] (or [1] for only_invar)
    bool has_norm_map = false;
  };

  // ---- FCRight2Left (metadata driven; 3L FCs are right->left, left_coefs=False) ----
  struct FC {
    std::string name;
    int n_out = 0;
    bool left_coefs = false;
    int n_funcs_left = 0, n_funcs_right = 0;
    int n_in_left = 0, n_in_right = 0, w_shape_left = 0, w_shape_right = 0;
    double norm_left = 1.0, norm_right = 1.0;
    std::vector<double> w_left;           // [n_out, n_in_left, w_shape_left] (empty if !left_coefs)
    std::vector<double> w_right;          // [n_out, n_in_right, w_shape_right]
    std::vector<double> norm_out_factor;  // [n_funcs_left]
    std::vector<int> w_tile_left, w_tile_right;    // [n_funcs_left], [n_funcs_right]
    std::vector<int> collect_to, collect_from;     // [n_funcs_left], [n_funcs_right]
  };

  // ---- EquivariantRMSNorm (eq1_norm, eq2_norm) ----
  struct EqNorm {
    std::string name;
    int center_l0 = 0, n_groups = 0, n_out = 0, M = 0;
    double eps = 0.0;
    std::vector<double> affine_weight;    // [n_groups, n_out]
    std::vector<double> degree_weights;   // [M]
    std::vector<double> l0_mask;          // [M]
    std::vector<int> expand_index;        // [M]
  };

  // ---- InvariantLayerRMSNorm (rho1/2/3_norm) ----
  struct RMSNorm {
    std::string name;
    int type = 0, n_out = 0;              // type: 0=full (scale len n_out), 1=only_nonlin (len n_out-1)
    std::vector<double> scale;            // [n_out] or [n_out-1]
  };

  // ---- Ordered inventories (order == the *_names index lists in the npz) ----
  std::vector<SPBF> spbf;                 // [A1, A2, A3]
  std::vector<CPL> prod;                  // [A1_2, A1_3, A1_4, A2_2, A2_3, A2_4, A3_2, A3_3, A3_4]
  std::vector<Reduce> reduce;             // 16 reduces (product reduces + eq*/rho*)
  std::vector<FC> fc;                     // [A1_2a, A1_2b, A2_2a, A2_2b, A3_2a, A3_2b]
  std::vector<EqNorm> eqnorm;             // [eq1_norm, eq2_norm]
  std::vector<RMSNorm> rmsnorm;           // [rho1_norm, rho2_norm, rho3_norm]

  // ---- Readout ----
  MLP energy_mlp;                         // (31->64) silu, (64->1) linear, no bias
  int energy_mlp_activation = 0;          // 0=silu, 1=tanh
  std::vector<std::string> out1_origins;  // [rho1_norm, rho2_norm, rho3_norm]
  int energy_mlp_origin_n_out = 0;        // per-origin density width (32)
  int energy_mlp_n_origins = 0;           // 3

  void load(const std::string &filepath);
};

// ============================================================================
// PairGRACE3LKokkos: KOKKOS pair style for GRACE-3L
//
// Compile-time capacity constants. These MUST agree with the exporter's
// authoritative CAPS_3L table (kokkos_export_3l.py). Named GRACE3L_<CAPS-key>
// so the mapping to that table is one-to-one.
// ============================================================================
static constexpr int GRACE3L_MAX_NRADMAX       = 64;     // SPBF.n_rad_max (A1/A2/A3)
static constexpr int GRACE3L_LMAX              = 4;      // SPBF._lmax / Lmax
static constexpr int GRACE3L_MAX_NRADBASIS     = 12;     // SPBF.n_rad_basis (Chebyshev)
static constexpr int GRACE3L_MAX_MLP_DIM       = 320;    // radial-MLP hidden / output width
static constexpr int GRACE3L_MAX_MLP_LAYERS    = 4;      // radial MLP depth
static constexpr int GRACE3L_MAX_SPBF_NFUNC    = 960;    // A3 CG-couple nfunc
static constexpr int GRACE3L_MAX_SPBF_NLMIND   = 64;     // A3 n_lm_indicator
static constexpr int GRACE3L_MAX_PROD_RANK     = 80;     // A*_2 cp_l rank
static constexpr int GRACE3L_MAX_PROD_NFUNC    = 1920;   // A2_3 product nfunc
static constexpr int GRACE3L_MAX_PROD_NIN      = 160;    // max(n_left,n_right) across products
static constexpr int GRACE3L_MAX_PROD_NGROUPS  = 16;     // n_groups_left/right (L2/L3 products)
static constexpr int GRACE3L_MAX_PROD_NCG      = 10400;  // A2_3 cp_l coupling terms
static constexpr int GRACE3L_MAX_REDUCE_NOUT   = 160;    // A*_2_red / A*_red n_out
static constexpr int GRACE3L_MAX_REDUCE_NIN    = 160;    // largest reduce collector input width
static constexpr int GRACE3L_MAX_REDUCE_NFUNCS = 64;     // largest reduce coupling_meta_data
static constexpr int GRACE3L_MAX_REDUCE_WSHAPE = 336;    // A2_3_red reducing-weight w_shape
static constexpr int GRACE3L_MAX_FC_NOUT       = 160;    // FC n_out / feature width
static constexpr int GRACE3L_MAX_FC_TILE       = 16;     // FC w_right tile dim (A2/A3 FCs)
static constexpr int GRACE3L_MAX_FC_NFUNCS     = 64;     // FC n_funcs_left/right accumulator cap (Task 2.6)
static constexpr int GRACE3L_MAX_RHO_NOUT      = 40;     // rho*_norm / rho* density width
static constexpr int GRACE3L_MAX_ENERGY_MLP_DIM = 80;    // readout MLP hidden
static constexpr int GRACE3L_MAX_ATOM_TYPES    = 112;    // n_elements / elem-dependent reduce dim

// Inventory-size ceilings for the fixed device-side arrays below. These bound
// the *count* of instructions (not their inner dims — those use GRACE3L_MAX_*
// above). GRACE-3L-OMAT has 3 SPBF / 9 products / 16 reduces / 6 FCs /
// 2 eqnorms / 3 rmsnorms and <=4 collectors per reduce; the values carry slack.
static constexpr int GRACE3L_N_SPBF          = 3;   // A1, A2, A3 (architectural)
static constexpr int GRACE3L_MAX_PRODS       = 12;  // cp_l products (9 used)
static constexpr int GRACE3L_MAX_REDUCES     = 20;  // FunctionReduceN (16 used)
static constexpr int GRACE3L_MAX_FCS         = 8;   // FCRight2Left (6 used)
static constexpr int GRACE3L_MAX_EQNORMS     = 4;   // EquivariantRMSNorm (2 used)
static constexpr int GRACE3L_MAX_RMSNORMS    = 4;   // InvariantLayerRMSNorm (3 used)
static constexpr int GRACE3L_MAX_REDUCE_INSTR = 6;  // collectors per reduce (4 used)
static constexpr int GRACE3L_MAX_REDUCE_ROUNDS = 64; // Opt 17: max cuBLAS fwd rounds (=max connections sharing one output f; 34 observed)
static constexpr int GRACE3L_UQ_MAX_RP_DIM = 136;    // bounds the FULL UQ feature dim D = rp_dim + n_density (z/f/delta stack arrays); 136 = 128 rp_dim + up to 4 density channels + slack; a larger model fails the D > MAX guard at load

template<class DeviceType, typename NNScalarT = float, typename GeomScalarT = double>
class PairGRACE3LKokkos : public Pair, public KokkosBase {
 public:
  // Kernel tags
  struct TagPackForwardComm{};
  struct TagUnpackForwardComm{};
  struct TagPackReverseComm{};
  struct TagUnpackReverseComm{};
  struct TagComputeNeigh{};
  struct TagComputeRadialBasis{};
  struct TagComputeMLPRadial{};      // generic radial MLP (per bond, spbf-indexed) -> d_R1_nl
  struct TagComputeA1{};             // Layer-1 scalar-indicator SPBF basis -> d_A1
  struct TagComputeA2{};             // Layer-2 equivariant-indicator SPBF basis -> d_A2
  struct TagComputeA3{};             // Layer-3 equivariant-indicator SPBF basis -> d_A3
  struct TagComputeSPBFFused{};      // Opt 6: fused equiv-SPBF (A2/A3) — one team per (ii,n); prod slice + CG couple in team shared
  struct TagCPLProject{};            // cp_l Pass A: per-lm U/V projection to rank axis
  struct TagCPLCouple{};             // cp_l Pass B: bilinear CG on rank axis + segment-sum
  struct TagCPLFused{};              // Opt 5: fused fwd cp_l (team/(a,r): project+couple in shared, out via shared atomics)
  struct TagCPLBwdFused{};           // Opt 5: fused bwd cp_l (team/(a,r): recompute-project + bwd-couple in shared)
  struct TagReduceN{};                // Task 2.6: generic FunctionReduceN (multi-instr, elem-dep)
  struct TagFCRight2Left{};           // Task 2.6: generic FCRight2Left (right->left scatter)
  struct TagEqNorm{};                 // Task 2.7: generic EquivariantRMSNorm (eq1_norm, eq2_norm)
  struct TagRMSNorm{};                // Task 2.7: generic InvariantLayerRMSNorm (rho*_norm)
  struct TagComputeMLPEnergy{};       // Task 2.11: 3-density readout (rho1/2/3_norm) -> per-atom energy
  struct TagComputeUQ_rho1{};         // basis-RP: rho1 (L1) block -> raw proj z (Phase 1; d_A1.. current) — forward only
  struct TagComputeUQ_rho2{};         // basis-RP: rho2 (L2) block += into z (Phase 2; d_A2_red.. current) — forward only
  struct TagComputeUQ{};              // basis-RP: rho3 (L3) block += z, normalize + density + GMM (gamma/sigma/gmm_cluster) in Phase 3 — forward only

  // ---- Task 3.1: backward (reverse-mode adjoint) tags ----
  struct TagReadoutBwd{};             // readout backward -> d_rho{1,2,3}_norm_adj (all equal)
  struct TagRMSNormBwd{};             // generic InvariantLayerRMSNorm backward (adj_out -> adj_in)
  struct TagReduceNBwd{};             // generic FunctionReduceN backward (per-instruction, accumulate)
  struct TagFCBwdLeft{};              // generic FCRight2Left backward, left path -> adj_left
  struct TagFCBwdRight{};             // generic FCRight2Left backward, right path -> adj_right (atomic)
  struct TagCPLBwdCouple{};           // generic cp_l backward Pass B (couple) -> adj_lproj/adj_rproj
  struct TagCPLBwdProject{};          // generic cp_l backward Pass A (project) -> adj_left/adj_right
  struct TagComputeAdjProd{};         // Opt 4: adj_prod = inv_avg*A_adj@cg_W^T -> d_adj_prod (output-parallel, shared by SPBF-bwd + ForceEquiv)
  struct TagA3SPBFBwd{};              // L3 equiv-SPBF backward: d_A3_adj -> d_eq2_norm_adj (message adj)
  struct TagEqNormBwd{};              // Task 3.2: generic EquivariantRMSNorm VJP (adj_out -> adj_in)
  struct TagA2SPBFBwd{};              // Task 3.2: L2 equiv-SPBF backward: d_A2_adj -> d_grad_I_global (eq1_norm message adj)
  struct TagSPBFBwdFused{};           // Opt 10: team+shared equiv-SPBF backward (A2/A3) — adj_prod row staged in shared once
  // ---- Task 3.3: radial r-derivative + per-bond force (geometry backward) tags ----
  struct TagComputeMLPRadialDeriv{}; // radial MLP forward + dR/dr -> d_R1_nl, d_DR1_nl, d_denv (per layer)
#ifdef KOKKOS_ENABLE_CUDA
  struct TagMLPAssembleVal{};        // Opt 15: build cuBLAS input X[M x n_in0] (value path), zero padding rows
  struct TagMLPAssembleDeriv{};      // Opt 15b: build X[M x n_in0] + dX[M x nradbase] + d_denv (deriv path)
  struct TagMLPActVal{};             // Opt 15: silu activation between GEMMs (value only)
  struct TagMLPActDeriv{};           // Opt 15b: silu (value) + silu'(s)*ds (deriv) between GEMMs
  struct TagMLPScatterR{};           // Opt 15: scatter LayoutRight R/DR[M x d3] -> LayoutLeft d_R1_nl/d_DR1_nl
#endif
  struct TagForceL1{};               // L1 scalar-SPBF geometry backward -> d_f_ij
  struct TagForceEquiv{};            // L2/L3 equiv-SPBF geometry backward -> d_f_ij
  struct TagForceEquivFused{};       // Opt 11: team-per-(ii,n) ForceEquiv — adj_prod row staged in shared, one bond per thread

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  using GeomScalar = GeomScalarT;
  using NNScalar = NNScalarT;

  PairGRACE3LKokkos(class LAMMPS *);
  ~PairGRACE3LKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  // UQ / extrapolation grade (forward only): expose GMM gamma/sigma triggers +
  // per-atom arrays to `fix pair` (same interface as GRACE-1L/2L).
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  // MPI communication (host-side fallback)
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  // Kokkos-native communication (device-side)
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d,
                                DAT::tdual_double_1d &, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &) override;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPackForwardComm, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagUnpackForwardComm, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPackReverseComm, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagUnpackReverseComm, const int&) const;

  // Kernel operators (Task 2.3: neighbor + Chebyshev radial + Y_lm precompute)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeNeigh, const typename Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeRadialBasis, const typename Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>::member_type& team) const;

  // Kernel operators (Task 2.4: L1 scalar-indicator SPBF basis A1)
  // Generic radial MLP: serves A1 (Task 2.4) and A2 (Task 2.8), spbf-indexed by
  // the member radial_mlp_spbf set before launch.
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPRadial, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeA1, const typename Kokkos::TeamPolicy<DeviceType, TagComputeA1>::member_type& team) const;

  // Kernel operator (Task 2.8: L2 equivariant-indicator SPBF basis A2)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeA2, const typename Kokkos::TeamPolicy<DeviceType, TagComputeA2>::member_type& team) const;

  // Kernel operator (Task 2.10: L3 equivariant-indicator SPBF basis A3)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeA3, const typename Kokkos::TeamPolicy<DeviceType, TagComputeA3>::member_type& team) const;

  // Kernel operator (Opt 6: fused equiv-SPBF for A2/A3 — team per (ii,n), prod slice + CG couple in shared)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeSPBFFused, const typename Kokkos::TeamPolicy<DeviceType, TagComputeSPBFFused>::member_type& team) const;

  // Kernel operators (Task 2.5: generic cp_l GeneralProductFunction product)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLProject, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLCouple, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLFused,
      const typename Kokkos::TeamPolicy<DeviceType, TagCPLFused>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLBwdFused,
      const typename Kokkos::TeamPolicy<DeviceType, TagCPLBwdFused>::member_type& team) const;

  // Kernel operators (Task 2.6: generic FunctionReduceN + FCRight2Left)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagReduceN, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagFCRight2Left, const int&) const;

  // Kernel operators (Task 2.7: generic EquivariantRMSNorm + InvariantLayerRMSNorm)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagEqNorm, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagRMSNorm, const int&) const;

  // Kernel operator (Task 2.11: 3-density readout -> per-atom energy)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPEnergy, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeUQ_rho1, const typename Kokkos::TeamPolicy<DeviceType, TagComputeUQ_rho1>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeUQ_rho2, const typename Kokkos::TeamPolicy<DeviceType, TagComputeUQ_rho2>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeUQ, const typename Kokkos::TeamPolicy<DeviceType, TagComputeUQ>::member_type& team) const;

  // Kernel operators (Task 3.1: readout + L3 backward / generic bwd routines)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagReadoutBwd, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagRMSNormBwd, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagReduceNBwd, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagFCBwdLeft, const int, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagFCBwdRight, const int, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLBwdCouple, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagCPLBwdProject, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeAdjProd, const int, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagA3SPBFBwd, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagEqNormBwd, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagA2SPBFBwd, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSPBFBwdFused, const typename Kokkos::TeamPolicy<DeviceType, TagSPBFBwdFused>::member_type& team) const;

  // ---- Task 3.3 force kernels ----
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPRadialDeriv, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadialDeriv>::member_type& team) const;
#ifdef KOKKOS_ENABLE_CUDA
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPAssembleVal, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPAssembleDeriv, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPActVal, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPActDeriv, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPScatterR, const int) const;
#endif
  KOKKOS_INLINE_FUNCTION
  void operator() (TagForceL1, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagForceEquiv, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagForceEquivFused, const typename Kokkos::TeamPolicy<DeviceType, TagForceEquivFused>::member_type& team) const;
  // Shared per-bond geometry-force chain (spec §4.1): given the per-lm weight
  // array w[nlm_y] (= dE/da_nl for this (ii,jj,n) up to the R*Y*env factor),
  // accumulates the Cartesian force f = Σ_{l,m} w*(Y*DRenv*rhat + DY*Renv/r)
  // into (f0,f1,f2). Renv = R*env, DRenv = dR/dr*env + R*denv/dr; DY is the
  // tangential (rhat-orthogonal) gradient of the real SH. Reused by L1/L2/L3.
  KOKKOS_INLINE_FUNCTION
  void bond_geom_force(const int ii, const int jj, const int n,
      const NNScalar* w, GeomScalar& f0, GeomScalar& f1, GeomScalar& f2) const;

 protected:
  int inum, maxneigh, chunk_size, chunk_offset;
  int host_flag;
  int eflag, vflag;
  int neighflag;
  int chunksize;

  // Architecture params
  int nelements, lmax, nradbase, embedding_size;
  int radial_basis_p;
  GeomScalar rcut;

  // Layer-specific radial dimensions
  int L1_nradmax, L2_nradmax, L3_nradmax;

  // I feature dims (drive comm buffer sizing). Zero-initialized; the task that
  // wires up real message comm must set these AND comm_forward/comm_reverse
  // (mirroring pair_grace_2l_kokkos) before any forward/reverse_comm call.
  int I_n_out = 0, I_n_funcs = 0;
  // Round-2 (L3) indicator dims (Task 2.10). The Layer-3 indicator is eq2_norm,
  // which is PARITY-DOUBLED vs eq1_norm: I2_n_out = eq2_norm.M (=50) is twice
  // I_n_out (=25), so round 2 needs its own buffer d_eq2_norm_global (below) —
  // d_I_global [nall,64,25] cannot hold a [nall,64,50] payload.
  int I2_n_out = 0, I2_n_funcs = 0;
  // Forward-comm stage selector (Task 2.10): 1 = round 1 (eq1_norm -> d_I_global,
  // I_n_funcs*I_n_out/atom); 2 = round 2 (eq2_norm -> d_eq2_norm_global,
  // I2_n_funcs*I2_n_out/atom). Set on host immediately before each
  // comm->forward_comm; captured by value into the pack/unpack tag-functors.
  int comm_stage = 1;

  // LAMMPS arrays
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  DAT::ttransform_kkacc_1d_9 k_cvatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  typename AT::t_kkacc_1d_9 d_cvatom;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;

  typedef Kokkos::DualView<KK_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq, k_scale;
  typedef Kokkos::View<KK_FLOAT**, DeviceType> t_fparams;
  t_fparams d_cutsq, d_scale;
  typename AT::t_int_1d d_map;

  int need_dup;
  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;
  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  DupScatterView<KK_ACC_FLOAT*[9], typename DAT::t_kkacc_1d_9::array_layout> dup_cvatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*[9], typename DAT::t_kkacc_1d_9::array_layout> ndup_cvatom;

  friend void pair_virial_fdotr_compute<PairGRACE3LKokkos>(PairGRACE3LKokkos*);

  // ---------- Kokkos View typedefs ----------
  typedef Kokkos::View<int*, DeviceType> t_int_1d;
  typedef Kokkos::View<int**, DeviceType> t_int_2d;
  typedef Kokkos::View<GeomScalar*, DeviceType> t_geom_1d;
  typedef Kokkos::View<GeomScalar**, DeviceType> t_geom_2d;
  typedef Kokkos::View<GeomScalar***, DeviceType> t_geom_3d;
  typedef Kokkos::View<GeomScalar**[3], DeviceType> t_geom_3d3;
  typedef Kokkos::View<GeomScalar****, DeviceType> t_geom_4d;
  typedef Kokkos::View<NNScalar*, DeviceType> t_nn_1d;
  typedef Kokkos::View<NNScalar**, DeviceType> t_nn_2d;
  typedef Kokkos::View<NNScalar***, DeviceType> t_nn_3d;
  typedef Kokkos::View<NNScalar***, Kokkos::LayoutRight, DeviceType> t_nn_3d_r; // row-major 3D
  typedef Kokkos::View<NNScalar****, DeviceType> t_nn_4d;

  // ---------- Per-chunk geometry arrays (Task 2.3: neighbor + Chebyshev
  // radial + Y_lm precompute). GeomScalar=double for ALL geometry. ----------
  t_int_1d d_ncount;              // [chunk_size]
  t_int_2d d_nearest;             // [chunk_size, maxneigh] local/ghost neighbor index
  t_geom_2d d_rnorms;             // [chunk_size, maxneigh] bond length r
  t_geom_3d3 d_rhats;             // [chunk_size, maxneigh, 3] unit vector (x_j - x_i)/r
  t_int_1d d_mu_i;                // [chunk_size] central-atom model species index
  t_int_2d d_mu_j;                // [chunk_size, maxneigh] neighbor model species index
  t_geom_3d d_radial_basis;       // [chunk_size, maxneigh, nradbase] RAW Chebyshev T_1..T_nradbase (NO envelope)
  t_geom_2d d_env;                // [chunk_size, maxneigh] polynomial cutoff envelope (separate scalar)
  // Per-bond real spherical harmonics, precomputed once per chunk from
  // d_rhats and reused by every downstream SPBF layer. Declared LayoutRight
  // so the innermost (lm) dim is contiguous for readers that sweep it.
  Kokkos::View<GeomScalar***, Kokkos::LayoutRight, DeviceType> d_Y_bond;  // [chunk_size, maxneigh, (lmax+1)^2]

  // ---------- Layer-1 NN intermediates (Task 2.4: scalar-indicator SPBF A1).
  // NNScalar (fp32 in the Mixed default); geometry feeding these kernels
  // (d_radial_basis, d_Y_bond, d_env) is cast down from GeomScalar. ----------
  t_nn_4d d_R1_nl;                // [chunk_size, maxneigh, L1_nradmax, lmax+1] (pre-l_tile-gather)
  t_nn_3d d_A1;                   // [chunk_size, L1_nradmax, (lmax+1)^2]
  t_nn_3d d_A1_2;                 // [chunk_size, rank, nfunc] — first cp_l product (Task 2.5)
  t_nn_3d d_A1_2_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A1_2) (Task 2.6)
  t_nn_3d d_A1_2a;                // [chunk_size, n_out, n_lm] — FCRight2Left(A1_2_red, A1) (Task 2.6)
  // ---- Remaining L1 chain (Task 2.7) ----
  t_nn_3d d_A1_3;                 // [chunk_size, rank, nfunc] — cp_l(A1_2a, A1)
  t_nn_3d d_A1_3_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A1_3)
  t_nn_3d d_A1_2b;                // [chunk_size, n_out, n_lm] — FCRight2Left(A1_2_red, A1)
  t_nn_3d d_A1_4;                 // [chunk_size, rank, nfunc] — cp_l(A1_2b, A1_2b)
  t_nn_3d d_A1_4_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A1_4)
  t_nn_3d d_eq1;                  // [chunk_size, 64, 25] — multi-instr elem-dep ReduceN([A1,A1_2_red,A1_3_red,A1_4_red])
  t_nn_3d d_eq1_norm;             // [chunk_size, 64, 25] — EquivariantRMSNorm(eq1); THE Layer-2 indicator
  t_nn_3d d_rho1;                 // [chunk_size, 32, 1] — multi-instr elem-dep only_invar ReduceN (same 4 inputs as eq1)
  t_nn_3d d_rho1_norm;            // [chunk_size, 32, 1] — InvariantLayerRMSNorm(rho1, only_nonlin)

  // ---------- Layer-2 NN intermediates (Task 2.8: equiv-indicator SPBF A2) ----------
  t_nn_3d d_A2;                   // [chunk_size, A2_n_rad_max=64, A2_nfunc=375] — L2 equiv-ind SPBF basis
  // Which SPBF layer the generic radial-MLP kernel (TagComputeMLPRadial) uses.
  // Set on host before each launch (0 = A1, idx_spbf_A2 = A2); captured by value
  // when *this is copied into the parallel_for functor.
  int radial_mlp_spbf = 0;

  // ---- Remaining L2 chain (Task 2.9). Mirrors the L1 chain exactly EXCEPT
  // every downstream product/FC/reduce here is built on A2_red (a reduced
  // A2), not on the raw SPBF output A2 itself (L1 used raw A1 directly). ----
  t_nn_3d d_A2_red;               // [chunk_size, n_out, n_funcs] — FunctionReduceN({A2})
  t_nn_3d d_A2_2;                 // [chunk_size, rank, nfunc] — cp_l(A2_red, A2_red)
  t_nn_3d d_A2_2_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A2_2)
  t_nn_3d d_A2_2a;                // [chunk_size, n_out, n_lm] — FCRight2Left(A2_2_red, A2_red)
  t_nn_3d d_A2_3;                 // [chunk_size, rank, nfunc] — cp_l(A2_2a, A2_red)
  t_nn_3d d_A2_3_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A2_3)
  t_nn_3d d_A2_2b;                // [chunk_size, n_out, n_lm] — FCRight2Left(A2_2_red, A2_red)
  t_nn_3d d_A2_4;                 // [chunk_size, rank, nfunc] — cp_l(A2_2b, A2_2b)
  t_nn_3d d_A2_4_red;             // [chunk_size, n_out, n_funcs] — FunctionReduceN(A2_4)
  t_nn_3d d_eq2;                  // [chunk_size, n_out, M] — multi-instr elem-dep ReduceN([A2_red,A2_2_red,A2_3_red,A2_4_red])
  t_nn_3d d_eq2_norm;             // [chunk_size, n_out, M] — EquivariantRMSNorm(eq2); THE Layer-3 indicator
  t_nn_3d d_rho2;                 // [chunk_size, n_out, n_lm] — multi-instr elem-dep ReduceN (same 4 inputs as eq2)
  t_nn_3d d_rho2_norm;            // [chunk_size, n_out, n_lm] — InvariantLayerRMSNorm(rho2, FULL branch)

  // ---------- Layer-3 NN intermediates (Task 2.10: equiv-indicator SPBF A3;
  // indicator = eq2_norm). TERMINAL layer: produces rho3/rho3_norm ONLY — there
  // is NO eq3/eq3_norm. The A3_red..rho3 chain is structurally identical to L2's
  // A2_red..rho2 chain (built on A3_red, a reduce of the raw SPBF output A3).
  // A3_3/A3_4 are Lmax=0 (out width ~10), so their dims come from the resolved
  // per-instruction metadata, NOT L2's 50-wide shapes. ----------
  t_nn_3d d_A3;                   // [chunk_size, A3_n_rad_max=64, A3_nfunc=750] — L3 equiv-ind SPBF basis
  t_nn_3d d_A3_red;              // [chunk_size, n_out, n_funcs] — FunctionReduceN({A3})
  t_nn_3d d_A3_2;               // [chunk_size, rank, nfunc] — cp_l(A3_red, A3_red)
  t_nn_3d d_A3_2_red;           // [chunk_size, n_out, n_funcs] — FunctionReduceN(A3_2)
  t_nn_3d d_A3_2a;              // [chunk_size, n_out, n_lm] — FCRight2Left(A3_2_red, A3_red)
  t_nn_3d d_A3_3;               // [chunk_size, rank, nfunc] — cp_l(A3_2a, A3_red) (Lmax=0)
  t_nn_3d d_A3_3_red;           // [chunk_size, n_out, n_funcs] — FunctionReduceN(A3_3)
  t_nn_3d d_A3_2b;              // [chunk_size, n_out, n_lm] — FCRight2Left(A3_2_red, A3_red)
  t_nn_3d d_A3_4;               // [chunk_size, rank, nfunc] — cp_l(A3_2b, A3_2b) (Lmax=0)
  t_nn_3d d_A3_4_red;           // [chunk_size, n_out, n_funcs] — FunctionReduceN(A3_4)
  t_nn_3d d_rho3;               // [chunk_size, 32, 1] — multi-instr elem-dep only_invar ReduceN({A3_red,A3_2_red,A3_3_red,A3_4_red})
  t_nn_3d d_rho3_norm;          // [chunk_size, 32, 1] — InvariantLayerRMSNorm(rho3, FULL branch)

  // ---------- Backward (Task 3.1) adjoint buffers (dE/dX). NNScalar to match
  // the golden fp32-autodiff adjoints (tol 1e-4). Chunk-local (single-chunk
  // regime, chunk_offset=0), mirroring each forward tensor's shape; d_eq2_norm_adj
  // is the ONLY global (nall) one (message adjoint scattered to owned+ghost). ----
  t_nn_3d d_rho1_norm_adj, d_rho2_norm_adj, d_rho3_norm_adj;  // readout bwd (all 3 equal)
  // Task 3.5b: d_rho1_norm_adj/d_rho2_norm_adj are PROMOTED to natom ([nall],
  // allocated in grow_global) so TagReadoutBwd (run over ALL P3 chunks) can seed
  // every owned atom's rho1/rho2 readout adjoint at (ii+chunk_offset), surviving
  // until P5/P4. The shared TagRMSNormBwd reads adj_out chunk-local, so P5/P4
  // gather the natom seed into these chunk-scratch _local buffers per chunk
  // (mirrors d_eq{1,2}_norm_adj_local). d_rho3_norm_adj stays chunk-scratch.
  t_nn_3d d_rho1_norm_adj_local, d_rho2_norm_adj_local;
  t_nn_3d d_rho1_adj, d_rho2_adj, d_rho3_adj;                 // invariant-RMSNorm bwd
  t_nn_3d d_A3_red_adj, d_A3_2_adj, d_A3_2_red_adj, d_A3_2a_adj;
  t_nn_3d d_A3_3_adj, d_A3_3_red_adj, d_A3_2b_adj, d_A3_4_adj, d_A3_4_red_adj;
  t_nn_3d d_A3_adj;             // [chunk_size, A3_n_rad_max, A3_nfunc] — GATE 2
  // Opt 4: shared adj_prod buffer [chunk_size, MAX_NRADMAX, nlm_y*MAX_NLMIND].
  // adj_prod[ii,n,p] = inv_avg*Σ_f A_adj(ii,n,f)*cg_W(p,f) (p = ly*n_lm_ind+lr),
  // computed ONCE per layer by TagComputeAdjProd (output-parallel — no per-thread
  // local array) and read by BOTH the SPBF-bwd (A2/A3) AND ForceEquiv consumers,
  // which previously each recomputed it into a register-spilling 1600-float local.
  t_nn_3d d_adj_prod;
  t_nn_3d d_eq2_norm_adj;       // [nall, I2_n_funcs, I2_n_out] — message adjoint (owned PARTIAL, pre-reverse-comm)
  t_nn_3d d_cpl_lproj_adj, d_cpl_rproj_adj;  // cp_l bwd rank-projection scratch (== d_cpl_lproj/rproj dims)
  // ---- Task 3.2: Layer-2 backward adjoint buffers (mirror the L3 d_A3_* set) ----
  t_nn_3d d_A2_red_adj, d_A2_2_adj, d_A2_2_red_adj, d_A2_2a_adj;
  t_nn_3d d_A2_3_adj, d_A2_3_red_adj, d_A2_2b_adj, d_A2_4_adj, d_A2_4_red_adj;
  t_nn_3d d_A2_adj;             // [chunk_size, A2_n_rad_max, A2_nfunc] — GATE 2 (L2)
  t_nn_3d d_eq2_adj;            // [chunk_size, eq2.n_out, eq2.M] — EquivariantRMSNorm bwd input adj
  // Chunk-local gather of the reverse-comm-completed d_eq2_norm_adj (owned rows,
  // gathered by d_ilist so the generic equiv bwd stays chunk-local). GATE 1 dump src.
  t_nn_3d d_eq2_norm_adj_local; // [chunk_size, I2_n_funcs, I2_n_out]
  // d_eq1_norm_adj (round-1 message adjoint, owned partial) REUSES d_grad_I_global
  // [nall, I_n_funcs, I_n_out] — the round-1 reverse buffer; Task 3.3 reverse-comms it.

  // ---- Task 3.3: Layer-1 backward adjoint buffers. UNLIKE L2/L3, the L1 chain
  // consumes the RAW SPBF basis A1 directly (there is no A1_red), so d_A1_adj is
  // the DIRECT fan-out target (accumulates from eq1[0]+rho1[0] + cp_l A1_2 (both)
  // + fc A1_2a/2b (right) + cp_l A1_3 (right)). ----
  t_nn_3d d_A1_2_adj, d_A1_2_red_adj, d_A1_2a_adj, d_A1_3_adj, d_A1_3_red_adj;
  t_nn_3d d_A1_2b_adj, d_A1_4_adj, d_A1_4_red_adj;
  t_nn_3d d_A1_adj;             // [chunk_size, L1_nradmax, (lmax+1)^2] — GATE 2 (L1)
  t_nn_3d d_eq1_adj;            // d_rho1_adj already declared with the rho{1,2,3}_adj set above
  // Chunk-local gather of the reverse-comm-completed round-1 message adjoint
  // (d_grad_I_global owned rows). GATE 1 dump src for eq1_norm.
  t_nn_3d d_eq1_norm_adj_local; // [chunk_size, I_n_funcs, I_n_out]

  // ---- Task 3.3: per-bond forces + radial r-derivative (geometry backward) ----
  t_geom_3d3 d_f_ij;            // [chunk_size, maxneigh, 3] per-bond force (F_i += f_ij, F_j -= f_ij)
  t_nn_4d d_DR1_nl;             // [chunk_size, maxneigh, L1_nradmax, lmax+1] dR/dr (shared scratch, per layer)
  t_geom_2d d_denv;             // [chunk_size, maxneigh] denv/dr (layer-independent)
  // per-call state for the equivariant force kernel (TagForceEquiv)
  int force_spbf_idx = 0;       // idx_spbf_A2 or idx_spbf_A3
  int force_ind_stage = 1;      // 1 = L2 (indicator d_I_global), 2 = L3 (indicator d_eq2_norm_global)
  int force_n_lm_ind = 0;       // indicator lm width for the current layer
  t_nn_3d force_A_adj;          // d_A2_adj or d_A3_adj (the basis adjoint being backpropped)

#ifdef KOKKOS_ENABLE_CUDA
  // ---- Opt 15: cuBLAS true-fp32 radial-MLP scratch + launch state ----
  cublasHandle_t cublas_handle = nullptr;
  typedef Kokkos::View<NNScalar**, Kokkos::LayoutRight, DeviceType> t_nn_2d_r; // row-major [M x dim]
  // Per-(spbf,layer) weights repacked TIGHTLY row-major [n_in x n_out] at load
  // time. The native mlp_W view is LayoutLeft (column-major) + padded, so a fixed
  // layer's (k,j) block is NOT contiguous -> cannot feed cuBLAS directly. These
  // packed copies give lda=n_out (col-major W^T leading dim). One-time, ~1 MB.
  t_nn_2d_r d_mlp_Wg[GRACE3L_N_SPBF][GRACE3L_MAX_MLP_LAYERS];
  t_nn_2d_r d_mlp_X;    // [M x n_in0=138] value input
  t_nn_2d_r d_mlp_dX;   // [M x nradbase=10] deriv input (only Chebyshev cols are r-dependent)
  t_nn_2d_r d_mlp_h0;   // [M x 256]       layer-0 value activation (raw in-place -> silu)
  t_nn_2d_r d_mlp_dh0;  // [M x 256]       layer-0 deriv activation (raw in-place -> silu'*ds)
  t_nn_2d_r d_mlp_h1;   // [M x 128]       layer-1 value activation
  t_nn_2d_r d_mlp_dh1;  // [M x 128]       layer-1 deriv activation
  t_nn_2d_r d_mlp_R;    // [M x n_out(last)=320] value output GEMM (LayoutRight); scattered
                        // into the LayoutLeft d_R1_nl by MLPScatterR
  t_nn_2d_r d_mlp_dR;   // [M x n_out(last)=320] deriv output GEMM; scattered into d_DR1_nl
  int mlp_M = 0;                 // current bonds = chunk_size * maxneigh
  int mlp_scatter_which = 0;     // MLPScatterR target: 0 = R->d_R1_nl, 1 = DR->d_DR1_nl
  int mlp_act_layer = 0;         // layer index for the activation kernel
  int mlp_act_n = 0;             // n_out width for the activation kernel
  NNScalar mlp_act_norm = 1;     // norm for the activation kernel
  int mlp_assemble_idx = 0;      // spbf idx for the assembly kernel
  // Opt 17: device pointer-array buffer for cublasSgemmBatched (ReduceN). One packed
  // buffer holds [A ptrs | B ptrs | C ptrs], each block sized reduceN_max_conn, filled
  // on host + deep_copied per call (base pointers change only on grow, but rebuilding
  // is cheap vs the GEMMs). uintptr_t storage reinterpreted as float** for cuBLAS.
  Kokkos::View<uintptr_t*, DeviceType> d_redptr;
  Kokkos::View<uintptr_t*, Kokkos::HostSpace> h_redptr;   // explicit HostSpace (portable;
                                                          // ::HostMirror alias trips some nvcc)
  int reduceN_max_conn = 0;
#endif

  // per-call backward launch state (mirrors the forward per-call state members)
  t_nn_3d rmsnorm_adj_in, rmsnorm_adj_out;                 // TagRMSNormBwd
  t_nn_3d eqnorm_adj_out, eqnorm_adj_in;                   // TagEqNormBwd
  t_nn_3d reduceN_adj_out, reduceN_adj_in;                 // TagReduceNBwd
  int reduceN_bwd_instr = 0;                               // current reduce instruction index
  t_nn_3d fc_adj_out, fc_adj_left, fc_adj_right;           // TagFCBwd{Left,Right}
  t_nn_3d cpl_adj_out, cpl_adj_left, cpl_adj_right;        // TagCPLBwd{Couple,Project}
  int cpl_bwd_nmax = 0;                                    // max(n_left,n_right) stride for TagCPLBwdProject
  int cpl_bwd_projected = 0;                               // lproj/rproj supplied by cuBLAS

  // ---------- Global arrays used by the comm subsystem (size = nall) ----------
  t_nn_3d d_I_global;          // [nall, I_n_funcs, I_n_out] — I features (round-1 indicator eq1_norm) for forward comm
  t_nn_3d d_grad_I_global;     // [nall, I_n_funcs, I_n_out] — grad_I for reverse comm
  // Round-2 (L3) indicator: eq2_norm forward-comm'd to ghosts (Task 2.10). Sized
  // [nall, I2_n_funcs, I2_n_out] — parity-doubled lm axis, distinct from d_I_global.
  t_nn_3d d_eq2_norm_global;   // [nall, I2_n_funcs, I2_n_out] — round-2 indicator (eq2_norm)

  // ---- Task 3.5a: natom promotions of d_rho1_norm/d_rho2_norm so the P3
  // readout can read atoms from EARLIER chunks after P1/P2's chunk-scratch
  // d_rho{1,2}_norm has been overwritten by later chunks (multi-chunk forward
  // correctness). UNLIKE d_I_global/d_eq2_norm_global above, these are PURELY
  // LOCAL (no MPI comm — the 3L readout only ever needs OWNED-atom densities),
  // so they are filled by a dedicated per-chunk scatter (ScatterRho1Norm/
  // ScatterRho2Norm) indexed by the dense list-order ii+chunk_offset, NOT
  // d_ilist(ii+chunk_offset). Dims mirror d_rho1_norm/d_rho2_norm's trailing
  // dims exactly; sized to nall in grow_global() (safe superset of inum).
  t_nn_3d d_rho1_norm_full;    // [nall, 32, 1] — natom copy of d_rho1_norm
  t_nn_3d d_rho2_norm_full;    // [nall, n_out, n_lm] — natom copy of d_rho2_norm

  // Host mirrors for MPI communication
  typename t_nn_3d::host_mirror_type h_I_global;
  typename t_nn_3d::host_mirror_type h_grad_I_global;
  typename t_nn_3d::host_mirror_type h_eq2_norm_global;
  typename t_nn_3d::host_mirror_type h_eq2_norm_adj;   // Task 3.2 reverse-comm host fallback

  // ==========================================================================
  // Constant device weight storage (populated by copy_weights_to_device).
  //
  // These are POD-of-Kokkos::View structs held in FIXED C-arrays so the whole
  // set is safe to capture inside a `*this` tag-functor (as GRACE-2L does):
  // Views/scalars/C-arrays-thereof are device-copyable; std::vector would NOT
  // be, so it is deliberately confined to the host GRACE3LModel. Downstream
  // kernels index these arrays positionally (order == the npz *_names lists,
  // also mirrored in grace_model->{spbf,prod,reduce,fc,eqnorm,rmsnorm}); use
  // grace_model to map an instruction name to its array index at launch time.
  //
  // NORM CONVENTION: raw weight + separate `norm` scalar (kernel applies once),
  // EXCEPT SPBF z_tr/z_proj which fold the norm into a per-element precompute.
  // ==========================================================================

  // --- SPBF (A1 scalar, A2/A3 equivariant) ---
  struct DeviceSPBF {
    int equivariant = 0;
    int n_rad_max = 0, n_rad_basis = 0, spbf_lmax = 0, Lmax = 0, p = 0;
    int nfunc = 0, n_lm_ind = 0;        // equivariant only
    int mlp_n_layers = 0, mlp_max_dim = 0;
    GeomScalar rcut = 0, inv_avg_n_neigh = 1;
    t_int_1d l_tile;                    // [nlm]
    // radial MLP: W padded [n_layers, max_in, max_out]; bias padded
    // [n_layers, max_out] (hidden rows filled, output row zero); norms/dims
    t_nn_3d mlp_W;
    t_nn_2d mlp_b;
    t_nn_1d mlp_norms;                  // [n_layers]
    t_int_1d mlp_dims;                  // [n_layers+1]
    // scalar mode (A1): per-element species modulation (norm folded once)
    t_nn_2d z_tr;                       // [n_elements, n_rad_max]
    // equivariant mode (A2/A3): per-element L=0 chem injection (norm folded once)
    t_nn_2d z_proj;                     // [n_elements, n_rad_max]
    t_nn_1d chem_l0_mask;              // [n_lm_ind]
    t_nn_2d cg_W;                      // [(lmax+1)^2 * n_lm_ind, nfunc]
    // Opt 9: sparse forms of cg_W. The CG-coupling matrix is ~99% zeros
    // (mean ~3-6 nnz per row/column), so the dense couple loops streamed
    // mostly zeros. Ascending index order within each column/row preserves
    // the dense loops' accumulation order exactly.
    t_int_1d cg_col_ptr;               // [nfunc+1]  CSC (SPBFFused couple)
    t_int_1d cg_row_idx;               // [nnz]
    t_nn_1d cg_col_val;                // [nnz]
    t_int_1d cg_row_ptr;               // [P+1]      CSR (ComputeAdjProd)
    t_int_1d cg_col_idx;               // [nnz]
    t_nn_1d cg_row_val;                // [nnz]
  };
  DeviceSPBF d_spbf[GRACE3L_N_SPBF];
  int n_spbf = 0;

  // --- cp_l products ---
  struct DeviceCPL {
    int rank = 0, nfunc = 0, n_left = 0, n_right = 0;
    int n_groups_left = 0, n_groups_right = 0, n_cg = 0;
    NNScalar norm_u = 1, norm_v = 1;
    t_nn_3d U;                          // [n_groups_left, rank, n_left]
    t_nn_3d V;                          // [n_groups_right, rank, n_right]
    t_int_1d group_left, group_right;   // [n_lm_left], [n_lm_right]
    t_int_1d left_ind, right_ind, m_sum_ind;   // [n_cg]
    // Opt 14: the three couple indices packed into ONE int32 per entry
    // (m_sum_ind<<16 | right_ind<<8 | left_ind) so the fused couple loops
    // stream 8 B/entry (packed+cg) instead of 16 B. Bounds checked at load.
    t_int_1d cg_packed;                 // [n_cg]
    t_nn_1d cg;                         // [n_cg]
    t_int_1d out_l, out_parity;         // [nfunc]
    // Opt 16: cuBLAS CPLBwdProject scratch. Per-lm-slab (w) group-gathered +
    // transposed U/V, tightly packed so each w-slab is a col-major [rank x n]
    // matrix (ldb=rank, w-stride = n*rank) for cublasSgemmStridedBatched.
    // Uw(w,n,r) = U(group_left(w), r, n); Vw(w,n,r) = V(group_right(w), r, n).
    t_nn_3d_r Uw;                       // [n_lm_left,  n_left,  rank]
    t_nn_3d_r Vw;                       // [n_lm_right, n_right, rank]
  };
  DeviceCPL d_prod[GRACE3L_MAX_PRODS];
  int n_prods = 0;

  // ---- cp_l generic-product scratch + per-call launch state (Task 2.5) ----
  // Rank-axis projections written by TagCPLProject, read by TagCPLCouple. Sized
  // once to the max rank / max n_lm over ALL products so the same two views
  // serve every compute_cp_l call (A1_2, A1_3, ... A3_2). Per-(atom,rank)
  // ownership in Pass B makes the segment-sum serial -> no atomics needed.
  t_nn_3d d_cpl_lproj;                 // [chunk, cpl_scratch_rank, cpl_scratch_nlm]
  t_nn_3d d_cpl_rproj;                 // [chunk, cpl_scratch_rank, cpl_scratch_nlm]
  int cpl_scratch_rank = 0, cpl_scratch_nlm = 0;
  // Captured by the TagCPL* functors for the current compute_cp_l call.
  DeviceCPL cpl_p;
  t_nn_3d cpl_in_left, cpl_in_right, cpl_out;
  int cpl_nlm = 0;                     // max(n_lm_left, n_lm_right) for current product

  // --- Opt 6: fused equiv-SPBF (A2/A3) launch state (captured by TagComputeSPBFFused) ---
  int spbf_fused_idx = -1;             // idx_spbf for the current fused equiv-SPBF launch
  t_nn_3d spbf_fused_indicator;        // per-neighbor indicator source: d_I_global (A2) or d_eq2_norm_global (A3)
  t_nn_3d spbf_fused_out;              // output basis: d_A2 or d_A3

  // --- Opt 10: fused equiv-SPBF backward launch state (captured by TagSPBFBwdFused) ---
  int spbf_bwd_idx = -1;               // idx_spbf_A2 or idx_spbf_A3
  t_nn_3d spbf_bwd_target;             // message-adjoint scatter target: d_grad_I_global (A2) or d_eq2_norm_adj (A3)
  enum { SPBF_BWD_TILE = 32 };         // neighbours staged per shared a_nl tile

  // --- FunctionReduceN ---
  struct DeviceReduce {
    int n_out = 0, n_funcs = 0, n_instr = 0;
    int only_invar = 0, elem_dep = 0;
    int n_in[GRACE3L_MAX_REDUCE_INSTR];
    int w_shape[GRACE3L_MAX_REDUCE_INSTR];
    int n_conn[GRACE3L_MAX_REDUCE_INSTR];
    NNScalar norm[GRACE3L_MAX_REDUCE_INSTR];
    // W stored uniformly 4D: [n_types, n_out, n_in, w_shape] with
    // n_types = elem_dep ? n_elements : 1  (kernel uses e = elem_dep?mu_i:0).
    t_nn_4d W[GRACE3L_MAX_REDUCE_INSTR];
    t_int_1d collect_ind[GRACE3L_MAX_REDUCE_INSTR];   // [n_conn]
    t_int_1d w_l_tile[GRACE3L_MAX_REDUCE_INSTR];      // [n_conn]
    t_int_1d total_sum_ind[GRACE3L_MAX_REDUCE_INSTR]; // [n_conn]
    t_nn_1d norm_map;                   // [n_funcs] (or [1] for only_invar)
    int has_norm_map = 0;
    // Opt 17: cuBLAS batched-GEMM scratch for elem-indep, non-only_invar reduces
    // (n_instr==1). Wsc[c] = norm*norm_map(tsi(c)) * W(0,:,:,wlt(c)), stored as a
    // col-major [n_out x n_in] block per connection (c-stride n_out*n_in, ld n_out).
    // Serves BOTH forward (rounds; alpha=1,OP_N/OP_T) and backward (single batch;
    // alpha=1,OP_N/OP_N,beta=1) — the norm folding is identical for both.
    t_nn_1d Wsc;                        // [n_conn * n_out * n_in]  (instr 0 only)
  };
  DeviceReduce d_reduce[GRACE3L_MAX_REDUCES];
  int n_reduces = 0;
  // Opt 17: host-side per-reduce metadata for the cuBLAS batched path (index-parallel
  // to d_reduce). Holds the host copies of the connection indices + the forward round
  // schedule used to build the batched-GEMM device pointer arrays each call.
  struct ReduceHost {
    bool cublas_ok = false;            // elem_dep==0 && only_invar==0 && n_instr==1
    int n_conn = 0, n_out = 0, n_in = 0;
    std::vector<int> ci, tsi;          // [n_conn] collect_ind / total_sum_ind (host)
    std::vector<int> fwd_order;        // [n_conn] connections grouped by round
    int nrounds = 0;
    int round_ptr[GRACE3L_MAX_REDUCE_ROUNDS + 1] = {0};
  };
  ReduceHost h_reduce[GRACE3L_MAX_REDUCES];

  // --- FCRight2Left ---
  struct DeviceFC {
    int n_out = 0, left_coefs = 0;
    int n_funcs_left = 0, n_funcs_right = 0;
    int n_in_left = 0, n_in_right = 0, w_shape_left = 0, w_shape_right = 0;
    NNScalar norm_left = 1, norm_right = 1;
    t_nn_3d w_left;                     // [n_out, n_in_left, w_shape_left] (empty if !left_coefs)
    t_nn_3d w_right;                    // [n_out, n_in_right, w_shape_right]
    t_int_1d w_tile_left, w_tile_right; // [n_funcs_left], [n_funcs_right]
    t_int_1d collect_to, collect_from;  // [n_funcs_left], [n_funcs_right]
    t_nn_1d norm_out_factor;            // [n_funcs_left]
  };
  DeviceFC d_fc[GRACE3L_MAX_FCS];
  int n_fcs = 0;

  // ---- FunctionReduceN generic-kernel scratch + per-call launch state (Task 2.6).
  // Mirrors the cp_l pattern (DeviceCPL/cpl_p above): caller supplies the
  // per-instruction input views by name-matching grace_model->reduce[i].instr[j]
  // .name at the call site; the device kernel itself never does string lookup. ----
  DeviceReduce reduceN_r;
  t_nn_3d reduceN_inputs[GRACE3L_MAX_REDUCE_INSTR];
  t_nn_3d reduceN_out;

  // ---- FCRight2Left per-call launch state (Task 2.6) ----
  DeviceFC fc_p;
  t_nn_3d fc_in_left, fc_in_right, fc_out;

  // --- EquivariantRMSNorm ---
  struct DeviceEqNorm {
    int center_l0 = 0, n_groups = 0, n_out = 0, M = 0;
    NNScalar eps = 0;
    t_nn_2d affine_weight;              // [n_groups, n_out]
    t_nn_1d degree_weights;            // [M]
    t_int_1d expand_index;             // [M]
    t_nn_1d l0_mask;                   // [M]
  };
  DeviceEqNorm d_eqnorm[GRACE3L_MAX_EQNORMS];
  int n_eqnorms = 0;

  // --- InvariantLayerRMSNorm ---
  struct DeviceRMSNorm {
    int type = 0, n_out = 0, scale_len = 0;
    t_nn_1d scale;                     // [scale_len]
  };
  DeviceRMSNorm d_rmsnorm[GRACE3L_MAX_RMSNORMS];
  int n_rmsnorms = 0;

  // ---- EquivariantRMSNorm per-call launch state (Task 2.7; same pattern as
  // reduceN_r/fc_p above: caller stashes metadata + in/out views, kernel is
  // pure tag-functor code with no string lookups). ----
  DeviceEqNorm eqnorm_p;
  t_nn_3d eqnorm_in, eqnorm_out;

  // ---- InvariantLayerRMSNorm per-call launch state (Task 2.7) ----
  DeviceRMSNorm rmsnorm_p;
  t_nn_3d rmsnorm_in, rmsnorm_out;

  // --- Chemical embedding / shifts / readout ---
  t_nn_2d d_chem_embed;                 // [n_elements, embedding_size]
  t_nn_1d d_shifts;                     // [n_elements]
  NNScalar output_scale = 1;

  // Energy readout MLP (no bias): W padded [n_layers, max_in, max_out]
  t_nn_3d d_energy_W;
  t_nn_1d d_energy_norms;               // [n_layers]
  t_int_1d d_energy_dims;               // [n_layers+1]
  int energy_n_layers = 0;
  int energy_max_dim = 0;
  int energy_activation = 0;            // 0=silu, 1=tanh

  // Per-atom energy (Task 2.11). Deliberately KK_ACC_FLOAT (fp64 by default),
  // NOT NNScalar/GeomScalar: the MLP math runs in NNScalar (matching TF's fp32
  // readout), but the accumulated per-atom result — and the eng_vdwl reduction
  // over it — must stay fp64 regardless of which PairGRACE3LKokkos<> NNScalar/
  // GeomScalar instantiation (Mixed or FP32) is active. Chunk-local (size
  // natom, like GRACE-2L's d_e_atom); indexed by ii exactly like rho3_norm.
  Kokkos::View<KK_ACC_FLOAT*, DeviceType> d_e_atom;

  // Per-element-pair cutoffs
  t_geom_2d d_bond_cutoff;              // [n_elements, n_elements]

  // Spherical harmonics precomputed coefficients (verbatim from GRACE-2L)
  t_geom_1d d_idx_sph, alm, blm, cl, dl;
  int idx_sph_max = 0;
  t_int_1d d_l_from_lm;                 // [(lmax+1)^2] -> l

  // ---------- Helper functions ----------
  void grow(int natom, int maxneigh);
  // Allocate the global (nall-sized) inter-layer comm arrays (Task 2.8). Sized
  // to I_n_funcs/I_n_out (set in init_style); mirrors GRACE-2L's grow_global.
  void grow_global(int nall);
  // Generic cp_l GeneralProductFunction product (Task 2.5): out[a,r,f] from
  // left/right [a,n,lm] for product p over the current chunk. Reused for every
  // A*_* cp_l product; caller supplies the per-product weights + in/out views.
  void compute_cp_l(const DeviceCPL &p, const t_nn_3d &in_left,
                    const t_nn_3d &in_right, const t_nn_3d &out);
  // Opt 6: fused equiv-SPBF (A2/A3) launcher — one team per (ii,n); the prod
  // slice + dense CG couple live in team shared (replaces the register-spilling
  // TagComputeA2/TagComputeA3). idx picks the SPBF weights; indicator is the
  // per-neighbor message source (d_I_global for A2, d_eq2_norm_global for A3);
  // out is the basis output (d_A2 or d_A3).
  void compute_spbf_equiv(int idx, const t_nn_3d &indicator, const t_nn_3d &out);
  // Generic FunctionReduceN (Task 2.6): out[a,k,f] from per-instruction gathered
  // inputs for reduce r over the current chunk. Handles multi-instruction sums,
  // elem_dep W indexing (W stored uniformly 4D [n_types,n_out,n_in,w_shape],
  // n_types=1 if !elem_dep, so the kernel always indexes e=elem_dep?mu_i:0 with
  // no branch on the W read), only_invar collapse (every connection accumulates
  // into f=0, ignoring total_sum_ind) vs the general scatter-by-total_sum_ind,
  // and an optional norm_map applied once at the end. Reused for every reduce
  // (A*_*_red, eq*, rho*, ...); caller supplies inputs[j] matching instr[j].name.
  void compute_reduceN(const DeviceReduce &r, const t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR],
                       const t_nn_3d &out);
  // Generic FCRight2Left (Task 2.6): out[a,k,lm] = left_t[lm,a,k] with right's
  // per-connection contribution scattered in via collect_from/collect_to.
  // Handles left_coefs=True (left projected via w_left/w_tile_left, like the
  // right branch) and left_coefs=False (left already at this FC's n_out, just
  // transposed) generically; always applies norm_out_factor (the loader
  // requires it unconditionally, matching the oracle's `if W.has(...)` being
  // always-true for this exporter). Reused for every FC (A*_2a/2b, ...);
  // caller supplies the left/right input views.
  void compute_fc(const DeviceFC &f, const t_nn_3d &in_left, const t_nn_3d &in_right,
                  const t_nn_3d &out);
  // Generic EquivariantRMSNorm (Task 2.7): x=[atoms,n_feat=n_out,n_lm=M].
  // center_l0 (per-(atom,lm) cross-feature mean, masked to l==0 slots) ->
  // ONE degree-balanced RMS scalar per atom (mean over feat of sum over lm of
  // dw[lm]*x^2) -> per-(feat,lm) affine scale = affine_weight[expand_index[lm], feat].
  // Reused for eq1_norm/eq2_norm.
  void compute_equiv_rms_norm(const DeviceEqNorm &e, const t_nn_3d &in, const t_nn_3d &out);
  // Generic InvariantLayerRMSNorm (Task 2.7): x=[atoms,n_out,n_lm]. Branches
  // on scale_len vs n_out exactly like the oracle: full (scale_len==n_out) ->
  // every channel scaled by rsqrt(mean_over_channel(x^2)+eps); only_nonlin
  // (scale_len==n_out-1) -> channel 0 passthrough, channels 1..n_out-1 scaled
  // by rsqrt(mean_over_the_nonlin_channels(x^2)+eps). eps=1e-10 fixed (not
  // loaded from the npz — matches the oracle's hardcoded F32(1e-10)). Reused
  // for rho1_norm/rho2_norm/rho3_norm.
  void compute_invariant_rms_norm(const DeviceRMSNorm &r, const t_nn_3d &in, const t_nn_3d &out);

  // ---- Task 3.1: backward (reverse-mode adjoint) routines ----
  // Readout backward: seeds d_rho{1,2,3}_norm_adj = dE/d(rho_k_norm) (all equal,
  // since the 3 densities enter the readout symmetrically). ch0=output_scale
  // (linear skip), ch1..31 from the MLP backward; every lm>0 slot = 0.
  void readout_bwd();
  // Generic InvariantLayerRMSNorm VJP (mirrors compute_invariant_rms_norm): given
  // the forward input `in`, forward output `out` and adj_out, writes adj_in
  // (assignment; the rho densities have a single consumer). Handles both branches
  // (full / only_nonlin) data-driven off scale_len, eps=1e-10.
  void compute_invariant_rms_norm_bwd(const DeviceRMSNorm &r, const t_nn_3d &in,
      const t_nn_3d &out, const t_nn_3d &adj_out, const t_nn_3d &adj_in);
  // Generic FunctionReduceN VJP (mirrors compute_reduceN): scatters adj_out into
  // each per-instruction input adjoint via the transpose contraction, applying
  // norm + norm_map. ACCUMULATES (+=) into adj_inputs[j] (inputs fan out to many
  // reduces/products) — caller must zero the adjoint buffers first.
  void compute_reduceN_bwd(const DeviceReduce &r, const t_nn_3d &adj_out,
      const t_nn_3d adj_inputs[GRACE3L_MAX_REDUCE_INSTR]);
  // Generic FCRight2Left VJP (mirrors compute_fc): ACCUMULATES (+=) into adj_left
  // (identity/left_coefs path) and adj_right (scatter via collect_from/to, atomic).
  void compute_fc_bwd(const DeviceFC &f, const t_nn_3d &adj_out,
      const t_nn_3d &adj_left, const t_nn_3d &adj_right);
  // Generic cp_l VJP (product rule; spec §4.6). Recomputes the forward Pass-A
  // projections (shared scratch is overwritten per call), backprops the CG couple
  // to the rank-projection adjoints, then the U/V projection to adj_left/adj_right.
  // ACCUMULATES (+=). Self-products pass the same buffer for adj_left==adj_right;
  // both contributions accumulate correctly (serial per (a,n,w) thread).
  void compute_cp_l_bwd(const DeviceCPL &p, const t_nn_3d &in_left,
      const t_nn_3d &in_right, const t_nn_3d &adj_out,
      const t_nn_3d &adj_left, const t_nn_3d &adj_right);
  // L3 equivariant-SPBF backward: scatters d_A3_adj into d_eq2_norm_adj (the L3
  // message adjoint) over every (owned atom i, radial n, neighbor j). j may be
  // owned OR ghost -> d_eq2_norm_adj is nall-sized, atomic-add scatter. The
  // species-injection (z_proj/chem_l0_mask) term is constant and contributes
  // nothing to the indicator adjoint. NOT reverse-comm'd here (Task 3.2).
  void compute_A3_spbf_bwd();
  // Opt 4: fill d_adj_prod for the given layer (spbf_idx = idx_spbf_A2 or _A3,
  // A_adj = d_A2_adj or d_A3_adj) BEFORE its SPBF-bwd + ForceEquiv consumers.
  void compute_adj_prod(int spbf_idx, const t_nn_3d &A_adj);
  // Task 3.2: Generic EquivariantRMSNorm VJP (mirrors compute_equiv_rms_norm).
  // Given the forward input `in` and adj_out, ASSIGNS adj_in (single consumer:
  // eq2 feeds only eq2_norm). Recomputes center_l0 mean + degree-balanced rms
  // from `in` (the forward caches neither). Reused for eq1_norm in Task 3.3.
  void compute_equiv_rms_norm_bwd(const DeviceEqNorm &e, const t_nn_3d &in,
      const t_nn_3d &adj_out, const t_nn_3d &adj_in);
  // Task 3.2: L2 equivariant-SPBF backward: scatters d_A2_adj into d_grad_I_global
  // (the round-1 message adjoint d_eq1_norm_adj, owned partial) over (owned atom i,
  // radial n, neighbor j; j owned OR ghost -> atomic scatter). The z_proj/
  // chem_l0_mask species-injection term is constant -> no indicator adj. NOT
  // reverse-comm'd here (Task 3.3). Requires d_R1_nl to hold A2's radial (caller
  // recomputes it — A3's forward overwrote the shared d_R1_nl scratch).
  void compute_A2_spbf_bwd();
  // Opt 10: shared launcher for the fused team+shared equiv-SPBF backward
  // (TagSPBFBwdFused). idx = idx_spbf_A2 or idx_spbf_A3; adj_target = the
  // nall-sized message-adjoint buffer the kernel atomic-add scatters into
  // (d_grad_I_global for A2, d_eq2_norm_adj for A3).
  void compute_spbf_bwd_fused(int idx, const t_nn_3d &adj_target);
  // Task 3.3: radial MLP forward + its r-derivative for the layer selected by
  // radial_mlp_spbf. Writes d_R1_nl (R, matching the forward), d_DR1_nl (dR/dr),
  // and d_denv (denv/dr, layer-independent). Bit-identical R to TagComputeMLPRadial.
  void compute_mlp_radial_deriv();
  // Opt 15: value-only radial MLP forward for layer `radial_mlp_spbf`, writing
  // d_R1_nl. Dispatches to cuBLAS SGEMM (CUDA) or the hand kernel (fallback).
  // Replaces the inline TagComputeMLPRadial launches in l{1,2}_forward_chunk.
  void compute_mlp_radial();
#ifdef KOKKOS_ENABLE_CUDA
  // Opt 15 helper: C[M x n_out] = alpha * In[M x n_in] * Wg (true-fp32 SGEMM,
  // CUBLAS_PEDANTIC_MATH). Row-major operands; ld* are the physical row strides
  // (allocated widths): ldW=n_out for the tight Wg, ldIn=In's alloc width,
  // ldC=C's alloc width (e.g. L1_nradmax*(lmax+1) when C aliases d_R1_nl).
  void mlp_sgemm(const NNScalar* Wg, const NNScalar* In, NNScalar* C,
      int n_out, int M, int n_in, int ldW, int ldIn, int ldC, NNScalar alpha);
  // Opt 15: one-time layer-0 transpose/leading-dim self-test vs a host reference
  // (aborts on mismatch). Called from init_style before the cuBLAS path is used.
  void mlp_cublas_selftest();
  // Opt 16 helper: batched projection SGEMM for TagCPLBwdProject. For each lm-slab
  // w in [0,batch): C(:,:,w)[M x n] += alpha * A(:,:,w)[M x rank] * Bw[w][rank x n].
  // A/C are LayoutLeft slabs (col-major, ld = extent(0), w-stride = extent(0)*extent(1));
  // Bw is the tightly-packed Uw/Vw (col-major [rank x n], ldb=rank, w-stride=n*rank).
  // beta=1 (accumulates into the shared adj buffer). True-fp32 (PEDANTIC handle).
  void cpl_project_sgemm_batched(const NNScalar* A, long strideA, int ldA,
      const NNScalar* Bw, int n, int rank, NNScalar* C, long strideC, int ldC,
      int M, int batch, NNScalar alpha);
  // Opt 17: cuBLAS batched-GEMM ReduceN forward/backward (elem-indep, non-only_invar).
  // Return true if handled via cuBLAS; false -> caller runs the hand-kernel fallback.
  bool compute_reduceN_cublas(const DeviceReduce &r, int idx,
      const t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR], const t_nn_3d &out);
  bool compute_reduceN_bwd_cublas(const DeviceReduce &r, int idx,
      const t_nn_3d &adj_out, const t_nn_3d &adj_in);
  // Opt 18 (fp32 only): cuBLAS U/V rank-projection for cp_l forward. Writes the global
  // lproj/rproj scratch (d_cpl_lproj/d_cpl_rproj) via two strided-batched GEMMs (reusing
  // the Opt-16 Uw/Vw weights), replacing the shared-memory projection in TagCPLFused;
  // the couple then reads lproj/rproj back (TagCPLCouple). fp64 + non-CUDA keep the fused
  // kernel (no SGEMM benefit, and un-fusing reintroduces the round-trip Opt-5 removed).
  void compute_cp_l_project_cublas(const DeviceCPL &p,
      const t_nn_3d &in_left, const t_nn_3d &in_right);
#endif
  // Task 3.3: per-bond geometry-force kernels. Each ACCUMULATES (atomic +=) into
  // d_f_ij; caller zeros d_f_ij once before invoking all three layers.
  void compute_force_L1();
  void compute_force_equiv(int spbf_idx, const t_nn_3d &A_adj, int ind_stage, int n_lm_ind);
  void copy_weights_to_device();
  void precompute_harmonics();
  void allocate();

 public:
  // Public because it contains an extended __host__ __device__ lambda; nvcc
  // forbids those inside private/protected members.
  // Fills d_Y_bond from d_rhats for the current chunk (uses member
  // chunk_size, maxneigh). Must be (re)called whenever d_rhats for the
  // current chunk changes before any kernel reads d_Y_bond.
  void compute_Y_bond_chunk();
  // Task 3.5b: forward-chunk recompute helpers (the L1/L2 per-chunk forward
  // bodies extracted from the P1/P2 loops, WITHOUT the global eq-scatter) and
  // the per-chunk force scatter (extracted from the old single-shot P6b).
  // P1 calls l1_forward_chunk() then the eq1_norm scatter; P5 (L1 backward)
  // re-runs l1_forward_chunk() to repopulate the chunk-scratch forward
  // intermediates. Likewise l2_forward_chunk() for P2/P4. scatter_forces_chunk()
  // scatters the current chunk's d_f_ij into f[] + virial. Public for the same
  // reason as compute_Y_bond_chunk() (they contain extended device lambdas).
  void l1_forward_chunk();
  void l2_forward_chunk();
  void scatter_forces_chunk();
 protected:

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template <typename scratch_type>
  KOKKOS_INLINE_FUNCTION
  int scratch_size_helper(int values_per_team) const;

  // Activation functions (used by the NN kernels added in later tasks)
  KOKKOS_INLINE_FUNCTION
  static NNScalar silu(NNScalar x) {
    return x / (NNScalar(1.0) + Kokkos::exp(-x));
  }

  KOKKOS_INLINE_FUNCTION
  static NNScalar tanh_act(NNScalar x) {
    return Kokkos::tanh(x);
  }

  // Energy-only mode
  bool debug_no_energy_only_calc = false;

  // Comm helper views for Kokkos-native communication
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d v_buf;
  int comm_first;  // for unpack: starting ghost index

  // Host-side model data
  GRACE3LModel *grace_model;
  std::vector<int> h_type2model;

  // Task 2.6: reduce/FC array indices resolved by name-match at init_style
  // time (index-parallel with d_reduce/d_fc — see copy_weights_to_device,
  // which fills d_reduce[i]/d_fc[i] in the same order as grace_model->reduce
  // /fc). Set to -1 until resolved; init_style() errors out if not found.
  int idx_reduce_A1_2_red = -1;
  int idx_fc_A1_2a = -1;

  // Task 2.7: remaining L1 prod/reduce/FC/eqnorm/rmsnorm indices, resolved by
  // NAME-match at init_style time (same convention as idx_reduce_A1_2_red).
  int idx_prod_A1_3 = -1;
  int idx_prod_A1_4 = -1;
  int idx_reduce_A1_3_red = -1;
  int idx_fc_A1_2b = -1;
  int idx_reduce_A1_4_red = -1;
  int idx_reduce_eq1 = -1;
  int idx_reduce_rho1 = -1;
  int idx_eqnorm_eq1_norm = -1;
  int idx_rmsnorm_rho1_norm = -1;
  // Task 2.8: SPBF layer index of "A2" (equivariant), resolved by name in
  // init_style (spbf order is [A1,A2,A3]; index 1, but resolved defensively).
  int idx_spbf_A2 = -1;
  // eq1/rho1 both collect the same 4 named L1 tensors {A1, A1_2_red,
  // A1_3_red, A1_4_red}; role[j] tells compute() which of those 4 views to
  // plug into inputs[j] for instruction j, resolved by matching
  // grace_model->reduce[idx]->instr[j].name (0=A1,1=A1_2_red,2=A1_3_red,
  // 3=A1_4_red) — NOT assumed positional, per the loaded instruction_names.
  int eq1_input_role[GRACE3L_MAX_REDUCE_INSTR];
  int rho1_input_role[GRACE3L_MAX_REDUCE_INSTR];

  // Task 2.9: remaining L2 prod/reduce/FC/eqnorm/rmsnorm indices, resolved by
  // NAME-match at init_style time (same convention as the Task 2.7 L1 idx_*
  // members above). idx_reduce_A2_red is the extra reduce L2 has that L1
  // does not (L1's products consume raw A1; L2's consume A2_red).
  int idx_reduce_A2_red = -1;
  int idx_prod_A2_2 = -1;
  int idx_reduce_A2_2_red = -1;
  int idx_fc_A2_2a = -1;
  int idx_prod_A2_3 = -1;
  int idx_reduce_A2_3_red = -1;
  int idx_fc_A2_2b = -1;
  int idx_prod_A2_4 = -1;
  int idx_reduce_A2_4_red = -1;
  int idx_reduce_eq2 = -1;
  int idx_reduce_rho2 = -1;
  int idx_eqnorm_eq2_norm = -1;
  int idx_rmsnorm_rho2_norm = -1;
  // eq2/rho2 both collect the same 4 named L2 tensors {A2_red, A2_2_red,
  // A2_3_red, A2_4_red}; role[j] tells compute() which of those 4 views to
  // plug into inputs[j] for instruction j, resolved by matching
  // grace_model->reduce[idx]->instr[j].name (0=A2_red,1=A2_2_red,2=A2_3_red,
  // 3=A2_4_red) — NOT assumed positional, per the loaded instruction_names.
  int eq2_input_role[GRACE3L_MAX_REDUCE_INSTR];
  int rho2_input_role[GRACE3L_MAX_REDUCE_INSTR];

  // Task 2.10: L3 SPBF + remaining L3 prod/reduce/FC/rmsnorm indices, resolved
  // by NAME-match at init_style (same convention as the Task 2.9 L2 idx_*
  // members). L3 is TERMINAL: it has A3_red (like A2_red) and the same
  // cp_l/reduce/fc chain, but produces ONLY rho3/rho3_norm — there is NO eq3.
  int idx_spbf_A3 = -1;
  int idx_reduce_A3_red = -1;
  int idx_prod_A3_2 = -1;
  int idx_reduce_A3_2_red = -1;
  int idx_fc_A3_2a = -1;
  int idx_prod_A3_3 = -1;
  int idx_reduce_A3_3_red = -1;
  int idx_fc_A3_2b = -1;
  int idx_prod_A3_4 = -1;
  int idx_reduce_A3_4_red = -1;
  int idx_reduce_rho3 = -1;
  int idx_rmsnorm_rho3_norm = -1;
  // rho3 collects the same 4 named L3 tensors {A3_red, A3_2_red, A3_3_red,
  // A3_4_red}; role[j] maps instruction j to one of those 4 views, resolved by
  // matching grace_model->reduce[idx_reduce_rho3].instr[j].name (NOT positional).
  // Only rho3 needs a role array here — there is no eq3.
  int rho3_input_role[GRACE3L_MAX_REDUCE_INSTR];

  // ---------- UQ / extrapolation grade (forward only) ----------
  // FP64 device storage regardless of NNScalar (numerical safety). The GMM feature
  // is the basis-RP projection of the three per-layer scalar reduces rho1/rho2/rho3
  // (R-row order [rho1, rho2, rho3]; density channels [full, rho1, rho2, rho3]).
  typedef Kokkos::View<double*, DeviceType>    t_uq_1d;
  // LayoutRight so the last (contiguous) index — split/reduced across team lanes in
  // the team-parallel ComputeUQ kernels — gives coalesced lane reads (rp_matrix over
  // rp_dim, centroids/inv_cov over the feature dims, d_uq_z over rp_dim). Transparent
  // to the logical-index copy_* loaders.
  typedef Kokkos::View<double**, Kokkos::LayoutRight, DeviceType>   t_uq_2d;
  typedef Kokkos::View<double***, Kokkos::LayoutRight, DeviceType>  t_uq_3d;
  typedef Kokkos::View<double****, Kokkos::LayoutRight, DeviceType> t_uq_4d;
  t_uq_3d d_uq_centroids;          // [E, Kmax, D]
  t_uq_4d d_uq_inv_cov;            // [E, Kmax, D, D]
  t_int_1d d_uq_n_clusters;        // [E]
  t_uq_2d d_uq_interp_thresholds;  // [E, Kmax]
  t_uq_2d d_uq_rp_matrix;          // [D_basis, rp_dim] projection R
  t_uq_2d d_uq_z;                  // [nlocal, rp_dim] per-atom RAW proj accumulated across rho1 (P1) + rho2 (P2), finished in P3
  t_uq_1d d_uq_n2_rho1;            // [nlocal] ||B^(rho1)||^2 from Phase 1 (block density carry)
  t_uq_1d d_uq_n2_rho2;            // [nlocal] ||B^(rho2)||^2 from Phase 2 (block density carry)
  t_uq_1d d_gamma, d_sigma, d_gmm_cluster;

  bool has_uq = false;
  int uq_Kmax = 0, uq_D = 0;          // uq_D = full GMM feature dim (rp_dim + n_density)
  int uq_rp_dim = 0;                  // projection output width (R cols; e.g. 128)
  int uq_n_density = 0;               // density channels appended after proj (uq_D - rp_dim)
  int uq_normalize = 0;               // 1 -> L2-normalize basis before projection (v6)
  double uq_density_scale = 1.0;      // density_scale
  int uq_d_basis = 0;                 // total basis width (rows of R)
  int uq_d_basis_rho1 = 0;            // rho1 block width (R rows [0, rho1))
  int uq_d_basis_rho2 = 0;            // rho2 block width (R rows [rho1, rho1+rho2))
  int uq_d_basis_rho3 = 0;            // rho3 block width (R rows [rho1+rho2, d_basis))
  int flag_compute_gamma = 0, flag_compute_atomic_sigma = 0, flag_compute_gmm_cluster = 0;
  int nmax_uq = 0;
  double *gamma = nullptr, *atomic_sigma = nullptr, *gmm_cluster = nullptr;
};

using PairGRACE3LKokkos_Mixed_Device = PairGRACE3LKokkos<LMPDeviceType, float, double>;
using PairGRACE3LKokkos_FP32_Device  = PairGRACE3LKokkos<LMPDeviceType, float, float>;

}    // namespace LAMMPS_NS

#endif
#endif
