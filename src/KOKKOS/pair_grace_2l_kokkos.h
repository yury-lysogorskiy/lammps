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
   GRACE-2L KOKKOS implementation
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/2l/kk,              PairGRACE2LKokkos_FP64_Device);
PairStyle(grace/2l/kk/device,       PairGRACE2LKokkos_FP64_Device);
PairStyle(grace/2l/kk/host,         PairGRACE2LKokkos_FP64_Host);
PairStyle(grace/2l/kk/mixed,        PairGRACE2LKokkos_Mixed_Device);
PairStyle(grace/2l/kk/mixed/device, PairGRACE2LKokkos_Mixed_Device);
PairStyle(grace/2l/kk/mixed/host,   PairGRACE2LKokkos_Mixed_Host);
PairStyle(grace/2l/kk/fp32,         PairGRACE2LKokkos_FP32_Device);
PairStyle(grace/2l/kk/fp32/device,  PairGRACE2LKokkos_FP32_Device);
PairStyle(grace/2l/kk/fp32/host,    PairGRACE2LKokkos_FP32_Host);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRACE_2L_KOKKOS_H
#define LMP_PAIR_GRACE_2L_KOKKOS_H

#include "pair.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "pair_kokkos.h"

// Opt-2L-15: cuBLAS true-fp32 SGEMM for the radial MLP (CUDA backend only).
#ifdef KOKKOS_ENABLE_CUDA
#include <cublas_v2.h>
#endif

namespace LAMMPS_NS {

// ============================================================================
// GRACE2LModel: Host-side weight container loaded from .npz file
// ============================================================================
struct GRACE2LModel {
  // Architecture metadata
  int n_elements = 0;
  int embedding_size = 0;
  int lmax = 0;
  int nradbase = 0;
  int radial_basis_p = 0;
  double rcut = 0.0;
  bool has_bond_specific_cutoff = false;
  std::vector<double> bond_cutoff_map;  // [n_elements * n_elements]

  // Chemical embedding [n_elements, embedding_size]
  std::vector<double> chem_embedding;

  // MLP layer struct (reused for all MLPs)
  struct MLPLayer {
    std::vector<double> W;
    double norm = 1.0;
    int n_in = 0, n_out = 0;
  };

  // ---- Layer 1 ----
  int L1_nradmax = 42;
  int L1_mlp_rad_n_layers = 0;
  bool L1_mlp_rad_has_chem_emb = false;
  std::vector<MLPLayer> L1_mlp_rad_layers;  // R

  // A indicator transform
  std::vector<double> A_lin_transform_W;  // [embedding_size, L1_nradmax]
  double A_lin_transform_norm = 1.0;
  double A_inv_avg_n_neigh = 1.0;

  // FC weights
  struct FCWeights {
    std::vector<double> w_left, w_right;
    std::vector<int> w_tile_left, w_tile_right;
    std::vector<int> collect_to, collect_from;
    std::vector<double> norm_out_factor;
    double norm_left = 1.0, norm_right = 1.0;
    int n_out = 0;
    bool left_coefs = true;
    int n_funcs_left = 0, n_funcs_right = 0;
    int w_shape_left = 0, w_shape_right = 0;
  };
  FCWeights fc_A1, fc_AA1, fc_AA2;
  int L1_fc_n_out = 42;

  // CG products
  struct CGProduct {
    std::vector<int> left_ind, right_ind, m_sum_ind;
    std::vector<double> cg_coeff;
    int n_output_funcs = 0;
    int n_cg_terms = 0;
  };
  CGProduct prod_AA, prod_AAA, prod_AAAA;

  // Generalized ReduceN weights (supports equivariant output)
  struct ReduceWeights {
    std::vector<double> W;  // element-dep: [n_elem, n_out, n_in, w_shape]; non-elem: [n_out, n_in, w_shape]
    double norm = 1.0;
    int n_in = 0, w_shape = 0;
    int n_connections = 0;         // len(collect_ind) — may differ from w_shape for equivariant
    std::vector<int> collect_ind;  // [n_connections] input function indices
    std::vector<int> w_l_tile;     // [n_connections] tile index into W's last dim
    std::vector<int> total_sum_ind; // [n_connections] output function index
  };

  // I1 reducer (equivariant, element-dependent)
  ReduceWeights I1_reduce_A, I1_reduce_AA, I1_reduce_AAA, I1_reduce_AAAA;
  int I1_n_out = 12;
  int I1_n_funcs = 4;
  bool I1_is_elem_dependent = true;

  // I reducer (equivariant, NOT element-dependent)
  ReduceWeights I_reduce_I1;
  int I_n_out = 32;
  int I_n_funcs = 4;
  bool I_is_elem_dependent = false;

  // rho reducer (scalar, element-dependent) — same structure as 1L
  ReduceWeights rho_reduce_A, rho_reduce_AA, rho_reduce_AAA, rho_reduce_AAAA;
  int rho_n_out = 17;
  bool rho_is_elem_dependent = true;

  // RMSNorm: I_nl_LN (only_nonlin: scale has n_out-1 elements, rho[0] passed through)
  std::vector<double> I_nl_LN_scale;  // [rho_n_out - 1] = [16]
  int I_nl_LN_type = 1;  // 1 = only_nonlin
  int I_nl_LN_n_out = 17;

  // ---- Layer 2 ----
  int L2_nradmax = 32;
  int L2_mlp_rad_n_layers = 0;
  std::vector<MLPLayer> L2_mlp_rad_layers;  // R1

  // B0 indicator transform
  std::vector<double> B0_lin_transform_W;  // [embedding_size, L2_nradmax]
  double B0_lin_transform_norm = 1.0;
  double B0_inv_avg_n_neigh = 1.0;

  // YI equivariant indicator CG coupling
  struct YICoupling {
    std::vector<int> lr_inds;     // [n_cg_terms * 2] flattened (y_idx, I_idx pairs)
    std::vector<int> m_sum_ind;   // [n_cg_terms] output function index
    std::vector<double> cg_coeff; // [n_cg_terms]
    int n_output_funcs = 0;
    int n_cg_terms = 0;
    int Lmax = 0;
    int n_rad_max = 0;
    double inv_avg_n_neigh = 1.0;
  } yi;

  // B reducer (equivariant, NOT element-dependent)
  ReduceWeights B_reduce_YI, B_reduce_B0;
  int B_n_out = 64;
  int B_n_funcs = 31;
  bool B_is_elem_dependent = false;

  // Layer 2 FC/Products
  FCWeights fc_B1, fc_BB1, fc_BB2;
  int L2_fc_n_out = 64;
  CGProduct prod_BB, prod_BBB, prod_BBBB;

  // I2 reducer (scalar, element-dependent)
  ReduceWeights I2_reduce_B, I2_reduce_BB, I2_reduce_BBB, I2_reduce_BBBB;
  int I2_n_out = 17;
  bool I2_is_elem_dependent = true;

  // RMSNorm: I_0_LN (full: scale has n_out elements, all scaled)
  std::vector<double> I_0_LN_scale;  // [I2_n_out] = [17]
  int I_0_LN_type = 0;  // 0 = full
  int I_0_LN_n_out = 17;

  // Energy MLP (input=16, tanh activation)
  std::vector<MLPLayer> energy_mlp_layers;
  int energy_mlp_activation = 1;  // 0=silu, 1=tanh

  // Per-element shifts
  std::vector<double> shift_values;

  // ConstantScaleShiftTarget output scale; 1.0 if absent from npz
  double output_scale = 1.0;

  // Element names
  std::vector<std::string> element_names;

  // ---- UQ / extrapolation-grade artifacts (optional; present iff has_uq) ----
  // Schema v6 (basis-RP, normalize+density): the GMM feature phi is built from the
  // invariant (l=0) energy-path B-basis as proj = (B/||B||) * R (L2-normalized) with
  // n_density log-norm density channels appended: phi = concat(proj, dens). R is
  // stored verbatim; gamma is the only UQ signal (no force-error model).
  bool has_uq = false;
  int uq_schema_version = 0;
  int uq_n_elements = 0;
  int uq_max_clusters = 0;            // Kmax
  int uq_feature_dim = 0;             // D = rp_dim + n_density (full GMM feature dim)
  int uq_rp_dim = 0;                  // random-projection output dim (R cols; e.g. 128, read from artifact, not enforced)
  int uq_n_density = 0;              // density channels appended after proj (D - rp_dim)
  int uq_normalize = 0;             // 1 -> L2-normalize basis before projection (v6)
  double uq_density_scale = 1.0;     // density_scale (=1.0 for v6)
  int uq_d_basis = 0;                 // concatenated invariant basis width (rows of R)
  std::vector<double> uq_centroids;          // [E*Kmax*D]
  std::vector<double> uq_inv_cov;            // [E*Kmax*D*D] (pre-inverted precision)
  std::vector<int>    uq_n_clusters;         // [E]
  std::vector<double> uq_interp_thresholds;  // [E*Kmax]
  std::vector<double> uq_rp_matrix;          // [D_basis*rp_dim] row-major projection R

  void load(const std::string &filepath);
};

// ============================================================================
// PairGRACE2LKokkos: KOKKOS pair style for GRACE-2L
// ============================================================================
static constexpr int GRACE2L_MAX_MLP_DIM = 210;      // max hidden/output dim for MLP layers
static constexpr int GRACE2L_MAX_MLP_HIDDEN = 128;    // max input/hidden dim for MLP radial buffers (NOT output)
static constexpr int GRACE2L_MAX_NRADMAX = 42;        // max radial basis functions
static constexpr int GRACE2L_MAX_MLP_LAYERS = 8;      // max MLP depth
static constexpr int GRACE2L_MAX_YI_OUT_FUNCS = 200;  // max YI output functions
static constexpr int GRACE2L_MAX_I1_FUNCS = 16;       // max I1 equivariant channels (I1_n_funcs)
static constexpr int GRACE2L_MAX_I_FUNCS = 16;        // max I equivariant channels (I_n_funcs); matches grad_buf[16]
static constexpr int GRACE2L_UQ_MAX_RP_DIM = 132;       // bounds the FULL UQ feature dim D = rp_dim + n_density (z/f/delta stack arrays); 132 = 128 rp_dim + up to 3 density channels + 1 slack; a larger model fails the D > MAX guard at load

template<class DeviceType, typename NNScalarT = double, typename GeomScalarT = double>
class PairGRACE2LKokkos : public Pair, public KokkosBase {
 public:
  // Kernel tags
  struct TagPackForwardComm{};
  struct TagUnpackForwardComm{};
  struct TagPackReverseComm{};
  struct TagUnpackReverseComm{};
  struct TagComputeNeigh{};
  struct TagComputeRadialBasis{};
  struct TagComputeMLPRadial_R{};   // Layer 1 radial MLP
  struct TagComputeMLPRadial_R1{};  // Layer 2 radial MLP
  struct TagComputeAi{};            // Layer 1 single particle basis
  struct TagComputeAi_B0{};         // Layer 2 single particle basis
  struct TagComputeYI{};            // Layer 2 equivariant indicator basis
  struct TagComputeMLPEnergy{};
  struct TagComputeDerivative_L1{}; // Layer 1 force computation
  struct TagComputeUQ_L1basis{};    // basis-RP: L1 part of z = b_inv*R, captured in Phase 1 (d_A.. current)
  struct TagComputeUQ{};            // basis-RP: add L2 part + GMM (gamma/sigma/gmm_cluster) in Phase 3 — forward only
#ifdef KOKKOS_ENABLE_CUDA
  // Opt-2L-15: batched cuBLAS radial-MLP support kernels (CUDA + NNScalar=float only).
  struct TagMLPAssemble{};          // pack row-major X[M x nradbase] + dX[M x nradbase]; zero padding rows
  struct TagMLPActDeriv{};          // fused value+deriv silu between GEMMs (2L MLP has no bias)
  struct TagMLPScatterR{};          // scatter LayoutRight [M x d3] -> LayoutLeft d_{R,DR}{,1}_nl(ii,jj,n,l)
  struct TagUQBmat{};    // gather dense Bmat[Db x subN] + block ||B||^2 (cuBLAS UQ; level 1=rho, 2=I2)
  struct TagUQFinish{};  // combine L1 carry + L2 proj -> normalize + density + GMM
#endif

  template<int NEIGHFLAG, int EVFLAG>
  struct TagComputeForce{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  using GeomScalar = GeomScalarT;
  using NNScalar = NNScalarT;

  PairGRACE2LKokkos(class LAMMPS *);
  ~PairGRACE2LKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  // MPI communication for I features (host-side fallback)
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

  // Kernel operators
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeNeigh, const typename Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeRadialBasis, const typename Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPRadial_R, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial_R>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPRadial_R1, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial_R1>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeAi, const typename Kokkos::TeamPolicy<DeviceType, TagComputeAi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeAi_B0, const typename Kokkos::TeamPolicy<DeviceType, TagComputeAi_B0>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeYI, const typename Kokkos::TeamPolicy<DeviceType, TagComputeYI>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPEnergy, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeDerivative_L1, const typename Kokkos::TeamPolicy<DeviceType, TagComputeDerivative_L1>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeUQ_L1basis, const typename Kokkos::TeamPolicy<DeviceType, TagComputeUQ_L1basis>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeUQ, const typename Kokkos::TeamPolicy<DeviceType, TagComputeUQ>::member_type& team) const;

#ifdef KOKKOS_ENABLE_CUDA
  KOKKOS_INLINE_FUNCTION
  void operator() (TagUQBmat, const typename Kokkos::TeamPolicy<DeviceType, TagUQBmat>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagUQFinish, const typename Kokkos::TeamPolicy<DeviceType, TagUQFinish>::member_type& team) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPAssemble, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPActDeriv, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagMLPScatterR, const int) const;
#endif

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeForce<NEIGHFLAG,EVFLAG>, const int& ii) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeForce<NEIGHFLAG,EVFLAG>, const int& ii, EV_FLOAT&) const;

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

  // Layer-specific dimensions
  int L1_nradmax, L2_nradmax;
  int L1_fc_n_out, L2_fc_n_out;
  int I1_n_out, I1_n_funcs;
  int I_n_out, I_n_funcs;
  int rho_n_out;
  int B_n_out, B_n_funcs;
  int I2_n_out;
  int n_yi_cg, n_yi_out_funcs, n_yi_pairs;

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

  friend void pair_virial_fdotr_compute<PairGRACE2LKokkos>(PairGRACE2LKokkos*);

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
  typedef Kokkos::View<NNScalar****, DeviceType> t_nn_4d;

  // ---------- Per-chunk intermediate arrays (Layer 1) ----------

  // Neighbor list
  t_int_1d d_ncount;
  t_int_2d d_nearest;
  t_geom_2d d_rnorms;
  t_geom_3d3 d_rhats;
  t_int_1d d_mu_i;
  t_int_2d d_mu_j;

  // Radial basis (shared between R and R1 — sized for max(L1,L2) nradmax)
  t_geom_3d d_radial_basis;    // [chunk_size, maxneigh, nradbase]
  t_geom_3d d_dradial_basis;   // [chunk_size, maxneigh, nradbase]

  // Layer 1 MLP radial
  t_nn_4d d_R_nl;              // [chunk_size, maxneigh, L1_nradmax, lmax+1]
  // Opt-2L-16 precision specialization (allocated per instantiation, see grow):
  //  - float (fp32/mixed): d_DR_nl holds the output-layer radial derivative,
  //    precomputed once (cuBLAS DR SGEMM / fused hand kernel) so ComputeDerivative
  //    just reads it (register relief). d_dh2_R unused.
  //  - fp64: d_dh2_R holds the raw last-hidden fwd-mode derivative; ComputeDerivative
  //    forms DR inline via the W2^T*dh2 matvec (the original scheme). Precomputing the
  //    wider DR to global regresses fp64 (no SGEMM to amortize it), so fp64 keeps the
  //    old bit-identical path. d_DR_nl unused.
  t_nn_4d d_DR_nl;             // [chunk_size, maxneigh, L1_nradmax, lmax+1] (float only)
  t_nn_3d d_dh2_R;             // [chunk_size, maxneigh, mlp_rad_R_max_dim]   (fp64 only)

  // Layer 1 basis
  t_nn_3d d_A;                 // [chunk_size, L1_nradmax, nlm]
  t_nn_3d d_A1;                // [chunk_size, L1_fc_n_out, n_funcs_A1]
  t_nn_3d d_AA;                // [chunk_size, L1_fc_n_out, n_funcs_AA]
  t_nn_3d d_AA1;
  t_nn_3d d_AAA;
  t_nn_3d d_AA2;
  t_nn_3d d_AAAA;

  // Layer 1 equivariant features
  t_nn_3d d_I1;                // [chunk_size, I1_n_funcs, I1_n_out]
  t_nn_2d d_rho;               // [chunk_size, rho_n_out]

  // Layer 1 energy
  t_geom_1d d_e_atom;          // [chunk_size]
  t_nn_2d d_I_nl_LN;           // [chunk_size, rho_n_out] (after RMSNorm)

  // Layer 1 per-bond force
  t_geom_3d3 d_f_ij;           // [chunk_size, maxneigh, 3]

  // ---------- Per-chunk intermediate arrays (Layer 2) ----------

  // Layer 2 MLP radial
  t_nn_4d d_R1_nl;             // [chunk_size, maxneigh, L2_nradmax, lmax+1]
  t_nn_4d d_DR1_nl;            // [chunk_size, maxneigh, L2_nradmax, lmax+1] (float only, see d_DR_nl)
  t_nn_3d d_dh2_R1;            // [chunk_size, maxneigh, mlp_rad_R1_max_dim]  (fp64 only, see d_dh2_R)

  // Per-bond precomputed real spherical harmonics — hoisted out of
  // ComputeAi/ComputeAi_B0/ComputeYI/ReverseYI_grad_I (each previously
  // recomputed plm+Y_vals per (ii,jj,n), redundant L1_nradmax/L2_nradmax×).
  // Declared LayoutRight so the inner-most (lm) dim is contiguous — readers
  // sweep idx=0..n_yi for fixed (ii,jj) and need a cache-friendly stride.
  Kokkos::View<GeomScalar***, Kokkos::LayoutRight, DeviceType> d_Y_bond; // [chunk_size, maxneigh, (lmax+1)^2]

  // Layer 2 basis
  t_nn_3d d_B0;                // [chunk_size, L2_nradmax, nlm]
  t_nn_3d d_YI_pair;           // [chunk_size, L2_nradmax, n_unique(Y_lm,I_idx)]
  t_nn_3d d_YI;                // [chunk_size, L2_nradmax, n_yi_out_funcs]
  t_nn_3d d_B;                 // [chunk_size, B_n_out, B_n_funcs]
  t_nn_3d d_B1;                // [chunk_size, L2_fc_n_out, B_n_funcs]
  t_nn_3d d_BB;                // [chunk_size, L2_fc_n_out, n_funcs_BB]
  t_nn_3d d_BB1;
  t_nn_3d d_BBB;
  t_nn_3d d_BB2;
  t_nn_3d d_BBBB;

  // Layer 2 scalar features
  t_nn_2d d_I2;                // [chunk_size, I2_n_out]
  t_nn_2d d_I_0_LN;            // [chunk_size, I2_n_out] (after RMSNorm)

  // Layer 2 per-bond force
  t_geom_3d3 d_f_ij_L2;       // [chunk_size, maxneigh, 3]

  // ---------- Global arrays (persist across chunks, size = nall) ----------
  t_nn_3d d_I_global;           // [nall, I_n_funcs, I_n_out] — I features for comm
  t_nn_2d d_I_nl_LN_global;    // [nall, rho_n_out] — I_nl_LN for energy
  t_nn_3d d_grad_I_global;     // [nall, I_n_funcs, I_n_out] — grad_I for backward
  t_nn_2d d_I_nl_LN_adj_global; // [nall, rho_n_out] — energy MLP adjoint for L1 backward

  // Host mirrors for MPI communication
  typename t_nn_3d::host_mirror_type h_I_global;
  typename t_nn_3d::host_mirror_type h_grad_I_global;

  // ---------- Adjoint views (Layer 1 backward) ----------
  t_nn_2d d_rho_adj;
  t_nn_2d d_I_nl_LN_adj;       // [chunk_size, rho_n_out]
  t_nn_3d d_I1_adj;             // [chunk_size, I1_n_funcs, I1_n_out]
  t_nn_3d d_AAAA_adj;
  t_nn_3d d_AA2_adj;
  t_nn_3d d_AAA_adj;
  t_nn_3d d_AA1_adj;
  t_nn_3d d_AA_adj;
  t_nn_3d d_A1_adj;
  t_nn_3d d_A_adj;

  // ---------- Adjoint views (Layer 2 backward) ----------
  t_nn_2d d_I_0_LN_adj;        // [chunk_size, I2_n_out]
  t_nn_2d d_I2_adj;             // [chunk_size, I2_n_out]
  t_nn_3d d_BBBB_adj;
  t_nn_3d d_BB2_adj;
  t_nn_3d d_BBB_adj;
  t_nn_3d d_BB1_adj;
  t_nn_3d d_BB_adj;
  t_nn_3d d_B1_adj;
  t_nn_3d d_B_adj;              // [chunk_size, B_n_out, B_n_funcs]
  t_nn_3d d_YI_adj;             // [chunk_size, L2_nradmax, n_yi_out_funcs]
  t_nn_3d d_B0_adj;             // [chunk_size, L2_nradmax, nlm]

  // ---------- Constant device arrays (weights) ----------

  // Chemical embedding
  t_nn_2d d_chem_embed;        // [n_elements, embedding_size]

  // === Layer 1 weights ===

  // A indicator transform
  t_nn_2d d_A_lin_W;           // [embedding_size, L1_nradmax]
  t_nn_2d d_z_tr_A;            // [n_elements, L1_nradmax] — precomputed z_tr per element
  NNScalar A_lin_norm;
  GeomScalar L1_inv_avg_n_neigh;

  // MLP radial R weights (generalized N-layer)
  t_nn_3d d_mlp_rad_R_W;       // [n_layers, max_in, max_out]
  t_nn_1d d_mlp_rad_R_norms;
  t_int_1d d_mlp_rad_R_dims;   // [n_layers+1]
  int mlp_rad_R_n_layers;
  int mlp_rad_R_max_dim;

  // FC weights (Layer 1)
  t_nn_3d d_fc_A1_wl, d_fc_A1_wr;
  t_int_1d d_fc_A1_wtl, d_fc_A1_wtr;
  t_int_1d d_fc_A1_ct, d_fc_A1_cf;
  t_nn_1d d_fc_A1_nof;
  NNScalar fc_A1_nl, fc_A1_nr;

  t_nn_3d d_fc_AA1_wl, d_fc_AA1_wr;
  t_int_1d d_fc_AA1_wtl, d_fc_AA1_wtr;
  t_int_1d d_fc_AA1_ct, d_fc_AA1_cf;
  t_nn_1d d_fc_AA1_nof;
  NNScalar fc_AA1_nl, fc_AA1_nr;

  t_nn_3d d_fc_AA2_wl, d_fc_AA2_wr;
  t_int_1d d_fc_AA2_wtl, d_fc_AA2_wtr;
  t_int_1d d_fc_AA2_ct, d_fc_AA2_cf;
  t_nn_1d d_fc_AA2_nof;
  NNScalar fc_AA2_nl, fc_AA2_nr;

  // CG product metadata (Layer 1)
  t_int_1d d_prod_AA_li, d_prod_AA_ri, d_prod_AA_si;
  t_nn_1d d_prod_AA_cg;
  t_int_1d d_prod_AAA_li, d_prod_AAA_ri, d_prod_AAA_si;
  t_nn_1d d_prod_AAA_cg;
  t_int_1d d_prod_AAAA_li, d_prod_AAAA_ri, d_prod_AAAA_si;
  t_nn_1d d_prod_AAAA_cg;
  int n_funcs_A1, n_funcs_AA, n_funcs_AA1, n_funcs_AAA, n_funcs_AA2, n_funcs_AAAA;
  int n_cg_AA, n_cg_AAA, n_cg_AAAA;

  // I1 ReduceN (equivariant, element-dependent): W[n_elem, I1_n_out, n_in, w_shape]
  t_nn_4d d_I1_reduce_A_W, d_I1_reduce_AA_W, d_I1_reduce_AAA_W, d_I1_reduce_AAAA_W;
  NNScalar I1_reduce_A_norm, I1_reduce_AA_norm, I1_reduce_AAA_norm, I1_reduce_AAAA_norm;
  t_int_1d d_I1_ci_A, d_I1_ci_AA, d_I1_ci_AAA, d_I1_ci_AAAA;     // collect_ind
  t_int_1d d_I1_wlt_A, d_I1_wlt_AA, d_I1_wlt_AAA, d_I1_wlt_AAAA; // w_l_tile
  t_int_1d d_I1_tsi_A, d_I1_tsi_AA, d_I1_tsi_AAA, d_I1_tsi_AAAA; // total_sum_ind
  int I1_nc_A, I1_nc_AA, I1_nc_AAA, I1_nc_AAAA; // n_connections

  // I ReduceN (equivariant, NOT element-dependent): W[I_n_out, n_in, w_shape]
  t_nn_3d d_I_reduce_I1_W;     // [I_n_out, I1_n_out, w_shape]
  NNScalar I_reduce_I1_norm;
  t_int_1d d_I_ci_I1, d_I_wlt_I1, d_I_tsi_I1;
  int I_nc_I1;

  // rho ReduceN (scalar, element-dependent): W[n_elem, rho_n_out, n_in, w_shape]
  t_nn_4d d_rho_reduce_A_W, d_rho_reduce_AA_W, d_rho_reduce_AAA_W, d_rho_reduce_AAAA_W;
  NNScalar rho_reduce_A_norm, rho_reduce_AA_norm, rho_reduce_AAA_norm, rho_reduce_AAAA_norm;
  t_int_1d d_rho_ci_A, d_rho_ci_AA, d_rho_ci_AAA, d_rho_ci_AAAA;
  t_int_1d d_rho_wlt_A, d_rho_wlt_AA, d_rho_wlt_AAA, d_rho_wlt_AAAA;

  // I_nl_LN RMSNorm scale
  t_nn_1d d_I_nl_LN_scale;     // [rho_n_out - 1]

  // === Layer 2 weights ===

  // B0 indicator transform
  t_nn_2d d_B0_lin_W;          // [embedding_size, L2_nradmax]
  t_nn_2d d_z_tr_B0;           // [n_elements, L2_nradmax] — precomputed z_tr per element
  NNScalar B0_lin_norm;
  GeomScalar L2_inv_avg_n_neigh;

  // MLP radial R1 weights
  t_nn_3d d_mlp_rad_R1_W;
  t_nn_1d d_mlp_rad_R1_norms;
  t_int_1d d_mlp_rad_R1_dims;
  int mlp_rad_R1_n_layers;
  int mlp_rad_R1_max_dim;

#ifdef KOKKOS_ENABLE_CUDA
  // ---- Opt-2L-15: cuBLAS true-fp32 radial-MLP scratch + launch state ----
  // Active only for the CUDA device instantiation with NNScalar=float (fp32 &
  // Mixed styles). The fp64 device style and all host/non-CUDA backends keep the
  // hand kernels (handle stays null -> compute_mlp_radial() takes the fallback).
  cublasHandle_t cublas_handle = nullptr;
  typedef Kokkos::View<NNScalar**, Kokkos::LayoutRight, DeviceType> t_nn_2d_r; // row-major [M x dim]
  // Per-(which,layer) weights repacked TIGHTLY row-major [n_in x n_out] at load
  // time (the native padded d_mlp_rad_R{,1}_W is LayoutLeft -> a fixed layer's
  // (k,j) block is not contiguous). which: 0 = R (L1), 1 = R1 (L2). lda = n_out.
  t_nn_2d_r d_mlp_Wg[2][GRACE2L_MAX_MLP_LAYERS];
  // Batched scratch, M = chunk_size * maxneigh bonds. X/dX are the (shared)
  // Chebyshev radial-basis inputs; h/dh the two hidden activations (+fwd-mode
  // derivative); R the tightly-packed LayoutRight output GEMM (scattered into
  // the LayoutLeft d_{R,DR}{,1}_nl). d_mlp_R is reused for both the value (h1@W2)
  // and the output-layer derivative (dh1@W2) GEMMs.
  t_nn_2d_r d_mlp_X, d_mlp_dX, d_mlp_h0, d_mlp_h1, d_mlp_dh0, d_mlp_dh1, d_mlp_R;
  int mlp_M = 0;                 // current bonds = chunk_size * maxneigh
  int mlp_which = 0;            // 0 = R (L1), 1 = R1 (L2) — selects weights/output view
  int mlp_scatter_deriv = 0;    // TagMLPScatterR target: 0 = d_{R,R1}_nl, 1 = d_{DR,DR1}_nl
  int mlp_act_layer = 0;        // layer index for the activation kernel
  int mlp_act_n = 0;            // n_out width for the activation kernel
  NNScalar mlp_act_norm = 1;    // per-layer norm for the activation kernel

  // ---- Opt-2L-17: cuBLAS true-fp32 block-sparse FC (forward) ----
  // Each FC output is out(ii,k,func) [=|+=] norm*nof(func) * Σ_n W(k,n,tile(func))
  // * in(ii,n,src(func)), i.e. a contraction over n with a per-function TILE-
  // selected weight (5-27 tiles). Reformulated as: gather in-columns grouped by
  // tile -> one true-fp32 SGEMM per tile (W_t[n_out x n_in] @ G[n_in x cs*|funcs_t|])
  // -> scatter. The scatter keeps the hand kernel's (ii,k)-parallel, function-order
  // accumulation (no atomics) so only the Σ_n order changes vs the hand kernel.
  // Active only for CUDA device + NNScalar=float; fp64/host keep the FC hand kernels.
  // FCCublasOp is public because it appears in the signature of the public method
  // fc_forward_cublas (nvcc's extended-lambda registration needs the type visible);
  // the op instances themselves stay protected.
 public:
  struct FCCublasOp {
    int n_out = 0, n_in_l = 0, n_in_r = 0, N_l = 0, N_r = 0, ntiles = 0;
    NNScalar nl = 1, nr = 1;
    t_nn_1d Wl_packed, Wr_packed;   // [ntiles * n_out * n_in] col-major per tile
    t_nn_1d nof;                    // [N_l] per-output-function norm factor
    t_int_1d L_src, L_pos;          // L_src[p]=func(lm) tile-sorted; L_pos[lm]=p
    t_int_1d R_src, R_pos, R_tgt;   // R_src[p]=cf tile-sorted; R_pos[idx]=p; R_tgt[idx]=ct
    std::vector<int> L_pstart, L_count, R_pstart, R_count;  // host, [ntiles]
    bool active = false;
  };
  // Opt-2L-19: element-independent ReduceN (B <- YI, B0). Each connection c maps
  // an input func fi(c) to an output func fo(c) with a tile: B(ii,k,fo) += norm *
  // Σ_n W(k,n,tile) * X(ii,n,fi). Same gather->per-tile SGEMM->scatter shape as
  // the FC collect-half (accumulate, no left/assign half, single scalar norm, no
  // nof). Public for the same reason as FCCublasOp (appears in the signature of
  // the public reduce_forward_cublas). Float + CUDA only.
  struct ReduceCublasOp {
    int n_out = 0, n_in = 0, N = 0, ntiles = 0;
    NNScalar norm = 1;
    t_nn_1d W_packed;               // [ntiles * n_out * n_in] col-major per tile
    t_int_1d src;                   // [N] fi at tile-sorted col p (gather)
    t_int_1d pos, tgt;              // [N] pos[c]=sorted col; tgt[c]=fo (orig c order)
    std::vector<int> pstart, count; // host, [ntiles]
    bool active = false;
  };
 protected:
  FCCublasOp fc_op_A1, fc_op_AA1, fc_op_AA2, fc_op_B1, fc_op_BB1, fc_op_BB2;
  ReduceCublasOp reduce_op_B_YI, reduce_op_B_B0;
  t_nn_1d d_fc_G, d_fc_C;         // shared gather / GEMM-output scratch (float)
  int fc_cublas_max_nin = 0, fc_cublas_max_nout = 0, fc_cublas_max_N = 0;
#endif

  // YI CG coupling (sorted by output function, then y_idx = lr_inds[t*2])
  t_int_2d d_yi_lr_inds;       // [n_cg_terms, 2]
  t_int_1d d_yi_m_sum_ind;     // [n_cg_terms]
  t_nn_1d d_yi_cg_coeff;       // [n_cg_terms]
  t_int_1d d_yi_lm_offset;     // [nlm] — start offset per lm in sorted CG terms
  t_int_1d d_yi_lm_count;      // [nlm] — count of CG terms per lm
  t_int_1d d_yi_out_offset;    // [n_yi_out_funcs] — start offset per output function
  t_int_1d d_yi_out_count;     // [n_yi_out_funcs] — count of CG terms per output function
  t_int_1d d_yi_pair_yidx;     // [n_yi_pairs] — unique Y lm index
  t_int_1d d_yi_pair_iidx;     // [n_yi_pairs] — unique I channel
  t_int_1d d_yi_pair_l;        // [n_yi_pairs] — l value for the Y lm index
  t_int_1d d_yi_term_pair;     // [n_cg_terms] — pair index for each sorted CG term
  t_int_1d d_yi_pair_term_offset; // [n_yi_pairs] — start offset per pair in pair-term list
  t_int_1d d_yi_pair_term_count;  // [n_yi_pairs] — count per pair
  t_int_1d d_yi_pair_term_out;    // [n_cg_terms] — output function per pair-sorted term
  t_nn_1d d_yi_pair_term_cg;      // [n_cg_terms] — CG coefficient per pair-sorted term
  t_int_1d d_yi_i_pair_offset;    // [I_n_funcs] — start offset per I channel in pair list
  t_int_1d d_yi_i_pair_count;     // [I_n_funcs] — count of unique pairs per I channel
  t_int_1d d_yi_i_pair_index;     // [n_yi_pairs] — pair index grouped by I channel
  t_int_1d d_yi_i_yidx;        // [n_cg_terms] — I-sorted Y lm index for reverse dI
  t_int_1d d_yi_i_out;         // [n_cg_terms] — I-sorted output function for reverse dI
  t_nn_1d d_yi_i_cg;           // [n_cg_terms] — I-sorted CG coefficient for reverse dI
  t_int_1d d_yi_i_offset;      // [I_n_funcs] — start offset per I channel
  t_int_1d d_yi_i_count;       // [I_n_funcs] — count of CG terms per I channel
  GeomScalar yi_inv_avg_n_neigh;

  // B ReduceN (equivariant, NOT element-dependent)
  t_nn_3d d_B_reduce_YI_W;     // [B_n_out, n_in, w_shape]
  t_nn_3d d_B_reduce_B0_W;
  NNScalar B_reduce_YI_norm, B_reduce_B0_norm;
  t_int_1d d_B_ci_YI, d_B_ci_B0;
  t_int_1d d_B_wlt_YI, d_B_wlt_B0;
  t_int_1d d_B_tsi_YI, d_B_tsi_B0;
  int B_nc_YI, B_nc_B0;

  // FC weights (Layer 2)
  t_nn_3d d_fc_B1_wl, d_fc_B1_wr;
  t_int_1d d_fc_B1_wtl, d_fc_B1_wtr;
  t_int_1d d_fc_B1_ct, d_fc_B1_cf;
  t_nn_1d d_fc_B1_nof;
  NNScalar fc_B1_nl, fc_B1_nr;

  t_nn_3d d_fc_BB1_wl, d_fc_BB1_wr;
  t_int_1d d_fc_BB1_wtl, d_fc_BB1_wtr;
  t_int_1d d_fc_BB1_ct, d_fc_BB1_cf;
  t_nn_1d d_fc_BB1_nof;
  NNScalar fc_BB1_nl, fc_BB1_nr;

  t_nn_3d d_fc_BB2_wl, d_fc_BB2_wr;
  t_int_1d d_fc_BB2_wtl, d_fc_BB2_wtr;
  t_int_1d d_fc_BB2_ct, d_fc_BB2_cf;
  t_nn_1d d_fc_BB2_nof;
  NNScalar fc_BB2_nl, fc_BB2_nr;

  // CG product metadata (Layer 2)
  t_int_1d d_prod_BB_li, d_prod_BB_ri, d_prod_BB_si;
  t_nn_1d d_prod_BB_cg;
  t_int_1d d_prod_BBB_li, d_prod_BBB_ri, d_prod_BBB_si;
  t_nn_1d d_prod_BBB_cg;
  t_int_1d d_prod_BBBB_li, d_prod_BBBB_ri, d_prod_BBBB_si;
  t_nn_1d d_prod_BBBB_cg;
  int n_funcs_B1, n_funcs_BB, n_funcs_BB1, n_funcs_BBB, n_funcs_BB2, n_funcs_BBBB;
  int n_cg_BB, n_cg_BBB, n_cg_BBBB;

  // I2 ReduceN (scalar, element-dependent)
  t_nn_4d d_I2_reduce_B_W, d_I2_reduce_BB_W, d_I2_reduce_BBB_W, d_I2_reduce_BBBB_W;
  NNScalar I2_reduce_B_norm, I2_reduce_BB_norm, I2_reduce_BBB_norm, I2_reduce_BBBB_norm;
  t_int_1d d_I2_ci_B, d_I2_ci_BB, d_I2_ci_BBB, d_I2_ci_BBBB;
  t_int_1d d_I2_wlt_B, d_I2_wlt_BB, d_I2_wlt_BBB, d_I2_wlt_BBBB;

  // I_0_LN RMSNorm scale
  t_nn_1d d_I_0_LN_scale;      // [I2_n_out]

  // === Energy MLP ===
  t_nn_3d d_energy_W;           // [n_layers, max_in, max_out]
  t_nn_1d d_energy_norms;
  t_int_1d d_energy_dims;
  int energy_n_layers;
  int energy_max_dim;
  int energy_activation;        // 0=silu, 1=tanh

  // Per-element shifts
  t_nn_1d d_shifts;

  // ConstantScaleShiftTarget output scale (broadcast scalar)
  NNScalar output_scale;

  // Per-element-pair cutoffs
  t_geom_2d d_bond_cutoff;

  // Spherical harmonics precomputed coefficients
  t_geom_1d d_idx_sph, alm, blm, cl, dl;
  int idx_sph_max;

  // l value for each lm index (for YI coupling): l_from_lm[lm] = l
  t_int_1d d_l_from_lm;        // [(lmax+1)^2]

  // ---------- Helper functions ----------
  void grow(int natom, int maxneigh);
  void grow_global(int nall);
  void copy_weights_to_device();
  void precompute_harmonics();
  void allocate();
  // Opt-2L-15: radial MLP forward (value R + output-layer radial derivative DR)
  // for `which` (0 = R/L1, 1 = R1/L2). Dispatches to the batched true-fp32
  // cuBLAS SGEMM chain (CUDA + NNScalar=float + 3-layer MLP) or the hand kernel
  // fallback (every other case). Bit-reproduces the hand kernel's d_R{,1}_nl and
  // d_DR{,1}_nl outputs. Replaces the inline TagComputeMLPRadial_R{,1} launches.
  void compute_mlp_radial(int which);
#ifdef KOKKOS_ENABLE_CUDA
  // Opt-2L-15 helper: C[M x n_out] = alpha * In[M x n_in] * Wg[n_in x n_out]
  // (true-fp32 SGEMM, CUBLAS_PEDANTIC_MATH). Row-major operands; ld* are physical
  // row strides (allocated widths).
  void mlp_sgemm(const NNScalar* Wg, const NNScalar* In, NNScalar* C,
      int n_out, int M, int n_in, int ldW, int ldIn, int ldC, NNScalar alpha);
  // Opt-2L-15 hard gate: verify the layer-0 SGEMM (transpose/leading-dim
  // convention + true-fp32 path) reproduces a host fp64 reference before use.
  void mlp_cublas_selftest();
  // Opt-2L-17: build the tile-sorted metadata + repacked per-tile weights for one
  // block-sparse FC (host-side, from the model FCWeights). Also tracks scratch maxes.
  void build_fc_cublas_op(FCCublasOp &op, const GRACE2LModel::FCWeights &fc);
  // Opt-2L-19: same, for one element-independent ReduceN source (from ReduceWeights).
  void build_reduce_cublas_op(ReduceCublasOp &op, const GRACE2LModel::ReduceWeights &rw,
      int n_out);
#endif

 public:
  // Public because it contains an extended __host__ __device__ lambda; nvcc
  // forbids those inside private/protected members.
  // Fills d_Y_bond from d_rhats for the current chunk (uses member chunk_size,
  // maxneigh). Phase 1 computes d_Y_bond once per chunk; Phase 2 / Phase 5
  // recompute branches must call this after L2_ComputeNeigh refreshes d_rhats,
  // otherwise d_Y_bond holds stale values from Phase 1's last chunk.
  void compute_Y_bond_chunk();
#ifdef KOKKOS_ENABLE_CUDA
  // Opt-2L-17: run one forward FC (left assign + right accumulate) via gather ->
  // per-tile true-fp32 SGEMM -> scatter, into out (float only). Public for the
  // same reason as compute_Y_bond_chunk (contains extended device lambdas).
  void fc_forward_cublas(FCCublasOp &op, const t_nn_3d &in_left,
      const t_nn_3d &in_right, const t_nn_3d &out);
  // Opt-2L-19: run one element-independent ReduceN source (X -> Bout) via gather ->
  // per-tile true-fp32 SGEMM -> scatter-accumulate (float only). zero_first zeroes
  // all nBf output funcs before the first source. Public for the same reason as
  // fc_forward_cublas (contains extended device lambdas).
  void reduce_forward_cublas(ReduceCublasOp &op, const t_nn_3d &X,
      const t_nn_3d &Bout, bool zero_first, int nBf);
  // cuBLAS DGEMM UQ projection (two phases, sharing the d_uq_z carry):
  //  _L1 (Phase 1, d_A.. live): gather the rho block -> DGEMM -> RAW L1 proj into
  //      d_uq_z, and ||B^(L1)||^2 into d_uq_n2L1.
  //  _L2 (Phase 3, d_B.. live): gather the I2 block -> DGEMM -> add to d_uq_z ->
  //      normalize by the combined ||B|| -> density -> GMM.
  void compute_uq_cublas_L1();
  void compute_uq_cublas_L2();
#endif
 protected:

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template <typename scratch_type>
  KOKKOS_INLINE_FUNCTION
  int scratch_size_helper(int values_per_team) const;

  // Activation functions
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

  // ---------- UQ / extrapolation grade (forward only) ----------
  // FP64 device storage regardless of NNScalar (numerical safety).
  typedef Kokkos::View<double*, DeviceType>    t_uq_1d;
  // LayoutRight so the last (contiguous) index — reduced/split across team lanes in
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
  t_uq_2d d_uq_z;                  // [nlocal, rp_dim] per-atom RAW L1 projection (set in Phase 1)
  t_uq_1d d_uq_n2L1;               // [nlocal] ||B^(L1)||^2 from Phase 1 (block0 density carry)
  t_uq_1d d_gamma, d_sigma, d_gmm_cluster;

#ifdef KOKKOS_ENABLE_CUDA
  // ---- cuBLAS DGEMM UQ projection scratch (CUDA; built lazily on 1st UQ call) ----
  // Two gather maps replay the walk of each block (rho / I2) so Bmat rows line up
  // with the corresponding R-row range (I2 = rows [0,d_basis_L2); rho = rows
  // [d_basis_L2,d_basis)). Bmat/proj/n2 are sub-block scratch reused across L1/L2.
  Kokkos::View<double*, DeviceType> d_uq_Bmat;    // [max(Db_L1,Db_L2) * subN] col-major
  Kokkos::View<double*, DeviceType> d_uq_proj;    // [rp_dim * subN] col-major
  Kokkos::View<double*, DeviceType> d_uq_n2;      // [subN] current-block ||B_L2||^2 carry to finish
  t_int_1d d_uq_mapL1_view, d_uq_mapL1_n, d_uq_mapL1_col;  // [d_basis_L1] rho-block gather map (built lazily Phase 1)
  t_int_1d d_uq_mapL2_view, d_uq_mapL2_n, d_uq_mapL2_col;  // [d_basis_L2] I2-block gather map (built lazily Phase 3)
  int uq_subN = 0;                 // Bmat sub-chunk width
  int uq_bmat_s0 = 0;              // current sub-block base (into the chunk)
  int uq_bmat_level = 0;           // 1 = rho/L1 gather, 2 = I2/L2 gather (selects views/map/Db)
#endif

  bool has_uq = false;
  int uq_Kmax = 0, uq_D = 0;          // uq_D = full GMM feature dim (rp_dim + n_density)
  int uq_rp_dim = 0;                  // projection output width (R cols; e.g. 128, read from artifact, not enforced)
  int uq_n_density = 0;               // density channels appended after proj (uq_D - rp_dim)
  int uq_normalize = 0;               // 1 -> L2-normalize basis before projection (v6)
  double uq_density_scale = 1.0;      // density_scale
  int uq_d_basis = 0, uq_d_basis_L1 = 0;  // total basis width; L1 width = R-row offset where L2 begins
  int flag_compute_gamma = 0, flag_compute_atomic_sigma = 0, flag_compute_gmm_cluster = 0;
  int nmax_uq = 0;
  double *gamma = nullptr, *atomic_sigma = nullptr, *gmm_cluster = nullptr;

  // Comm helper views for Kokkos-native communication
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d v_buf;
  int comm_first;  // for unpack: starting ghost index

  // Host-side model data
  GRACE2LModel *grace_model;
  std::vector<int> h_type2model;
};

using PairGRACE2LKokkos_FP64_Device  = PairGRACE2LKokkos<LMPDeviceType, double>;
using PairGRACE2LKokkos_FP64_Host    = PairGRACE2LKokkos<LMPHostType,   double>;
using PairGRACE2LKokkos_Mixed_Device = PairGRACE2LKokkos<LMPDeviceType, float>;
using PairGRACE2LKokkos_Mixed_Host   = PairGRACE2LKokkos<LMPHostType,   float>;
using PairGRACE2LKokkos_FP32_Device  = PairGRACE2LKokkos<LMPDeviceType, float, float>;
using PairGRACE2LKokkos_FP32_Host    = PairGRACE2LKokkos<LMPHostType,   float, float>;

}    // namespace LAMMPS_NS

#endif
#endif
