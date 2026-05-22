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
   GRACE-1L KOKKOS implementation
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/1l/cpu/kk,              PairGRACE1LCPUKokkos_FP64_Device);
PairStyle(grace/1l/cpu/kk/device,       PairGRACE1LCPUKokkos_FP64_Device);
PairStyle(grace/1l/cpu/kk/host,         PairGRACE1LCPUKokkos_FP64_Host);
PairStyle(grace/1l/cpu/kk/mixed,        PairGRACE1LCPUKokkos_Mixed_Device);
PairStyle(grace/1l/cpu/kk/mixed/device, PairGRACE1LCPUKokkos_Mixed_Device);
PairStyle(grace/1l/cpu/kk/mixed/host,   PairGRACE1LCPUKokkos_Mixed_Host);
PairStyle(grace/1l/cpu/kk/fp32,         PairGRACE1LCPUKokkos_FP32_Device);
PairStyle(grace/1l/cpu/kk/fp32/device,  PairGRACE1LCPUKokkos_FP32_Device);
PairStyle(grace/1l/cpu/kk/fp32/host,    PairGRACE1LCPUKokkos_FP32_Host);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRACE_1L_CPU_KOKKOS_H
#define LMP_PAIR_GRACE_1L_CPU_KOKKOS_H

#include "pair.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

// ============================================================================
// GRACE1LCPUModel: Host-side weight container loaded from .npz file
// ============================================================================
struct GRACE1LCPUModel {
  // Architecture metadata
  int n_elements = 0;
  int embedding_size = 0;
  int lmax = 0;
  int nradbase = 0;
  int nradmax = 0;
  int radial_basis_p = 0;
  double rcut = 0.0;          // max cutoff (or single cutoff if uniform)
  bool has_bond_specific_cutoff = false;
  std::vector<double> bond_cutoff_map;  // [n_elements * n_elements] per-pair cutoffs

  // MLP radial function
  int mlp_rad_n_layers = 0;
  bool mlp_rad_has_chem_emb = false;

  // Chemical embedding [n_elements, embedding_size]
  std::vector<double> chem_embedding;

  // A indicator transform [embedding_size, nradmax]
  std::vector<double> A_lin_transform_W;
  double A_lin_transform_norm = 1.0;
  double A_inv_avg_n_neigh = 1.0;

  // MLP radial weights: layers[i] has W[n_in, n_out] and norm
  struct MLPLayer {
    std::vector<double> W;
    double norm = 1.0;
    int n_in = 0, n_out = 0;
  };
  std::vector<MLPLayer> mlp_rad_layers;

  // FC weights (A1, AA1, AA2)
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

  // Product CG metadata (AA, AAA, AAAA)
  struct CGProduct {
    std::vector<int> left_ind, right_ind, m_sum_ind;
    std::vector<double> cg_coeff;
    int n_output_funcs = 0;
    int n_cg_terms = 0;
  };
  CGProduct prod_AA, prod_AAA, prod_AAAA;

  // FunctionReduceN (rho)
  struct ReduceWeights {
    std::vector<double> W;  // [n_elements, n_out, n_in, w_shape]
    double norm = 1.0;
    int n_in = 0, w_shape = 0;
    std::vector<int> collect_ind;
    std::vector<int> w_l_tile;
    std::vector<int> total_sum_ind;
  };
  int rho_n_out = 0;
  bool rho_is_elem_dependent = true;
  ReduceWeights reduce_A, reduce_AA, reduce_AAA, reduce_AAAA;

  // Energy MLP
  std::vector<MLPLayer> energy_mlp_layers;
  int energy_mlp_activation = 0;  // 0=silu, 1=tanh

  // Per-element shifts
  std::vector<double> shift_values;

  // ConstantScaleShiftTarget output scale; 1.0 if absent from npz
  double output_scale = 1.0;

  // Element names (for type mapping)
  std::vector<std::string> element_names;

  void load(const std::string &filepath);
};

// ============================================================================
// PairGRACE1LCPUKokkos: KOKKOS pair style for GRACE-1L
// ============================================================================
// Compile-time max dimensions for stack arrays in CUDA kernels
static constexpr int GRACE1L_CPU_MAX_MLP_DIM = 64;   // max hidden dim for MLP layers
static constexpr int GRACE1L_CPU_MAX_NRADMAX = 32;    // max radial basis functions
static constexpr int GRACE1L_CPU_MAX_MLP_LAYERS = 8;  // max MLP depth
static constexpr int GRACE1L_CPU_MAX_LMAX = 4;
static constexpr int GRACE1L_CPU_MAX_PLM =
    (GRACE1L_CPU_MAX_LMAX + 1) * (GRACE1L_CPU_MAX_LMAX + 2) / 2;

template<class DeviceType, typename NNScalarT = double, typename GeomScalarT = double>
class PairGRACE1LCPUKokkos : public Pair {
 public:
  // Kernel tags
  struct TagComputeNeigh{};
  struct TagComputeRadialBasis{};
  struct TagComputeMLPRadial{};
  struct TagComputeAi{};
  // FC, Product, ReduceN dispatched via inline lambdas in compute()
  struct TagComputeMLPEnergy{};
  struct TagComputeDerivative{};

  // Force computation tags (Stage 2)
  template<int NEIGHFLAG, int EVFLAG>
  struct TagComputeForce{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  // Precision typedefs
  using GeomScalar = GeomScalarT;
  using NNScalar = NNScalarT;

  PairGRACE1LCPUKokkos(class LAMMPS *);
  ~PairGRACE1LCPUKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;

  // Kernel operators (to be implemented)
  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeNeigh, const typename Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeRadialBasis, const typename Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPRadial, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeAi, const int ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeMLPEnergy, const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeDerivative, const typename Kokkos::TeamPolicy<DeviceType, TagComputeDerivative>::member_type& team) const;

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

  // Architecture params (copied from GRACE1LCPUModel for device access)
  int nelements, lmax, nradmax, nradbase, embedding_size;
  int radial_basis_p;
  GeomScalar rcut;

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

  friend void pair_virial_fdotr_compute<PairGRACE1LCPUKokkos>(PairGRACE1LCPUKokkos*);

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

  // ---------- Per-chunk intermediate arrays ----------

  // Neighbor list (short, filtered by cutoff)
  t_int_1d d_ncount;           // [chunk_size] neighbor count per atom
  t_int_2d d_nearest;          // [chunk_size, maxneigh] neighbor indices
  t_geom_2d d_rnorms;          // [chunk_size, maxneigh] distances
  t_geom_3d3 d_rhats;          // [chunk_size, maxneigh, 3] unit vectors
  t_int_1d d_mu_i;             // [chunk_size] element type of center atom
  t_int_2d d_mu_j;             // [chunk_size, maxneigh] element type of neighbor

  // Radial basis g_k(r) and derivatives
  t_geom_3d d_radial_basis;    // [chunk_size, maxneigh, nradbase]
  t_geom_3d d_dradial_basis;   // [chunk_size, maxneigh, nradbase] (for forces)

  // MLP radial output and derivative intermediate
  t_nn_4d d_R_nl;              // [chunk_size, maxneigh, nradmax, lmax+1]
  t_nn_3d d_dh2;               // [chunk_size, maxneigh, mlp_hidden] (MLP layer 1 derivative)

  // Single particle basis A
  t_nn_3d d_A;                 // [chunk_size, (lmax+1)^2, nradmax]

  // FC outputs
  t_nn_3d d_A1;                // [chunk_size, n_funcs_A1, n_out_fc]
  t_nn_3d d_AA;                // [chunk_size, n_funcs_AA, n_out_fc]
  t_nn_3d d_AA1;               // [chunk_size, n_funcs_AA1, n_out_fc]
  t_nn_3d d_AAA;               // [chunk_size, n_funcs_AAA, n_out_fc]
  t_nn_3d d_AA2;               // [chunk_size, n_funcs_AA2, n_out_fc]
  t_nn_3d d_AAAA;              // [chunk_size, n_funcs_AAAA, n_out_fc]

  // Density and energy
  t_nn_2d d_rho;               // [chunk_size, rho_n_out]
  t_geom_1d d_e_atom;          // [chunk_size]

  // Per-bond force (for Stage 2)
  t_geom_3d3 d_f_ij;           // [chunk_size, maxneigh, 3]

  // Adjoint views for backward pass (Stage 2: Forces)
  t_nn_2d d_rho_adj;            // [chunk_size, rho_n_out]
  t_nn_3d d_AAAA_adj;           // [chunk_size, fc_n_out, n_funcs_AAAA]
  t_nn_3d d_AA2_adj;            // [chunk_size, fc_n_out, n_funcs_AA]
  t_nn_3d d_AAA_adj;            // [chunk_size, fc_n_out, n_funcs_AAA]
  t_nn_3d d_AA1_adj;            // [chunk_size, fc_n_out, n_funcs_AA]
  t_nn_3d d_AA_adj;             // [chunk_size, fc_n_out, n_funcs_AA]
  t_nn_3d d_A1_adj;             // [chunk_size, fc_n_out, n_funcs_A1]
  t_nn_3d d_A_adj;              // [chunk_size, nradmax, nlm] ("weights")

  // ---------- Constant device arrays (weights) ----------

  // Chemical embedding
  t_nn_2d d_chem_embed;        // [n_elements, embedding_size]

  // A indicator transform
  t_nn_2d d_A_lin_W;           // [embedding_size, nradmax]
  NNScalar A_lin_norm;
  GeomScalar inv_avg_n_neigh;

  // Precomputed z_tr[mu,n] = A_lin_norm * inv_avg_n_neigh * Σ_e A_lin_W[e,n] * chem_embed[mu,e]
  t_nn_2d d_z_tr;              // [n_elements, nradmax]

  // MLP radial weights (generalized N-layer)
  t_nn_3d d_mlp_rad_W;          // [n_layers, max_in, max_out] (padded)
  t_nn_1d d_mlp_rad_norms;      // [n_layers]
  t_int_1d d_mlp_rad_dims;      // [n_layers+1] — dims[0]=input, dims[i+1]=output of layer i
  int mlp_rad_n_layers;
  int mlp_rad_max_dim;           // max(all layer dims) for scratch sizing

  // FC weights
  t_nn_3d d_fc_A1_wl, d_fc_A1_wr;
  t_nn_3d d_fc_AA1_wl, d_fc_AA1_wr;
  t_nn_3d d_fc_AA2_wl, d_fc_AA2_wr;
  t_int_1d d_fc_A1_wtl, d_fc_A1_wtr;
  t_int_1d d_fc_AA1_wtl, d_fc_AA1_wtr;
  t_int_1d d_fc_AA2_wtl, d_fc_AA2_wtr;
  t_int_1d d_fc_A1_ct, d_fc_A1_cf;
  t_int_1d d_fc_AA1_ct, d_fc_AA1_cf;
  t_int_1d d_fc_AA2_ct, d_fc_AA2_cf;
  t_nn_1d d_fc_A1_nof, d_fc_AA1_nof, d_fc_AA2_nof;
  NNScalar fc_A1_nl, fc_A1_nr, fc_AA1_nl, fc_AA1_nr, fc_AA2_nl, fc_AA2_nr;
  int n_funcs_A1, n_funcs_AA, n_funcs_AA1, n_funcs_AAA, n_funcs_AA2, n_funcs_AAAA;
  int fc_n_out;  // 64 for all FC layers

  // Forward FC gather index: for each target lm, list of (src, tile) from right part
  // CSR format: entries [rg_start[lm], rg_start[lm+1]) index into rg_src/rg_tile
  t_int_1d d_fc_A1_rg_start, d_fc_A1_rg_src, d_fc_A1_rg_tile;
  t_int_1d d_fc_AA1_rg_start, d_fc_AA1_rg_src, d_fc_AA1_rg_tile;
  t_int_1d d_fc_AA2_rg_start, d_fc_AA2_rg_src, d_fc_AA2_rg_tile;

  // Reverse FC gather index: for each source lm, list of (tgt, tile) for right part backward
  t_int_1d d_fc_A1_rev_start, d_fc_A1_rev_src, d_fc_A1_rev_tile;
  t_int_1d d_fc_AA1_rev_start, d_fc_AA1_rev_src, d_fc_AA1_rev_tile;
  t_int_1d d_fc_AA2_rev_start, d_fc_AA2_rev_src, d_fc_AA2_rev_tile;

  // CG product metadata
  t_int_1d d_prod_AA_li, d_prod_AA_ri, d_prod_AA_si;
  t_nn_1d d_prod_AA_cg;
  t_int_1d d_prod_AAA_li, d_prod_AAA_ri, d_prod_AAA_si;
  t_nn_1d d_prod_AAA_cg;
  t_int_1d d_prod_AAAA_li, d_prod_AAAA_ri, d_prod_AAAA_si;
  t_nn_1d d_prod_AAAA_cg;
  int n_cg_AA, n_cg_AAA, n_cg_AAAA;

  // CG product gather indices (forward): for each output func, list of contributing CG terms
  // CSR: gs[n_out+1] start offsets, gli/gri/gcg[n_terms] sorted by output
  t_int_1d d_prod_AA_gs, d_prod_AA_gli, d_prod_AA_gri;  t_nn_1d d_prod_AA_gcg;
  t_int_1d d_prod_AAA_gs, d_prod_AAA_gli, d_prod_AAA_gri;  t_nn_1d d_prod_AAA_gcg;
  t_int_1d d_prod_AAAA_gs, d_prod_AAAA_gli, d_prod_AAAA_gri;  t_nn_1d d_prod_AAAA_gcg;
  // Reverse self-product gather: for each adjoint func f, entries where li==f or ri==f
  // Stores (oa_index, fwd_index, cg_coeff) — fwd_index is the "other" index (ri if li==f, li if ri==f)
  t_int_1d d_rprod_AA_gs, d_rprod_AA_goa, d_rprod_AA_gfwd;  t_nn_1d d_rprod_AA_gcg;
  t_int_1d d_rprod_AAAA_gs, d_rprod_AAAA_goa, d_rprod_AAAA_gfwd;  t_nn_1d d_rprod_AAAA_gcg;
  // Reverse asymmetric-product gather for AAA (AA1⊗A1): separate left/right CSRs
  t_int_1d d_rprod_AAA_left_gs, d_rprod_AAA_left_goa, d_rprod_AAA_left_gfwd;  t_nn_1d d_rprod_AAA_left_gcg;
  t_int_1d d_rprod_AAA_right_gs, d_rprod_AAA_right_goa, d_rprod_AAA_right_gfwd;  t_nn_1d d_rprod_AAA_right_gcg;

  // FunctionReduceN weights
  t_nn_4d d_reduce_A_W;        // [n_elements, rho_n_out, n_in, w_shape]
  t_nn_4d d_reduce_AA_W;
  t_nn_4d d_reduce_AAA_W;
  t_nn_4d d_reduce_AAAA_W;
  NNScalar reduce_A_norm, reduce_AA_norm, reduce_AAA_norm, reduce_AAAA_norm;
  t_int_1d d_reduce_A_ci, d_reduce_AA_ci, d_reduce_AAA_ci, d_reduce_AAAA_ci;
  int rho_n_out;

  // Energy MLP (generalized N-layer)
  t_nn_3d d_energy_W;           // [n_layers, max_in, max_out] (padded)
  t_nn_1d d_energy_norms;       // [n_layers]
  t_int_1d d_energy_dims;       // [n_layers+1]
  int energy_n_layers;
  int energy_max_dim;
  int energy_activation;        // 0=silu, 1=tanh (device-side mirror)

  // Per-element shifts
  t_nn_1d d_shifts;            // [n_elements]

  // ConstantScaleShiftTarget output scale (broadcast scalar)
  NNScalar output_scale;

  // Per-element-pair cutoffs
  t_geom_2d d_bond_cutoff;     // [n_elements, n_elements]

  // Spherical harmonics precomputed coefficients
  t_geom_1d d_idx_sph, alm, blm, cl, dl;
  int idx_sph_max;

  // ---------- Helper functions ----------
  void grow(int natom, int maxneigh);
  void copy_weights_to_device();
  void precompute_harmonics();
  void allocate();

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template <typename scratch_type>
  KOKKOS_INLINE_FUNCTION
  int scratch_size_helper(int values_per_team) const;

  // silu activation: x * sigmoid(x)
  KOKKOS_INLINE_FUNCTION
  static NNScalar silu(NNScalar x) {
    return x / (NNScalar(1.0) + Kokkos::exp(-x));
  }

  KOKKOS_INLINE_FUNCTION
  static NNScalar tanh_act(NNScalar x) {
    return Kokkos::tanh(x);
  }

  // Energy-only mode (skips backward pass / forces)
  int flag_compute_energy_only = 0;
  bool debug_no_energy_only_calc = false;

  // Host-side model data
  GRACE1LCPUModel *grace_model;
  std::vector<int> h_type2model;  // LAMMPS type -> model element index (for init_one)
};

using PairGRACE1LCPUKokkos_FP64_Device  = PairGRACE1LCPUKokkos<LMPDeviceType, double>;
using PairGRACE1LCPUKokkos_FP64_Host    = PairGRACE1LCPUKokkos<LMPHostType,   double>;
using PairGRACE1LCPUKokkos_Mixed_Device = PairGRACE1LCPUKokkos<LMPDeviceType, float>;
using PairGRACE1LCPUKokkos_Mixed_Host   = PairGRACE1LCPUKokkos<LMPHostType,   float>;
using PairGRACE1LCPUKokkos_FP32_Device  = PairGRACE1LCPUKokkos<LMPDeviceType, float, float>;
using PairGRACE1LCPUKokkos_FP32_Host    = PairGRACE1LCPUKokkos<LMPHostType,   float, float>;

}    // namespace LAMMPS_NS

#endif
#endif
