// clang-format off
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

   Task 2.1: class scaffolding, style registration and compile-time caps.
   Task 2.2: full weight loader (GRACE3LModel::load + copy_weights_to_device).
   Task 2.3: forward geometry front-end — neighbor list build, Chebyshev
   radial basis, polynomial envelope and real spherical harmonics Y_lm,
   chunked exactly like the validated GRACE-2L KOKKOS style. compute() still
   carries no NN/energy/force physics beyond this geometry precompute; that
   lands in later Phase 2/3 tasks.
------------------------------------------------------------------------- */

#include "pair_grace_3l_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"
#include "neigh_list_kokkos.h"

#include "cnpy/cnpy.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace LAMMPS_NS;

// Spherical harmonics constants (used by compute_Y_bond_chunk + the Task 3.3
// geometry-force plm/dplm recurrence; verbatim from the validated GRACE-2L
// KOKKOS style, where Y00 == 1.0 so the plm normalization matches d_Y_bond).
static constexpr double sq2   = 1.4142135623730950488;
static constexpr double sq3   = 1.7320508075688772935;   // sqrt(3)
static constexpr double sq3o2 = 1.2247448713915890491;   // sqrt(3/2)

// ======================================================================
// FindMaxNumNeighs helper
// ======================================================================
template<class DeviceType>
struct FindMaxNumNeighs3L {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;
  FindMaxNumNeighs3L(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs3L() {k_list.copymode = 1;}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& maxneigh) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (maxneigh < num_neighs) maxneigh = num_neighs;
  }
};

// ======================================================================
// npz helpers (generic; shared shape with 1L/2L loaders)
// ======================================================================

static std::vector<double> npz_get_double(const cnpy::npz_t &npz, const std::string &key) {
  auto it = npz.find(key);
  if (it == npz.end()) return {};
  const auto &arr = it->second;
  if (arr.word_size == 4) {                       // float32 -> upcast to double
    const float *p = arr.data<float>();
    return std::vector<double>(p, p + arr.num_vals);
  }
  const double *p = arr.data<double>();            // float64
  return std::vector<double>(p, p + arr.num_vals);
}

static std::vector<int> npz_get_int(const cnpy::npz_t &npz, const std::string &key) {
  auto it = npz.find(key);
  if (it == npz.end()) return {};
  const auto &arr = it->second;
  if (arr.word_size == 4) {
    const int32_t *p = arr.data<int32_t>();
    return std::vector<int>(p, p + arr.num_vals);
  } else {
    const int64_t *p = arr.data<int64_t>();
    return std::vector<int>(p, p + arr.num_vals);
  }
}

static int npz_get_int_scalar(const cnpy::npz_t &npz, const std::string &key, int default_val = 0) {
  auto v = npz_get_int(npz, key);
  return v.empty() ? default_val : v[0];
}

static double npz_get_double_scalar(const cnpy::npz_t &npz, const std::string &key, double default_val = 0.0) {
  auto v = npz_get_double(npz, key);
  return v.empty() ? default_val : v[0];
}

// require-form: throw if a key the pair style cannot run without is absent from
// the .npz. Catches the silent-zero failure mode where an optional npz_get_*
// returns an empty vector that would become zero-filled device storage (and
// hence zero forces). Read keys by EXACT name only — never prefix/startswith
// matching: e.g. the reduce `rho1`'s `rho1_norm_map` shares a prefix with the
// RMS-norm `rho1_norm`'s `rho1_norm_scale`; they are DISTINCT keys.
static std::vector<double> npz_require_double(const cnpy::npz_t &npz, const std::string &key) {
  auto v = npz_get_double(npz, key);
  if (v.empty())
    throw std::runtime_error("GRACE-3L/KK: required weight '" + key +
                             "' missing from .npz file (empty load -> zero forces). "
                             "Re-export the model with a current kokkos_export_3l.");
  return v;
}

static std::vector<int> npz_require_int(const cnpy::npz_t &npz, const std::string &key) {
  auto v = npz_get_int(npz, key);
  if (v.empty())
    throw std::runtime_error("GRACE-3L/KK: required index array '" + key +
                             "' missing from .npz file. "
                             "Re-export the model with a current kokkos_export_3l.");
  return v;
}

// shape helper: dims of a stored array (empty vector if absent).
static std::vector<size_t> npz_shape(const cnpy::npz_t &npz, const std::string &key) {
  auto it = npz.find(key);
  if (it == npz.end()) return {};
  return it->second.shape;
}

// decode an S-dtype names index array (element_names, *_names) to strings.
static std::vector<std::string> npz_names(const cnpy::npz_t &npz, const std::string &key) {
  std::vector<std::string> out;
  auto it = npz.find(key);
  if (it == npz.end()) return out;
  const auto &arr = it->second;
  const int n = arr.shape.empty() ? 0 : (int) arr.shape[0];
  const int wsize = arr.word_size;
  const char *raw = arr.data<char>();
  out.resize(n);
  for (int i = 0; i < n; i++) {
    std::string s(raw + (size_t) i * wsize, wsize);
    while (!s.empty() && s.back() == '\0') s.pop_back();
    out[i] = s;
  }
  return out;
}

// ======================================================================
// GRACE3LModel::load — minimal loader (globals + cutoffs + element names).
// Per-layer weight tensors are loaded by later Phase 2/3 tasks; the empty
// compute() needs no weights, but coeff()/init_one() need this metadata.
// ======================================================================

void GRACE3LModel::load(const std::string &filepath) {
  cnpy::npz_t npz = cnpy::npz_load(filepath);

  // ---- Global metadata ----
  n_elements = npz_get_int_scalar(npz, "n_elements", 0);
  embedding_size = npz_get_int_scalar(npz, "embedding_size", 64);
  nradbase = npz_get_int_scalar(npz, "nradbase", 10);
  // No bare `lmax` key in the 3L npz — the angular lmax is SPBF _lmax (A1_lmax).
  lmax = npz_get_int_scalar(npz, "lmax", npz_get_int_scalar(npz, "A1_lmax", 4));
  // Global radial_basis_p (=16) is NOT the envelope order; geometry uses the
  // SPBF p (=A1_p=5). Keep the global one only for the record; radial_basis_p
  // (the envelope order actually used) is overwritten from A1_p below.
  radial_basis_p_global = npz_get_int_scalar(npz, "radial_basis_p", 0);
  radial_basis_p = radial_basis_p_global;   // provisional; set to SPBF p below
  rcut = npz_get_double_scalar(npz, "rcut", 6.0);
  has_bond_specific_cutoff = npz_get_int_scalar(npz, "has_bond_specific_cutoff", 0) != 0;
  bond_cutoff_map = npz_get_double(npz, "bond_cutoff_map");
  if (bond_cutoff_map.empty())
    bond_cutoff_map.resize((size_t) n_elements * n_elements, rcut);

  chem_embedding = npz_require_double(npz, "chem_embedding");
  shift_values = npz_require_double(npz, "shift_values");
  output_scale = npz_get_double_scalar(npz, "output_scale", 1.0);  // absent -> 1.0
  element_names = npz_names(npz, "element_names");

  // ---- Helper: generic MLP (radial MLP has bias on hidden layers; energy MLP
  //      has no bias). W{i} raw + scalar norm{i}; b{i} (hidden only) if present. ----
  auto load_mlp = [&](const std::string &prefix, MLP &mlp) {
    mlp.n_layers = npz_get_int_scalar(npz, prefix + "_n_layers", 0);
    mlp.W.resize(mlp.n_layers);
    mlp.b.resize(mlp.n_layers);
    mlp.norm.resize(mlp.n_layers, 1.0);
    mlp.n_in.resize(mlp.n_layers, 0);
    mlp.n_out.resize(mlp.n_layers, 0);
    mlp.has_bias = false;
    for (int i = 0; i < mlp.n_layers; i++) {
      const std::string wk = prefix + "_W" + std::to_string(i);
      mlp.W[i] = npz_require_double(npz, wk);
      mlp.norm[i] = npz_get_double_scalar(npz, prefix + "_norm" + std::to_string(i), 1.0);
      auto sh = npz_shape(npz, wk);
      if (sh.size() == 2) { mlp.n_in[i] = (int) sh[0]; mlp.n_out[i] = (int) sh[1]; }
      // bias b{i} is stored [1, n_out] on hidden layers only (output has none).
      auto b = npz_get_double(npz, prefix + "_b" + std::to_string(i));
      if (!b.empty()) { mlp.b[i] = std::move(b); mlp.has_bias = true; }
    }
  };

  // ---- SPBF (A1 scalar, A2/A3 equivariant) ----
  auto spbf_names = npz_names(npz, "spbf_names");
  spbf.clear();
  spbf.reserve(spbf_names.size());
  for (const auto &nm : spbf_names) {
    SPBF s;
    s.name = nm;
    s.equivariant = npz_get_int_scalar(npz, nm + "_equivariant", 0) != 0;
    s.n_rad_max = npz_get_int_scalar(npz, nm + "_n_rad_max", 0);
    s.n_rad_basis = npz_get_int_scalar(npz, nm + "_n_rad_basis", 0);
    s.lmax = npz_get_int_scalar(npz, nm + "_lmax", 4);
    s.Lmax = npz_get_int_scalar(npz, nm + "_Lmax", 4);
    s.p = npz_get_int_scalar(npz, nm + "_p", 5);
    s.rcut = npz_get_double_scalar(npz, nm + "_rcut", rcut);
    s.inv_avg_n_neigh = npz_get_double_scalar(npz, nm + "_inv_avg_n_neigh", 1.0);
    s.l_tile = npz_require_int(npz, nm + "_l_tile");
    load_mlp(nm + "_mlp", s.mlp);
    if (!s.equivariant) {
      s.lin_transform_W = npz_require_double(npz, nm + "_lin_transform_W");
      s.lin_transform_norm = npz_get_double_scalar(npz, nm + "_lin_transform_norm", 1.0);
    } else {
      s.chem_linear_W = npz_require_double(npz, nm + "_chem_linear_W");
      s.chem_linear_norm = npz_get_double_scalar(npz, nm + "_chem_linear_norm", 1.0);
      s.chem_l0_mask = npz_require_double(npz, nm + "_chem_l0_mask");
      s.cg_W = npz_require_double(npz, nm + "_cg_W");
      s.nfunc = npz_get_int_scalar(npz, nm + "_nfunc", 0);
      s.n_lm_ind = npz_get_int_scalar(npz, nm + "_n_lm_ind", 0);
    }
    spbf.push_back(std::move(s));
  }
  // Effective radial ceilings + envelope order from the SPBF layers.
  if (spbf.size() >= 1) { L1_nradmax = spbf[0].n_rad_max; radial_basis_p = spbf[0].p; }
  if (spbf.size() >= 2) L2_nradmax = spbf[1].n_rad_max;
  if (spbf.size() >= 3) L3_nradmax = spbf[2].n_rad_max;

  // ---- cp_l products ----
  auto prod_names = npz_names(npz, "prod_names");
  prod.clear();
  prod.reserve(prod_names.size());
  for (const auto &nm : prod_names) {
    CPL c;
    c.name = nm;
    c.rank = npz_get_int_scalar(npz, nm + "_rank", 0);
    c.nfunc = npz_get_int_scalar(npz, nm + "_nfunc", 0);
    c.n_left = npz_get_int_scalar(npz, nm + "_n_left", 0);
    c.n_right = npz_get_int_scalar(npz, nm + "_n_right", 0);
    c.n_groups_left = npz_get_int_scalar(npz, nm + "_n_groups_left", 0);
    c.n_groups_right = npz_get_int_scalar(npz, nm + "_n_groups_right", 0);
    c.norm_u = npz_get_double_scalar(npz, nm + "_norm_u", 1.0);
    c.norm_v = npz_get_double_scalar(npz, nm + "_norm_v", 1.0);
    c.U = npz_require_double(npz, nm + "_U");
    c.V = npz_require_double(npz, nm + "_V");
    c.cg = npz_require_double(npz, nm + "_cg");
    c.n_cg = (int) c.cg.size();
    c.group_left = npz_require_int(npz, nm + "_group_left");
    c.group_right = npz_require_int(npz, nm + "_group_right");
    c.left_ind = npz_require_int(npz, nm + "_left_ind");
    c.right_ind = npz_require_int(npz, nm + "_right_ind");
    c.m_sum_ind = npz_require_int(npz, nm + "_m_sum_ind");
    c.out_l = npz_require_int(npz, nm + "_out_l");
    c.out_parity = npz_require_int(npz, nm + "_out_parity");
    prod.push_back(std::move(c));
  }

  // ---- FunctionReduceN (product reduces + eq*/rho*) ----
  auto reduce_names = npz_names(npz, "reduce_names");
  reduce.clear();
  reduce.reserve(reduce_names.size());
  for (const auto &nm : reduce_names) {
    Reduce r;
    r.name = nm;
    r.n_out = npz_get_int_scalar(npz, nm + "_n_out", 0);
    r.n_funcs = npz_get_int_scalar(npz, nm + "_n_funcs", 0);
    r.only_invar = npz_get_int_scalar(npz, nm + "_only_invar", 0) != 0;
    r.elem_dep = npz_get_int_scalar(npz, nm + "_is_central_atom_type_dependent", 0) != 0;
    auto instr_names = npz_names(npz, nm + "_instruction_names");
    r.instr.reserve(instr_names.size());
    for (const auto &cn : instr_names) {
      ReduceInstr ri;
      ri.name = cn;
      ri.W = npz_require_double(npz, nm + "_reduce_" + cn + "_W");
      ri.norm = npz_get_double_scalar(npz, nm + "_reduce_" + cn + "_norm", 1.0);
      ri.n_in = npz_get_int_scalar(npz, nm + "_n_in_" + cn, 0);
      ri.w_shape = npz_get_int_scalar(npz, nm + "_w_shape_" + cn, 0);
      ri.collect_ind = npz_require_int(npz, nm + "_collect_ind_" + cn);
      ri.w_l_tile = npz_require_int(npz, nm + "_w_l_tile_" + cn);
      ri.total_sum_ind = npz_require_int(npz, nm + "_total_sum_ind_" + cn);  // [n_conn,1] flat
      ri.n_conn = (int) ri.collect_ind.size();
      r.instr.push_back(std::move(ri));
    }
    // EXACT key — `{name}_norm_map` (NOT the RMS-norm `{name}_norm` footgun).
    r.norm_map = npz_get_double(npz, nm + "_norm_map");
    r.has_norm_map = !r.norm_map.empty();
    reduce.push_back(std::move(r));
  }

  // ---- FCRight2Left (keys prefixed fc_<name>; index list stores bare name) ----
  auto fc_names = npz_names(npz, "fc_names");
  fc.clear();
  fc.reserve(fc_names.size());
  for (const auto &nm : fc_names) {
    FC f;
    f.name = nm;
    const std::string p = "fc_" + nm;
    f.n_out = npz_get_int_scalar(npz, p + "_n_out", 0);
    f.left_coefs = npz_get_int_scalar(npz, p + "_left_coefs", 0) != 0;
    f.w_right = npz_require_double(npz, p + "_w_right");
    f.norm_right = npz_get_double_scalar(npz, p + "_norm_right", 1.0);
    f.w_tile_left = npz_require_int(npz, p + "_w_tile_left");
    f.w_tile_right = npz_require_int(npz, p + "_w_tile_right");
    f.collect_to = npz_require_int(npz, p + "_collect_to");
    f.collect_from = npz_require_int(npz, p + "_collect_from");
    f.norm_out_factor = npz_require_double(npz, p + "_norm_out_factor");  // [n_funcs_left,1,1] flat
    f.n_funcs_left = (int) f.w_tile_left.size();
    f.n_funcs_right = (int) f.collect_from.size();
    auto wr = npz_shape(npz, p + "_w_right");
    if (wr.size() == 3) { f.n_in_right = (int) wr[1]; f.w_shape_right = (int) wr[2]; }
    if (f.left_coefs) {   // 3L FCs are right->left (left_coefs=False); guard anyway
      f.w_left = npz_get_double(npz, p + "_w_left");
      f.norm_left = npz_get_double_scalar(npz, p + "_norm_left", 1.0);
      auto wl = npz_shape(npz, p + "_w_left");
      if (wl.size() == 3) { f.n_in_left = (int) wl[1]; f.w_shape_left = (int) wl[2]; }
    }
    fc.push_back(std::move(f));
  }

  // ---- EquivariantRMSNorm (eq1_norm, eq2_norm) ----
  auto eqnorm_names = npz_names(npz, "eqnorm_names");
  eqnorm.clear();
  eqnorm.reserve(eqnorm_names.size());
  for (const auto &nm : eqnorm_names) {
    EqNorm e;
    e.name = nm;
    e.affine_weight = npz_require_double(npz, nm + "_affine_weight");
    e.degree_weights = npz_require_double(npz, nm + "_degree_weights");
    e.l0_mask = npz_require_double(npz, nm + "_l0_mask");
    e.expand_index = npz_require_int(npz, nm + "_expand_index");
    e.eps = npz_get_double_scalar(npz, nm + "_eps", 1e-8);
    e.center_l0 = npz_get_int_scalar(npz, nm + "_center_l0", 1);
    e.M = (int) e.degree_weights.size();
    auto aw = npz_shape(npz, nm + "_affine_weight");
    if (aw.size() == 2) { e.n_groups = (int) aw[0]; e.n_out = (int) aw[1]; }
    eqnorm.push_back(std::move(e));
  }

  // ---- InvariantLayerRMSNorm (rho1/2/3_norm). EXACT keys — the reduce
  //      `rho1`'s `rho1_norm_map` is a different array from `rho1_norm_scale`. ----
  auto rmsnorm_names = npz_names(npz, "rmsnorm_names");
  rmsnorm.clear();
  rmsnorm.reserve(rmsnorm_names.size());
  for (const auto &nm : rmsnorm_names) {
    RMSNorm r;
    r.name = nm;
    r.scale = npz_require_double(npz, nm + "_scale");
    r.type = npz_get_int_scalar(npz, nm + "_type", 0);
    r.n_out = npz_get_int_scalar(npz, nm + "_n_out", (int) r.scale.size());
    rmsnorm.push_back(std::move(r));
  }

  // ---- Readout energy MLP (no bias) + 3-origin metadata ----
  load_mlp("energy_mlp", energy_mlp);
  energy_mlp_activation = npz_get_int_scalar(npz, "energy_mlp_activation", 0);  // 0=silu
  out1_origins = npz_names(npz, "out1_origins");
  energy_mlp_origin_n_out = npz_get_int_scalar(npz, "energy_mlp_origin_n_out", 0);
  energy_mlp_n_origins = npz_get_int_scalar(npz, "energy_mlp_n_origins", (int) out1_origins.size());

  // ---- UQ / extrapolation-grade artifacts (optional; schema v6 basis-RP) ----
  // Exclusive to uqv6 (schema 6: L2-normalize + log-norm density channels); older
  // v3/v4 artifacts are rejected at the gate so there is no silent miscompute. The
  // 3L invariant basis is the concatenation of the three per-layer scalar reduces
  // rho1 | rho2 | rho3 (sorted by reduce name), so uq_n_density == 4 ([full, rho1,
  // rho2, rho3]) — one more than 2L's 3.
  uq_schema_version = npz_get_int_scalar(npz, "uq_schema_version", 0);
  if (uq_schema_version != 0) {
    const int uq_schema_supported = 6;
    if (uq_schema_version != uq_schema_supported)
      throw std::runtime_error("GRACE-3L/KK: UQ schema v" +
                               std::to_string(uq_schema_version) + " unsupported by this GRACE "
                               "pair style build (expects v" + std::to_string(uq_schema_supported) +
                               ", uqv6 basis-RP normalize+density). Re-export the UQ artifacts at "
                               "schema v6 with a current grace_utils, or update the binary.");
    has_uq = true;
    uq_n_elements = npz_get_int_scalar(npz, "uq_n_elements", n_elements);
    uq_max_clusters = npz_get_int_scalar(npz, "uq_max_clusters", 0);
    uq_feature_dim = npz_get_int_scalar(npz, "uq_feature_dim", 0);
    uq_rp_dim = npz_get_int_scalar(npz, "uq_rp_dim", 0);
    // v6 feature-build flags (confirmatory: must match the v6 profile below).
    uq_normalize = npz_get_int_scalar(npz, "uq_rp_normalize", 0);
    const int uq_add_density = npz_get_int_scalar(npz, "uq_rp_add_density_channel", 0);
    uq_density_scale = npz_get_double_scalar(npz, "uq_rp_density_scale", 1.0);
    const int uq_transform = npz_get_int_scalar(npz, "uq_feature_transform", 0);
    uq_centroids = npz_require_double(npz, "uq_centroids");
    uq_inv_cov = npz_require_double(npz, "uq_inv_cov");
    uq_n_clusters = npz_require_int(npz, "uq_n_clusters");
    uq_interp_thresholds = npz_require_double(npz, "uq_interp_thresholds");
    uq_rp_matrix = npz_require_double(npz, "uq_rp_matrix");   // [D_basis, rp_dim] row-major
    // Confirm the feature-build flags agree with schema v6; a mismatch means the
    // artifact is corrupt/mislabeled (this build implements normalize+density only,
    // identity transform — no asinh, which would break the linear stream-accumulate).
    if (uq_normalize != 1 || uq_add_density != 1 || uq_transform != 0)
      throw std::runtime_error("GRACE-3L/KK UQ: schema v6 expects normalize=1, "
          "add_density_channel=1, feature_transform=0 (identity); got normalize=" +
          std::to_string(uq_normalize) + ", add_density=" + std::to_string(uq_add_density) +
          ", transform=" + std::to_string(uq_transform) + ". Artifact mislabeled or unsupported variant.");
    // Validate the dense arrays against the declared (E, Kmax, D) and derive
    // D_basis from R; copy_*d index these flat with no bounds check, so a
    // stale/mismatched export must fail here rather than read out of bounds.
    if (uq_rp_dim <= 0 || uq_feature_dim < uq_rp_dim)
      throw std::runtime_error("GRACE-3L/KK UQ: uq_feature_dim (" + std::to_string(uq_feature_dim) +
          ") must be >= uq_rp_dim (" + std::to_string(uq_rp_dim) + ")");
    uq_n_density = uq_feature_dim - uq_rp_dim;   // = 1 + n_blocks (=4 for 3L: full + rho1/rho2/rho3)
    // The 3L UQ kernel writes exactly the [full, rho1, rho2, rho3] density channels
    // (RP+0..RP+3); a higher count would leave f[RP+4..) uninitialized yet read by
    // the GMM loops.
    if (uq_n_density > 4)
      throw std::runtime_error("GRACE-3L/KK UQ: uq_n_density (" + std::to_string(uq_n_density) +
          ") > 4; the 3L UQ kernel only writes [full, rho1, rho2, rho3] density channels. "
          "Artifact has an unexpected invariant-block count.");
    if (uq_rp_matrix.size() % (size_t)uq_rp_dim != 0)
      throw std::runtime_error("GRACE-3L/KK UQ: uq_rp_matrix size (" +
          std::to_string(uq_rp_matrix.size()) + ") not divisible by rp_dim " + std::to_string(uq_rp_dim));
    uq_d_basis = (int)(uq_rp_matrix.size() / (size_t)uq_rp_dim);
    // Guard against a transposed/mis-strided R: the flat-size divisibility check
    // above is blind to a [rp_dim, D_basis] layout (identical element count), which
    // copy_2d would then read with the wrong stride -> silently wrong gamma. Verify
    // the stored 2D column count is rp_dim.
    {
      auto rp_it = npz.find("uq_rp_matrix");
      if (rp_it != npz.end() && rp_it->second.shape.size() == 2 &&
          (int)rp_it->second.shape[1] != uq_rp_dim)
        throw std::runtime_error("GRACE-3L/KK UQ: uq_rp_matrix has shape [" +
            std::to_string(rp_it->second.shape[0]) + ", " + std::to_string(rp_it->second.shape[1]) +
            "]; expected column count == uq_rp_dim (" + std::to_string(uq_rp_dim) +
            "). R may be transposed or mislabeled.");
    }
    const size_t E = uq_n_elements, K = uq_max_clusters, D = uq_feature_dim;
    auto uq_chk = [](const char *name, size_t got, size_t want) {
      if (got != want)
        throw std::runtime_error(std::string("GRACE-3L/KK UQ: ") + name + " has " +
            std::to_string(got) + " elements, expected " + std::to_string(want));
    };
    uq_chk("uq_centroids",         uq_centroids.size(),         E * K * D);
    uq_chk("uq_inv_cov",           uq_inv_cov.size(),           E * K * D * D);
    uq_chk("uq_n_clusters",        uq_n_clusters.size(),        E);
    uq_chk("uq_interp_thresholds", uq_interp_thresholds.size(), E * K);
    // uqv6 artifacts must cover ALL model elements with real GMM clusters; a zero
    // here means an incomplete export, which we reject rather than silently emit a
    // sentinel gamma for that element at runtime.
    for (size_t e = 0; e < uq_n_clusters.size(); e++)
      if (uq_n_clusters[e] <= 0)
        throw std::runtime_error("GRACE-3L/KK UQ: element index " + std::to_string(e) +
            " has no GMM clusters; uqv6 artifacts must cover all elements. Re-export the UQ "
            "artifacts with a current grace_utils.");
    // Defensive finite-value check: a partially-populated / mis-inverted export can
    // leave NaN/Inf in the GMM arrays (e.g. a singular covariance), which would
    // silently poison the Mahalanobis sigma at runtime. Reject at load instead.
    auto uq_finite = [](const char *name, const std::vector<double> &v) {
      for (double x : v)
        if (!std::isfinite(x))
          throw std::runtime_error(std::string("GRACE-3L/KK UQ: ") + name +
              " contains a non-finite value (NaN/Inf); re-export the UQ artifacts.");
    };
    uq_finite("uq_centroids",         uq_centroids);
    uq_finite("uq_inv_cov",           uq_inv_cov);
    uq_finite("uq_interp_thresholds", uq_interp_thresholds);
    uq_finite("uq_rp_matrix",         uq_rp_matrix);
  }
}

// ======================================================================
// Constructor / Destructor
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::PairGRACE3LKokkos(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  kokkosable = 1;
  reverse_comm_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  host_flag = (execution_space == Host);

  grace_model = new GRACE3LModel();
  // The 3L backward has a single-chunk fast path that avoids recomputing the
  // L1/L2 forward intermediates.  Use the largest cap validated to fit each
  // GPU memory class while retaining the user-facing chunksize override.
  chunksize = 4096;
#ifdef KOKKOS_ENABLE_CUDA
  int cuda_device = 0;
  cudaDeviceProp cuda_prop{};
  if (cudaGetDevice(&cuda_device) == cudaSuccess &&
      cudaGetDeviceProperties(&cuda_prop, cuda_device) == cudaSuccess) {
    constexpr size_t GiB = size_t(1) << 30;
    if (cuda_prop.totalGlobalMem >= 70 * GiB)
      chunksize = 16384;
    else if (cuda_prop.totalGlobalMem >= 40 * GiB)
      chunksize = 8192;
  }
#endif
  no_virial_fdotr_compute = 1;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::~PairGRACE3LKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_eatom, eatom);
  memoryKK->destroy_kokkos(k_vatom, vatom);
  memoryKK->destroy_kokkos(k_cvatom, cvatom);
  memory->destroy(gamma);
  memory->destroy(atomic_sigma);
  memory->destroy(gmm_cluster);
  delete grace_model;
#ifdef KOKKOS_ENABLE_CUDA
  if (cublas_handle) { cublasDestroy(cublas_handle); cublas_handle = nullptr; }
#endif
}

// ======================================================================
// settings
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::settings(int narg, char **arg)
{
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "chunksize") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style grace/3l/kk chunksize", error);
      chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (chunksize <= 0)
        error->all(FLERR, "pair_style grace/3l/kk chunksize must be > 0 (got {})", chunksize);
      iarg += 2;
    } else if (strcmp(arg[iarg], "debug_no_energy_only_calc") == 0) {
      debug_no_energy_only_calc = true;
      iarg += 1;
    } else {
      error->all(FLERR, "Unknown pair_style grace/3l/kk keyword: {}", arg[iarg]);
    }
  }
}

// ======================================================================
// extract: expose UQ GMM triggers (global scalar flags) to `fix pair`
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void *PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::extract(const char *str, int &dim)
{
  dim = 0;   // every value exposed here is a global scalar (fix pair triggers, flags)
  if (strcmp(str, "debug_no_energy_only_calc") == 0) return (void *) &debug_no_energy_only_calc;
  if (has_uq) {
    if (strcmp(str, "gamma_flag") == 0) return (void *) &flag_compute_gamma;
    if (strcmp(str, "atomic_sigma_flag") == 0) return (void *) &flag_compute_atomic_sigma;
    if (strcmp(str, "gmm_cluster_flag") == 0) return (void *) &flag_compute_gmm_cluster;
  }
  return nullptr;
}

// ======================================================================
// extract_peratom: expose per-atom UQ arrays to `fix pair`
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void *PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::extract_peratom(const char *str, int &ncol)
{
  ncol = 0;
  if (!has_uq) return nullptr;
  if (strcmp(str, "gamma") == 0) return (void *) gamma;
  if (strcmp(str, "atomic_sigma") == 0) return (void *) atomic_sigma;
  if (strcmp(str, "gmm_cluster") == 0) return (void *) gmm_cluster;
  return nullptr;
}

// ======================================================================
// allocate
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  map = new int[n];

  MemKK::realloc_kokkos(d_map, "grace_3l:map", n);
  MemKK::realloc_kokkos(k_cutsq, "grace_3l:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();
}

// ======================================================================
// coeff
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients");

  map_element2type(narg - 3, arg + 3);
  std::string weights_path = arg[2];
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-3L/KK] Loading weights from {}\n", weights_path);

  grace_model->load(weights_path);

  // Publish the UQ availability flag now (at pair_coeff time) so `fix pair` can
  // extract the gamma/atomic_sigma/gmm_cluster triggers, which it queries at
  // fix-creation time — BEFORE init_style()/copy_weights_to_device() stage the
  // GMM arrays onto the device. The device arrays + dims are set in init_style().
  has_uq = grace_model->has_uq;

  if (comm->me == 0) {
    auto &m = *grace_model;
    // Loaded-architecture summary (validation-gate anchor line).
    int emlp_in  = m.energy_mlp.n_layers > 0 ? m.energy_mlp.n_in.front()  : 0;
    int emlp_hid = m.energy_mlp.n_layers > 1 ? m.energy_mlp.n_out.front() : 0;
    int emlp_out = m.energy_mlp.n_layers > 0 ? m.energy_mlp.n_out.back()  : 0;
    double a1_norm0 = (!m.spbf.empty() && !m.spbf[0].mlp.norm.empty()) ? m.spbf[0].mlp.norm[0] : 0.0;
    double a1_inv_avg = m.spbf.empty() ? 0.0 : m.spbf[0].inv_avg_n_neigh;
    utils::logmesg(lmp,
      "[GRACE-3L/KK] GRACE-3L loaded: n_elements={} lmax={} nradbase={} radial_p={} rcut={} | "
      "SPBF={} products={} reduces={} FCs={} eqnorms={} rmsnorms={} | "
      "energy_mlp[{}->{}->{}] act={} origins={} | "
      "A1_inv_avg_n_neigh={} A1_mlp_norm0={} output_scale={}\n",
      m.n_elements, m.lmax, m.nradbase, m.radial_basis_p, m.rcut,
      m.spbf.size(), m.prod.size(), m.reduce.size(), m.fc.size(), m.eqnorm.size(), m.rmsnorm.size(),
      emlp_in, emlp_hid, emlp_out, m.energy_mlp_activation, m.energy_mlp_n_origins,
      a1_inv_avg, a1_norm0, m.output_scale);
  }

  // Build element mapping
  auto h_map = Kokkos::create_mirror_view(d_map);
  int nuser = narg - 3;
  for (int i = 1; i <= atom->ntypes; i++) {
    if (map[i] < 0 || map[i] >= nuser) {
      h_map(i) = -1;
      continue;
    }
    std::string user_elem = arg[3 + map[i]];
    int model_idx = -1;
    for (int j = 0; j < (int)grace_model->element_names.size(); j++) {
      if (grace_model->element_names[j] == user_elem) {
        model_idx = j;
        break;
      }
    }
    if (model_idx < 0)
      error->all(FLERR, "Element '{}' not found in GRACE-3L model element list", user_elem);
    h_map(i) = model_idx;
    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE-3L/KK] LAMMPS type {} -> {} (model index {})\n", i, user_elem, model_idx);
  }
  Kokkos::deep_copy(d_map, h_map);

  h_type2model.resize(atom->ntypes + 1, -1);
  for (int i = 1; i <= atom->ntypes; i++)
    h_type2model[i] = h_map(i);
}

// ======================================================================
// init_style
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace/3l/kk requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace/3l/kk requires newton pair on");

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR, "Must use half neighbor list style with pair grace/3l/kk");

  auto &m = *grace_model;

  // Copy scalar architecture params
  nelements = m.n_elements;
  lmax = m.lmax;
  nradbase = m.nradbase;
  embedding_size = m.embedding_size;
  radial_basis_p = m.radial_basis_p;
  rcut = m.rcut;

  L1_nradmax = m.L1_nradmax;
  L2_nradmax = m.L2_nradmax;
  L3_nradmax = m.L3_nradmax;

  // ---- Validate loaded dims against compile-time caps / stack bounds ----
  {
    const int nlm = (lmax + 1) * (lmax + 1);
    const int plm_size = (lmax + 1) * (lmax + 2) / 2;
    if (nlm > 25)
      error->all(FLERR, "GRACE-3L/KK: (lmax+1)^2={} exceeds stack bound 25 (lmax={})", nlm, lmax);
    if (plm_size > 15)
      error->all(FLERR, "GRACE-3L/KK: plm size {} exceeds stack bound 15 (lmax={})", plm_size, lmax);
    if (lmax > GRACE3L_LMAX)
      error->all(FLERR, "GRACE-3L/KK: lmax={} exceeds GRACE3L_LMAX={}", lmax, GRACE3L_LMAX);
    if (nelements > GRACE3L_MAX_ATOM_TYPES)
      error->all(FLERR, "GRACE-3L/KK: n_elements={} exceeds GRACE3L_MAX_ATOM_TYPES={}", nelements, GRACE3L_MAX_ATOM_TYPES);

    if ((int) m.spbf.size() > GRACE3L_N_SPBF)
      error->all(FLERR, "GRACE-3L/KK: #SPBF={} exceeds GRACE3L_N_SPBF={}", (int) m.spbf.size(), GRACE3L_N_SPBF);
    if ((int) m.prod.size() > GRACE3L_MAX_PRODS)
      error->all(FLERR, "GRACE-3L/KK: #products={} exceeds GRACE3L_MAX_PRODS={}", (int) m.prod.size(), GRACE3L_MAX_PRODS);
    if ((int) m.reduce.size() > GRACE3L_MAX_REDUCES)
      error->all(FLERR, "GRACE-3L/KK: #reduces={} exceeds GRACE3L_MAX_REDUCES={}", (int) m.reduce.size(), GRACE3L_MAX_REDUCES);
    if ((int) m.fc.size() > GRACE3L_MAX_FCS)
      error->all(FLERR, "GRACE-3L/KK: #FCs={} exceeds GRACE3L_MAX_FCS={}", (int) m.fc.size(), GRACE3L_MAX_FCS);
    if ((int) m.eqnorm.size() > GRACE3L_MAX_EQNORMS)
      error->all(FLERR, "GRACE-3L/KK: #eqnorms={} exceeds GRACE3L_MAX_EQNORMS={}", (int) m.eqnorm.size(), GRACE3L_MAX_EQNORMS);
    // EquivariantRMSNorm lm-axis M (25 for eq1_norm, 50 for eq2_norm) shares the
    // downstream SPBF indicator lm-axis; both are bounded by the same cap. This
    // guards the per-thread l0mean/colsum stack arrays in TagEqNorm/TagEqNormBwd.
    for (const auto &e : m.eqnorm)
      if (e.M > GRACE3L_MAX_SPBF_NLMIND)
        error->all(FLERR, "GRACE-3L/KK: {}.M={} exceeds GRACE3L_MAX_SPBF_NLMIND={}", e.name, e.M, GRACE3L_MAX_SPBF_NLMIND);
    if ((int) m.rmsnorm.size() > GRACE3L_MAX_RMSNORMS)
      error->all(FLERR, "GRACE-3L/KK: #rmsnorms={} exceeds GRACE3L_MAX_RMSNORMS={}", (int) m.rmsnorm.size(), GRACE3L_MAX_RMSNORMS);

    for (const auto &s : m.spbf) {
      if (s.n_rad_max > GRACE3L_MAX_NRADMAX)
        error->all(FLERR, "GRACE-3L/KK: {}.n_rad_max={} exceeds GRACE3L_MAX_NRADMAX={}", s.name, s.n_rad_max, GRACE3L_MAX_NRADMAX);
      if (s.n_rad_basis > GRACE3L_MAX_NRADBASIS)
        error->all(FLERR, "GRACE-3L/KK: {}.n_rad_basis={} exceeds GRACE3L_MAX_NRADBASIS={}", s.name, s.n_rad_basis, GRACE3L_MAX_NRADBASIS);
      if (s.mlp.n_layers > GRACE3L_MAX_MLP_LAYERS)
        error->all(FLERR, "GRACE-3L/KK: {}.mlp_n_layers={} exceeds GRACE3L_MAX_MLP_LAYERS={}", s.name, s.mlp.n_layers, GRACE3L_MAX_MLP_LAYERS);
      for (int i = 0; i < s.mlp.n_layers; i++)
        if (s.mlp.n_out[i] > GRACE3L_MAX_MLP_DIM || s.mlp.n_in[i] > GRACE3L_MAX_MLP_DIM)
          error->all(FLERR, "GRACE-3L/KK: {}.mlp layer {} dim exceeds GRACE3L_MAX_MLP_DIM={}", s.name, i, GRACE3L_MAX_MLP_DIM);
      if (s.equivariant) {
        if (s.nfunc > GRACE3L_MAX_SPBF_NFUNC)
          error->all(FLERR, "GRACE-3L/KK: {}.nfunc={} exceeds GRACE3L_MAX_SPBF_NFUNC={}", s.name, s.nfunc, GRACE3L_MAX_SPBF_NFUNC);
        if (s.n_lm_ind > GRACE3L_MAX_SPBF_NLMIND)
          error->all(FLERR, "GRACE-3L/KK: {}.n_lm_ind={} exceeds GRACE3L_MAX_SPBF_NLMIND={}", s.name, s.n_lm_ind, GRACE3L_MAX_SPBF_NLMIND);
      }
    }
    for (const auto &c : m.prod) {
      if (c.rank > GRACE3L_MAX_PROD_RANK)
        error->all(FLERR, "GRACE-3L/KK: {}.rank={} exceeds GRACE3L_MAX_PROD_RANK={}", c.name, c.rank, GRACE3L_MAX_PROD_RANK);
      if (c.nfunc > GRACE3L_MAX_PROD_NFUNC)
        error->all(FLERR, "GRACE-3L/KK: {}.nfunc={} exceeds GRACE3L_MAX_PROD_NFUNC={}", c.name, c.nfunc, GRACE3L_MAX_PROD_NFUNC);
      if (c.n_left > GRACE3L_MAX_PROD_NIN || c.n_right > GRACE3L_MAX_PROD_NIN)
        error->all(FLERR, "GRACE-3L/KK: {} n_left/n_right exceeds GRACE3L_MAX_PROD_NIN={}", c.name, GRACE3L_MAX_PROD_NIN);
      if (c.n_groups_left > GRACE3L_MAX_PROD_NGROUPS || c.n_groups_right > GRACE3L_MAX_PROD_NGROUPS)
        error->all(FLERR, "GRACE-3L/KK: {} n_groups exceeds GRACE3L_MAX_PROD_NGROUPS={}", c.name, GRACE3L_MAX_PROD_NGROUPS);
      if (c.n_cg > GRACE3L_MAX_PROD_NCG)
        error->all(FLERR, "GRACE-3L/KK: {}.n_cg={} exceeds GRACE3L_MAX_PROD_NCG={}", c.name, c.n_cg, GRACE3L_MAX_PROD_NCG);
    }
    for (const auto &r : m.reduce) {
      if ((int) r.instr.size() > GRACE3L_MAX_REDUCE_INSTR)
        error->all(FLERR, "GRACE-3L/KK: {} has {} collectors > GRACE3L_MAX_REDUCE_INSTR={}", r.name, (int) r.instr.size(), GRACE3L_MAX_REDUCE_INSTR);
      if (r.n_out > GRACE3L_MAX_REDUCE_NOUT)
        error->all(FLERR, "GRACE-3L/KK: {}.n_out={} exceeds GRACE3L_MAX_REDUCE_NOUT={}", r.name, r.n_out, GRACE3L_MAX_REDUCE_NOUT);
      if (r.n_funcs > GRACE3L_MAX_REDUCE_NFUNCS)
        error->all(FLERR, "GRACE-3L/KK: {}.n_funcs={} exceeds GRACE3L_MAX_REDUCE_NFUNCS={}", r.name, r.n_funcs, GRACE3L_MAX_REDUCE_NFUNCS);
      for (const auto &ri : r.instr) {
        if (ri.n_in > GRACE3L_MAX_REDUCE_NIN)
          error->all(FLERR, "GRACE-3L/KK: {}[{}].n_in={} exceeds GRACE3L_MAX_REDUCE_NIN={}", r.name, ri.name, ri.n_in, GRACE3L_MAX_REDUCE_NIN);
        if (ri.w_shape > GRACE3L_MAX_REDUCE_WSHAPE)
          error->all(FLERR, "GRACE-3L/KK: {}[{}].w_shape={} exceeds GRACE3L_MAX_REDUCE_WSHAPE={}", r.name, ri.name, ri.w_shape, GRACE3L_MAX_REDUCE_WSHAPE);
      }
    }
    for (const auto &f : m.fc) {
      if (f.n_out > GRACE3L_MAX_FC_NOUT)
        error->all(FLERR, "GRACE-3L/KK: fc_{}.n_out={} exceeds GRACE3L_MAX_FC_NOUT={}", f.name, f.n_out, GRACE3L_MAX_FC_NOUT);
      if (f.w_shape_right > GRACE3L_MAX_FC_TILE)
        error->all(FLERR, "GRACE-3L/KK: fc_{} tile={} exceeds GRACE3L_MAX_FC_TILE={}", f.name, f.w_shape_right, GRACE3L_MAX_FC_TILE);
      // TagFCRight2Left (Task 2.6) accumulates into a GRACE3L_MAX_FC_NFUNCS-sized
      // stack array indexed by n_funcs_left/right; catch an overflowing model here
      // rather than silently corrupting the stack on device.
      if (f.n_funcs_left > GRACE3L_MAX_FC_NFUNCS)
        error->all(FLERR, "GRACE-3L/KK: fc_{}.n_funcs_left={} exceeds GRACE3L_MAX_FC_NFUNCS={}", f.name, f.n_funcs_left, GRACE3L_MAX_FC_NFUNCS);
      if (f.n_funcs_right > GRACE3L_MAX_FC_NFUNCS)
        error->all(FLERR, "GRACE-3L/KK: fc_{}.n_funcs_right={} exceeds GRACE3L_MAX_FC_NFUNCS={}", f.name, f.n_funcs_right, GRACE3L_MAX_FC_NFUNCS);
    }
    for (const auto &r : m.rmsnorm)
      if (r.n_out > GRACE3L_MAX_RHO_NOUT)
        error->all(FLERR, "GRACE-3L/KK: {}.n_out={} exceeds GRACE3L_MAX_RHO_NOUT={}", r.name, r.n_out, GRACE3L_MAX_RHO_NOUT);
    if (m.energy_mlp.n_layers > GRACE3L_MAX_MLP_LAYERS)
      error->all(FLERR, "GRACE-3L/KK: energy_mlp n_layers={} exceeds GRACE3L_MAX_MLP_LAYERS={}", m.energy_mlp.n_layers, GRACE3L_MAX_MLP_LAYERS);
    // TagComputeMLPEnergy/TagReadoutBwd stack arrays (ebuf_a/ebuf_b, adj_a/adj_b,
    // pre_act_buf) are sized GRACE3L_MAX_ENERGY_MLP_DIM; energy_max_dim itself is
    // only computed later in copy_weights_to_device, so validate every energy-MLP
    // layer width here (before any weights are staged to device).
    for (int i = 0; i < m.energy_mlp.n_layers; i++)
      if (m.energy_mlp.n_out[i] > GRACE3L_MAX_ENERGY_MLP_DIM || m.energy_mlp.n_in[i] > GRACE3L_MAX_ENERGY_MLP_DIM)
        error->all(FLERR, "GRACE-3L/KK: energy_mlp layer {} dim exceeds GRACE3L_MAX_ENERGY_MLP_DIM={}", i, GRACE3L_MAX_ENERGY_MLP_DIM);
  }

  // ---- Spherical harmonics precompute (GeomScalar; verbatim from GRACE-2L) ----
  MemKK::realloc_kokkos(d_idx_sph, "g3l:idx_sph", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(alm, "g3l:alm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(blm, "g3l:blm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(cl, "g3l:cl", lmax + 1);
  MemKK::realloc_kokkos(dl, "g3l:dl", lmax + 1);
  precompute_harmonics();

  // l_from_lm lookup: for each lm index, store the corresponding l
  {
    int nlm = (lmax + 1) * (lmax + 1);
    MemKK::realloc_kokkos(d_l_from_lm, "g3l:l_from_lm", nlm);
    auto h_lfm = Kokkos::create_mirror_view(d_l_from_lm);
    for (int l = 0; l <= lmax; l++)
      for (int mm = -l; mm <= l; mm++)
        h_lfm(l * (l + 1) + mm) = l;
    Kokkos::deep_copy(d_l_from_lm, h_lfm);
  }

  // ---- Copy all weights to device ----
  copy_weights_to_device();

  // ---- Task 2.6: locate A1_2_red / fc_A1_2a by name (index-parallel with
  // d_reduce/d_fc, filled in the same loop order by copy_weights_to_device
  // above). Match by NAME (not position) so a re-ordered npz still resolves. ----
  idx_reduce_A1_2_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A1_2_red") { idx_reduce_A1_2_red = i; break; }
  if (idx_reduce_A1_2_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_2_red' not found in loaded model");
  if (m.reduce[idx_reduce_A1_2_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A1_2_red].instr[0].name != "A1_2")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_2_red' instruction layout unexpected "
               "(expected single instruction named 'A1_2'; compute() wires d_A1_2 "
               "into inputs[0] positionally)");

  idx_fc_A1_2a = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A1_2a") { idx_fc_A1_2a = i; break; }
  if (idx_fc_A1_2a < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A1_2a' not found in loaded model");

  // ---- Task 2.7: locate the remaining L1 chain by NAME (index-parallel with
  // d_prod/d_reduce/d_fc/d_eqnorm/d_rmsnorm, filled in the same order as
  // grace_model->{prod,reduce,fc,eqnorm,rmsnorm} by copy_weights_to_device). ----
  idx_prod_A1_3 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A1_3") { idx_prod_A1_3 = i; break; }
  if (idx_prod_A1_3 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A1_3' not found in loaded model");

  idx_prod_A1_4 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A1_4") { idx_prod_A1_4 = i; break; }
  if (idx_prod_A1_4 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A1_4' not found in loaded model");

  idx_reduce_A1_3_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A1_3_red") { idx_reduce_A1_3_red = i; break; }
  if (idx_reduce_A1_3_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_3_red' not found in loaded model");
  if (m.reduce[idx_reduce_A1_3_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A1_3_red].instr[0].name != "A1_3")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_3_red' instruction layout unexpected "
               "(expected single instruction named 'A1_3')");

  idx_fc_A1_2b = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A1_2b") { idx_fc_A1_2b = i; break; }
  if (idx_fc_A1_2b < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A1_2b' not found in loaded model");

  idx_reduce_A1_4_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A1_4_red") { idx_reduce_A1_4_red = i; break; }
  if (idx_reduce_A1_4_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_4_red' not found in loaded model");
  if (m.reduce[idx_reduce_A1_4_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A1_4_red].instr[0].name != "A1_4")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A1_4_red' instruction layout unexpected "
               "(expected single instruction named 'A1_4')");

  // eq1 / rho1: multi-instruction (4 collectors), elem-dependent reduces over
  // the SAME 4 named L1 tensors. Resolve role[j] by matching instr[j].name —
  // do NOT assume the npz's instruction_names order matches this enumeration.
  auto resolve_l1_roles = [&](const GRACE3LModel::Reduce &r, int roles[GRACE3L_MAX_REDUCE_INSTR]) {
    static const char *l1_names[4] = {"A1", "A1_2_red", "A1_3_red", "A1_4_red"};
    if ((int) r.instr.size() != 4)
      error->all(FLERR, "GRACE-3L/KK: reduce '{}' expected 4 L1 collectors, found {}",
                 r.name, (int) r.instr.size());
    for (int j = 0; j < 4; j++) {
      int role = -1;
      for (int k = 0; k < 4; k++)
        if (r.instr[j].name == l1_names[k]) { role = k; break; }
      if (role < 0)
        error->all(FLERR, "GRACE-3L/KK: reduce '{}' instruction[{}] name '{}' is not one "
                   "of the expected L1 inputs {{A1,A1_2_red,A1_3_red,A1_4_red}}",
                   r.name, j, r.instr[j].name);
      roles[j] = role;
    }
  };

  idx_reduce_eq1 = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "eq1") { idx_reduce_eq1 = i; break; }
  if (idx_reduce_eq1 < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'eq1' not found in loaded model");
  resolve_l1_roles(m.reduce[idx_reduce_eq1], eq1_input_role);

  idx_reduce_rho1 = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "rho1") { idx_reduce_rho1 = i; break; }
  if (idx_reduce_rho1 < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'rho1' not found in loaded model");
  resolve_l1_roles(m.reduce[idx_reduce_rho1], rho1_input_role);

  idx_eqnorm_eq1_norm = -1;
  for (int i = 0; i < (int) m.eqnorm.size(); i++)
    if (m.eqnorm[i].name == "eq1_norm") { idx_eqnorm_eq1_norm = i; break; }
  if (idx_eqnorm_eq1_norm < 0)
    error->all(FLERR, "GRACE-3L/KK: eqnorm 'eq1_norm' not found in loaded model");

  idx_rmsnorm_rho1_norm = -1;
  for (int i = 0; i < (int) m.rmsnorm.size(); i++)
    if (m.rmsnorm[i].name == "rho1_norm") { idx_rmsnorm_rho1_norm = i; break; }
  if (idx_rmsnorm_rho1_norm < 0)
    error->all(FLERR, "GRACE-3L/KK: rmsnorm 'rho1_norm' not found in loaded model");

  // ---- Task 2.8: locate the L2 SPBF "A2" (equivariant) by NAME (defensive;
  // spbf order is [A1,A2,A3]). Used by the radial-MLP + A2-assembly kernels. ----
  idx_spbf_A2 = -1;
  for (int i = 0; i < (int) m.spbf.size(); i++)
    if (m.spbf[i].name == "A2") { idx_spbf_A2 = i; break; }
  if (idx_spbf_A2 < 0)
    error->all(FLERR, "GRACE-3L/KK: SPBF 'A2' not found in loaded model");
  if (!m.spbf[idx_spbf_A2].equivariant)
    error->all(FLERR, "GRACE-3L/KK: SPBF 'A2' is not equivariant (unexpected schema)");

  // ---- Task 2.9: locate the remaining L2 chain by NAME (index-parallel with
  // d_prod/d_reduce/d_fc/d_eqnorm/d_rmsnorm, same convention as the Task 2.7
  // L1 block above). L2 differs from L1 in ONE structural way: every
  // downstream product/FC here is built on A2_red (a reduce of the raw SPBF
  // output A2), not on A2 itself — so there is one extra reduce (A2_red)
  // that L1 does not have. ----
  idx_reduce_A2_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A2_red") { idx_reduce_A2_red = i; break; }
  if (idx_reduce_A2_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_red' not found in loaded model");
  if (m.reduce[idx_reduce_A2_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A2_red].instr[0].name != "A2")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_red' instruction layout unexpected "
               "(expected single instruction named 'A2'; compute() wires d_A2 "
               "into inputs[0] positionally)");

  idx_prod_A2_2 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A2_2") { idx_prod_A2_2 = i; break; }
  if (idx_prod_A2_2 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A2_2' not found in loaded model");

  idx_reduce_A2_2_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A2_2_red") { idx_reduce_A2_2_red = i; break; }
  if (idx_reduce_A2_2_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_2_red' not found in loaded model");
  if (m.reduce[idx_reduce_A2_2_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A2_2_red].instr[0].name != "A2_2")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_2_red' instruction layout unexpected "
               "(expected single instruction named 'A2_2')");

  idx_fc_A2_2a = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A2_2a") { idx_fc_A2_2a = i; break; }
  if (idx_fc_A2_2a < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A2_2a' not found in loaded model");

  idx_prod_A2_3 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A2_3") { idx_prod_A2_3 = i; break; }
  if (idx_prod_A2_3 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A2_3' not found in loaded model");

  idx_reduce_A2_3_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A2_3_red") { idx_reduce_A2_3_red = i; break; }
  if (idx_reduce_A2_3_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_3_red' not found in loaded model");
  if (m.reduce[idx_reduce_A2_3_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A2_3_red].instr[0].name != "A2_3")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_3_red' instruction layout unexpected "
               "(expected single instruction named 'A2_3')");

  idx_fc_A2_2b = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A2_2b") { idx_fc_A2_2b = i; break; }
  if (idx_fc_A2_2b < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A2_2b' not found in loaded model");

  idx_prod_A2_4 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A2_4") { idx_prod_A2_4 = i; break; }
  if (idx_prod_A2_4 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A2_4' not found in loaded model");

  idx_reduce_A2_4_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A2_4_red") { idx_reduce_A2_4_red = i; break; }
  if (idx_reduce_A2_4_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_4_red' not found in loaded model");
  if (m.reduce[idx_reduce_A2_4_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A2_4_red].instr[0].name != "A2_4")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A2_4_red' instruction layout unexpected "
               "(expected single instruction named 'A2_4')");

  // eq2 / rho2: multi-instruction (4 collectors), elem-dependent reduces over
  // the SAME 4 named L2 tensors. Resolve role[j] by matching instr[j].name —
  // do NOT assume the npz's instruction_names order matches this enumeration.
  auto resolve_l2_roles = [&](const GRACE3LModel::Reduce &r, int roles[GRACE3L_MAX_REDUCE_INSTR]) {
    static const char *l2_names[4] = {"A2_red", "A2_2_red", "A2_3_red", "A2_4_red"};
    if ((int) r.instr.size() != 4)
      error->all(FLERR, "GRACE-3L/KK: reduce '{}' expected 4 L2 collectors, found {}",
                 r.name, (int) r.instr.size());
    for (int j = 0; j < 4; j++) {
      int role = -1;
      for (int k = 0; k < 4; k++)
        if (r.instr[j].name == l2_names[k]) { role = k; break; }
      if (role < 0)
        error->all(FLERR, "GRACE-3L/KK: reduce '{}' instruction[{}] name '{}' is not one "
                   "of the expected L2 inputs {{A2_red,A2_2_red,A2_3_red,A2_4_red}}",
                   r.name, j, r.instr[j].name);
      roles[j] = role;
    }
  };

  idx_reduce_eq2 = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "eq2") { idx_reduce_eq2 = i; break; }
  if (idx_reduce_eq2 < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'eq2' not found in loaded model");
  resolve_l2_roles(m.reduce[idx_reduce_eq2], eq2_input_role);

  idx_reduce_rho2 = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "rho2") { idx_reduce_rho2 = i; break; }
  if (idx_reduce_rho2 < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'rho2' not found in loaded model");
  resolve_l2_roles(m.reduce[idx_reduce_rho2], rho2_input_role);

  idx_eqnorm_eq2_norm = -1;
  for (int i = 0; i < (int) m.eqnorm.size(); i++)
    if (m.eqnorm[i].name == "eq2_norm") { idx_eqnorm_eq2_norm = i; break; }
  if (idx_eqnorm_eq2_norm < 0)
    error->all(FLERR, "GRACE-3L/KK: eqnorm 'eq2_norm' not found in loaded model");

  idx_rmsnorm_rho2_norm = -1;
  for (int i = 0; i < (int) m.rmsnorm.size(); i++)
    if (m.rmsnorm[i].name == "rho2_norm") { idx_rmsnorm_rho2_norm = i; break; }
  if (idx_rmsnorm_rho2_norm < 0)
    error->all(FLERR, "GRACE-3L/KK: rmsnorm 'rho2_norm' not found in loaded model");

  // ---- Task 2.10: locate the L3 SPBF "A3" (equivariant) by NAME (defensive;
  // spbf order is [A1,A2,A3]). Its indicator is eq2_norm (parity-doubled). ----
  idx_spbf_A3 = -1;
  for (int i = 0; i < (int) m.spbf.size(); i++)
    if (m.spbf[i].name == "A3") { idx_spbf_A3 = i; break; }
  if (idx_spbf_A3 < 0)
    error->all(FLERR, "GRACE-3L/KK: SPBF 'A3' not found in loaded model");
  if (!m.spbf[idx_spbf_A3].equivariant)
    error->all(FLERR, "GRACE-3L/KK: SPBF 'A3' is not equivariant (unexpected schema)");

  // ---- Guard: d_R1_nl/d_DR1_nl are allocated to L1_nradmax (A1's n_rad_max)
  // but reused as shared radial scratch by A2/A3 in TagComputeMLPRadial /
  // TagComputeMLPRadialDeriv. Safe only because every layer currently shares
  // n_rad_max; nothing else enforces that, so a model with a larger L2/L3
  // nradmax would silently overflow the shared scratch on device. ----
  if (m.spbf[idx_spbf_A2].n_rad_max > L1_nradmax)
    error->all(FLERR, "GRACE-3L/KK: A2.n_rad_max={} exceeds shared radial scratch "
               "size L1_nradmax={} (d_R1_nl/d_DR1_nl sized to L1_nradmax)",
               m.spbf[idx_spbf_A2].n_rad_max, L1_nradmax);
  if (m.spbf[idx_spbf_A3].n_rad_max > L1_nradmax)
    error->all(FLERR, "GRACE-3L/KK: A3.n_rad_max={} exceeds shared radial scratch "
               "size L1_nradmax={} (d_R1_nl/d_DR1_nl sized to L1_nradmax)",
               m.spbf[idx_spbf_A3].n_rad_max, L1_nradmax);

  // ---- Task 2.10: remaining L3 chain by NAME. Structurally identical to L2
  // (built on A3_red, a reduce of the raw SPBF output A3), EXCEPT L3 is
  // TERMINAL: it produces ONLY rho3/rho3_norm — there is NO eq3. ----
  idx_reduce_A3_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A3_red") { idx_reduce_A3_red = i; break; }
  if (idx_reduce_A3_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_red' not found in loaded model");
  if (m.reduce[idx_reduce_A3_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A3_red].instr[0].name != "A3")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_red' instruction layout unexpected "
               "(expected single instruction named 'A3'; compute() wires d_A3 "
               "into inputs[0] positionally)");

  idx_prod_A3_2 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A3_2") { idx_prod_A3_2 = i; break; }
  if (idx_prod_A3_2 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A3_2' not found in loaded model");

  idx_reduce_A3_2_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A3_2_red") { idx_reduce_A3_2_red = i; break; }
  if (idx_reduce_A3_2_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_2_red' not found in loaded model");
  if (m.reduce[idx_reduce_A3_2_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A3_2_red].instr[0].name != "A3_2")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_2_red' instruction layout unexpected "
               "(expected single instruction named 'A3_2')");

  idx_fc_A3_2a = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A3_2a") { idx_fc_A3_2a = i; break; }
  if (idx_fc_A3_2a < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A3_2a' not found in loaded model");

  idx_prod_A3_3 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A3_3") { idx_prod_A3_3 = i; break; }
  if (idx_prod_A3_3 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A3_3' not found in loaded model");

  idx_reduce_A3_3_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A3_3_red") { idx_reduce_A3_3_red = i; break; }
  if (idx_reduce_A3_3_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_3_red' not found in loaded model");
  if (m.reduce[idx_reduce_A3_3_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A3_3_red].instr[0].name != "A3_3")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_3_red' instruction layout unexpected "
               "(expected single instruction named 'A3_3')");

  idx_fc_A3_2b = -1;
  for (int i = 0; i < (int) m.fc.size(); i++)
    if (m.fc[i].name == "A3_2b") { idx_fc_A3_2b = i; break; }
  if (idx_fc_A3_2b < 0)
    error->all(FLERR, "GRACE-3L/KK: FC 'A3_2b' not found in loaded model");

  idx_prod_A3_4 = -1;
  for (int i = 0; i < (int) m.prod.size(); i++)
    if (m.prod[i].name == "A3_4") { idx_prod_A3_4 = i; break; }
  if (idx_prod_A3_4 < 0)
    error->all(FLERR, "GRACE-3L/KK: product 'A3_4' not found in loaded model");

  idx_reduce_A3_4_red = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "A3_4_red") { idx_reduce_A3_4_red = i; break; }
  if (idx_reduce_A3_4_red < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_4_red' not found in loaded model");
  if (m.reduce[idx_reduce_A3_4_red].instr.size() != 1 ||
      m.reduce[idx_reduce_A3_4_red].instr[0].name != "A3_4")
    error->all(FLERR, "GRACE-3L/KK: reduce 'A3_4_red' instruction layout unexpected "
               "(expected single instruction named 'A3_4')");

  // rho3: multi-instruction (4 collectors), elem-dependent, only_invar reduce
  // over the SAME 4 named L3 tensors. Resolve role[j] by matching instr[j].name
  // — do NOT assume the npz's instruction_names order matches this enumeration.
  // (L3 is terminal: no eq3, so only rho3 needs a role array.)
  auto resolve_l3_roles = [&](const GRACE3LModel::Reduce &r, int roles[GRACE3L_MAX_REDUCE_INSTR]) {
    static const char *l3_names[4] = {"A3_red", "A3_2_red", "A3_3_red", "A3_4_red"};
    if ((int) r.instr.size() != 4)
      error->all(FLERR, "GRACE-3L/KK: reduce '{}' expected 4 L3 collectors, found {}",
                 r.name, (int) r.instr.size());
    for (int j = 0; j < 4; j++) {
      int role = -1;
      for (int k = 0; k < 4; k++)
        if (r.instr[j].name == l3_names[k]) { role = k; break; }
      if (role < 0)
        error->all(FLERR, "GRACE-3L/KK: reduce '{}' instruction[{}] name '{}' is not one "
                   "of the expected L3 inputs {{A3_red,A3_2_red,A3_3_red,A3_4_red}}",
                   r.name, j, r.instr[j].name);
      roles[j] = role;
    }
  };

  idx_reduce_rho3 = -1;
  for (int i = 0; i < (int) m.reduce.size(); i++)
    if (m.reduce[i].name == "rho3") { idx_reduce_rho3 = i; break; }
  if (idx_reduce_rho3 < 0)
    error->all(FLERR, "GRACE-3L/KK: reduce 'rho3' not found in loaded model");
  resolve_l3_roles(m.reduce[idx_reduce_rho3], rho3_input_role);

  idx_rmsnorm_rho3_norm = -1;
  for (int i = 0; i < (int) m.rmsnorm.size(); i++)
    if (m.rmsnorm[i].name == "rho3_norm") { idx_rmsnorm_rho3_norm = i; break; }
  if (idx_rmsnorm_rho3_norm < 0)
    error->all(FLERR, "GRACE-3L/KK: rmsnorm 'rho3_norm' not found in loaded model");

  // ---- UQ basis-RP: pin the per-layer scalar-basis block widths. The GMM feature
  // is the projection of the invariant (l=0) products entering rho1|rho2|rho3, read
  // in the exact (source, n, collect) order TagReduceN iterates them, so R's row
  // count D_basis must equal the summed per-block widths [Σ n_in·n_conn per source].
  // R-row layout is [rho1, rho2, rho3] (sorted by reduce name). This guards against
  // an artifact/model mismatch. ----
  if (m.has_uq) {
    auto block_width = [&](int idx_reduce) {
      const auto &r = m.reduce[idx_reduce];
      int w = 0;
      for (const auto &ins : r.instr) w += ins.n_in * ins.n_conn;
      return w;
    };
    uq_d_basis_rho1 = block_width(idx_reduce_rho1);
    uq_d_basis_rho2 = block_width(idx_reduce_rho2);
    uq_d_basis_rho3 = block_width(idx_reduce_rho3);
    if (uq_d_basis_rho1 + uq_d_basis_rho2 + uq_d_basis_rho3 != m.uq_d_basis)
      error->all(FLERR, "GRACE-3L/KK UQ: scalar-basis width rho1({})+rho2({})+rho3({}) "
                 "!= uq_rp_matrix rows ({})",
                 uq_d_basis_rho1, uq_d_basis_rho2, uq_d_basis_rho3, m.uq_d_basis);
  }

  // ---- Task 2.8/2.10: inter-layer message comm sizing. Round 1 forward-comms
  // the Layer-2 indicator eq1_norm [nall, n_out=64, M=25] (I_n_funcs*I_n_out =
  // 1600/atom); round 2 forward-comms the Layer-3 indicator eq2_norm
  // [nall, n_out=64, M=50] (I2_n_funcs*I2_n_out = 3200/atom, PARITY-DOUBLED).
  // comm_forward is the MAX of the two payloads so the shared comm buffer fits
  // both stages; comm_stage selects which at each forward_comm. comm_reverse is
  // ALSO the MAX (Task 3.2): round-1 reverse-comms grad_I [I_n_funcs*I_n_out=1600]
  // (Task 3.3), round-2 reverse-comms the eq2_norm message adj [I2_n_funcs*
  // I2_n_out=3200]; the shared reverse buffer must fit the larger stage. ----
  I_n_funcs = m.eqnorm[idx_eqnorm_eq1_norm].n_out;
  I_n_out   = m.eqnorm[idx_eqnorm_eq1_norm].M;
  I2_n_funcs = m.eqnorm[idx_eqnorm_eq2_norm].n_out;
  I2_n_out   = m.eqnorm[idx_eqnorm_eq2_norm].M;
  comm_forward = std::max(I_n_funcs * I_n_out, I2_n_funcs * I2_n_out);
  comm_reverse = std::max(I_n_funcs * I_n_out, I2_n_funcs * I2_n_out);

  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-3L/KK] Initialization complete: weights on device "
                   "(SPBF={} products={} reduces={} FCs={} eqnorms={} rmsnorms={}) | "
                   "comm_forward={} (round1 I_n_funcs={} I_n_out={}; "
                   "round2 I2_n_funcs={} I2_n_out={}).\n",
                   n_spbf, n_prods, n_reduces, n_fcs, n_eqnorms, n_rmsnorms,
                   comm_forward, I_n_funcs, I_n_out, I2_n_funcs, I2_n_out);

#ifdef KOKKOS_ENABLE_CUDA
  // ---- Opt 15: create the cuBLAS handle once, bind to the SAME stream the
  // Kokkos parallel_fors launch on (the default Kokkos::Cuda instance), and force
  // TRUE fp32 (CUBLAS_PEDANTIC_MATH disables the Ampere+ tf32 default). Same-stream
  // binding makes the GEMMs and the surrounding Kokkos assembly/activation kernels
  // execute in issue order with NO race and NO explicit fences. ----
  if (!cublas_handle) {
    if (cublasCreate(&cublas_handle) != CUBLAS_STATUS_SUCCESS)
      error->all(FLERR, "GRACE-3L/KK: cublasCreate failed");
    cublasSetStream(cublas_handle, DeviceType().cuda_stream());
    cublasSetMathMode(cublas_handle, CUBLAS_PEDANTIC_MATH);
    mlp_cublas_selftest();   // tiny-case transpose/leading-dim gate (aborts on mismatch)
  }
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-3L/KK] radial MLP: true-fp32 cuBLAS SGEMM enabled (CUDA)\n");
#else
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-3L/KK] radial MLP: hand-written Kokkos kernels "
                   "(cuBLAS is CUDA-only; default on non-CUDA backends, e.g. HIP/AMD)\n");
#endif
}

// ======================================================================
// init_one
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
double PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::init_one(int i, int j)
{
  double rcut_ij = rcut;
  if (grace_model->has_bond_specific_cutoff) {
    int mu_i = h_type2model[i];
    int mu_j = h_type2model[j];
    if (mu_i >= 0 && mu_j >= 0)
      rcut_ij = grace_model->bond_cutoff_map[mu_i * nelements + mu_j];
  }
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = rcut_ij * rcut_ij;
  k_cutsq.modify_host();
  return rcut_ij;
}

// ======================================================================
// precompute_harmonics — GeomScalar recurrence coefficients for the real
// spherical harmonics. COPIED VERBATIM from pair_grace_2l_kokkos (same
// alm/blm/cl/dl/idx_sph convention, consumed by the geometry kernels).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::precompute_harmonics()
{
  auto h_idx_sph = Kokkos::create_mirror_view(d_idx_sph);
  auto h_alm = Kokkos::create_mirror_view(alm);
  auto h_blm = Kokkos::create_mirror_view(blm);
  auto h_cl = Kokkos::create_mirror_view(cl);
  auto h_dl = Kokkos::create_mirror_view(dl);

  Kokkos::deep_copy(h_idx_sph, -1);

  int idx_sph = 0;
  for (int m_val = 0; m_val <= lmax; m_val++) {
    const double msq = m_val * m_val;
    for (int l = m_val; l <= lmax; l++) {
      const int idx = l * (l + 1) + m_val;
      h_idx_sph(idx) = idx_sph;
      double a = 0.0, b = 0.0;
      if (l > 1 && m_val < l - 1) {
        const double lsq = l * l;
        const double ld = 2 * l;
        const double l1 = (4 * lsq - 1);
        const double l2 = lsq - ld + 1;
        a = sqrt(double(l1) / double(lsq - msq));
        b = -sqrt(double(l2 - msq) / double(4 * l2 - 1));
      }
      h_alm(idx_sph) = a;
      h_blm(idx_sph) = b;
      idx_sph++;
    }
  }
  idx_sph_max = idx_sph;

  for (int l = 1; l <= lmax; l++) {
    h_cl(l) = -sqrt(1.0 + 0.5 / double(l));
    h_dl(l) = sqrt(double(2 * (l - 1) + 3));
  }

  Kokkos::deep_copy(d_idx_sph, h_idx_sph);
  Kokkos::deep_copy(alm, h_alm);
  Kokkos::deep_copy(blm, h_blm);
  Kokkos::deep_copy(cl, h_cl);
  Kokkos::deep_copy(dl, h_dl);
}

// ======================================================================
// copy_weights_to_device — stage every host GRACE3LModel tensor into the
// constant device views. RAW weight + separate norm scalar (kernel applies
// once), EXCEPT the SPBF per-element species tables (z_tr / z_proj) which
// fold the layer norm in once at precompute time.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::copy_weights_to_device()
{
  auto &m = *grace_model;

  // ---- generic copy helpers (dst value type drives the double->NN/Geom cast) ----
  auto copy_1d_int = [](const std::vector<int> &src, auto &dst, const char *label) {
    MemKK::realloc_kokkos(dst, label, std::max<size_t>(src.size(), 1));
    auto h = Kokkos::create_mirror_view(dst);
    for (size_t i = 0; i < src.size(); i++) h(i) = src[i];
    Kokkos::deep_copy(dst, h);
  };
  auto copy_1d = [](const std::vector<double> &src, auto &dst, const char *label) {
    MemKK::realloc_kokkos(dst, label, std::max<size_t>(src.size(), 1));
    auto h = Kokkos::create_mirror_view(dst);
    for (size_t i = 0; i < src.size(); i++) h(i) = src[i];
    Kokkos::deep_copy(dst, h);
  };
  auto copy_2d = [this](const std::vector<double> &src, auto &dst, int d0, int d1, const char *label) {
    if (src.size() != (size_t) d0 * d1)
      error->all(FLERR, "GRACE-3L/KK: {} size {} does not match expected {}x{}={} "
                 "(desynced weight export?)", label, src.size(), d0, d1, (size_t) d0 * d1);
    MemKK::realloc_kokkos(dst, label, d0, d1);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        h(i, j) = src[(size_t) i * d1 + j];
    Kokkos::deep_copy(dst, h);
  };
  auto copy_3d = [this](const std::vector<double> &src, auto &dst, int d0, int d1, int d2, const char *label) {
    if (src.size() != (size_t) d0 * d1 * d2)
      error->all(FLERR, "GRACE-3L/KK: {} size {} does not match expected {}x{}x{}={} "
                 "(desynced weight export?)", label, src.size(), d0, d1, d2, (size_t) d0 * d1 * d2);
    MemKK::realloc_kokkos(dst, label, d0, d1, d2);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        for (int k = 0; k < d2; k++)
          h(i, j, k) = src[((size_t) i * d1 + j) * d2 + k];
    Kokkos::deep_copy(dst, h);
  };
  auto copy_4d = [this](const std::vector<double> &src, auto &dst, int d0, int d1, int d2, int d3, const char *label) {
    if (src.size() != (size_t) d0 * d1 * d2 * d3)
      error->all(FLERR, "GRACE-3L/KK: {} size {} does not match expected {}x{}x{}x{}={} "
                 "(desynced weight export?)", label, src.size(), d0, d1, d2, d3, (size_t) d0 * d1 * d2 * d3);
    MemKK::realloc_kokkos(dst, label, d0, d1, d2, d3);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        for (int k = 0; k < d2; k++)
          for (int l = 0; l < d3; l++)
            h(i, j, k, l) = src[(((size_t) i * d1 + j) * d2 + k) * d3 + l];
    Kokkos::deep_copy(dst, h);
  };

  // ---- Chemical embedding / shifts / bond cutoffs / output scale ----
  copy_2d(m.chem_embedding, d_chem_embed, m.n_elements, m.embedding_size, "g3l:chem_embed");
  copy_1d(m.shift_values, d_shifts, "g3l:shifts");
  copy_2d(m.bond_cutoff_map, d_bond_cutoff, m.n_elements, m.n_elements, "g3l:bond_cutoff");
  output_scale = (NNScalar) m.output_scale;

  // ---- SPBF layers (radial MLP w/ hidden bias; species modulation; CG couple) ----
  n_spbf = (int) m.spbf.size();
  for (int si = 0; si < n_spbf; si++) {
    const auto &s = m.spbf[si];
    DeviceSPBF &ds = d_spbf[si];
    ds.equivariant = s.equivariant ? 1 : 0;
    ds.n_rad_max = s.n_rad_max; ds.n_rad_basis = s.n_rad_basis;
    ds.spbf_lmax = s.lmax; ds.Lmax = s.Lmax; ds.p = s.p;
    ds.nfunc = s.nfunc; ds.n_lm_ind = s.n_lm_ind;
    ds.rcut = (GeomScalar) s.rcut;
    ds.inv_avg_n_neigh = (GeomScalar) s.inv_avg_n_neigh;
    copy_1d_int(s.l_tile, ds.l_tile, "g3l:spbf_l_tile");

    // radial MLP (padded), with hidden-layer biases (output row zero)
    const auto &mlp = s.mlp;
    int nl = mlp.n_layers, max_in = 0, max_out = 0;
    for (int i = 0; i < nl; i++) { max_in = std::max(max_in, mlp.n_in[i]); max_out = std::max(max_out, mlp.n_out[i]); }
    ds.mlp_n_layers = nl;
    ds.mlp_max_dim = std::max(max_in, max_out);
    std::string sp = "g3l:" + s.name + "_mlp";
    MemKK::realloc_kokkos(ds.mlp_W, (sp + "_W").c_str(), std::max(nl, 1), std::max(max_in, 1), std::max(max_out, 1));
    MemKK::realloc_kokkos(ds.mlp_b, (sp + "_b").c_str(), std::max(nl, 1), std::max(max_out, 1));
    auto hW = Kokkos::create_mirror_view(ds.mlp_W);
    auto hb = Kokkos::create_mirror_view(ds.mlp_b);
    Kokkos::deep_copy(hW, (NNScalar) 0); Kokkos::deep_copy(hb, (NNScalar) 0);
    for (int l = 0; l < nl; l++) {
      for (int i = 0; i < mlp.n_in[l]; i++)
        for (int j = 0; j < mlp.n_out[l]; j++)
          hW(l, i, j) = (NNScalar) mlp.W[l][(size_t) i * mlp.n_out[l] + j];
      if (!mlp.b[l].empty())
        for (int j = 0; j < mlp.n_out[l] && j < (int) mlp.b[l].size(); j++)
          hb(l, j) = (NNScalar) mlp.b[l][j];
    }
    Kokkos::deep_copy(ds.mlp_W, hW); Kokkos::deep_copy(ds.mlp_b, hb);
#ifdef KOKKOS_ENABLE_CUDA
    // Opt 15: pack each layer's weights TIGHTLY row-major [n_in x n_out] for cuBLAS
    // (the padded LayoutLeft ds.mlp_W is not a contiguous per-layer block). lda=n_out.
    for (int l = 0; l < nl; l++) {
      const int ni = mlp.n_in[l], no = mlp.n_out[l];
      MemKK::realloc_kokkos(d_mlp_Wg[si][l], (sp + "_Wg").c_str(),
                            std::max(ni, 1), std::max(no, 1));
      auto hWg = Kokkos::create_mirror_view(d_mlp_Wg[si][l]);
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < no; j++)
          hWg(i, j) = (NNScalar) mlp.W[l][(size_t) i * no + j];
      Kokkos::deep_copy(d_mlp_Wg[si][l], hWg);
    }
#endif
    MemKK::realloc_kokkos(ds.mlp_norms, (sp + "_norms").c_str(), std::max(nl, 1));
    auto hn = Kokkos::create_mirror_view(ds.mlp_norms);
    for (int i = 0; i < nl; i++) hn(i) = (NNScalar) mlp.norm[i];
    Kokkos::deep_copy(ds.mlp_norms, hn);
    MemKK::realloc_kokkos(ds.mlp_dims, (sp + "_dims").c_str(), nl + 1);
    auto hd = Kokkos::create_mirror_view(ds.mlp_dims);
    hd(0) = mlp.n_in.empty() ? 0 : mlp.n_in[0];
    for (int i = 0; i < nl; i++) hd(i + 1) = mlp.n_out[i];
    Kokkos::deep_copy(ds.mlp_dims, hd);

    // Per-element species modulation — norm FOLDED here (once).
    if (!s.equivariant) {
      // A1: z_tr[mu][n] = norm * sum_e chem_embed[mu,e] * lin_transform_W[e,n]
      MemKK::realloc_kokkos(ds.z_tr, ("g3l:" + s.name + "_z_tr").c_str(), m.n_elements, s.n_rad_max);
      auto h = Kokkos::create_mirror_view(ds.z_tr);
      for (int mu = 0; mu < m.n_elements; mu++)
        for (int n = 0; n < s.n_rad_max; n++) {
          double v = 0.0;
          for (int e = 0; e < m.embedding_size; e++)
            v += m.chem_embedding[(size_t) mu * m.embedding_size + e] *
                 s.lin_transform_W[(size_t) e * s.n_rad_max + n];
          h(mu, n) = (NNScalar) (v * s.lin_transform_norm);
        }
      Kokkos::deep_copy(ds.z_tr, h);
    } else {
      // A2/A3: z_proj[mu][n] = norm * sum_e chem_embed[mu,e] * chem_linear_W[e,n]
      const int cn_out = (int) (s.chem_linear_W.size() / std::max(m.embedding_size, 1));
      MemKK::realloc_kokkos(ds.z_proj, ("g3l:" + s.name + "_z_proj").c_str(), m.n_elements, cn_out);
      auto h = Kokkos::create_mirror_view(ds.z_proj);
      for (int mu = 0; mu < m.n_elements; mu++)
        for (int n = 0; n < cn_out; n++) {
          double v = 0.0;
          for (int e = 0; e < m.embedding_size; e++)
            v += m.chem_embedding[(size_t) mu * m.embedding_size + e] *
                 s.chem_linear_W[(size_t) e * cn_out + n];
          h(mu, n) = (NNScalar) (v * s.chem_linear_norm);
        }
      Kokkos::deep_copy(ds.z_proj, h);
      copy_1d(s.chem_l0_mask, ds.chem_l0_mask, "g3l:spbf_chem_l0_mask");
      // dense CG couple: [(lmax+1)^2 * n_lm_ind, nfunc]
      const int cg_rows = (int) (s.cg_W.size() / std::max(s.nfunc, 1));
      copy_2d(s.cg_W, ds.cg_W, cg_rows, s.nfunc, "g3l:spbf_cg_W");
      // Opt 9: build sparse CSC/CSR forms of cg_W. CG selection rules make the
      // matrix ~99% zeros (mean ~3-6 nnz per row/column), so the dense couple
      // loops in TagComputeSPBFFused (per-f over P) and TagComputeAdjProd
      // (per-p over f) wasted ~100-200x compute + L2 bandwidth streaming zeros.
      // Ascending index order within each column/row preserves the dense loops'
      // accumulation order exactly (skipped terms are exact +0 contributions),
      // so results are numerically identical. Worst case (a dense cg_W) the
      // sparse loop only adds index-load overhead — still correct.
      {
        const int P = cg_rows, NF = s.nfunc;
        std::vector<int> cptr(NF + 1, 0), ridx;
        std::vector<double> cval;
        std::vector<int> rptr(P + 1, 0), cidx;
        std::vector<double> rval;
        for (int f = 0; f < NF; f++) {
          cptr[f] = (int) ridx.size();
          for (int p = 0; p < P; p++) {
            const double w = s.cg_W[(size_t) p * NF + f];
            if (w != 0.0) { ridx.push_back(p); cval.push_back(w); }
          }
        }
        cptr[NF] = (int) ridx.size();
        for (int p = 0; p < P; p++) {
          rptr[p] = (int) cidx.size();
          for (int f = 0; f < NF; f++) {
            const double w = s.cg_W[(size_t) p * NF + f];
            if (w != 0.0) { cidx.push_back(f); rval.push_back(w); }
          }
        }
        rptr[P] = (int) cidx.size();
        copy_1d_int(cptr, ds.cg_col_ptr, "g3l:spbf_cg_col_ptr");
        copy_1d_int(ridx, ds.cg_row_idx, "g3l:spbf_cg_row_idx");
        copy_1d(cval, ds.cg_col_val, "g3l:spbf_cg_col_val");
        copy_1d_int(rptr, ds.cg_row_ptr, "g3l:spbf_cg_row_ptr");
        copy_1d_int(cidx, ds.cg_col_idx, "g3l:spbf_cg_col_idx");
        copy_1d(rval, ds.cg_row_val, "g3l:spbf_cg_row_val");
      }
    }
  }

  // ---- cp_l products ----
  n_prods = (int) m.prod.size();
  for (int pi = 0; pi < n_prods; pi++) {
    const auto &c = m.prod[pi];
    DeviceCPL &dc = d_prod[pi];
    dc.rank = c.rank; dc.nfunc = c.nfunc; dc.n_left = c.n_left; dc.n_right = c.n_right;
    dc.n_groups_left = c.n_groups_left; dc.n_groups_right = c.n_groups_right; dc.n_cg = c.n_cg;
    dc.norm_u = (NNScalar) c.norm_u; dc.norm_v = (NNScalar) c.norm_v;
    copy_3d(c.U, dc.U, c.n_groups_left, c.rank, c.n_left, "g3l:prod_U");
    copy_3d(c.V, dc.V, c.n_groups_right, c.rank, c.n_right, "g3l:prod_V");
    copy_1d(c.cg, dc.cg, "g3l:prod_cg");
    copy_1d_int(c.group_left, dc.group_left, "g3l:prod_gl");
    copy_1d_int(c.group_right, dc.group_right, "g3l:prod_gr");
    copy_1d_int(c.left_ind, dc.left_ind, "g3l:prod_li");
    copy_1d_int(c.right_ind, dc.right_ind, "g3l:prod_ri");
    copy_1d_int(c.m_sum_ind, dc.m_sum_ind, "g3l:prod_si");
    // Opt 14: pack (m_sum_ind, right_ind, left_ind) into one int32 per entry so
    // the fused couple loops stream half the index bytes. Hard-check the packing
    // bounds (li/ri < 256, m < 65536) — far beyond any plausible cp_l dims
    // (nlm <= (Lmax+1)^2, nfunc <= a few thousand).
    {
      std::vector<int> packed(std::max(c.n_cg, 1), 0);
      for (int i = 0; i < c.n_cg; i++) {
        const int li = c.left_ind[i], ri = c.right_ind[i], mi = c.m_sum_ind[i];
        if (li < 0 || li > 255 || ri < 0 || ri > 255 || mi < 0 || mi > 65535)
          error->all(FLERR, "GRACE-3L/KK: cp_l couple index out of packing range "
                     "(li={}, ri={}, m={})", li, ri, mi);
        packed[i] = (mi << 16) | (ri << 8) | li;
      }
      copy_1d_int(packed, dc.cg_packed, "g3l:prod_cg_packed");
    }
    copy_1d_int(c.out_l, dc.out_l, "g3l:prod_out_l");
    copy_1d_int(c.out_parity, dc.out_parity, "g3l:prod_out_parity");
#ifdef KOKKOS_ENABLE_CUDA
    // Opt 16: pre-gather + transpose U/V into cuBLAS-friendly per-lm-slab weights
    // for the batched CPLBwdProject GEMM. Uw(w,n,r) = U(group_left(w), r, n) so each
    // w-slab is a tight col-major [rank x n_left] matrix (ldb=rank, w-stride=n_left*rank).
    // Same for Vw. One-time; ~a few MB total across all products.
    {
      const int nlm_l = (int) c.group_left.size();
      const int nlm_r = (int) c.group_right.size();
      MemKK::realloc_kokkos(dc.Uw, "g3l:prod_Uw", std::max(nlm_l, 1),
                            std::max(c.n_left, 1), std::max(c.rank, 1));
      MemKK::realloc_kokkos(dc.Vw, "g3l:prod_Vw", std::max(nlm_r, 1),
                            std::max(c.n_right, 1), std::max(c.rank, 1));
      auto hUw = Kokkos::create_mirror_view(dc.Uw);
      auto hVw = Kokkos::create_mirror_view(dc.Vw);
      for (int w = 0; w < nlm_l; w++) {
        const int g = c.group_left[w];
        for (int n = 0; n < c.n_left; n++)
          for (int r = 0; r < c.rank; r++)
            hUw(w, n, r) = (NNScalar) c.U[((size_t) g * c.rank + r) * c.n_left + n];
      }
      for (int w = 0; w < nlm_r; w++) {
        const int g = c.group_right[w];
        for (int n = 0; n < c.n_right; n++)
          for (int r = 0; r < c.rank; r++)
            hVw(w, n, r) = (NNScalar) c.V[((size_t) g * c.rank + r) * c.n_right + n];
      }
      Kokkos::deep_copy(dc.Uw, hUw);
      Kokkos::deep_copy(dc.Vw, hVw);
    }
#endif
  }

  // cp_l scratch caps: max rank / max n_lm over ALL products so the two shared
  // projection scratch views (d_cpl_lproj/rproj) serve every compute_cp_l call.
  cpl_scratch_rank = 0;
  cpl_scratch_nlm = 0;
  for (int pi = 0; pi < n_prods; pi++) {
    const auto &c = m.prod[pi];
    cpl_scratch_rank = std::max(cpl_scratch_rank, c.rank);
    cpl_scratch_nlm = std::max(cpl_scratch_nlm,
                               std::max((int) c.group_left.size(), (int) c.group_right.size()));
  }

  // ---- FunctionReduceN ----
  n_reduces = (int) m.reduce.size();
  for (int ri = 0; ri < n_reduces; ri++) {
    const auto &r = m.reduce[ri];
    DeviceReduce &dr = d_reduce[ri];
    dr.n_out = r.n_out; dr.n_funcs = r.n_funcs;
    dr.only_invar = r.only_invar ? 1 : 0;
    dr.elem_dep = r.elem_dep ? 1 : 0;
    dr.n_instr = (int) r.instr.size();
    const int n_types = r.elem_dep ? m.n_elements : 1;
    for (int ci = 0; ci < dr.n_instr; ci++) {
      const auto &ins = r.instr[ci];
      dr.n_in[ci] = ins.n_in; dr.w_shape[ci] = ins.w_shape; dr.n_conn[ci] = ins.n_conn;
      dr.norm[ci] = (NNScalar) ins.norm;
      // W stored uniformly 4D [n_types, n_out, n_in, w_shape] (n_types=1 if !elem_dep).
      copy_4d(ins.W, dr.W[ci], n_types, r.n_out, ins.n_in, ins.w_shape, "g3l:reduce_W");
      copy_1d_int(ins.collect_ind, dr.collect_ind[ci], "g3l:reduce_ci");
      copy_1d_int(ins.w_l_tile, dr.w_l_tile[ci], "g3l:reduce_wlt");
      copy_1d_int(ins.total_sum_ind, dr.total_sum_ind[ci], "g3l:reduce_tsi");
    }
    dr.has_norm_map = r.has_norm_map ? 1 : 0;
    if (r.has_norm_map) copy_1d(r.norm_map, dr.norm_map, "g3l:reduce_norm_map");
  }

#ifdef KOKKOS_ENABLE_CUDA
  // ---- Opt 17: cuBLAS batched-GEMM scratch for the element-independent, non-only_invar,
  //      single-instruction product reduces (A*_*_red). Fold norm*norm_map(f) into a
  //      per-connection col-major [n_out x n_in] weight block Wsc[c] (serves BOTH the
  //      forward rounds and the backward single batch), and precompute the forward round
  //      schedule (connections grouped so every round writes DISTINCT output f -> the
  //      beta=1 batched accumulate into `out` is race-free within a round). collect_ind must
  //      be a bijection over the input lm-columns so the backward writes DISTINCT adj_in
  //      slabs -> a single beta=1 batch with no scatter; this is verified per reduce below
  //      and cuBLAS is disabled (hand-kernel fallback) for any reduce that violates it. ----
  reduceN_max_conn = 0;
  for (int ri = 0; ri < n_reduces; ri++) {
    const auto &r = m.reduce[ri];
    ReduceHost &hr = h_reduce[ri];
    hr.cublas_ok = (!r.elem_dep && !r.only_invar && r.instr.size() == 1);
    if (!hr.cublas_ok) continue;
    const auto &ins = r.instr[0];
    // The cuBLAS backward accumulates into adj_in slabs indexed by collect_ind in a single
    // beta=1 batch, which is only race-free if collect_ind is a bijection (distinct value
    // per connection). Verify it here; if any index repeats, disable cuBLAS for this reduce
    // and fall back to the hand kernel (which scatters safely).
    {
      std::unordered_set<int> seen_ci;
      bool ci_unique = true;
      for (int c = 0; c < ins.n_conn; c++)
        if (!seen_ci.insert(ins.collect_ind[c]).second) { ci_unique = false; break; }
      if (!ci_unique) { hr.cublas_ok = false; continue; }
    }
    const int n_out = r.n_out, n_in = ins.n_in, n_conn = ins.n_conn, w_shape = ins.w_shape;
    const int n_f = r.n_funcs;
    hr.n_conn = n_conn; hr.n_out = n_out; hr.n_in = n_in;
    hr.ci = ins.collect_ind;                 // [n_conn]
    hr.tsi = ins.total_sum_ind;              // [n_conn]
    // per-connection scaled weight blocks: Wsc[c](k,n) = norm*norm_map(tsi[c])*W(0,k,n,tile)
    // stored col-major [n_out x n_in] (k fastest), c-stride = n_out*n_in.
    MemKK::realloc_kokkos(d_reduce[ri].Wsc, "g3l:reduce_Wsc", (size_t) n_conn * n_out * n_in);
    auto hWsc = Kokkos::create_mirror_view(d_reduce[ri].Wsc);
    for (int c = 0; c < n_conn; c++) {
      const int tile = ins.w_l_tile[c];
      const int f = hr.tsi[c];
      const double sc = ins.norm * (r.has_norm_map ? r.norm_map[f] : 1.0);
      const size_t base = (size_t) c * n_out * n_in;
      for (int n = 0; n < n_in; n++)
        for (int k = 0; k < n_out; k++)
          hWsc(base + k + (size_t) n_out * n) =
              (NNScalar) (sc * ins.W[((size_t) k * n_in + n) * w_shape + tile]);
    }
    Kokkos::deep_copy(d_reduce[ri].Wsc, hWsc);
    // forward round schedule: round rr collects the rr-th connection of every output f.
    std::vector<std::vector<int>> byf(n_f);
    for (int c = 0; c < n_conn; c++) byf[hr.tsi[c]].push_back(c);
    int maxr = 0;
    for (int f = 0; f < n_f; f++) maxr = std::max(maxr, (int) byf[f].size());
    if (maxr > GRACE3L_MAX_REDUCE_ROUNDS)
      error->all(FLERR, "GRACE-3L/KK: reduce {} needs {} cuBLAS fwd rounds > cap {}",
                 r.name, maxr, GRACE3L_MAX_REDUCE_ROUNDS);
    hr.fwd_order.clear();
    hr.round_ptr[0] = 0;
    for (int rr = 0; rr < maxr; rr++) {
      for (int f = 0; f < n_f; f++)
        if ((int) byf[f].size() > rr) hr.fwd_order.push_back(byf[f][rr]);
      hr.round_ptr[rr + 1] = (int) hr.fwd_order.size();
    }
    hr.nrounds = maxr;
    reduceN_max_conn = std::max(reduceN_max_conn, n_conn);
  }
  if (reduceN_max_conn > 0) {
    MemKK::realloc_kokkos(d_redptr, "g3l:reduceN_ptrbuf", (size_t) 3 * reduceN_max_conn);
    h_redptr = Kokkos::View<uintptr_t*, Kokkos::HostSpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "g3l:reduceN_hptr"),
        (size_t) 3 * reduceN_max_conn);
  }
#endif

  // ---- FCRight2Left ----
  n_fcs = (int) m.fc.size();
  for (int fi = 0; fi < n_fcs; fi++) {
    const auto &f = m.fc[fi];
    DeviceFC &df = d_fc[fi];
    df.n_out = f.n_out; df.left_coefs = f.left_coefs ? 1 : 0;
    df.n_funcs_left = f.n_funcs_left; df.n_funcs_right = f.n_funcs_right;
    df.n_in_left = f.n_in_left; df.n_in_right = f.n_in_right;
    df.w_shape_left = f.w_shape_left; df.w_shape_right = f.w_shape_right;
    df.norm_left = (NNScalar) f.norm_left; df.norm_right = (NNScalar) f.norm_right;
    copy_3d(f.w_right, df.w_right, f.n_out, f.n_in_right, f.w_shape_right, "g3l:fc_wr");
    if (f.left_coefs && !f.w_left.empty())
      copy_3d(f.w_left, df.w_left, f.n_out, f.n_in_left, f.w_shape_left, "g3l:fc_wl");
    copy_1d_int(f.w_tile_left, df.w_tile_left, "g3l:fc_wtl");
    copy_1d_int(f.w_tile_right, df.w_tile_right, "g3l:fc_wtr");
    copy_1d_int(f.collect_to, df.collect_to, "g3l:fc_ct");
    copy_1d_int(f.collect_from, df.collect_from, "g3l:fc_cf");
    copy_1d(f.norm_out_factor, df.norm_out_factor, "g3l:fc_nof");
  }

  // ---- EquivariantRMSNorm ----
  n_eqnorms = (int) m.eqnorm.size();
  for (int ei = 0; ei < n_eqnorms; ei++) {
    const auto &e = m.eqnorm[ei];
    DeviceEqNorm &de = d_eqnorm[ei];
    de.center_l0 = e.center_l0; de.n_groups = e.n_groups; de.n_out = e.n_out; de.M = e.M;
    de.eps = (NNScalar) e.eps;
    copy_2d(e.affine_weight, de.affine_weight, e.n_groups, e.n_out, "g3l:eqnorm_aff");
    copy_1d(e.degree_weights, de.degree_weights, "g3l:eqnorm_dw");
    copy_1d(e.l0_mask, de.l0_mask, "g3l:eqnorm_l0");
    copy_1d_int(e.expand_index, de.expand_index, "g3l:eqnorm_exp");
  }

  // ---- InvariantLayerRMSNorm ----
  n_rmsnorms = (int) m.rmsnorm.size();
  for (int ni = 0; ni < n_rmsnorms; ni++) {
    const auto &r = m.rmsnorm[ni];
    DeviceRMSNorm &dn = d_rmsnorm[ni];
    dn.type = r.type; dn.n_out = r.n_out; dn.scale_len = (int) r.scale.size();
    copy_1d(r.scale, dn.scale, "g3l:rmsnorm_scale");
  }

  // ---- Readout energy MLP (no bias) ----
  {
    const auto &mlp = m.energy_mlp;
    energy_n_layers = mlp.n_layers;
    energy_activation = m.energy_mlp_activation;
    int max_in = 0, max_out = 0;
    for (int i = 0; i < energy_n_layers; i++) { max_in = std::max(max_in, mlp.n_in[i]); max_out = std::max(max_out, mlp.n_out[i]); }
    energy_max_dim = std::max(max_in, max_out);
    MemKK::realloc_kokkos(d_energy_W, "g3l:energy_W", std::max(energy_n_layers, 1), std::max(max_in, 1), std::max(max_out, 1));
    auto hW = Kokkos::create_mirror_view(d_energy_W);
    Kokkos::deep_copy(hW, (NNScalar) 0);
    for (int l = 0; l < energy_n_layers; l++)
      for (int i = 0; i < mlp.n_in[l]; i++)
        for (int j = 0; j < mlp.n_out[l]; j++)
          hW(l, i, j) = (NNScalar) mlp.W[l][(size_t) i * mlp.n_out[l] + j];
    Kokkos::deep_copy(d_energy_W, hW);
    MemKK::realloc_kokkos(d_energy_norms, "g3l:energy_norms", std::max(energy_n_layers, 1));
    auto hn = Kokkos::create_mirror_view(d_energy_norms);
    for (int i = 0; i < energy_n_layers; i++) hn(i) = (NNScalar) mlp.norm[i];
    Kokkos::deep_copy(d_energy_norms, hn);
    MemKK::realloc_kokkos(d_energy_dims, "g3l:energy_dims", energy_n_layers + 1);
    auto hd = Kokkos::create_mirror_view(d_energy_dims);
    hd(0) = mlp.n_in.empty() ? 0 : mlp.n_in[0];
    for (int i = 0; i < energy_n_layers; i++) hd(i + 1) = mlp.n_out[i];
    Kokkos::deep_copy(d_energy_dims, hd);
  }

  // ---- UQ / extrapolation-grade artifacts (optional; schema v6 basis-RP) ----
  // Stage the GMM arrays + R onto the device; set the scalar members. The
  // per-block basis-width validation (needs idx_reduce_rho{1,2,3}) happens in
  // init_style() after those reduce indices are resolved.
  if (m.has_uq) {
    const int E = m.uq_n_elements;
    const int K = m.uq_max_clusters;
    const int D = m.uq_feature_dim;        // = rp_dim + n_density
    const int RP = m.uq_rp_dim;            // projection width (R cols)
    if (E != m.n_elements)
      error->all(FLERR, "GRACE-3L/KK UQ: uq_n_elements ({}) != model n_elements ({})",
                 E, m.n_elements);
    if (D > GRACE3L_UQ_MAX_RP_DIM)
      error->all(FLERR, "GRACE-3L/KK UQ: feature_dim ({}) exceeds GRACE3L_UQ_MAX_RP_DIM ({})",
                 D, GRACE3L_UQ_MAX_RP_DIM);
    has_uq = true;
    uq_Kmax = K;
    uq_D = D;
    uq_rp_dim = RP;
    uq_n_density = m.uq_n_density;
    uq_normalize = m.uq_normalize;
    uq_density_scale = m.uq_density_scale;
    uq_d_basis = m.uq_d_basis;
    copy_3d(m.uq_centroids, d_uq_centroids, E, K, D, "g3l:uq_centroids");
    copy_4d(m.uq_inv_cov, d_uq_inv_cov, E, K, D, D, "g3l:uq_inv_cov");
    copy_1d_int(m.uq_n_clusters, d_uq_n_clusters, "g3l:uq_n_clusters");
    copy_2d(m.uq_interp_thresholds, d_uq_interp_thresholds, E, K, "g3l:uq_interp_thr");
    copy_2d(m.uq_rp_matrix, d_uq_rp_matrix, m.uq_d_basis, RP, "g3l:uq_rp_matrix");
  }
}

// ======================================================================
// grow: allocate per-chunk geometry arrays (Task 2.3). Later Phase 2/3 tasks
// extend this with per-atom / per-bond NN intermediate arrays as they add
// the kernels that need them.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::grow(int natom, int maxneigh_in)
{
  if ((int) d_ncount.extent(0) < natom) {
    MemKK::realloc_kokkos(d_ncount, "g3l:ncount", natom);
    MemKK::realloc_kokkos(d_mu_i, "g3l:mu_i", natom);
    MemKK::realloc_kokkos(d_A1, "g3l:A1", natom, L1_nradmax, (lmax + 1) * (lmax + 1));
    // cp_l product output A1_2 [atoms, rank, nfunc] + shared projection scratch.
    MemKK::realloc_kokkos(d_A1_2, "g3l:A1_2", natom,
                          std::max(d_prod[0].rank, 1), std::max(d_prod[0].nfunc, 1));
    MemKK::realloc_kokkos(d_cpl_lproj, "g3l:cpl_lproj", natom,
                          std::max(cpl_scratch_rank, 1), std::max(cpl_scratch_nlm, 1));
    MemKK::realloc_kokkos(d_cpl_rproj, "g3l:cpl_rproj", natom,
                          std::max(cpl_scratch_rank, 1), std::max(cpl_scratch_nlm, 1));
    // FunctionReduceN / FCRight2Left outputs (Task 2.6).
    MemKK::realloc_kokkos(d_A1_2_red, "g3l:A1_2_red", natom,
                          std::max(d_reduce[idx_reduce_A1_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A1_2a, "g3l:A1_2a", natom,
                          std::max(d_fc[idx_fc_A1_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A1_2a].n_funcs_left, 1));
    // Remaining L1 chain (Task 2.7).
    MemKK::realloc_kokkos(d_A1_3, "g3l:A1_3", natom,
                          std::max(d_prod[idx_prod_A1_3].rank, 1),
                          std::max(d_prod[idx_prod_A1_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A1_3_red, "g3l:A1_3_red", natom,
                          std::max(d_reduce[idx_reduce_A1_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A1_2b, "g3l:A1_2b", natom,
                          std::max(d_fc[idx_fc_A1_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A1_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A1_4, "g3l:A1_4", natom,
                          std::max(d_prod[idx_prod_A1_4].rank, 1),
                          std::max(d_prod[idx_prod_A1_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A1_4_red, "g3l:A1_4_red", natom,
                          std::max(d_reduce[idx_reduce_A1_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_eq1, "g3l:eq1", natom,
                          std::max(d_reduce[idx_reduce_eq1].n_out, 1),
                          std::max(d_reduce[idx_reduce_eq1].n_funcs, 1));
    MemKK::realloc_kokkos(d_eq1_norm, "g3l:eq1_norm", natom,
                          std::max(d_eqnorm[idx_eqnorm_eq1_norm].n_out, 1),
                          std::max(d_eqnorm[idx_eqnorm_eq1_norm].M, 1));
    MemKK::realloc_kokkos(d_rho1, "g3l:rho1", natom,
                          std::max(d_reduce[idx_reduce_rho1].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho1].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho1_norm, "g3l:rho1_norm", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho1_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho1].n_funcs, 1));
    // Task 2.8: L2 equiv-ind SPBF basis A2 [atoms, A2_n_rad_max, A2_nfunc].
    MemKK::realloc_kokkos(d_A2, "g3l:A2", natom,
                          std::max(d_spbf[idx_spbf_A2].n_rad_max, 1),
                          std::max(d_spbf[idx_spbf_A2].nfunc, 1));
    // Task 2.9: remaining L2 chain (identical structure to L1's Task 2.6/2.7
    // allocations, EXCEPT everything downstream of A2 here is built on
    // A2_red — the one extra reduce L2 has that L1 does not).
    MemKK::realloc_kokkos(d_A2_red, "g3l:A2_red", natom,
                          std::max(d_reduce[idx_reduce_A2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2, "g3l:A2_2", natom,
                          std::max(d_prod[idx_prod_A2_2].rank, 1),
                          std::max(d_prod[idx_prod_A2_2].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_2_red, "g3l:A2_2_red", natom,
                          std::max(d_reduce[idx_reduce_A2_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2a, "g3l:A2_2a", natom,
                          std::max(d_fc[idx_fc_A2_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A2_2a].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A2_3, "g3l:A2_3", natom,
                          std::max(d_prod[idx_prod_A2_3].rank, 1),
                          std::max(d_prod[idx_prod_A2_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_3_red, "g3l:A2_3_red", natom,
                          std::max(d_reduce[idx_reduce_A2_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2b, "g3l:A2_2b", natom,
                          std::max(d_fc[idx_fc_A2_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A2_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A2_4, "g3l:A2_4", natom,
                          std::max(d_prod[idx_prod_A2_4].rank, 1),
                          std::max(d_prod[idx_prod_A2_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_4_red, "g3l:A2_4_red", natom,
                          std::max(d_reduce[idx_reduce_A2_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_eq2, "g3l:eq2", natom,
                          std::max(d_reduce[idx_reduce_eq2].n_out, 1),
                          std::max(d_reduce[idx_reduce_eq2].n_funcs, 1));
    MemKK::realloc_kokkos(d_eq2_norm, "g3l:eq2_norm", natom,
                          std::max(d_eqnorm[idx_eqnorm_eq2_norm].n_out, 1),
                          std::max(d_eqnorm[idx_eqnorm_eq2_norm].M, 1));
    MemKK::realloc_kokkos(d_rho2, "g3l:rho2", natom,
                          std::max(d_reduce[idx_reduce_rho2].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho2].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho2_norm, "g3l:rho2_norm", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho2_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho2].n_funcs, 1));
    // Task 2.10: L3 equiv-ind SPBF basis A3 [atoms, A3_n_rad_max, A3_nfunc] +
    // the terminal L3 chain. Identical structure to L2's Task 2.9 allocations
    // (everything downstream of A3 is built on A3_red) EXCEPT there is NO eq3;
    // A3_3/A3_4 are Lmax=0 so their dims come from the resolved metadata.
    MemKK::realloc_kokkos(d_A3, "g3l:A3", natom,
                          std::max(d_spbf[idx_spbf_A3].n_rad_max, 1),
                          std::max(d_spbf[idx_spbf_A3].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_red, "g3l:A3_red", natom,
                          std::max(d_reduce[idx_reduce_A3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2, "g3l:A3_2", natom,
                          std::max(d_prod[idx_prod_A3_2].rank, 1),
                          std::max(d_prod[idx_prod_A3_2].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_2_red, "g3l:A3_2_red", natom,
                          std::max(d_reduce[idx_reduce_A3_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2a, "g3l:A3_2a", natom,
                          std::max(d_fc[idx_fc_A3_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A3_2a].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A3_3, "g3l:A3_3", natom,
                          std::max(d_prod[idx_prod_A3_3].rank, 1),
                          std::max(d_prod[idx_prod_A3_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_3_red, "g3l:A3_3_red", natom,
                          std::max(d_reduce[idx_reduce_A3_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2b, "g3l:A3_2b", natom,
                          std::max(d_fc[idx_fc_A3_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A3_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A3_4, "g3l:A3_4", natom,
                          std::max(d_prod[idx_prod_A3_4].rank, 1),
                          std::max(d_prod[idx_prod_A3_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_4_red, "g3l:A3_4_red", natom,
                          std::max(d_reduce[idx_reduce_A3_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho3, "g3l:rho3", natom,
                          std::max(d_reduce[idx_reduce_rho3].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho3].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho3_norm, "g3l:rho3_norm", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho3_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho3].n_funcs, 1));
    // Task 2.11: per-atom readout energy (fp64), chunk-local like rho3_norm.
    MemKK::realloc_kokkos(d_e_atom, "g3l:e_atom", natom);

    // ---- Task 3.1: backward adjoint buffers (each mirrors its forward tensor's
    // shape). rho{1,2}_(norm_)adj are produced for Tasks 3.2/3.3 reuse; only the
    // rho3 -> L3 -> A3 chain is exercised/gated here. ----
    // Task 3.5b: d_rho1_norm_adj/d_rho2_norm_adj are now natom ([nall]) and
    // allocated in grow_global(); here we allocate only their chunk-scratch
    // gather targets (fed to the shared TagRMSNormBwd per chunk in P5/P4).
    MemKK::realloc_kokkos(d_rho1_norm_adj_local, "g3l:rho1_norm_adj_local", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho1_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho1].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho2_norm_adj_local, "g3l:rho2_norm_adj_local", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho2_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho2].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho3_norm_adj, "g3l:rho3_norm_adj", natom,
                          std::max(d_rmsnorm[idx_rmsnorm_rho3_norm].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho3].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho1_adj, "g3l:rho1_adj", natom,
                          std::max(d_reduce[idx_reduce_rho1].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho1].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho2_adj, "g3l:rho2_adj", natom,
                          std::max(d_reduce[idx_reduce_rho2].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho2].n_funcs, 1));
    MemKK::realloc_kokkos(d_rho3_adj, "g3l:rho3_adj", natom,
                          std::max(d_reduce[idx_reduce_rho3].n_out, 1),
                          std::max(d_reduce[idx_reduce_rho3].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_red_adj, "g3l:A3_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2_adj, "g3l:A3_2_adj", natom,
                          std::max(d_prod[idx_prod_A3_2].rank, 1),
                          std::max(d_prod[idx_prod_A3_2].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_2_red_adj, "g3l:A3_2_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A3_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2a_adj, "g3l:A3_2a_adj", natom,
                          std::max(d_fc[idx_fc_A3_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A3_2a].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A3_3_adj, "g3l:A3_3_adj", natom,
                          std::max(d_prod[idx_prod_A3_3].rank, 1),
                          std::max(d_prod[idx_prod_A3_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_3_red_adj, "g3l:A3_3_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A3_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_2b_adj, "g3l:A3_2b_adj", natom,
                          std::max(d_fc[idx_fc_A3_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A3_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A3_4_adj, "g3l:A3_4_adj", natom,
                          std::max(d_prod[idx_prod_A3_4].rank, 1),
                          std::max(d_prod[idx_prod_A3_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A3_4_red_adj, "g3l:A3_4_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A3_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A3_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A3_adj, "g3l:A3_adj", natom,
                          std::max(d_spbf[idx_spbf_A3].n_rad_max, 1),
                          std::max(d_spbf[idx_spbf_A3].nfunc, 1));
    // Opt 4: shared adj_prod buffer. Sized to the MAX over the equiv layers
    // (A2/A3) of n_rad_max and P = nlm_y*n_lm_ind, so one buffer serves both.
    {
      const int nlm_y = (lmax + 1) * (lmax + 1);
      const int adjp_nrad = std::max({d_spbf[idx_spbf_A2].n_rad_max,
                                      d_spbf[idx_spbf_A3].n_rad_max, 1});
      const int adjp_P = std::max({nlm_y * d_spbf[idx_spbf_A2].n_lm_ind,
                                   nlm_y * d_spbf[idx_spbf_A3].n_lm_ind, 1});
      MemKK::realloc_kokkos(d_adj_prod, "g3l:adj_prod", natom, adjp_nrad, adjp_P);
    }
    // cp_l backward rank-projection scratch (same size as the forward scratch).
    MemKK::realloc_kokkos(d_cpl_lproj_adj, "g3l:cpl_lproj_adj", natom,
                          std::max(cpl_scratch_rank, 1), std::max(cpl_scratch_nlm, 1));
    MemKK::realloc_kokkos(d_cpl_rproj_adj, "g3l:cpl_rproj_adj", natom,
                          std::max(cpl_scratch_rank, 1), std::max(cpl_scratch_nlm, 1));

    // ---- Task 3.2: Layer-2 backward adjoint buffers (mirror the L3 d_A3_* set;
    // each matches its forward L2 tensor's shape). d_rho2_(norm_)adj already
    // allocated above. ----
    MemKK::realloc_kokkos(d_A2_red_adj, "g3l:A2_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2_adj, "g3l:A2_2_adj", natom,
                          std::max(d_prod[idx_prod_A2_2].rank, 1),
                          std::max(d_prod[idx_prod_A2_2].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_2_red_adj, "g3l:A2_2_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A2_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2a_adj, "g3l:A2_2a_adj", natom,
                          std::max(d_fc[idx_fc_A2_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A2_2a].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A2_3_adj, "g3l:A2_3_adj", natom,
                          std::max(d_prod[idx_prod_A2_3].rank, 1),
                          std::max(d_prod[idx_prod_A2_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_3_red_adj, "g3l:A2_3_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A2_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_2b_adj, "g3l:A2_2b_adj", natom,
                          std::max(d_fc[idx_fc_A2_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A2_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A2_4_adj, "g3l:A2_4_adj", natom,
                          std::max(d_prod[idx_prod_A2_4].rank, 1),
                          std::max(d_prod[idx_prod_A2_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A2_4_red_adj, "g3l:A2_4_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A2_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A2_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A2_adj, "g3l:A2_adj", natom,
                          std::max(d_spbf[idx_spbf_A2].n_rad_max, 1),
                          std::max(d_spbf[idx_spbf_A2].nfunc, 1));
    MemKK::realloc_kokkos(d_eq2_adj, "g3l:eq2_adj", natom,
                          std::max(d_reduce[idx_reduce_eq2].n_out, 1),
                          std::max(d_reduce[idx_reduce_eq2].n_funcs, 1));
    // Chunk-local gather of the reverse-comm-completed message adjoint.
    MemKK::realloc_kokkos(d_eq2_norm_adj_local, "g3l:eq2_norm_adj_local", natom,
                          std::max(I2_n_funcs, 1), std::max(I2_n_out, 1));

    // ---- Task 3.3: Layer-1 backward adjoint buffers (mirror the L1 forward
    // tensor shapes). d_A1_adj is the direct fan-out target (raw A1). ----
    MemKK::realloc_kokkos(d_A1_2_adj, "g3l:A1_2_adj", natom,
                          std::max(d_prod[0].rank, 1), std::max(d_prod[0].nfunc, 1));
    MemKK::realloc_kokkos(d_A1_2_red_adj, "g3l:A1_2_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A1_2_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_2_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A1_2a_adj, "g3l:A1_2a_adj", natom,
                          std::max(d_fc[idx_fc_A1_2a].n_out, 1),
                          std::max(d_fc[idx_fc_A1_2a].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A1_3_adj, "g3l:A1_3_adj", natom,
                          std::max(d_prod[idx_prod_A1_3].rank, 1),
                          std::max(d_prod[idx_prod_A1_3].nfunc, 1));
    MemKK::realloc_kokkos(d_A1_3_red_adj, "g3l:A1_3_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A1_3_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_3_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A1_2b_adj, "g3l:A1_2b_adj", natom,
                          std::max(d_fc[idx_fc_A1_2b].n_out, 1),
                          std::max(d_fc[idx_fc_A1_2b].n_funcs_left, 1));
    MemKK::realloc_kokkos(d_A1_4_adj, "g3l:A1_4_adj", natom,
                          std::max(d_prod[idx_prod_A1_4].rank, 1),
                          std::max(d_prod[idx_prod_A1_4].nfunc, 1));
    MemKK::realloc_kokkos(d_A1_4_red_adj, "g3l:A1_4_red_adj", natom,
                          std::max(d_reduce[idx_reduce_A1_4_red].n_out, 1),
                          std::max(d_reduce[idx_reduce_A1_4_red].n_funcs, 1));
    MemKK::realloc_kokkos(d_A1_adj, "g3l:A1_adj", natom,
                          L1_nradmax, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(d_eq1_adj, "g3l:eq1_adj", natom,
                          std::max(d_reduce[idx_reduce_eq1].n_out, 1),
                          std::max(d_reduce[idx_reduce_eq1].n_funcs, 1));
    // d_rho1_adj already allocated with the rho{1,2,3}_adj set (Task 3.1).
    MemKK::realloc_kokkos(d_eq1_norm_adj_local, "g3l:eq1_norm_adj_local", natom,
                          std::max(I_n_funcs, 1), std::max(I_n_out, 1));
  }

  if ((int) d_radial_basis.extent(0) < natom || (int) d_radial_basis.extent(1) < maxneigh_in) {
    MemKK::realloc_kokkos(d_nearest, "g3l:nearest", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_rnorms, "g3l:rnorms", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_rhats, "g3l:rhats", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_mu_j, "g3l:mu_j", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_radial_basis, "g3l:radial_basis", natom, maxneigh_in, nradbase);
    MemKK::realloc_kokkos(d_env, "g3l:env", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_Y_bond, "g3l:Y_bond", natom, maxneigh_in, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(d_R1_nl, "g3l:R1_nl", natom, maxneigh_in, L1_nradmax, lmax + 1);
    // Task 3.3: radial r-derivative scratch (same shape as d_R1_nl, reused per
    // layer), the per-bond envelope derivative, and the per-bond force buffer.
    MemKK::realloc_kokkos(d_DR1_nl, "g3l:DR1_nl", natom, maxneigh_in, L1_nradmax, lmax + 1);
    MemKK::realloc_kokkos(d_denv, "g3l:denv", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_f_ij, "g3l:f_ij", natom, maxneigh_in);
#ifdef KOKKOS_ENABLE_CUDA
    // Opt 15: batched cuBLAS radial-MLP scratch (value path). M = natom*maxneigh
    // bonds; hidden widths = max n_out of layers 0/1 across all radial SPBFs.
    {
      const int Mrows = natom * maxneigh_in;
      int nin0 = 1, mh0 = 1, mh1 = 1, mout = 1;
      for (int si = 0; si < n_spbf; si++) {
        const auto &mlp = grace_model->spbf[si].mlp;
        const int nl = mlp.n_layers;
        if (nl >= 1) { nin0 = std::max(nin0, mlp.n_in[0]); mh0 = std::max(mh0, mlp.n_out[0]); }
        if (nl >= 2) { mh1 = std::max(mh1, mlp.n_out[1]); }
        if (nl >= 1) { mout = std::max(mout, mlp.n_out[nl - 1]); }
      }
      MemKK::realloc_kokkos(d_mlp_X,  "g3l:mlp_X",  std::max(Mrows, 1), nin0);
      MemKK::realloc_kokkos(d_mlp_h0, "g3l:mlp_h0", std::max(Mrows, 1), mh0);
      MemKK::realloc_kokkos(d_mlp_h1, "g3l:mlp_h1", std::max(Mrows, 1), mh1);
      MemKK::realloc_kokkos(d_mlp_R,  "g3l:mlp_R",  std::max(Mrows, 1), mout);
      // Opt 15b deriv-path scratch: dX is Chebyshev-only (nradbase cols).
      MemKK::realloc_kokkos(d_mlp_dX,  "g3l:mlp_dX",  std::max(Mrows, 1), std::max(nradbase, 1));
      MemKK::realloc_kokkos(d_mlp_dh0, "g3l:mlp_dh0", std::max(Mrows, 1), mh0);
      MemKK::realloc_kokkos(d_mlp_dh1, "g3l:mlp_dh1", std::max(Mrows, 1), mh1);
      MemKK::realloc_kokkos(d_mlp_dR,  "g3l:mlp_dR",  std::max(Mrows, 1), mout);
    }
#endif
  }
}

// ======================================================================
// grow_global (Task 2.8): allocate the global (nall-sized) inter-layer comm
// arrays for the Layer-2 indicator (eq1_norm) + its host mirror. Mirrors
// GRACE-2L's grow_global; d_grad_I_global is sized too so the (unused-in-this-
// task) reverse-comm path stays well-formed. Sized to I_n_funcs/I_n_out set
// in init_style.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::grow_global(int nall)
{
  if ((int) d_I_global.extent(0) < nall) {
    MemKK::realloc_kokkos(d_I_global, "g3l:I_global", nall, I_n_funcs, I_n_out);
    MemKK::realloc_kokkos(d_grad_I_global, "g3l:grad_I_global", nall, I_n_funcs, I_n_out);
    h_I_global = Kokkos::create_mirror_view(d_I_global);
    h_grad_I_global = Kokkos::create_mirror_view(d_grad_I_global);
  }
  // Round-2 (L3) indicator buffer (Task 2.10): eq2_norm is parity-doubled
  // (I2_n_out=50 vs I_n_out=25), so it needs its own [nall, I2_n_funcs, I2_n_out]
  // array — d_I_global cannot hold it. Separate extent check (distinct dims).
  if ((int) d_eq2_norm_global.extent(0) < nall) {
    MemKK::realloc_kokkos(d_eq2_norm_global, "g3l:eq2_norm_global", nall, I2_n_funcs, I2_n_out);
    h_eq2_norm_global = Kokkos::create_mirror_view(d_eq2_norm_global);
  }
  // Task 3.1: L3 message adjoint d_eq2_norm_adj — same [nall, I2_n_funcs, I2_n_out]
  // sizing as d_eq2_norm_global (scattered by atomic_add to owned+ghost neighbors).
  if ((int) d_eq2_norm_adj.extent(0) < nall) {
    MemKK::realloc_kokkos(d_eq2_norm_adj, "g3l:eq2_norm_adj", nall, I2_n_funcs, I2_n_out);
    // Task 3.2: host mirror for the reverse-comm host fallback path.
    h_eq2_norm_adj = Kokkos::create_mirror_view(d_eq2_norm_adj);
  }
  // Task 3.5a: natom promotions of d_rho1_norm/d_rho2_norm (see .h comment).
  // grow() above always runs before grow_global() in compute(), so
  // d_rho1_norm/d_rho2_norm already carry their correct (natom-independent)
  // trailing dims by this point — reuse them rather than re-deriving from the
  // model metadata indices.
  if ((int) d_rho1_norm_full.extent(0) < nall) {
    MemKK::realloc_kokkos(d_rho1_norm_full, "g3l:rho1_norm_full", nall,
                          (int) d_rho1_norm.extent(1), (int) d_rho1_norm.extent(2));
  }
  if ((int) d_rho2_norm_full.extent(0) < nall) {
    MemKK::realloc_kokkos(d_rho2_norm_full, "g3l:rho2_norm_full", nall,
                          (int) d_rho2_norm.extent(1), (int) d_rho2_norm.extent(2));
  }
  // Task 3.5b: natom promotions of the rho1/rho2 readout-adjoint SEEDS. Written
  // by TagReadoutBwd at (ii+chunk_offset) across ALL P3 chunks (every owned row
  // fully ASSIGNED -> no pre-zero), then gathered natom->_local per chunk in
  // P5/P4. Trailing dims mirror d_rho{1,2}_norm (== d_rho3_norm_adj's dims, which
  // TagReadoutBwd uses to bound the write loop). Sized to nall (safe superset).
  if ((int) d_rho1_norm_adj.extent(0) < nall) {
    MemKK::realloc_kokkos(d_rho1_norm_adj, "g3l:rho1_norm_adj", nall,
                          (int) d_rho1_norm.extent(1), (int) d_rho1_norm.extent(2));
  }
  if ((int) d_rho2_norm_adj.extent(0) < nall) {
    MemKK::realloc_kokkos(d_rho2_norm_adj, "g3l:rho2_norm_adj", nall,
                          (int) d_rho2_norm.extent(1), (int) d_rho2_norm.extent(2));
  }
}

// ======================================================================
// compute_Y_bond_chunk: fill d_Y_bond from d_rhats for the current chunk.
// COPIED VERBATIM (module renamed) from the validated GRACE-2L KOKKOS style —
// same alm/blm/cl/dl/idx_sph recurrence, same real-SH sign/ordering
// convention. Must be (re)called whenever d_rhats for the current chunk is
// refreshed, before any kernel reads d_Y_bond.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_Y_bond_chunk()
{
  auto Y_bond = d_Y_bond;
  auto rhats = d_rhats;
  auto nc_v = d_ncount;
  auto idx_s = d_idx_sph;
  auto alm_v = alm; auto blm_v = blm; auto cl_v = cl; auto dl_v = dl;
  int lmax_v = lmax;
  GeomScalar sq2_v = sq2;
  int cs = chunk_size;
  int mn = maxneigh;

  Kokkos::parallel_for("ComputeYBond",
    Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, mn}),
    KOKKOS_LAMBDA(const int ii, const int jj) {
      if (jj >= nc_v(ii)) return;
      const GeomScalar rx = rhats(ii, jj, 0);
      const GeomScalar ry = rhats(ii, jj, 1);
      const GeomScalar rz = rhats(ii, jj, 2);

      GeomScalar plm_arr[15];
      plm_arr[0] = GeomScalar(1.0);
      if (lmax_v > 0) {
        plm_arr[1] = GeomScalar(1.7320508075688772935) * rz;
        plm_arr[2] = GeomScalar(-1.2247448713915890491);
        for (int ll = 2; ll <= lmax_v; ll++) {
          for (int mm = 0; mm < ll - 1; mm++) {
            const int idx = ll*(ll+1)/2 + mm;
            const int i1 = (ll-1)*ll/2 + mm, i2 = (ll-2)*(ll-1)/2 + mm;
            const int ai = idx_s(ll*(ll+1) + mm);
            plm_arr[idx] = alm_v(ai) * (rz * plm_arr[i1] + blm_v(ai) * plm_arr[i2]);
          }
          { int idx = ll*(ll+1)/2+ll-1;
            plm_arr[idx] = dl_v(ll) * plm_arr[(ll-1)*ll/2+ll-1] * rz; }
          { int idx = ll*(ll+1)/2+ll;
            plm_arr[idx] = cl_v(ll) * plm_arr[(ll-1)*ll/2+ll-1]; }
        }
      }

      const GeomScalar phase_re = rx, phase_im = ry;
      for (int l = 0; l <= lmax_v; l++) {
        Y_bond(ii, jj, l*(l+1)) = plm_arr[l*(l+1)/2];
        GeomScalar pm_re = phase_re, pm_im = phase_im;
        for (int m = 1; m <= l; m++) {
          const int fac = (m % 2 == 0) ? 1 : -1;
          const GeomScalar pv = plm_arr[l*(l+1)/2 + m];
          Y_bond(ii, jj, l*(l+1)+m) = sq2_v * fac * pm_re * pv;
          Y_bond(ii, jj, l*(l+1)-m) = sq2_v * fac * pm_im * pv;
          const GeomScalar t_re = pm_re * phase_re - pm_im * phase_im;
          const GeomScalar t_im = pm_re * phase_im + pm_im * phase_re;
          pm_re = t_re; pm_im = t_im;
        }
      }
    });
}

// ======================================================================
// Kernel: ComputeNeigh — half-list -> per-atom filtered neighbor arrays.
// COPIED VERBATIM (module renamed) from the validated GRACE-2L KOKKOS style.
// d_rhats stores the unit vector (x_j - x_i)/r — i.e. delx = x_i - x_j
// negated — matching the oracle's bond_vector = r_j - r_i convention.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeNeigh,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type(i);
  const GeomScalar xtmp = x(i,0);
  const GeomScalar ytmp = x(i,1);
  const GeomScalar ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * maxneigh;
  int* inside = (int*)team.team_shmem().get_shmem(team.team_size() * maxneigh * sizeof(int), 0) + scratch_shift;

  int ncount = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, jnum),
      [&] (const int jj, int& count) {
    int j = d_neighbors(i, jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const GeomScalar delx = xtmp - x(j,0);
    const GeomScalar dely = ytmp - x(j,1);
    const GeomScalar delz = ztmp - x(j,2);
    const GeomScalar rsq = delx*delx + dely*dely + delz*delz;
    inside[jj] = -1;
    if (rsq < d_cutsq(itype, jtype)) {
      inside[jj] = 1;
      count++;
    }
  }, ncount);

  d_ncount(ii) = ncount;
  d_mu_i(ii) = d_map(itype);

  Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, jnum),
      [&] (const int jj, int& offset, bool final) {
    if (inside[jj] < 0) return;
    if (final) {
      int j = d_neighbors(i, jj);
      j &= NEIGHMASK;
      const GeomScalar delx = xtmp - x(j,0);
      const GeomScalar dely = ytmp - x(j,1);
      const GeomScalar delz = ztmp - x(j,2);
      const GeomScalar rsq = delx*delx + dely*dely + delz*delz;
      // Oracle's BondLength.frwrd regularizes: r = sqrt(sum(bond_vector^2) + 1e-10).
      // Deviation from the 2L port (which omits this) — required to hit the
      // 1e-12 gate tolerance; 2L never needed it at its looser tolerance.
      const GeomScalar r = Kokkos::sqrt(rsq + GeomScalar(1e-10));
      const GeomScalar rinv = GeomScalar(1.0) / r;
      d_mu_j(ii, offset) = d_map(type(j));
      d_rnorms(ii, offset) = r;
      d_rhats(ii, offset, 0) = -delx * rinv;
      d_rhats(ii, offset, 1) = -dely * rinv;
      d_rhats(ii, offset, 2) = -delz * rinv;
      d_nearest(ii, offset) = j;
    }
    offset++;
  });
}

// ======================================================================
// Kernel: ComputeRadialBasis — Chebyshev (kind=1) basis T_1..T_nradbase,
// RAW (no envelope folded in) + a separate polynomial cutoff envelope d_env.
// Adapted from the validated GRACE-2L KOKKOS style: 2L multiplies the
// envelope into d_radial_basis and also tracks its r-derivative; the 3L
// oracle (numpy_forward_3l.chebyshev_basis/envelope) keeps them SEPARATE and
// this task does not yet need the derivative (added by the force-kernel
// task), so both are dropped here.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeRadialBasis,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (jj >= d_ncount(ii)) return;

  const GeomScalar r = d_rnorms(ii, jj);
  const GeomScalar rcut_ij = d_bond_cutoff(d_mu_i(ii), d_mu_j(ii, jj));
  const GeomScalar x_norm = r / rcut_ij;

  const int p = radial_basis_p;
  GeomScalar xp = GeomScalar(1.0);
  for (int ip = 0; ip < p; ip++) xp *= x_norm;
  const GeomScalar xp1 = xp * x_norm, xp2 = xp1 * x_norm;
  const GeomScalar pp1 = GeomScalar(p) * (p + 1);
  const GeomScalar pp2 = GeomScalar(p) * (p + 2);
  const GeomScalar p1p2 = GeomScalar(p + 1) * (p + 2);
  const GeomScalar fcut = GeomScalar(1.0) - GeomScalar(0.5) * p1p2 * xp + pp2 * xp1 - GeomScalar(0.5) * pp1 * xp2;
  d_env(ii, jj) = fcut;

  const GeomScalar x_cheb = GeomScalar(2.0) * x_norm - GeomScalar(1.0);

  GeomScalar T_prev = GeomScalar(1.0), T_curr = x_cheb;
  d_radial_basis(ii, jj, 0) = T_curr;

  for (int k = 1; k < nradbase; k++) {
    const GeomScalar T_next = GeomScalar(2.0) * x_cheb * T_curr - T_prev;
    T_prev = T_curr;
    T_curr = T_next;
    d_radial_basis(ii, jj, k) = T_curr;
  }
}

// ======================================================================
// Kernel: TagComputeMLPRadial (Task 2.4 A1 / Task 2.8 A2) — generic radial
// MLP, one bond per thread, for the SPBF layer selected by the member
// radial_mlp_spbf (0=A1, idx_spbf_A2=A2). Reproduces
// numpy_forward_3l.radial_mlp EXACTLY:
//   x = concat([basis32, z_table[mu_i], z_table[mu_j]])         (width 138)
//   x = x @ (W_i*norm_i); hidden layers (i<n_layers-1) add b_i*norm_i, silu;
//   output layer (i==n_layers-1) is linear, no bias.
// Output reshaped [n_rad_max, lmax+1] (row-major n*[lmax+1]+l) and stored
// RAW (pre l_tile-gather) in d_R1_nl; the A*-assembly kernel gathers by
// l_tile when it forms a_nl. All arithmetic in NNScalar (fp32 Mixed
// default); d_radial_basis (GeomScalar) is cast down on read, d_chem_embed
// is already NNScalar. Team-policy (ii,jj) decode copied from
// TagComputeRadialBasis (Task 2.3). The A1/A2 MLPs share this code path but
// use their own weights (d_spbf[radial_mlp_spbf].mlp_*) and n_rad_max; the
// chem-embedding concat input is spbf-independent.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeMLPRadial,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (jj >= d_ncount(ii)) return;

  const DeviceSPBF &sp = d_spbf[radial_mlp_spbf];   // A1 or A2 (set on host)
  const int n_layers = sp.mlp_n_layers;
  int mlp_dims[GRACE3L_MAX_MLP_LAYERS + 1];
  NNScalar mlp_norms[GRACE3L_MAX_MLP_LAYERS];
  for (int i = 0; i <= n_layers; i++) mlp_dims[i] = sp.mlp_dims(i);
  for (int i = 0; i < n_layers; i++) mlp_norms[i] = sp.mlp_norms(i);

  NNScalar buf_a[GRACE3L_MAX_MLP_DIM], buf_b[GRACE3L_MAX_MLP_DIM];
  NNScalar* h_cur = buf_a; NNScalar* h_nxt = buf_b;

  // Input: concat([basis32, z_i32, z_j32]) — order pinned by the oracle.
  const int mu_i = d_mu_i(ii);
  const int mu_j = d_mu_j(ii, jj);
  int k = 0;
  for (int b = 0; b < nradbase; b++) h_cur[k++] = (NNScalar) d_radial_basis(ii, jj, b);
  for (int e = 0; e < embedding_size; e++) h_cur[k++] = d_chem_embed(mu_i, e);
  for (int e = 0; e < embedding_size; e++) h_cur[k++] = d_chem_embed(mu_j, e);

  for (int layer = 0; layer < n_layers - 1; layer++) {
    const int nin = mlp_dims[layer], nout = mlp_dims[layer + 1];
    const NNScalar norm = mlp_norms[layer];
    for (int j = 0; j < nout; j++) {
      NNScalar s = 0.0;
      for (int kk = 0; kk < nin; kk++)
        s += sp.mlp_W(layer, kk, j) * h_cur[kk];
      s = s * norm + sp.mlp_b(layer, j) * norm;
      h_nxt[j] = silu(s);
    }
    NNScalar* tmp = h_cur; h_cur = h_nxt; h_nxt = tmp;
  }

  // Output layer: linear, no bias.
  const int last_layer = n_layers - 1;
  const int n_last_hidden = mlp_dims[last_layer];
  const NNScalar out_norm = mlp_norms[last_layer];
  const int lmax1 = lmax + 1;
  for (int n = 0; n < sp.n_rad_max; n++) {
    for (int l = 0; l < lmax1; l++) {
      const int j = n * lmax1 + l;
      NNScalar sum = 0.0;
      for (int kk = 0; kk < n_last_hidden; kk++)
        sum += sp.mlp_W(last_layer, kk, j) * h_cur[kk];
      d_R1_nl(ii, jj, n, l) = sum * out_norm;
    }
  }
}

// ======================================================================
// Kernel: ComputeA1 (Task 2.4) — Layer-1 scalar-indicator SPBF assembly,
// per (central-atom ii, radial-channel n). Reproduces
// numpy_forward_3l.spbf_scalar EXACTLY:
//   a_nl[lm] = R1_nl[n, l_tile[lm]] * Y_bond[lm] * z_tr[mu_j][n] * env
//   A1[ii,n,lm] = inv_avg_n_neigh * Sum_{j in neigh(ii)} a_nl[lm]
// z_tr is the PRECOMPUTED per-element species table (norm folded once at
// copy_weights_to_device); env/Y_bond (GeomScalar) are cast down to
// NNScalar on read. One thread owns the whole (ii,n,:) row -> no atomics,
// zero-init folded into the same kernel. Team-policy (ii,n) decode mirrors
// TagComputeRadialBasis's (ii,jj) decode, dividing by L1_nradmax instead of
// maxneigh.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeA1,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeA1>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int n = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (n >= L1_nradmax) return;

  const DeviceSPBF &sp = d_spbf[0];   // A1
  const int nlm = (lmax + 1) * (lmax + 1);
  NNScalar acc[25];
  for (int lm = 0; lm < nlm; lm++) acc[lm] = NNScalar(0.0);

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const NNScalar z_n = sp.z_tr(mu_j, n);
    const NNScalar env = (NNScalar) d_env(ii, jj);
    for (int lm = 0; lm < nlm; lm++) {
      const int l = sp.l_tile(lm);
      NNScalar a = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, lm);
      a *= z_n;
      a *= env;
      acc[lm] += a;
    }
  }
  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;
  for (int lm = 0; lm < nlm; lm++)
    d_A1(ii, n, lm) = acc[lm] * inv_avg;
}

// ======================================================================
// Kernel: TagComputeA2 (Task 2.8) — Layer-2 equivariant-indicator SPBF
// assembly, per (central-atom ii, radial-channel n). Reproduces
// numpy_forward_3l.spbf_equivariant EXACTLY for the A2 layer:
//   a_nl[lm_y] = R1_nl[n, l_tile[lm_y]] * Y_bond[lm_y] * env       (no z_tr)
//   bond_I[lm_ind] = indicator[ind_j][n, lm_ind]
//                    + z_proj[mu_j][n] * chem_l0_mask[lm_ind]      (L0 chem inject)
//   prod[lm_y, lm_ind] = inv_avg_n_neigh * Sum_j a_nl[lm_y] * bond_I[lm_ind]
//   A2[n, f] = Sum_p prod_flat[p] * cg_W[p, f]   (p = lm_y * n_lm_ind + lm_ind)
// indicator[ind_j] is read from the forward-comm'd global array d_I_global
// (layout [nall, n=I_n_funcs, lm_ind=I_n_out]); ind_j = d_nearest(ii,jj) is the
// LAMMPS local index of the neighbor (owned OR ghost) — ghost rows are filled
// by comm->forward_comm in compute() before this kernel runs.
// z_proj is the PRECOMPUTED per-element chem projection (chem_linear folded once
// at copy_weights_to_device). One thread owns the whole (ii,n,:) A2 row -> the
// per-atom 4-D product slice for this (ii,n) is a private [nlm_y, n_lm_ind]
// stack accumulator, so no atomics. Team-policy (ii,n) decode mirrors
// TagComputeA1, dividing by n_rad_max (A2's radial ceiling).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeA2,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeA2>::member_type& team) const
{
  const DeviceSPBF &sp = d_spbf[idx_spbf_A2];   // A2 (equivariant)
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int n = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (n >= sp.n_rad_max) return;

  const int nlm_y = (lmax + 1) * (lmax + 1);    // Y lm axis (<=25)
  const int n_lm_ind = sp.n_lm_ind;             // indicator lm axis (<=GRACE3L_MAX_SPBF_NLMIND)
  const int nfunc = sp.nfunc;                   // CG-couple output width
  const int P = nlm_y * n_lm_ind;               // flattened prod / cg_W rows

  // per-(ii,n) private 4-D product slice prod[lm_y, lm_ind], row-major.
  NNScalar prod_acc[25 * GRACE3L_MAX_SPBF_NLMIND];
  for (int p = 0; p < P; p++) prod_acc[p] = NNScalar(0.0);

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const int j_global = d_nearest(ii, jj);     // neighbor local index (owned OR ghost)
    const NNScalar env = (NNScalar) d_env(ii, jj);

    NNScalar a_nl[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int l = sp.l_tile(ly);
      a_nl[ly] = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, ly) * env;
    }

    NNScalar bond_I[GRACE3L_MAX_SPBF_NLMIND];
    const NNScalar zpn = sp.z_proj(mu_j, n);
    for (int lr = 0; lr < n_lm_ind; lr++)
      bond_I[lr] = d_I_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);

    for (int ly = 0; ly < nlm_y; ly++) {
      const NNScalar a = a_nl[ly];
      const int base = ly * n_lm_ind;
      for (int lr = 0; lr < n_lm_ind; lr++)
        prod_acc[base + lr] += a * bond_I[lr];
    }
  }

  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;
  for (int p = 0; p < P; p++) prod_acc[p] *= inv_avg;

  // dense CG couple: A2[n,f] = Sum_p prod_flat[p] * cg_W[p,f].
  for (int f = 0; f < nfunc; f++) {
    NNScalar s = NNScalar(0.0);
    for (int p = 0; p < P; p++)
      s += prod_acc[p] * sp.cg_W(p, f);
    d_A2(ii, n, f) = s;
  }
}

// ======================================================================
// Kernel: TagComputeA3 (Task 2.10) — Layer-3 equivariant-indicator SPBF
// assembly. IDENTICAL to TagComputeA2 except: (1) uses the A3 SPBF weights;
// (2) reads the indicator from the round-2 global buffer d_eq2_norm_global
// (n_lm_ind = A3.n_lm_ind = 50, parity-doubled vs A2's 25); (3) writes d_A3.
// The stack arrays are unchanged: n_lm_ind=50 <= GRACE3L_MAX_SPBF_NLMIND=64 and
// P = 25*50 = 1250 <= 25*GRACE3L_MAX_SPBF_NLMIND = 1600. Ghost rows of
// d_eq2_norm_global are filled by round-2 comm->forward_comm before this runs.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeA3,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeA3>::member_type& team) const
{
  const DeviceSPBF &sp = d_spbf[idx_spbf_A3];   // A3 (equivariant)
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int n = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (n >= sp.n_rad_max) return;

  const int nlm_y = (lmax + 1) * (lmax + 1);    // Y lm axis (<=25)
  const int n_lm_ind = sp.n_lm_ind;             // indicator lm axis (=50 for A3, <=GRACE3L_MAX_SPBF_NLMIND)
  const int nfunc = sp.nfunc;                   // CG-couple output width
  const int P = nlm_y * n_lm_ind;               // flattened prod / cg_W rows (=1250 for A3)

  // per-(ii,n) private 4-D product slice prod[lm_y, lm_ind], row-major.
  NNScalar prod_acc[25 * GRACE3L_MAX_SPBF_NLMIND];
  for (int p = 0; p < P; p++) prod_acc[p] = NNScalar(0.0);

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const int j_global = d_nearest(ii, jj);     // neighbor local index (owned OR ghost)
    const NNScalar env = (NNScalar) d_env(ii, jj);

    NNScalar a_nl[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int l = sp.l_tile(ly);
      a_nl[ly] = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, ly) * env;
    }

    NNScalar bond_I[GRACE3L_MAX_SPBF_NLMIND];
    const NNScalar zpn = sp.z_proj(mu_j, n);
    for (int lr = 0; lr < n_lm_ind; lr++)
      bond_I[lr] = d_eq2_norm_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);

    for (int ly = 0; ly < nlm_y; ly++) {
      const NNScalar a = a_nl[ly];
      const int base = ly * n_lm_ind;
      for (int lr = 0; lr < n_lm_ind; lr++)
        prod_acc[base + lr] += a * bond_I[lr];
    }
  }

  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;
  for (int p = 0; p < P; p++) prod_acc[p] *= inv_avg;

  // dense CG couple: A3[n,f] = Sum_p prod_flat[p] * cg_W[p,f].
  for (int f = 0; f < nfunc; f++) {
    NNScalar s = NNScalar(0.0);
    for (int p = 0; p < P; p++)
      s += prod_acc[p] * sp.cg_W(p, f);
    d_A3(ii, n, f) = s;
  }
}

// ======================================================================
// Kernel: TagComputeSPBFFused (Opt 6) — team+shared replacement for
// TagComputeA2/TagComputeA3. The old kernels gave one thread the whole (ii,n)
// row, so the [nlm_y, n_lm_ind] product slice (P=1250 for A3) lived in a
// per-thread stack array — REG:72 STACK:6768, i.e. a massive local-memory spill
// that made these the #1 hotspot (ComputeA3 ~25.6%). Here ONE TEAM owns one
// (ii,n): the product slice prod_sh[P] lives in TEAM SHARED, filled cooperatively
// over neighbours (a_nl/bond_I staged in shared per neighbour so global loads stay
// minimal — no atomics: each thread owns disjoint prod slots p), then the dense
// CG couple A[n,f]=Σ_p prod_sh[p]*cg_W[p,f] is split across the team (one f per
// thread). Parameterized by launch state so A2 (indicator d_I_global, out d_A2)
// and A3 (indicator d_eq2_norm_global, out d_A3) share one kernel.
// Shared layout: [ prod_sh(P) | anl_sh(nlm_y) | bondI_sh(n_lm_ind) ] NNScalars.
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeSPBFFused,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeSPBFFused>::member_type& team) const
{
  const DeviceSPBF &sp = d_spbf[spbf_fused_idx];
  const int ii = team.league_rank() % chunk_size;   // league = chunk_size * n_rad_max
  const int n  = team.league_rank() / chunk_size;    // -> ii<chunk_size, n<n_rad_max always

  const int nlm_y = (lmax + 1) * (lmax + 1);         // Y lm axis (<=25)
  const int n_lm_ind = sp.n_lm_ind;                  // indicator lm axis (25 for A2, 50 for A3)
  const int nfunc = sp.nfunc;                        // CG-couple output width
  const int P = nlm_y * n_lm_ind;                    // flattened prod / cg_W rows

  NNScalar* sh = (NNScalar*) team.team_shmem().get_shmem((P + nlm_y + n_lm_ind) * sizeof(NNScalar), 0);
  NNScalar* prod_sh  = sh;
  NNScalar* anl_sh   = sh + P;
  NNScalar* bondI_sh = sh + P + nlm_y;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, P), [&] (const int p) {
    prod_sh[p] = NNScalar(0.0);
  });
  team.team_barrier();

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const int j_global = d_nearest(ii, jj);          // neighbor local index (owned OR ghost)
    const NNScalar env = (NNScalar) d_env(ii, jj);
    const NNScalar zpn = sp.z_proj(mu_j, n);
    // Stage this neighbour's a_nl[ly] and bond_I[lr] into shared (disjoint regions).
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlm_y), [&] (const int ly) {
      const int l = sp.l_tile(ly);
      anl_sh[ly] = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, ly) * env;
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n_lm_ind), [&] (const int lr) {
      bondI_sh[lr] = spbf_fused_indicator(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);
    });
    team.team_barrier();
    // Outer-product accumulate into the shared prod slice; slot p=(ly,lr) is owned
    // by exactly one thread across the loop, so no atomics are needed.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, P), [&] (const int p) {
      const int ly = p / n_lm_ind;
      const int lr = p - ly * n_lm_ind;
      prod_sh[p] += anl_sh[ly] * bondI_sh[lr];
    });
    team.team_barrier();
  }

  // Sparse CG couple (Opt 9), one output function f per thread. cg_W is ~99%
  // zeros, so iterate only column f's nonzeros (CSC, ascending p = the dense
  // loop's accumulation order — numerically identical). inv_avg is folded into
  // the per-f result (A[n,f] = inv_avg·Σ_p prod·cg), sparing a whole
  // prod-scaling pass + barrier. All prod_sh writes are complete (barrier at
  // the end of the last neighbour iteration), so no extra barrier before this read.
  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nfunc), [&] (const int f) {
    NNScalar s = NNScalar(0.0);
    const int e1 = sp.cg_col_ptr(f + 1);
    for (int e = sp.cg_col_ptr(f); e < e1; e++)
      s += prod_sh[sp.cg_row_idx(e)] * sp.cg_col_val(e);
    spbf_fused_out(ii, n, f) = s * inv_avg;
  });
}

// ======================================================================
// compute_spbf_equiv (Opt 6) — launcher for the fused equiv-SPBF (A2/A3).
// Sets the per-call launch state captured by TagComputeSPBFFused, then launches
// one team per (ii,n) over the chunk with the prod slice + staging buffers in
// team shared. Replaces the TagComputeA2/TagComputeA3 launches.
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_spbf_equiv(
    int idx, const t_nn_3d &indicator, const t_nn_3d &out)
{
  spbf_fused_idx = idx;
  spbf_fused_indicator = indicator;
  spbf_fused_out = out;

  const int nradmax = d_spbf[idx].n_rad_max;
  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int n_lm_ind = d_spbf[idx].n_lm_ind;
  const int P = nlm_y * n_lm_ind;

  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 128;
  const int league = chunk_size * nradmax;
  int scratch_size = scratch_size_helper<NNScalar>(P + nlm_y + n_lm_ind);
  check_team_size_for<TagComputeSPBFFused>(league, team_size, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeSPBFFused>(league, team_size, vector_length);
  policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
  Kokkos::parallel_for("ComputeSPBFFused", policy, *this);
}

// ======================================================================
// Kernel: TagCPLProject (Task 2.5) — cp_l Pass A. Projects each input onto the
// rank axis per lm, reproducing the oracle einsum("wrn,anw->arw", U[gl], left)
// * norm_u  (w=lm, n=feature, r=rank):
//   lproj[a,r,w] = norm_u * Sum_n U[group_left[w], r, n] * left[a, n, w]
//   rproj[a,r,w] = norm_v * Sum_n V[group_right[w], r, n] * right[a, n, w]
// Flat index over (atom a, rank r, lm w); one thread owns one (a,r,w) triple,
// so the scratch write is exclusive (no atomics). n_lm_left may differ from
// n_lm_right, so each projection is guarded by its own group_*.extent(0).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLProject,
    const int& idx) const
{
  const int rank = cpl_p.rank;
  const int nlm = cpl_nlm;
  const int a = idx / (rank * nlm);
  const int rem = idx - a * (rank * nlm);
  const int r = rem / nlm;
  const int w = rem - r * nlm;

  const int n_lm_left = (int) cpl_p.group_left.extent(0);
  if (w < n_lm_left) {
    const int gl = cpl_p.group_left(w);
    NNScalar s = NNScalar(0.0);
    const int nL = cpl_p.n_left;
    for (int nn = 0; nn < nL; nn++)
      s += cpl_p.U(gl, r, nn) * cpl_in_left(a, nn, w);
    d_cpl_lproj(a, r, w) = s * cpl_p.norm_u;
  }

  const int n_lm_right = (int) cpl_p.group_right.extent(0);
  if (w < n_lm_right) {
    const int gr = cpl_p.group_right(w);
    NNScalar s = NNScalar(0.0);
    const int nR = cpl_p.n_right;
    for (int nn = 0; nn < nR; nn++)
      s += cpl_p.V(gr, r, nn) * cpl_in_right(a, nn, w);
    d_cpl_rproj(a, r, w) = s * cpl_p.norm_v;
  }
}

// ======================================================================
// Kernel: TagCPLCouple (Task 2.5) — cp_l Pass B. Bilinear Clebsch-Gordan
// coupling on the rank axis + segment-sum over m_sum_ind:
//   out[a,r,:] = 0
//   for c: out[a,r,m_sum_ind[c]] += lproj[a,r,left_ind[c]]
//                                   * rproj[a,r,right_ind[c]] * cg[c]
// Flat index over (atom a, rank r); one thread owns the whole out[a,r,:] row,
// so the segment-sum accumulation is serial within a thread -> no atomics.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLCouple,
    const int& idx) const
{
  const int rank = cpl_p.rank;
  // t_nn_3d is LayoutLeft on CUDA, so atom is the unit-stride dimension.
  // Map adjacent threads to adjacent atoms for coalesced lproj/rproj reads and
  // output writes. Each (a,r) still has one owner and scans CG entries in the
  // same order, so arithmetic and determinism are unchanged.
  const int a = idx % chunk_size;
  const int r = idx / chunk_size;

  const int nfunc = cpl_p.nfunc;
  for (int fcn = 0; fcn < nfunc; fcn++)
    cpl_out(a, r, fcn) = NNScalar(0.0);

  const int n_cg = cpl_p.n_cg;
  for (int c = 0; c < n_cg; c++) {
    const NNScalar v = d_cpl_lproj(a, r, cpl_p.left_ind(c)) *
                       d_cpl_rproj(a, r, cpl_p.right_ind(c)) *
                       cpl_p.cg(c);
    cpl_out(a, r, cpl_p.m_sum_ind(c)) += v;
  }
}

// ======================================================================
// Kernel: TagCPLFused (Opt 5) — fused forward cp_l. ONE TEAM per (atom a, rank r).
// The two-pass Project+Couple, previously two global kernels round-tripping lproj/rproj
// [chunk,rank,nlm] through DRAM and doing an O(n_cg) global read-modify-write into the
// out[a,r,:] row, is fused: lproj/rproj and the out row live in TEAM SHARED memory. The
// couple's segment-sum uses SHARED atomics (fast) and out is written to global ONCE.
// Diagnosis (PROFILE.md): CPL kernels are memory-bandwidth-bound (REG:40/STACK:0, high
// occupancy), so cutting the lproj/rproj round-trip + the n_cg out-RMW is the lever.
// Shared layout: [ lproj_sh(nlm) | rproj_sh(nlm) | out_sh(nfunc) ] NNScalars.
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLFused,
    const typename Kokkos::TeamPolicy<DeviceType, TagCPLFused>::member_type& team) const
{
  const int rank = cpl_p.rank;
  const int a = team.league_rank() / rank;
  const int r = team.league_rank() - a * rank;
  const int nlm = cpl_nlm;
  const int nfunc = cpl_p.nfunc;
  const int n_lm_left = (int) cpl_p.group_left.extent(0);
  const int n_lm_right = (int) cpl_p.group_right.extent(0);
  const int nL = cpl_p.n_left;
  const int nR = cpl_p.n_right;

  NNScalar* sh = (NNScalar*) team.team_shmem().get_shmem((2 * nlm + nfunc) * sizeof(NNScalar), 0);
  NNScalar* lproj_sh = sh;
  NNScalar* rproj_sh = sh + nlm;
  NNScalar* out_sh   = sh + 2 * nlm;

  // Pass A (project) into shared: lproj_sh[w] = norm_u*Σ_n U[gl(w),r,n]*left[a,n,w].
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlm), [&] (const int w) {
    NNScalar sl = NNScalar(0.0);
    if (w < n_lm_left) {
      const int gl = cpl_p.group_left(w);
      for (int nn = 0; nn < nL; nn++) sl += cpl_p.U(gl, r, nn) * cpl_in_left(a, nn, w);
      sl *= cpl_p.norm_u;
    }
    lproj_sh[w] = sl;
    NNScalar sr = NNScalar(0.0);
    if (w < n_lm_right) {
      const int gr = cpl_p.group_right(w);
      for (int nn = 0; nn < nR; nn++) sr += cpl_p.V(gr, r, nn) * cpl_in_right(a, nn, w);
      sr *= cpl_p.norm_v;
    }
    rproj_sh[w] = sr;
  });
  // Zero the out row in shared.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nfunc), [&] (const int fcn) {
    out_sh[fcn] = NNScalar(0.0);
  });
  team.team_barrier();

  // Pass B (couple): out_sh[m_sum_ind[c]] += lproj_sh[left_ind[c]]*rproj_sh[right_ind[c]]*cg[c]
  // via shared atomics (segment-sum over the CG connections, split across the team).
  // Opt 14: the three indices come packed in one int32 (8 B/entry with cg vs 16 B).
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, cpl_p.n_cg), [&] (const int c) {
    const int pk = cpl_p.cg_packed(c);
    const NNScalar v = lproj_sh[pk & 255] * rproj_sh[(pk >> 8) & 255] * cpl_p.cg(c);
    Kokkos::atomic_add(&out_sh[pk >> 16], v);
  });
  team.team_barrier();

  // Write the out row to global ONCE.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nfunc), [&] (const int fcn) {
    cpl_out(a, r, fcn) = out_sh[fcn];
  });
}

// ======================================================================
// compute_cp_l (Task 2.5, Opt 5) — generic cp_l GeneralProductFunction product.
// Sets the per-call launch state (weights + in/out views) captured by the TagCPL*
// functors, then launches the fused team kernel (one team per (a,r)) over the chunk.
// Reused verbatim by every A*_* cp_l product (A1_2, A1_3, ... A3_2).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_cp_l(
    const DeviceCPL &p, const t_nn_3d &in_left, const t_nn_3d &in_right,
    const t_nn_3d &out)
{
  cpl_p = p;
  cpl_in_left = in_left;
  cpl_in_right = in_right;
  cpl_out = out;
  const int n_lm_left = (int) p.group_left.extent(0);
  const int n_lm_right = (int) p.group_right.extent(0);
  cpl_nlm = (n_lm_left > n_lm_right) ? n_lm_left : n_lm_right;

#ifdef KOKKOS_ENABLE_CUDA
  // Opt 18 (fp32 only): un-fuse. cuBLAS U/V rank-projection into the global lproj/rproj
  // scratch, then the couple reads them back. The projection is ~78% of the fused cp_l
  // cost and is the mirror of the Opt-16 backward projection GEMM (reuses p.Uw/p.Vw).
  // fp64 + non-CUDA keep the fused shared-memory kernel below (no SGEMM benefit, and
  // un-fusing reintroduces the lproj/rproj round-trip that Opt-5 removed).
  if (cublas_handle && sizeof(NNScalar) == 4) {
    compute_cp_l_project_cublas(p, in_left, in_right);
    Kokkos::parallel_for("CPLCouple",
        Kokkos::RangePolicy<DeviceType, TagCPLCouple>(0, chunk_size * p.rank), *this);
    return;
  }
#endif

  // Fused Project+Couple, one team per (a,r); lproj/rproj + out row in team shared.
  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 128;
  const int league = chunk_size * p.rank;
  int scratch_size = scratch_size_helper<NNScalar>(2 * cpl_nlm + p.nfunc);
  check_team_size_for<TagCPLFused>(league, team_size, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagCPLFused>(league, team_size, vector_length);
  policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
  Kokkos::parallel_for("CPLFused", policy, *this);
}

#ifdef KOKKOS_ENABLE_CUDA
// Opt 18: cuBLAS cp_l forward U/V rank-projection into the global lproj/rproj scratch.
//   lproj(a,r,w) = norm_u * Σ_n in_left(a,n,w) * U(gl(w),r,n)   (batch over w = lm)
//   rproj(a,r,w) = norm_v * Σ_n in_right(a,n,w) * V(gr(w),r,n)
// The LayoutLeft w-slabs in_left(:,:,w)/lproj(:,:,w) are column-major [chunk x n]/[chunk x
// rank] (ld=extent(0), w-stride=extent(0)*extent(1)); p.Uw[w] is the pre-gathered col-major
// [rank x n_left] weight (ldb=rank, w-stride=n_left*rank) built in Opt-16. OP_N/OP_T so
// op(Uw)=[n_left x rank] gives C(a,r)=Σ_n in_left(a,n)*Uw(r,n). beta=0 (full overwrite).
// True fp32 (PEDANTIC handle).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_cp_l_project_cublas(
    const DeviceCPL &p, const t_nn_3d &in_left, const t_nn_3d &in_right)
{
  const int n_lm_left  = (int) p.group_left.extent(0);
  const int n_lm_right = (int) p.group_right.extent(0);
  const float beta = 0.0f;
  if (n_lm_left > 0 && p.n_left > 0) {
    const float alpha = (float) p.norm_u;
    cublasStatus_t st = cublasSgemmStridedBatched(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_T,
        chunk_size, p.rank, p.n_left, &alpha,
        (const float*) in_left.data(), (int) in_left.extent(0),
        (long long) in_left.extent(0) * in_left.extent(1),
        (const float*) p.Uw.data(), p.rank, (long long) p.n_left * p.rank,
        &beta,
        (float*) d_cpl_lproj.data(), (int) d_cpl_lproj.extent(0),
        (long long) d_cpl_lproj.extent(0) * d_cpl_lproj.extent(1),
        n_lm_left);
    if (st != CUBLAS_STATUS_SUCCESS)
      error->all(FLERR, "GRACE-3L/KK: cp_l fwd left projection cublas failed (status {})", (int) st);
  }
  if (n_lm_right > 0 && p.n_right > 0) {
    const float alpha = (float) p.norm_v;
    cublasStatus_t st = cublasSgemmStridedBatched(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_T,
        chunk_size, p.rank, p.n_right, &alpha,
        (const float*) in_right.data(), (int) in_right.extent(0),
        (long long) in_right.extent(0) * in_right.extent(1),
        (const float*) p.Vw.data(), p.rank, (long long) p.n_right * p.rank,
        &beta,
        (float*) d_cpl_rproj.data(), (int) d_cpl_rproj.extent(0),
        (long long) d_cpl_rproj.extent(0) * d_cpl_rproj.extent(1),
        n_lm_right);
    if (st != CUBLAS_STATUS_SUCCESS)
      error->all(FLERR, "GRACE-3L/KK: cp_l fwd right projection cublas failed (status {})", (int) st);
  }
}
#endif

// ======================================================================
// Kernel: TagReduceN (Task 2.6) — generic FunctionReduceN. Reproduces
// numpy_forward_3l.reduce_n EXACTLY for one reduce r over the current chunk:
//   for each instruction ins with n_conn connections:
//     for each connection c: gather A[a, n, collect_ind[c]] over n in
//     [0,n_in), dot with W[e, k, n, w_l_tile[c]] (e = elem_dep ? mu_i : 0;
//     W stored uniformly 4D so no branch on the read), * norm;
//     only_invar -> accumulate into f=0 (ignore total_sum_ind);
//     else       -> accumulate into f=total_sum_ind[c].
//   collection *= norm_map (if present).
// One thread owns the whole out[a,k,:] row (flat index over atom a, output
// channel k) -> the per-instruction/per-connection accumulation is serial,
// no atomics needed (same ownership pattern as TagCPLCouple).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagReduceN,
    const int& idx) const
{
  const int n_out = reduceN_r.n_out;
  const int a = idx % chunk_size;
  const int k = idx / chunk_size;

  const bool only_invar = reduceN_r.only_invar != 0;
  const bool elem_dep = reduceN_r.elem_dep != 0;
  const int shape0 = only_invar ? 1 : reduceN_r.n_funcs;
  const int e = elem_dep ? d_mu_i(a) : 0;

  NNScalar acc[GRACE3L_MAX_REDUCE_NFUNCS];
  for (int f = 0; f < shape0; f++) acc[f] = NNScalar(0.0);

  for (int ins = 0; ins < reduceN_r.n_instr; ins++) {
    const auto &Wv = reduceN_r.W[ins];
    const auto &ci = reduceN_r.collect_ind[ins];
    const auto &wlt = reduceN_r.w_l_tile[ins];
    const auto &tsi = reduceN_r.total_sum_ind[ins];
    const int n_in = reduceN_r.n_in[ins];
    const int n_conn = reduceN_r.n_conn[ins];
    const NNScalar norm = reduceN_r.norm[ins];
    const auto &A = reduceN_inputs[ins];

    for (int c = 0; c < n_conn; c++) {
      const int lm = ci(c);
      const int tile = wlt(c);
      NNScalar s = NNScalar(0.0);
      for (int n = 0; n < n_in; n++)
        s += Wv(e, k, n, tile) * A(a, n, lm);
      s *= norm;
      if (only_invar) acc[0] += s;
      else acc[tsi(c)] += s;
    }
  }

  if (reduceN_r.has_norm_map)
    for (int f = 0; f < shape0; f++) acc[f] *= reduceN_r.norm_map(f);

  for (int f = 0; f < shape0; f++) reduceN_out(a, k, f) = acc[f];
}

// ======================================================================
// compute_reduceN (Task 2.6) — generic FunctionReduceN. Sets the per-call
// launch state (reduce metadata + per-instruction input views + output view)
// captured by TagReduceN, then launches one flat pass over (atom, n_out) for
// the current chunk. Reused by every A*_*_red / eq* / rho* reduce; the caller
// is responsible for matching inputs[j] to grace_model->reduce[i].instr[j]
// .name (device kernels have no strings to do that lookup themselves).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_reduceN(
    const DeviceReduce &r, const t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR],
    const t_nn_3d &out)
{
#ifdef KOKKOS_ENABLE_CUDA
  if (cublas_handle && sizeof(NNScalar) == 4) {
    const DeviceReduce *rbase = d_reduce;
    const int idx = (int) (&r - rbase);
    if (idx >= 0 && idx < n_reduces && h_reduce[idx].cublas_ok &&
        compute_reduceN_cublas(r, idx, inputs, out))
      return;
  }
#endif
  reduceN_r = r;
  for (int i = 0; i < GRACE3L_MAX_REDUCE_INSTR; i++) reduceN_inputs[i] = inputs[i];
  reduceN_out = out;
  Kokkos::parallel_for("ReduceN",
      Kokkos::RangePolicy<DeviceType, TagReduceN>(0, chunk_size * r.n_out), *this);
}

#ifdef KOKKOS_ENABLE_CUDA
// Opt 17: cuBLAS batched-GEMM ReduceN forward. out(a,k,f) = Σ_{c:tsi(c)=f} Σ_n
// A(a,n,ci(c)) * Wsc[c](k,n), where Wsc already carries norm*norm_map(f). Zero `out`,
// then issue one cublasSgemmBatched per round (each round has DISTINCT f -> the beta=1
// accumulate into disjoint out(:,:,f) slabs is race-free). A/out LayoutLeft slabs are
// column-major (ld=extent(0), slab stride=extent(0)*extent(1)); the batch selects the
// gathered slabs via a device pointer array. OP_N/OP_T (K=n_in). True fp32 (PEDANTIC).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
bool PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_reduceN_cublas(
    const DeviceReduce &r, int idx, const t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR],
    const t_nn_3d &out)
{
  const ReduceHost &hr = h_reduce[idx];
  const t_nn_3d &A = inputs[0];
  const int nc = hr.n_conn, n_out = hr.n_out, n_in = hr.n_in;
  if (nc <= 0) return true;
  Kokkos::deep_copy(out, NNScalar(0));   // forward is an assignment (some f have no conn)
  const long strideA = (long) A.extent(0) * A.extent(1);
  const long strideC = (long) out.extent(0) * out.extent(1);
  const uintptr_t baseA = (uintptr_t) A.data();
  const uintptr_t baseW = (uintptr_t) r.Wsc.data();
  const uintptr_t baseC = (uintptr_t) out.data();
  const size_t szf = sizeof(NNScalar);
  for (int i = 0; i < nc; i++) {
    const int c = hr.fwd_order[i];
    h_redptr(i)          = baseA + (uintptr_t) (strideA * hr.ci[c]) * szf;
    h_redptr(nc + i)     = baseW + (uintptr_t) ((size_t) c * n_out * n_in) * szf;
    h_redptr(2 * nc + i) = baseC + (uintptr_t) (strideC * hr.tsi[c]) * szf;
  }
  Kokkos::deep_copy(d_redptr, h_redptr);
  const float** Aarr = (const float**) d_redptr.data();
  const float** Barr = (const float**) (d_redptr.data() + nc);
  float** Carr       = (float**)       (d_redptr.data() + 2 * nc);
  const float alpha = 1.0f, beta = 1.0f;
  const int ldA = (int) A.extent(0), ldC = (int) out.extent(0);
  for (int rr = 0; rr < hr.nrounds; rr++) {
    const int off = hr.round_ptr[rr];
    const int bs = hr.round_ptr[rr + 1] - off;
    if (bs <= 0) continue;
    cublasStatus_t st = cublasSgemmBatched(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_T,
        chunk_size, n_out, n_in, &alpha,
        Aarr + off, ldA, Barr + off, n_out, &beta,
        Carr + off, ldC, bs);
    if (st != CUBLAS_STATUS_SUCCESS)
      error->all(FLERR, "GRACE-3L/KK: reduceN fwd cublasSgemmBatched failed (status {})", (int) st);
  }
  return true;
}
#endif

// ======================================================================
// Kernel: TagFCRight2Left (Task 2.6) — generic FCRight2Left. Reproduces
// numpy_forward_3l.fc_right2left EXACTLY for one FC f over the current chunk:
//   left_t[lm,a,k] = left_coefs ? norm_left * Sum_n w_left[k,n,w_tile_left[lm]]
//                                            * left[a,n,lm]     (lm in [0,n_funcs_left))
//                               : left[a,k,lm]      (left already at this FC's n_out)
//   for each connection c in [0,n_funcs_right):
//     right_t = norm_right * Sum_n w_right[k,n,w_tile_right[c]] * right[a,n,collect_from[c]]
//     left_t[collect_to[c]] += right_t
//   left_t *= norm_out_factor   (loader requires it unconditionally -> always applied)
// One thread owns the whole out[a,k,:] row (flat index over atom a, output
// channel k) -> the scatter-add over connections is serial, no atomics needed
// (same ownership pattern as TagReduceN / TagCPLCouple).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagFCRight2Left,
    const int& idx) const
{
  const int n_out = fc_p.n_out;
  const int a = idx % chunk_size;
  const int k = idx / chunk_size;

  const int n_funcs_left = fc_p.n_funcs_left;
  NNScalar acc[GRACE3L_MAX_FC_NFUNCS];

  if (fc_p.left_coefs) {
    const int n_in_left = fc_p.n_in_left;
    for (int w = 0; w < n_funcs_left; w++) {
      const int tile = fc_p.w_tile_left(w);
      NNScalar s = NNScalar(0.0);
      for (int n = 0; n < n_in_left; n++)
        s += fc_p.w_left(k, n, tile) * fc_in_left(a, n, w);
      acc[w] = s * fc_p.norm_left;
    }
  } else {
    for (int w = 0; w < n_funcs_left; w++)
      acc[w] = fc_in_left(a, k, w);
  }

  const int n_funcs_right = fc_p.n_funcs_right;
  const int n_in_right = fc_p.n_in_right;
  for (int c = 0; c < n_funcs_right; c++) {
    const int src = fc_p.collect_from(c);
    const int dst = fc_p.collect_to(c);
    const int tile = fc_p.w_tile_right(c);
    NNScalar s = NNScalar(0.0);
    for (int n = 0; n < n_in_right; n++)
      s += fc_p.w_right(k, n, tile) * fc_in_right(a, n, src);
    acc[dst] += s * fc_p.norm_right;
  }

  for (int w = 0; w < n_funcs_left; w++)
    acc[w] *= fc_p.norm_out_factor(w);

  for (int w = 0; w < n_funcs_left; w++)
    fc_out(a, k, w) = acc[w];
}

// ======================================================================
// compute_fc (Task 2.6) — generic FCRight2Left. Sets the per-call launch
// state (FC metadata + left/right input views + output view) captured by
// TagFCRight2Left, then launches one flat pass over (atom, n_out) for the
// current chunk. Reused by every A*_2a/2b FC.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_fc(
    const DeviceFC &f, const t_nn_3d &in_left, const t_nn_3d &in_right,
    const t_nn_3d &out)
{
  fc_p = f;
  fc_in_left = in_left;
  fc_in_right = in_right;
  fc_out = out;
  Kokkos::parallel_for("FCRight2Left",
      Kokkos::RangePolicy<DeviceType, TagFCRight2Left>(0, chunk_size * f.n_out), *this);
}

// ======================================================================
// Kernel: TagEqNorm (Task 2.7) — generic EquivariantRMSNorm. Reproduces
// numpy_forward_3l.equivariant_rms_norm EXACTLY for one eqnorm e over the
// current chunk. x = eqnorm_in[a, feat, lm] with feat in [0,n_out=n_feat),
// lm in [0,M):
//   center_l0: l0_mean[lm] = mean_over_feat(x[a,:,lm]) * l0_mask[lm]
//              x <- x - l0_mean[lm]                      (broadcast over feat)
//   norm[a]  = mean_over_feat( sum_over_lm( dw[lm] * x^2 ) )   -- ONE scalar/atom
//   rms[a]   = 1/sqrt(norm[a] + eps)
//   scale[feat,lm] = affine_weight[expand_index[lm], feat]
//   out[a,feat,lm] = x[a,feat,lm] * rms[a] * scale[feat,lm]
// One thread owns an entire atom (all feat*lm entries never touched by any
// other thread) -> three sequential passes over the same (feat,lm) grid
// (mean, sum-of-squares, write), re-reading eqnorm_in each time instead of
// staging the centered x in a [feat,lm]-sized stack buffer (up to 64*25 —
// too large to keep as registers). l0_mean[lm] (<=25 lm) is the only
// intermediate kept in a stack array, matching the existing acc[25] idiom
// used by TagComputeA1 for the same lm-axis bound.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagEqNorm,
    const int& a) const
{
  const int n_feat = eqnorm_p.n_out;
  const int M = eqnorm_p.M;
  const bool center_l0 = eqnorm_p.center_l0 != 0;
  const NNScalar eps = eqnorm_p.eps;

  NNScalar l0mean[GRACE3L_MAX_SPBF_NLMIND];   // M (<=64) validated at load
  for (int lm = 0; lm < M; lm++) l0mean[lm] = NNScalar(0.0);

  if (center_l0) {
    for (int lm = 0; lm < M; lm++) {
      NNScalar s = NNScalar(0.0);
      for (int f = 0; f < n_feat; f++) s += eqnorm_in(a, f, lm);
      l0mean[lm] = (s / (NNScalar) n_feat) * eqnorm_p.l0_mask(lm);
    }
  }

  NNScalar sumsq = NNScalar(0.0);
  for (int f = 0; f < n_feat; f++)
    for (int lm = 0; lm < M; lm++) {
      const NNScalar xc = eqnorm_in(a, f, lm) - l0mean[lm];
      sumsq += xc * xc * eqnorm_p.degree_weights(lm);
    }
  const NNScalar norm = sumsq / (NNScalar) n_feat;
  const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(norm + eps);

  for (int f = 0; f < n_feat; f++)
    for (int lm = 0; lm < M; lm++) {
      const NNScalar xc = eqnorm_in(a, f, lm) - l0mean[lm];
      const int grp = eqnorm_p.expand_index(lm);
      const NNScalar scale = eqnorm_p.affine_weight(grp, f);
      eqnorm_out(a, f, lm) = xc * rms * scale;
    }
}

// ======================================================================
// compute_equiv_rms_norm (Task 2.7) — generic EquivariantRMSNorm. Sets the
// per-call launch state captured by TagEqNorm, then launches one thread per
// atom for the current chunk. Reused for eq1_norm/eq2_norm.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_equiv_rms_norm(
    const DeviceEqNorm &e, const t_nn_3d &in, const t_nn_3d &out)
{
  eqnorm_p = e;
  eqnorm_in = in;
  eqnorm_out = out;
  Kokkos::parallel_for("EqNorm",
      Kokkos::RangePolicy<DeviceType, TagEqNorm>(0, chunk_size), *this);
}

// ======================================================================
// Kernel: TagEqNormBwd (Task 3.2) — generic EquivariantRMSNorm VJP. Analytic
// transpose of TagEqNorm. Forward (per atom a; x=eqnorm_in[a,f,lm], f<n_feat,
// lm<M):
//   l0mean[lm] = mean_f(x[a,:,lm]) * l0_mask[lm]     (center_l0)
//   xc         = x - l0mean
//   S[a]       = (1/F) Σ_{f,lm} dw[lm]*xc^2          (F=n_feat, ONE scalar/atom)
//   rms[a]     = 1/sqrt(S[a]+eps)
//   out        = xc*rms*scale[f,lm],  scale[f,lm]=affine_weight[expand_index[lm],f]
// VJP w.r.t x (affine is a constant weight, no weight-grad needed):
//   adj_S   = -0.5*rms^3 * Σ_{f,lm} adj_out*xc*scale
//   adj_xc  = adj_out*rms*scale + adj_S*(2/F)*dw[lm]*xc
//   adj_x   = adj_xc - l0_mask[lm]*(1/F)*Σ_f adj_xc     (center_l0 transpose)
// One thread/atom (owns all f*lm). Recomputes xc/rms from the forward input
// (not cached). ASSIGNS adj_in (eq2 has a single consumer: eq2_norm). Stack
// arrays sized to the parity-doubled lm ceiling 50 (eq2_norm M=50).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagEqNormBwd,
    const int& a) const
{
  const int n_feat = eqnorm_p.n_out;
  const int M = eqnorm_p.M;
  const bool center_l0 = eqnorm_p.center_l0 != 0;
  const NNScalar eps = eqnorm_p.eps;
  const NNScalar invF = NNScalar(1.0) / (NNScalar) n_feat;

  NNScalar l0mean[GRACE3L_MAX_SPBF_NLMIND];   // M (<=64) validated at load
  for (int lm = 0; lm < M; lm++) l0mean[lm] = NNScalar(0.0);
  if (center_l0) {
    for (int lm = 0; lm < M; lm++) {
      NNScalar s = NNScalar(0.0);
      for (int f = 0; f < n_feat; f++) s += eqnorm_in(a, f, lm);
      l0mean[lm] = s * invF * eqnorm_p.l0_mask(lm);
    }
  }

  // Recompute S/rms and adj_S = -0.5*rms^3 * Σ adj_out*xc*scale.
  NNScalar sumsq = NNScalar(0.0);
  NNScalar adot = NNScalar(0.0);
  for (int f = 0; f < n_feat; f++)
    for (int lm = 0; lm < M; lm++) {
      const NNScalar xc = eqnorm_in(a, f, lm) - l0mean[lm];
      sumsq += xc * xc * eqnorm_p.degree_weights(lm);
      const NNScalar scale = eqnorm_p.affine_weight(eqnorm_p.expand_index(lm), f);
      adot += eqnorm_adj_out(a, f, lm) * xc * scale;
    }
  const NNScalar S = sumsq * invF;
  const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(S + eps);
  const NNScalar adj_S = NNScalar(-0.5) * rms * rms * rms * adot;

  // colsum[lm] = Σ_f adj_xc[a,f,lm] (only needed for the center_l0 transpose).
  NNScalar colsum[GRACE3L_MAX_SPBF_NLMIND];   // M (<=64) validated at load
  for (int lm = 0; lm < M; lm++) colsum[lm] = NNScalar(0.0);
  if (center_l0) {
    for (int lm = 0; lm < M; lm++) {
      const NNScalar dw = eqnorm_p.degree_weights(lm);
      NNScalar cs = NNScalar(0.0);
      for (int f = 0; f < n_feat; f++) {
        const NNScalar xc = eqnorm_in(a, f, lm) - l0mean[lm];
        const NNScalar scale = eqnorm_p.affine_weight(eqnorm_p.expand_index(lm), f);
        const NNScalar adj_xc = eqnorm_adj_out(a, f, lm) * rms * scale
                              + adj_S * (NNScalar(2.0) * invF) * dw * xc;
        cs += adj_xc;
      }
      colsum[lm] = cs;
    }
  }

  for (int f = 0; f < n_feat; f++)
    for (int lm = 0; lm < M; lm++) {
      const NNScalar dw = eqnorm_p.degree_weights(lm);
      const NNScalar xc = eqnorm_in(a, f, lm) - l0mean[lm];
      const NNScalar scale = eqnorm_p.affine_weight(eqnorm_p.expand_index(lm), f);
      const NNScalar adj_xc = eqnorm_adj_out(a, f, lm) * rms * scale
                            + adj_S * (NNScalar(2.0) * invF) * dw * xc;
      // center_l0 transpose (colsum==0 when !center_l0 -> pure identity).
      eqnorm_adj_in(a, f, lm) = adj_xc - eqnorm_p.l0_mask(lm) * invF * colsum[lm];
    }
}

// compute_equiv_rms_norm_bwd (Task 3.2) — sets the per-call state captured by
// TagEqNormBwd, then launches one thread per atom. ASSIGNS adj_in.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_equiv_rms_norm_bwd(
    const DeviceEqNorm &e, const t_nn_3d &in, const t_nn_3d &adj_out,
    const t_nn_3d &adj_in)
{
  eqnorm_p = e;
  eqnorm_in = in;
  eqnorm_adj_out = adj_out;
  eqnorm_adj_in = adj_in;
  Kokkos::parallel_for("EqNormBwd",
      Kokkos::RangePolicy<DeviceType, TagEqNormBwd>(0, chunk_size), *this);
}

// ======================================================================
// Kernel: TagRMSNorm (Task 2.7) — generic InvariantLayerRMSNorm. Reproduces
// numpy_forward_3l.invariant_rms_norm EXACTLY for one rmsnorm r over the
// current chunk, at one (atom, lm) column x[a,:,lm]:
//   full        (scale_len==n_out):   rms = 1/sqrt(mean_over_ch(x^2)+eps);
//                                      out[ch] = x[ch]*rms*scale[ch]
//   only_nonlin (scale_len==n_out-1): out[0] = x[0] (passthrough); rms over
//                                      channels [1,n_out) only;
//                                      out[ch] = x[ch]*rms*scale[ch-1], ch>=1
// eps=1e-10 fixed (oracle hardcodes F32(1e-10), not an npz-loaded value).
// One thread owns the whole column x[a,:,lm] (flat index over atom a, lm) ->
// no atomics needed (same ownership pattern as TagReduceN/TagFCRight2Left).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagRMSNorm,
    const int& idx) const
{
  const int n_lm = (int) rmsnorm_in.extent(2);
  const int a = idx / n_lm;
  const int lm = idx - a * n_lm;

  const int n_out = rmsnorm_p.n_out;
  const NNScalar eps = NNScalar(1e-10);

  if (rmsnorm_p.scale_len == n_out) {
    // full
    NNScalar sumsq = NNScalar(0.0);
    for (int c = 0; c < n_out; c++) {
      const NNScalar v = rmsnorm_in(a, c, lm);
      sumsq += v * v;
    }
    const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(sumsq / (NNScalar) n_out + eps);
    for (int c = 0; c < n_out; c++)
      rmsnorm_out(a, c, lm) = rmsnorm_in(a, c, lm) * rms * rmsnorm_p.scale(c);
  } else {
    // only_nonlin: scale_len == n_out - 1; channel 0 passthrough.
    rmsnorm_out(a, 0, lm) = rmsnorm_in(a, 0, lm);
    const int n_nl = n_out - 1;
    NNScalar sumsq = NNScalar(0.0);
    for (int c = 1; c < n_out; c++) {
      const NNScalar v = rmsnorm_in(a, c, lm);
      sumsq += v * v;
    }
    const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(sumsq / (NNScalar) n_nl + eps);
    for (int c = 1; c < n_out; c++)
      rmsnorm_out(a, c, lm) = rmsnorm_in(a, c, lm) * rms * rmsnorm_p.scale(c - 1);
  }
}

// ======================================================================
// compute_invariant_rms_norm (Task 2.7) — generic InvariantLayerRMSNorm.
// Sets the per-call launch state captured by TagRMSNorm, then launches one
// flat pass over (atom, lm) for the current chunk. Reused for
// rho1_norm/rho2_norm/rho3_norm.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_invariant_rms_norm(
    const DeviceRMSNorm &r, const t_nn_3d &in, const t_nn_3d &out)
{
  rmsnorm_p = r;
  rmsnorm_in = in;
  rmsnorm_out = out;
  const int n_lm = (int) in.extent(2);
  Kokkos::parallel_for("RMSNorm",
      Kokkos::RangePolicy<DeviceType, TagRMSNorm>(0, chunk_size * n_lm), *this);
}

// ======================================================================
// Kernel: ComputeMLPEnergy (Task 2.11) — 3-density readout -> per-atom energy.
// Oracle: numpy_forward_3l.py readout() (out1_origins = [rho1_norm, rho2_norm,
// rho3_norm]). For each origin, s = origin[:, :, 0] (the l=0 slice, LAST axis
// index 0); channel 0 is a weightless linear skip, channels 1..31 feed the
// MLP(31->64 silu ->1). ae_nn = e_mlp + lin; E_atom = ae_nn*output_scale +
// shift[species]. rho3_norm is a chunk-scratch buffer indexed by the SAME
// chunk-local ii used by the L3 rho kernels (written+read within this SAME P3
// chunk, so that is safe). rho1_norm/rho2_norm are instead read from the natom
// buffers d_rho{1,2}_norm_full at ii+chunk_offset (Task 3.5a) — the chunk-local
// d_rho1_norm/d_rho2_norm were written back in P1/P2 and, under multi-chunk,
// hold only the LAST chunk by the time P3 runs, so reading them chunk-local
// here would silently pick up the wrong atoms' densities for every chunk but
// the last. Backward (readout_bwd, TagReadoutBwd below) is UNCHANGED — it
// still reads d_rho1_norm/d_rho2_norm chunk-local, which is fine only in the
// single-chunk regime and stays a known limitation until Task 3.5b.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeMLPEnergy,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>::member_type& team) const
{
  const int ii = team.league_rank();
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  // Linear (weightless) skip: channel 0, l=0 slice, summed over the 3 origins.
  const NNScalar linear_term = d_rho1_norm_full(ii + chunk_offset, 0, 0) +
                                d_rho2_norm_full(ii + chunk_offset, 0, 0) +
                                d_rho3_norm(ii, 0, 0);

  NNScalar ebuf_a[GRACE3L_MAX_ENERGY_MLP_DIM], ebuf_b[GRACE3L_MAX_ENERGY_MLP_DIM];
  NNScalar* h_cur = ebuf_a;
  NNScalar* h_nxt = ebuf_b;

  // MLP input: channels 1..31, l=0 slice, summed over the 3 origins.
  const int nin0 = d_energy_dims(0);
  for (int c = 0; c < nin0; c++)
    h_cur[c] = d_rho1_norm_full(ii + chunk_offset, c + 1, 0) +
               d_rho2_norm_full(ii + chunk_offset, c + 1, 0) +
               d_rho3_norm(ii, c + 1, 0);

  // Hidden layers with activation (3L uses silu: energy_activation == 0).
  for (int layer = 0; layer < energy_n_layers - 1; layer++) {
    const int nin = d_energy_dims(layer);
    const int nout = d_energy_dims(layer + 1);
    const NNScalar norm = d_energy_norms(layer);
    for (int j = 0; j < nout; j++) {
      NNScalar sum = 0.0;
      for (int k = 0; k < nin; k++)
        sum += d_energy_W(layer, k, j) * h_cur[k];
      sum *= norm;
      if (energy_activation == 1)
        h_nxt[j] = tanh_act(sum);
      else
        h_nxt[j] = silu(sum);
    }
    NNScalar* tmp = h_cur; h_cur = h_nxt; h_nxt = tmp;
  }

  // Output layer (no activation).
  const int last = energy_n_layers - 1;
  const int nin_last = d_energy_dims(last);
  const int nout_last = d_energy_dims(last + 1);
  NNScalar e = 0.0;
  for (int o = 0; o < nout_last; o++) {
    NNScalar sum = 0.0;
    for (int k = 0; k < nin_last; k++)
      sum += d_energy_W(last, k, o) * h_cur[k];
    e += sum * d_energy_norms(last);
  }

  // E_atom = (MLP_out + linear_term) * output_scale + shift[species]; the MLP
  // math runs in NNScalar (matching TF's fp32 readout), the stored per-atom
  // result is cast to fp64 (d_e_atom is KK_ACC_FLOAT regardless of NNScalar).
  d_e_atom(ii) = (KK_ACC_FLOAT)((e + linear_term) * output_scale + d_shifts(mu_i));
}

// ======================================================================
// UQ basis-RP kernels (forward only). The GMM feature is the projection of the
// invariant (l=0) products entering the three per-layer scalar reduces
// rho1|rho2|rho3, read in the EXACT (source, n, collect) order TagReduceN
// iterates them (n outer, connection inner), so R's row cursor advances in the
// concatenation order [rho1, rho2, rho3]. Because feature_transform is identity
// and normalize is a scalar divide, the projection proj = B·R is linear and is
// stream-accumulated across the three layer phases into d_uq_z; the final kernel
// normalizes by ||B|| = sqrt(Σ n2) and appends the [full, rho1, rho2, rho3]
// density channels before the GMM. FP64 throughout for numerical safety.
// ======================================================================

// Stream one rho block's sources into the RAW projection ACC and its block norm² N2.
// For each gathered product b_inv = A(ii, n, ci(c)) (n outer, connection inner —
// TagReduceN order), accumulate ACC[d] += b_inv*R[B,d] over the RP columns, add
// b_inv² into N2, and advance the R-row cursor B. Expands against each kernel's
// local RP/ii, the d_uq_rp_matrix member, and the reduce metadata.
#define GRACE3L_UQ_ACC_BLOCK(REDUCE, VIEWS, ROLE, ACC, N2, B)                       \
    { const DeviceReduce &_r = (REDUCE);                                            \
      for (int _j = 0; _j < _r.n_instr; _j++) {                                     \
        const t_nn_3d &_A = (VIEWS)[(ROLE)[_j]];                                    \
        const int _nin = _r.n_in[_j];                                              \
        const int _nc = _r.n_conn[_j];                                             \
        const auto &_ci = _r.collect_ind[_j];                                      \
        for (int _n = 0; _n < _nin; _n++)                                          \
          for (int _c = 0; _c < _nc; _c++) {                                       \
            const double _bv = (double)(_A)(ii, _n, _ci(_c));                       \
            (N2) += _bv * _bv;                                                      \
            for (int _d = 0; _d < RP; _d++) (ACC)[_d] += _bv * d_uq_rp_matrix((B), _d); \
            (B)++;                                                                  \
          } } }

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeUQ_rho1, const int& ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int RP = uq_rp_dim;
  double z[GRACE3L_UQ_MAX_RP_DIM];
  for (int d = 0; d < RP; d++) z[d] = 0.0;
  double n2_rho1 = 0.0;
  // rho1 block occupies R-rows [0, d_basis_rho1); cursor starts at 0.
  int b = 0;
  t_nn_3d views[4] = {d_A1, d_A1_2_red, d_A1_3_red, d_A1_4_red};
  GRACE3L_UQ_ACC_BLOCK(d_reduce[idx_reduce_rho1], views, rho1_input_role, z, n2_rho1, b);
  for (int d = 0; d < RP; d++) d_uq_z(i, d) = z[d];
  d_uq_n2_rho1(i) = n2_rho1;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeUQ_rho2, const int& ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int RP = uq_rp_dim;
  double z[GRACE3L_UQ_MAX_RP_DIM];
  for (int d = 0; d < RP; d++) z[d] = d_uq_z(i, d);   // carry rho1 projection
  double n2_rho2 = 0.0;
  // rho2 block occupies R-rows [d_basis_rho1, d_basis_rho1 + d_basis_rho2).
  int b = uq_d_basis_rho1;
  t_nn_3d views[4] = {d_A2_red, d_A2_2_red, d_A2_3_red, d_A2_4_red};
  GRACE3L_UQ_ACC_BLOCK(d_reduce[idx_reduce_rho2], views, rho2_input_role, z, n2_rho2, b);
  for (int d = 0; d < RP; d++) d_uq_z(i, d) = z[d];
  d_uq_n2_rho2(i) = n2_rho2;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeUQ, const int& ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));
  const int D = uq_D;          // full feature dim (rp_dim + n_density)
  const int RP = uq_rp_dim;    // projection width (R cols)

  // RAW projection proj = B·R: rho1 (P1) + rho2 (P2) from d_uq_z + rho3 here.
  double f[GRACE3L_UQ_MAX_RP_DIM];
  for (int d = 0; d < RP; d++) f[d] = d_uq_z(i, d);
  const double n2_rho1 = d_uq_n2_rho1(i);
  const double n2_rho2 = d_uq_n2_rho2(i);
  double n2_rho3 = 0.0;
  // rho3 block occupies R-rows [d_basis_rho1 + d_basis_rho2, d_basis).
  int b = uq_d_basis_rho1 + uq_d_basis_rho2;
  t_nn_3d views[4] = {d_A3_red, d_A3_2_red, d_A3_3_red, d_A3_4_red};
  GRACE3L_UQ_ACC_BLOCK(d_reduce[idx_reduce_rho3], views, rho3_input_role, f, n2_rho3, b);
  #undef GRACE3L_UQ_ACC_BLOCK

  // --- L2-normalize the projection: proj = (B/||B||)·R = (B·R)/||B|| ---
  const double n2_full = n2_rho1 + n2_rho2 + n2_rho3;
  const double nrm_full = Kokkos::sqrt(n2_full);   // ||B||, reused below
  if (uq_normalize) {
    const double inv = 1.0 / (nrm_full + 1.0e-12);
    for (int d = 0; d < RP; d++) f[d] *= inv;
  }

  // --- Density channels appended after the projection: [full, rho1, rho2, rho3] ---
  // Block order is the sorted-by-reduce-name order (rho1 < rho2 < rho3), matching
  // the R-row layout and the reference centroid density block.
  if (uq_n_density > 0) f[RP + 0] = Kokkos::log(nrm_full + 1.0e-12) * uq_density_scale;
  if (uq_n_density > 1) f[RP + 1] = Kokkos::log(Kokkos::sqrt(n2_rho1) + 1.0e-12) * uq_density_scale;
  if (uq_n_density > 2) f[RP + 2] = Kokkos::log(Kokkos::sqrt(n2_rho2) + 1.0e-12) * uq_density_scale;
  if (uq_n_density > 3) f[RP + 3] = Kokkos::log(Kokkos::sqrt(n2_rho3) + 1.0e-12) * uq_density_scale;

  // --- Cluster assignment: nearest centroid (squared Euclidean) ---
  // ncl > 0 for every element is guaranteed at load time (uqv6 artifacts cover all elements).
  const int ncl = d_uq_n_clusters(mu_i);
  int kstar = 0;
  double best = 1.0e300;
  for (int k = 0; k < ncl; k++) {
    double d2 = 0.0;
    for (int p = 0; p < D; p++) {
      const double diff = f[p] - d_uq_centroids(mu_i, k, p);
      d2 += diff * diff;
    }
    if (d2 < best) { best = d2; kstar = k; }
  }

  // --- Mahalanobis distance to assigned cluster ---
  double delta[GRACE3L_UQ_MAX_RP_DIM];
  for (int p = 0; p < D; p++) delta[p] = f[p] - d_uq_centroids(mu_i, kstar, p);
  double sig2 = 0.0;
  for (int p = 0; p < D; p++) {
    double acc = 0.0;
    for (int q = 0; q < D; q++)
      acc += d_uq_inv_cov(mu_i, kstar, p, q) * delta[q];
    sig2 += delta[p] * acc;
  }
  const double sigma = Kokkos::sqrt(sig2 + 1.0e-8);

  double t = d_uq_interp_thresholds(mu_i, kstar);
  if (t < 1.0e-10) t = 1.0e-10;

  d_sigma(i) = sigma;
  d_gamma(i) = sigma / t;
  d_gmm_cluster(i) = (double) kstar;
}

// ======================================================================
// ============ Task 3.1: readout + Layer-3 backward (adjoints) ==========
// GENERIC reverse-mode VJP kernels + launchers (reused by Tasks 3.2/3.3 for
// L2/L1). Each is the analytic transpose of the matching forward routine.
// ======================================================================

// ---- TagReadoutBwd: VJP of TagComputeMLPEnergy. Seed dE/dE_atom=1. Recomputes
// the hidden pre-activations from the SUMMED 3-density input, backprops through
// the energy MLP, then writes d_rho{1,2,3}_norm_adj (all EQUAL — the 3 densities
// enter symmetrically): ch0 = output_scale (linear skip), ch (c+1) = adj_input[c],
// every lm>0 slot = 0. Mirrors GRACE-2L ReverseMLPEnergy generalized to 3 origins.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagReadoutBwd,
    const int& ii) const
{
  if (ii >= chunk_size) return;

  NNScalar pre_act_buf[GRACE3L_MAX_MLP_LAYERS * GRACE3L_MAX_ENERGY_MLP_DIM];
  NNScalar adj_a[GRACE3L_MAX_ENERGY_MLP_DIM], adj_b[GRACE3L_MAX_ENERGY_MLP_DIM];
  NNScalar* pre_act = pre_act_buf;
  NNScalar* adj_cur = adj_a;
  NNScalar* adj_nxt = adj_b;

  const int emax = energy_max_dim;
  const int nin0 = d_energy_dims(0);

  // Forward recompute: layer-0 pre-activations from the summed 3-density input.
  {
    const int nout = d_energy_dims(1);
    for (int j = 0; j < nout; j++) {
      NNScalar s = NNScalar(0.0);
      for (int k = 0; k < nin0; k++) {
        const NNScalar h = d_rho1_norm_full(ii + chunk_offset, k + 1, 0) + d_rho2_norm_full(ii + chunk_offset, k + 1, 0) + d_rho3_norm(ii, k + 1, 0);
        s += d_energy_W(0, k, j) * h;
      }
      pre_act[0 * emax + j] = s * d_energy_norms(0);
    }
  }
  for (int layer = 1; layer < energy_n_layers - 1; layer++) {
    const int nin = d_energy_dims(layer);
    const int nout = d_energy_dims(layer + 1);
    for (int j = 0; j < nout; j++) {
      NNScalar s = NNScalar(0.0);
      for (int k = 0; k < nin; k++) {
        const NNScalar pa = pre_act[(layer - 1) * emax + k];
        const NNScalar act = (energy_activation == 1) ? tanh_act(pa) : silu(pa);
        s += d_energy_W(layer, k, j) * act;
      }
      pre_act[layer * emax + j] = s * d_energy_norms(layer);
    }
  }

  // Backward: output layer, seeded by output_scale (E = (e+lin)*output_scale+shift).
  {
    const int last = energy_n_layers - 1;
    const int nin = d_energy_dims(last);
    const int nout = d_energy_dims(last + 1);
    for (int k = 0; k < nin; k++) {
      NNScalar v = NNScalar(0.0);
      for (int o = 0; o < nout; o++) v += d_energy_W(last, k, o) * d_energy_norms(last);
      adj_cur[k] = v * output_scale;
    }
  }
  for (int layer = energy_n_layers - 2; layer >= 0; layer--) {
    const int nin = d_energy_dims(layer);
    const int nout = d_energy_dims(layer + 1);
    for (int j = 0; j < nout; j++) {
      const NNScalar pa = pre_act[layer * emax + j];
      if (energy_activation == 1) {
        const NNScalar t = tanh_act(pa);
        adj_cur[j] *= (NNScalar(1.0) - t * t);            // tanh'
      } else {
        const NNScalar sig = NNScalar(1.0) / (NNScalar(1.0) + Kokkos::exp(-pa));
        adj_cur[j] *= sig * (NNScalar(1.0) + pa * (NNScalar(1.0) - sig));  // silu'
      }
    }
    for (int k = 0; k < nin; k++) {
      NNScalar v = NNScalar(0.0);
      for (int j = 0; j < nout; j++) v += adj_cur[j] * d_energy_W(layer, k, j);
      adj_nxt[k] = v * d_energy_norms(layer);
    }
    NNScalar* tmp = adj_cur; adj_cur = adj_nxt; adj_nxt = tmp;
  }

  // Write the shared readout adjoint into all 3 densities (identical). Zero every
  // channel/lm first (n_lm==1 for rho, but stay general), then set ch0 + ch1..nin0.
  const int nch = (int) d_rho3_norm_adj.extent(1);
  const int n_lm = (int) d_rho3_norm_adj.extent(2);
  // Task 3.5b: rho1/rho2 readout adjoints are natom-promoted -> write the dense
  // list-order row (ii+chunk_offset) so they survive until P5/P4; rho3 stays
  // chunk-local (ii, consumed within THIS P3 chunk).
  const int io = ii + chunk_offset;
  for (int c = 0; c < nch; c++)
    for (int lm = 0; lm < n_lm; lm++) {
      d_rho1_norm_adj(io, c, lm) = NNScalar(0.0);
      d_rho2_norm_adj(io, c, lm) = NNScalar(0.0);
      d_rho3_norm_adj(ii, c, lm) = NNScalar(0.0);
    }
  d_rho1_norm_adj(io, 0, 0) = output_scale;
  d_rho2_norm_adj(io, 0, 0) = output_scale;
  d_rho3_norm_adj(ii, 0, 0) = output_scale;
  for (int c = 0; c < nin0; c++) {
    d_rho1_norm_adj(io, c + 1, 0) = adj_cur[c];
    d_rho2_norm_adj(io, c + 1, 0) = adj_cur[c];
    d_rho3_norm_adj(ii, c + 1, 0) = adj_cur[c];
  }
}

// ---- TagRMSNormBwd: VJP of TagRMSNorm (InvariantLayerRMSNorm). Per (atom, lm)
// column, recompute rms from the forward input, dot = Σ_c adj_out·out (forward
// output), and dx[c] = rms*(adj_out[c]*scale[c] - x[c]*dot*rms/N). Both branches
// (full / only_nonlin, ch0 passthrough) selected by scale_len. eps=1e-10 (3L).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagRMSNormBwd,
    const int& idx) const
{
  const int n_lm = (int) rmsnorm_adj_in.extent(2);
  const int a = idx / n_lm;
  const int lm = idx - a * n_lm;

  const int n_out = rmsnorm_p.n_out;
  const NNScalar eps = NNScalar(1e-10);

  if (rmsnorm_p.scale_len == n_out) {
    // full
    NNScalar sumsq = NNScalar(0.0);
    for (int c = 0; c < n_out; c++) {
      const NNScalar v = rmsnorm_in(a, c, lm);
      sumsq += v * v;
    }
    const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(sumsq / (NNScalar) n_out + eps);
    NNScalar dot = NNScalar(0.0);
    for (int c = 0; c < n_out; c++)
      dot += rmsnorm_adj_out(a, c, lm) * rmsnorm_out(a, c, lm);
    for (int c = 0; c < n_out; c++) {
      const NNScalar x = rmsnorm_in(a, c, lm);
      rmsnorm_adj_in(a, c, lm) =
          rms * (rmsnorm_adj_out(a, c, lm) * rmsnorm_p.scale(c) - x * dot * rms / (NNScalar) n_out);
    }
  } else {
    // only_nonlin: channel 0 passthrough, RMS over channels [1,n_out).
    rmsnorm_adj_in(a, 0, lm) = rmsnorm_adj_out(a, 0, lm);
    const int n_nl = n_out - 1;
    NNScalar sumsq = NNScalar(0.0);
    for (int c = 1; c < n_out; c++) {
      const NNScalar v = rmsnorm_in(a, c, lm);
      sumsq += v * v;
    }
    const NNScalar rms = NNScalar(1.0) / Kokkos::sqrt(sumsq / (NNScalar) n_nl + eps);
    NNScalar dot = NNScalar(0.0);
    for (int c = 1; c < n_out; c++)
      dot += rmsnorm_adj_out(a, c, lm) * rmsnorm_out(a, c, lm);
    for (int c = 1; c < n_out; c++) {
      const NNScalar x = rmsnorm_in(a, c, lm);
      rmsnorm_adj_in(a, c, lm) =
          rms * (rmsnorm_adj_out(a, c, lm) * rmsnorm_p.scale(c - 1) - x * dot * rms / (NNScalar) n_nl);
    }
  }
}

// ---- TagReduceNBwd: VJP of TagReduceN, for ONE instruction (reduceN_bwd_instr).
// One thread per (a, n_in). adj_input[a,n,collect_ind[c]] += norm*norm_map[dst]*
// Σ_k adj_out[a,k,dst]*W[e,k,n,w_l_tile[c]], accumulated serially over c (same-fi
// collisions handled within the owning thread; += onto the pre-zeroed buffer to
// accumulate a tensor's fan-out across reduces/products).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagReduceNBwd,
    const int a, const int n) const
{
  const int ins = reduceN_bwd_instr;
  const bool only_invar = reduceN_r.only_invar != 0;
  const bool elem_dep = reduceN_r.elem_dep != 0;
  const int n_out = reduceN_r.n_out;
  const int e = elem_dep ? d_mu_i(a) : 0;

  const auto &Wv = reduceN_r.W[ins];
  const auto &ci = reduceN_r.collect_ind[ins];
  const auto &wlt = reduceN_r.w_l_tile[ins];
  const auto &tsi = reduceN_r.total_sum_ind[ins];
  const int n_conn = reduceN_r.n_conn[ins];
  const NNScalar norm = reduceN_r.norm[ins];
  const bool hnm = reduceN_r.has_norm_map;

  for (int c = 0; c < n_conn; c++) {
    const int fi = ci(c);
    const int tile = wlt(c);
    const int dst = only_invar ? 0 : tsi(c);
    NNScalar v = NNScalar(0.0);
    for (int k = 0; k < n_out; k++)
      v += reduceN_adj_out(a, k, dst) * Wv(e, k, n, tile);
    NNScalar factor = norm;
    if (hnm) factor *= reduceN_r.norm_map(dst);
    reduceN_adj_in(a, n, fi) += v * factor;
  }
}

// ---- TagFCBwdLeft: VJP of the FCRight2Left left path. left_coefs=false (identity,
// i==output channel k): adj_left[a,k,w] += adj_out[a,k,w]*norm_out_factor[w].
// left_coefs=true (i==input feature n): adj_left[a,n,w] += norm_left*norm_out_factor[w]
// * Σ_k w_left[k,n,w_tile_left[w]]*adj_out[a,k,w]. Owning (a,i,w) thread -> no atomic.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagFCBwdLeft,
    const int a, const int i, const int w) const
{
  const NNScalar nof = fc_p.norm_out_factor(w);
  if (!fc_p.left_coefs) {
    fc_adj_left(a, i, w) += fc_adj_out(a, i, w) * nof;
  } else {
    const int n_out = fc_p.n_out;
    const int tile = fc_p.w_tile_left(w);
    NNScalar s = NNScalar(0.0);
    for (int k = 0; k < n_out; k++)
      s += fc_p.w_left(k, i, tile) * fc_adj_out(a, k, w);
    fc_adj_left(a, i, w) += s * fc_p.norm_left * nof;
  }
}

// ---- TagFCBwdRight: VJP of the FCRight2Left right path. For connection c
// (src=collect_from[c], dst=collect_to[c]): adj_right[a,n,src] += norm_right*
// norm_out_factor[dst]*Σ_k w_right[k,n,w_tile_right[c]]*adj_out[a,k,dst]. Multiple
// c may share src -> atomic_add.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagFCBwdRight,
    const int a, const int c, const int n) const
{
  const int src = fc_p.collect_from(c);
  const int dst = fc_p.collect_to(c);
  const int tile = fc_p.w_tile_right(c);
  const NNScalar nof = fc_p.norm_out_factor(dst);
  const int n_out = fc_p.n_out;
  NNScalar s = NNScalar(0.0);
  for (int k = 0; k < n_out; k++)
    s += fc_p.w_right(k, n, tile) * fc_adj_out(a, k, dst);
  Kokkos::atomic_add(&fc_adj_right(a, n, src), s * fc_p.norm_right * nof);
}

// ---- TagCPLBwdCouple: VJP of the cp_l CG couple (Pass B). One thread per (a,r);
// recomputed lproj/rproj live in d_cpl_lproj/d_cpl_rproj (forward Pass A re-run
// just before). adj_p[c]=adj_out[a,r,m_sum_ind[c]]*cg[c]; adj_lproj[left_ind[c]] +=
// adj_p*rproj[right_ind[c]]; adj_rproj[right_ind[c]] += adj_p*lproj[left_ind[c]].
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLBwdCouple,
    const int& idx) const
{
  const int rank = cpl_p.rank;
  const int a = idx / rank;
  const int r = idx - a * rank;
  const int nlm = cpl_nlm;

  for (int w = 0; w < nlm; w++) {
    d_cpl_lproj_adj(a, r, w) = NNScalar(0.0);
    d_cpl_rproj_adj(a, r, w) = NNScalar(0.0);
  }
  const int n_cg = cpl_p.n_cg;
  for (int c = 0; c < n_cg; c++) {
    const int li = cpl_p.left_ind(c);
    const int ri = cpl_p.right_ind(c);
    const int fsum = cpl_p.m_sum_ind(c);
    const NNScalar d = cpl_adj_out(a, r, fsum) * cpl_p.cg(c);
    d_cpl_lproj_adj(a, r, li) += d * d_cpl_rproj(a, r, ri);
    d_cpl_rproj_adj(a, r, ri) += d * d_cpl_lproj(a, r, li);
  }
}

// ---- TagCPLBwdProject: VJP of the cp_l U/V projection (Pass A). One thread per
// (a,n,w). adj_left[a,n,w] += norm_u*Σ_r U[gl[w],r,n]*adj_lproj[a,r,w]; adj_right
// analogous with V/gr/adj_rproj. For SELF-products adj_left==adj_right -> the two
// += accumulate into the one buffer (serial per (a,n,w) thread).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLBwdProject,
    const int& idx) const
{
  const int nmax = cpl_bwd_nmax;
  const int nlm = cpl_nlm;
  const int rank = cpl_p.rank;
  const int a = idx / (nmax * nlm);
  const int rem = idx - a * (nmax * nlm);
  const int n = rem / nlm;
  const int w = rem - n * nlm;

  const int n_lm_left = (int) cpl_p.group_left.extent(0);
  if (w < n_lm_left && n < cpl_p.n_left) {
    const int gl = cpl_p.group_left(w);
    NNScalar s = NNScalar(0.0);
    for (int r = 0; r < rank; r++)
      s += cpl_p.U(gl, r, n) * d_cpl_lproj_adj(a, r, w);
    cpl_adj_left(a, n, w) += s * cpl_p.norm_u;
  }
  const int n_lm_right = (int) cpl_p.group_right.extent(0);
  if (w < n_lm_right && n < cpl_p.n_right) {
    const int gr = cpl_p.group_right(w);
    NNScalar s = NNScalar(0.0);
    for (int r = 0; r < rank; r++)
      s += cpl_p.V(gr, r, n) * d_cpl_rproj_adj(a, r, w);
    cpl_adj_right(a, n, w) += s * cpl_p.norm_v;
  }
}

// ---- TagComputeAdjProd (Opt 4): compute adj_prod = inv_avg*Σ_f A_adj(ii,n,f)*cg_W(p,f)
// into the shared d_adj_prod buffer, ONE thread per OUTPUT element (ii,n,p). Unlike the
// previous per-(ii,n) computation (a register-spilling 1600-float local array in EACH of
// A{2,3}SPBFBwd and ForceEquiv), this holds a single accumulator per thread — no spill —
// and both consumers then just READ d_adj_prod. Opt 9: the f-reduction walks only row
// p's ~3-6 nonzeros (CSR) instead of the dense nfunc range (cg_W is ~99% zeros).
// force_spbf_idx/force_A_adj set by the launcher.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeAdjProd,
    const int ii, const int n, const int p) const
{
  const DeviceSPBF &sp = d_spbf[force_spbf_idx];
  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;
  // Opt 9: iterate only row p's nonzeros (CSR, ascending f = the dense loop's
  // accumulation order — numerically identical; cg_W is ~99% zeros).
  NNScalar s = NNScalar(0.0);
  const int e1 = sp.cg_row_ptr(p + 1);
  for (int e = sp.cg_row_ptr(p); e < e1; e++)
    s += force_A_adj(ii, n, sp.cg_col_idx(e)) * sp.cg_row_val(e);
  d_adj_prod(ii, n, p) = s * inv_avg;
}

// ---- TagA3SPBFBwd: VJP of TagComputeA3 down to the L3 indicator eq2_norm. One
// thread per (owned atom ii, radial n). adj_prod[p] = inv_avg*Σ_f d_A3_adj[ii,n,f]
// *cg_W[p,f]; then for each neighbor j and lr: d_eq2_norm_adj[j,n,lr] += Σ_ly
// a_nl(ii,j,n)[ly]*adj_prod[ly*n_lm_ind+lr] (atomic; j owned OR ghost). The
// species-injection term (z_proj·chem_l0_mask) is constant -> no indicator adj.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagA3SPBFBwd,
    const int ii, const int n) const
{
  const DeviceSPBF &sp = d_spbf[idx_spbf_A3];
  if (n >= sp.n_rad_max) return;

  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int n_lm_ind = sp.n_lm_ind;
  // Opt 4: adj_prod precomputed into d_adj_prod by compute_adj_prod(idx_spbf_A3).

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int j_global = d_nearest(ii, jj);
    const NNScalar env = (NNScalar) d_env(ii, jj);
    NNScalar a_nl[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int l = sp.l_tile(ly);
      a_nl[ly] = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, ly) * env;
    }
    for (int lr = 0; lr < n_lm_ind; lr++) {
      NNScalar s = NNScalar(0.0);
      for (int ly = 0; ly < nlm_y; ly++)
        s += a_nl[ly] * d_adj_prod(ii, n, ly * n_lm_ind + lr);
      Kokkos::atomic_add(&d_eq2_norm_adj(j_global, n, lr), s);
    }
  }
}

// ---------------------- backward launchers ----------------------

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::readout_bwd()
{
  Kokkos::parallel_for("ReadoutBwd",
      Kokkos::RangePolicy<DeviceType, TagReadoutBwd>(0, chunk_size), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_invariant_rms_norm_bwd(
    const DeviceRMSNorm &r, const t_nn_3d &in, const t_nn_3d &out,
    const t_nn_3d &adj_out, const t_nn_3d &adj_in)
{
  rmsnorm_p = r;
  rmsnorm_in = in;
  rmsnorm_out = out;
  rmsnorm_adj_out = adj_out;
  rmsnorm_adj_in = adj_in;
  const int n_lm = (int) in.extent(2);
  Kokkos::parallel_for("RMSNormBwd",
      Kokkos::RangePolicy<DeviceType, TagRMSNormBwd>(0, chunk_size * n_lm), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_reduceN_bwd(
    const DeviceReduce &r, const t_nn_3d &adj_out,
    const t_nn_3d adj_inputs[GRACE3L_MAX_REDUCE_INSTR])
{
#ifdef KOKKOS_ENABLE_CUDA
  if (cublas_handle && sizeof(NNScalar) == 4) {
    const DeviceReduce *rbase = d_reduce;
    const int idx = (int) (&r - rbase);
    if (idx >= 0 && idx < n_reduces && h_reduce[idx].cublas_ok &&
        compute_reduceN_bwd_cublas(r, idx, adj_out, adj_inputs[0]))
      return;
  }
#endif
  reduceN_r = r;
  reduceN_adj_out = adj_out;
  for (int ins = 0; ins < r.n_instr; ins++) {
    reduceN_bwd_instr = ins;
    reduceN_adj_in = adj_inputs[ins];
    const int n_in = r.n_in[ins];
    Kokkos::parallel_for("ReduceNBwd",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>, TagReduceNBwd>({0, 0}, {chunk_size, n_in}), *this);
  }
}

#ifdef KOKKOS_ENABLE_CUDA
// Opt 17: cuBLAS batched-GEMM ReduceN backward (VJP). adj_in(a,n,ci(c)) += Σ_k
// adj_out(a,k,tsi(c)) * Wsc[c](k,n), Wsc carrying norm*norm_map(tsi(c)) (same fold as
// forward). collect_ind is a bijection (enforced at setup: cublas_ok is cleared for any
// reduce with duplicate collect_ind) -> every connection writes a DISTINCT adj_in slab,
// so a single beta=1 batch accumulates into the (pre-zeroed, fan-in) adj_in with no
// scatter and no race. OP_N/OP_N (K=n_out). True fp32 (PEDANTIC handle).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
bool PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_reduceN_bwd_cublas(
    const DeviceReduce &r, int idx, const t_nn_3d &adj_out, const t_nn_3d &adj_in)
{
  const ReduceHost &hr = h_reduce[idx];
  const int nc = hr.n_conn, n_out = hr.n_out, n_in = hr.n_in;
  if (nc <= 0) return true;
  const long strideAo = (long) adj_out.extent(0) * adj_out.extent(1);
  const long strideAi = (long) adj_in.extent(0) * adj_in.extent(1);
  const uintptr_t baseAo = (uintptr_t) adj_out.data();
  const uintptr_t baseW  = (uintptr_t) r.Wsc.data();
  const uintptr_t baseAi = (uintptr_t) adj_in.data();
  const size_t szf = sizeof(NNScalar);
  for (int c = 0; c < nc; c++) {
    h_redptr(c)          = baseAo + (uintptr_t) (strideAo * hr.tsi[c]) * szf;
    h_redptr(nc + c)     = baseW  + (uintptr_t) ((size_t) c * n_out * n_in) * szf;
    h_redptr(2 * nc + c) = baseAi + (uintptr_t) (strideAi * hr.ci[c]) * szf;
  }
  Kokkos::deep_copy(d_redptr, h_redptr);
  const float** Aarr = (const float**) d_redptr.data();
  const float** Barr = (const float**) (d_redptr.data() + nc);
  float** Carr       = (float**)       (d_redptr.data() + 2 * nc);
  const float alpha = 1.0f, beta = 1.0f;
  const int ldAo = (int) adj_out.extent(0), ldAi = (int) adj_in.extent(0);
  cublasStatus_t st = cublasSgemmBatched(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
      chunk_size, n_in, n_out, &alpha,
      Aarr, ldAo, Barr, n_out, &beta,
      Carr, ldAi, nc);
  if (st != CUBLAS_STATUS_SUCCESS)
    error->all(FLERR, "GRACE-3L/KK: reduceN bwd cublasSgemmBatched failed (status {})", (int) st);
  return true;
}
#endif

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_fc_bwd(
    const DeviceFC &f, const t_nn_3d &adj_out, const t_nn_3d &adj_left,
    const t_nn_3d &adj_right)
{
  fc_p = f;
  fc_adj_out = adj_out;
  fc_adj_left = adj_left;
  fc_adj_right = adj_right;
  const int dim1_left = f.left_coefs ? f.n_in_left : f.n_out;
  Kokkos::parallel_for("FCBwdLeft",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>, TagFCBwdLeft>(
          {0, 0, 0}, {chunk_size, dim1_left, f.n_funcs_left}), *this);
  Kokkos::parallel_for("FCBwdRight",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>, TagFCBwdRight>(
          {0, 0, 0}, {chunk_size, f.n_funcs_right, f.n_in_right}), *this);
}

// ======================================================================
// Kernel: TagCPLBwdFused (Opt 5) — fused backward cp_l passes 1+2. ONE TEAM per (a,r).
// Replaces the separate CPLBwdRecomputeProject (global lproj/rproj) + CPLBwdCouple
// (n_cg re-reads of the big adj_out[a,r,:] row + n_cg global RMW into the small adj
// rows). Here lproj/rproj are RECOMPUTED into shared, adj_out[a,r,:] is loaded ONCE into
// shared, and the couple-VJP accumulates lproj_adj/rproj_adj in shared (shared atomics),
// written to global ONCE. (Pass 3, TagCPLBwdProject, stays separate — it reduces over r.)
// Shared: [ lproj_sh(nlm) | rproj_sh(nlm) | lproj_adj_sh(nlm) | rproj_adj_sh(nlm) | adj_out_sh(nfunc) ].
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagCPLBwdFused,
    const typename Kokkos::TeamPolicy<DeviceType, TagCPLBwdFused>::member_type& team) const
{
  const int rank = cpl_p.rank;
  const int a = team.league_rank() / rank;
  const int r = team.league_rank() - a * rank;
  const int nlm = cpl_nlm;
  const int nfunc = cpl_p.nfunc;
  const int n_lm_left = (int) cpl_p.group_left.extent(0);
  const int n_lm_right = (int) cpl_p.group_right.extent(0);
  const int nL = cpl_p.n_left;
  const int nR = cpl_p.n_right;

  NNScalar* sh = (NNScalar*) team.team_shmem().get_shmem((4 * nlm + nfunc) * sizeof(NNScalar), 0);
  NNScalar* lproj_sh     = sh;
  NNScalar* rproj_sh     = sh + nlm;
  NNScalar* lproj_adj_sh = sh + 2 * nlm;
  NNScalar* rproj_adj_sh = sh + 3 * nlm;
  NNScalar* adj_out_sh   = sh + 4 * nlm;

  // Load the cuBLAS projections when available; otherwise recompute forward
  // Pass A into shared. The VJP below remains fused/shared in both cases.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlm), [&] (const int w) {
    NNScalar sl = NNScalar(0.0);
    if (w < n_lm_left) {
      if (cpl_bwd_projected) {
        sl = d_cpl_lproj(a, r, w);
      } else {
        const int gl = cpl_p.group_left(w);
        for (int nn = 0; nn < nL; nn++) sl += cpl_p.U(gl, r, nn) * cpl_in_left(a, nn, w);
        sl *= cpl_p.norm_u;
      }
    }
    lproj_sh[w] = sl; lproj_adj_sh[w] = NNScalar(0.0);
    NNScalar sr = NNScalar(0.0);
    if (w < n_lm_right) {
      if (cpl_bwd_projected) {
        sr = d_cpl_rproj(a, r, w);
      } else {
        const int gr = cpl_p.group_right(w);
        for (int nn = 0; nn < nR; nn++) sr += cpl_p.V(gr, r, nn) * cpl_in_right(a, nn, w);
        sr *= cpl_p.norm_v;
      }
    }
    rproj_sh[w] = sr; rproj_adj_sh[w] = NNScalar(0.0);
  });
  // Load the adj_out[a,r,:] row into shared once.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nfunc), [&] (const int fcn) {
    adj_out_sh[fcn] = cpl_adj_out(a, r, fcn);
  });
  team.team_barrier();

  // Couple VJP: adj_p = adj_out_sh[m_sum_ind[c]]*cg[c];
  //   lproj_adj_sh[left_ind[c]] += adj_p*rproj_sh[right_ind[c]];
  //   rproj_adj_sh[right_ind[c]] += adj_p*lproj_sh[left_ind[c]]  (shared atomics).
  // Opt 14: the three indices come packed in one int32 (8 B/entry with cg vs 16 B).
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, cpl_p.n_cg), [&] (const int c) {
    const int pk = cpl_p.cg_packed(c);
    const int li = pk & 255;
    const int ri = (pk >> 8) & 255;
    const NNScalar d = adj_out_sh[pk >> 16] * cpl_p.cg(c);
    Kokkos::atomic_add(&lproj_adj_sh[li], d * rproj_sh[ri]);
    Kokkos::atomic_add(&rproj_adj_sh[ri], d * lproj_sh[li]);
  });
  team.team_barrier();

  // Write the adj rank-projections to global ONCE (consumed by TagCPLBwdProject).
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlm), [&] (const int w) {
    d_cpl_lproj_adj(a, r, w) = lproj_adj_sh[w];
    d_cpl_rproj_adj(a, r, w) = rproj_adj_sh[w];
  });
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_cp_l_bwd(
    const DeviceCPL &p, const t_nn_3d &in_left, const t_nn_3d &in_right,
    const t_nn_3d &adj_out, const t_nn_3d &adj_left, const t_nn_3d &adj_right)
{
  cpl_p = p;
  cpl_in_left = in_left;
  cpl_in_right = in_right;
  cpl_adj_out = adj_out;
  cpl_adj_left = adj_left;
  cpl_adj_right = adj_right;
  const int n_lm_left = (int) p.group_left.extent(0);
  const int n_lm_right = (int) p.group_right.extent(0);
  cpl_nlm = (n_lm_left > n_lm_right) ? n_lm_left : n_lm_right;
  cpl_bwd_nmax = (p.n_left > p.n_right) ? p.n_left : p.n_right;
  cpl_bwd_projected = 0;

  // 1+2) fused recompute-project + couple-VJP (team per (a,r); lproj/rproj + adj_out row
  //      + adj accumulators in shared) -> d_cpl_lproj_adj/d_cpl_rproj_adj.
#ifdef KOKKOS_ENABLE_CUDA
  // Use the same true-fp32 cuBLAS projections as the forward path, but keep the
  // VJP fused in shared memory. This avoids the hand projection's strided
  // LayoutLeft traversal without reintroducing global scattered VJP updates.
  if (cublas_handle && sizeof(NNScalar) == 4) {
    compute_cp_l_project_cublas(p, in_left, in_right);
    cpl_bwd_projected = 1;
  }
#endif
  {
    int team_size = 1, vector_length = 1;
    if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 128;
    const int league = chunk_size * p.rank;
    int scratch_size = scratch_size_helper<NNScalar>(4 * cpl_nlm + p.nfunc);
    check_team_size_for<TagCPLBwdFused>(league, team_size, vector_length);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagCPLBwdFused>(league, team_size, vector_length);
    policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for("CPLBwdFused", policy, *this);
  }
  // 3) project bwd -> adj_left/adj_right (accumulate; reduces over r, stays separate).
#ifdef KOKKOS_ENABLE_CUDA
  // Opt 16: two batched col-major GEMMs (left with Uw, right with Vw), each beta=1
  // into the shared adj buffer. The LayoutLeft slabs A=lproj_adj(:,:,w) / C=adj(:,:,w)
  // are col-major with ld=extent(0), w-stride=extent(0)*extent(1) -> no scatter needed.
  // For SELF-products cpl_adj_left aliases cpl_adj_right; the two GEMMs run sequentially
  // on the SAME (default-instance) stream, so the beta=1 accumulation is race-free.
  if (cublas_handle && sizeof(NNScalar) == 4) {
    const long strideA = (long) d_cpl_lproj_adj.extent(0) * d_cpl_lproj_adj.extent(1);
    const int  ldA     = (int)  d_cpl_lproj_adj.extent(0);
    cpl_project_sgemm_batched(
        d_cpl_lproj_adj.data(), strideA, ldA,
        p.Uw.data(), p.n_left, p.rank,
        cpl_adj_left.data(),
        (long) cpl_adj_left.extent(0) * cpl_adj_left.extent(1),
        (int)  cpl_adj_left.extent(0),
        chunk_size, n_lm_left, p.norm_u);
    cpl_project_sgemm_batched(
        d_cpl_rproj_adj.data(), strideA, ldA,
        p.Vw.data(), p.n_right, p.rank,
        cpl_adj_right.data(),
        (long) cpl_adj_right.extent(0) * cpl_adj_right.extent(1),
        (int)  cpl_adj_right.extent(0),
        chunk_size, n_lm_right, p.norm_v);
    return;
  }
#endif
  Kokkos::parallel_for("CPLBwdProject",
      Kokkos::RangePolicy<DeviceType, TagCPLBwdProject>(0, chunk_size * cpl_bwd_nmax * cpl_nlm), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_adj_prod(
    int spbf_idx, const t_nn_3d &A_adj)
{
  force_spbf_idx = spbf_idx;
  force_A_adj = A_adj;
  const int nradmax = d_spbf[spbf_idx].n_rad_max;
  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int P = nlm_y * d_spbf[spbf_idx].n_lm_ind;
  Kokkos::parallel_for("ComputeAdjProd",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>, TagComputeAdjProd>(
          {0, 0, 0}, {chunk_size, nradmax, P}), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_A3_spbf_bwd()
{
  // Opt 10: fused team+shared version (old TagA3SPBFBwd kept for reference).
  compute_spbf_bwd_fused(idx_spbf_A3, d_eq2_norm_adj);
}

// ---- TagA2SPBFBwd (Task 3.2): VJP of TagComputeA2 down to the L2 indicator
// eq1_norm. IDENTICAL structure to TagA3SPBFBwd except: reads d_A2_adj, uses the
// A2 SPBF weights (idx_spbf_A2, n_lm_ind=25), and scatters into d_grad_I_global
// [nall, I_n_funcs=64, I_n_out=25] (the round-1 reverse buffer = d_eq1_norm_adj,
// owned partial). Requires d_R1_nl to hold A2's radial (caller recomputes it).
// The z_proj/chem_l0_mask species-injection term is constant -> no indicator adj.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagA2SPBFBwd,
    const int ii, const int n) const
{
  const DeviceSPBF &sp = d_spbf[idx_spbf_A2];
  if (n >= sp.n_rad_max) return;

  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int n_lm_ind = sp.n_lm_ind;
  // Opt 4: adj_prod precomputed into d_adj_prod by compute_adj_prod(idx_spbf_A2).

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int j_global = d_nearest(ii, jj);
    const NNScalar env = (NNScalar) d_env(ii, jj);
    NNScalar a_nl[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int l = sp.l_tile(ly);
      a_nl[ly] = d_R1_nl(ii, jj, n, l) * (NNScalar) d_Y_bond(ii, jj, ly) * env;
    }
    for (int lr = 0; lr < n_lm_ind; lr++) {
      NNScalar s = NNScalar(0.0);
      for (int ly = 0; ly < nlm_y; ly++)
        s += a_nl[ly] * d_adj_prod(ii, n, ly * n_lm_ind + lr);
      Kokkos::atomic_add(&d_grad_I_global(j_global, n, lr), s);
    }
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_A2_spbf_bwd()
{
  // Opt 10: fused team+shared version (old TagA2SPBFBwd kept for reference).
  compute_spbf_bwd_fused(idx_spbf_A2, d_grad_I_global);
}

// ---- TagSPBFBwdFused (Opt 10): team+shared equiv-SPBF backward, shared by A2
// (target d_grad_I_global) and A3 (target d_eq2_norm_adj). The old per-(ii,n)
// THREAD kernels re-read the full d_adj_prod row (P = nlm_y*n_lm_ind, 5 KB for
// A3) from global once PER NEIGHBOUR (~nct times) and walked the jj/lr/ly loops
// serially (~nct*n_lm_ind*nlm_y MACs per thread). Here ONE TEAM owns one (ii,n):
// the adj_prod row is staged into TEAM SHARED once, neighbours are processed in
// tiles of SPBF_BWD_TILE with each tile's a_nl[tn,nlm_y] staged cooperatively,
// and the (jj,lr) scatter pairs are split across the team. The atomic scatter
// into the nall-sized target is unchanged (same count, same nondeterministic
// accumulation order as before); the inner ly sum keeps ascending order.
// Shared layout: [ adjp_sh(P) | anl_sh(SPBF_BWD_TILE*nlm_y) ] NNScalars.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagSPBFBwdFused,
    const typename Kokkos::TeamPolicy<DeviceType, TagSPBFBwdFused>::member_type& team) const
{
  const DeviceSPBF &sp = d_spbf[spbf_bwd_idx];
  const int ii = team.league_rank() % chunk_size;   // league = chunk_size * n_rad_max
  const int n  = team.league_rank() / chunk_size;
  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int nli = sp.n_lm_ind;
  const int P = nlm_y * nli;

  NNScalar* sh = (NNScalar*) team.team_shmem().get_shmem((P + (int) SPBF_BWD_TILE * nlm_y) * sizeof(NNScalar), 0);
  NNScalar* adjp_sh = sh;                 // [P] this (ii,n)'s adj_prod row
  NNScalar* anl_sh  = sh + P;             // [tn, nlm_y] staged tile of neighbour a_nl

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, P), [&] (const int p) {
    adjp_sh[p] = d_adj_prod(ii, n, p);
  });
  team.team_barrier();

  const int nct = d_ncount(ii);
  for (int t0 = 0; t0 < nct; t0 += (int) SPBF_BWD_TILE) {
    const int tn = (nct - t0 < (int) SPBF_BWD_TILE) ? (nct - t0) : (int) SPBF_BWD_TILE;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, tn * nlm_y), [&] (const int q) {
      const int j = q / nlm_y, ly = q - j * nlm_y;
      const int l = sp.l_tile(ly);
      anl_sh[q] = d_R1_nl(ii, t0 + j, n, l) * (NNScalar) d_Y_bond(ii, t0 + j, ly)
                  * (NNScalar) d_env(ii, t0 + j);
    });
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, tn * nli), [&] (const int q) {
      const int j = q / nli, lr = q - j * nli;
      const NNScalar* an = anl_sh + j * nlm_y;
      NNScalar s = NNScalar(0.0);
      for (int ly = 0; ly < nlm_y; ly++)
        s += an[ly] * adjp_sh[ly * nli + lr];
      Kokkos::atomic_add(&spbf_bwd_target(d_nearest(ii, t0 + j), n, lr), s);
    });
    team.team_barrier();
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_spbf_bwd_fused(
    int idx, const t_nn_3d &adj_target)
{
  spbf_bwd_idx = idx;
  spbf_bwd_target = adj_target;
  const int nradmax = d_spbf[idx].n_rad_max;
  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int P = nlm_y * d_spbf[idx].n_lm_ind;

  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 128;
  const int league = chunk_size * nradmax;
  int scratch_size = scratch_size_helper<NNScalar>(P + (int) SPBF_BWD_TILE * nlm_y);
  check_team_size_for<TagSPBFBwdFused>(league, team_size, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagSPBFBwdFused>(league, team_size, vector_length);
  policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
  Kokkos::parallel_for("SPBFBwdFused", policy, *this);
}

// ======================================================================
// Task 3.3: geometry backward (per-bond forces).
// ======================================================================

// ---- TagComputeMLPRadialDeriv: radial MLP forward + its r-derivative, one bond
// per thread, for the layer selected by radial_mlp_spbf. Reproduces
// TagComputeMLPRadial's value path EXACTLY (bit-identical R via the same raw
// Chebyshev basis read from d_radial_basis) and forward-propagates the r-derivative
// (only the Chebyshev inputs depend on r). Writes d_R1_nl, d_DR1_nl (=dR/dr) and
// the layer-independent envelope derivative d_denv (=denv/dr). silu'(s) =
// sig*(1+s*(1-sig)).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeMLPRadialDeriv,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadialDeriv>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  if (jj >= d_ncount(ii)) return;

  const DeviceSPBF &sp = d_spbf[radial_mlp_spbf];
  const int n_layers = sp.mlp_n_layers;
  int mlp_dims[GRACE3L_MAX_MLP_LAYERS + 1];
  NNScalar mlp_norms[GRACE3L_MAX_MLP_LAYERS];
  for (int i = 0; i <= n_layers; i++) mlp_dims[i] = sp.mlp_dims(i);
  for (int i = 0; i < n_layers; i++) mlp_norms[i] = sp.mlp_norms(i);

  // Geometry: bond length + per-element cutoff (== TagComputeRadialBasis).
  const int mu_i = d_mu_i(ii);
  const int mu_j = d_mu_j(ii, jj);
  const GeomScalar r = d_rnorms(ii, jj);
  const GeomScalar rcut_ij = d_bond_cutoff(mu_i, mu_j);
  const GeomScalar x_norm = r / rcut_ij;
  const GeomScalar x_cheb = GeomScalar(2.0) * x_norm - GeomScalar(1.0);
  const GeomScalar dxcheb_dr = GeomScalar(2.0) / rcut_ij;   // d x_cheb / dr

  // Envelope derivative denv/dr (env = 1 - 0.5*p1p2*x^p + pp2*x^{p+1} - 0.5*pp1*x^{p+2}).
  {
    const int p = radial_basis_p;
    GeomScalar xp = GeomScalar(1.0);
    for (int ip = 0; ip < p; ip++) xp *= x_norm;             // x^p
    const GeomScalar xpm1 = xp / x_norm;                     // x^{p-1} (r>0)
    const GeomScalar xp1  = xp * x_norm;                     // x^{p+1}
    const GeomScalar pp1  = GeomScalar(p) * (p + 1);
    const GeomScalar pp2  = GeomScalar(p) * (p + 2);
    const GeomScalar p1p2 = GeomScalar(p + 1) * (p + 2);
    const GeomScalar denv_dx = GeomScalar(-0.5) * p1p2 * GeomScalar(p) * xpm1
                             + pp2 * GeomScalar(p + 1) * xp
                             - GeomScalar(0.5) * pp1 * GeomScalar(p + 2) * xp1;
    d_denv(ii, jj) = denv_dx / rcut_ij;                      // d x_norm / dr = 1/rcut
  }

  NNScalar buf_a[GRACE3L_MAX_MLP_DIM], buf_b[GRACE3L_MAX_MLP_DIM];
  NNScalar dbuf_a[GRACE3L_MAX_MLP_DIM], dbuf_b[GRACE3L_MAX_MLP_DIM];
  NNScalar* h_cur = buf_a;  NNScalar* h_nxt = buf_b;
  NNScalar* dh_cur = dbuf_a; NNScalar* dh_nxt = dbuf_b;

  // Input value path: read the raw Chebyshev basis (bit-identical to the forward).
  for (int b = 0; b < nradbase; b++) h_cur[b] = (NNScalar) d_radial_basis(ii, jj, b);
  // Input derivative path for the Chebyshev basis: basis[b] = T_{b+1}(x_cheb),
  // dT_1/dx=1, dT_{k+1}/dx = 2*T_k + 2*x*dT_k - dT_{k-1}; * dxcheb_dr for d/dr.
  {
    GeomScalar dT_prev = GeomScalar(0.0);   // dT_0/dx
    GeomScalar dT_curr = GeomScalar(1.0);   // dT_1/dx
    dh_cur[0] = (NNScalar)(dT_curr * dxcheb_dr);
    for (int kk = 1; kk < nradbase; kk++) {
      const GeomScalar T_kk = (GeomScalar) d_radial_basis(ii, jj, kk - 1);  // T_kk
      const GeomScalar dT_next = GeomScalar(2.0) * T_kk
                               + GeomScalar(2.0) * x_cheb * dT_curr - dT_prev;
      dh_cur[kk] = (NNScalar)(dT_next * dxcheb_dr);
      dT_prev = dT_curr; dT_curr = dT_next;
    }
  }
  int k = nradbase;
  for (int e = 0; e < embedding_size; e++) { h_cur[k] = d_chem_embed(mu_i, e); dh_cur[k] = NNScalar(0.0); k++; }
  for (int e = 0; e < embedding_size; e++) { h_cur[k] = d_chem_embed(mu_j, e); dh_cur[k] = NNScalar(0.0); k++; }

  for (int layer = 0; layer < n_layers - 1; layer++) {
    const int nin = mlp_dims[layer], nout = mlp_dims[layer + 1];
    const NNScalar norm = mlp_norms[layer];
    for (int j = 0; j < nout; j++) {
      NNScalar s = 0.0, ds = 0.0;
      for (int kk = 0; kk < nin; kk++) {
        const NNScalar Wv = sp.mlp_W(layer, kk, j);
        s  += Wv * h_cur[kk];
        ds += Wv * dh_cur[kk];
      }
      s  = s * norm + sp.mlp_b(layer, j) * norm;
      ds = ds * norm;
      const NNScalar sig = NNScalar(1.0) / (NNScalar(1.0) + Kokkos::exp(-s));
      h_nxt[j]  = s * sig;                                          // silu(s)
      dh_nxt[j] = sig * (NNScalar(1.0) + s * (NNScalar(1.0) - sig)) * ds;  // silu'(s)*ds
    }
    NNScalar* t; t = h_cur; h_cur = h_nxt; h_nxt = t;
    t = dh_cur; dh_cur = dh_nxt; dh_nxt = t;
  }

  const int last_layer = n_layers - 1;
  const int n_last_hidden = mlp_dims[last_layer];
  const NNScalar out_norm = mlp_norms[last_layer];
  const int lmax1 = lmax + 1;
  for (int n = 0; n < sp.n_rad_max; n++) {
    for (int l = 0; l < lmax1; l++) {
      const int j = n * lmax1 + l;
      NNScalar R = 0.0, DR = 0.0;
      for (int kk = 0; kk < n_last_hidden; kk++) {
        const NNScalar Wv = sp.mlp_W(last_layer, kk, j);
        R  += Wv * h_cur[kk];
        DR += Wv * dh_cur[kk];
      }
      d_R1_nl(ii, jj, n, l)  = R * out_norm;
      d_DR1_nl(ii, jj, n, l) = DR * out_norm;
    }
  }
}

#ifdef KOKKOS_ENABLE_CUDA
// ======================================================================
// Opt 15 (value path): support kernels for the batched cuBLAS radial MLP.
// TagMLPAssembleVal builds the tightly-packed row-major input matrix
// d_mlp_X[M x n_in0], one bond per thread. Bond b decodes to
// (ii = b / bond_stride, jj = b % bond_stride) with bond_stride =
// d_radial_basis.extent(1) (== allocated maxneigh). Padding rows (jj >=
// ncount) are zeroed so the batched SGEMM stays finite; their outputs are
// ignored downstream (consumers loop jj < ncount). Input order is pinned to
// the oracle: [Chebyshev basis, chem_embed(mu_i), chem_embed(mu_j)].
// TagMLPActVal applies the per-layer affine (s = raw*norm + bias*norm) and
// silu in place on d_mlp_h{0,1} (raw = X*W from the GEMM, no norm folded).
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagMLPAssembleVal,
    const int b) const
{
  const int bond_stride = (int) d_radial_basis.extent(1);
  const int ii = b / bond_stride;
  const int jj = b - ii * bond_stride;
  const int n_in0 = nradbase + 2 * embedding_size;
  if (jj >= d_ncount(ii)) {
    for (int c = 0; c < n_in0; c++) d_mlp_X(b, c) = NNScalar(0.0);
    return;
  }
  const int mu_i = d_mu_i(ii);
  const int mu_j = d_mu_j(ii, jj);
  int k = 0;
  for (int c = 0; c < nradbase; c++)     d_mlp_X(b, k++) = (NNScalar) d_radial_basis(ii, jj, c);
  for (int e = 0; e < embedding_size; e++) d_mlp_X(b, k++) = d_chem_embed(mu_i, e);
  for (int e = 0; e < embedding_size; e++) d_mlp_X(b, k++) = d_chem_embed(mu_j, e);
}

// TagMLPAssembleDeriv (Opt 15b): value input d_mlp_X[b, :n_in0] (identical to
// AssembleVal) + deriv input d_mlp_dX[b, :nradbase] (Chebyshev dT/dr recurrence;
// embeddings are r-independent so excluded, giving the K=nradbase input sub-GEMM)
// + the per-bond envelope derivative d_denv (layer-independent). Reproduces the
// hand kernel TagComputeMLPRadialDeriv's input/envelope math EXACTLY. Padding
// bonds (jj>=ncount) zero X and dX (finite GEMM) and skip d_denv (unused).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagMLPAssembleDeriv,
    const int b) const
{
  const int bond_stride = (int) d_radial_basis.extent(1);
  const int ii = b / bond_stride;
  const int jj = b - ii * bond_stride;
  const int n_in0 = nradbase + 2 * embedding_size;
  if (jj >= d_ncount(ii)) {
    for (int c = 0; c < n_in0; c++)   d_mlp_X(b, c)  = NNScalar(0.0);
    for (int c = 0; c < nradbase; c++) d_mlp_dX(b, c) = NNScalar(0.0);
    return;
  }
  const int mu_i = d_mu_i(ii);
  const int mu_j = d_mu_j(ii, jj);

  // Geometry (== TagComputeRadialBasis / hand deriv kernel).
  const GeomScalar r = d_rnorms(ii, jj);
  const GeomScalar rcut_ij = d_bond_cutoff(mu_i, mu_j);
  const GeomScalar x_norm = r / rcut_ij;
  const GeomScalar x_cheb = GeomScalar(2.0) * x_norm - GeomScalar(1.0);
  const GeomScalar dxcheb_dr = GeomScalar(2.0) / rcut_ij;

  // Envelope derivative denv/dr.
  {
    const int p = radial_basis_p;
    GeomScalar xp = GeomScalar(1.0);
    for (int ip = 0; ip < p; ip++) xp *= x_norm;
    const GeomScalar xpm1 = xp / x_norm;
    const GeomScalar xp1  = xp * x_norm;
    const GeomScalar pp1  = GeomScalar(p) * (p + 1);
    const GeomScalar pp2  = GeomScalar(p) * (p + 2);
    const GeomScalar p1p2 = GeomScalar(p + 1) * (p + 2);
    const GeomScalar denv_dx = GeomScalar(-0.5) * p1p2 * GeomScalar(p) * xpm1
                             + pp2 * GeomScalar(p + 1) * xp
                             - GeomScalar(0.5) * pp1 * GeomScalar(p + 2) * xp1;
    d_denv(ii, jj) = denv_dx / rcut_ij;
  }

  // Value input: [Chebyshev basis, chem_embed(mu_i), chem_embed(mu_j)].
  int k = 0;
  for (int c = 0; c < nradbase; c++)     d_mlp_X(b, k++) = (NNScalar) d_radial_basis(ii, jj, c);
  for (int e = 0; e < embedding_size; e++) d_mlp_X(b, k++) = d_chem_embed(mu_i, e);
  for (int e = 0; e < embedding_size; e++) d_mlp_X(b, k++) = d_chem_embed(mu_j, e);

  // Deriv input: dT_{b+1}/dr for the Chebyshev basis (dT_1/dx=1;
  // dT_{k+1}/dx = 2*T_k + 2*x*dT_k - dT_{k-1}); *dxcheb_dr for d/dr.
  {
    GeomScalar dT_prev = GeomScalar(0.0);   // dT_0/dx
    GeomScalar dT_curr = GeomScalar(1.0);   // dT_1/dx
    d_mlp_dX(b, 0) = (NNScalar)(dT_curr * dxcheb_dr);
    for (int kk = 1; kk < nradbase; kk++) {
      const GeomScalar T_kk = (GeomScalar) d_radial_basis(ii, jj, kk - 1);  // T_kk
      const GeomScalar dT_next = GeomScalar(2.0) * T_kk
                               + GeomScalar(2.0) * x_cheb * dT_curr - dT_prev;
      d_mlp_dX(b, kk) = (NNScalar)(dT_next * dxcheb_dr);
      dT_prev = dT_curr; dT_curr = dT_next;
    }
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagMLPActVal,
    const int idx) const
{
  const int nout = mlp_act_n;
  const int b = idx / nout;
  const int j = idx - b * nout;
  const DeviceSPBF &sp = d_spbf[mlp_assemble_idx];
  const NNScalar norm = mlp_act_norm;
  if (mlp_act_layer == 0) {
    const NNScalar s = d_mlp_h0(b, j) * norm + sp.mlp_b(0, j) * norm;
    d_mlp_h0(b, j) = silu(s);
  } else {
    const NNScalar s = d_mlp_h1(b, j) * norm + sp.mlp_b(1, j) * norm;
    d_mlp_h1(b, j) = silu(s);
  }
}

// TagMLPActDeriv (Opt 15b): fused value+deriv activation. Reads the raw value
// GEMM output (d_mlp_h*) and raw deriv GEMM output (d_mlp_dh*), forms
// s = raw*norm + bias*norm and ds = raw_deriv*norm, then writes silu(s) and
// silu'(s)*ds in place. silu'(s) = sig*(1 + s*(1-sig)) -- identical to the hand
// kernel TagComputeMLPRadialDeriv.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagMLPActDeriv,
    const int idx) const
{
  const int nout = mlp_act_n;
  const int b = idx / nout;
  const int j = idx - b * nout;
  const DeviceSPBF &sp = d_spbf[mlp_assemble_idx];
  const NNScalar norm = mlp_act_norm;
  if (mlp_act_layer == 0) {
    const NNScalar s  = d_mlp_h0(b, j) * norm + sp.mlp_b(0, j) * norm;
    const NNScalar ds = d_mlp_dh0(b, j) * norm;
    const NNScalar sig = NNScalar(1.0) / (NNScalar(1.0) + Kokkos::exp(-s));
    d_mlp_h0(b, j)  = s * sig;
    d_mlp_dh0(b, j) = sig * (NNScalar(1.0) + s * (NNScalar(1.0) - sig)) * ds;
  } else {
    const NNScalar s  = d_mlp_h1(b, j) * norm + sp.mlp_b(1, j) * norm;
    const NNScalar ds = d_mlp_dh1(b, j) * norm;
    const NNScalar sig = NNScalar(1.0) / (NNScalar(1.0) + Kokkos::exp(-s));
    d_mlp_h1(b, j)  = s * sig;
    d_mlp_dh1(b, j) = sig * (NNScalar(1.0) + s * (NNScalar(1.0) - sig)) * ds;
  }
}

// TagMLPScatterR: write the LayoutRight output-GEMM result into the (LayoutLeft)
// 4D target via the View operator() (layout-agnostic). mlp_scatter_which selects
// value (d_mlp_R -> d_R1_nl) vs deriv (d_mlp_dR -> d_DR1_nl). Padding bonds
// (jj >= ncount) are skipped so the untouched entries match the hand kernel.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagMLPScatterR,
    const int idx) const
{
  const int lm1 = lmax + 1;
  const int D3  = d_spbf[mlp_assemble_idx].n_rad_max * lm1;
  const int bs  = (int) d_radial_basis.extent(1);
  const int b   = idx / D3;
  const int j   = idx - b * D3;
  const int ii  = b / bs;
  const int jj  = b - ii * bs;
  if (jj >= d_ncount(ii)) return;
  if (mlp_scatter_which == 0) d_R1_nl(ii, jj, j / lm1, j % lm1)  = d_mlp_R(b, j);
  else                        d_DR1_nl(ii, jj, j / lm1, j % lm1) = d_mlp_dR(b, j);
}
#endif

// ---- bond_geom_force: shared per-bond force chain (spec §4.1). Adapted VERBATIM
// from the validated GRACE-2L ComputeDerivative kernels (Y00==1 plm/dplm
// recurrence, tangential-gradient DY projection, f += w*(Y*DR*rhat + DY*R/r)),
// with the 3L twist that the envelope is a SEPARATE factor: Renv = R*env,
// DRenv = dR/dr*env + R*denv/dr. Geometry math is GeomScalar (fp64).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::bond_geom_force(
    const int ii, const int jj, const int n, const NNScalar* w,
    GeomScalar& f0, GeomScalar& f1, GeomScalar& f2) const
{
  const GeomScalar rinv = GeomScalar(1.0) / d_rnorms(ii, jj);
  const GeomScalar rx = d_rhats(ii, jj, 0);
  const GeomScalar ry = d_rhats(ii, jj, 1);
  const GeomScalar rz = d_rhats(ii, jj, 2);
  const GeomScalar env  = d_env(ii, jj);
  const GeomScalar denv = d_denv(ii, jj);

  const GeomScalar Y00_v = GeomScalar(1.0), sq3_v = sq3, sq3o2_v = sq3o2, sq2_v = sq2;

  GeomScalar plm[15];
  GeomScalar dplm[15];
  plm[0] = Y00_v; dplm[0] = GeomScalar(0.0);
  if (lmax > 0) {
    plm[1] = Y00_v * sq3_v * rz; dplm[1] = Y00_v * sq3_v;
    plm[2] = -sq3o2_v * Y00_v;   dplm[2] = GeomScalar(0.0);
    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        const int idx = l*(l+1)/2 + m;
        const int i1 = (l-1)*l/2 + m, i2 = (l-2)*(l-1)/2 + m;
        const int ai = d_idx_sph(l*(l+1) + m);
        const GeomScalar a = alm(ai), b = blm(ai);
        plm[idx]  = a * (rz * plm[i1] + b * plm[i2]);
        dplm[idx] = a * (plm[i1] + rz * dplm[i1] + b * dplm[i2]);
      }
      { const int idx = l*(l+1)/2+l-1, prev = (l-1)*l/2+l-1;
        const GeomScalar t = dl(l) * plm[prev];
        plm[idx] = t * rz; dplm[idx] = t; }
      { const int idx = l*(l+1)/2+l, prev = (l-1)*l/2+l-1;
        plm[idx] = cl(l) * plm[prev]; dplm[idx] = GeomScalar(0.0); }
    }
  }

  const GeomScalar phase_re = rx, phase_im = ry;
  for (int l = 0; l <= lmax; l++) {
    const GeomScalar R_raw  = (GeomScalar) d_R1_nl(ii, jj, n, l);
    const GeomScalar DR_raw = (GeomScalar) d_DR1_nl(ii, jj, n, l);
    const GeomScalar Renv  = R_raw * env;
    const GeomScalar DRenv = DR_raw * env + R_raw * denv;   // d(R*env)/dr
    const GeomScalar R_over_r = Renv * rinv;

    // m = 0
    {
      const GeomScalar Y = plm[l*(l+1)/2];
      const GeomScalar dp = dplm[l*(l+1)/2];
      const GeomScalar rdy = dp * rz;
      const GeomScalar DY_x = -rdy * rx, DY_y = -rdy * ry, DY_z = dp - rdy * rz;
      const GeomScalar wv = (GeomScalar) w[l*(l+1)];
      const GeomScalar YDR = Y * DRenv;
      f0 += wv * (YDR * rx + DY_x * R_over_r);
      f1 += wv * (YDR * ry + DY_y * R_over_r);
      f2 += wv * (YDR * rz + DY_z * R_over_r);
    }

    // m >= 1
    GeomScalar pm_re = phase_re, pm_im = phase_im;
    for (int m = 1; m <= l; m++) {
      const int fac = (m % 2 == 0) ? 1 : -1;
      const GeomScalar pv = plm[l*(l+1)/2+m], dpv = dplm[l*(l+1)/2+m];
      const GeomScalar ylm_re = pm_re * pv, ylm_im = pm_im * pv;
      const GeomScalar rYp = sq2_v * fac * ylm_re, rYn = sq2_v * fac * ylm_im;

      GeomScalar dyx_re, dyx_im, dyy_re, dyy_im;
      const GeomScalar dyz_re = dpv * pm_re, dyz_im = dpv * pm_im;
      if (m == 1) {
        dyx_re = pv; dyx_im = GeomScalar(0.0); dyy_re = GeomScalar(0.0); dyy_im = pv;
      } else {
        const GeomScalar s2 = rx*rx + ry*ry;
        if (s2 > GeomScalar(1e-14)) {
          const GeomScalar is2 = GeomScalar(1.0)/s2;
          const GeomScalar p1r = (pm_re*rx + pm_im*ry)*is2;
          const GeomScalar p1i = (pm_im*rx - pm_re*ry)*is2;
          const GeomScalar mp = GeomScalar(m) * pv;
          dyx_re = mp*p1r; dyx_im = mp*p1i;
          dyy_re = -dyx_im; dyy_im = dyx_re;
        } else { dyx_re = dyx_im = dyy_re = dyy_im = GeomScalar(0.0); }
      }
      const GeomScalar rdy_re = rx*dyx_re + ry*dyy_re + rz*dyz_re;
      const GeomScalar rdy_im = rx*dyx_im + ry*dyy_im + rz*dyz_im;
      const GeomScalar Dpx_re = dyx_re - rdy_re*rx, Dpx_im = dyx_im - rdy_im*rx;
      const GeomScalar Dpy_re = dyy_re - rdy_re*ry, Dpy_im = dyy_im - rdy_im*ry;
      const GeomScalar Dpz_re = dyz_re - rdy_re*rz, Dpz_im = dyz_im - rdy_im*rz;
      const GeomScalar DYpx = sq2_v*fac*Dpx_re, DYpy = sq2_v*fac*Dpy_re, DYpz = sq2_v*fac*Dpz_re;
      const GeomScalar DYnx = sq2_v*fac*Dpx_im, DYny = sq2_v*fac*Dpy_im, DYnz = sq2_v*fac*Dpz_im;

      const int Ap = l*(l+1)+m, An = l*(l+1)-m;
      const GeomScalar wp = (GeomScalar) w[Ap];
      const GeomScalar wn = (GeomScalar) w[An];

      { const GeomScalar YDR = rYp * DRenv;
        f0 += wp * (YDR*rx + DYpx*R_over_r);
        f1 += wp * (YDR*ry + DYpy*R_over_r);
        f2 += wp * (YDR*rz + DYpz*R_over_r); }
      { const GeomScalar YDR = rYn * DRenv;
        f0 += wn * (YDR*rx + DYnx*R_over_r);
        f1 += wn * (YDR*ry + DYny*R_over_r);
        f2 += wn * (YDR*rz + DYnz*R_over_r); }

      const GeomScalar t_re = pm_re*phase_re - pm_im*phase_im;
      const GeomScalar t_im = pm_re*phase_im + pm_im*phase_re;
      pm_re = t_re; pm_im = t_im;
    }
  }
}

// ---- TagForceL1: L1 scalar-SPBF geometry backward. One thread per (owned atom
// ii, radial n); loops neighbors jj. Per-lm weight w[lm] = d_A1_adj*z_tr*inv_avg
// (a_nl = R*Y*z_tr*env for L1). Accumulates into d_f_ij (atomic; multiple n
// threads write the same bond).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagForceL1,
    const int ii, const int n) const
{
  const DeviceSPBF &sp = d_spbf[0];   // A1
  if (n >= sp.n_rad_max) return;
  const int nlm = (lmax + 1) * (lmax + 1);
  const NNScalar inv_avg = (NNScalar) sp.inv_avg_n_neigh;

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const NNScalar z_n = sp.z_tr(mu_j, n);
    NNScalar w[25];
    for (int lm = 0; lm < nlm; lm++)
      w[lm] = d_A1_adj(ii, n, lm) * z_n * inv_avg;
    GeomScalar f0 = GeomScalar(0.0), f1 = GeomScalar(0.0), f2 = GeomScalar(0.0);
    bond_geom_force(ii, jj, n, w, f0, f1, f2);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 0), f0);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 1), f1);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 2), f2);
  }
}

// ---- TagForceEquiv: L2/L3 equivariant-SPBF geometry backward. One thread per
// (owned atom ii, radial n); loops neighbors jj. adj_prod[p] = inv_avg*Σ_f
// force_A_adj*cg_W (same as compute_A{2,3}_spbf_bwd), the per-lm weight is
// w[ly] = Σ_lr adj_prod[ly*n_lm_ind+lr]*bond_I[lr] with bond_I the FULL indicator
// (indicator_global[j] + z_proj injection — a_nl multiplies both). Accumulates
// into d_f_ij (atomic).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagForceEquiv,
    const int ii, const int n) const
{
  const DeviceSPBF &sp = d_spbf[force_spbf_idx];
  if (n >= sp.n_rad_max) return;

  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int n_lm_ind = force_n_lm_ind;
  // Opt 4: adj_prod precomputed into d_adj_prod by compute_adj_prod(force_spbf_idx).

  const int nct = d_ncount(ii);
  for (int jj = 0; jj < nct; jj++) {
    const int j_global = d_nearest(ii, jj);
    const int mu_j = d_mu_j(ii, jj);
    const NNScalar zpn = sp.z_proj(mu_j, n);

    NNScalar bond_I[GRACE3L_MAX_SPBF_NLMIND];
    if (force_ind_stage == 2) {
      for (int lr = 0; lr < n_lm_ind; lr++)
        bond_I[lr] = d_eq2_norm_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);
    } else {
      for (int lr = 0; lr < n_lm_ind; lr++)
        bond_I[lr] = d_I_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);
    }

    NNScalar w[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int base = ly * n_lm_ind;
      NNScalar s = NNScalar(0.0);
      for (int lr = 0; lr < n_lm_ind; lr++)
        s += d_adj_prod(ii, n, base + lr) * bond_I[lr];
      w[ly] = s;
    }
    GeomScalar f0 = GeomScalar(0.0), f1 = GeomScalar(0.0), f2 = GeomScalar(0.0);
    bond_geom_force(ii, jj, n, w, f0, f1, f2);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 0), f0);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 1), f1);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 2), f2);
  }
}

// ---- TagForceEquivFused (Opt 11): team-per-(ii,n) version of TagForceEquiv.
// The old kernel gave one THREAD the whole (ii,n) row and re-read the full
// d_adj_prod row (P = nlm_y*n_lm_ind, 5 KB for A3) from global once per
// neighbour (~nct times). Here ONE TEAM owns one (ii,n): the adj_prod row is
// staged into TEAM SHARED once, and the bonds are split across the team (one
// bond per thread, strided). The per-bond body — bond_I gather, w[ly] couple
// (ascending lr, unchanged accumulation order), bond_geom_force, d_f_ij atomic
// — is VERBATIM from TagForceEquiv, so per-bond arithmetic is identical.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagForceEquivFused,
    const typename Kokkos::TeamPolicy<DeviceType, TagForceEquivFused>::member_type& team) const
{
  const DeviceSPBF &sp = d_spbf[force_spbf_idx];
  const int ii = team.league_rank() % chunk_size;   // league = chunk_size * n_rad_max
  const int n  = team.league_rank() / chunk_size;

  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int n_lm_ind = force_n_lm_ind;
  const int P = nlm_y * n_lm_ind;

  NNScalar* adjp_sh = (NNScalar*) team.team_shmem().get_shmem(P * sizeof(NNScalar), 0);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, P), [&] (const int p) {
    adjp_sh[p] = d_adj_prod(ii, n, p);
  });
  team.team_barrier();

  const int nct = d_ncount(ii);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nct), [&] (const int jj) {
    const int j_global = d_nearest(ii, jj);
    const int mu_j = d_mu_j(ii, jj);
    const NNScalar zpn = sp.z_proj(mu_j, n);

    NNScalar bond_I[GRACE3L_MAX_SPBF_NLMIND];
    if (force_ind_stage == 2) {
      for (int lr = 0; lr < n_lm_ind; lr++)
        bond_I[lr] = d_eq2_norm_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);
    } else {
      for (int lr = 0; lr < n_lm_ind; lr++)
        bond_I[lr] = d_I_global(j_global, n, lr) + zpn * sp.chem_l0_mask(lr);
    }

    NNScalar w[25];
    for (int ly = 0; ly < nlm_y; ly++) {
      const int base = ly * n_lm_ind;
      NNScalar s = NNScalar(0.0);
      for (int lr = 0; lr < n_lm_ind; lr++)
        s += adjp_sh[base + lr] * bond_I[lr];
      w[ly] = s;
    }
    GeomScalar f0 = GeomScalar(0.0), f1 = GeomScalar(0.0), f2 = GeomScalar(0.0);
    bond_geom_force(ii, jj, n, w, f0, f1, f2);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 0), f0);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 1), f1);
    Kokkos::atomic_add(&d_f_ij(ii, jj, 2), f2);
  });
}

// ======================================================================
// Opt 15: value-only radial MLP forward for layer `radial_mlp_spbf`, writing
// d_R1_nl. On CUDA with a live handle and the expected 3-layer architecture,
// runs the batched true-fp32 cuBLAS path (assemble -> SGEMM/act x2 -> output
// SGEMM straight into d_R1_nl); otherwise falls back to the hand kernel
// TagComputeMLPRadial. Row-major operands throughout (see mlp_sgemm).
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_mlp_radial()
{
#ifdef KOKKOS_ENABLE_CUDA
  const auto &hmlp = grace_model->spbf[radial_mlp_spbf].mlp;
  if (cublas_handle && hmlp.n_layers == 3 && sizeof(NNScalar) == 4) {
    const int si = radial_mlp_spbf;
    const int bond_stride = (int) d_R1_nl.extent(1);   // == d_radial_basis.extent(1)
    const int M = chunk_size * bond_stride;
    mlp_M = M;
    const int n_in0 = hmlp.n_in[0];        // 138
    const int d1    = hmlp.n_out[0];       // 256
    const int d2    = hmlp.n_out[1];       // 128
    const int d3    = hmlp.n_out[2];       // n_rad_max*(lmax+1)

    // NOTE: no fences between these steps. The cuBLAS handle is bound (in
    // init_style) to DeviceType().cuda_stream() -- the SAME stream the Kokkos
    // parallel_fors below launch on -- so assemble -> GEMM -> act -> GEMM -> ...
    // -> scatter execute in strict issue order on one stream (race-free).

    // 1. assemble X[M x n_in0]
    mlp_assemble_idx = si;
    Kokkos::parallel_for("MLPAssembleVal",
        Kokkos::RangePolicy<DeviceType, TagMLPAssembleVal>(0, M), *this);

    // 2. h0 = X * W0                              (raw; norm+bias in Act)
    mlp_sgemm(d_mlp_Wg[si][0].data(), d_mlp_X.data(), d_mlp_h0.data(),
        d1, M, n_in0, (int)d_mlp_Wg[si][0].extent(1), (int)d_mlp_X.extent(1),
        (int)d_mlp_h0.extent(1), NNScalar(1.0));
    // 3. Act0: silu(h0*norm0 + b0*norm0)
    mlp_act_layer = 0; mlp_act_n = d1; mlp_act_norm = (NNScalar) hmlp.norm[0];
    Kokkos::parallel_for("MLPActVal0",
        Kokkos::RangePolicy<DeviceType, TagMLPActVal>(0, M * d1), *this);

    // 4. h1 = h0 * W1
    mlp_sgemm(d_mlp_Wg[si][1].data(), d_mlp_h0.data(), d_mlp_h1.data(),
        d2, M, d1, (int)d_mlp_Wg[si][1].extent(1), (int)d_mlp_h0.extent(1),
        (int)d_mlp_h1.extent(1), NNScalar(1.0));
    // 5. Act1: silu(h1*norm1 + b1*norm1)
    mlp_act_layer = 1; mlp_act_n = d2; mlp_act_norm = (NNScalar) hmlp.norm[1];
    Kokkos::parallel_for("MLPActVal1",
        Kokkos::RangePolicy<DeviceType, TagMLPActVal>(0, M * d2), *this);

    // 6. R = out_norm * (h1 * W2) into the LayoutRight scratch d_mlp_R
    //    (output layer: linear, no bias); scattered into d_R1_nl below.
    mlp_sgemm(d_mlp_Wg[si][2].data(), d_mlp_h1.data(), d_mlp_R.data(),
        d3, M, d2, (int)d_mlp_Wg[si][2].extent(1), (int)d_mlp_h1.extent(1),
        (int)d_mlp_R.extent(1), (NNScalar) hmlp.norm[2]);
    // d_R1_nl is a plain t_nn_4d (LayoutLeft on CUDA), so the (n,l) block of a
    // bond is NOT contiguous -> cuBLAS cannot write it directly. Scatter the
    // LayoutRight GEMM output d_mlp_R[b,j] into d_R1_nl(ii,jj, j/(lmax+1),
    // j%(lmax+1)) with a layout-agnostic kernel (skip padding bonds to match
    // the hand kernel, which returns early for jj>=ncount).
    mlp_scatter_which = 0;   // R -> d_R1_nl
    Kokkos::parallel_for("MLPScatterR",
        Kokkos::RangePolicy<DeviceType, TagMLPScatterR>(0, M * d3), *this);
    return;
  }
#endif
  // Fallback: original hand kernel.
  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 32;
  int ts = team_size;
  check_team_size_for<TagComputeMLPRadial>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>(
      ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
  Kokkos::parallel_for("ComputeMLPRadial", policy, *this);
}

#ifdef KOKKOS_ENABLE_CUDA
// Opt 15: true-fp32 batched SGEMM. Row-major C[M x n_out] = alpha*In[M x n_in]*Wg[n_in x n_out].
// cuBLAS is column-major, so this equals the col-major product
// C^T[n_out x M] = alpha * Wg^T[n_out x n_in] * In^T[n_in x M], i.e. OP_N/OP_N with
// leading dims = the row-major operands' physical row strides. CUBLAS_PEDANTIC_MATH
// (set on the handle) forces true fp32 (no tf32).
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::mlp_sgemm(
    const NNScalar* Wg, const NNScalar* In, NNScalar* C,
    int n_out, int M, int n_in, int ldW, int ldIn, int ldC, NNScalar alpha)
{
  const float a = (float) alpha, beta = 0.0f;
  cublasStatus_t st = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
      n_out, M, n_in, &a,
      (const float*) Wg, ldW,
      (const float*) In, ldIn,
      &beta, (float*) C, ldC);
  if (st != CUBLAS_STATUS_SUCCESS)
    error->all(FLERR, "GRACE-3L/KK: cublasSgemm failed (status {})", (int) st);
}

// Opt 16: batched projection SGEMM for CPLBwdProject. Standard col-major batched GEMM
// (NOT the row-major transpose trick): the LayoutLeft slabs A=lproj_adj(:,:,w) and
// C=adj(:,:,w) are already col-major [M x rank] / [M x n] with leading dim = extent(0)
// and per-batch (w) stride = extent(0)*extent(1). Bw (Uw/Vw) is packed col-major
// [rank x n], ldb=rank, w-stride=n*rank. beta=1 accumulates into the shared adj buffer
// (which already holds contributions from earlier backward stages). PEDANTIC = true fp32.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::cpl_project_sgemm_batched(
    const NNScalar* A, long strideA, int ldA,
    const NNScalar* Bw, int n, int rank, NNScalar* C, long strideC, int ldC,
    int M, int batch, NNScalar alpha)
{
  if (batch <= 0 || n <= 0 || M <= 0 || rank <= 0) return;
  const float a = (float) alpha, beta = 1.0f;
  cublasStatus_t st = cublasSgemmStridedBatched(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
      M, n, rank, &a,
      (const float*) A,  ldA,  (long long) strideA,
      (const float*) Bw, rank, (long long) n * rank,
      &beta, (float*) C, ldC,  (long long) strideC,
      batch);
  if (st != CUBLAS_STATUS_SUCCESS)
    error->all(FLERR, "GRACE-3L/KK: cublasSgemmStridedBatched failed (status {})", (int) st);
}

// Opt 15 hard gate: verify the layer-0 SGEMM (transpose/leading-dim convention +
// true-fp32 path) reproduces a host fp64 reference on a tiny deterministic batch
// before the full chain is ever used. Aborts on gross mismatch (transpose bug ->
// O(1) error); fp32 rounding over K~138 stays well under the 1e-4 gate.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::mlp_cublas_selftest()
{
  const auto &mlp = grace_model->spbf[0].mlp;   // A1
  if (mlp.n_layers < 1) return;
  const int n_in = mlp.n_in[0], n_out = mlp.n_out[0];
  const int Mt = 4;
  t_nn_2d_r X("g3l:selftest_X", Mt, n_in), C("g3l:selftest_C", Mt, n_out);
  auto hX = Kokkos::create_mirror_view(X);
  for (int b = 0; b < Mt; b++)
    for (int k = 0; k < n_in; k++)
      hX(b, k) = (NNScalar)(0.01 * (double)(((b * 7 + k * 3) % 13) - 6));
  Kokkos::deep_copy(X, hX);

  mlp_sgemm(d_mlp_Wg[0][0].data(), X.data(), C.data(), n_out, Mt, n_in,
      (int)d_mlp_Wg[0][0].extent(1), (int)X.extent(1), (int)C.extent(1), NNScalar(1.0));
  Kokkos::fence();

  auto hC = Kokkos::create_mirror_view(C);
  Kokkos::deep_copy(hC, C);
  double maxrel = 0.0;
  for (int b = 0; b < Mt; b++)
    for (int j = 0; j < n_out; j++) {
      double ref = 0.0;
      for (int k = 0; k < n_in; k++)
        ref += (double) hX(b, k) * mlp.W[0][(size_t) k * n_out + j];
      const double rel = std::fabs((double) hC(b, j) - ref) / std::max(1e-6, std::fabs(ref));
      maxrel = std::max(maxrel, rel);
    }
  // Threshold sized for the fp32 cublasSgemm this selftest exercises: correct
  // fp32-GEMM rounding is ~a few e-4 for realistic MLP widths, while a
  // transpose/leading-dim convention bug (what this gate guards against) yields
  // O(0.1-1) error. 7e-4 passes the former and still catches the latter.
  if (maxrel > 7e-4)
    error->all(FLERR, "GRACE-3L/KK: cuBLAS layer-0 selftest FAILED (max rel err {})", maxrel);
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-3L/KK] cuBLAS radial-MLP selftest OK (layer-0 max rel {:.2e})\n", maxrel);
}
#endif

// ---------------------- force launchers ----------------------

// Opt 15b: radial MLP forward value + r-derivative for layer `radial_mlp_spbf`,
// writing d_R1_nl (R), d_DR1_nl (dR/dr) and d_denv. On CUDA with a live handle
// and the expected 3-layer architecture, runs the batched true-fp32 cuBLAS path
// (6 SGEMMs: 3 value + 3 deriv, deriv layer-0 is a K=nradbase sub-GEMM since only
// the Chebyshev inputs are r-dependent); otherwise the hand kernel
// TagComputeMLPRadialDeriv. Same-stream serialization => race-free, no fences.
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_mlp_radial_deriv()
{
#ifdef KOKKOS_ENABLE_CUDA
  const auto &hmlp = grace_model->spbf[radial_mlp_spbf].mlp;
  if (cublas_handle && hmlp.n_layers == 3 && sizeof(NNScalar) == 4) {
    const int si = radial_mlp_spbf;
    const int bond_stride = (int) d_R1_nl.extent(1);
    const int M = chunk_size * bond_stride;
    mlp_M = M;
    const int n_in0 = hmlp.n_in[0];   // 138
    const int d1    = hmlp.n_out[0];  // 256
    const int d2    = hmlp.n_out[1];  // 128
    const int d3    = hmlp.n_out[2];  // n_rad_max*(lmax+1)
    mlp_assemble_idx = si;

    // 1. assemble X[M x n_in0], dX[M x nradbase], d_denv
    Kokkos::parallel_for("MLPAssembleDeriv",
        Kokkos::RangePolicy<DeviceType, TagMLPAssembleDeriv>(0, M), *this);

    // 2. layer 0: value h0 = X*W0 (K=n_in0); deriv dh0 = dX*W0[0:nradbase] (K=nradbase)
    mlp_sgemm(d_mlp_Wg[si][0].data(), d_mlp_X.data(),  d_mlp_h0.data(),
        d1, M, n_in0,    (int)d_mlp_Wg[si][0].extent(1), (int)d_mlp_X.extent(1),
        (int)d_mlp_h0.extent(1),  NNScalar(1.0));
    mlp_sgemm(d_mlp_Wg[si][0].data(), d_mlp_dX.data(), d_mlp_dh0.data(),
        d1, M, nradbase, (int)d_mlp_Wg[si][0].extent(1), (int)d_mlp_dX.extent(1),
        (int)d_mlp_dh0.extent(1), NNScalar(1.0));
    // 3. act-deriv layer 0
    mlp_act_layer = 0; mlp_act_n = d1; mlp_act_norm = (NNScalar) hmlp.norm[0];
    Kokkos::parallel_for("MLPActDeriv0",
        Kokkos::RangePolicy<DeviceType, TagMLPActDeriv>(0, M * d1), *this);

    // 4. layer 1: value h1 = h0*W1; deriv dh1 = dh0*W1
    mlp_sgemm(d_mlp_Wg[si][1].data(), d_mlp_h0.data(),  d_mlp_h1.data(),
        d2, M, d1, (int)d_mlp_Wg[si][1].extent(1), (int)d_mlp_h0.extent(1),
        (int)d_mlp_h1.extent(1),  NNScalar(1.0));
    mlp_sgemm(d_mlp_Wg[si][1].data(), d_mlp_dh0.data(), d_mlp_dh1.data(),
        d2, M, d1, (int)d_mlp_Wg[si][1].extent(1), (int)d_mlp_dh0.extent(1),
        (int)d_mlp_dh1.extent(1), NNScalar(1.0));
    // 5. act-deriv layer 1
    mlp_act_layer = 1; mlp_act_n = d2; mlp_act_norm = (NNScalar) hmlp.norm[1];
    Kokkos::parallel_for("MLPActDeriv1",
        Kokkos::RangePolicy<DeviceType, TagMLPActDeriv>(0, M * d2), *this);

    // 6. output layer (linear, no bias): value R=h1*W2*out_norm -> d_mlp_R -> d_R1_nl;
    //    deriv DR=dh1*W2*out_norm -> d_mlp_dR -> d_DR1_nl.
    const NNScalar out_norm = (NNScalar) hmlp.norm[2];
    mlp_sgemm(d_mlp_Wg[si][2].data(), d_mlp_h1.data(),  d_mlp_R.data(),
        d3, M, d2, (int)d_mlp_Wg[si][2].extent(1), (int)d_mlp_h1.extent(1),
        (int)d_mlp_R.extent(1),  out_norm);
    mlp_sgemm(d_mlp_Wg[si][2].data(), d_mlp_dh1.data(), d_mlp_dR.data(),
        d3, M, d2, (int)d_mlp_Wg[si][2].extent(1), (int)d_mlp_dh1.extent(1),
        (int)d_mlp_dR.extent(1), out_norm);
    mlp_scatter_which = 0;   // R -> d_R1_nl
    Kokkos::parallel_for("MLPScatterR",
        Kokkos::RangePolicy<DeviceType, TagMLPScatterR>(0, M * d3), *this);
    mlp_scatter_which = 1;   // DR -> d_DR1_nl
    Kokkos::parallel_for("MLPScatterDR",
        Kokkos::RangePolicy<DeviceType, TagMLPScatterR>(0, M * d3), *this);
    return;
  }
#endif
  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 32;
  check_team_size_for<TagComputeMLPRadialDeriv>(((chunk_size+team_size-1)/team_size)*maxneigh, team_size, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadialDeriv>(
      ((chunk_size+team_size-1)/team_size)*maxneigh, team_size, vector_length);
  Kokkos::parallel_for("ComputeMLPRadialDeriv", policy, *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_force_L1()
{
  const int nradmax = d_spbf[0].n_rad_max;
  Kokkos::parallel_for("ForceL1",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>, TagForceL1>({0, 0}, {chunk_size, nradmax}), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute_force_equiv(
    int spbf_idx, const t_nn_3d &A_adj, int ind_stage, int n_lm_ind)
{
  force_spbf_idx = spbf_idx;
  force_A_adj = A_adj;
  force_ind_stage = ind_stage;
  force_n_lm_ind = n_lm_ind;
  const int nradmax = d_spbf[spbf_idx].n_rad_max;
  // Opt 11: team-per-(ii,n) with the adj_prod row staged in team shared (old
  // per-thread TagForceEquiv kept for reference). team_size 64: the per-bond
  // body is register-heavy (bond_I[64] + w[25] + fp64 Plm/dPlm recurrence).
  const int nlm_y = (lmax + 1) * (lmax + 1);
  const int P = nlm_y * n_lm_ind;
  int team_size = 1, vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1) team_size = 64;
  const int league = chunk_size * nradmax;
  int scratch_size = scratch_size_helper<NNScalar>(P);
  check_team_size_for<TagForceEquivFused>(league, team_size, vector_length);
  auto policy = Kokkos::TeamPolicy<DeviceType, TagForceEquivFused>(league, team_size, vector_length);
  policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
  Kokkos::parallel_for("ForceEquivFused", policy, *this);
}

// ======================================================================
// MPI communication (host-side fallback) — buffers I / grad_I features.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
int PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  // Stage 2 (Task 2.10) packs the parity-doubled round-2 indicator eq2_norm from
  // h_eq2_norm_global; stage 1 packs eq1_norm from h_I_global. Same buffer layout.
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  auto &G = s2 ? h_eq2_norm_global : h_I_global;
  const int nmax = (int)G.extent(0);
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    if (j >= 0 && j < nmax) {
      for (int f = 0; f < nf; f++)
        for (int k = 0; k < no; k++)
          buf[m++] = G(j, f, k);
    } else {
      for (int f = 0; f < nf; f++)
        for (int k = 0; k < no; k++)
          buf[m++] = 0.0;
    }
  }
  return m;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  auto &G = s2 ? h_eq2_norm_global : h_I_global;
  const int nmax = (int)G.extent(0);
  for (int i = 0; i < n; i++) {
    const int j = first + i;
    if (j >= 0 && j < nmax) {
      for (int f = 0; f < nf; f++)
        for (int k = 0; k < no; k++)
          G(j, f, k) = buf[m++];
    } else {
      m += nf * no;  // skip
    }
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
int PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  // Stage-aware (Task 3.2), mirroring the forward pack: stage 1 reverse-comms the
  // round-1 grad_I from h_grad_I_global (Task 3.3); stage 2 reverse-comms the
  // parity-doubled round-2 eq2_norm message adjoint from h_eq2_norm_adj.
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  auto &G = s2 ? h_eq2_norm_adj : h_grad_I_global;
  for (int i = 0; i < n; i++) {
    const int j = first + i;
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        buf[m++] = G(j, f, k);
  }
  return m;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  auto &G = s2 ? h_eq2_norm_adj : h_grad_I_global;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        G(j, f, k) += buf[m++];
  }
}

// ======================================================================
// Kokkos-native communication: pack/unpack on device
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
int PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist_in,
                                                             DAT::tdual_double_1d &buf,
                                                             int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist_in.view<DeviceType>();
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for("PackForwardComm",
    Kokkos::RangePolicy<DeviceType, TagPackForwardComm>(0, n), *this);
  return n * (comm_stage == 2 ? I2_n_funcs * I2_n_out : I_n_funcs * I_n_out);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagPackForwardComm, const int& i) const
{
  const int j = d_sendlist(i);
  // comm_stage (captured by value in *this) selects round-1 (eq1_norm ->
  // d_I_global) vs round-2 (eq2_norm -> d_eq2_norm_global, parity-doubled lm).
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  int m = i * nf * no;
  if (s2) {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        v_buf(m++) = d_eq2_norm_global(j, f, k);
  } else {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        v_buf(m++) = d_I_global(j, f, k);
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::unpack_forward_comm_kokkos(int n, int first_in,
                                                                DAT::tdual_double_1d &buf)
{
  comm_first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for("UnpackForwardComm",
    Kokkos::RangePolicy<DeviceType, TagUnpackForwardComm>(0, n), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagUnpackForwardComm, const int& i) const
{
  const int j = comm_first + i;
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  int m = i * nf * no;
  if (s2) {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        d_eq2_norm_global(j, f, k) = v_buf(m++);
  } else {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        d_I_global(j, f, k) = v_buf(m++);
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
int PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::pack_reverse_comm_kokkos(int n, int first_in,
                                                              DAT::tdual_double_1d &buf)
{
  comm_first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for("PackReverseComm",
    Kokkos::RangePolicy<DeviceType, TagPackReverseComm>(0, n), *this);
  return n * (comm_stage == 2 ? I2_n_funcs * I2_n_out : I_n_funcs * I_n_out);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagPackReverseComm, const int& i) const
{
  const int j = comm_first + i;
  // Stage 2 (Task 3.2): pack the round-2 eq2_norm message adjoint (d_eq2_norm_adj,
  // parity-doubled lm); stage 1: pack round-1 grad_I (d_grad_I_global, Task 3.3).
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  int m = i * nf * no;
  if (s2) {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        v_buf(m++) = d_eq2_norm_adj(j, f, k);
  } else {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        v_buf(m++) = d_grad_I_global(j, f, k);
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist_in,
                                                                 DAT::tdual_double_1d &buf)
{
  d_sendlist = k_sendlist_in.view<DeviceType>();
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for("UnpackReverseComm",
    Kokkos::RangePolicy<DeviceType, TagUnpackReverseComm>(0, n), *this);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagUnpackReverseComm, const int& i) const
{
  const int j = d_sendlist(i);
  const bool s2 = (comm_stage == 2);
  const int nf = s2 ? I2_n_funcs : I_n_funcs;
  const int no = s2 ? I2_n_out   : I_n_out;
  int m = i * nf * no;
  if (s2) {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        Kokkos::atomic_add(&d_eq2_norm_adj(j, f, k), v_buf(m++));
  } else {
    for (int f = 0; f < nf; f++)
      for (int k = 0; k < no; k++)
        Kokkos::atomic_add(&d_grad_I_global(j, f, k), v_buf(m++));
  }
}

// ======================================================================
// Task 3.5b: forward-chunk recompute helpers. l1_forward_chunk() runs the full
// Layer-1 forward body (geometry -> A1 chain -> eq1_norm/rho1_norm + the natom
// ScatterRho1Norm) for the CURRENT chunk (member chunk_offset/chunk_size), but
// STOPS BEFORE the global eq1_norm scatter (P1 does that; the backward recompute
// in P5 must NOT touch d_I_global). l2_forward_chunk() is the Layer-2 analogue
// (through rho2_norm + ScatterRho2Norm, excluding the eq2_norm global scatter).
// Because the inter-layer indicators d_I_global/d_eq2_norm_global are already
// ghost-complete (forward_comm ran in P1.5/P2.5), re-running these bodies in
// P4/P5 reproduces bit-identical per-chunk forward intermediates.
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::l1_forward_chunk()
{
  int team_size = 1;
  int vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1)
    team_size = 32;

  // ---- Stage 1: ComputeNeigh ----
  {
    check_team_size_for<TagComputeNeigh>(chunk_size, team_size, vector_length);
    int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>(chunk_size, team_size, vector_length);
    policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for("ComputeNeigh", policy, *this);
  }

  // ---- Stage 1b: ComputeYBond — precompute real spherical harmonics per bond ----
  compute_Y_bond_chunk();

  // ---- Stage 2: ComputeRadialBasis — Chebyshev basis (raw) + envelope ----
  {
    int ts = team_size;
    check_team_size_for<TagComputeRadialBasis>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>(
        ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
    Kokkos::parallel_for("ComputeRadialBasis", policy, *this);
  }

  // ---- Stage 3: ComputeMLPRadial — Layer-1 radial MLP (Task 2.4) ----
  radial_mlp_spbf = 0;   // A1
  compute_mlp_radial();  // Opt 15: cuBLAS SGEMM path (or hand-kernel fallback)

  // ---- Stage 4: ComputeA1 — Layer-1 scalar-indicator SPBF basis (Task 2.4) ----
  {
    int ts = team_size;
    check_team_size_for<TagComputeA1>(((chunk_size+ts-1)/ts)*L1_nradmax, ts, vector_length);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeA1>(
        ((chunk_size+ts-1)/ts)*L1_nradmax, ts, vector_length);
    Kokkos::parallel_for("ComputeA1", policy, *this);
  }

  // ---- Stage 5: cp_l product A1_2 = cp_l(A1, A1) (Task 2.5) ----
  compute_cp_l(d_prod[0], d_A1, d_A1, d_A1_2);

  // ---- Stage 6: FunctionReduceN A1_2_red = reduce_n({"A1_2": A1_2}) (Task 2.6) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A1_2;   // matches reduce["A1_2_red"].instr[0].name == "A1_2" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A1_2_red], inputs, d_A1_2_red);
  }

  // ---- Stage 7: FCRight2Left A1_2a = fc_right2left(left=A1_2_red, right=A1) (Task 2.6) ----
  compute_fc(d_fc[idx_fc_A1_2a], d_A1_2_red, d_A1, d_A1_2a);

  // ---- Stage 8: cp_l product A1_3 = cp_l(A1_2a, A1) (Task 2.7) ----
  compute_cp_l(d_prod[idx_prod_A1_3], d_A1_2a, d_A1, d_A1_3);

  // ---- Stage 9: FunctionReduceN A1_3_red = reduce_n({"A1_3": A1_3}) (Task 2.7) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A1_3;   // matches reduce["A1_3_red"].instr[0].name == "A1_3" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A1_3_red], inputs, d_A1_3_red);
  }

  // ---- Stage 10: FCRight2Left A1_2b = fc_right2left(left=A1_2_red, right=A1) (Task 2.7) ----
  compute_fc(d_fc[idx_fc_A1_2b], d_A1_2_red, d_A1, d_A1_2b);

  // ---- Stage 11: cp_l product A1_4 = cp_l(A1_2b, A1_2b) (Task 2.7) ----
  compute_cp_l(d_prod[idx_prod_A1_4], d_A1_2b, d_A1_2b, d_A1_4);

  // ---- Stage 12: FunctionReduceN A1_4_red = reduce_n({"A1_4": A1_4}) (Task 2.7) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A1_4;   // matches reduce["A1_4_red"].instr[0].name == "A1_4" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A1_4_red], inputs, d_A1_4_red);
  }

  // ---- Stage 13: FunctionReduceN eq1 = reduce_n({A1,A1_2_red,A1_3_red,A1_4_red}) ----
  {
    t_nn_3d l1_views[4] = {d_A1, d_A1_2_red, d_A1_3_red, d_A1_4_red};
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    const int n_instr = d_reduce[idx_reduce_eq1].n_instr;
    for (int j = 0; j < n_instr; j++) inputs[j] = l1_views[eq1_input_role[j]];
    compute_reduceN(d_reduce[idx_reduce_eq1], inputs, d_eq1);
  }

  // ---- Stage 14: EquivariantRMSNorm eq1_norm = equivariant_rms_norm(eq1) ----
  compute_equiv_rms_norm(d_eqnorm[idx_eqnorm_eq1_norm], d_eq1, d_eq1_norm);
  // NOTE: Stage 14b (ScatterEq1Norm -> d_I_global) is done by the P1 caller ONLY.

  // ---- Stage 15: FunctionReduceN rho1 = reduce_n({A1,A1_2_red,A1_3_red,A1_4_red}) ----
  {
    t_nn_3d l1_views[4] = {d_A1, d_A1_2_red, d_A1_3_red, d_A1_4_red};
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    const int n_instr = d_reduce[idx_reduce_rho1].n_instr;
    for (int j = 0; j < n_instr; j++) inputs[j] = l1_views[rho1_input_role[j]];
    compute_reduceN(d_reduce[idx_reduce_rho1], inputs, d_rho1);
  }

  // ---- Stage 16: InvariantLayerRMSNorm rho1_norm = invariant_rms_norm(rho1) ----
  compute_invariant_rms_norm(d_rmsnorm[idx_rmsnorm_rho1_norm], d_rho1, d_rho1_norm);

  // ---- Stage 16b (Task 3.5a): scatter rho1_norm into the natom buffer
  // d_rho1_norm_full (dense list-order index ii+chunk_offset; local-only, no comm). ----
  {
    auto rho1n = d_rho1_norm;
    auto full = d_rho1_norm_full;
    const int co = chunk_offset;
    const int cs = chunk_size;
    const int nc = (int) rho1n.extent(1), nl = (int) rho1n.extent(2);
    Kokkos::parallel_for("ScatterRho1Norm",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nc, nl}),
      KOKKOS_LAMBDA(const int ii, const int c, const int lm) {
        full(ii + co, c, lm) = rho1n(ii, c, lm);
      });
  }
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::l2_forward_chunk()
{
  int team_size = 1;
  int vector_length = 1;
  if (Kokkos::DefaultExecutionSpace().concurrency() > 1)
    team_size = 32;

  // ---- L2 Stage 0: recompute neighbor data for this chunk ----
  {
    check_team_size_for<TagComputeNeigh>(chunk_size, team_size, vector_length);
    int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>(chunk_size, team_size, vector_length);
    policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for("L2_ComputeNeigh", policy, *this);
  }

  // ---- L2 Stage 0b: ComputeYBond ----
  compute_Y_bond_chunk();

  // ---- L2 Stage 0c: ComputeRadialBasis ----
  {
    int ts = team_size;
    check_team_size_for<TagComputeRadialBasis>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
    auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>(
        ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
    Kokkos::parallel_for("L2_ComputeRadialBasis", policy, *this);
  }

  // ---- L2 Stage 1: ComputeMLPRadial — A2 radial MLP (spbf idx_spbf_A2) ----
  radial_mlp_spbf = idx_spbf_A2;
  compute_mlp_radial();  // Opt 15: cuBLAS SGEMM path (or hand-kernel fallback)

  // ---- L2 Stage 2: ComputeA2 — equiv-indicator SPBF basis (Task 2.8) ----
  // Opt 6: fused team+shared kernel (indicator = d_I_global, out = d_A2).
  compute_spbf_equiv(idx_spbf_A2, d_I_global, d_A2);

  // ---- L2 Stage 3: FunctionReduceN A2_red = reduce_n({"A2": A2}) (Task 2.9) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A2;   // matches reduce["A2_red"].instr[0].name == "A2" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A2_red], inputs, d_A2_red);
  }

  // ---- L2 Stage 4: cp_l product A2_2 = cp_l(A2_red, A2_red) (Task 2.9) ----
  compute_cp_l(d_prod[idx_prod_A2_2], d_A2_red, d_A2_red, d_A2_2);

  // ---- L2 Stage 5: FunctionReduceN A2_2_red = reduce_n({"A2_2": A2_2}) (Task 2.9) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A2_2;   // matches reduce["A2_2_red"].instr[0].name == "A2_2" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A2_2_red], inputs, d_A2_2_red);
  }

  // ---- L2 Stage 6: FCRight2Left A2_2a = fc_right2left(left=A2_2_red, right=A2_red) (Task 2.9) ----
  compute_fc(d_fc[idx_fc_A2_2a], d_A2_2_red, d_A2_red, d_A2_2a);

  // ---- L2 Stage 7: cp_l product A2_3 = cp_l(A2_2a, A2_red) (Task 2.9) ----
  compute_cp_l(d_prod[idx_prod_A2_3], d_A2_2a, d_A2_red, d_A2_3);

  // ---- L2 Stage 8: FunctionReduceN A2_3_red = reduce_n({"A2_3": A2_3}) (Task 2.9) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A2_3;   // matches reduce["A2_3_red"].instr[0].name == "A2_3" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A2_3_red], inputs, d_A2_3_red);
  }

  // ---- L2 Stage 9: FCRight2Left A2_2b = fc_right2left(left=A2_2_red, right=A2_red) (Task 2.9) ----
  compute_fc(d_fc[idx_fc_A2_2b], d_A2_2_red, d_A2_red, d_A2_2b);

  // ---- L2 Stage 10: cp_l product A2_4 = cp_l(A2_2b, A2_2b) (Task 2.9) ----
  compute_cp_l(d_prod[idx_prod_A2_4], d_A2_2b, d_A2_2b, d_A2_4);

  // ---- L2 Stage 11: FunctionReduceN A2_4_red = reduce_n({"A2_4": A2_4}) (Task 2.9) ----
  {
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    inputs[0] = d_A2_4;   // matches reduce["A2_4_red"].instr[0].name == "A2_4" (checked in init_style)
    compute_reduceN(d_reduce[idx_reduce_A2_4_red], inputs, d_A2_4_red);
  }

  // ---- L2 Stage 12: FunctionReduceN eq2 = reduce_n({A2_red,A2_2_red,A2_3_red,A2_4_red}) ----
  {
    t_nn_3d l2_views[4] = {d_A2_red, d_A2_2_red, d_A2_3_red, d_A2_4_red};
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    const int n_instr = d_reduce[idx_reduce_eq2].n_instr;
    for (int j = 0; j < n_instr; j++) inputs[j] = l2_views[eq2_input_role[j]];
    compute_reduceN(d_reduce[idx_reduce_eq2], inputs, d_eq2);
  }

  // ---- L2 Stage 13: EquivariantRMSNorm eq2_norm = equivariant_rms_norm(eq2) ----
  compute_equiv_rms_norm(d_eqnorm[idx_eqnorm_eq2_norm], d_eq2, d_eq2_norm);
  // NOTE: Stage 13b (ScatterEq2Norm -> d_eq2_norm_global) is done by the P2 caller ONLY.

  // ---- L2 Stage 14: FunctionReduceN rho2 = reduce_n({A2_red,A2_2_red,A2_3_red,A2_4_red}) ----
  {
    t_nn_3d l2_views[4] = {d_A2_red, d_A2_2_red, d_A2_3_red, d_A2_4_red};
    t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
    const int n_instr = d_reduce[idx_reduce_rho2].n_instr;
    for (int j = 0; j < n_instr; j++) inputs[j] = l2_views[rho2_input_role[j]];
    compute_reduceN(d_reduce[idx_reduce_rho2], inputs, d_rho2);
  }

  // ---- L2 Stage 15: InvariantLayerRMSNorm rho2_norm = invariant_rms_norm(rho2) ----
  compute_invariant_rms_norm(d_rmsnorm[idx_rmsnorm_rho2_norm], d_rho2, d_rho2_norm);

  // ---- L2 Stage 15b (Task 3.5a): scatter rho2_norm into the natom buffer
  // d_rho2_norm_full (dense list-order index ii+chunk_offset; local-only, no comm). ----
  {
    auto rho2n = d_rho2_norm;
    auto full = d_rho2_norm_full;
    const int co = chunk_offset;
    const int cs = chunk_size;
    const int nc = (int) rho2n.extent(1), nl = (int) rho2n.extent(2);
    Kokkos::parallel_for("ScatterRho2Norm",
      Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nc, nl}),
      KOKKOS_LAMBDA(const int ii, const int c, const int lm) {
        full(ii + co, c, lm) = rho2n(ii, c, lm);
      });
  }
}

// ======================================================================
// Task 3.5b: scatter_forces_chunk — scatter the CURRENT chunk's per-bond force
// buffer d_f_ij into the global force array f[] (atomic; F_i += f_ij, F_j -=
// f_ij) and accumulate the per-bond virial into virial[] (+=) / d_vatom /
// d_cvatom. Extracted VERBATIM from the old single-shot PHASE 6b. Called once
// per layer per chunk (after that layer's force kernel), so each layer's
// partial force accumulates into f[]/virial[] independently -> identical total.
// ======================================================================
template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::scatter_forces_chunk()
{
  auto fij = d_f_ij; auto nc = d_ncount; auto near = d_nearest;
  auto fout = f; auto il = d_ilist;
  auto rh = d_rhats; auto rn = d_rnorms;
  int co = chunk_offset, cs = chunk_size;
  bool do_v = (vflag_global != 0);
  bool do_va = (vflag_atom != 0);
  bool do_cva = (cvflag_atom != 0);
  auto va = d_vatom; auto cva = d_cvatom;

  EV_FLOAT ev_force;
  Kokkos::parallel_reduce("ComputeForce3L", Kokkos::RangePolicy<DeviceType>(0, cs),
    KOKKOS_LAMBDA(const int ii, EV_FLOAT& ev) {
      const int i = il[ii + co];
      const int nct = nc(ii);
      for (int jj = 0; jj < nct; jj++) {
        const int j = near(ii, jj);
        const GeomScalar fx = fij(ii, jj, 0);
        const GeomScalar fy = fij(ii, jj, 1);
        const GeomScalar fz = fij(ii, jj, 2);
        Kokkos::atomic_add(&fout(i, 0), fx);
        Kokkos::atomic_add(&fout(i, 1), fy);
        Kokkos::atomic_add(&fout(i, 2), fz);
        Kokkos::atomic_add(&fout(j, 0), -fx);
        Kokkos::atomic_add(&fout(j, 1), -fy);
        Kokkos::atomic_add(&fout(j, 2), -fz);
        if (do_v || do_va || do_cva) {
          const GeomScalar delx = -rh(ii, jj, 0) * rn(ii, jj);
          const GeomScalar dely = -rh(ii, jj, 1) * rn(ii, jj);
          const GeomScalar delz = -rh(ii, jj, 2) * rn(ii, jj);
          const GeomScalar v0 = delx*fx, v1 = dely*fy, v2 = delz*fz;
          const GeomScalar v3 = delx*fy, v4 = delx*fz, v5 = dely*fz;
          if (do_v) { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
          if (do_va) {
            Kokkos::atomic_add(&va(i,0),(KK_FLOAT)(GeomScalar(0.5)*v0)); Kokkos::atomic_add(&va(i,1),(KK_FLOAT)(GeomScalar(0.5)*v1));
            Kokkos::atomic_add(&va(i,2),(KK_FLOAT)(GeomScalar(0.5)*v2)); Kokkos::atomic_add(&va(i,3),(KK_FLOAT)(GeomScalar(0.5)*v3));
            Kokkos::atomic_add(&va(i,4),(KK_FLOAT)(GeomScalar(0.5)*v4)); Kokkos::atomic_add(&va(i,5),(KK_FLOAT)(GeomScalar(0.5)*v5));
            Kokkos::atomic_add(&va(j,0),(KK_FLOAT)(GeomScalar(0.5)*v0)); Kokkos::atomic_add(&va(j,1),(KK_FLOAT)(GeomScalar(0.5)*v1));
            Kokkos::atomic_add(&va(j,2),(KK_FLOAT)(GeomScalar(0.5)*v2)); Kokkos::atomic_add(&va(j,3),(KK_FLOAT)(GeomScalar(0.5)*v3));
            Kokkos::atomic_add(&va(j,4),(KK_FLOAT)(GeomScalar(0.5)*v4)); Kokkos::atomic_add(&va(j,5),(KK_FLOAT)(GeomScalar(0.5)*v5));
          }
          if (do_cva) {
            const GeomScalar v6=dely*fx, v7=delz*fx, v8=delz*fy;
            Kokkos::atomic_add(&cva(i,0),(KK_FLOAT)(GeomScalar(0.5)*v0)); Kokkos::atomic_add(&cva(i,1),(KK_FLOAT)(GeomScalar(0.5)*v1));
            Kokkos::atomic_add(&cva(i,2),(KK_FLOAT)(GeomScalar(0.5)*v2)); Kokkos::atomic_add(&cva(i,3),(KK_FLOAT)(GeomScalar(0.5)*v3));
            Kokkos::atomic_add(&cva(i,4),(KK_FLOAT)(GeomScalar(0.5)*v4)); Kokkos::atomic_add(&cva(i,5),(KK_FLOAT)(GeomScalar(0.5)*v5));
            Kokkos::atomic_add(&cva(i,6),(KK_FLOAT)(GeomScalar(0.5)*v6)); Kokkos::atomic_add(&cva(i,7),(KK_FLOAT)(GeomScalar(0.5)*v7));
            Kokkos::atomic_add(&cva(i,8),(KK_FLOAT)(GeomScalar(0.5)*v8));
            Kokkos::atomic_add(&cva(j,0),(KK_FLOAT)(GeomScalar(0.5)*v0)); Kokkos::atomic_add(&cva(j,1),(KK_FLOAT)(GeomScalar(0.5)*v1));
            Kokkos::atomic_add(&cva(j,2),(KK_FLOAT)(GeomScalar(0.5)*v2)); Kokkos::atomic_add(&cva(j,3),(KK_FLOAT)(GeomScalar(0.5)*v3));
            Kokkos::atomic_add(&cva(j,4),(KK_FLOAT)(GeomScalar(0.5)*v4)); Kokkos::atomic_add(&cva(j,5),(KK_FLOAT)(GeomScalar(0.5)*v5));
            Kokkos::atomic_add(&cva(j,6),(KK_FLOAT)(GeomScalar(0.5)*v6)); Kokkos::atomic_add(&cva(j,7),(KK_FLOAT)(GeomScalar(0.5)*v7));
            Kokkos::atomic_add(&cva(j,8),(KK_FLOAT)(GeomScalar(0.5)*v8));
          }
        }
      }
    }, ev_force);
  if (vflag_global) {
    virial[0]+=ev_force.v[0]; virial[1]+=ev_force.v[1]; virial[2]+=ev_force.v[2];
    virial[3]+=ev_force.v[3]; virial[4]+=ev_force.v[4]; virial[5]+=ev_force.v[5];
  }
}

// ======================================================================
// compute: full forward BASIS chain (L1/L2/L3) + 3-density readout (Task
// 2.11) -> per-atom/total energy. NO forces yet: d_e_atom/eng_vdwl/d_eatom
// are the only physics contributions; the force kernels land in Phase 3.
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag, vflag, 0);

  // eflag_only is populated by ev_init above; the energy-only fast path skips the
  // entire backward (L3 Stage 15 + PHASES 3.5/4/4.5/5 + per-chunk force scatter)
  // when the caller requests energy without forces (e.g. MC fixes, fix_numdiff,
  // compute_fep). UQ is forward-only and still runs. `debug_no_energy_only_calc`
  // (pair_style keyword) forces the full backward for validation/benchmarking.
  const bool do_energy_only = eflag_only && !debug_no_energy_only_calc;

  // (Re)allocate and zero the per-atom accumulators this style owns.
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
    Kokkos::deep_copy(d_eatom, 0.0);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
    Kokkos::deep_copy(d_vatom, 0.0);
  }
  if (cvflag_atom) {
    memoryKK->destroy_kokkos(k_cvatom, cvatom);
    memoryKK->create_kokkos(k_cvatom, cvatom, maxcvatom, "pair:cvatom");
    d_cvatom = k_cvatom.view<DeviceType>();
    Kokkos::deep_copy(d_cvatom, 0.0);
  }

  copymode = 1;

  atomKK->sync(execution_space, X_MASK | F_MASK | TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  maxneigh = 0;
  Kokkos::parallel_reduce("grace_3l::find_maxneigh", inum,
      FindMaxNumNeighs3L<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  grow(chunk_size, maxneigh);

  // Global (nall-sized) arrays for the inter-layer message comm (Task 2.8).
  const int nall = atom->nlocal + atom->nghost;
  grow_global(nall);
  // Zero the whole indicator array: owned rows are overwritten by the L1 loop,
  // ghost rows by comm->forward_comm; zeroing guards any atom left untouched.
  Kokkos::deep_copy(d_I_global, NNScalar(0.0));
  // Round-2 (L3) indicator buffer: owned rows written by the L2 loop's eq2_norm
  // scatter, ghost rows by the round-2 forward_comm (Task 2.10). Zero first.
  Kokkos::deep_copy(d_eq2_norm_global, NNScalar(0.0));

  // ---- UQ / extrapolation grade: allocate per-atom outputs up front. The
  // basis-RP feature is assembled in three parts across the forward: rho1 in
  // Phase 1, rho2 in Phase 2 (both accumulate the RAW proj into d_uq_z), and
  // rho3 + the GMM in Phase 3. All UQ kernels are any_uq-gated, so the energy
  // path is byte-for-byte unchanged when UQ is not requested. ----
  const bool any_uq = has_uq && (flag_compute_gamma || flag_compute_atomic_sigma || flag_compute_gmm_cluster);
  if (any_uq) {
    const int nlocal = atom->nlocal;
    if (nlocal > nmax_uq) {
      memory->destroy(gamma);          memory->create(gamma, nlocal, "grace:gamma");
      memory->destroy(atomic_sigma);   memory->create(atomic_sigma, nlocal, "grace:atomic_sigma");
      memory->destroy(gmm_cluster);    memory->create(gmm_cluster, nlocal, "grace:gmm_cluster");
      nmax_uq = nlocal;
    }
    // Also reallocate if d_uq_z's projection width no longer matches uq_rp_dim: a
    // mid-session pair_coeff that loads a model with a different rp_dim (while nlocal
    // has not grown) would otherwise leave d_uq_z sized to the OLD width, and the
    // rho1 kernel would write d_uq_z(i,d) out of bounds.
    if ((int)d_gamma.extent(0) < nlocal || (int)d_uq_z.extent(1) != uq_rp_dim) {
      MemKK::realloc_kokkos(d_gamma, "g3l:d_gamma", nlocal);
      MemKK::realloc_kokkos(d_sigma, "g3l:d_sigma", nlocal);
      MemKK::realloc_kokkos(d_gmm_cluster, "g3l:d_gmm_cluster", nlocal);
      MemKK::realloc_kokkos(d_uq_z, "g3l:d_uq_z", nlocal, uq_rp_dim);  // per-atom RAW proj (rho1+rho2 carry)
      MemKK::realloc_kokkos(d_uq_n2_rho1, "g3l:d_uq_n2_rho1", nlocal); // ||B^(rho1)||^2 carry
      MemKK::realloc_kokkos(d_uq_n2_rho2, "g3l:d_uq_n2_rho2", nlocal); // ||B^(rho2)||^2 carry
    }
  }

  // ============ PHASE 1: Layer-1 forward (chunked, ALL owned atoms) ============
  // Runs the full L1 pipeline (geometry -> A1 chain -> eq1_norm/rho1_norm) for
  // every owned atom, scattering each atom's eq1_norm into the global indicator
  // array d_I_global at its LAMMPS local index. After all chunks, forward_comm
  // fills the ghost rows; then Phase 2 computes A2 reading indicator[ind_j].
  chunk_offset = 0;
  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    // ---- L1 forward chunk body (Stages 1-16, incl. ScatterRho1Norm), EXCEPT
    // the Stage-14b eq1_norm global scatter which P1 does below (and P5's
    // backward recompute must NOT do). See l1_forward_chunk(). ----
    l1_forward_chunk();

    // ---- Stage 14b: scatter eq1_norm into the global indicator array (Task 2.8).
    // d_eq1_norm is chunk-indexed (ii); write it to d_I_global at the atom's
    // LAMMPS local index i=d_ilist(ii+chunk_offset) so forward_comm can then
    // fill the ghost rows. Layout matches: (n_out, M) == (I_n_funcs, I_n_out). ----
    {
      auto eq1n = d_eq1_norm;
      auto Ig = d_I_global;
      auto il = d_ilist;
      const int co = chunk_offset;
      const int nf = I_n_funcs, no = I_n_out, cs = chunk_size;
      Kokkos::parallel_for("ScatterEq1Norm",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nf, no}),
        KOKKOS_LAMBDA(const int ii, const int fdim, const int k) {
          Ig(il(ii + co), fdim, k) = eq1n(ii, fdim, k);
        });
    }

    // ---- UQ basis-RP (Phase 1): project the rho1 (L1) scalar basis into d_uq_z
    // for this chunk (d_A1/d_A1_2_red/d_A1_3_red/d_A1_4_red are current here). ----
    if (any_uq)
      Kokkos::parallel_for("ComputeUQ_rho1",
          Kokkos::RangePolicy<DeviceType, TagComputeUQ_rho1>(0, chunk_size), *this);

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  // ============ PHASE 1.5: forward-comm eq1_norm to ghost atoms (Task 2.8) ====
  // Every owned atom's eq1_norm now lives in d_I_global; forward_comm copies
  // owned rows into their ghost images so the A2 kernel can read the neighbor
  // indicator at ghost ind_j. When legacy pair comm is active (GPU-aware MPI
  // off / multi-proc host path), LAMMPS pack/unpack reads h_I_global, so sync
  // device->host before and host->device after; otherwise the kokkos-native
  // pack/unpack operates on d_I_global directly. (Mirrors GRACE-2L.)
  comm_stage = 1;   // round 1: pack/unpack eq1_norm from d_I_global (Task 2.10)
  if (lmp->kokkos->forward_pair_comm_legacy) {
    Kokkos::deep_copy(h_I_global, d_I_global);
    comm->forward_comm(this);
    Kokkos::deep_copy(d_I_global, h_I_global);
  } else {
    comm->forward_comm(this);
  }

  // ============ PHASE 2: Layer-2 forward (chunked) — A2 ============
  // The per-chunk geometry arrays (d_ncount/d_nearest/d_rhats/d_radial_basis/
  // d_Y_bond/d_R1_nl) hold the LAST L1 chunk's data, so each L2 chunk must
  // recompute geometry + the A2 radial MLP before assembling A2.
  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    // ---- L2 forward chunk body (Stages 0-15, incl. ScatterRho2Norm), EXCEPT
    // the Stage-13b eq2_norm global scatter which P2 does below (and P4's
    // backward recompute must NOT do). See l2_forward_chunk(). ----
    l2_forward_chunk();

    // ---- L2 Stage 13b: scatter eq2_norm into the round-2 global indicator
    // array (Task 2.10). d_eq2_norm is chunk-indexed (ii); write it to
    // d_eq2_norm_global at the atom's LAMMPS local index i=d_ilist(ii+chunk_offset)
    // so the round-2 forward_comm (after all L2 chunks) can fill the ghost rows.
    // Layout matches: (n_out, M) == (I2_n_funcs, I2_n_out). Mirrors Stage 14b. ----
    {
      auto eq2n = d_eq2_norm;
      auto Ig = d_eq2_norm_global;
      auto il = d_ilist;
      const int co = chunk_offset;
      const int nf = I2_n_funcs, no = I2_n_out, cs = chunk_size;
      Kokkos::parallel_for("ScatterEq2Norm",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nf, no}),
        KOKKOS_LAMBDA(const int ii, const int fdim, const int k) {
          Ig(il(ii + co), fdim, k) = eq2n(ii, fdim, k);
        });
    }

    // ---- UQ basis-RP (Phase 2): accumulate the rho2 (L2) scalar basis into d_uq_z
    // (d_A2_red/d_A2_2_red/d_A2_3_red/d_A2_4_red are current here). ----
    if (any_uq)
      Kokkos::parallel_for("ComputeUQ_rho2",
          Kokkos::RangePolicy<DeviceType, TagComputeUQ_rho2>(0, chunk_size), *this);

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  // ============ PHASE 2.5: forward-comm eq2_norm to ghost atoms (Task 2.10) ====
  // Every owned atom's eq2_norm now lives in d_eq2_norm_global; forward_comm
  // copies owned rows into their ghost images so the A3 kernel can read the
  // neighbor indicator at ghost ind_j. comm_stage=2 selects the round-2 buffer +
  // parity-doubled dims (I2_n_funcs*I2_n_out=3200/atom) in the pack/unpack (host
  // AND kokkos-native) paths. Mirrors the round-1 forward_comm block above.
  comm_stage = 2;
  if (lmp->kokkos->forward_pair_comm_legacy) {
    Kokkos::deep_copy(h_eq2_norm_global, d_eq2_norm_global);
    comm->forward_comm(this);
    Kokkos::deep_copy(d_eq2_norm_global, h_eq2_norm_global);
  } else {
    comm->forward_comm(this);
  }

  // ============ PHASE 3: Layer-3 forward (chunked) — A3 chain (TERMINAL) ======
  // The per-chunk geometry arrays hold the LAST L2 chunk's data, so each L3
  // chunk must recompute geometry + the A3 radial MLP before assembling A3.
  // A3's indicator is eq2_norm (in d_eq2_norm_global). L3 is TERMINAL: the
  // A3_red..rho3 chain mirrors L2's A2_red..rho2, but ends at rho3/rho3_norm —
  // there is NO eq3.
  // Task 3.5b: zero the L3 message adjoint d_eq2_norm_adj ONCE, BEFORE the P3
  // loop. It is an [nall] buffer that compute_A3_spbf_bwd atomic_add-scatters
  // across ALL chunks; zeroing it per-chunk (its old spot, inside the loop)
  // would wipe earlier chunks' contributions under multi-chunk.
  if (!do_energy_only) Kokkos::deep_copy(d_eq2_norm_adj, NNScalar(0.0));
  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    int team_size = 1;
    int vector_length = 1;
    if (Kokkos::DefaultExecutionSpace().concurrency() > 1)
      team_size = 32;

    // ---- L3 Stage 0: recompute neighbor data for this chunk ----
    {
      check_team_size_for<TagComputeNeigh>(chunk_size, team_size, vector_length);
      int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeNeigh>(chunk_size, team_size, vector_length);
      policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("L3_ComputeNeigh", policy, *this);
    }

    // ---- L3 Stage 0b: ComputeYBond ----
    compute_Y_bond_chunk();

    // ---- L3 Stage 0c: ComputeRadialBasis ----
    {
      int ts = team_size;
      check_team_size_for<TagComputeRadialBasis>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>(
          ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      Kokkos::parallel_for("L3_ComputeRadialBasis", policy, *this);
    }

    // ---- L3 Stage 1: A3 radial MLP (spbf idx_spbf_A3). Opt 13: run the DERIV
    // variant here instead of the plain forward — it reproduces the forward's
    // d_R1_nl bit-identically AND fills d_DR1_nl/d_denv, which Stage 15's force
    // path needs; nothing between Stage 1 and Stage 15 writes those buffers
    // (SPBFFused/reduce/cp_l/fc only READ d_R1_nl), so the Stage-15 radial
    // recompute is dropped. Saves one full radial-MLP forward per eval. ----
    radial_mlp_spbf = idx_spbf_A3;
    compute_mlp_radial_deriv();

    // ---- L3 Stage 2: ComputeA3 — equiv-indicator SPBF basis (indicator=eq2_norm) ----
    // Opt 6: fused team+shared kernel (indicator = d_eq2_norm_global, out = d_A3).
    compute_spbf_equiv(idx_spbf_A3, d_eq2_norm_global, d_A3);

    // ---- L3 Stage 3: FunctionReduceN A3_red = reduce_n({"A3": A3}). L3's
    // downstream products/FCs are ALL built on this reduced A3_red (like L2). ----
    {
      t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
      inputs[0] = d_A3;   // matches reduce["A3_red"].instr[0].name == "A3" (checked in init_style)
      compute_reduceN(d_reduce[idx_reduce_A3_red], inputs, d_A3_red);
    }

    // ---- L3 Stage 4: cp_l product A3_2 = cp_l(A3_red, A3_red) ----
    compute_cp_l(d_prod[idx_prod_A3_2], d_A3_red, d_A3_red, d_A3_2);

    // ---- L3 Stage 5: FunctionReduceN A3_2_red = reduce_n({"A3_2": A3_2}) ----
    {
      t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
      inputs[0] = d_A3_2;   // matches reduce["A3_2_red"].instr[0].name == "A3_2"
      compute_reduceN(d_reduce[idx_reduce_A3_2_red], inputs, d_A3_2_red);
    }

    // ---- L3 Stage 6: FCRight2Left A3_2a = fc_right2left(left=A3_2_red, right=A3_red) ----
    compute_fc(d_fc[idx_fc_A3_2a], d_A3_2_red, d_A3_red, d_A3_2a);

    // ---- L3 Stage 7: cp_l product A3_3 = cp_l(A3_2a, A3_red) (Lmax=0) ----
    compute_cp_l(d_prod[idx_prod_A3_3], d_A3_2a, d_A3_red, d_A3_3);

    // ---- L3 Stage 8: FunctionReduceN A3_3_red = reduce_n({"A3_3": A3_3}) ----
    {
      t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
      inputs[0] = d_A3_3;   // matches reduce["A3_3_red"].instr[0].name == "A3_3"
      compute_reduceN(d_reduce[idx_reduce_A3_3_red], inputs, d_A3_3_red);
    }

    // ---- L3 Stage 9: FCRight2Left A3_2b = fc_right2left(left=A3_2_red, right=A3_red) ----
    compute_fc(d_fc[idx_fc_A3_2b], d_A3_2_red, d_A3_red, d_A3_2b);

    // ---- L3 Stage 10: cp_l product A3_4 = cp_l(A3_2b, A3_2b) (Lmax=0) ----
    compute_cp_l(d_prod[idx_prod_A3_4], d_A3_2b, d_A3_2b, d_A3_4);

    // ---- L3 Stage 11: FunctionReduceN A3_4_red = reduce_n({"A3_4": A3_4}) ----
    {
      t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
      inputs[0] = d_A3_4;   // matches reduce["A3_4_red"].instr[0].name == "A3_4"
      compute_reduceN(d_reduce[idx_reduce_A3_4_red], inputs, d_A3_4_red);
    }

    // ---- L3 Stage 12: FunctionReduceN rho3 = reduce_n({A3_red,A3_2_red,A3_3_red,A3_4_red})
    // (multi-instruction + elem-dependent + only_invar). inputs[j] matched to
    // instruction j's NAME via rho3_input_role (resolved in init_style). There is
    // NO eq3 — rho3 is the ONLY terminal reduce of the L3 collectors. ----
    {
      t_nn_3d l3_views[4] = {d_A3_red, d_A3_2_red, d_A3_3_red, d_A3_4_red};
      t_nn_3d inputs[GRACE3L_MAX_REDUCE_INSTR];
      const int n_instr = d_reduce[idx_reduce_rho3].n_instr;
      for (int j = 0; j < n_instr; j++) inputs[j] = l3_views[rho3_input_role[j]];
      compute_reduceN(d_reduce[idx_reduce_rho3], inputs, d_rho3);
    }

    // ---- UQ basis-RP (Phase 3, final): accumulate the rho3 (L3) scalar basis into
    // d_uq_z (d_A3_red/d_A3_2_red/d_A3_3_red/d_A3_4_red are current here), then
    // L2-normalize, append density channels [full, rho1, rho2, rho3], and run the
    // GMM (nearest centroid, Mahalanobis sigma, gamma). Forward only. Placed here so
    // it reads the L3 forward views before the L3 backward (Stage 15) begins. ----
    if (any_uq)
      Kokkos::parallel_for("ComputeUQ",
          Kokkos::RangePolicy<DeviceType, TagComputeUQ>(0, chunk_size), *this);

    // ---- L3 Stage 13: InvariantLayerRMSNorm rho3_norm = invariant_rms_norm(rho3)
    // (FULL branch, same as rho2_norm — data-driven off scale_len). END of the
    // forward feature chain before readout. ----
    compute_invariant_rms_norm(d_rmsnorm[idx_rmsnorm_rho3_norm], d_rho3, d_rho3_norm);

    // ---- L3 Stage 14 (Task 2.11): Readout — 3-density sum -> energy MLP
    // (31->64 silu ->1) + lin skip + shift -> per-atom energy. rho3_norm is
    // read chunk-local (ii) since it was just written above for THIS chunk;
    // rho1_norm/rho2_norm are read from the natom buffers d_rho{1,2}_norm_full
    // (ii+chunk_offset, Task 3.5a) because their chunk-scratch originals
    // (d_rho1_norm/d_rho2_norm, written back in P1/P2) hold only the LAST
    // chunk under multi-chunk by the time P3 runs. ----
    {
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>(chunk_size, 1, 1);
      Kokkos::parallel_for("ComputeMLPEnergy", policy, *this);
    }

    if (eflag_global) {
      auto atom_energies = d_e_atom;
      KK_ACC_FLOAT energy_partial = 0.0;
      Kokkos::parallel_reduce("ComputeEnergySum", chunk_size,
          KOKKOS_LAMBDA(const int ii, KK_ACC_FLOAT& update) {
          update += atom_energies(ii);
      }, energy_partial);
      eng_vdwl += energy_partial;
    }

    if (eflag_atom) {
      auto e_chunk = d_e_atom;
      auto e_lammps = d_eatom;
      auto il = d_ilist;
      int co = chunk_offset;
      Kokkos::parallel_for("WriteEatom", Kokkos::RangePolicy<DeviceType>(0, chunk_size),
        KOKKOS_LAMBDA(const int ii) {
          e_lammps(il[ii + co]) += e_chunk(ii);
        });
    }

    // ---- L3 Stage 15 (Task 3.1): readout + Layer-3 backward (adjoints -> A3
    // basis + the L3 message adjoint d_eq2_norm_adj). Skipped in the energy-only
    // fast path (forces not requested); otherwise runs per chunk, reusing this
    // chunk's forward buffers. Zero every ACCUMULATING adjoint buffer, then walk
    // the L3 DAG in reverse topological order, ACCUMULATING (+=) where a tensor
    // fans out to multiple consumers (see task-3.1 DAG). ----
    if (!do_energy_only) {
      Kokkos::deep_copy(d_A3_red_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_2_red_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_3_red_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_4_red_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_2a_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_2b_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_2_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_3_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_4_adj, NNScalar(0.0));
      Kokkos::deep_copy(d_A3_adj, NNScalar(0.0));
      // NOTE (Task 3.5b): d_eq2_norm_adj is zeroed ONCE before the P3 loop (it
      // is [nall] and accumulated across chunks) — NOT here per chunk.

      // (a) readout bwd -> d_rho{1,2,3}_norm_adj (all EQUAL; GATE 1 = rho3_norm).
      readout_bwd();

      // (b) rho3_norm bwd (InvariantLayerRMSNorm, FULL branch) -> d_rho3_adj.
      compute_invariant_rms_norm_bwd(d_rmsnorm[idx_rmsnorm_rho3_norm],
          d_rho3, d_rho3_norm, d_rho3_norm_adj, d_rho3_adj);

      // (c) rho3 reduce bwd (only_invar, elem-dep, 4 instr) -> SEEDS the 4 *_red adj.
      {
        t_nn_3d l3_adj[4] = {d_A3_red_adj, d_A3_2_red_adj, d_A3_3_red_adj, d_A3_4_red_adj};
        t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR];
        const int n_instr = d_reduce[idx_reduce_rho3].n_instr;
        for (int j = 0; j < n_instr; j++) adj_in[j] = l3_adj[rho3_input_role[j]];
        compute_reduceN_bwd(d_reduce[idx_reduce_rho3], d_rho3_adj, adj_in);
      }

      // (d) A3_4_red -> A3_4.
      { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A3_4_adj;
        compute_reduceN_bwd(d_reduce[idx_reduce_A3_4_red], d_A3_4_red_adj, adj_in); }
      // (e) A3_4 = cp_l(A3_2b, A3_2b) self-product -> d_A3_2b_adj (both operands).
      compute_cp_l_bwd(d_prod[idx_prod_A3_4], d_A3_2b, d_A3_2b, d_A3_4_adj,
                       d_A3_2b_adj, d_A3_2b_adj);
      // (f) A3_2b = fc(A3_2_red, A3_red) -> d_A3_2_red_adj +=, d_A3_red_adj +=.
      compute_fc_bwd(d_fc[idx_fc_A3_2b], d_A3_2b_adj, d_A3_2_red_adj, d_A3_red_adj);
      // (g) A3_3_red -> A3_3.
      { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A3_3_adj;
        compute_reduceN_bwd(d_reduce[idx_reduce_A3_3_red], d_A3_3_red_adj, adj_in); }
      // (h) A3_3 = cp_l(A3_2a, A3_red) -> d_A3_2a_adj +=, d_A3_red_adj +=.
      compute_cp_l_bwd(d_prod[idx_prod_A3_3], d_A3_2a, d_A3_red, d_A3_3_adj,
                       d_A3_2a_adj, d_A3_red_adj);
      // (i) A3_2a = fc(A3_2_red, A3_red) -> d_A3_2_red_adj +=, d_A3_red_adj +=.
      compute_fc_bwd(d_fc[idx_fc_A3_2a], d_A3_2a_adj, d_A3_2_red_adj, d_A3_red_adj);
      // (j) A3_2_red -> A3_2  (d_A3_2_red_adj is now COMPLETE).
      { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A3_2_adj;
        compute_reduceN_bwd(d_reduce[idx_reduce_A3_2_red], d_A3_2_red_adj, adj_in); }
      // (k) A3_2 = cp_l(A3_red, A3_red) self-product -> d_A3_red_adj += (both operands).
      compute_cp_l_bwd(d_prod[idx_prod_A3_2], d_A3_red, d_A3_red, d_A3_2_adj,
                       d_A3_red_adj, d_A3_red_adj);
      // (l) A3_red -> A3  (d_A3_red_adj is now COMPLETE) => d_A3_adj  [GATE 2].
      { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A3_adj;
        compute_reduceN_bwd(d_reduce[idx_reduce_A3_red], d_A3_red_adj, adj_in); }
      // (m) A3 equiv-SPBF bwd -> d_eq2_norm_adj (L3 message adj; NOT reverse-comm'd
      // here — its gate is Task 3.2, after reverse-comm completes ghost contribs).
      // Opt 4: compute adj_prod ONCE here; A3SPBFBwd AND the later ForceEquiv(A3)
      // both read d_adj_prod (not overwritten between them in P3).
      compute_adj_prod(idx_spbf_A3, d_A3_adj);
      compute_A3_spbf_bwd();

      // ---- Task 3.5b: fold the L3 geometry force into P3 (per chunk). d_A3_adj
      // is complete for THIS chunk; A3's radial + dR/dr (d_R1_nl/d_DR1_nl/d_denv)
      // are still resident from Stage 1's deriv run (Opt 13 — nothing in Stages
      // 2-14 writes them), so no radial recompute here. Zero the per-bond force
      // buffer, run the equiv force kernel, and scatter into f[]/virial. f[]
      // (atomic_add) and virial[] (+=) accumulate across the 3 layers/phases and
      // all chunks -> identical total. ----
      Kokkos::deep_copy(d_f_ij, GeomScalar(0.0));
      compute_force_equiv(idx_spbf_A3, d_A3_adj, /*ind_stage=*/2, d_spbf[idx_spbf_A3].n_lm_ind);
      scatter_forces_chunk();
    }

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  // ---- Energy-only fast path: everything below (PHASES 3.5/4/4.5/5) is pure
  // force computation; skip it all when forces are not requested. ----
  if (!do_energy_only) {

  // ============ PHASE 3.5 (Task 3.2): reverse-comm the L2 message adjoint ======
  // compute_A3_spbf_bwd scattered d_eq2_norm_adj into owned+ghost rows (the owned
  // PARTIAL); reverse_comm ADDS every ghost image back into its owning rank's row
  // so the OWNED rows are now COMPLETE (== golden adj_eq2_norm after the tag fold).
  // comm_stage=2 selects the parity-doubled eq2_norm buffer + dims (I2_n_funcs*
  // I2_n_out=3200/atom) in all 4 reverse pack/unpack methods. Mirrors the round-1
  // forward_comm branch structure (device vs legacy-host fallback). Do NOT re-zero
  // d_eq2_norm_adj — the scatter partial must survive into the reverse-comm. ----
  comm_stage = 2;
  if (lmp->kokkos->reverse_pair_comm_legacy) {
    Kokkos::deep_copy(h_eq2_norm_adj, d_eq2_norm_adj);
    comm->reverse_comm(this);
    Kokkos::deep_copy(d_eq2_norm_adj, h_eq2_norm_adj);
  } else {
    comm->reverse_comm(this);
  }

  // ============ PHASE 4 (Task 3.2/3.5b): Layer-2 backward (CHUNKED) ============
  // Each chunk RE-RUNS the L2 forward body (l2_forward_chunk) to repopulate the
  // chunk-scratch L2 forward intermediates (d_A2*, d_eq2, d_rho2, ...) + geometry
  // + A2's radial (d_R1_nl) for THIS chunk (correct because d_I_global is already
  // ghost-complete), gathers the reverse-comm-completed seeds, walks the L2 DAG
  // -> d_A2_adj [GATE 2] -> A2-SPBF bwd -> d_grad_I_global (= d_eq1_norm_adj,
  // owned partial), then FOLDS the L2 geometry force into f[]/virial.
  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  // Perf (single-chunk fast path): when chunksize >= inum there is exactly ONE
  // chunk, so the L2 forward chunk-scratch (d_A2*, d_eq2, d_rho2, d_rho2_norm) is
  // still RESIDENT from PHASE 2 — PHASE 3 (L3) writes only L3-specific buffers and
  // the shared radial scratch d_R1_nl. So the L2 forward RECOMPUTE is pure overhead
  // here; skip it and restore only d_R1_nl/d_DR1_nl/d_denv (the L3-clobbered shared
  // scratch) via the reordered compute_mlp_radial_deriv below. Multi-chunk (a later
  // chunk overwrites the earlier chunk's scratch) still recomputes.
  const bool single_chunk = (chunk_size >= inum);
  // Task 3.5b: zero the round-1 message adjoint ONCE before the loop ([nall]
  // buffer, atomic_add-scattered by compute_A2_spbf_bwd across ALL chunks).
  Kokkos::deep_copy(d_grad_I_global, NNScalar(0.0));
  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    // (i) recompute the L2 forward intermediates + geometry + A2 radial for THIS
    // chunk (leaves d_R1_nl = A2's radial; and rho2_norm bwd reads this chunk's
    // d_rho2/d_rho2_norm). SKIPPED for single_chunk — buffers already resident.
    if (!single_chunk) l2_forward_chunk();

    // (ii) Gather the reverse-comm-completed L2 message adjoint from LAMMPS-local
    // rows (d_eq2_norm_adj[d_ilist(ii+co)]) into chunk-local rows (comm'd buffer
    // -> keyed by the LAMMPS local index).
    {
      auto src = d_eq2_norm_adj;
      auto dst = d_eq2_norm_adj_local;
      auto il = d_ilist;
      const int co = chunk_offset;
      const int nf = I2_n_funcs, no = I2_n_out, cs = chunk_size;
      Kokkos::parallel_for("GatherEq2NormAdj",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nf, no}),
        KOKKOS_LAMBDA(const int ii, const int f, const int k) {
          dst(ii, f, k) = src(il(ii + co), f, k);
        });
    }
    // Gather the natom rho2 readout-adjoint seed (dense list-order ii+co; NOT
    // d_ilist — this buffer is local-only) into the chunk-local _local buffer
    // that the shared TagRMSNormBwd reads.
    {
      auto src = d_rho2_norm_adj;
      auto dst = d_rho2_norm_adj_local;
      const int co = chunk_offset;
      const int cs = chunk_size;
      const int nc = (int) dst.extent(1), nl = (int) dst.extent(2);
      Kokkos::parallel_for("GatherRho2NormAdj",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nc, nl}),
        KOKKOS_LAMBDA(const int ii, const int c, const int lm) {
          dst(ii, c, lm) = src(ii + co, c, lm);
        });
    }

    // (iii) Zero every ACCUMULATING adjoint buffer (the 4 *_red seeds, the
    // product/FC adjoints, d_A2_adj written by reduceN_bwd +=). d_eq2_adj /
    // d_rho2_adj are ASSIGNED (single consumer). d_grad_I_global is zeroed ONCE
    // before the loop (above), NOT here.
    Kokkos::deep_copy(d_A2_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_2_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_3_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_4_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_2a_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_2b_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_2_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_3_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_4_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A2_adj, NNScalar(0.0));

    // (iv) L2 backward body (verbatim).
    // (a) eq2_norm bwd (EquivariantRMSNorm VJP) -> d_eq2_adj (ASSIGN).
    compute_equiv_rms_norm_bwd(d_eqnorm[idx_eqnorm_eq2_norm], d_eq2,
        d_eq2_norm_adj_local, d_eq2_adj);
    // (b) eq2 reduce bwd -> SEED the 4 *_red adj (by eq2_input_role, ACCUMULATE).
    {
      t_nn_3d l2_adj[4] = {d_A2_red_adj, d_A2_2_red_adj, d_A2_3_red_adj, d_A2_4_red_adj};
      t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR];
      const int n_instr = d_reduce[idx_reduce_eq2].n_instr;
      for (int j = 0; j < n_instr; j++) adj_in[j] = l2_adj[eq2_input_role[j]];
      compute_reduceN_bwd(d_reduce[idx_reduce_eq2], d_eq2_adj, adj_in);
    }
    // (c) rho2_norm bwd (InvariantLayerRMSNorm, FULL) -> d_rho2_adj (ASSIGN).
    //     seed = d_rho2_norm_adj_local (gathered from the natom d_rho2_norm_adj).
    compute_invariant_rms_norm_bwd(d_rmsnorm[idx_rmsnorm_rho2_norm],
        d_rho2, d_rho2_norm, d_rho2_norm_adj_local, d_rho2_adj);
    // (d) rho2 reduce bwd -> ACCUMULATE into the SAME 4 *_red adj (by rho2_input_role).
    {
      t_nn_3d l2_adj[4] = {d_A2_red_adj, d_A2_2_red_adj, d_A2_3_red_adj, d_A2_4_red_adj};
      t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR];
      const int n_instr = d_reduce[idx_reduce_rho2].n_instr;
      for (int j = 0; j < n_instr; j++) adj_in[j] = l2_adj[rho2_input_role[j]];
      compute_reduceN_bwd(d_reduce[idx_reduce_rho2], d_rho2_adj, adj_in);
    }

    // ---- L2 DAG walk (reverse topological; mirrors the L3 walk exactly) ----
    // (e) A2_4_red -> A2_4.
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A2_4_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A2_4_red], d_A2_4_red_adj, adj_in); }
    // (f) A2_4 = cp_l(A2_2b, A2_2b) self-product -> d_A2_2b_adj (both operands).
    compute_cp_l_bwd(d_prod[idx_prod_A2_4], d_A2_2b, d_A2_2b, d_A2_4_adj,
                     d_A2_2b_adj, d_A2_2b_adj);
    // (g) A2_2b = fc(A2_2_red, A2_red) -> d_A2_2_red_adj +=, d_A2_red_adj +=.
    compute_fc_bwd(d_fc[idx_fc_A2_2b], d_A2_2b_adj, d_A2_2_red_adj, d_A2_red_adj);
    // (h) A2_3_red -> A2_3.
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A2_3_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A2_3_red], d_A2_3_red_adj, adj_in); }
    // (i) A2_3 = cp_l(A2_2a, A2_red) -> d_A2_2a_adj +=, d_A2_red_adj +=.
    compute_cp_l_bwd(d_prod[idx_prod_A2_3], d_A2_2a, d_A2_red, d_A2_3_adj,
                     d_A2_2a_adj, d_A2_red_adj);
    // (j) A2_2a = fc(A2_2_red, A2_red) -> d_A2_2_red_adj +=, d_A2_red_adj +=.
    compute_fc_bwd(d_fc[idx_fc_A2_2a], d_A2_2a_adj, d_A2_2_red_adj, d_A2_red_adj);
    // (k) A2_2_red -> A2_2  (d_A2_2_red_adj now COMPLETE: eq2+rho2 seeds + fc_2b/2a left).
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A2_2_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A2_2_red], d_A2_2_red_adj, adj_in); }
    // (l) A2_2 = cp_l(A2_red, A2_red) self-product -> d_A2_red_adj += (both operands).
    compute_cp_l_bwd(d_prod[idx_prod_A2_2], d_A2_red, d_A2_red, d_A2_2_adj,
                     d_A2_red_adj, d_A2_red_adj);
    // (m) A2_red -> A2  (d_A2_red_adj now COMPLETE) => d_A2_adj  [GATE 2].
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A2_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A2_red], d_A2_red_adj, adj_in); }

    // (n) Restore A2's radial d_R1_nl/d_DR1_nl (+ d_denv) for THIS chunk.
    // compute_mlp_radial_deriv recomputes the radial straight from geometry
    // (d_rnorms), so it repairs the shared d_R1_nl that L3 clobbered — the ONLY
    // thing the single_chunk path needs restored. Reordered BEFORE
    // compute_A2_spbf_bwd (which reads d_R1_nl = A2's radial) and
    // compute_force_equiv. (Multi-chunk: l2_forward_chunk already set d_R1_nl = A2,
    // so this recompute is identical and the reorder is a no-op there.)
    radial_mlp_spbf = idx_spbf_A2;
    compute_mlp_radial_deriv();

    // (o) A2 equiv-SPBF bwd -> d_grad_I_global (= d_eq1_norm_adj, owned partial).
    // Opt 4: compute adj_prod ONCE here; A2SPBFBwd AND the later ForceEquiv(A2)
    // both read d_adj_prod (not overwritten between them in P4).
    compute_adj_prod(idx_spbf_A2, d_A2_adj);
    compute_A2_spbf_bwd();

    // (v) fold the L2 geometry force into f[]/virial for THIS chunk.
    Kokkos::deep_copy(d_f_ij, GeomScalar(0.0));
    compute_force_equiv(idx_spbf_A2, d_A2_adj, /*ind_stage=*/1, d_spbf[idx_spbf_A2].n_lm_ind);
    scatter_forces_chunk();

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  // ============ PHASE 4.5 (Task 3.3): reverse-comm the L1 message adjoint ======
  // compute_A2_spbf_bwd scattered d_grad_I_global (= d_eq1_norm_adj) into owned+
  // ghost rows (owned PARTIAL); reverse_comm (comm_stage=1 -> the round-1 grad_I
  // branch in all 4 reverse pack/unpack methods) ADDS every ghost image back into
  // its owning rank's row, so the OWNED rows are COMPLETE (== golden adj_eq1_norm).
  // Do NOT re-zero d_grad_I_global — the scatter partial must survive. ----
  comm_stage = 1;
  if (lmp->kokkos->reverse_pair_comm_legacy) {
    Kokkos::deep_copy(h_grad_I_global, d_grad_I_global);
    comm->reverse_comm(this);
    Kokkos::deep_copy(d_grad_I_global, h_grad_I_global);
  } else {
    comm->reverse_comm(this);
  }

  // ============ PHASE 5 (Task 3.3/3.5b): Layer-1 backward (CHUNKED) ============
  // Each chunk RE-RUNS the L1 forward body (l1_forward_chunk) to repopulate the
  // chunk-scratch L1 forward buffers (d_A1, d_A1_2, ..., d_eq1, d_rho1,
  // d_rho1_norm) + geometry + A1's radial (d_R1_nl) for THIS chunk. UNLIKE L2/L3,
  // the L1 chain consumes RAW A1 (there is no A1_red), so d_A1_adj is the DIRECT
  // fan-out target. Dual-seed the 4 l1_input adjoints from BOTH eq1 (equiv) and
  // rho1 (invariant), walk the L1 DAG -> d_A1_adj [GATE 2], then FOLD the L1
  // geometry force into f[]/virial. (P5 previously had NO chunk loop.)
  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  // Perf (single-chunk fast path): the L1 forward chunk-scratch (d_A1*, d_eq1,
  // d_rho1, d_rho1_norm) is still RESIDENT from PHASE 1 — nothing between P1 and
  // here writes the L1-specific buffers (P2/P3/P4 use L2/L3 buffers + the shared
  // radial scratch). Skip the L1 forward recompute; compute_mlp_radial_deriv below
  // (already before compute_force_L1) restores the shared d_R1_nl/d_DR1_nl/d_denv.
  // Reuses `single_chunk` from PHASE 4 (same condition: chunksize >= inum).
  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    // (i) recompute the L1 forward intermediates + geometry + A1 radial for THIS
    // chunk. SKIPPED for single_chunk — buffers already resident from PHASE 1.
    if (!single_chunk) l1_forward_chunk();

    // (ii) Gather the reverse-comm-completed round-1 message adjoint from
    // LAMMPS-local rows (d_grad_I_global[d_ilist(ii+co)]) into chunk-local rows.
    {
      auto src = d_grad_I_global;
      auto dst = d_eq1_norm_adj_local;
      auto il = d_ilist;
      const int co = chunk_offset;
      const int nf = I_n_funcs, no = I_n_out, cs = chunk_size;
      Kokkos::parallel_for("GatherEq1NormAdj",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nf, no}),
        KOKKOS_LAMBDA(const int ii, const int f, const int kk) {
          dst(ii, f, kk) = src(il(ii + co), f, kk);
        });
    }
    // Gather the natom rho1 readout-adjoint seed (dense list-order ii+co; NOT
    // d_ilist) into the chunk-local _local buffer the shared TagRMSNormBwd reads.
    {
      auto src = d_rho1_norm_adj;
      auto dst = d_rho1_norm_adj_local;
      const int co = chunk_offset;
      const int cs = chunk_size;
      const int nc = (int) dst.extent(1), nl = (int) dst.extent(2);
      Kokkos::parallel_for("GatherRho1NormAdj",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nc, nl}),
        KOKKOS_LAMBDA(const int ii, const int c, const int lm) {
          dst(ii, c, lm) = src(ii + co, c, lm);
        });
    }

    // Zero every ACCUMULATING L1 adjoint buffer (bwd routines +=). d_eq1_adj /
    // d_rho1_adj are ASSIGNED by the rms-norm bwds (single consumer).
    Kokkos::deep_copy(d_A1_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_2_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_2_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_2a_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_3_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_3_red_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_2b_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_4_adj, NNScalar(0.0));
    Kokkos::deep_copy(d_A1_4_red_adj, NNScalar(0.0));

    // (1) eq1_norm bwd (EquivariantRMSNorm VJP) -> d_eq1_adj (ASSIGN).
    compute_equiv_rms_norm_bwd(d_eqnorm[idx_eqnorm_eq1_norm], d_eq1,
        d_eq1_norm_adj_local, d_eq1_adj);
    // (2) eq1 reduce bwd -> SEED the 4 l1_input adj (by eq1_input_role; role 0 = raw A1).
    {
      t_nn_3d l1_adj[4] = {d_A1_adj, d_A1_2_red_adj, d_A1_3_red_adj, d_A1_4_red_adj};
      t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR];
      const int n_instr = d_reduce[idx_reduce_eq1].n_instr;
      for (int j = 0; j < n_instr; j++) adj_in[j] = l1_adj[eq1_input_role[j]];
      compute_reduceN_bwd(d_reduce[idx_reduce_eq1], d_eq1_adj, adj_in);
    }
    // (3) rho1_norm bwd (InvariantLayerRMSNorm, only_nonlin) -> d_rho1_adj (ASSIGN).
    //     seed = d_rho1_norm_adj_local (gathered from the natom d_rho1_norm_adj).
    compute_invariant_rms_norm_bwd(d_rmsnorm[idx_rmsnorm_rho1_norm],
        d_rho1, d_rho1_norm, d_rho1_norm_adj_local, d_rho1_adj);
    // (4) rho1 reduce bwd -> ACCUMULATE into the SAME 4 (by rho1_input_role).
    {
      t_nn_3d l1_adj[4] = {d_A1_adj, d_A1_2_red_adj, d_A1_3_red_adj, d_A1_4_red_adj};
      t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR];
      const int n_instr = d_reduce[idx_reduce_rho1].n_instr;
      for (int j = 0; j < n_instr; j++) adj_in[j] = l1_adj[rho1_input_role[j]];
      compute_reduceN_bwd(d_reduce[idx_reduce_rho1], d_rho1_adj, adj_in);
    }

    // ---- L1 DAG walk (reverse topological) ----
    // (d) A1_4_red -> A1_4.
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A1_4_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A1_4_red], d_A1_4_red_adj, adj_in); }
    // (e) A1_4 = cp_l(A1_2b, A1_2b) self-product -> d_A1_2b_adj (both operands).
    compute_cp_l_bwd(d_prod[idx_prod_A1_4], d_A1_2b, d_A1_2b, d_A1_4_adj,
                     d_A1_2b_adj, d_A1_2b_adj);
    // (f) A1_2b = fc(A1_2_red, A1) -> d_A1_2_red_adj +=, d_A1_adj += (right = raw A1).
    compute_fc_bwd(d_fc[idx_fc_A1_2b], d_A1_2b_adj, d_A1_2_red_adj, d_A1_adj);
    // (g) A1_3_red -> A1_3.
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A1_3_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A1_3_red], d_A1_3_red_adj, adj_in); }
    // (h) A1_3 = cp_l(A1_2a, A1) -> d_A1_2a_adj +=, d_A1_adj += (right = raw A1).
    compute_cp_l_bwd(d_prod[idx_prod_A1_3], d_A1_2a, d_A1, d_A1_3_adj,
                     d_A1_2a_adj, d_A1_adj);
    // (i) A1_2a = fc(A1_2_red, A1) -> d_A1_2_red_adj +=, d_A1_adj += (right = raw A1).
    compute_fc_bwd(d_fc[idx_fc_A1_2a], d_A1_2a_adj, d_A1_2_red_adj, d_A1_adj);
    // (j) A1_2_red -> A1_2  (d_A1_2_red_adj now COMPLETE: eq1+rho1 seeds + fc_2b/2a left).
    { t_nn_3d adj_in[GRACE3L_MAX_REDUCE_INSTR]; adj_in[0] = d_A1_2_adj;
      compute_reduceN_bwd(d_reduce[idx_reduce_A1_2_red], d_A1_2_red_adj, adj_in); }
    // (k) A1_2 = cp_l(A1, A1) self-product -> d_A1_adj += (both operands).
    //     d_A1_adj now COMPLETE: eq1[0]+rho1[0] + fc_2b(R) + cp_l_3(R) + fc_2a(R) + cp_l_2(both).
    compute_cp_l_bwd(d_prod[0], d_A1, d_A1, d_A1_2_adj, d_A1_adj, d_A1_adj);

    // (v) fold the L1 geometry force into f[]/virial for THIS chunk. compute_force_L1
    // needs d_R1_nl/d_DR1_nl = A1's radial + dR/dr, so recompute the radial deriv
    // (l1_forward_chunk left d_R1_nl = A1's radial value, but not d_DR1_nl).
    radial_mlp_spbf = 0;
    compute_mlp_radial_deriv();
    Kokkos::deep_copy(d_f_ij, GeomScalar(0.0));
    compute_force_L1();
    scatter_forces_chunk();

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  } // end !do_energy_only (backward PHASES 3.5/4/4.5/5)

  // NOTE (Task 3.5b): the old single-shot PHASE 6 (all 3 layers' geometry force
  // assembled into one d_f_ij) and PHASE 6b (single scatter of d_f_ij -> f[] +
  // virial) are DELETED. Each layer's geometry force is now folded into its
  // owning backward phase — L3 -> PHASE 3, L2 -> PHASE 4, L1 -> PHASE 5 — and
  // scattered per chunk via scatter_forces_chunk(). f[] (atomic_add) and virial[]
  // (+=) accumulate over all 3 layers and all chunks to the identical total.

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
  if (cvflag_atom) {
    k_cvatom.template modify<DeviceType>();
    k_cvatom.sync_host();
  }

  // ---- Copy per-atom UQ results device -> host arrays (indexed by local atom) ----
  if (any_uq) {
    const int nlocal = atom->nlocal;
    auto copy_out = [&](double *host, t_uq_1d &dview) {
      Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> hwrap(host, nlocal);
      Kokkos::deep_copy(hwrap, Kokkos::subview(dview, Kokkos::make_pair(0, nlocal)));
    };
    copy_out(gamma, d_gamma);
    copy_out(atomic_sigma, d_sigma);
    copy_out(gmm_cluster, d_gmm_cluster);
  }

  atomKK->modified(execution_space, F_MASK);
  copymode = 0;
}

// ======================================================================
// Utility functions (COPIED VERBATIM from the validated GRACE-2L KOKKOS
// style, module renamed).
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<class TagStyle>
void PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::check_team_size_for(int inum_in, int &team_size, int vector_length) {
  int team_size_max;
  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum_in,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<typename scratch_type>
KOKKOS_INLINE_FUNCTION
int PairGRACE3LKokkos<DeviceType, NNScalarT, GeomScalarT>::scratch_size_helper(int values_per_team) const {
  typedef Kokkos::View<scratch_type*, typename DeviceType::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ScratchViewType;
  return ScratchViewType::shmem_size(values_per_team);
}

// ======================================================================
// Template instantiations
// ======================================================================

namespace LAMMPS_NS {
// Mixed (NN=float, geometry=double) — default
template class PairGRACE3LKokkos<LMPDeviceType, float, double>;

// FP32 (all float)
template class PairGRACE3LKokkos<LMPDeviceType, float, float>;
}
