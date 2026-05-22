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
   GRACE-1L KOKKOS implementation
------------------------------------------------------------------------- */

#include "pair_grace_1l_cpu_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"
#include "neigh_list_kokkos.h"

#include "cnpy/cnpy.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

// Spherical harmonics constants
static constexpr double Y00 = 1.0;
static constexpr double sq3 = 1.7320508075688772935;
static constexpr double sq3o2 = 1.2247448713915890491;
static constexpr double sq2 = 1.4142135623730950488;

// ======================================================================
// FindMaxNumNeighs helper
// ======================================================================
template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;
  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& maxneigh) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (maxneigh < num_neighs) maxneigh = num_neighs;
  }
};

// ======================================================================
// GRACE1LCPUModel: load weights from .npz
// ======================================================================

static std::vector<double> npz_get_double(const cnpy::npz_t &npz, const std::string &key) {
  auto it = npz.find(key);
  if (it == npz.end()) return {};
  const auto &arr = it->second;
  const double *p = arr.data<double>();
  return std::vector<double>(p, p + arr.num_vals);
}

static std::vector<int> npz_get_int(const cnpy::npz_t &npz, const std::string &key) {
  auto it = npz.find(key);
  if (it == npz.end()) return {};
  const auto &arr = it->second;
  // Handle both int32 and int64
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

void GRACE1LCPUModel::load(const std::string &filepath) {
  cnpy::npz_t npz = cnpy::npz_load(filepath);

  // Metadata
  n_elements = npz_get_int_scalar(npz, "n_elements", 89);
  embedding_size = npz_get_int_scalar(npz, "embedding_size", 128);
  lmax = npz_get_int_scalar(npz, "lmax", 4);
  nradbase = npz_get_int_scalar(npz, "nradbase", 10);
  nradmax = npz_get_int_scalar(npz, "nradmax", 32);
  radial_basis_p = npz_get_int_scalar(npz, "radial_basis_p", 16);
  rcut = npz_get_double_scalar(npz, "rcut", 6.0);
  has_bond_specific_cutoff = npz_get_int_scalar(npz, "has_bond_specific_cutoff", 0) != 0;
  bond_cutoff_map = npz_get_double(npz, "bond_cutoff_map");
  // If no bond_cutoff_map provided, fill uniformly with rcut
  if (bond_cutoff_map.empty()) {
    bond_cutoff_map.resize(n_elements * n_elements, rcut);
  }
  mlp_rad_has_chem_emb = npz_get_int_scalar(npz, "mlp_rad_has_chem_emb", 0) != 0;

  // Chemical embedding
  chem_embedding = npz_get_double(npz, "chem_embedding");

  // A indicator transform
  A_lin_transform_W = npz_get_double(npz, "A_lin_transform_W");
  A_lin_transform_norm = npz_get_double_scalar(npz, "A_lin_transform_norm", 1.0);
  A_inv_avg_n_neigh = npz_get_double_scalar(npz, "A_inv_avg_n_neigh", 1.0);

  // MLP radial layers
  mlp_rad_n_layers = npz_get_int_scalar(npz, "mlp_rad_n_layers", 3);
  mlp_rad_layers.resize(mlp_rad_n_layers);
  for (int i = 0; i < mlp_rad_n_layers; i++) {
    auto &layer = mlp_rad_layers[i];
    layer.W = npz_get_double(npz, "mlp_rad_W" + std::to_string(i));
    layer.norm = npz_get_double_scalar(npz, "mlp_rad_norm" + std::to_string(i), 1.0);
    // Infer dimensions from weight matrix
    auto it = npz.find("mlp_rad_W" + std::to_string(i));
    if (it != npz.end() && it->second.shape.size() == 2) {
      layer.n_in = it->second.shape[0];
      layer.n_out = it->second.shape[1];
    }
  }

  // Helper to load FC weights
  auto load_fc = [&](const std::string &prefix, FCWeights &fc) {
    fc.w_left = npz_get_double(npz, prefix + "_w_left");
    fc.w_right = npz_get_double(npz, prefix + "_w_right");
    fc.w_tile_left = npz_get_int(npz, prefix + "_w_tile_left");
    fc.w_tile_right = npz_get_int(npz, prefix + "_w_tile_right");
    fc.collect_to = npz_get_int(npz, prefix + "_collect_to");
    fc.collect_from = npz_get_int(npz, prefix + "_collect_from");
    fc.norm_out_factor = npz_get_double(npz, prefix + "_norm_out_factor");
    fc.norm_left = npz_get_double_scalar(npz, prefix + "_norm_left", 1.0);
    fc.norm_right = npz_get_double_scalar(npz, prefix + "_norm_right", 1.0);
    fc.n_out = npz_get_int_scalar(npz, prefix + "_n_out", 64);
    fc.left_coefs = npz_get_int_scalar(npz, prefix + "_left_coefs", 1) != 0;
    fc.n_funcs_left = (int)fc.w_tile_left.size();
    fc.n_funcs_right = (int)fc.collect_from.size();
    // Infer w_shape from weight array
    auto it = npz.find(prefix + "_w_left");
    if (it != npz.end() && it->second.shape.size() == 3) {
      fc.w_shape_left = it->second.shape[2];
      fc.w_shape_right = fc.w_shape_left; // usually same
    }
    it = npz.find(prefix + "_w_right");
    if (it != npz.end() && it->second.shape.size() == 3) {
      fc.w_shape_right = it->second.shape[2];
    }
  };

  load_fc("fc_A1", fc_A1);
  load_fc("fc_AA1", fc_AA1);
  load_fc("fc_AA2", fc_AA2);

  // Helper to load CG product
  auto load_prod = [&](const std::string &prefix, CGProduct &prod) {
    prod.left_ind = npz_get_int(npz, prefix + "_left_ind");
    prod.right_ind = npz_get_int(npz, prefix + "_right_ind");
    prod.m_sum_ind = npz_get_int(npz, prefix + "_m_sum_ind");
    prod.cg_coeff = npz_get_double(npz, prefix + "_cg_coeff");
    prod.n_output_funcs = npz_get_int_scalar(npz, prefix + "_n_output_funcs", 0);
    prod.n_cg_terms = npz_get_int_scalar(npz, prefix + "_n_cg_terms", 0);
  };

  load_prod("prod_AA", prod_AA);
  load_prod("prod_AAA", prod_AAA);
  load_prod("prod_AAAA", prod_AAAA);

  // FunctionReduceN
  rho_n_out = npz_get_int_scalar(npz, "rho_n_out", 17);
  rho_is_elem_dependent = npz_get_int_scalar(npz, "rho_is_central_atom_type_dependent", 1) != 0;

  auto load_reduce = [&](const std::string &instr_name, ReduceWeights &rw) {
    rw.W = npz_get_double(npz, "rho_reduce_" + instr_name + "_W");
    rw.norm = npz_get_double_scalar(npz, "rho_reduce_" + instr_name + "_norm", 1.0);
    rw.collect_ind = npz_get_int(npz, "rho_collect_ind_" + instr_name);
    rw.w_l_tile = npz_get_int(npz, "rho_w_l_tile_" + instr_name);
    rw.n_in = npz_get_int_scalar(npz, "rho_n_in_" + instr_name, 0);
    rw.w_shape = npz_get_int_scalar(npz, "rho_w_shape_" + instr_name, 0);
  };

  load_reduce("A", reduce_A);
  load_reduce("AA", reduce_AA);
  load_reduce("AAA", reduce_AAA);
  load_reduce("AAAA", reduce_AAAA);

  // Energy MLP
  int energy_n_layers = npz_get_int_scalar(npz, "energy_mlp_n_layers", 2);
  energy_mlp_layers.resize(energy_n_layers);
  for (int i = 0; i < energy_n_layers; i++) {
    auto &layer = energy_mlp_layers[i];
    std::string wkey = "energy_mlp_W" + std::to_string(i);
    std::string nkey = "energy_mlp_norm" + std::to_string(i);
    layer.W = npz_get_double(npz, wkey);
    layer.norm = npz_get_double_scalar(npz, nkey, 1.0);
    auto wit = npz.find(wkey);
    if (wit != npz.end() && wit->second.shape.size() == 2) {
      layer.n_in = wit->second.shape[0];
      layer.n_out = wit->second.shape[1];
    }
  }
  energy_mlp_activation = npz_get_int_scalar(npz, "energy_mlp_activation", 0);  // 0=silu

  // Shifts
  shift_values = npz_get_double(npz, "shift_values");
  output_scale = npz_get_double_scalar(npz, "output_scale", 1.0);

  // Element names (stored as S2 = 2-byte ASCII strings)
  auto it_en = npz.find("element_names");
  if (it_en != npz.end()) {
    const auto &arr = it_en->second;
    int n = arr.shape[0];
    int wsize = arr.word_size;  // 2 for S2
    const char *raw = arr.data<char>();
    element_names.resize(n);
    for (int i = 0; i < n; i++) {
      std::string s(raw + i * wsize, wsize);
      // Trim trailing null bytes
      while (!s.empty() && s.back() == '\0') s.pop_back();
      element_names[i] = s;
    }
  }
}

// ======================================================================
// PairGRACE1LCPUKokkos: Constructor / Destructor
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::PairGRACE1LCPUKokkos(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  host_flag = (execution_space == Host);

  grace_model = new GRACE1LCPUModel();
  chunksize = 4096;
  no_virial_fdotr_compute = 1;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::~PairGRACE1LCPUKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_eatom, eatom);
  memoryKK->destroy_kokkos(k_vatom, vatom);
  memoryKK->destroy_kokkos(k_cvatom, cvatom);
  delete grace_model;
}

// ======================================================================
// settings: parse pair_style arguments
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::settings(int narg, char **arg)
{
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "chunksize") == 0) {
      chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "debug_no_energy_only_calc") == 0) {
      debug_no_energy_only_calc = true;
      iarg += 1;
    } else {
      error->all(FLERR, "Unknown pair_style grace/1l/cpu/kk keyword: {}", arg[iarg]);
    }
  }
}

// ======================================================================
// extract: expose internal pointers to fixes/computes
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void *PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::extract(const char *str, int & /*dim*/)
{
  if (strcmp(str, "compute_energy_only") == 0) return (void *) &flag_compute_energy_only;
  if (strcmp(str, "debug_no_energy_only_calc") == 0) return (void *) &debug_no_energy_only_calc;
  return nullptr;
}

// ======================================================================
// allocate
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  map = new int[n];

  MemKK::realloc_kokkos(d_map, "grace_1l:map", n);
  MemKK::realloc_kokkos(k_cutsq, "grace_1l:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();
}

// ======================================================================
// coeff: load model weights
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients");

  // pair_coeff * * weights.npz elem1 elem2 ...
  map_element2type(narg - 3, arg + 3);

  std::string weights_path = arg[2];
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE-1L/KK] Loading weights from {}\n", weights_path);

  grace_model->load(weights_path);

  if (comm->me == 0) {
    utils::logmesg(lmp, "[GRACE-1L/KK] n_elements={}, lmax={}, nradbase={}, nradmax={}, rcut={}, bond_specific_cutoff={}\n",
                   grace_model->n_elements, grace_model->lmax,
                   grace_model->nradbase, grace_model->nradmax, grace_model->rcut,
                   grace_model->has_bond_specific_cutoff ? "yes" : "no");
  }

  // Build element mapping: LAMMPS atom type -> model's internal element index
  // map_element2type sets map[i] = index into the user-provided element list (0-based)
  // We need to convert that to the model's internal element index using element_names
  auto h_map = Kokkos::create_mirror_view(d_map);
  int nuser = narg - 3;  // number of element names provided by user
  for (int i = 1; i <= atom->ntypes; i++) {
    if (map[i] < 0 || map[i] >= nuser) {
      h_map(i) = -1;
      continue;
    }
    // Get the user-provided element name for this type
    std::string user_elem = arg[3 + map[i]];
    // Find it in the model's element_names list
    int model_idx = -1;
    for (int j = 0; j < (int)grace_model->element_names.size(); j++) {
      if (grace_model->element_names[j] == user_elem) {
        model_idx = j;
        break;
      }
    }
    if (model_idx < 0)
      error->all(FLERR, "Element '{}' not found in GRACE-1L model element list", user_elem);
    h_map(i) = model_idx;
    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE-1L/KK] LAMMPS type {} -> {} (model index {})\n", i, user_elem, model_idx);
  }
  Kokkos::deep_copy(d_map, h_map);

  // Store host-side type -> model index mapping for init_one()
  h_type2model.resize(atom->ntypes + 1, -1);
  for (int i = 1; i <= atom->ntypes; i++)
    h_type2model[i] = h_map(i);
}

// ======================================================================
// init_style: copy weights to device
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace/1l/cpu/kk requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace/1l/cpu/kk requires newton pair on");

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  // No restriction on neighflag: pair style requested a full neighbor list above,
  // and the kernels iterate that list directly. Force application is per-bond on
  // the central (local) atom and its neighbor j; ghost-j writes are reverse-comm'd
  // back to the owner via standard pair-style force reverse_comm.

  // Copy scalar architecture params
  nelements = grace_model->n_elements;
  lmax = grace_model->lmax;
  nradmax = grace_model->nradmax;
  nradbase = grace_model->nradbase;
  embedding_size = grace_model->embedding_size;
  radial_basis_p = grace_model->radial_basis_p;
  rcut = grace_model->rcut;
  rho_n_out = grace_model->rho_n_out;
  fc_n_out = grace_model->fc_A1.n_out;  // 64

  n_funcs_A1 = grace_model->fc_A1.n_funcs_left;   // 25
  n_funcs_AA = grace_model->prod_AA.n_output_funcs; // 141
  n_funcs_AA1 = n_funcs_AA;  // FC AA1 has same # output funcs
  n_funcs_AAA = grace_model->prod_AAA.n_output_funcs; // 14
  n_funcs_AA2 = n_funcs_AA;
  n_funcs_AAAA = grace_model->prod_AAAA.n_output_funcs; // 69
  n_cg_AA = grace_model->prod_AA.n_cg_terms;
  n_cg_AAA = grace_model->prod_AAA.n_cg_terms;
  n_cg_AAAA = grace_model->prod_AAAA.n_cg_terms;

  // Precompute spherical harmonics coefficients
  MemKK::realloc_kokkos(d_idx_sph, "grace_1l:idx_sph", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(alm, "grace_1l:alm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(blm, "grace_1l:blm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(cl, "grace_1l:cl", lmax + 1);
  MemKK::realloc_kokkos(dl, "grace_1l:dl", lmax + 1);
  precompute_harmonics();

  // Copy all weights to device
  copy_weights_to_device();
}

// ======================================================================
// init_one: set cutoff for pair (i,j)
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
double PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::init_one(int i, int j)
{
  double rcut_ij = rcut;  // default: global cutoff
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
// precompute_harmonics: same as GRACE-FS
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::precompute_harmonics()
{
  auto h_idx_sph = Kokkos::create_mirror_view(d_idx_sph);
  auto h_alm = Kokkos::create_mirror_view(alm);
  auto h_blm = Kokkos::create_mirror_view(blm);
  auto h_cl = Kokkos::create_mirror_view(cl);
  auto h_dl = Kokkos::create_mirror_view(dl);

  Kokkos::deep_copy(h_idx_sph, -1);

  int idx_sph = 0;
  for (int m = 0; m <= lmax; m++) {
    const double msq = m * m;
    for (int l = m; l <= lmax; l++) {
      const int idx = l * (l + 1) + m;
      h_idx_sph(idx) = idx_sph;
      double a = 0.0, b = 0.0;
      if (l > 1 && m < l - 1) {
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
// copy_weights_to_device: transfer all model weights to Kokkos views
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::copy_weights_to_device()
{
  auto &m = *grace_model;

  // Helper: copy flat vector to 2D Kokkos view
  auto copy_2d = [](const std::vector<double> &src, auto &dst, int d0, int d1, const char *label) {
    MemKK::realloc_kokkos(dst, label, d0, d1);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        h(i, j) = src[i * d1 + j];
    Kokkos::deep_copy(dst, h);
  };

  auto copy_3d = [](const std::vector<double> &src, auto &dst, int d0, int d1, int d2, const char *label) {
    MemKK::realloc_kokkos(dst, label, d0, d1, d2);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        for (int k = 0; k < d2; k++)
          h(i, j, k) = src[(i * d1 + j) * d2 + k];
    Kokkos::deep_copy(dst, h);
  };

  auto copy_4d = [](const std::vector<double> &src, auto &dst, int d0, int d1, int d2, int d3, const char *label) {
    MemKK::realloc_kokkos(dst, label, d0, d1, d2, d3);
    auto h = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        for (int k = 0; k < d2; k++)
          for (int l = 0; l < d3; l++)
            h(i, j, k, l) = src[((i * d1 + j) * d2 + k) * d3 + l];
    Kokkos::deep_copy(dst, h);
  };

  auto copy_1d_int = [](const std::vector<int> &src, auto &dst, const char *label) {
    MemKK::realloc_kokkos(dst, label, src.size());
    auto h = Kokkos::create_mirror_view(dst);
    for (size_t i = 0; i < src.size(); i++) h(i) = src[i];
    Kokkos::deep_copy(dst, h);
  };

  auto copy_1d_double = [](const std::vector<double> &src, auto &dst, const char *label) {
    MemKK::realloc_kokkos(dst, label, src.size());
    auto h = Kokkos::create_mirror_view(dst);
    for (size_t i = 0; i < src.size(); i++) h(i) = src[i];
    Kokkos::deep_copy(dst, h);
  };

  // Chemical embedding [n_elements, embedding_size]
  copy_2d(m.chem_embedding, d_chem_embed, m.n_elements, m.embedding_size, "g1l:chem_embed");

  // A indicator transform [embedding_size, nradmax]
  copy_2d(m.A_lin_transform_W, d_A_lin_W, m.embedding_size, m.nradmax, "g1l:A_lin_W");
  A_lin_norm = (NNScalar) m.A_lin_transform_norm;
  inv_avg_n_neigh = (GeomScalar) m.A_inv_avg_n_neigh;

  // Precompute z_tr[mu, n] = (A_lin_norm * inv_avg_n_neigh) * Σ_e A_lin_W[e, n] * chem_embed[mu, e]
  {
    MemKK::realloc_kokkos(d_z_tr, "g1l:z_tr", m.n_elements, m.nradmax);
    auto h_z_tr = Kokkos::create_mirror_view(d_z_tr);
    NNScalar z_norm = A_lin_norm * inv_avg_n_neigh;
    for (int mu = 0; mu < m.n_elements; mu++)
      for (int n = 0; n < m.nradmax; n++) {
        NNScalar val = 0.0;
        for (int e = 0; e < m.embedding_size; e++)
          val += m.A_lin_transform_W[e * m.nradmax + n] * m.chem_embedding[mu * m.embedding_size + e];
        h_z_tr(mu, n) = val * z_norm;
      }
    Kokkos::deep_copy(d_z_tr, h_z_tr);
  }

  // MLP radial weights (generalized N-layer)
  {
    mlp_rad_n_layers = m.mlp_rad_n_layers;
    // Compute dims array: dims[0]=input of layer 0, dims[i+1]=output of layer i
    std::vector<int> dims(mlp_rad_n_layers + 1);
    dims[0] = m.mlp_rad_layers[0].n_in;
    for (int i = 0; i < mlp_rad_n_layers; i++)
      dims[i + 1] = m.mlp_rad_layers[i].n_out;
    // Compute max hidden dim for stack array sizing (exclude output layer dim)
    mlp_rad_max_dim = 0;
    for (int i = 0; i < mlp_rad_n_layers; i++)
      if (dims[i] > mlp_rad_max_dim) mlp_rad_max_dim = dims[i];
    // Validate compile-time limits
    if (mlp_rad_n_layers > GRACE1L_CPU_MAX_MLP_LAYERS)
      error->all(FLERR, "GRACE-1L: mlp_rad_n_layers ({}) exceeds GRACE1L_CPU_MAX_MLP_LAYERS ({})",
                 mlp_rad_n_layers, GRACE1L_CPU_MAX_MLP_LAYERS);
    if (mlp_rad_max_dim > GRACE1L_CPU_MAX_MLP_DIM)
      error->all(FLERR, "GRACE-1L: mlp hidden dim ({}) exceeds GRACE1L_CPU_MAX_MLP_DIM ({})",
                 mlp_rad_max_dim, GRACE1L_CPU_MAX_MLP_DIM);
    if (nradmax > GRACE1L_CPU_MAX_NRADMAX)
      error->all(FLERR, "GRACE-1L: nradmax ({}) exceeds GRACE1L_CPU_MAX_NRADMAX ({})",
                 nradmax, GRACE1L_CPU_MAX_NRADMAX);
    // Find max_in and max_out for weight tensor padding
    int max_in = 0, max_out = 0;
    for (int i = 0; i < mlp_rad_n_layers; i++) {
      if (m.mlp_rad_layers[i].n_in > max_in) max_in = m.mlp_rad_layers[i].n_in;
      if (m.mlp_rad_layers[i].n_out > max_out) max_out = m.mlp_rad_layers[i].n_out;
    }
    // Allocate and fill padded 3D weight tensor [n_layers, max_in, max_out]
    MemKK::realloc_kokkos(d_mlp_rad_W, "g1l:mlp_rad_W", mlp_rad_n_layers, max_in, max_out);
    auto h_W = Kokkos::create_mirror_view(d_mlp_rad_W);
    Kokkos::deep_copy(h_W, 0.0);
    for (int layer = 0; layer < mlp_rad_n_layers; layer++) {
      auto &L = m.mlp_rad_layers[layer];
      for (int i = 0; i < L.n_in; i++)
        for (int j = 0; j < L.n_out; j++)
          h_W(layer, i, j) = L.W[i * L.n_out + j];
    }
    Kokkos::deep_copy(d_mlp_rad_W, h_W);
    // Norms
    MemKK::realloc_kokkos(d_mlp_rad_norms, "g1l:mlp_rad_norms", mlp_rad_n_layers);
    auto h_norms = Kokkos::create_mirror_view(d_mlp_rad_norms);
    for (int i = 0; i < mlp_rad_n_layers; i++)
      h_norms(i) = (NNScalar) m.mlp_rad_layers[i].norm;
    Kokkos::deep_copy(d_mlp_rad_norms, h_norms);
    // Dims
    MemKK::realloc_kokkos(d_mlp_rad_dims, "g1l:mlp_rad_dims", mlp_rad_n_layers + 1);
    auto h_dims = Kokkos::create_mirror_view(d_mlp_rad_dims);
    for (int i = 0; i <= mlp_rad_n_layers; i++)
      h_dims(i) = dims[i];
    Kokkos::deep_copy(d_mlp_rad_dims, h_dims);
  }

  // FC weights
  auto copy_fc = [&](const GRACE1LCPUModel::FCWeights &fc, const char *pfx,
                      auto &dwl, auto &dwr, auto &dwtl, auto &dwtr,
                      auto &dct, auto &dcf, auto &dnof,
                      NNScalar &nl, NNScalar &nr) {
    std::string s(pfx);
    copy_3d(fc.w_left, dwl, fc.n_out, fc.w_left.size() / (fc.n_out * fc.w_shape_left), fc.w_shape_left, (s+"_wl").c_str());
    copy_3d(fc.w_right, dwr, fc.n_out, fc.w_right.size() / (fc.n_out * fc.w_shape_right), fc.w_shape_right, (s+"_wr").c_str());
    copy_1d_int(fc.w_tile_left, dwtl, (s+"_wtl").c_str());
    copy_1d_int(fc.w_tile_right, dwtr, (s+"_wtr").c_str());
    copy_1d_int(fc.collect_to, dct, (s+"_ct").c_str());
    copy_1d_int(fc.collect_from, dcf, (s+"_cf").c_str());
    // norm_out_factor is [n_funcs, 1, 1] - flatten to 1d
    std::vector<double> nof_flat(fc.n_funcs_left);
    for (int i = 0; i < fc.n_funcs_left && i * 1 < (int)fc.norm_out_factor.size(); i++)
      nof_flat[i] = fc.norm_out_factor[i];  // stride is 1*1=1 in flattened [n,1,1]
    copy_1d_double(nof_flat, dnof, (s+"_nof").c_str());
    nl = (NNScalar) fc.norm_left;
    nr = (NNScalar) fc.norm_right;
  };

  copy_fc(m.fc_A1, "g1l:fc_A1", d_fc_A1_wl, d_fc_A1_wr, d_fc_A1_wtl, d_fc_A1_wtr,
          d_fc_A1_ct, d_fc_A1_cf, d_fc_A1_nof, fc_A1_nl, fc_A1_nr);
  copy_fc(m.fc_AA1, "g1l:fc_AA1", d_fc_AA1_wl, d_fc_AA1_wr, d_fc_AA1_wtl, d_fc_AA1_wtr,
          d_fc_AA1_ct, d_fc_AA1_cf, d_fc_AA1_nof, fc_AA1_nl, fc_AA1_nr);
  copy_fc(m.fc_AA2, "g1l:fc_AA2", d_fc_AA2_wl, d_fc_AA2_wr, d_fc_AA2_wtl, d_fc_AA2_wtr,
          d_fc_AA2_ct, d_fc_AA2_cf, d_fc_AA2_nof, fc_AA2_nl, fc_AA2_nr);

  // Build gather indices for fused FC kernels (forward: group by target, reverse: group by source)
  auto build_gather = [&](const GRACE1LCPUModel::FCWeights &fc, int n_tgt,
                          auto &rg_start, auto &rg_src, auto &rg_tile,
                          auto &rev_start, auto &rev_src, auto &rev_tile,
                          const char *pfx) {
    int ncol = (int)fc.collect_to.size();
    int n_src_max = 0;
    for (int i = 0; i < ncol; i++)
      if (fc.collect_from[i] + 1 > n_src_max) n_src_max = fc.collect_from[i] + 1;
    std::string s(pfx);

    // Forward gather: group scatter entries by target (ct)
    {
      std::vector<int> count(n_tgt, 0);
      for (int i = 0; i < ncol; i++) count[fc.collect_to[i]]++;
      std::vector<int> start(n_tgt + 1, 0);
      for (int i = 0; i < n_tgt; i++) start[i + 1] = start[i] + count[i];
      std::vector<int> src(ncol), tile(ncol);
      std::vector<int> pos(n_tgt, 0);
      for (int i = 0; i < ncol; i++) {
        int t = fc.collect_to[i];
        int idx = start[t] + pos[t]++;
        src[idx] = fc.collect_from[i];
        tile[idx] = fc.w_tile_right[i];
      }
      copy_1d_int(start, rg_start, (s + "_rg_start").c_str());
      copy_1d_int(src, rg_src, (s + "_rg_src").c_str());
      copy_1d_int(tile, rg_tile, (s + "_rg_tile").c_str());
    }
    // Reverse gather: group scatter entries by source (cf)
    {
      std::vector<int> count(n_src_max, 0);
      for (int i = 0; i < ncol; i++) count[fc.collect_from[i]]++;
      std::vector<int> start(n_src_max + 1, 0);
      for (int i = 0; i < n_src_max; i++) start[i + 1] = start[i] + count[i];
      std::vector<int> tgt(ncol), tile(ncol);
      std::vector<int> pos(n_src_max, 0);
      for (int i = 0; i < ncol; i++) {
        int s2 = fc.collect_from[i];
        int idx = start[s2] + pos[s2]++;
        tgt[idx] = fc.collect_to[i];
        tile[idx] = fc.w_tile_right[i];
      }
      copy_1d_int(start, rev_start, (s + "_rev_start").c_str());
      copy_1d_int(tgt, rev_src, (s + "_rev_tgt").c_str());
      copy_1d_int(tile, rev_tile, (s + "_rev_tile").c_str());
    }
  };
  build_gather(m.fc_A1, n_funcs_A1,
    d_fc_A1_rg_start, d_fc_A1_rg_src, d_fc_A1_rg_tile,
    d_fc_A1_rev_start, d_fc_A1_rev_src, d_fc_A1_rev_tile, "g1l:fc_A1");
  build_gather(m.fc_AA1, n_funcs_AA,
    d_fc_AA1_rg_start, d_fc_AA1_rg_src, d_fc_AA1_rg_tile,
    d_fc_AA1_rev_start, d_fc_AA1_rev_src, d_fc_AA1_rev_tile, "g1l:fc_AA1");
  build_gather(m.fc_AA2, n_funcs_AA,
    d_fc_AA2_rg_start, d_fc_AA2_rg_src, d_fc_AA2_rg_tile,
    d_fc_AA2_rev_start, d_fc_AA2_rev_src, d_fc_AA2_rev_tile, "g1l:fc_AA2");

  // CG product metadata
  auto copy_prod = [&](const GRACE1LCPUModel::CGProduct &prod, const char *pfx,
                        auto &dli, auto &dri, auto &dsi, auto &dcg) {
    std::string s(pfx);
    copy_1d_int(prod.left_ind, dli, (s+"_li").c_str());
    copy_1d_int(prod.right_ind, dri, (s+"_ri").c_str());
    copy_1d_int(prod.m_sum_ind, dsi, (s+"_si").c_str());
    copy_1d_double(prod.cg_coeff, dcg, (s+"_cg").c_str());
  };

  copy_prod(m.prod_AA, "g1l:prod_AA", d_prod_AA_li, d_prod_AA_ri, d_prod_AA_si, d_prod_AA_cg);
  copy_prod(m.prod_AAA, "g1l:prod_AAA", d_prod_AAA_li, d_prod_AAA_ri, d_prod_AAA_si, d_prod_AAA_cg);
  copy_prod(m.prod_AAAA, "g1l:prod_AAAA", d_prod_AAAA_li, d_prod_AAAA_ri, d_prod_AAAA_si, d_prod_AAAA_cg);

  // Build forward product gather indices (group CG terms by output function)
  auto build_prod_gather = [&](const GRACE1LCPUModel::CGProduct &prod,
                               auto &gs, auto &gli, auto &gri, auto &gcg, const char *pfx) {
    int nt = prod.n_cg_terms, no = prod.n_output_funcs;
    std::string s(pfx);
    std::vector<int> count(no, 0);
    for (int t = 0; t < nt; t++) count[prod.m_sum_ind[t]]++;
    std::vector<int> start(no + 1, 0);
    for (int i = 0; i < no; i++) start[i + 1] = start[i] + count[i];
    std::vector<int> li(nt), ri(nt); std::vector<double> cg(nt);
    std::vector<int> pos(no, 0);
    for (int t = 0; t < nt; t++) {
      int o = prod.m_sum_ind[t];
      int idx = start[o] + pos[o]++;
      li[idx] = prod.left_ind[t];
      ri[idx] = prod.right_ind[t];
      cg[idx] = prod.cg_coeff[t];
    }
    copy_1d_int(start, gs, (s + "_gs").c_str());
    copy_1d_int(li, gli, (s + "_gli").c_str());
    copy_1d_int(ri, gri, (s + "_gri").c_str());
    copy_1d_double(cg, gcg, (s + "_gcg").c_str());
  };
  build_prod_gather(m.prod_AA, d_prod_AA_gs, d_prod_AA_gli, d_prod_AA_gri, d_prod_AA_gcg, "g1l:prod_AA");
  build_prod_gather(m.prod_AAA, d_prod_AAA_gs, d_prod_AAA_gli, d_prod_AAA_gri, d_prod_AAA_gcg, "g1l:prod_AAA");
  build_prod_gather(m.prod_AAAA, d_prod_AAAA_gs, d_prod_AAAA_gli, d_prod_AAAA_gri, d_prod_AAAA_gcg, "g1l:prod_AAAA");

  // Build reverse self-product gather: for each adjoint func f, entries where li==f or ri==f
  auto build_rprod_gather = [&](const GRACE1LCPUModel::CGProduct &prod, int n_adj,
                                auto &gs, auto &goa, auto &gfwd, auto &gcg, const char *pfx) {
    int nt = prod.n_cg_terms;
    std::string s(pfx);
    // Count entries per adjoint function
    std::vector<int> count(n_adj, 0);
    for (int t = 0; t < nt; t++) { count[prod.left_ind[t]]++; count[prod.right_ind[t]]++; }
    std::vector<int> start(n_adj + 1, 0);
    for (int i = 0; i < n_adj; i++) start[i + 1] = start[i] + count[i];
    int total = start[n_adj];
    std::vector<int> oa(total), fwd(total); std::vector<double> cg(total);
    std::vector<int> pos(n_adj, 0);
    for (int t = 0; t < nt; t++) {
      // as left: adj[li] += cg * oa[si] * fwd[ri]
      int f = prod.left_ind[t];
      int idx = start[f] + pos[f]++;
      oa[idx] = prod.m_sum_ind[t]; fwd[idx] = prod.right_ind[t]; cg[idx] = prod.cg_coeff[t];
      // as right: adj[ri] += cg * oa[si] * fwd[li]
      f = prod.right_ind[t];
      idx = start[f] + pos[f]++;
      oa[idx] = prod.m_sum_ind[t]; fwd[idx] = prod.left_ind[t]; cg[idx] = prod.cg_coeff[t];
    }
    copy_1d_int(start, gs, (s + "_gs").c_str());
    copy_1d_int(oa, goa, (s + "_goa").c_str());
    copy_1d_int(fwd, gfwd, (s + "_gfwd").c_str());
    copy_1d_double(cg, gcg, (s + "_gcg").c_str());
  };
  build_rprod_gather(m.prod_AA, n_funcs_A1,
    d_rprod_AA_gs, d_rprod_AA_goa, d_rprod_AA_gfwd, d_rprod_AA_gcg, "g1l:rprod_AA");
  build_rprod_gather(m.prod_AAAA, n_funcs_AA,
    d_rprod_AAAA_gs, d_rprod_AAAA_goa, d_rprod_AAAA_gfwd, d_rprod_AAAA_gcg, "g1l:rprod_AAAA");

  // Build reverse gather for asymmetric AAA product (AA1⊗A1 → AAA): separate left/right CSRs
  {
    const auto &prod = m.prod_AAA;
    int nt = prod.n_cg_terms;
    // --- Left gather (adj indexed by left_ind space: [0, n_funcs_AA)) ---
    {
      std::vector<int> count(n_funcs_AA, 0);
      for (int t = 0; t < nt; t++) count[prod.left_ind[t]]++;
      std::vector<int> start(n_funcs_AA + 1, 0);
      for (int i = 0; i < n_funcs_AA; i++) start[i + 1] = start[i] + count[i];
      int total = start[n_funcs_AA];
      std::vector<int> oa(total), fwd(total); std::vector<double> cg(total);
      std::vector<int> pos(n_funcs_AA, 0);
      for (int t = 0; t < nt; t++) {
        int f = prod.left_ind[t], idx = start[f] + pos[f]++;
        oa[idx] = prod.m_sum_ind[t]; fwd[idx] = prod.right_ind[t]; cg[idx] = prod.cg_coeff[t];
      }
      copy_1d_int(start, d_rprod_AAA_left_gs,  "g1l:rprod_AAA_left_gs");
      copy_1d_int(oa,    d_rprod_AAA_left_goa, "g1l:rprod_AAA_left_goa");
      copy_1d_int(fwd,   d_rprod_AAA_left_gfwd,"g1l:rprod_AAA_left_gfwd");
      copy_1d_double(cg, d_rprod_AAA_left_gcg, "g1l:rprod_AAA_left_gcg");
    }
    // --- Right gather (adj indexed by right_ind space: [0, n_funcs_A1)) ---
    {
      std::vector<int> count(n_funcs_A1, 0);
      for (int t = 0; t < nt; t++) count[prod.right_ind[t]]++;
      std::vector<int> start(n_funcs_A1 + 1, 0);
      for (int i = 0; i < n_funcs_A1; i++) start[i + 1] = start[i] + count[i];
      int total = start[n_funcs_A1];
      std::vector<int> oa(total), fwd(total); std::vector<double> cg(total);
      std::vector<int> pos(n_funcs_A1, 0);
      for (int t = 0; t < nt; t++) {
        int f = prod.right_ind[t], idx = start[f] + pos[f]++;
        oa[idx] = prod.m_sum_ind[t]; fwd[idx] = prod.left_ind[t]; cg[idx] = prod.cg_coeff[t];
      }
      copy_1d_int(start, d_rprod_AAA_right_gs,  "g1l:rprod_AAA_right_gs");
      copy_1d_int(oa,    d_rprod_AAA_right_goa, "g1l:rprod_AAA_right_goa");
      copy_1d_int(fwd,   d_rprod_AAA_right_gfwd,"g1l:rprod_AAA_right_gfwd");
      copy_1d_double(cg, d_rprod_AAA_right_gcg, "g1l:rprod_AAA_right_gcg");
    }
  }

  // FunctionReduceN weights [n_elements, rho_n_out, n_in, w_shape]
  auto copy_reduce = [&](const GRACE1LCPUModel::ReduceWeights &rw, const char *pfx,
                          auto &dw, NNScalar &norm, auto &dci) {
    std::string s(pfx);
    copy_4d(rw.W, dw, m.n_elements, m.rho_n_out, rw.n_in, rw.w_shape, (s+"_W").c_str());
    norm = (NNScalar) rw.norm;
    copy_1d_int(rw.collect_ind, dci, (s+"_ci").c_str());
  };

  copy_reduce(m.reduce_A, "g1l:red_A", d_reduce_A_W, reduce_A_norm, d_reduce_A_ci);
  copy_reduce(m.reduce_AA, "g1l:red_AA", d_reduce_AA_W, reduce_AA_norm, d_reduce_AA_ci);
  copy_reduce(m.reduce_AAA, "g1l:red_AAA", d_reduce_AAA_W, reduce_AAA_norm, d_reduce_AAA_ci);
  copy_reduce(m.reduce_AAAA, "g1l:red_AAAA", d_reduce_AAAA_W, reduce_AAAA_norm, d_reduce_AAAA_ci);

  // Energy MLP (generalized N-layer)
  {
    energy_n_layers = (int) m.energy_mlp_layers.size();
    std::vector<int> edims(energy_n_layers + 1);
    edims[0] = m.energy_mlp_layers[0].n_in;
    for (int i = 0; i < energy_n_layers; i++)
      edims[i + 1] = m.energy_mlp_layers[i].n_out;
    energy_max_dim = 0;
    for (int i = 0; i <= energy_n_layers; i++)
      if (edims[i] > energy_max_dim) energy_max_dim = edims[i];
    if (energy_n_layers > GRACE1L_CPU_MAX_MLP_LAYERS)
      error->all(FLERR, "GRACE-1L: energy_n_layers ({}) exceeds GRACE1L_CPU_MAX_MLP_LAYERS ({})",
                 energy_n_layers, GRACE1L_CPU_MAX_MLP_LAYERS);
    if (energy_max_dim > GRACE1L_CPU_MAX_MLP_DIM)
      error->all(FLERR, "GRACE-1L: energy hidden dim ({}) exceeds GRACE1L_CPU_MAX_MLP_DIM ({})",
                 energy_max_dim, GRACE1L_CPU_MAX_MLP_DIM);
    if (m.energy_mlp_activation != 0 && m.energy_mlp_activation != 1)
      error->all(FLERR, "GRACE-1L/KK: energy_mlp_activation={} not supported "
                        "(only 0=silu, 1=tanh)", m.energy_mlp_activation);
    energy_activation = m.energy_mlp_activation;
    int emax_in = 0, emax_out = 0;
    for (int i = 0; i < energy_n_layers; i++) {
      if (m.energy_mlp_layers[i].n_in > emax_in) emax_in = m.energy_mlp_layers[i].n_in;
      if (m.energy_mlp_layers[i].n_out > emax_out) emax_out = m.energy_mlp_layers[i].n_out;
    }
    MemKK::realloc_kokkos(d_energy_W, "g1l:energy_W", energy_n_layers, emax_in, emax_out);
    auto h_eW = Kokkos::create_mirror_view(d_energy_W);
    Kokkos::deep_copy(h_eW, 0.0);
    for (int layer = 0; layer < energy_n_layers; layer++) {
      auto &L = m.energy_mlp_layers[layer];
      for (int i = 0; i < L.n_in; i++)
        for (int j = 0; j < L.n_out; j++)
          h_eW(layer, i, j) = L.W[i * L.n_out + j];
    }
    Kokkos::deep_copy(d_energy_W, h_eW);
    MemKK::realloc_kokkos(d_energy_norms, "g1l:energy_norms", energy_n_layers);
    auto h_enorms = Kokkos::create_mirror_view(d_energy_norms);
    for (int i = 0; i < energy_n_layers; i++)
      h_enorms(i) = (NNScalar) m.energy_mlp_layers[i].norm;
    Kokkos::deep_copy(d_energy_norms, h_enorms);
    MemKK::realloc_kokkos(d_energy_dims, "g1l:energy_dims", energy_n_layers + 1);
    auto h_edims = Kokkos::create_mirror_view(d_energy_dims);
    for (int i = 0; i <= energy_n_layers; i++)
      h_edims(i) = edims[i];
    Kokkos::deep_copy(d_energy_dims, h_edims);
  }

  // Shifts
  copy_1d_double(m.shift_values, d_shifts, "g1l:shifts");
  output_scale = static_cast<NNScalar>(m.output_scale);

  // Bond-specific cutoff map [nelements, nelements]
  copy_2d(m.bond_cutoff_map, d_bond_cutoff, nelements, nelements, "g1l:bond_cutoff");
}

// ======================================================================
// grow: allocate per-chunk intermediate arrays
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::grow(int natom, int maxneigh_in)
{
  if ((int)d_A.extent(0) < natom) {
    const int nlm = (lmax + 1) * (lmax + 1);  // 25
    // Layout matches Python: [batch, n_channels, n_funcs]

    MemKK::realloc_kokkos(d_ncount, "g1l:ncount", natom);
    MemKK::realloc_kokkos(d_mu_i, "g1l:mu_i", natom);
    MemKK::realloc_kokkos(d_e_atom, "g1l:e_atom", natom);
    MemKK::realloc_kokkos(d_A, "g1l:A", natom, nradmax, nlm);
    MemKK::realloc_kokkos(d_A1, "g1l:A1", natom, fc_n_out, n_funcs_A1);
    MemKK::realloc_kokkos(d_AA, "g1l:AA", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_AA1, "g1l:AA1", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_AAA, "g1l:AAA", natom, fc_n_out, n_funcs_AAA);
    MemKK::realloc_kokkos(d_AA2, "g1l:AA2", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_AAAA, "g1l:AAAA", natom, fc_n_out, n_funcs_AAAA);
    MemKK::realloc_kokkos(d_rho, "g1l:rho", natom, rho_n_out);

    // Adjoint views for backward pass
    MemKK::realloc_kokkos(d_rho_adj, "g1l:rho_adj", natom, rho_n_out);
    MemKK::realloc_kokkos(d_AAAA_adj, "g1l:AAAA_adj", natom, fc_n_out, n_funcs_AAAA);
    MemKK::realloc_kokkos(d_AA2_adj, "g1l:AA2_adj", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_AAA_adj, "g1l:AAA_adj", natom, fc_n_out, n_funcs_AAA);
    MemKK::realloc_kokkos(d_AA1_adj, "g1l:AA1_adj", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_AA_adj, "g1l:AA_adj", natom, fc_n_out, n_funcs_AA);
    MemKK::realloc_kokkos(d_A1_adj, "g1l:A1_adj", natom, fc_n_out, n_funcs_A1);
    MemKK::realloc_kokkos(d_A_adj, "g1l:A_adj", natom, nradmax, nlm);
  }

  if ((int)d_radial_basis.extent(0) < natom || (int)d_radial_basis.extent(1) < maxneigh_in) {
    MemKK::realloc_kokkos(d_nearest, "g1l:nearest", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_rnorms, "g1l:rnorms", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_rhats, "g1l:rhats", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_mu_j, "g1l:mu_j", natom, maxneigh_in);
    MemKK::realloc_kokkos(d_radial_basis, "g1l:radial_basis", natom, maxneigh_in, nradbase);
    MemKK::realloc_kokkos(d_dradial_basis, "g1l:dradial_basis", natom, maxneigh_in, nradbase);
    MemKK::realloc_kokkos(d_R_nl, "g1l:R_nl", natom, maxneigh_in, nradmax, lmax + 1);
    MemKK::realloc_kokkos(d_dh2, "g1l:dh2", natom, maxneigh_in, mlp_rad_max_dim);
    MemKK::realloc_kokkos(d_f_ij, "g1l:f_ij", natom, maxneigh_in);
  }
}

// ======================================================================
// compute: main chunked pipeline
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  bool do_energy_only = flag_compute_energy_only && !debug_no_energy_only_calc;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag, vflag, 0);

  // Reallocate per-atom arrays if necessary
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }
  if (cvflag_atom) {
    memoryKK->destroy_kokkos(k_cvatom, cvatom);
    memoryKK->create_kokkos(k_cvatom, cvatom, maxcvatom, "pair:cvatom");
    d_cvatom = k_cvatom.view<DeviceType>();
  }

  copymode = 1;
  if (!force->newton_pair)
    error->all(FLERR, "PairGRACE1LCPUKokkos requires 'newton on'");

  atomKK->sync(execution_space, X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  // Find max neighbors
  maxneigh = 0;
  Kokkos::parallel_reduce("grace_1l::find_maxneigh", inum,
      FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  chunk_size = MIN(chunksize, inum);
  chunk_offset = 0;
  grow(chunk_size, maxneigh);

  EV_FLOAT ev;

  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    // Zero per-chunk arrays
    Kokkos::deep_copy(d_A, 0.0);
    Kokkos::deep_copy(d_rho, 0.0);

    EV_FLOAT ev_tmp;

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

    // ---- Stage 2: ComputeRadialBasis ----
    {
      int ts = team_size;
      check_team_size_for<TagComputeRadialBasis>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>(
          ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      Kokkos::parallel_for("ComputeRadialBasis", policy, *this);
    }

    // ---- Stage 3: ComputeMLPRadial ----
    {
      int ts = team_size;
      check_team_size_for<TagComputeMLPRadial>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>(
          ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      Kokkos::parallel_for("ComputeMLPRadial", policy, *this);
    }

    // ---- Stage 4: ComputeAi (CPU sandbox: per-atom, non-atomic) ----
    {
      Kokkos::parallel_for("ComputeAi",
        Kokkos::RangePolicy<DeviceType, TagComputeAi>(0, chunk_size), *this);
    }

    // ---- Stage 5: FC A1 = FC(A, A) — fused Left+Right (gather, no atomics) ----
    {
      auto left_in = d_A;  auto right_in = d_A;
      auto out = d_A1;
      auto wl = d_fc_A1_wl;  auto wr = d_fc_A1_wr;
      auto wtl = d_fc_A1_wtl;
      auto rgs = d_fc_A1_rg_start; auto rg_src = d_fc_A1_rg_src; auto rg_tile = d_fc_A1_rg_tile;
      auto nof = d_fc_A1_nof;
      NNScalar nl = fc_A1_nl, nr = fc_A1_nr;
      int n_lm = n_funcs_A1, cs = chunk_size, no = fc_n_out;

      Kokkos::parallel_for("FC_A1",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, no, n_lm}),
        KOKKOS_LAMBDA(const int ii, const int k, const int lm) {
          int nin_l = left_in.extent(1);
          NNScalar sum_l = 0.0;
          int tile_l = wtl(lm);
          for (int n = 0; n < nin_l; n++)
            sum_l += wl(k, n, tile_l) * left_in(ii, n, lm);
          NNScalar result = sum_l * nl;
          int nin_r = right_in.extent(1);
          for (int g = rgs(lm); g < rgs(lm + 1); g++) {
            int src = rg_src(g), tile_r = rg_tile(g);
            NNScalar sum_r = 0.0;
            for (int n = 0; n < nin_r; n++)
              sum_r += wr(k, n, tile_r) * right_in(ii, n, src);
            result += sum_r * nr;
          }
          out(ii, k, lm) = result * nof(lm);
        });
    }

    // ---- Stage 6: Product AA = A1 ⊗ A1 (CG-coupled, gather, no atomics) ----
    {
      auto left = d_A1; auto right = d_A1; auto out = d_AA;
      auto gs = d_prod_AA_gs; auto gli = d_prod_AA_gli;
      auto gri = d_prod_AA_gri; auto gcg = d_prod_AA_gcg;
      int n_out = n_funcs_AA, n_ch = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("Product_AA",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_out}),
        KOKKOS_LAMBDA(const int ii, const int o) {
          for (int k = 0; k < n_ch; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(o); g < gs(o + 1); g++)
              sum += gcg(g) * left(ii, k, gli(g)) * right(ii, k, gri(g));
            out(ii, k, o) = sum;
          }
        });
    }

    // ---- Stage 7: FC AA1 = FC(AA, A) — fused Left+Right (gather, no atomics) ----
    {
      auto left_in = d_AA;  auto right_in = d_A;
      auto out = d_AA1;
      auto wl = d_fc_AA1_wl; auto wr = d_fc_AA1_wr;
      auto wtl = d_fc_AA1_wtl;
      auto rgs = d_fc_AA1_rg_start; auto rg_src = d_fc_AA1_rg_src; auto rg_tile = d_fc_AA1_rg_tile;
      auto nof = d_fc_AA1_nof;
      NNScalar nl = fc_AA1_nl, nr = fc_AA1_nr;
      int nlm = n_funcs_AA, cs = chunk_size, no = fc_n_out;

      Kokkos::parallel_for("FC_AA1",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, no, nlm}),
        KOKKOS_LAMBDA(const int ii, const int k, const int lm) {
          int nin_l = left_in.extent(1);
          NNScalar sum_l = 0.0;
          int tile_l = wtl(lm);
          for (int n = 0; n < nin_l; n++)
            sum_l += wl(k, n, tile_l) * left_in(ii, n, lm);
          NNScalar result = sum_l * nl;
          int nin_r = right_in.extent(1);
          for (int g = rgs(lm); g < rgs(lm + 1); g++) {
            int src = rg_src(g), tile_r = rg_tile(g);
            NNScalar sum_r = 0.0;
            for (int n = 0; n < nin_r; n++)
              sum_r += wr(k, n, tile_r) * right_in(ii, n, src);
            result += sum_r * nr;
          }
          out(ii, k, lm) = result * nof(lm);
        });
    }

    // ---- Stage 8: Product AAA = AA1 ⊗ A1 (gather, no atomics) ----
    {
      auto left = d_AA1; auto right = d_A1; auto out = d_AAA;
      auto gs = d_prod_AAA_gs; auto gli = d_prod_AAA_gli;
      auto gri = d_prod_AAA_gri; auto gcg = d_prod_AAA_gcg;
      int n_out = n_funcs_AAA, n_ch = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("Product_AAA",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_out}),
        KOKKOS_LAMBDA(const int ii, const int o) {
          for (int k = 0; k < n_ch; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(o); g < gs(o + 1); g++)
              sum += gcg(g) * left(ii, k, gli(g)) * right(ii, k, gri(g));
            out(ii, k, o) = sum;
          }
        });
    }

    // ---- Stage 9: FC AA2 = FC(AA, A) — fused Left+Right (gather, no atomics) ----
    {
      auto left_in = d_AA;  auto right_in = d_A;
      auto out = d_AA2;
      auto wl = d_fc_AA2_wl; auto wr = d_fc_AA2_wr;
      auto wtl = d_fc_AA2_wtl;
      auto rgs = d_fc_AA2_rg_start; auto rg_src = d_fc_AA2_rg_src; auto rg_tile = d_fc_AA2_rg_tile;
      auto nof = d_fc_AA2_nof;
      NNScalar nl = fc_AA2_nl, nr = fc_AA2_nr;
      int n_lm = n_funcs_AA, cs = chunk_size, no = fc_n_out;

      Kokkos::parallel_for("FC_AA2",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, no, n_lm}),
        KOKKOS_LAMBDA(const int ii, const int k, const int lm) {
          int nin_l = left_in.extent(1);
          NNScalar sum_l = 0.0;
          int tile_l = wtl(lm);
          for (int n = 0; n < nin_l; n++)
            sum_l += wl(k, n, tile_l) * left_in(ii, n, lm);
          NNScalar result = sum_l * nl;
          int nin_r = right_in.extent(1);
          for (int g = rgs(lm); g < rgs(lm + 1); g++) {
            int src = rg_src(g), tile_r = rg_tile(g);
            NNScalar sum_r = 0.0;
            for (int n = 0; n < nin_r; n++)
              sum_r += wr(k, n, tile_r) * right_in(ii, n, src);
            result += sum_r * nr;
          }
          out(ii, k, lm) = result * nof(lm);
        });
    }

    // ---- Stage 10: Product AAAA = AA2 ⊗ AA2 (gather, no atomics) ----
    {
      auto left = d_AA2; auto right = d_AA2; auto out = d_AAAA;
      auto gs = d_prod_AAAA_gs; auto gli = d_prod_AAAA_gli;
      auto gri = d_prod_AAAA_gri; auto gcg = d_prod_AAAA_gcg;
      int n_out = n_funcs_AAAA, n_ch = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("Product_AAAA",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_out}),
        KOKKOS_LAMBDA(const int ii, const int o) {
          for (int k = 0; k < n_ch; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(o); g < gs(o + 1); g++)
              sum += gcg(g) * left(ii, k, gli(g)) * right(ii, k, gri(g));
            out(ii, k, o) = sum;
          }
        });
    }

    // ---- Stage 11: ComputeReduceN ----
    // Parallelize over (atom, rho_channel) for better GPU utilization
    {
      auto rho_out = d_rho;
      auto il = d_ilist; auto mp = d_map; auto tp = type;
      int co = chunk_offset, cs = chunk_size, nrho = rho_n_out;

      // Reduce from A: W[89,17,32,1] * A[ii,32,25]
      {
        auto W = d_reduce_A_W; auto ci = d_reduce_A_ci; auto X = d_A;
        NNScalar norm = reduce_A_norm;
        Kokkos::parallel_for("ReduceN_A",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, nrho}),
          KOKKOS_LAMBDA(const int ii, const int k) {
            const int mu_i = mp(tp(il[ii + co]));
            const int nin = W.extent(2), ws = W.extent(3);
            NNScalar val = 0.0;
            for (int w = 0; w < ws; w++) {
              const int fi = ci(w);
              for (int n = 0; n < nin; n++)
                val += W(mu_i, k, n, w) * X(ii, n, fi);
            }
            rho_out(ii, k) += val * norm;
          });
      }
      // Reduce from AA: W[89,17,64,5] * AA[ii,64,141]
      {
        auto W = d_reduce_AA_W; auto ci = d_reduce_AA_ci; auto X = d_AA;
        NNScalar norm = reduce_AA_norm;
        Kokkos::parallel_for("ReduceN_AA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, nrho}),
          KOKKOS_LAMBDA(const int ii, const int k) {
            const int mu_i = mp(tp(il[ii + co]));
            const int nin = W.extent(2), ws = W.extent(3);
            NNScalar val = 0.0;
            for (int w = 0; w < ws; w++) {
              const int fi = ci(w);
              for (int n = 0; n < nin; n++)
                val += W(mu_i, k, n, w) * X(ii, n, fi);
            }
            rho_out(ii, k) += val * norm;
          });
      }
      // Reduce from AAA: W[89,17,64,14] * AAA[ii,64,14]
      {
        auto W = d_reduce_AAA_W; auto ci = d_reduce_AAA_ci; auto X = d_AAA;
        NNScalar norm = reduce_AAA_norm;
        Kokkos::parallel_for("ReduceN_AAA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, nrho}),
          KOKKOS_LAMBDA(const int ii, const int k) {
            const int mu_i = mp(tp(il[ii + co]));
            const int nin = W.extent(2), ws = W.extent(3);
            NNScalar val = 0.0;
            for (int w = 0; w < ws; w++) {
              const int fi = ci(w);
              for (int n = 0; n < nin; n++)
                val += W(mu_i, k, n, w) * X(ii, n, fi);
            }
            rho_out(ii, k) += val * norm;
          });
      }
      // Reduce from AAAA: W[89,17,64,69] * AAAA[ii,64,69]
      {
        auto W = d_reduce_AAAA_W; auto ci = d_reduce_AAAA_ci; auto X = d_AAAA;
        NNScalar norm = reduce_AAAA_norm;
        Kokkos::parallel_for("ReduceN_AAAA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, nrho}),
          KOKKOS_LAMBDA(const int ii, const int k) {
            const int mu_i = mp(tp(il[ii + co]));
            const int nin = W.extent(2), ws = W.extent(3);
            NNScalar val = 0.0;
            for (int w = 0; w < ws; w++) {
              const int fi = ci(w);
              for (int n = 0; n < nin; n++)
                val += W(mu_i, k, n, w) * X(ii, n, fi);
            }
            rho_out(ii, k) += val * norm;
          });
      }
    }

    // ---- Stage 12: ComputeMLPEnergy ----
    {
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>(chunk_size, 1, 1);
      Kokkos::parallel_for("ComputeMLPEnergy", policy, *this);
    }

    // Accumulate energy
    if (eflag_global) {
      auto atom_energies = d_e_atom;
      GeomScalar energy_partial = 0.0;
      Kokkos::parallel_reduce("ComputeEnergySum", chunk_size,
          KOKKOS_LAMBDA(const int ii, GeomScalar& update) {
          update += atom_energies(ii);
      }, energy_partial);
      eng_vdwl += energy_partial;
    }

    // Write per-atom energy to LAMMPS array
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

    // ============ BACKWARD PASS (Forces) ============
    if (!do_energy_only) {

    // Zero adjoint arrays that are only partially written by RevReduceN.
    // d_AAA_adj and d_AAAA_adj are fully written (ws == n_funcs, ci injective),
    // so their deep_copy is not needed.
    Kokkos::deep_copy(d_A_adj, 0.0);
    Kokkos::deep_copy(d_AA_adj, 0.0);

    // ---- B1: ReverseMLP_Energy → d_rho_adj ----
    {
      auto rho_in = d_rho;
      auto rho_adj = d_rho_adj;
      auto eW = d_energy_W;
      auto enorms = d_energy_norms;
      auto edims = d_energy_dims;
      int en_layers = energy_n_layers;
      int emax = energy_max_dim;
      int cs = chunk_size;
      int eact = energy_activation;
      const NNScalar os = output_scale;
      Kokkos::parallel_for("ReverseMLPEnergy", Kokkos::RangePolicy<DeviceType>(0, cs),
        KOKKOS_LAMBDA(const int ii) {
          // Stack arrays: pre_act for hidden layers, adj ping-pong
          NNScalar pre_act_buf[GRACE1L_CPU_MAX_MLP_LAYERS * GRACE1L_CPU_MAX_MLP_DIM];
          NNScalar adj_a[GRACE1L_CPU_MAX_MLP_DIM], adj_b[GRACE1L_CPU_MAX_MLP_DIM];
          NNScalar* pre_act = pre_act_buf;
          NNScalar* adj_cur = adj_a;
          NNScalar* adj_nxt = adj_b;

          rho_adj(ii, 0) = os;  // linear term, scaled by output_scale

          // --- Forward pass: recompute all hidden pre-activations ---
          // Input to energy MLP is rho[1:] (skip rho[0])
          // Layer 0: pre_act[0] = norm0 * W0 @ rho[1:]
          {
            const int nin = edims(0);
            const int nout = edims(1);
            for (int j = 0; j < nout; j++) {
              NNScalar sum = 0.0;
              for (int k = 0; k < nin; k++)
                sum += eW(0, k, j) * rho_in(ii, k + 1);
              pre_act[0 * emax + j] = sum * enorms(0);
            }
          }
          // Hidden layers 1..en_layers-2: pre_act[layer] = norm * W @ act(pre_act[layer-1])
          for (int layer = 1; layer < en_layers - 1; layer++) {
            const int nin = edims(layer);
            const int nout = edims(layer + 1);
            for (int j = 0; j < nout; j++) {
              NNScalar sum = 0.0;
              for (int k = 0; k < nin; k++) {
                NNScalar pa = pre_act[(layer - 1) * emax + k];
                NNScalar act;
                if (eact == 1) {
                  act = Kokkos::tanh(pa);
                } else {
                  NNScalar sig = 1.0 / (1.0 + Kokkos::exp(-pa));
                  act = pa * sig;
                }
                sum += eW(layer, k, j) * act;
              }
              pre_act[layer * emax + j] = sum * enorms(layer);
            }
          }

          // --- Backward pass ---
          // Last layer (en_layers-1): output layer, no activation
          // adj for last hidden = W_last[:,0] * norm_last (output is scalar, index 0)
          {
            const int last = en_layers - 1;
            const int nin = edims(last);
            const int nout = edims(last + 1);
            // adj_cur = adjoint at input of last layer = silu(pre_act[last-1])
            // d_output/d_silu_input[k] = norm_last * W_last[k,o] for each output o
            // Since output is scalar (nout=1 typically), sum over outputs.
            // Multiply seed by output_scale so the chain rule scales all forces.
            for (int k = 0; k < nin; k++) {
              NNScalar val = 0.0;
              for (int o = 0; o < nout; o++)
                val += eW(last, k, o) * enorms(last);
              adj_cur[k] = val * os;
            }
          }

          // Backward through hidden layers (en_layers-2 down to 0)
          for (int layer = en_layers - 2; layer >= 0; layer--) {
            const int nin = edims(layer);
            const int nout = edims(layer + 1);
            // adj_cur[j] is adjoint at output of layer = act(pre_act[layer])
            // Apply act': d_pre_act = adj_cur * act'(pre_act)
            for (int j = 0; j < nout; j++) {
              NNScalar pa = pre_act[layer * emax + j];
              if (eact == 1) {
                NNScalar t = Kokkos::tanh(pa);
                adj_cur[j] *= (1.0 - t * t);
              } else {
                NNScalar sig = 1.0 / (1.0 + Kokkos::exp(-pa));
                adj_cur[j] *= sig * (1.0 + pa * (1.0 - sig));
              }
            }
            // Propagate through weight matrix: adj_nxt[k] = norm * Σ_j adj_cur[j] * W[layer,k,j]
            for (int k = 0; k < nin; k++) {
              NNScalar val = 0.0;
              for (int j = 0; j < nout; j++)
                val += adj_cur[j] * eW(layer, k, j);
              adj_nxt[k] = val * enorms(layer);
            }
            // Swap for next iteration
            NNScalar* tmp = adj_cur; adj_cur = adj_nxt; adj_nxt = tmp;
          }
          // adj_cur now holds d_output/d_rho[k+1] for k=0..edims(0)-1
          const int nin0 = edims(0);
          for (int k = 0; k < nin0; k++)
            rho_adj(ii, k + 1) = adj_cur[k];
        });
    }

    // ---- B2: ReverseReduceN → d_A_adj, d_AA_adj, d_AAA_adj, d_AAAA_adj ----
    // Parallelize over (atom, w) for each reduction block
    {
      auto rho_adj_v = d_rho_adj;
      auto il = d_ilist; auto tp = type; auto mp = d_map;
      int co = chunk_offset, cs = chunk_size, n_rho = rho_n_out;

      // All collect_ind arrays are verified injective (unique fi per w), so
      // each (ii, n, fi) is written by exactly one thread — no atomic needed.
      // RevReduceN_AAA/AAAA: ci is a full bijection over all funcs → deep_copy not needed.
      // RevReduceN_A/AA: ci is injective but partial → deep_copy zeros the rest (done above).

      // Reverse A
      {
        auto W = d_reduce_A_W; auto ci = d_reduce_A_ci; auto X = d_A_adj;
        NNScalar norm = reduce_A_norm;
        int ws = W.extent(3);
        Kokkos::parallel_for("RevReduceN_A",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, ws}),
          KOKKOS_LAMBDA(const int ii, const int w) {
            const int mu_i = mp(tp(il[ii + co]));
            const int fi = ci(w), nin = W.extent(2);
            for (int n = 0; n < nin; n++) {
              NNScalar val = 0.0;
              for (int k = 0; k < n_rho; k++)
                val += rho_adj_v(ii, k) * W(mu_i, k, n, w);
              X(ii, n, fi) = val * norm;
            }
          });
      }
      // Reverse AA
      {
        auto W = d_reduce_AA_W; auto ci = d_reduce_AA_ci; auto X = d_AA_adj;
        NNScalar norm = reduce_AA_norm;
        int ws = W.extent(3);
        Kokkos::parallel_for("RevReduceN_AA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, ws}),
          KOKKOS_LAMBDA(const int ii, const int w) {
            const int mu_i = mp(tp(il[ii + co]));
            const int fi = ci(w), nin = W.extent(2);
            for (int n = 0; n < nin; n++) {
              NNScalar val = 0.0;
              for (int k = 0; k < n_rho; k++)
                val += rho_adj_v(ii, k) * W(mu_i, k, n, w);
              X(ii, n, fi) = val * norm;
            }
          });
      }
      // Reverse AAA (ci is full bijection → no prior deep_copy needed)
      {
        auto W = d_reduce_AAA_W; auto ci = d_reduce_AAA_ci; auto X = d_AAA_adj;
        NNScalar norm = reduce_AAA_norm;
        int ws = W.extent(3);
        Kokkos::parallel_for("RevReduceN_AAA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, ws}),
          KOKKOS_LAMBDA(const int ii, const int w) {
            const int mu_i = mp(tp(il[ii + co]));
            const int fi = ci(w), nin = W.extent(2);
            for (int n = 0; n < nin; n++) {
              NNScalar val = 0.0;
              for (int k = 0; k < n_rho; k++)
                val += rho_adj_v(ii, k) * W(mu_i, k, n, w);
              X(ii, n, fi) = val * norm;
            }
          });
      }
      // Reverse AAAA (ci is full bijection → no prior deep_copy needed)
      {
        auto W = d_reduce_AAAA_W; auto ci = d_reduce_AAAA_ci; auto X = d_AAAA_adj;
        NNScalar norm = reduce_AAAA_norm;
        int ws = W.extent(3);
        Kokkos::parallel_for("RevReduceN_AAAA",
          Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, ws}),
          KOKKOS_LAMBDA(const int ii, const int w) {
            const int mu_i = mp(tp(il[ii + co]));
            const int fi = ci(w), nin = W.extent(2);
            for (int n = 0; n < nin; n++) {
              NNScalar val = 0.0;
              for (int k = 0; k < n_rho; k++)
                val += rho_adj_v(ii, k) * W(mu_i, k, n, w);
              X(ii, n, fi) = val * norm;
            }
          });
      }
    }

    // ---- B3: ReverseProduct_AAAA (AA2 ⊗ AA2, gather, no atomics) → d_AA2_adj ----
    {
      auto fwd = d_AA2; auto oa = d_AAAA_adj; auto adj = d_AA2_adj;
      auto gs = d_rprod_AAAA_gs; auto goa = d_rprod_AAAA_goa;
      auto gfwd = d_rprod_AAAA_gfwd; auto gcg = d_rprod_AAAA_gcg;
      int n_adj = n_funcs_AA, nc = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("RevProd_AAAA",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_adj}),
        KOKKOS_LAMBDA(const int ii, const int f) {
          for (int k = 0; k < nc; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(f); g < gs(f + 1); g++)
              sum += gcg(g) * oa(ii, k, goa(g)) * fwd(ii, k, gfwd(g));
            adj(ii, k, f) = sum;
          }
        });
    }

    // ---- B4: ReverseFC_AA2 → d_AA_adj +=, d_A_adj += ----
    {
      auto oa = d_AA2_adj; auto la = d_AA_adj; auto ra = d_A_adj;
      auto wl = d_fc_AA2_wl; auto wr = d_fc_AA2_wr;
      auto wtl = d_fc_AA2_wtl; auto wtr = d_fc_AA2_wtr;
      auto ct = d_fc_AA2_ct; auto cf = d_fc_AA2_cf; auto nof = d_fc_AA2_nof;
      NNScalar nl = fc_AA2_nl, nr = fc_AA2_nr;
      int no = fc_n_out, nlm = n_funcs_AA, ncol = (int)d_fc_AA2_cf.extent(0), cs = chunk_size;

      // Left part: la[ii,n,lm] += nof[lm]*nl * Σ_k wl[k,n,tile] * oa[ii,k,lm]
      { int nil = (int)la.extent(1);
      Kokkos::parallel_for("RevFC_AA2_L",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nil, nlm}),
        KOKKOS_LAMBDA(const int ii, const int n, const int lm) {
          NNScalar nf = nof(lm) * nl;
          int tile = wtl(lm);
          NNScalar val = 0.0;
          for (int k = 0; k < no; k++) val += wl(k, n, tile) * oa(ii, k, lm);
          la(ii, n, lm) += val * nf;
        }); }
      // Right part: gather by source — ra[ii,n,src] += Σ_{g} nof[tgt]*nr * Σ_k wr[k,n,tile] * oa[ii,k,tgt]
      { int nir = (int)ra.extent(1);
        auto revs = d_fc_AA2_rev_start; auto rev_tgt = d_fc_AA2_rev_src; auto rev_tile = d_fc_AA2_rev_tile;
        int n_src = (int)revs.extent(0) - 1;
      Kokkos::parallel_for("RevFC_AA2_R",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nir, n_src}),
        KOKKOS_LAMBDA(const int ii, const int n, const int src) {
          NNScalar val = 0.0;
          for (int g = revs(src); g < revs(src + 1); g++) {
            int tl = rev_tgt(g), tile = rev_tile(g);
            NNScalar gval = 0.0;
            for (int k = 0; k < no; k++) gval += wr(k, n, tile) * oa(ii, k, tl);
            val += gval * nof(tl);
          }
          ra(ii, n, src) += val * nr;
        }); }
    }

    // ---- B5: ReverseProduct_AAA (AA1 ⊗ A1, gather, no atomics) → d_AA1_adj, d_A1_adj ----
    {
      auto fwd_r = d_A1; auto oa = d_AAA_adj; auto adj = d_AA1_adj;
      auto gs = d_rprod_AAA_left_gs; auto goa = d_rprod_AAA_left_goa;
      auto gfwd = d_rprod_AAA_left_gfwd; auto gcg = d_rprod_AAA_left_gcg;
      int n_adj = n_funcs_AA, nc = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("RevProd_AAA_L",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_adj}),
        KOKKOS_LAMBDA(const int ii, const int f) {
          for (int k = 0; k < nc; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(f); g < gs(f + 1); g++)
              sum += gcg(g) * oa(ii, k, goa(g)) * fwd_r(ii, k, gfwd(g));
            adj(ii, k, f) = sum;
          }
        });
    }
    {
      auto fwd_l = d_AA1; auto oa = d_AAA_adj; auto adj = d_A1_adj;
      auto gs = d_rprod_AAA_right_gs; auto goa = d_rprod_AAA_right_goa;
      auto gfwd = d_rprod_AAA_right_gfwd; auto gcg = d_rprod_AAA_right_gcg;
      int n_adj = n_funcs_A1, nc = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("RevProd_AAA_R",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_adj}),
        KOKKOS_LAMBDA(const int ii, const int f) {
          for (int k = 0; k < nc; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(f); g < gs(f + 1); g++)
              sum += gcg(g) * oa(ii, k, goa(g)) * fwd_l(ii, k, gfwd(g));
            adj(ii, k, f) = sum;
          }
        });
    }

    // ---- B6: ReverseFC_AA1 → d_AA_adj +=, d_A_adj += ----
    {
      auto oa = d_AA1_adj; auto la = d_AA_adj; auto ra = d_A_adj;
      auto wl = d_fc_AA1_wl; auto wr = d_fc_AA1_wr;
      auto wtl = d_fc_AA1_wtl; auto wtr = d_fc_AA1_wtr;
      auto ct = d_fc_AA1_ct; auto cf = d_fc_AA1_cf; auto nof = d_fc_AA1_nof;
      NNScalar nl = fc_AA1_nl, nr = fc_AA1_nr;
      int no = fc_n_out, nlm = n_funcs_AA, ncol = (int)d_fc_AA1_cf.extent(0), cs = chunk_size;

      // Left part
      { int nil = (int)la.extent(1);
      Kokkos::parallel_for("RevFC_AA1_L",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nil, nlm}),
        KOKKOS_LAMBDA(const int ii, const int n, const int lm) {
          NNScalar nf = nof(lm) * nl;
          int tile = wtl(lm);
          NNScalar val = 0.0;
          for (int k = 0; k < no; k++) val += wl(k, n, tile) * oa(ii, k, lm);
          la(ii, n, lm) += val * nf;
        }); }
      // Right part: gather by source
      { int nir = (int)ra.extent(1);
        auto revs = d_fc_AA1_rev_start; auto rev_tgt = d_fc_AA1_rev_src; auto rev_tile = d_fc_AA1_rev_tile;
        int n_src = (int)revs.extent(0) - 1;
      Kokkos::parallel_for("RevFC_AA1_R",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, nir, n_src}),
        KOKKOS_LAMBDA(const int ii, const int n, const int src) {
          NNScalar val = 0.0;
          for (int g = revs(src); g < revs(src + 1); g++) {
            int tl = rev_tgt(g), tile = rev_tile(g);
            NNScalar gval = 0.0;
            for (int k = 0; k < no; k++) gval += wr(k, n, tile) * oa(ii, k, tl);
            val += gval * nof(tl);
          }
          ra(ii, n, src) += val * nr;
        }); }
    }

    // ---- B7: ReverseProduct_AA (A1 ⊗ A1, gather, no atomics) → d_A1_adj += ----
    {
      auto fwd = d_A1; auto oa = d_AA_adj; auto adj = d_A1_adj;
      auto gs = d_rprod_AA_gs; auto goa_v = d_rprod_AA_goa;
      auto gfwd_v = d_rprod_AA_gfwd; auto gcg = d_rprod_AA_gcg;
      int n_adj = n_funcs_A1, nc = fc_n_out, cs = chunk_size;

      Kokkos::parallel_for("RevProd_AA",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>({0,0}, {cs, n_adj}),
        KOKKOS_LAMBDA(const int ii, const int f) {
          for (int k = 0; k < nc; k++) {
            NNScalar sum = 0.0;
            for (int g = gs(f); g < gs(f + 1); g++)
              sum += gcg(g) * oa(ii, k, goa_v(g)) * fwd(ii, k, gfwd_v(g));
            adj(ii, k, f) += sum;
          }
        });
    }

    // ---- B8: ReverseFC_A1 → d_A_adj += (both left and right inputs are A) ----
    {
      auto oa = d_A1_adj; auto aa = d_A_adj;
      auto wl = d_fc_A1_wl; auto wr = d_fc_A1_wr;
      auto wtl = d_fc_A1_wtl; auto wtr = d_fc_A1_wtr;
      auto ct = d_fc_A1_ct; auto cf = d_fc_A1_cf; auto nof = d_fc_A1_nof;
      NNScalar nl = fc_A1_nl, nr = fc_A1_nr;
      int no = fc_n_out, nlm = n_funcs_A1, ncol = (int)d_fc_A1_cf.extent(0), cs = chunk_size;

      // Left part
      { int ni = (int)aa.extent(1);
      Kokkos::parallel_for("RevFC_A1_L",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, ni, nlm}),
        KOKKOS_LAMBDA(const int ii, const int n, const int lm) {
          NNScalar nf = nof(lm) * nl;
          int tile = wtl(lm);
          NNScalar val = 0.0;
          for (int k = 0; k < no; k++) val += wl(k, n, tile) * oa(ii, k, lm);
          aa(ii, n, lm) += val * nf;
        }); }
      // Right part: gather by source
      { int ni = (int)aa.extent(1);
        auto revs = d_fc_A1_rev_start; auto rev_tgt = d_fc_A1_rev_src; auto rev_tile = d_fc_A1_rev_tile;
        int n_src = (int)revs.extent(0) - 1;
      Kokkos::parallel_for("RevFC_A1_R",
        Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<3>>({0,0,0}, {cs, ni, n_src}),
        KOKKOS_LAMBDA(const int ii, const int n, const int src) {
          NNScalar val = 0.0;
          for (int g = revs(src); g < revs(src + 1); g++) {
            int tl = rev_tgt(g), tile = rev_tile(g);
            NNScalar gval = 0.0;
            for (int k = 0; k < no; k++) gval += wr(k, n, tile) * oa(ii, k, tl);
            val += gval * nof(tl);
          }
          aa(ii, n, src) += val * nr;
        }); }
    }

    // ---- B9: ComputeDerivative (per-bond forces from d_A_adj) ----
    {
      int ts = team_size;
      check_team_size_for<TagComputeDerivative>(((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      int plm_size = (lmax + 1) * (lmax + 2) / 2;
      // Scratch: plm + dplm (z_tr removed — precomputed as d_z_tr)
      int scratch_size = scratch_size_helper<GeomScalar>(2 * plm_size);
      auto policy = Kokkos::TeamPolicy<DeviceType, TagComputeDerivative>(
          ((chunk_size+ts-1)/ts)*maxneigh, ts, vector_length);
      policy = policy.set_scratch_size(0, Kokkos::PerThread(scratch_size));
      Kokkos::parallel_for("ComputeDerivative", policy, *this);
    }

    // ---- B10: ComputeForce (accumulate f_ij to atom forces + virial) ----
    {
      auto fij = d_f_ij; auto nc = d_ncount; auto near = d_nearest;
      auto fout = f; auto il = d_ilist;
      auto rh = d_rhats; auto rn = d_rnorms;
      int co = chunk_offset, cs = chunk_size;
      bool do_v = (vflag_global != 0);
      bool do_va = (vflag_atom != 0);
      bool do_cva = (cvflag_atom != 0);
      auto va = d_vatom;
      auto cva = d_cvatom;

      EV_FLOAT ev_force;
      Kokkos::parallel_reduce("ComputeForce", Kokkos::RangePolicy<DeviceType>(0, cs),
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
              const GeomScalar v0 = delx * fx, v1 = dely * fy, v2 = delz * fz;
              const GeomScalar v3 = delx * fy, v4 = delx * fz, v5 = dely * fz;
              if (do_v) {
                ev.v[0] += v0; ev.v[1] += v1; ev.v[2] += v2;
                ev.v[3] += v3; ev.v[4] += v4; ev.v[5] += v5;
              }
              if (do_va) {
                // Half to each atom (half neighbor list)
                const GeomScalar half = GeomScalar(0.5);
                Kokkos::atomic_add(&va(i, 0), (KK_FLOAT)(half * v0));
                Kokkos::atomic_add(&va(i, 1), (KK_FLOAT)(half * v1));
                Kokkos::atomic_add(&va(i, 2), (KK_FLOAT)(half * v2));
                Kokkos::atomic_add(&va(i, 3), (KK_FLOAT)(half * v3));
                Kokkos::atomic_add(&va(i, 4), (KK_FLOAT)(half * v4));
                Kokkos::atomic_add(&va(i, 5), (KK_FLOAT)(half * v5));
                Kokkos::atomic_add(&va(j, 0), (KK_FLOAT)(half * v0));
                Kokkos::atomic_add(&va(j, 1), (KK_FLOAT)(half * v1));
                Kokkos::atomic_add(&va(j, 2), (KK_FLOAT)(half * v2));
                Kokkos::atomic_add(&va(j, 3), (KK_FLOAT)(half * v3));
                Kokkos::atomic_add(&va(j, 4), (KK_FLOAT)(half * v4));
                Kokkos::atomic_add(&va(j, 5), (KK_FLOAT)(half * v5));
              }
              if (do_cva) {
                const GeomScalar half = GeomScalar(0.5);
                const GeomScalar v6 = dely * fx, v7 = delz * fx, v8 = delz * fy;
                Kokkos::atomic_add(&cva(i, 0), (KK_FLOAT)(half * v0));
                Kokkos::atomic_add(&cva(i, 1), (KK_FLOAT)(half * v1));
                Kokkos::atomic_add(&cva(i, 2), (KK_FLOAT)(half * v2));
                Kokkos::atomic_add(&cva(i, 3), (KK_FLOAT)(half * v3));
                Kokkos::atomic_add(&cva(i, 4), (KK_FLOAT)(half * v4));
                Kokkos::atomic_add(&cva(i, 5), (KK_FLOAT)(half * v5));
                Kokkos::atomic_add(&cva(i, 6), (KK_FLOAT)(half * v6));
                Kokkos::atomic_add(&cva(i, 7), (KK_FLOAT)(half * v7));
                Kokkos::atomic_add(&cva(i, 8), (KK_FLOAT)(half * v8));
                Kokkos::atomic_add(&cva(j, 0), (KK_FLOAT)(half * v0));
                Kokkos::atomic_add(&cva(j, 1), (KK_FLOAT)(half * v1));
                Kokkos::atomic_add(&cva(j, 2), (KK_FLOAT)(half * v2));
                Kokkos::atomic_add(&cva(j, 3), (KK_FLOAT)(half * v3));
                Kokkos::atomic_add(&cva(j, 4), (KK_FLOAT)(half * v4));
                Kokkos::atomic_add(&cva(j, 5), (KK_FLOAT)(half * v5));
                Kokkos::atomic_add(&cva(j, 6), (KK_FLOAT)(half * v6));
                Kokkos::atomic_add(&cva(j, 7), (KK_FLOAT)(half * v7));
                Kokkos::atomic_add(&cva(j, 8), (KK_FLOAT)(half * v8));
              }
            }
          }
        }, ev_force);

      if (vflag_global) {
        virial[0] += ev_force.v[0]; virial[1] += ev_force.v[1]; virial[2] += ev_force.v[2];
        virial[3] += ev_force.v[3]; virial[4] += ev_force.v[4]; virial[5] += ev_force.v[5];
      }
    }

    } // end if (!do_energy_only)

    Kokkos::fence();
    chunk_offset += chunk_size;
  }

  if (vflag_fdotr && !do_energy_only) pair_virial_fdotr_compute(this);

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

  atomKK->modified(execution_space, F_MASK);
  copymode = 0;
}

// ======================================================================
// Kernel: ComputeNeigh (reused from GRACE-FS)
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeNeigh,
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
      const GeomScalar r = Kokkos::sqrt(rsq);
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
// Kernel: ComputeRadialBasis - Chebyshev polynomials + cutoff
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeRadialBasis,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeRadialBasis>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const GeomScalar r = d_rnorms(ii, jj);
  // Per-bond cutoff: look up from bond_cutoff_map[mu_i, mu_j]
  const GeomScalar rcut_ij = d_bond_cutoff(d_mu_i(ii), d_mu_j(ii, jj));
  const GeomScalar x_norm = r / rcut_ij;  // x in [0, 1]

  // Smooth polynomial envelope cutoff (p-order):
  // fcut(x) = 1 - (p+1)(p+2)/2 * x^p + p*(p+2) * x^(p+1) - p*(p+1)/2 * x^(p+2)
  const int p = radial_basis_p;
  GeomScalar xp = GeomScalar(1.0);
  for (int ip = 0; ip < p; ip++) xp *= x_norm;
  const GeomScalar xp1 = xp * x_norm;     // x^(p+1)
  const GeomScalar xp2 = xp1 * x_norm;    // x^(p+2)
  const GeomScalar pp1 = GeomScalar(p) * (p + 1);
  const GeomScalar pp2 = GeomScalar(p) * (p + 2);
  const GeomScalar p1p2 = GeomScalar(p + 1) * (p + 2);
  const GeomScalar fcut = GeomScalar(1.0) - GeomScalar(0.5) * p1p2 * xp + pp2 * xp1 - GeomScalar(0.5) * pp1 * xp2;

  // Derivative of fcut w.r.t. r:
  // dfcut/dr = (1/rcut_ij) * [-p*(p+1)*(p+2)/2 * x^(p-1)
  //            + p*(p+1)*(p+2) * x^p - p*(p+1)*(p+2)/2 * x^(p+1)]
  GeomScalar dfcut = GeomScalar(0.0);
  if (x_norm > GeomScalar(1e-14)) {
    const GeomScalar xp_m1 = xp / x_norm;  // x^(p-1)
    const GeomScalar coeff = GeomScalar(p) * (p + 1) * (p + 2) * GeomScalar(0.5);
    dfcut = coeff * (-xp_m1 + GeomScalar(2.0) * xp - xp1) / rcut_ij;
  }

  // Chebyshev variable: x_cheb = 2*(r/rcut_ij) - 1 (maps [0, rcut_ij] to [-1, 1])
  const GeomScalar x_cheb = GeomScalar(2.0) * x_norm - GeomScalar(1.0);
  const GeomScalar dx_cheb_dr = GeomScalar(2.0) / rcut_ij;

  // Chebyshev polynomials T_1..T_{nradbase} (skip T_0, use kind=1)
  // T_0 = 1, T_1 = x, T_k = 2*x*T_{k-1} - T_{k-2}
  GeomScalar T_prev = GeomScalar(1.0);    // T_0
  GeomScalar T_curr = x_cheb;             // T_1
  GeomScalar dT_prev = GeomScalar(0.0);   // dT_0/dx = 0
  GeomScalar dT_curr = GeomScalar(1.0);   // dT_1/dx = 1

  // Store T_1 as first basis function (k=0 in output)
  d_radial_basis(ii, jj, 0) = T_curr * fcut;
  d_dradial_basis(ii, jj, 0) = dT_curr * dx_cheb_dr * fcut + T_curr * dfcut;

  for (int k = 1; k < nradbase; k++) {
    // T_{k+1} = 2*x*T_k - T_{k-1}
    const GeomScalar T_next = GeomScalar(2.0) * x_cheb * T_curr - T_prev;
    const GeomScalar dT_next = GeomScalar(2.0) * (T_curr + x_cheb * dT_curr) - dT_prev;
    T_prev = T_curr;
    dT_prev = dT_curr;
    T_curr = T_next;
    dT_curr = dT_next;

    // g_k = T_{k+1} * fcut  (output index k stores T_{k+1})
    d_radial_basis(ii, jj, k) = T_curr * fcut;
    d_dradial_basis(ii, jj, k) = dT_curr * dx_cheb_dr * fcut + T_curr * dfcut;
  }
}

// ======================================================================
// Kernel: ComputeMLPRadial - MLP(g_k) -> R_nl
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeMLPRadial,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPRadial>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  // Cache layer dims/norms in registers (avoid repeated device view reads)
  const int n_mlp_layers = mlp_rad_n_layers;
  int mlp_dims[GRACE1L_CPU_MAX_MLP_LAYERS + 1];
  NNScalar mlp_norms[GRACE1L_CPU_MAX_MLP_LAYERS];
  for (int i = 0; i <= n_mlp_layers; i++) mlp_dims[i] = d_mlp_rad_dims(i);
  for (int i = 0; i < n_mlp_layers; i++) mlp_norms[i] = d_mlp_rad_norms(i);

  NNScalar buf_a[GRACE1L_CPU_MAX_MLP_DIM], buf_b[GRACE1L_CPU_MAX_MLP_DIM];
  NNScalar dbuf_a[GRACE1L_CPU_MAX_MLP_DIM], dbuf_b[GRACE1L_CPU_MAX_MLP_DIM];
  NNScalar* h_cur = buf_a;  NNScalar* h_nxt = buf_b;
  NNScalar* dh_cur = dbuf_a; NNScalar* dh_nxt = dbuf_b;

  const int n_input = mlp_dims[0];
  for (int k = 0; k < n_input; k++) {
    h_cur[k] = d_radial_basis(ii, jj, k);
    dh_cur[k] = d_dradial_basis(ii, jj, k);
  }

  // Hidden layers: 0 .. n_layers-2 (all with silu activation)
  for (int layer = 0; layer < n_mlp_layers - 1; layer++) {
    const int nin  = mlp_dims[layer];
    const int nout = mlp_dims[layer + 1];
    const NNScalar norm = mlp_norms[layer];
    for (int j = 0; j < nout; j++) {
      NNScalar s = 0.0, ds = 0.0;
      for (int k = 0; k < nin; k++) {
        NNScalar w = d_mlp_rad_W(layer, k, j);
        s  += w * h_cur[k];
        ds += w * dh_cur[k];
      }
      s *= norm; ds *= norm;
      NNScalar sig = 1.0 / (1.0 + Kokkos::exp(-s));
      h_nxt[j]  = s * sig;
      dh_nxt[j] = sig * (1.0 + s * (1.0 - sig)) * ds;
    }
    // Swap cur/nxt
    NNScalar* tmp;
    tmp = h_cur;  h_cur  = h_nxt;  h_nxt  = tmp;
    tmp = dh_cur; dh_cur = dh_nxt; dh_nxt = tmp;
  }

  // Store dh_cur (derivative at last hidden layer output) for ComputeDerivative
  const int n_last_hidden = mlp_dims[n_mlp_layers - 1];
  for (int j = 0; j < n_last_hidden; j++)
    d_dh2(ii, jj, j) = dh_cur[j];

  // Output layer (last layer, no activation): fused into R_nl computation
  const int last_layer = n_mlp_layers - 1;
  const NNScalar out_norm = mlp_norms[last_layer];
  for (int n = 0; n < nradmax; n++) {
    for (int l = 0; l <= lmax; l++) {
      const int j = n * (lmax + 1) + l;
      NNScalar sum = 0.0;
      for (int k = 0; k < n_last_hidden; k++)
        sum += d_mlp_rad_W(last_layer, k, j) * h_cur[k];
      d_R_nl(ii, jj, n, l) = sum * out_norm;
    }
  }
}

// ======================================================================
// Kernel: ComputeAi - A[i,lm,n] = sum_j R_nl * Y_lm * Z_tr
// Adapted from GRACE-FS ComputeAi with R_nl from MLP and Z_tr from embedding
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeAi,
    const int ii) const
{
  if (ii >= chunk_size) return;
  const int ncount = d_ncount(ii);

  const GeomScalar Y00_v = Y00;
  const GeomScalar sq3_v = sq3;
  const GeomScalar sq3o2_v = sq3o2;
  const GeomScalar sq2_v = sq2;

  GeomScalar plm[GRACE1L_CPU_MAX_PLM];

  for (int jj = 0; jj < ncount; jj++) {
    const int mu_j = d_mu_j(ii, jj);
    const GeomScalar rx = d_rhats(ii, jj, 0);
    const GeomScalar ry = d_rhats(ii, jj, 1);
    const GeomScalar rz = d_rhats(ii, jj, 2);

    plm[0] = Y00_v;
    if (lmax > 0) {
      plm[1] = Y00_v * sq3_v * rz;
      plm[2] = -sq3o2_v * Y00_v;
      for (int l = 2; l <= lmax; l++) {
        for (int m = 0; m < l - 1; m++) {
          const int idx = l * (l + 1) / 2 + m;
          const int idx_l1_m = (l - 1) * l / 2 + m;
          const int idx_l2_m = (l - 2) * (l - 1) / 2 + m;
          const int alm_idx = d_idx_sph(l * (l + 1) + m);
          plm[idx] = alm(alm_idx) * (rz * plm[idx_l1_m] + blm(alm_idx) * plm[idx_l2_m]);
        }
        {
          const int idx = l * (l + 1) / 2 + l - 1;
          plm[idx] = dl(l) * plm[(l - 1) * l / 2 + l - 1] * rz;
        }
        {
          const int idx = l * (l + 1) / 2 + l;
          plm[idx] = cl(l) * plm[(l - 1) * l / 2 + l - 1];
        }
      }
    }

    const GeomScalar phase_re = rx;
    const GeomScalar phase_im = ry;

    // m = 0 (non-atomic: this thread owns ii)
    for (int l = 0; l <= lmax; l++) {
      const GeomScalar Y = plm[l * (l + 1) / 2];
      const int A_idx = l * (l + 1);
      for (int n = 0; n < nradmax; n++) {
        d_A(ii, n, A_idx) +=
            (NNScalar)(d_R_nl(ii, jj, n, l) * Y) * d_z_tr(mu_j, n);
      }
    }

    GeomScalar phasem_re = phase_re;
    GeomScalar phasem_im = phase_im;

    for (int m = 1; m <= lmax; m++) {
      if (m >= 2) {
        const GeomScalar tmp_re = phasem_re * phase_re - phasem_im * phase_im;
        const GeomScalar tmp_im = phasem_re * phase_im + phasem_im * phase_re;
        phasem_re = tmp_re;
        phasem_im = tmp_im;
      }
      const int factor = (m % 2 == 0) ? 1 : -1;

      for (int l = m; l <= lmax; l++) {
        const GeomScalar plm_val = plm[l * (l + 1) / 2 + m];
        const GeomScalar ylm_re = phasem_re * plm_val;
        const GeomScalar ylm_im = phasem_im * plm_val;

        const GeomScalar real_Y_pos = sq2_v * factor * ylm_re;
        const GeomScalar real_Y_neg = sq2_v * factor * ylm_im;

        const int A_idx_pos = l * (l + 1) + m;
        const int A_idx_neg = l * (l + 1) - m;

        for (int n = 0; n < nradmax; n++) {
          const NNScalar R = d_R_nl(ii, jj, n, l);
          const NNScalar zn = d_z_tr(mu_j, n);
          d_A(ii, n, A_idx_pos) += R * (NNScalar)real_Y_pos * zn;
          d_A(ii, n, A_idx_neg) += R * (NNScalar)real_Y_neg * zn;
        }
      }
    }
  }
}

// ======================================================================
// Kernel: ComputeMLPEnergy - E = MLP(rho[1:]) + linear(rho[0]) + shift
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeMLPEnergy,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeMLPEnergy>::member_type& team) const
{
  const int ii = team.league_rank();
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  // LinMLPOut2ScalarTarget: input is rho (17 channels)
  // rho[0] is the "linear" channel (added directly)
  // rho[1:] go through MLP
  const NNScalar linear_term = d_rho(ii, 0);

  NNScalar ebuf_a[GRACE1L_CPU_MAX_MLP_DIM], ebuf_b[GRACE1L_CPU_MAX_MLP_DIM];
  NNScalar* h_cur = ebuf_a;
  NNScalar* h_nxt = ebuf_b;

  // Initialize h_cur = rho[1:] (skip rho[0])
  const int nin0 = d_energy_dims(0);
  for (int k = 0; k < nin0; k++)
    h_cur[k] = d_rho(ii, k + 1);

  // Hidden layers: 0 .. energy_n_layers-2 (silu or tanh activation)
  for (int layer = 0; layer < energy_n_layers - 1; layer++) {
    const int nin  = d_energy_dims(layer);
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

  // Output layer (last layer, no activation)
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

  // Total: E_atom = (MLP_out + linear_term) * output_scale + shift
  d_e_atom(ii) = (e + linear_term) * output_scale + d_shifts(mu_i);
}

// ======================================================================
// Kernel: ComputeDerivative - per-bond forces from d_A_adj
// Uses forward-mode AD through MLP for dR_nl/dr + SH derivatives
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeDerivative,
    const typename Kokkos::TeamPolicy<DeviceType, TagComputeDerivative>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const int mu_j = d_mu_j(ii, jj);
  const GeomScalar rinv = GeomScalar(1.0) / d_rnorms(ii, jj);
  const GeomScalar rx = d_rhats(ii, jj, 0);
  const GeomScalar ry = d_rhats(ii, jj, 1);
  const GeomScalar rz = d_rhats(ii, jj, 2);

  // Local copies of file-scope constants to avoid double promotion
  const GeomScalar Y00_v = Y00;
  const GeomScalar sq3_v = sq3;
  const GeomScalar sq3o2_v = sq3o2;
  const GeomScalar sq2_v = sq2;

  // ---- Scratch: plm + dplm (z_tr precomputed in d_z_tr) ----
  int plm_size = (lmax + 1) * (lmax + 2) / 2;
  GeomScalar* scratch_geom = (GeomScalar*) team.thread_scratch(0).get_shmem(
      2 * plm_size * sizeof(GeomScalar));
  GeomScalar* plm = scratch_geom;
  GeomScalar* dplm = scratch_geom + plm_size;

  // dh2 was precomputed in ComputeMLPRadial
  const int n_last_hidden = d_mlp_rad_dims(mlp_rad_n_layers - 1);
  const int last_layer = mlp_rad_n_layers - 1;
  const NNScalar norm_last = d_mlp_rad_norms(last_layer);

  plm[0] = Y00_v; dplm[0] = GeomScalar(0.0);
  if (lmax > 0) {
    plm[1] = Y00_v * sq3_v * rz;     dplm[1] = Y00_v * sq3_v;
    plm[2] = -sq3o2_v * Y00_v;       dplm[2] = GeomScalar(0.0);
    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        const int idx = l*(l+1)/2 + m;
        const int i1 = (l-1)*l/2 + m, i2 = (l-2)*(l-1)/2 + m;
        const int ai = d_idx_sph(l*(l+1) + m);
        const GeomScalar a = alm(ai), b = blm(ai);
        plm[idx] = a * (rz * plm[i1] + b * plm[i2]);
        dplm[idx] = a * (plm[i1] + rz * dplm[i1] + b * dplm[i2]);
      }
      { const int idx = l*(l+1)/2+l-1, prev = (l-1)*l/2+l-1;
        const GeomScalar t = dl(l) * plm[prev];
        plm[idx] = t * rz;  dplm[idx] = t; }
      { const int idx = l*(l+1)/2+l, prev = (l-1)*l/2+l-1;
        plm[idx] = cl(l) * plm[prev];  dplm[idx] = 0.0; }
    }
  }

  // ---- Cache dh2[k] (read 5*32=160 times per (ii,jj) in DR_n inner loop) ----
  // and z_tr[mu_j, n] (read ~25 times per (ii,jj,n) across the (l,m,n) loops).
  // Same hoist pattern as 2L Derivative_L1/L2 (commit ca2acd6567, -2.2% FP64).
  NNScalar dh2_cache[GRACE1L_CPU_MAX_MLP_DIM];
  for (int k = 0; k < n_last_hidden; k++) dh2_cache[k] = d_dh2(ii, jj, k);
  NNScalar zn_cache[GRACE1L_CPU_MAX_NRADMAX];
  for (int n = 0; n < nradmax; n++) zn_cache[n] = d_z_tr(mu_j, n);

  // ---- Accumulate per-bond force ----
  // Restructured loop: (l, m) outermost, n innermost
  // Angular derivatives depend only on (l, m, rhat) — computed once per (l,m)
  // DR_n[n] and R_n[n] precomputed per l value
  GeomScalar f_ji[3] = {GeomScalar(0.0), GeomScalar(0.0), GeomScalar(0.0)};
  const GeomScalar phase_re = rx, phase_im = ry;

  NNScalar DR_n[GRACE1L_CPU_MAX_NRADMAX];
  GeomScalar R_n[GRACE1L_CPU_MAX_NRADMAX];

  for (int l = 0; l <= lmax; l++) {
    // Precompute DR_n[n] and R_n[n] for this l
    for (int n = 0; n < nradmax; n++) {
      R_n[n] = d_R_nl(ii, jj, n, l);
      const int j = n * (lmax + 1) + l;
      NNScalar dr = 0.0;
      for (int k = 0; k < n_last_hidden; k++)
        dr += d_mlp_rad_W(last_layer, k, j) * dh2_cache[k];
      DR_n[n] = dr * norm_last;
    }

    // m = 0: angular derivatives computed once
    {
      const GeomScalar Y = plm[l*(l+1)/2];
      const GeomScalar dp = dplm[l*(l+1)/2];
      const GeomScalar rdy = dp * rz;
      const GeomScalar DY_x = -rdy * rx, DY_y = -rdy * ry, DY_z = dp - rdy * rz;
      const int A_idx = l * (l + 1);

      for (int n = 0; n < nradmax; n++) {
        const NNScalar w = d_A_adj(ii, n, A_idx) * zn_cache[n];
        const GeomScalar R_over_r = R_n[n] * rinv;
        const GeomScalar YDR = Y * DR_n[n];
        f_ji[0] += w * (YDR * rx + DY_x * R_over_r);
        f_ji[1] += w * (YDR * ry + DY_y * R_over_r);
        f_ji[2] += w * (YDR * rz + DY_z * R_over_r);
      }
    }

    // m >= 1: angular derivatives computed once per m, then n-inner loop
    GeomScalar pm_re = phase_re, pm_im = phase_im;
    for (int m = 1; m <= l; m++) {
      const int fac = (m % 2 == 0) ? 1 : -1;
      const GeomScalar pv = plm[l*(l+1)/2+m], dpv = dplm[l*(l+1)/2+m];
      const GeomScalar ylm_re = pm_re * pv, ylm_im = pm_im * pv;
      const GeomScalar rYp = sq2_v * fac * ylm_re, rYn = sq2_v * fac * ylm_im;

      // Angular derivatives (independent of n — computed once)
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
          const GeomScalar mp = (GeomScalar)m * pv;
          dyx_re = mp*p1r; dyx_im = mp*p1i;
          dyy_re = -dyx_im; dyy_im = dyx_re;
        } else {
          dyx_re = dyx_im = dyy_re = dyy_im = GeomScalar(0.0);
        }
      }

      const GeomScalar rdy_re = rx*dyx_re + ry*dyy_re + rz*dyz_re;
      const GeomScalar rdy_im = rx*dyx_im + ry*dyy_im + rz*dyz_im;
      const GeomScalar Dpx_re = dyx_re - rdy_re*rx, Dpx_im = dyx_im - rdy_im*rx;
      const GeomScalar Dpy_re = dyy_re - rdy_re*ry, Dpy_im = dyy_im - rdy_im*ry;
      const GeomScalar Dpz_re = dyz_re - rdy_re*rz, Dpz_im = dyz_im - rdy_im*rz;
      const GeomScalar DYpx = sq2_v*fac*Dpx_re, DYpy = sq2_v*fac*Dpy_re, DYpz = sq2_v*fac*Dpz_re;
      const GeomScalar DYnx = sq2_v*fac*Dpx_im, DYny = sq2_v*fac*Dpy_im, DYnz = sq2_v*fac*Dpz_im;

      const int Ap = l*(l+1)+m, An = l*(l+1)-m;

      for (int n = 0; n < nradmax; n++) {
        const NNScalar zn = zn_cache[n];
        const NNScalar wp = d_A_adj(ii, n, Ap) * zn;
        const NNScalar wn = d_A_adj(ii, n, An) * zn;
        const GeomScalar R_over_r = R_n[n] * rinv;
        {
          const GeomScalar YDR = rYp * DR_n[n];
          f_ji[0] += wp * (YDR*rx + DYpx*R_over_r);
          f_ji[1] += wp * (YDR*ry + DYpy*R_over_r);
          f_ji[2] += wp * (YDR*rz + DYpz*R_over_r);
        }
        {
          const GeomScalar YDR = rYn * DR_n[n];
          f_ji[0] += wn * (YDR*rx + DYnx*R_over_r);
          f_ji[1] += wn * (YDR*ry + DYny*R_over_r);
          f_ji[2] += wn * (YDR*rz + DYnz*R_over_r);
        }
      }

      const GeomScalar t_re = pm_re*phase_re - pm_im*phase_im;
      const GeomScalar t_im = pm_re*phase_im + pm_im*phase_re;
      pm_re = t_re; pm_im = t_im;
    }
  }

  d_f_ij(ii, jj, 0) = f_ji[0];
  d_f_ij(ii, jj, 1) = f_ji[1];
  d_f_ij(ii, jj, 2) = f_ji[2];
}

// ======================================================================
// Kernel: ComputeForce (placeholder - force done via lambda in compute())
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeForce<NEIGHFLAG,EVFLAG>, const int& /*ii*/) const
{
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::operator() (TagComputeForce<NEIGHFLAG,EVFLAG>, const int& /*ii*/, EV_FLOAT& /*ev*/) const
{
}

// ======================================================================
// Utility functions
// ======================================================================

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<class TagStyle>
void PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::check_team_size_for(int inum_in, int &team_size, int vector_length) {
  int team_size_max;
  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum_in,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
template<typename scratch_type>
KOKKOS_INLINE_FUNCTION
int PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::scratch_size_helper(int values_per_team) const {
  typedef Kokkos::View<scratch_type*, typename DeviceType::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ScratchViewType;
  return ScratchViewType::shmem_size(values_per_team);
}

template<class DeviceType, typename NNScalarT, typename GeomScalarT>
double PairGRACE1LCPUKokkos<DeviceType, NNScalarT, GeomScalarT>::memory_usage()
{
  double bytes = 0;
  bytes += MemKK::memory_usage(d_A);
  bytes += MemKK::memory_usage(d_R_nl);
  bytes += MemKK::memory_usage(d_radial_basis);
  bytes += MemKK::memory_usage(d_e_atom);
  bytes += MemKK::memory_usage(d_rho);
  return bytes;
}

// ======================================================================
// Template instantiations
// ======================================================================

namespace LAMMPS_NS {
// FP64 (all double)
template class PairGRACE1LCPUKokkos<LMPDeviceType, double>;
#ifdef LMP_KOKKOS_GPU
template class PairGRACE1LCPUKokkos<LMPHostType, double>;
#endif

// Mixed (NN=float, geometry=double)
template class PairGRACE1LCPUKokkos<LMPDeviceType, float>;
#ifdef LMP_KOKKOS_GPU
template class PairGRACE1LCPUKokkos<LMPHostType, float>;
#endif

// FP32 (all float)
template class PairGRACE1LCPUKokkos<LMPDeviceType, float, float>;
#ifdef LMP_KOKKOS_GPU
template class PairGRACE1LCPUKokkos<LMPHostType, float, float>;
#endif
}
