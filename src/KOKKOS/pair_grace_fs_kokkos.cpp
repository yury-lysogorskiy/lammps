// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Yury Lysogorskiy (ICAMS)
------------------------------------------------------------------------- */

#include "pair_grace_fs_kokkos.h"

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

#include "ace-evaluator/ace_version.h"
#include "ace-evaluator/ace_radial.h"
#include "ace/grace_fs_evaluator.h"

#include <cstring>

// ACEImpl for GRACE FS
namespace LAMMPS_NS {
  struct ACEImpl {
    GRACEFSBasisSet *basis_set;
    GRACEFSBEvaluator *ace;
  };
} // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

enum{FS,FS_SHIFTEDSCALED};

// Constants for spherical harmonics
static constexpr double Y00_kk = 1.0; // 0.28209479177387814347403972578;
static constexpr double sq3_kk = 1.7320508075688772935;
static constexpr double sq3o2_kk = 1.2247448713915890491;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGRACEFSKokkos<DeviceType>::PairGRACEFSKokkos(LAMMPS *lmp) : PairGRACEFS(lmp)
{
  respa_enable = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  host_flag = (execution_space == Host);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGRACEFSKokkos<DeviceType>::~PairGRACEFSKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
  deallocate_views_of_views();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::deallocate_views_of_views()
{
  if (k_splines_gk.view_host().data()) {
    for (int i = 0; i < nelements; i++) {
      k_splines_gk.view_host()(i).deallocate();
      k_splines_rnl.view_host()(i).deallocate();
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::allocate()
{
  PairGRACEFS::allocate();
  int n = atom->ntypes + 1;
  MemKK::realloc_kokkos(d_map, "grace_fs:map", n);
  MemKK::realloc_kokkos(k_cutsq, "grace_fs:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();
  MemKK::realloc_kokkos(k_scale, "grace_fs:scale", n, n);
  d_scale = k_scale.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairGRACEFS::coeff(narg,arg);

  auto h_map = Kokkos::create_mirror_view(d_map);
  for (int i = 1; i <= atom->ntypes; i++)
    h_map(i) = map[i];
  Kokkos::deep_copy(d_map,h_map);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::init_style()
{
  if (host_flag) {
    PairGRACEFS::init_style();
    return;
  }

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace/fs/kk requires atom IDs");
  if (force->newton_pair == 0) {
    if (comm->me == 0) error->warning(FLERR, "Pair style grace/fs/kk requires newton pair on. Auto-enabling 'newton on'.");
    force->newton_pair = 1;
  }

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair grace/fs/kk");

  if (request_extrapolation) {
      flag_compute_extrapolation_grade = 0;
      if (extrapolation_grade_gamma == nullptr) {
          int nmax = atom->nmax;
          memory->create(extrapolation_grade_gamma, nmax, "grace_fs/atom:gamma");
          memset(extrapolation_grade_gamma, 0, nmax * sizeof(*extrapolation_grade_gamma));
      }
  }

  auto basis_set = aceimpl->basis_set;
  nelements = basis_set->nelements;
  lmax = basis_set->lmax;
  nradmax = basis_set->nradmax;
  nradbase = basis_set->nradbase;
  nnorm = basis_set->nnorm;
  energy_scale = basis_set->scale;
  energy_shift = basis_set->shift;

  // spherical harmonics
  MemKK::realloc_kokkos(d_idx_sph, "grace_fs:idx_sph", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(alm, "grace_fs:alm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(blm, "grace_fs:blm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(cl, "grace_fs:cl", lmax + 1);
  MemKK::realloc_kokkos(dl, "grace_fs:dl", lmax + 1);

  pre_compute_harmonics(lmax);
  copy_pertype();
  copy_splines();
  copy_radial_Z();
  copy_tilde();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairGRACEFSKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairGRACEFS::init_one(i,j);
  k_scale.view_host()(i,j) = k_scale.view_host()(j,i) = scale[i][j];
  k_scale.modify_host();
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::grow(int natom, int maxneigh)
{
  auto basis_set = aceimpl->basis_set;

  if ((int)A.extent(0) < natom) {
    MemKK::realloc_kokkos(A_sph, "grace_fs:A_sph", natom, (lmax + 1) * (lmax + 1), nradmax);
    MemKK::realloc_kokkos(A, "grace_fs:A", natom, (lmax + 1) * (lmax + 1), nradmax);
    MemKK::realloc_kokkos(A_list, "grace_fs:A_list", natom, idx_ms_combs_max * basis_set->rankmax);
    MemKK::realloc_kokkos(A_forward_prod, "grace_fs:A_forward_prod", natom, idx_ms_combs_max * (basis_set->rankmax + 1));
    MemKK::realloc_kokkos(e_atom, "grace_fs:e_atom", natom);
    MemKK::realloc_kokkos(rhos, "grace_fs:rhos", natom, basis_set->ndensitymax);
    MemKK::realloc_kokkos(dF_drho, "grace_fs:dF_drho", natom, basis_set->ndensitymax);
    MemKK::realloc_kokkos(weights, "grace_fs:weights", natom, (lmax + 1) * (lmax + 1), nradmax);
    MemKK::realloc_kokkos(dB_flatten, "grace_fs:dB_flatten", natom, idx_ms_combs_max * basis_set->rankmax);
    MemKK::realloc_kokkos(projections, "grace_fs:projections", natom, total_num_functions_max);
    MemKK::realloc_kokkos(d_gamma, "grace_fs:gamma", natom);
  }

  if (((int)fr.extent(0) < natom) || ((int)fr.extent(1) < maxneigh)) {
    MemKK::realloc_kokkos(fr, "grace_fs:fr", natom, maxneigh, lmax + 1, nradmax);
    MemKK::realloc_kokkos(dfr, "grace_fs:dfr", natom, maxneigh, lmax + 1, nradmax);
    MemKK::realloc_kokkos(gr, "grace_fs:gr", natom, maxneigh, nradbase);
    MemKK::realloc_kokkos(dgr, "grace_fs:dgr", natom, maxneigh, nradbase);
    const int max_num_functions = MAX(nradbase, nradmax*(lmax + 1));
    MemKK::realloc_kokkos(d_values, "grace_fs:d_values", natom, maxneigh, max_num_functions);
    MemKK::realloc_kokkos(d_derivatives, "grace_fs:d_derivatives", natom, maxneigh, max_num_functions);
    MemKK::realloc_kokkos(d_ncount, "grace_fs:ncount", natom);
    MemKK::realloc_kokkos(d_mu, "grace_fs:mu", natom, maxneigh);
    MemKK::realloc_kokkos(d_rhats, "grace_fs:rhats", natom, maxneigh);
    MemKK::realloc_kokkos(d_rnorms, "grace_fs:rnorms", natom, maxneigh);
    MemKK::realloc_kokkos(d_nearest, "grace_fs:nearest", natom, maxneigh);
    MemKK::realloc_kokkos(f_ij, "grace_fs:f_ij", natom, maxneigh);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::copy_pertype()
{
  auto basis_set = aceimpl->basis_set;

  MemKK::realloc_kokkos(d_E0vals, "grace_fs:E0vals", nelements);
  MemKK::realloc_kokkos(d_ndensity, "grace_fs:ndensity", nelements);
  MemKK::realloc_kokkos(d_npoti, "grace_fs:npoti", nelements);

  auto h_E0vals = Kokkos::create_mirror_view(d_E0vals);
  auto h_ndensity = Kokkos::create_mirror_view(d_ndensity);
  auto h_npoti = Kokkos::create_mirror_view(d_npoti);

  for (int n = 0; n < nelements; n++) {
    h_E0vals(n) = basis_set->E0_shift[n];
    h_ndensity(n) = basis_set->embedding_specifications.ndensity;
    string npoti = basis_set->embedding_specifications.type;
    if (npoti == "FinnisSinclair")
      h_npoti(n) = FS;
    else if (npoti == "FinnisSinclairShiftedScaled")
      h_npoti(n) = FS_SHIFTEDSCALED;
  }

  Kokkos::deep_copy(d_E0vals, h_E0vals);
  Kokkos::deep_copy(d_ndensity, h_ndensity);
  Kokkos::deep_copy(d_npoti, h_npoti);

  MemKK::realloc_kokkos(d_wpre, "grace_fs:wpre", nelements, basis_set->ndensitymax);
  MemKK::realloc_kokkos(d_mexp, "grace_fs:mexp", nelements, basis_set->ndensitymax);

  auto h_wpre = Kokkos::create_mirror_view(d_wpre);
  auto h_mexp = Kokkos::create_mirror_view(d_mexp);

  for (int n = 0; n < nelements; n++) {
    const int ndensity = basis_set->embedding_specifications.ndensity;
    for (int p = 0; p < ndensity; p++) {
      h_wpre(n, p) = basis_set->embedding_specifications.FS_parameters[p * 2 + 0];
      h_mexp(n, p) = basis_set->embedding_specifications.FS_parameters[p * 2 + 1];
    }
  }

  Kokkos::deep_copy(d_wpre, h_wpre);
  Kokkos::deep_copy(d_mexp, h_mexp);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::copy_splines()
{
  auto basis_set = aceimpl->basis_set;
  deallocate_views_of_views();

  k_splines_gk = Kokkos::DualView<SplineInterpolatorKokkos*, DeviceType>("grace_fs:splines_gk", nelements);
  k_splines_rnl = Kokkos::DualView<SplineInterpolatorKokkos*, DeviceType>("grace_fs:splines_rnl", nelements);

  for (int i = 0; i < nelements; i++) {
    k_splines_gk.view_host()(i) = basis_set->radial_functions.splines_gk;
    k_splines_rnl.view_host()(i) = basis_set->radial_functions.splines_rnl;
  }

  k_splines_gk.modify_host();
  k_splines_rnl.modify_host();
  k_splines_gk.sync_device();
  k_splines_rnl.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::copy_radial_Z()
{
  auto basis_set = aceimpl->basis_set;

  MemKK::realloc_kokkos(d_Z, "grace_fs:Z", nelements, nradmax);
  auto h_Z = Kokkos::create_mirror_view(d_Z);

  for (int mu_j = 0; mu_j < nelements; mu_j++) {
    for (int n = 0; n < nradmax; n++) {
      h_Z(mu_j, n) = basis_set->radial_functions.Z(mu_j, n);
    }
  }

  Kokkos::deep_copy(d_Z, h_Z);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::copy_tilde()
{
  auto basis_set = aceimpl->basis_set;

  idx_ms_combs_max = 0;
  total_num_functions_max = 0;

  MemKK::realloc_kokkos(d_idx_ms_combs_count, "grace_fs:idx_ms_combs_count", nelements);
  auto h_idx_ms_combs_count = Kokkos::create_mirror_view(d_idx_ms_combs_count);

  MemKK::realloc_kokkos(d_total_basis_size, "grace_fs:total_basis_size", nelements);
  auto h_total_basis_size = Kokkos::create_mirror_view(d_total_basis_size);

  for (int mu = 0; mu < nelements; mu++) {
    int idx_ms_combs = 0;
    const int total_basis_size = basis_set->basis[mu].size();

    for (int idx_func = 0; idx_func < total_basis_size; ++idx_func) {
      auto &func = basis_set->basis[mu][idx_func];
      for (int ms_ind = 0; ms_ind < func.num_ms_combs; ++ms_ind)
        idx_ms_combs++;
    }
    h_idx_ms_combs_count(mu) = idx_ms_combs;
    idx_ms_combs_max = MAX(idx_ms_combs_max, idx_ms_combs);
    total_num_functions_max = MAX(total_num_functions_max, total_basis_size);
    h_total_basis_size(mu) = total_basis_size;
  }

  Kokkos::deep_copy(d_idx_ms_combs_count, h_idx_ms_combs_count);
  Kokkos::deep_copy(d_total_basis_size, h_total_basis_size);

  MemKK::realloc_kokkos(d_rank, "grace_fs:rank", nelements, total_num_functions_max);
  MemKK::realloc_kokkos(d_num_ms_combs, "grace_fs:num_ms_combs", nelements, total_num_functions_max);
  MemKK::realloc_kokkos(d_idx_funcs, "grace_fs:idx_funcs", nelements, idx_ms_combs_max);
  MemKK::realloc_kokkos(d_ns, "grace_fs:ns", nelements, total_num_functions_max);
  MemKK::realloc_kokkos(d_ls, "grace_fs:ls", nelements, total_num_functions_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ms_combs, "grace_fs:ms_combs", nelements, idx_ms_combs_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_gen_cgs, "grace_fs:gen_cgs", nelements, idx_ms_combs_max);
  MemKK::realloc_kokkos(d_coeffs, "grace_fs:coeffs", nelements, total_num_functions_max, basis_set->ndensitymax);

  auto h_rank = Kokkos::create_mirror_view(d_rank);
  auto h_num_ms_combs = Kokkos::create_mirror_view(d_num_ms_combs);
  auto h_idx_funcs = Kokkos::create_mirror_view(d_idx_funcs);
  auto h_ns = Kokkos::create_mirror_view(d_ns);
  auto h_ls = Kokkos::create_mirror_view(d_ls);
  auto h_ms_combs = Kokkos::create_mirror_view(d_ms_combs);
  auto h_gen_cgs = Kokkos::create_mirror_view(d_gen_cgs);
  auto h_coeffs = Kokkos::create_mirror_view(d_coeffs);

  const int ndensity = basis_set->embedding_specifications.ndensity;

  for (int mu = 0; mu < nelements; mu++) {
    const int total_basis_size = basis_set->basis[mu].size();
    int idx_ms_combs = 0;

    for (int idx_func = 0; idx_func < total_basis_size; ++idx_func) {
      auto &func = basis_set->basis[mu][idx_func];
      const int rank = func.rank;

      h_rank(mu, idx_func) = rank;
      h_num_ms_combs(mu, idx_func) = func.num_ms_combs;
      h_ns(mu, idx_func) = func.ns;

      for (int t = 0; t < rank; t++)
        h_ls(mu, idx_func, t) = func.ls(t);

      for (int p = 0; p < ndensity; ++p)
        h_coeffs(mu, idx_func, p) = func.coeff(p);

      for (int ms_ind = 0; ms_ind < func.num_ms_combs; ++ms_ind) {
        auto ms = &func.ms_combs(ms_ind * rank);
        for (int t = 0; t < rank; t++)
          h_ms_combs(mu, idx_ms_combs, t) = ms[t];

        h_gen_cgs(mu, idx_ms_combs) = func.gen_cgs(ms_ind);
        h_idx_funcs(mu, idx_ms_combs) = idx_func;
        idx_ms_combs++;
      }
    }
  }

  Kokkos::deep_copy(d_rank, h_rank);
  Kokkos::deep_copy(d_num_ms_combs, h_num_ms_combs);
  Kokkos::deep_copy(d_idx_funcs, h_idx_funcs);
  Kokkos::deep_copy(d_ns, h_ns);
  Kokkos::deep_copy(d_ls, h_ls);
  Kokkos::deep_copy(d_ms_combs, h_ms_combs);
  Kokkos::deep_copy(d_gen_cgs, h_gen_cgs);
  Kokkos::deep_copy(d_coeffs, h_coeffs);

  // ASI for extrapolation grade
  // find max dimension for ASI
  int max_asi_dim = total_num_functions_max; // Use the max basis size found
  
  // Allocate d_ASI
  // We use a temporary non-const view to copy data, then assign to const d_ASI (if d_ASI is const in header)
  // Header defines: tc_ace_3d d_ASI; which is const.
  // We need to verify if we can assign a non-const view to a const view member? Yes.
  // But to allocate we need the non-const view.
  t_ace_3d d_ASI_temp("grace_fs:ASI", nelements, max_asi_dim, max_asi_dim);
  auto h_ASI = Kokkos::create_mirror_view(d_ASI_temp);
  Kokkos::deep_copy(h_ASI, 0.0);

  for (int mu = 0; mu < nelements; mu++) {
     if (aceimpl->ace->A_active_set_inv.count(mu) > 0) {
        const auto &A_as_inv = aceimpl->ace->A_active_set_inv.at(mu);
        int rows = A_as_inv.get_dim(0);
        int cols = A_as_inv.get_dim(1); // Should match basis size
        
        // Check bounds
        if (rows > max_asi_dim || cols > max_asi_dim) {
            // Should not happen if total_basis_size matches
        }

        // Copy and transpose: h_ASI(mu, col, row) = A_as_inv(row, col)
        // We want d_ASI(mu, k, j) where k sums with projections.
        // A_as_inv(i, k) in CPU code means i=gamma_idx, k=basis_idx.
        // So CPU: gamma[i] += proj[k] * A(i, k).
        // We want: gamma[i] += proj[k] * d_ASI(mu, k, i).
        // So d_ASI(mu, k, i) = A(i, k).
        // So we copy A(row, col) to h_ASI(mu, col, row).
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                h_ASI(mu, c, r) = A_as_inv(r, c);
            }
        }
     }
  }
  Kokkos::deep_copy(d_ASI_temp, h_ASI);
  d_ASI = d_ASI_temp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::pre_compute_harmonics(int lmax)
{
  auto h_idx_sph = Kokkos::create_mirror_view(d_idx_sph);
  auto h_alm = Kokkos::create_mirror_view(alm);
  auto h_blm = Kokkos::create_mirror_view(blm);
  auto h_cl = Kokkos::create_mirror_view(cl);
  auto h_dl = Kokkos::create_mirror_view(dl);

  Kokkos::deep_copy(h_idx_sph,-1);

  int idx_sph = 0;
  for (int m = 0; m <= lmax; m++) {
    const double msq = m * m;
    for (int l = m; l <= lmax; l++) {
      const int idx = l * (l + 1) + m;
      h_idx_sph(idx) = idx_sph;

      double a = 0.0, b = 0.0;
      // CPU condition: for (MS_TYPE m = 0; m < l - 1; m++) in loop l=1..lmax
      // This means alm(l, m) is computed when m < l - 1, i.e., l > m + 1
      // Also l must be > 1 for the recursion to work
      if (l > 1 && m < l - 1) {
        const double lsq = l * l;
        const double ld = 2 * l;
        const double l1 = (4 * lsq - 1);
        const double l2 = lsq - ld + 1;
        a = sqrt((double(l1)) / (double(lsq - msq)));
        b = -sqrt((double(l2 - msq)) / (double(4 * l2 - 1)));
      }
      h_alm(idx_sph) = a;
      h_blm(idx_sph) = b;
      idx_sph++;
    }
  }
  idx_sph_max = idx_sph;

  for (int l = 1; l <= lmax; l++) {
    h_cl(l) = -sqrt(1.0 + 0.5 / (double(l)));
    h_dl(l) = sqrt(double(2 * (l - 1) + 3));
  }

  Kokkos::deep_copy(d_idx_sph, h_idx_sph);
  Kokkos::deep_copy(alm, h_alm);
  Kokkos::deep_copy(blm, h_blm);
  Kokkos::deep_copy(cl, h_cl);
  Kokkos::deep_copy(dl, h_dl);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (host_flag) {
    atomKK->sync(Host,X_MASK|TYPE_MASK);
    PairGRACEFS::compute(eflag_in,vflag_in);
    atomKK->modified(Host,F_MASK);
    return;
  }

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  if (flag_compute_extrapolation_grade && atom->nlocal > nmax) {
    memory->destroy(extrapolation_grade_gamma);
    nmax = atom->nlocal;
    memory->create(extrapolation_grade_gamma, nmax, "grace_fs/atom:gamma");
    memset(extrapolation_grade_gamma, 0, nmax * sizeof(*extrapolation_grade_gamma));
  }

  copymode = 1;
  if (!force->newton_pair)
    error->all(FLERR,"PairGRACEFSKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_scale.template sync<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  int vector_length_default = 1;
  int team_size_default = 1;
  if(Kokkos::DefaultExecutionSpace::concurrency() > 1)
    team_size_default = 32;

  k_splines_gk.sync_device();
  k_splines_rnl.sync_device();
  d_splines_gk = k_splines_gk.view_device();
  d_splines_rnl = k_splines_rnl.view_device();

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  maxneigh = 0;
  Kokkos::parallel_reduce("grace_fs::find_maxneigh", inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  chunk_size = MIN(chunksize,inum);
  chunk_offset = 0;

  grow(chunk_size, maxneigh);

  EV_FLOAT ev;

  while (chunk_offset < inum) {
    Kokkos::deep_copy(weights, 0.0);
    Kokkos::deep_copy(A, 0.0);
    Kokkos::deep_copy(rhos, 0.0);
    Kokkos::deep_copy(projections, 0.0);
    Kokkos::deep_copy(d_gamma, 0.0);

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairGRACEFSComputeNeigh>(chunk_size,team_size,vector_length);
      int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
      typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeNeigh> policy_neigh(chunk_size,team_size,vector_length);
      policy_neigh = policy_neigh.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
    }

    //ComputeRadial
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairGRACEFSComputeRadial>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeRadial> policy_radial(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeRadial",policy_radial,*this);
    }

    //ComputeAi
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairGRACEFSComputeAi>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      int plm_size = (lmax + 1) * (lmax + 2) / 2;
      int scratch_size = scratch_size_helper<KK_FLOAT>(plm_size);
      typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeAi> policy_ai(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      policy_ai = policy_ai.set_scratch_size(0, Kokkos::PerThread(scratch_size));
      Kokkos::parallel_for("ComputeAi",policy_ai,*this);
    }

    //ComputeRho
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeRho> policy_rho(0, chunk_size * idx_ms_combs_max);
      Kokkos::parallel_for("ComputeRho",policy_rho,*this);
    }

    //ComputeFS
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeFS> policy_fs(0,chunk_size);
      Kokkos::parallel_for("ComputeFS",policy_fs,*this);
    }
    

    //ComputeGamma[OPTIONAL]
    if (flag_compute_extrapolation_grade) {
      typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeGamma> policy_gamma(0,chunk_size);
      Kokkos::parallel_for("ComputeGamma",policy_gamma,*this);
    }

    //ComputeWeights
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeWeights> policy_weights(0,chunk_size * idx_ms_combs_max);
      Kokkos::parallel_for("ComputeWeights",policy_weights,*this);
    }

    //ComputeDerivative
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairGRACEFSComputeDerivative>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      int plm_size = (lmax + 1) * (lmax + 2) / 2;
      int scratch_size = scratch_size_helper<KK_FLOAT>(2 * plm_size);
      typename Kokkos::TeamPolicy<DeviceType,TagPairGRACEFSComputeDerivative> policy_derivative(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      policy_derivative = policy_derivative.set_scratch_size(0, Kokkos::PerThread(scratch_size));
      Kokkos::parallel_for("ComputeDerivative",policy_derivative,*this);
    }

    //ComputeForce
    {
      if (neighflag == HALF) {
        if (evflag) {
          typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeForce<HALF,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce("ComputeForce", policy_force, *this, ev_tmp);
        } else {
          typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeForce<HALF,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for("ComputeForce", policy_force, *this);
        }
      } else if (neighflag == HALFTHREAD) {
        if (evflag) {
          typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeForce<HALFTHREAD,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce("ComputeForce", policy_force, *this, ev_tmp);
        } else {
          typename Kokkos::RangePolicy<DeviceType,TagPairGRACEFSComputeForce<HALFTHREAD,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for("ComputeForce", policy_force, *this);
        }
      }
    }

    if (eflag_global) eng_vdwl += ev_tmp.evdwl;
    if (vflag_global) {
      for (int k = 0; k < 6; k++) virial[k] += ev_tmp.v[k];
    }

    if (flag_compute_extrapolation_grade){
      h_gamma = Kokkos::create_mirror_view(d_gamma);
      Kokkos::deep_copy(h_gamma, d_gamma);
      memcpy(extrapolation_grade_gamma+chunk_offset, (void *) h_gamma.data(), sizeof(double)*chunk_size);
    }

    Kokkos::fence();

    chunk_offset += chunk_size;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  atomKK->modified(execution_space,F_MASK);
  copymode = 0;

  if (need_dup) {
    dup_f = {};
    dup_vatom = {};
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::evaluate_splines(const int ii, const int jj, KK_FLOAT r, int /*nradbase_c*/, int /*nradial_c*/, int mu_i, int mu_j) const
{
  // Radial functions for B basis (rnl)
  d_splines_rnl(0).calcSplines(ii, jj, r, d_values, d_derivatives);
  for (int ll = 0; ll < (int)fr.extent(2); ll++) {
    for (int kk = 0; kk < (int)fr.extent(3); kk++) {
      const KK_FLOAT Z_val = d_Z(mu_j, kk);
      const int flatten = kk * (lmax + 1) + ll;
      fr(ii, jj, ll, kk) = d_values(ii, jj, flatten) * Z_val * nnorm;
      dfr(ii, jj, ll, kk) = d_derivatives(ii, jj, flatten) * Z_val * nnorm;
    }
  }

  // Radial functions for density (gk)
  d_splines_gk(0).calcSplines(ii, jj, r, d_values, d_derivatives);
  for (int kk = 0; kk < (int)gr.extent(2); kk++) {
    gr(ii, jj, kk) = d_values(ii, jj, kk);
    dgr(ii, jj, kk) = d_derivatives(ii, jj, kk);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type(i);
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  // scratch memory for filtering neighbors
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * maxneigh;
  int* inside = (int*)team.team_shmem().get_shmem(team.team_size() * maxneigh * sizeof(int), 0) + scratch_shift;

  int ncount = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& count) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    inside[jj] = -1;
    if (rsq < d_cutsq(itype,jtype)) {
     inside[jj] = 1;
     count++;
    }
  },ncount);

  d_ncount(ii) = ncount;

  Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& offset, bool final) {
    if (inside[jj] < 0) return;
    if (final) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const KK_FLOAT delx = xtmp - x(j,0);
      const KK_FLOAT dely = ytmp - x(j,1);
      const KK_FLOAT delz = ztmp - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      const KK_FLOAT r = sqrt(rsq);
      const KK_FLOAT rinv = 1.0/r;
      const int mu_j = d_map(type(j));
      d_mu(ii,offset) = mu_j;
      d_rnorms(ii,offset) = r;
      d_rhats(ii,offset,0) = -delx*rinv;
      d_rhats(ii,offset,1) = -dely*rinv;
      d_rhats(ii,offset,2) = -delz*rinv;
      d_nearest(ii,offset) = j;
    }
    offset++;
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeRadial, const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeRadial>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const KK_FLOAT r_norm = d_rnorms(ii, jj);
  const int mu_i = d_map(type(i));
  const int mu_j = d_mu(ii, jj);

  // Note: nradbase and nradmax are member variables of the class, accessible inside operator
  evaluate_splines(ii, jj, r_norm, 0, 0, mu_i, mu_j);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeAi, const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeAi>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const int mu_j = d_mu(ii, jj);

  // Get unit vector components
  const KK_FLOAT rx = d_rhats(ii, jj, 0);
  const KK_FLOAT ry = d_rhats(ii, jj, 1);
  const KK_FLOAT rz = d_rhats(ii, jj, 2);

  static constexpr KK_FLOAT sq2 = 1.4142135623730950488;

  int plm_size = (lmax + 1) * (lmax + 2) / 2;
  KK_FLOAT* plm = (KK_FLOAT*) team.thread_scratch(0).get_shmem(this->scratch_size_helper<KK_FLOAT>(plm_size));

  // =====================================================
  // STEP 1: Compute associated Legendre polynomials (barplm)
  // Following CPU ace_spherical_cart.cpp::compute_barplm exactly
  // =====================================================

  // l=0, m=0
  plm[0] = Y00_kk;

  if (lmax > 0) {
    // l=1, m=0: index = 1
    plm[1] = Y00_kk * sq3_kk * rz;

    // l=1, m=1: index = 2
    plm[2] = -sq3o2_kk * Y00_kk;

    // loop l = 2, lmax
    for (int l = 2; l <= lmax; l++) {
      // m = 0 to l-2
      for (int m = 0; m < l - 1; m++) {
        const int idx = l * (l + 1) / 2 + m;
        const int idx_l1_m = (l - 1) * l / 2 + m;
        const int idx_l2_m = (l - 2) * (l - 1) / 2 + m;
        // alm, blm are stored as alm(l, m) -> flat index l*(lmax+1)+m
        // But in pre_compute_harmonics they are stored by idx_sph
        // Need to get the correct alm/blm for this (l, m)
        const int alm_idx = d_idx_sph(l * (l + 1) + m);
        const KK_FLOAT a = alm(alm_idx);
        const KK_FLOAT b = blm(alm_idx);
        plm[idx] = a * (rz * plm[idx_l1_m] + b * plm[idx_l2_m]);
      }
      // m = l-1
      {
        const int idx = l * (l + 1) / 2 + l - 1;
        const int idx_l1_l1 = (l - 1) * l / 2 + l - 1;
        const KK_FLOAT t = dl(l) * plm[idx_l1_l1];
        plm[idx] = t * rz;
      }
      // m = l
      {
        const int idx = l * (l + 1) / 2 + l;
        const int idx_l1_l1 = (l - 1) * l / 2 + l - 1;
        plm[idx] = cl(l) * plm[idx_l1_l1];
      }
    }
  }

  // =====================================================
  // STEP 2: Compute complex Ylm and convert to real Ylm
  // Following CPU ace_spherical_cart.cpp::compute_real_ylm
  // =====================================================

  // Complex phase = (rx, ry)
  const KK_FLOAT phase_re = rx;
  const KK_FLOAT phase_im = ry;

  // m = 0: real_ylm(l, 0) = plm(l, 0)
  for (int l = 0; l <= lmax; l++) {
    const int plm_idx = l * (l + 1) / 2 + 0;
    const KK_FLOAT real_ylm_0 = plm[plm_idx];
    const int A_idx = l * (l + 1) + 0; // m = 0

    for (int n = 0; n < nradmax; n++) {
      Kokkos::atomic_add(&A(ii, A_idx, n), fr(ii, jj, l, n) * real_ylm_0);
    }
  }

  // m = 1: ylm = phase * plm
  if (lmax >= 1) {
    for (int l = 1; l <= lmax; l++) {
      const int plm_idx = l * (l + 1) / 2 + 1;
      const KK_FLOAT plm_val = plm[plm_idx];

      // ylm.real = phase.real * plm, ylm.imag = phase.imag * plm
      const KK_FLOAT ylm_re = phase_re * plm_val;
      const KK_FLOAT ylm_im = phase_im * plm_val;

      // factor = (-1)^m
      const int factor = -1; // m=1

      // real_ylm(l, +m) = sq2 * factor * ylm.real
      // real_ylm(l, -m) = sq2 * factor * ylm.imag
      const KK_FLOAT real_ylm_pos = sq2 * factor * ylm_re;
      const KK_FLOAT real_ylm_neg = sq2 * factor * ylm_im;

      const int A_idx_pos = l * (l + 1) + 1;  // m = +1
      const int A_idx_neg = l * (l + 1) - 1;  // m = -1

      for (int n = 0; n < nradmax; n++) {
        const KK_FLOAT R = fr(ii, jj, l, n);
        Kokkos::atomic_add(&A(ii, A_idx_pos, n), R * real_ylm_pos);
        Kokkos::atomic_add(&A(ii, A_idx_neg, n), R * real_ylm_neg);
      }
    }
  }

  // m >= 2: phasem = phase^m
  KK_FLOAT phasem_re = phase_re;
  KK_FLOAT phasem_im = phase_im;

  for (int m = 2; m <= lmax; m++) {
    // Update phasem = phasem * phase
    const KK_FLOAT tmp_re = phasem_re * phase_re - phasem_im * phase_im;
    const KK_FLOAT tmp_im = phasem_re * phase_im + phasem_im * phase_re;
    phasem_re = tmp_re;
    phasem_im = tmp_im;

    // factor = (-1)^m
    const int factor = (m % 2 == 0) ? 1 : -1;

    for (int l = m; l <= lmax; l++) {
      const int plm_idx = l * (l + 1) / 2 + m;
      const KK_FLOAT plm_val = plm[plm_idx];

      // ylm = phasem * plm
      const KK_FLOAT ylm_re = phasem_re * plm_val;
      const KK_FLOAT ylm_im = phasem_im * plm_val;

      // real_ylm(l, +m) = sq2 * factor * ylm.real
      // real_ylm(l, -m) = sq2 * factor * ylm.imag
      const KK_FLOAT real_ylm_pos = sq2 * factor * ylm_re;
      const KK_FLOAT real_ylm_neg = sq2 * factor * ylm_im;

      const int A_idx_pos = l * (l + 1) + m;
      const int A_idx_neg = l * (l + 1) - m;

      for (int n = 0; n < nradmax; n++) {
        const KK_FLOAT R = fr(ii, jj, l, n);
        Kokkos::atomic_add(&A(ii, A_idx_pos, n), R * real_ylm_pos);
        Kokkos::atomic_add(&A(ii, A_idx_neg, n), R * real_ylm_neg);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeRho, const int& iter) const
{
  // Decode ii and idx from iter
  const int ii = iter / idx_ms_combs_max;
  const int idx = iter % idx_ms_combs_max;
  if (ii >= chunk_size) return;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));
  if (idx >= d_idx_ms_combs_count(mu_i)) return;

  const int idx_func = d_idx_funcs(mu_i, idx);
  const int rank = d_rank(mu_i, idx_func);
  const int ns = d_ns(mu_i, idx_func);
  const int ndensity = d_ndensity(mu_i);

  // Compute product of A values
  KK_FLOAT val = 1.0;

  for (int t = 0; t < rank; t++) {
    const int l = d_ls(mu_i, idx_func, t);
    const int m = d_ms_combs(mu_i, idx, t);
    const int idx_sph = l * (l + 1) + m;
    KK_FLOAT a_val = A(ii, idx_sph, ns - 1);
    val *= a_val;
  }

  val *= d_gen_cgs(mu_i, idx);

  // Accumulate to rhos
  for (int p = 0; p < ndensity; p++) {
    const KK_FLOAT coeff = d_coeffs(mu_i, idx_func, p);
    Kokkos::atomic_add(&rhos(ii, p), val * coeff);
  }

  if (flag_compute_extrapolation_grade) {
     Kokkos::atomic_add(&projections(ii, idx_func), val);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeFS, const int& ii) const
{
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  // Compute FS embedding (evdwl is computed internally, stored in e_atom(ii) temporarily)
  KK_FLOAT evdwl = 0.0;
  FS_values_and_derivatives(ii, evdwl, mu_i);

  // Apply scale, shift, and E0_shift like CPU at line 664:
  // e_atom = basis_set.scale * evdwl + basis_set.shift + basis_set.E0_shift.at(mu_i);
  e_atom(ii) = energy_scale * evdwl + energy_shift + d_E0vals(mu_i);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeGamma, const int& ii) const
{
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));
  
  // Initialize projections and gamma
  // d_gamma(ii) = 0.0; // Already init in compute loop
  // projections is reused, need to zero out?
  // Parallel_for "ComputeGamma" is per atom.
  // projections size: [natom, total_num_functions_max]
  
  const int total_basis = d_total_basis_size(mu_i);
  const int ms_combs_count = d_idx_ms_combs_count(mu_i);
  
  // Zero out projections for this atom
  // 1. Compute basis function values (B) and store in projections
  // REMOVED: Now computed in ComputeRho

  // 2. Compute Gamma = max | projections * ASI |
  // ASI is transposed in d_ASI(mu, k, j) where k is basis idx, j is gamma idx.
  // d_ASI shape: [nelements, max_dim, max_dim]
  // We iterate j (gamma component) and sum over k (basis).
  // Need to know number of gamma components (rows of ASI). 
  // We don't distinctly store 'rows of ASI' on device.
  // Assume it's same as basis size? ASI is Inverse of Active Set, so square matrix.
  // rows = cols = total_basis (approx).
  // Safest to iterate up to total_basis.
  
  KK_FLOAT max_gamma = 0.0;
  
  for (int j = 0; j < total_basis; j++) {
      KK_FLOAT current_gamma = 0.0;
      for (int k = 0; k < total_basis; k++) {
          current_gamma += projections(ii, k) * d_ASI(mu_i, k, j);
      }
      if (abs(current_gamma) > max_gamma) max_gamma = abs(current_gamma);
  }
  
  d_gamma(ii) = max_gamma;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeWeights, const int& iter) const
{
  const int ii = iter / idx_ms_combs_max;
  const int idx = iter % idx_ms_combs_max;
  if (ii >= chunk_size) return;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));
  if (idx >= d_idx_ms_combs_count(mu_i)) return;

  const int idx_func = d_idx_funcs(mu_i, idx);
  const int rank = d_rank(mu_i, idx_func);
  const int ns = d_ns(mu_i, idx_func);
  const int ndensity = d_ndensity(mu_i);

  // Compute dB/dA for each A in the product
  for (int t = 0; t < rank; t++) {
    const int l_t = d_ls(mu_i, idx_func, t);
    const int m_t = d_ms_combs(mu_i, idx, t);
    const int idx_sph_t = l_t * (l_t + 1) + m_t;

    // Compute product of A's except t-th
    KK_FLOAT val = 1.0;

    for (int s = 0; s < rank; s++) {
      if (s == t) continue;
      const int l_s = d_ls(mu_i, idx_func, s);
      const int m_s = d_ms_combs(mu_i, idx, s);
      const int idx_sph_s = l_s * (l_s + 1) + m_s;
      KK_FLOAT a_val = A(ii, idx_sph_s, ns - 1);
      val *= a_val;
    }

    val *= d_gen_cgs(mu_i, idx);

    // Accumulate weights with dF/drho
    for (int p = 0; p < ndensity; p++) {
      const KK_FLOAT coeff = d_coeffs(mu_i, idx_func, p);
      const KK_FLOAT dF = dF_drho(ii, p);
      const KK_FLOAT weight_contrib = dF * coeff * val;
      Kokkos::atomic_add(&weights(ii, idx_sph_t, ns - 1), weight_contrib);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeDerivative, const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeDerivative>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const KK_FLOAT rx = d_rhats(ii, jj, 0);
  const KK_FLOAT ry = d_rhats(ii, jj, 1);
  const KK_FLOAT rz = d_rhats(ii, jj, 2);
  const KK_FLOAT rinv = 1.0 / d_rnorms(ii, jj);  // Fix: rinv should be 1/r, not r

  KK_ACC_FLOAT f_ji[3];
  f_ji[0] = f_ji[1] = f_ji[2] = 0.0;

  static constexpr KK_FLOAT sq2 = 1.4142135623730950488;

  // =====================================================
  // STEP 1: Compute Plm and dPlm (associated Legendre polynomials)
  // =====================================================
  int plm_size = (lmax + 1) * (lmax + 2) / 2;
  KK_FLOAT* scratch = (KK_FLOAT*) team.thread_scratch(0).get_shmem(this->scratch_size_helper<KK_FLOAT>(2 * plm_size));
  KK_FLOAT* plm = scratch;
  KK_FLOAT* dplm = scratch + plm_size;

  plm[0] = Y00_kk;
  dplm[0] = 0.0;

  if (lmax > 0) {
    plm[1] = Y00_kk * sq3_kk * rz;
    dplm[1] = Y00_kk * sq3_kk;
    plm[2] = -sq3o2_kk * Y00_kk;
    dplm[2] = 0.0;

    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        const int idx = l * (l + 1) / 2 + m;
        const int idx_l1_m = (l - 1) * l / 2 + m;
        const int idx_l2_m = (l - 2) * (l - 1) / 2 + m;
        const int alm_idx = d_idx_sph(l * (l + 1) + m);
        const KK_FLOAT a = alm(alm_idx);
        const KK_FLOAT b = blm(alm_idx);
        plm[idx] = a * (rz * plm[idx_l1_m] + b * plm[idx_l2_m]);
        dplm[idx] = a * (plm[idx_l1_m] + rz * dplm[idx_l1_m] + b * dplm[idx_l2_m]);
      }
      {
        const int idx = l * (l + 1) / 2 + l - 1;
        const int idx_l1_l1 = (l - 1) * l / 2 + l - 1;
        const KK_FLOAT t = dl(l) * plm[idx_l1_l1];
        plm[idx] = t * rz;
        dplm[idx] = t;
      }
      {
        const int idx = l * (l + 1) / 2 + l;
        const int idx_l1_l1 = (l - 1) * l / 2 + l - 1;
        plm[idx] = cl(l) * plm[idx_l1_l1];
        dplm[idx] = 0.0;
      }
    }
  }

  // =====================================================
  // STEP 2: Compute real Y_lm and DY_lm, then accumulate forces
  // Following CPU: for each (n, l, m=-l..+l), w * grad(R * Y)
  // =====================================================

  const KK_FLOAT phase_re = rx;
  const KK_FLOAT phase_im = ry;

  // Loop over (n, l) first, then m from -l to +l
  for (int n = 0; n < nradmax; n++) {
    for (int l = 0; l <= lmax; l++) {
      const KK_FLOAT R = fr(ii, jj, l, n);
      const KK_FLOAT DR = dfr(ii, jj, l, n);
      const KK_FLOAT R_over_r = R * rinv;

      // m = 0: real_ylm = plm, real_dylm = projection of dplm
      {
        const int plm_idx = l * (l + 1) / 2 + 0;
        const KK_FLOAT Y = plm[plm_idx];
        const KK_FLOAT dplm_val = dplm[plm_idx];

        // DY for m=0: dY/dr_hat projected to Cartesian
        const KK_FLOAT rdy = dplm_val * rz;
        const KK_FLOAT DY_x = -rdy * rx;
        const KK_FLOAT DY_y = -rdy * ry;
        const KK_FLOAT DY_z = dplm_val - rdy * rz;

        const int A_idx = l * (l + 1) + 0;
        const KK_FLOAT w = weights(ii, A_idx, n);
        if (w != 0.0) {
          const KK_FLOAT Y_DR = Y * DR;
          f_ji[0] += w * (Y_DR * rx + DY_x * R_over_r);
          f_ji[1] += w * (Y_DR * ry + DY_y * R_over_r);
          f_ji[2] += w * (Y_DR * rz + DY_z * R_over_r);
        }
      }

      // m >= 1: compute complex Ylm, convert to real, and derivatives
      KK_FLOAT phasem_re = phase_re;
      KK_FLOAT phasem_im = phase_im;

      for (int m = 1; m <= l; m++) {
        const int factor = (m % 2 == 0) ? 1 : -1;
        const KK_FLOAT m_kk = (KK_FLOAT)m;

        const int plm_idx = l * (l + 1) / 2 + m;
        const KK_FLOAT plm_val = plm[plm_idx];
        const KK_FLOAT dplm_val = dplm[plm_idx];

        // Complex ylm = phasem * plm
        const KK_FLOAT ylm_re = phasem_re * plm_val;
        const KK_FLOAT ylm_im = phasem_im * plm_val;

        // Real spherical harmonics
        const KK_FLOAT real_Y_pos = sq2 * factor * ylm_re;  // Y(l, +m)
        const KK_FLOAT real_Y_neg = sq2 * factor * ylm_im;  // Y(l, -m)

        // Complex derivatives - match CPU formula exactly
        // For m=1: dyx = (plm, 0), dyy = (0, plm)
        // For m>=2: dyx = m * phase^(m-1) * plm, dyy = i * dyx
        // dyz = phase^m * dplm for all m>=1
        
        KK_FLOAT dyx_re, dyx_im, dyy_re, dyy_im;
        const KK_FLOAT dyz_re = dplm_val * phasem_re;
        const KK_FLOAT dyz_im = dplm_val * phasem_im;
        
        if (m == 1) {
            // m=1: dyx = (plm, 0), dyy = (0, plm)
            dyx_re = plm_val;
            dyx_im = 0.0;
            dyy_re = 0.0;
            dyy_im = plm_val;
        } else {
            // m>=2: dyx = m * phase^(m-1) * plm, dyy = i * dyx
            // phasem = phase^m, so phase^(m-1) = phasem / phase
            // We compute m * phase^(m-1) * plm using: mphasem1 = m * phase^(m-1)
            // mphasem1_re = (phasem_re * phase_re + phasem_im * phase_im) / |phase|^2 * m
            // But since |phase|^2 = rx^2 + ry^2 = s2_safe, and we need phase^(m-1):
            // phase^(m-1) = phasem * conj(phase) / |phase|^2
            const KK_FLOAT s2 = rx * rx + ry * ry;
            if (s2 > 1e-12) {
                const KK_FLOAT inv_s2 = 1.0 / s2;
                // phase^(m-1) = phasem * conj(phase) / s2
                const KK_FLOAT pm1_re = (phasem_re * rx + phasem_im * ry) * inv_s2;
                const KK_FLOAT pm1_im = (phasem_im * rx - phasem_re * ry) * inv_s2;
                // mphasem1 = m * phase^(m-1) * plm
                const KK_FLOAT m_plm = m_kk * plm_val;
                dyx_re = m_plm * pm1_re;
                dyx_im = m_plm * pm1_im;
                // dyy = i * dyx = (-dyx_im, dyx_re)
                dyy_re = -dyx_im;
                dyy_im = dyx_re;
            } else {
                // Pole case: s2 ≈ 0
                dyx_re = dyx_im = dyy_re = dyy_im = 0.0;
            }
        }
        
        // rdy = rx * dyx + ry * dyy + rz * dyz (full projection!)
        const KK_FLOAT rdy_re = rx * dyx_re + ry * dyy_re + rz * dyz_re;
        const KK_FLOAT rdy_im = rx * dyx_im + ry * dyy_im + rz * dyz_im;
        
        // dylm = d? - rdy * r? for each component
        const KK_FLOAT dylm_x_re = dyx_re - rdy_re * rx;
        const KK_FLOAT dylm_x_im = dyx_im - rdy_im * rx;
        const KK_FLOAT dylm_y_re = dyy_re - rdy_re * ry;
        const KK_FLOAT dylm_y_im = dyy_im - rdy_im * ry;
        const KK_FLOAT dylm_z_re = dyz_re - rdy_re * rz;
        const KK_FLOAT dylm_z_im = dyz_im - rdy_im * rz;

        // Real DY
        const KK_FLOAT DY_pos_x = sq2 * factor * dylm_x_re;
        const KK_FLOAT DY_pos_y = sq2 * factor * dylm_y_re;
        const KK_FLOAT DY_pos_z = sq2 * factor * dylm_z_re;
        const KK_FLOAT DY_neg_x = sq2 * factor * dylm_x_im;
        const KK_FLOAT DY_neg_y = sq2 * factor * dylm_y_im;
        const KK_FLOAT DY_neg_z = sq2 * factor * dylm_z_im;

        // +m contribution
        const int A_idx_pos = l * (l + 1) + m;
        const KK_FLOAT w_pos = weights(ii, A_idx_pos, n);
          if (w_pos != 0.0) {
            const KK_FLOAT Y_DR = real_Y_pos * DR;
            f_ji[0] += w_pos * (Y_DR * rx + DY_pos_x * R_over_r);
            f_ji[1] += w_pos * (Y_DR * ry + DY_pos_y * R_over_r);
            f_ji[2] += w_pos * (Y_DR * rz + DY_pos_z * R_over_r);
          }

        // -m contribution
        const int A_idx_neg = l * (l + 1) - m;
        const KK_FLOAT w_neg = weights(ii, A_idx_neg, n);
        if (w_neg != 0.0) {
          const KK_FLOAT Y_DR = real_Y_neg * DR;
          f_ji[0] += w_neg * (Y_DR * rx + DY_neg_x * R_over_r);
          f_ji[1] += w_neg * (Y_DR * ry + DY_neg_y * R_over_r);
          f_ji[2] += w_neg * (Y_DR * rz + DY_neg_z * R_over_r);
        }

        // Update phasem for next m
        const KK_FLOAT tmp_re = phasem_re * phase_re - phasem_im * phase_im;
        const KK_FLOAT tmp_im = phasem_re * phase_im + phasem_im * phase_re;
        phasem_re = tmp_re;
        phasem_im = tmp_im;
      }
    }

  }  // end n loop

  f_ij(ii, jj, 0) = f_ji[0];
  f_ij(ii, jj, 1) = f_ji[1];
  f_ij(ii, jj, 2) = f_ji[2];
}



template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeForce<NEIGHFLAG,EVFLAG>, const int& ii) const
{
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int ncount = d_ncount(ii);

  for (int jj = 0; jj < ncount; jj++) {
    const int j = d_nearest(ii, jj);
    const KK_FLOAT fx = f_ij(ii, jj, 0);
    const KK_FLOAT fy = f_ij(ii, jj, 1);
    const KK_FLOAT fz = f_ij(ii, jj, 2);

    if (NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) {
      Kokkos::atomic_add(&f(i, 0), fx * energy_scale);
      Kokkos::atomic_add(&f(i, 1), fy * energy_scale);
      Kokkos::atomic_add(&f(i, 2), fz * energy_scale);
      Kokkos::atomic_add(&f(j, 0), -fx * energy_scale);
      Kokkos::atomic_add(&f(j, 1), -fy * energy_scale);
      Kokkos::atomic_add(&f(j, 2), -fz * energy_scale);
    } else {
      Kokkos::atomic_add(&f(i, 0), fx * energy_scale);
      Kokkos::atomic_add(&f(i, 1), fy * energy_scale);
      Kokkos::atomic_add(&f(i, 2), fz * energy_scale);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::operator() (TagPairGRACEFSComputeForce<NEIGHFLAG,EVFLAG>, const int& ii, EV_FLOAT& ev) const
{
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];
  const int ncount = d_ncount(ii);

  for (int jj = 0; jj < ncount; jj++) {
    const int j = d_nearest(ii, jj);
    const KK_FLOAT fx = f_ij(ii, jj, 0);
    const KK_FLOAT fy = f_ij(ii, jj, 1);
    const KK_FLOAT fz = f_ij(ii, jj, 2);

    if (NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) {
      Kokkos::atomic_add(&f(i, 0), fx * energy_scale);
      Kokkos::atomic_add(&f(i, 1), fy * energy_scale);
      Kokkos::atomic_add(&f(i, 2), fz * energy_scale);
      Kokkos::atomic_add(&f(j, 0), -fx * energy_scale);
      Kokkos::atomic_add(&f(j, 1), -fy * energy_scale);
      Kokkos::atomic_add(&f(j, 2), -fz * energy_scale);
    } else {
      // FULL neighbor list logic
      Kokkos::atomic_add(&f(i, 0), fx * energy_scale);
      Kokkos::atomic_add(&f(i, 1), fy * energy_scale);
      Kokkos::atomic_add(&f(i, 2), fz * energy_scale);
    }

    if (EVFLAG) {
      const KK_FLOAT rx = d_rhats(ii, jj, 0);
      const KK_FLOAT ry = d_rhats(ii, jj, 1);
      const KK_FLOAT rz = d_rhats(ii, jj, 2);
      const KK_FLOAT r = d_rnorms(ii, jj);
      const KK_FLOAT delx = -rx * r;
      const KK_FLOAT dely = -ry * r;
      const KK_FLOAT delz = -rz * r;
      const KK_FLOAT v0 = delx * fx * energy_scale;
      const KK_FLOAT v1 = dely * fy * energy_scale;
      const KK_FLOAT v2 = delz * fz * energy_scale;
      const KK_FLOAT v3 = delx * fy * energy_scale;
      const KK_FLOAT v4 = delx * fz * energy_scale;
      const KK_FLOAT v5 = dely * fz * energy_scale;

      if (NEIGHFLAG == HALF) {
        ev.v[0] += v0;
        ev.v[1] += v1;
        ev.v[2] += v2;
        ev.v[3] += v3;
        ev.v[4] += v4;
        ev.v[5] += v5;
      } else {
        ev.v[0] += 0.5 * v0;
        ev.v[1] += 0.5 * v1;
        ev.v[2] += 0.5 * v2;
        ev.v[3] += 0.5 * v3;
        ev.v[4] += 0.5 * v4;
        ev.v[5] += 0.5 * v5;
      }
    }
  }

  if (EVFLAG) {
    ev.evdwl += e_atom(ii);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::Fexp(const KK_FLOAT x, const KK_FLOAT m, KK_FLOAT &F, KK_FLOAT &DF) const
{
  const KK_FLOAT w = 1.e6;
  const KK_FLOAT eps = 1e-10;
  const KK_FLOAT lambda = pow(1.0 / w, m - 1.0);
  if (abs(x) > eps) {
    KK_FLOAT g;
    const KK_FLOAT a = abs(x);
    const KK_FLOAT am = pow(a, m);
    const KK_FLOAT w3x3 = pow(w * a, 3);
    const KK_FLOAT sign_factor = (signbit(x) ? -1 : 1);
    if (w3x3 > 30.0)
      g = 0.0;
    else
      g = exp(-w3x3);
    const KK_FLOAT omg = 1.0 - g;
    F = sign_factor * (omg * am + lambda * g * a);
    const KK_FLOAT dg = -3.0 * w * w * w * a * a * g;
    DF = m * pow(a, m - 1.0) * omg - am * dg + lambda * dg * a + lambda * g;
  } else {
    F = lambda * x;
    DF = lambda;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::FexpShiftedScaled(const KK_FLOAT rho, const KK_FLOAT mexp, KK_FLOAT &F, KK_FLOAT &DF) const
{
  const KK_FLOAT eps = 1e-10;
  if (abs(mexp - 1.0) < eps) {
    F = rho;
    DF = 1;
  } else {
    const KK_FLOAT a = abs(rho);
    const KK_FLOAT exprho = exp(-a);
    const KK_FLOAT nx = 1. / mexp;
    const KK_FLOAT xoff = pow(nx, (nx / (1.0 - nx))) * exprho;
    const KK_FLOAT yoff = pow(nx, (1 / (1.0 - nx))) * exprho;
    const KK_FLOAT sign_factor = (signbit(rho) ? -1 : 1);
    F = sign_factor * (pow(xoff + a, mexp) - yoff);
    DF = yoff + mexp * (-xoff + 1.0) * pow(xoff + a, mexp - 1.);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::FS_values_and_derivatives(const int ii, KK_FLOAT &evdwl, const int mu_i) const
{
  KK_FLOAT F, DF = 0;
  int npoti = d_npoti(mu_i);
  int ndensity = d_ndensity(mu_i);
  for (int p = 0; p < ndensity; p++) {
    const KK_FLOAT wpre = d_wpre(mu_i, p);
    const KK_FLOAT mexp = d_mexp(mu_i, p);
    if (npoti == FS)
      Fexp(rhos(ii, p), mexp, F, DF);
    else if (npoti == FS_SHIFTEDSCALED)
      FexpShiftedScaled(rhos(ii, p), mexp, F, DF);
    evdwl += F * wpre;
    dF_drho(ii, p) = DF * wpre;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairGRACEFSKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;
  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairGRACEFSKokkos<DeviceType>::check_team_size_reduce(int inum, int &team_size, int vector_length) {
  int team_size_max;
  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<typename scratch_type>
KOKKOS_INLINE_FUNCTION
int PairGRACEFSKokkos<DeviceType>::scratch_size_helper(int values_per_team) const {
  typedef Kokkos::View<scratch_type*, typename DeviceType::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;
  return ScratchViewType::shmem_size(values_per_team);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGRACEFSKokkos<DeviceType>::SplineInterpolatorKokkos::operator=(const SplineInterpolator &spline) {
  cutoff = spline.cutoff;
  deltaSplineBins = spline.deltaSplineBins;
  ntot = spline.ntot;
  nlut = spline.nlut;
  invrscalelookup = spline.invrscalelookup;
  rscalelookup = spline.rscalelookup;
  num_of_functions = spline.num_of_functions;

  lookupTable = t_ace_3d4_lr("lookupTable", ntot+1, num_of_functions);
  auto h_lookupTable = Kokkos::create_mirror_view(lookupTable);
  for (int i = 0; i < ntot+1; i++)
    for (int j = 0; j < num_of_functions; j++)
      for (int k = 0; k < 4; k++)
        h_lookupTable(i, j, k) = spline.lookupTable(i, j, k);
  Kokkos::deep_copy(lookupTable, h_lookupTable);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGRACEFSKokkos<DeviceType>::SplineInterpolatorKokkos::calcSplines(const int ii, const int jj, const KK_FLOAT r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const
{
  KK_FLOAT wl, wl2, wl3, w2l1, w3l2;
  KK_FLOAT c[4];
  KK_FLOAT x = r * rscalelookup;
  int nl = static_cast<int>(floor(x));

  if (nl <= 0)
    Kokkos::abort("Encountered very small distance. Stopping.");

  if (nl < nlut) {
    wl = x - KK_FLOAT(nl);
    wl2 = wl * wl;
    wl3 = wl2 * wl;
    w2l1 = 2.0 * wl;
    w3l2 = 3.0 * wl2;
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      for (int idx = 0; idx < 4; idx++)
        c[idx] = lookupTable(nl, func_id, idx);
      d_values(ii, jj, func_id) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
      d_derivatives(ii, jj, func_id) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;
    }
  } else { // fill with zeroes
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      d_values(ii, jj, func_id) = 0.0;
      d_derivatives(ii, jj, func_id) = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairGRACEFSKokkos<DeviceType>::memory_usage()

{
  double bytes = 0;
  bytes += MemKK::memory_usage(A);
  bytes += MemKK::memory_usage(e_atom);
  bytes += MemKK::memory_usage(rhos);
  bytes += MemKK::memory_usage(weights);
  bytes += MemKK::memory_usage(fr);
  bytes += MemKK::memory_usage(dfr);
  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairGRACEFSKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGRACEFSKokkos<LMPHostType>;
#endif
}
