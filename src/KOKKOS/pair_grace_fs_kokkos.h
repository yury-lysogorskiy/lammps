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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/fs/kk,PairGRACEFSKokkos<LMPDeviceType>);
PairStyle(grace/fs/kk/device,PairGRACEFSKokkos<LMPDeviceType>);
PairStyle(grace/fs/kk/host,PairGRACEFSKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRACE_FS_KOKKOS_H
#define LMP_PAIR_GRACE_FS_KOKKOS_H

#include "pair_grace_fs.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

class SplineInterpolator;

namespace LAMMPS_NS {

template<class DeviceType>
class PairGRACEFSKokkos : public PairGRACEFS {
 public:
  // Tags for Kokkos kernels
  struct TagPairGRACEFSComputeNeigh{};
  struct TagPairGRACEFSComputeRadial{};
  struct TagPairGRACEFSComputeAi{};
  struct TagPairGRACEFSComputeRho{};
  struct TagPairGRACEFSComputeFS{};
  struct TagPairGRACEFSComputeGamma{};
  struct TagPairGRACEFSComputeWeights{};
  struct TagPairGRACEFSComputeDerivative{};

  template<int NEIGHFLAG, int EVFLAG>
  struct TagPairGRACEFSComputeForce{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  using complex = SNAComplex<KK_FLOAT>;

  PairGRACEFSKokkos(class LAMMPS *);
  ~PairGRACEFSKokkos() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  // Kokkos kernel operators
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeRadial,const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeRadial>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeAi,const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeAi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeRho,const int& iter) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeFS,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeGamma, const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeWeights,const int& iter) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeDerivative,const typename Kokkos::TeamPolicy<DeviceType, TagPairGRACEFSComputeDerivative>::member_type& team) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairGRACEFSComputeForce<NEIGHFLAG,EVFLAG>,const int& ii, EV_FLOAT&) const;

 protected:
  int inum, maxneigh, chunk_size, chunk_offset, idx_ms_combs_max, total_num_functions_max, idx_sph_max;
  int host_flag;

  int eflag, vflag;

  int neighflag;
  int nelements, lmax, nradmax, nradbase;
  KK_FLOAT nnorm;  // normalization factor
  KK_FLOAT energy_scale;  // scale factor for energy output
  KK_FLOAT energy_shift;  // shift factor for energy output

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

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

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  friend void pair_virial_fdotr_compute<PairGRACEFSKokkos>(PairGRACEFSKokkos*);

  void grow(int, int);
  void copy_pertype();
  void copy_splines();
  void copy_radial_Z();
  void copy_tilde();
  void allocate() override;
  void precompute_harmonics();
  double memory_usage() override;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &fx, const KK_FLOAT &fy, const KK_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  void Fexp(const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void FexpShiftedScaled(const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void FS_values_and_derivatives(const int, KK_FLOAT&, const int) const;

  KOKKOS_INLINE_FUNCTION
  void evaluate_splines(const int, const int, KK_FLOAT, int, int, int, int) const;

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template<class TagStyle>
  void check_team_size_reduce(int, int&, int);

  template <typename scratch_type>
  KOKKOS_INLINE_FUNCTION
  int scratch_size_helper(int values_per_team) const;

  // Kokkos view typedefs
  typedef Kokkos::View<int*, DeviceType> t_ace_1i;
  typedef Kokkos::View<int**, DeviceType> t_ace_2i;
  typedef Kokkos::View<int**, Kokkos::LayoutRight, DeviceType> t_ace_2i_lr;
  typedef Kokkos::View<int***, DeviceType> t_ace_3i;
  typedef Kokkos::View<int***, Kokkos::LayoutRight, DeviceType> t_ace_3i_lr;
  typedef Kokkos::View<KK_FLOAT*, DeviceType> t_ace_1d;
  typedef Kokkos::View<KK_FLOAT**, DeviceType> t_ace_2d;
  typedef Kokkos::View<KK_FLOAT**, Kokkos::LayoutRight, DeviceType> t_ace_2d_lr;
  typedef Kokkos::View<KK_FLOAT*[3], DeviceType> t_ace_2d3;
  typedef Kokkos::View<KK_FLOAT***, DeviceType> t_ace_3d;
  typedef Kokkos::View<const KK_FLOAT***, DeviceType> tc_ace_3d;
  typedef Kokkos::View<KK_FLOAT**[3], DeviceType> t_ace_3d3;
  typedef Kokkos::View<KK_FLOAT**[4], DeviceType> t_ace_3d4;
  typedef Kokkos::View<KK_FLOAT**[4], Kokkos::LayoutRight, DeviceType> t_ace_3d4_lr;
  typedef Kokkos::View<KK_FLOAT****, DeviceType> t_ace_4d;
  typedef Kokkos::View<complex*, DeviceType> t_ace_1c;
  typedef Kokkos::View<complex**, DeviceType> t_ace_2c;
  typedef Kokkos::View<complex***, DeviceType> t_ace_3c;
  typedef Kokkos::View<complex**[3], DeviceType> t_ace_3c3;

  typedef typename Kokkos::View<KK_FLOAT*, DeviceType>::HostMirror th_ace_1d;

  // GRACE FS specific: A arrays (simplified, no element indexing)
  t_ace_3c A;          // [natom, (lmax+1)^2, nradmax]
  t_ace_3c A_sph;      // [natom, idx_sph_max, nradmax]

  t_ace_2c A_list;     // [natom, idx_ms_combs_max * rankmax]
  t_ace_2c A_forward_prod;

  t_ace_3c weights;     // [natom, idx_sph_max, nradmax]

  t_ace_1d e_atom;
  t_ace_2d rhos;
  t_ace_2d dF_drho;

  t_ace_2c dB_flatten;

  // radial functions
  t_ace_4d fr;
  t_ace_4d dfr;
  t_ace_3d gr;
  t_ace_3d dgr;
  t_ace_3d d_values;
  t_ace_3d d_derivatives;

  // Z coefficients for element-radial weighting
  t_ace_2d d_Z;  // [nelements, nradmax]

  // inverted active set for extrapolation grades
  tc_ace_3d d_ASI;
  t_ace_2d projections;
  t_ace_1d d_gamma;
  th_ace_1d h_gamma;

  // Spherical Harmonics
  void pre_compute_harmonics(int);

  t_ace_1d d_idx_sph;
  t_ace_1d alm;
  t_ace_1d blm;
  t_ace_1d cl;
  t_ace_1d dl;

  // short neigh list
  t_ace_1i d_ncount;
  t_ace_2d d_mu;
  t_ace_2d d_rnorms;
  t_ace_3d3 d_rhats;
  t_ace_2i d_nearest;

  // per-type
  t_ace_1i d_ndensity;
  t_ace_1i d_npoti;
  t_ace_1d d_E0vals;
  t_ace_2d_lr d_wpre;
  t_ace_2d_lr d_mexp;
  t_ace_1d d_shift;
  t_ace_1d d_scale_factor;

  // tilde - GRACE FS uses scalar ns per function
  t_ace_1i d_idx_ms_combs_count;
  t_ace_1i d_total_basis_size;
  t_ace_2i_lr d_rank;
  t_ace_2i_lr d_num_ms_combs;
  t_ace_2i_lr d_idx_funcs;
  t_ace_2i_lr d_ns;        // scalar ns per function (not per rank)
  t_ace_3i_lr d_ls;
  t_ace_3i_lr d_ms_combs;
  t_ace_2d d_gen_cgs;
  t_ace_3d d_coeffs;

  t_ace_3d3 f_ij;

  void deallocate_views_of_views();

 public:
  struct SplineInterpolatorKokkos {
    int ntot, nlut, num_of_functions;
    KK_FLOAT cutoff, deltaSplineBins, invrscalelookup, rscalelookup;

    t_ace_3d4_lr lookupTable;

    void operator=(const SplineInterpolator &spline);

    void deallocate() {
      lookupTable = t_ace_3d4_lr();
    }

    KK_FLOAT memory_usage() {
      return lookupTable.span() * sizeof(typename decltype(lookupTable)::value_type);
    }

    KOKKOS_INLINE_FUNCTION
    void calcSplines(const int ii, const int jj, const KK_FLOAT r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const;
  };

  Kokkos::DualView<SplineInterpolatorKokkos*, DeviceType> k_splines_gk;
  Kokkos::DualView<SplineInterpolatorKokkos*, DeviceType> k_splines_rnl;
  typename Kokkos::View<SplineInterpolatorKokkos*, DeviceType> d_splines_gk;
  typename Kokkos::View<SplineInterpolatorKokkos*, DeviceType> d_splines_rnl;

};
}    // namespace LAMMPS_NS

#endif
#endif
