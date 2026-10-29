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
   grace/runtime/kk and pace/runtime/kk: LAMMPS glue over the libgrace Kokkos engines (grace::kk, grace-runtime) of
   GRACE-1L / 2L / 3L and PACE, picked at pair_coeff from the model file. GRACE: bond geometry in fp64, the network
   and the radial tables in fp32, energies, forces and virials accumulated in fp64; PACE in fp64. The engine holds the
   kernels; this class only moves LAMMPS data in and out.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/runtime/kk,        PairGRACERuntimeKokkos_Device);
PairStyle(grace/runtime/kk/device, PairGRACERuntimeKokkos_Device);
PairStyle(pace/runtime/kk,         PairPACERuntimeKokkos_Device);
PairStyle(pace/runtime/kk/device,  PairPACERuntimeKokkos_Device);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRACE_RUNTIME_KOKKOS_H
#define LMP_PAIR_GRACE_RUNTIME_KOKKOS_H

#include "pair.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "grace_decouple.h"
#include "grace_runtime_settings.h"
#include "grace_uq.h"

#include "grace/kk/grace_1l_fused_kk.h"
#include "grace/kk/grace_2l_fused_kk.h"
#include "grace/kk/grace_3l_fused_kk.h"
#include "grace/kk/pace_kk.h"

#include <string>
#include <variant>
#include <vector>

namespace LAMMPS_NS {

// the LAMMPS view types for the libgrace engines
template<class DeviceType> struct GraceLammpsTraits {
  using AT = ArrayTypes<DeviceType>;
  using kk_float = KK_FLOAT;
  using kk_acc_float = KK_ACC_FLOAT;
  using ev_type = EV_FLOAT;
  using kk_device = typename KKDevice<DeviceType>::value;
  using layout_f = typename DAT::t_kkacc_1d_3::array_layout;
  using layout_vatom = typename DAT::t_kkacc_1d_6::array_layout;
  using layout_cvatom = typename DAT::t_kkacc_1d_9::array_layout;
  using t_kkfloat_1d_3_lr_randomread = typename AT::t_kkfloat_1d_3_lr_randomread;
  using t_kkacc_1d_3 = typename AT::t_kkacc_1d_3;
  using t_kkacc_1d = typename AT::t_kkacc_1d;
  using t_kkacc_1d_6 = typename AT::t_kkacc_1d_6;
  using t_kkacc_1d_9 = typename AT::t_kkacc_1d_9;
  using t_int_1d = typename AT::t_int_1d;
  using t_double_1d = typename AT::t_double_1d;
  using t_int_1d_randomread = typename AT::t_int_1d_randomread;
  using t_neighbors_2d = typename AT::t_neighbors_2d;
  static constexpr int neighmask = NEIGHMASK;
};

// one engine per model kind, held by value; the layer messages of 2L (one) and 3L (two, comm_layer) go to the ghost
// atoms through LAMMPS comm (forward_comm / reverse_comm of this pair style), on the device (KokkosBase) or, with the
// legacy pair comm or a non-Kokkos comm, through the host
template<class DeviceType> class PairGRACERuntimeKokkos : public Pair, public KokkosBase {
 public:
  using Traits = GraceLammpsTraits<DeviceType>;
  using E1 = grace::kk::Grace1LFusedKK<DeviceType, Traits>;
  using E2 = grace::kk::Grace2LFusedKK<DeviceType, Traits>;
  using E3 = grace::kk::Grace3LFusedKK<DeviceType, Traits>;
  using EP = grace::kk::PaceKK<DeviceType, Traits>;

  PairGRACERuntimeKokkos(class LAMMPS *, bool pace = false);
  ~PairGRACERuntimeKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void modify_params(int, char **) override;
  int modify_param(int, char **) override;
  double init_one(int, int) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &) override;

 protected:
  const bool pace;           // pace/runtime/kk: PACE models only; grace/runtime/kk: GRACE-1L / 2L / 3L only
  const std::string name;    // style name for messages
  std::variant<std::monostate, E1, E2, E3, EP> eng;
  GraceRuntimeSettings set;
  GraceRuntime::ModelFile mf;              // the model of the last pair_coeff
  std::vector<std::string> coeff_args;    // its arguments (a new pair_style replays them; empty before pair_coeff)
  std::vector<int> type2mu;    // LAMMPS type -> model element (-1 unmapped)
  bool has_uq = false;
  int ranks_per_device = 0;    // processes on this rank's GPU (chunksize auto), counted at the first init_style
  int comm_layer = 0;          // 3L: layer of the exchange in progress
  GraceUQOut uq;
  GraceGlue::Decouple dc;      // pair_modify decouple/4d (GRACE): two-group coupling, lambda, dE/dlambda
  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  DAT::ttransform_kkacc_1d_9 k_cvatom;

  struct Exchange2 : E2::Exchange {
    PairGRACERuntimeKokkos *p = nullptr;
    void forward() override;
    void reverse() override;
  } exchange2;
  struct Exchange3 : E3::Exchange3 {
    PairGRACERuntimeKokkos *p = nullptr;
    void forward(int layer) override;
    void reverse(int layer) override;
  } exchange3;

  void allocate();
  void build();
  template<class F> void with_engine(F &&);    // f(engine) for the engine of the model; nothing before pair_coeff
  template<class E> void load(const std::string &path, const std::string &asi);
  template<class E> void run(E &);
  void prepare_peratom();
  void finish_peratom();
  // one exchange of the layer messages through LAMMPS comm (forward or reverse); to_host / from_host move the rows
  // between device and host around a host pack / unpack
  template<class ToHost, class FromHost> void pair_exchange(bool reverse, ToHost &&to_host, FromHost &&from_host);
};

template<class DeviceType> class PairPACERuntimeKokkos : public PairGRACERuntimeKokkos<DeviceType> {
 public:
  PairPACERuntimeKokkos(class LAMMPS *lmp) : PairGRACERuntimeKokkos<DeviceType>(lmp, true) {}
};

using PairGRACERuntimeKokkos_Device = PairGRACERuntimeKokkos<LMPDeviceType>;
using PairPACERuntimeKokkos_Device = PairPACERuntimeKokkos<LMPDeviceType>;

}    // namespace LAMMPS_NS

// the GRACE engines are instantiated in grace_runtime_kokkos_{1l,2l,3l}.cpp, compiled in parallel with the style (PACE:
// header-only, with the style)
#ifndef GRACE_KK_EXTERN_1L
#error "grace-runtime too old for this LAMMPS (no GRACE_KK_EXTERN_*): use the pinned revision"
#endif
GRACE_KK_EXTERN_1L(LMPDeviceType, LAMMPS_NS::GraceLammpsTraits<LMPDeviceType>);
GRACE_KK_EXTERN_2L(LMPDeviceType, LAMMPS_NS::GraceLammpsTraits<LMPDeviceType>);
GRACE_KK_EXTERN_3L(LMPDeviceType, LAMMPS_NS::GraceLammpsTraits<LMPDeviceType>);

#endif
#endif
#endif    // NO_GRACE_RUNTIME
