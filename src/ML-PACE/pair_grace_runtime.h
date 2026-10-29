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
   grace/runtime and pace/runtime (CPU, OpenMP + MPI): the libgrace CPU drivers (grace-runtime) of GRACE-1L / 2L / 3L
   (kokkos.npz) and PACE (pacemaker B-basis .yaml or its npz), picked at pair_coeff from the model file. The KOKKOS
   styles grace/runtime/kk and pace/runtime/kk are registered from KOKKOS/pair_grace_runtime_kokkos.h.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/runtime,PairGRACERuntime);
PairStyle(pace/runtime,PairPACERuntime);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_RUNTIME_H
#define LMP_PAIR_GRACE_RUNTIME_H

#include "grace_decouple.h"
#include "grace_runtime_settings.h"
#include "grace_uq.h"
#include "pair.h"

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace grace {
class Grace1LFused;
class Grace2LFused;
class Grace3LFused;
class PaceCpu;
}    // namespace grace

namespace LAMMPS_NS {

// glue over the libgrace CPU drivers: neighbour list, the layer-message comm of 2L / 3L, per-atom energy / virial,
// UQ for fix pair, pair_modify decouple/4d (GRACE), settings and pair_coeff
class PairGRACERuntime : public Pair {
 public:
  PairGRACERuntime(class LAMMPS *, bool pace = false);
  ~PairGRACERuntime() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void modify_params(int, char **) override;
  int modify_param(int, char **) override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  int comm_layer = 0;    // message layer of the pending forward / reverse comm (3L: 0 I1, 1 I2; 2L: 0)

 protected:
  const bool pace;            // pace/runtime: PACE models only; grace/runtime: GRACE-1L / 2L / 3L only
  const std::string name;     // style name for messages
  std::variant<std::monostate, std::unique_ptr<grace::Grace1LFused>, std::unique_ptr<grace::Grace2LFused>,
               std::unique_ptr<grace::Grace3LFused>, std::unique_ptr<grace::PaceCpu>>
      drv;
  GraceRuntimeSettings set;
  GraceRuntime::ModelFile mf;              // the model of the last pair_coeff
  std::vector<std::string> coeff_args;    // its arguments (a new pair_style replays them; empty before pair_coeff)
  bool isa_warned = false;                // the -march warning once
  int nthreads = 1;
  GraceUQOut uq;               // UQ (gamma, atomic_sigma, gmm_cluster; PACE gamma only) for fix pair
  GraceGlue::Decouple dc;      // pair_modify decouple/4d (GRACE): two-group coupling, lambda, dE/dlambda
  bool has_uq = false;         // the model's UQ head passed enable_uq() (PACE: an active set)

  void allocate();
  void build();
  template<class F> void with_driver(F &&);    // f(driver) for the driver of the model; nothing before pair_coeff
};

class PairPACERuntime : public PairGRACERuntime {
 public:
  PairPACERuntime(class LAMMPS *lmp) : PairGRACERuntime(lmp, true) {}
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif    // NO_GRACE_RUNTIME
