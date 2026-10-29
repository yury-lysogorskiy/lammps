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
   Shared by grace/runtime and pace/runtime (CPU and KOKKOS): the pair_style keywords (pair_style takes the union;
   pair_coeff, once the model kind is known (grace::model_kind), rejects a keyword that does not apply to it; a keyword
   of the other backend that only tunes performance is ignored with a warning, so a CPU input runs unchanged with
   -sf kk), the model file of pair_coeff, the type map, the libgrace log callback and the rank count per node / GPU.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#ifndef LMP_GRACE_RUNTIME_SETTINGS_H
#define LMP_GRACE_RUNTIME_SETTINGS_H

#include "grace/model_kind.h"

#include <mpi.h>
#include <string>
#include <vector>

namespace LAMMPS_NS {
class LAMMPS;

// the bit of a model kind in the kind masks of the keyword table
constexpr unsigned grace_kind_bit(grace::ModelKind kind)
{
  return 1u << static_cast<int>(kind);
}

struct GraceRuntimeSettings {
  static constexpr grace::ModelKind kinds_all[] = {grace::ModelKind::Grace1L, grace::ModelKind::Grace2L,
                                                   grace::ModelKind::Grace3L, grace::ModelKind::Pace};
  enum : unsigned {    // model kinds
    K1L = grace_kind_bit(grace::ModelKind::Grace1L),
    K2L = grace_kind_bit(grace::ModelKind::Grace2L),
    K3L = grace_kind_bit(grace::ModelKind::Grace3L),
    KPACE = grace_kind_bit(grace::ModelKind::Pace),
    KALL = K1L | K2L | K3L | KPACE
  };
  enum : unsigned { CPU = 1, KK = 2 };    // backends
  enum Id { SPLINE, SPLINE_LAYERS, SPLINE_MB, STORE_FWD, FC_FOLD, STORE_L1, BLOCK, FUSE_RADIAL, NOSHAPE, CHUNKSIZE,
            DEBUG_NO_ENERGY_ONLY_CALC };
  struct Keyword {
    Id id;
    const char *name;
    int nvalues;
    unsigned kinds, backends;
    bool perf_only;    // tunes speed or memory only: ignored with a warning on the other backend
  };
  static constexpr Keyword keywords[] = {
      {SPLINE, "spline", 1, KALL, CPU | KK, false},
      {SPLINE_LAYERS, "spline_layers", 3, K3L, CPU | KK, false},
      {SPLINE_MB, "spline_mb", 1, K3L, CPU | KK, false},
      {STORE_FWD, "store_fwd", 1, K3L, CPU | KK, false},
      {FC_FOLD, "fc_fold", 1, K3L, CPU | KK, false},
      {STORE_L1, "store_l1", 1, K2L, CPU | KK, false},
      {BLOCK, "block", 1, K3L, CPU, true},
      {FUSE_RADIAL, "fuse_radial", 1, K3L, KK, true},
      {NOSHAPE, "noshape", 0, K1L | K2L, KK, true},
      {CHUNKSIZE, "chunksize", 1, KALL, KK, true},
      {DEBUG_NO_ENERGY_ONLY_CALC, "debug_no_energy_only_calc", 0, KALL, KK, false},
  };

  int spline = 1024;                   // intervals of the radial tables (3L: layer 3, layers 1-2 half); PACE: 0 none
  int spline_layer[3] = {0, 0, 0};     // 3L spline_layers N1 N2 N3, 0: from spline
  double spline_mb = -1.0;             // 3L radial table budget in MB, < 0: auto
  int store_fwd = -1;                  // 3L: 1 yes, 0 no, -1 auto (the CPU keeps them: yes)
  bool fc_fold = true;                 // 3L
  bool store_l1 = true;                // 2L
  int block = 16;                      // 3L CPU
  bool fuse_radial = true;             // 3L KK
  bool noshape = false;                // 1L / 2L KK
  int chunksize = -1;                  // KK: > 0 atoms per chunk, 0 one chunk, -1 auto
  bool debug_no_energy_only_calc = false;    // KK
  std::vector<const Keyword *> given;        // the keywords of pair_style that apply to this backend

  // pair_style: defaults, then the keywords (syntax and range); backend CPU or KK
  void parse(LAMMPS *lmp, const std::string &style, unsigned backend, int narg, char **arg);
  // pair_coeff: every keyword given must apply to the model kind
  void check_kind(LAMMPS *lmp, const std::string &style, grace::ModelKind kind, const std::string &path) const;
};

namespace GraceRuntime {

// the model file of pair_coeff * * model [model.asi] El ...: its kind, which must match the style (pace: PACE only,
// else GRACE only) and the keywords of pair_style; nskip: arguments before the elements beyond the model (.asi)
struct ModelFile {
  grace::ModelKind kind = grace::ModelKind::Grace1L;
  std::string path, asi;
  int nskip = 0;
};
ModelFile model_file(LAMMPS *lmp, const std::string &style, bool pace, const GraceRuntimeSettings &set, int narg,
                     char **arg);
// LAMMPS type -> model element (map[] indexes the deduplicated elements[]; -1: NULL of pair_style hybrid, unmapped)
std::vector<int> type2mu(LAMMPS *lmp, const std::string &style, const int *map, char *const *elements,
                         const std::vector<std::string> &names);
// libgrace log callback (user: the LAMMPS instance), rank 0: info lines to the log, warnings through Error::warning
void log(void *user, int level, const char *msg);
// the ranks of this node (shared-memory communicator) with the same key (e.g. the GPU's device key; "": every rank of
// the node); serial build (STUBS): 1
int ranks_sharing(MPI_Comm world, const std::string &key);

}    // namespace GraceRuntime

}    // namespace LAMMPS_NS

#endif
#endif    // NO_GRACE_RUNTIME
