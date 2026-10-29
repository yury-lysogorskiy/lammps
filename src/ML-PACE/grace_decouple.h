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

// 4D decoupling (Pablo Pena Cano) of the libgrace-backed grace/runtime(/kk) styles:
//   pair_modify decouple/4d <group1> <group2> [eq]     eq 1..4 (default 4), see grace-runtime grace/coupling.h
//   pair_modify decouple/4d off
// lambda is a global scalar of the style (extract("lambda"), set by fix adapt, also inside rerun), dE/dlambda the
// first extra quantity (compute pair ... -> c_ID[1], summed over the ranks). Lambda 0: the bonds are unchanged.

#ifndef LMP_GRACE_DECOUPLE_H
#define LMP_GRACE_DECOUPLE_H

#include "atom.h"
#include "error.h"
#include "group.h"
#include "lammps.h"
#include "pair.h"
#include "utils.h"

#include <cstring>
#include <mpi.h>
#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace GraceGlue {

struct Decouple {
  int eq = 0;                // shift equation 1..4, 0: off
  double lambda = 0.0;       // fix adapt target
  std::string group1, group2;
  int bit1 = 0, bit2 = 0;    // group bitmasks, from init()
  double pvector[1] = {0.0};    // the pair style's extra quantity: dE/dlambda of this rank's atoms

  // one pair_modify keyword at arg[0]: decouple/4d with its arguments applied, returns the number of args consumed
  // (0: not decouple/4d). Pair::modify_param of the styles, so the keyword also works in pair_modify pair <style> ...
  int modify_param(LAMMPS *lmp, const char *style, int narg, char **arg)
  {
    if (narg < 1 || strcmp(arg[0], "decouple/4d") != 0) return 0;
    if (narg < 2) utils::missing_cmd_args(FLERR, "pair_modify decouple/4d", lmp->error);
    if (strcmp(arg[1], "off") == 0) {
      eq = 0;
      return 2;
    }
    if (narg < 3) utils::missing_cmd_args(FLERR, "pair_modify decouple/4d", lmp->error);
    group1 = arg[1];
    group2 = arg[2];
    if (group1 == group2) lmp->error->all(FLERR, "pair {}: decouple/4d needs two different groups", style);
    eq = 4;
    if (narg > 3 && utils::is_integer(arg[3])) {
      eq = utils::inumeric(FLERR, arg[3], false, lmp);
      if (eq < 1 || eq > 4) lmp->error->all(FLERR, "pair {}: decouple/4d eq must be 1..4, not {}", style, eq);
      return 4;
    }
    return 3;
  }

  // the style's pair_modify: its own keywords here, the others to Pair::modify_params (as PairGRACE::modify_params)
  static void modify_params(Pair *pair, int narg, char **arg)
  {
    std::vector<char *> rest;
    for (int iarg = 0; iarg < narg;) {
      const int n = pair->modify_param(narg - iarg, &arg[iarg]);
      if (n > 0) iarg += n;
      else rest.push_back(arg[iarg++]);
    }
    if (!rest.empty()) pair->Pair::modify_params((int) rest.size(), rest.data());
  }

  // init_style: group bitmasks (the groups may be defined after pair_modify); no atom may be in both groups
  void init(LAMMPS *lmp, const char *style)
  {
    if (!eq) return;
    const int g1 = lmp->group->find(group1), g2 = lmp->group->find(group2);
    if (g1 < 0) lmp->error->all(FLERR, "pair {}: decouple/4d group {} does not exist", style, group1);
    if (g2 < 0) lmp->error->all(FLERR, "pair {}: decouple/4d group {} does not exist", style, group2);
    bit1 = lmp->group->bitmask[g1];
    bit2 = lmp->group->bitmask[g2];
    const int *mask = lmp->atom->mask;
    int both = 0, sum = 0;
    for (int i = 0; i < lmp->atom->nlocal; i++)
      if ((mask[i] & bit1) && (mask[i] & bit2)) both++;
    MPI_Allreduce(&both, &sum, 1, MPI_INT, MPI_SUM, lmp->world);
    if (sum)
      lmp->error->all(FLERR, "pair {}: decouple/4d: {} atoms are in both groups {} and {}", style, sum, group1,
                      group2);
  }

  void *extract(const char *str, int &dim)
  {
    if (strcmp(str, "lambda") != 0) return nullptr;
    dim = 0;
    return &lambda;
  }
};

}    // namespace GraceGlue
}    // namespace LAMMPS_NS

#endif
