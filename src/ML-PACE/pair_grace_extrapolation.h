/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_TF
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/extrapolation,PairGRACEExtrapolation);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_EXTRAPOLATION_H
#define LMP_PAIR_GRACE_EXTRAPOLATION_H

#include "pair_grace.h"

namespace LAMMPS_NS {

class PairGRACEExtrapolation : public PairGRACE {
 public:
  PairGRACEExtrapolation(class LAMMPS *);

  void settings(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif    //#ifndef NO_GRACE_TF
