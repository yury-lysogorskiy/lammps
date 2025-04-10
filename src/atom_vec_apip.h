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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(apip,AtomVecApip);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_APIP_H
#define LMP_ATOM_VEC_APIP_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecApip : public AtomVec {
 public:
  AtomVecApip(class LAMMPS *);

  void grow_pointers();

  void force_clear(int, size_t) override;
  void data_atom_post(int) override;

 protected:
  double *lambda, *lambda_input, *lambda_const, *lambda_input_ta, *e_simple, *e_complex,
      **f_const_lambda, **f_dyn_lambda;
  int *lambda_required;
};

#ifndef LMP_ATOM_LAMBDA_REQUIRED
#define LMP_ATOM_LAMBDA_REQUIRED
namespace LambdaRequired {

  enum {
    UNKNOWN = 0,
    SIMPLE = 1 << 0,
    NO_SIMPLE = 1 << 1,
    COMPLEX = 1 << 2,
    NO_COMPLEX = 1 << 3,
  };
}    // namespace LambdaRequired
#endif

}    // namespace LAMMPS_NS

#endif
#endif
