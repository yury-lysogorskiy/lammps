/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_apip.h"

#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecApip::AtomVecApip(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  double *lambda, *lambda_input, *lambda_input_ta, *lambda_const, *e_simple, *e_complex,
      **f_const_lambda, **f_dyn_lambda;
  int *lambda_required;
  atom->lambda_flag = 1;
  atom->lambda_input_flag = 1;
  atom->lambda_input_ta_flag = 1;
  atom->lambda_const_flag = 1;
  atom->lambda_required_flag = 1;
  atom->e_simple_flag = 1;
  atom->e_complex_flag = 1;
  atom->f_const_lambda_flag = 1;
  atom->f_dyn_lambda_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  // The full list of fields is in atom_vec.cpp
  fields_copy = {"lambda", "lambda_required", "lambda_input", "lambda_input_ta", "lambda_const"};
  fields_comm = {"lambda", "lambda_required", "lambda_input_ta", "lambda_const"};
  fields_comm_vel = {};
  fields_border = {"lambda", "lambda_required", "lambda_input_ta", "lambda_const"};
  fields_border_vel = {};
  fields_exchange = {"lambda", "lambda_required", "lambda_input_ta", "lambda_const"};
  fields_restart = {"lambda", "lambda_required", "lambda_input", "lambda_input_ta", "lambda_const"};
  fields_create = {};
  fields_grow = {
      "lambda",   "lambda_required", "lambda_input",   "lambda_input_ta", "lambda_const",
      "e_simple", "e_complex",       "f_const_lambda", "f_dyn_lambda"};    // allocates memory
  fields_reverse = {"f_const_lambda",
                    "f_dyn_lambda"};    // communication of force after calculation
  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecApip::grow_pointers()
{
  lambda = atom->lambda;
  lambda_required = atom->lambda_required;
  lambda_input = atom->lambda_input;
  lambda_input_ta = atom->lambda_input_ta;
  lambda_const = atom->lambda_const;
  e_simple = atom->e_simple;
  e_complex = atom->e_complex;
  f_const_lambda = atom->f_const_lambda;
  f_dyn_lambda = atom->f_dyn_lambda;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecApip::data_atom_post(int ilocal)
{
  lambda[ilocal] = 0;
  lambda_const[ilocal] = 0;
  lambda_required[ilocal] = LambdaRequired::UNKNOWN;
  lambda_input[ilocal] = 0;
  lambda_input_ta[ilocal] = 0;
  e_simple[ilocal] = 0;
  e_complex[ilocal] = 0;
  f_const_lambda[ilocal][0] = 0;
  f_const_lambda[ilocal][1] = 0;
  f_const_lambda[ilocal][2] = 0;
  f_dyn_lambda[ilocal][0] = 0;
  f_dyn_lambda[ilocal][1] = 0;
  f_dyn_lambda[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom n
   natoms = # of atoms to clear
   nbytes = natoms * sizeof(double)
   requires forceclearflag = 1 to be called
------------------------------------------------------------------------- */

void AtomVecApip::force_clear(int n, size_t nbytes)
{
  memset(&f_const_lambda[n][0], 0, 3 * nbytes);
  memset(&f_dyn_lambda[n][0], 0, 3 * nbytes);
  memset(&lambda_required[n], 0, nbytes / sizeof(double) * sizeof(int));
}
