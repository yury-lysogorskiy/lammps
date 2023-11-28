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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(set2,Set2);
// clang-format on
#else

#ifndef LMP_SET2_H
#define LMP_SET2_H

#include "command.h"

namespace LAMMPS_NS {

class Set2 : public Command {
 public:
  Set2(class LAMMPS *lmp) : Command(lmp){};
  ~Set2();
  
  void command(int, char **) override;

  void process_args(int, int, char **);
  void selection(int);
  void invoke_actions();

 private: 
  int caller;           // SETCOMMAND or FIXSET

  // params for atom selection

  int style;
  char *id;
  int nlo,nhi;
  bigint nlobig,nhibig;
  int groupbit;
  class Region *region;

  // atom property setting params for keyword/value pairs

  int ivalue, newtype, count, index_custom, icol_custom;
  int ximage, yimage, zimage, ximageflag, yimageflag, zimageflag;
  int cc_index;
  bigint nsubset;
  double dvalue, xvalue, yvalue, zvalue, wvalue, fraction;
  int discflag;

  // one Action = one keyword/value pair
  
  struct Action {
    int keyword;
    int varflag;
    int varflag1, varflag2, varflag3, varflag4;
    int ivar1, ivar2, ivar3, ivar4;
    int ivalue1;
    tagint tvalue1;
    double dvalue1,dvalue2,dvalue3,dvalue4;
  };
  
  int naction,maxaction;
  Action *actions;

  typedef void (Set2::*FnPtrPack)();
  FnPtrPack *invoke_choice;    // list of ptrs to invoke functions

  // storage for evaluated variables
  
  double *vec1, *vec2, *vec3, *vec4;

  // flag for selected atoms
  
  int *select;

  // private functions
  
  void setrandom(int);
  void topology(int);
  void varparse(const char *, int);
 
  // customize by adding a process method prototype

  void process_angle(int &, int, char **);
  void process_angmom(int &, int, char **);
  void process_bond(int &, int, char **);
  void process_cc(int &, int, char **);
  void process_charge(int &, int, char **);
  void process_density(int &, int, char **);
  void process_diameter(int &, int, char **);
  void process_dihedral(int &, int, char **);
  void process_dipole(int &, int, char **);
  void process_dipole_random(int &, int, char **);
  void process_dpd_theta(int &, int, char **);
  void process_edpd_cv(int &, int, char **);
  void process_edpd_temp(int &, int, char **);
  void process_epsilon(int &, int, char **);
  void process_image(int &, int, char **);
  void process_improper(int &, int, char **);
  void process_length(int &, int, char **);
  void process_mass(int &, int, char **);
  void process_mol(int &, int, char **);
  void process_omega(int &, int, char **);
  void process_quat(int &, int, char **);
  void process_quat_random(int &, int, char **);
  void process_radius_election(int &, int, char **);
  void process_shape(int &, int, char **);
  void process_smd_contact_radius(int &, int, char **);
  void process_smd_mass_density(int &, int, char **);
  void process_sph_cv(int &, int, char **);
  void process_sph_e(int &, int, char **);
  void process_sph_rho(int &, int, char **);
  void process_spin_atom(int &, int, char **);
  void process_spin_atom_random(int &, int, char **);
  void process_spin_electron(int &, int, char **);
  void process_temperature(int &, int, char **);
  void process_theta(int &, int, char **);
  void process_theta_random(int &, int, char **);
  void process_tri(int &, int, char **);
  void process_type(int &, int, char **);
  void process_type_fraction(int &, int, char **);
  void process_type_ratio(int &, int, char **);
  void process_type_subset(int &, int, char **);
  void process_volume(int &, int, char **);
  void process_vx(int &, int, char **);
  void process_vy(int &, int, char **);
  void process_vz(int &, int, char **);
  void process_x(int &, int, char **);
  void process_y(int &, int, char **);
  void process_z(int &, int, char **);

  void process_custom(int &, int, char **);

  // customize by adding an invoke method prototype

  void invoke_angle();
  void invoke_angmom();
  void invoke_bond();
  void invoke_cc();
  void invoke_charge();
  void invoke_density();
  void invoke_diameter();
  void invoke_dihedral();
  void invoke_dipole();
  void invoke_dipole_random();
  void invoke_dpd_theta();
  void invoke_edpd_cv();
  void invoke_edpd_temp();
  void invoke_epsilon();
  void invoke_image();
  void invoke_improper();
  void invoke_length();
  void invoke_mass();
  void invoke_mol();
  void invoke_omega();
  void invoke_quat();
  void invoke_quat_random();
  void invoke_radius_election();
  void invoke_shape();
  void invoke_smd_contact_radius();
  void invoke_smd_mass_density();
  void invoke_sph_cv();
  void invoke_sph_e();
  void invoke_sph_rho();
  void invoke_spin_atom();
  void invoke_spin_atom_random();
  void invoke_spin_electron();
  void invoke_temperature();
  void invoke_theta();
  void invoke_theta_random();
  void invoke_tri();
  void invoke_type();
  void invoke_type_fraction();
  void invoke_type_ratio();
  void invoke_type_subset();
  void invoke_volume();
  void invoke_vx();
  void invoke_vy();
  void invoke_vz();
  void invoke_x();
  void invoke_y();
  void invoke_z();

  void invoke_custom();
};

}    // namespace LAMMPS_NS

#endif
#endif

