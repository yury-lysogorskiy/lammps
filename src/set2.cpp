// clang-format off
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

#include "set2.h"

#include "arg_info.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "random_park.h"
#include "region.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{SET,FIXSET};

enum{ATOM_SELECT,MOL_SELECT,TYPE_SELECT,GROUP_SELECT,REGION_SELECT};

enum{TYPE,TYPE_FRACTION,TYPE_RATIO,TYPE_SUBSET,
     MOLECULE,X,Y,Z,VX,VY,VZ,CHARGE,MASS,SHAPE,LENGTH,TRI,
     DIPOLE,DIPOLE_RANDOM,SPIN_ATOM,SPIN_RANDOM,SPIN_ELECTRON,RADIUS_ELECTRON,
     QUAT,QUAT_RANDOM,THETA,THETA_RANDOM,ANGMOM,OMEGA,TEMPERATURE,
     DIAMETER,DENSITY,VOLUME,IMAGE,BOND,ANGLE,DIHEDRAL,IMPROPER,
     SPH_E,SPH_CV,SPH_RHO,EDPD_TEMP,EDPD_CV,CC,SMD_MASS_DENSITY,
     SMD_CONTACT_RADIUS,DPDTHETA,EPSILON,IVEC,DVEC,IARRAY,DARRAY};

#define BIG INT_MAX
#define DELTA 4

/* ---------------------------------------------------------------------- */

void Set2::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Set command on system without atoms");

  select = nullptr;

  process_args(SET,narg-3,&arg[3]);

  if (comm->me == 0) utils::logmesg(lmp,"Setting atom values ...\n");

  selection(atom->nlocal);
  
  // loop over list of actions to reset attributes

  invoke_actions();

  // statistics
  // for CC option, include species index

  int origarg;
  
  int count,allcount;
  MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
  
  if (comm->me == 0) {
    if (strcmp(arg[origarg],"cc") == 0)
      utils::logmesg(lmp,"  {} settings made for {} index {}\n",
                     allcount,arg[origarg],arg[origarg+1]);
    else
      utils::logmesg(lmp,"  {} settings made for {}\n",
                     allcount,arg[origarg]);
  }

  // clean up

  memory->destroy(select);
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom coordinates
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set2::process_args(int caller_flag, int narg, char **arg)
{
  caller = caller_flag;
  
  if (narg < 3) error->all(FLERR,"Illegal set command");

  // style and ID info

  id = utils::strdup(arg[1]);

  if (strcmp(arg[0],"atom") == 0) {
    style = ATOM_SELECT;
    if (atom->tag_enable == 0)
      error->all(FLERR,"Cannot use set atom with no atom IDs defined");
    utils::bounds(FLERR,id,1,MAXTAGINT,nlobig,nhibig,error);

  } else if (strcmp(arg[0],"mol") == 0) {
    style = MOL_SELECT;
    if (atom->molecule_flag == 0)
      error->all(FLERR,"Cannot use set mol with no molecule IDs defined");
    utils::bounds(FLERR,id,1,MAXTAGINT,nlobig,nhibig,error);
    
  } else if (strcmp(arg[0],"type") == 0) {
    style = TYPE_SELECT;
    if (char *typestr = utils::expand_type(FLERR,id,Atom::ATOM,lmp)) {
      delete [] id;
      id = typestr;
    }
    utils::bounds(FLERR,id,1,atom->ntypes,nlo,nhi,error);

  } else if (strcmp(arg[0],"group") == 0) {
    style = GROUP_SELECT;
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR,"Could not find set group ID {}", id);
    groupbit = group->bitmask[igroup];

  } else if (strcmp(arg[0],"region") == 0) {
    style = REGION_SELECT;
    region = domain->get_region_by_id(id);
    if (!region) error->all(FLERR,"Set region {} does not exist", id);

  } else error->all(FLERR,"Unknown set command style: {}", arg[0]);

  delete [] id;
  
  // loop over keyword/value pairs to create list of actions
  // one action = keyword/value pair

  naction = maxaction = 0;
  actions = nullptr;
  
  int iarg = 2;
  while (iarg < narg) {
    if (naction == maxaction) {
      maxaction += DELTA;
      actions = (Action *) memory->srealloc(actions,maxaction*sizeof(Action),"set:actions");
      invoke_choice = (FnPtrPack *)
        memory->srealloc(invoke_choice,maxaction*sizeof(FnPtrPack),"set:invoke_choice");
    }

    if (strcmp(arg[iarg],"angle") == 0) {
      process_angle(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_angle;
    } else if (strcmp(arg[iarg],"angmom") == 0) {
      process_angmom(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_angmom;
    } else if (strcmp(arg[iarg],"bond") == 0) {
      process_bond(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_bond;
    } else if (strcmp(arg[iarg],"cc") == 0) {
      process_cc(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_cc;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      process_charge(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_charge;
    } else if (strcmp(arg[iarg],"density") == 0 ||(strcmp(arg[iarg],"density/disc") == 0)) {
      process_density(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_density;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      process_diameter(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_diameter;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      process_dihedral(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_dihedral;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      process_dipole(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_dipole;
    } else if (strcmp(arg[iarg],"dipole/random") == 0) {
      process_dipole_random(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_dipole_random;
    } else if (strcmp(arg[iarg],"dpd/theta") == 0) {
      process_dpd_theta(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_dpd_theta;
    } else if (strcmp(arg[iarg],"edpd/cv") == 0) {
      process_edpd_cv(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_edpd_cv;
    } else if (strcmp(arg[iarg],"edpd/temp") == 0) {
      process_edpd_temp(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_edpd_temp;
    } else if (strcmp(arg[iarg],"epsilon") == 0) {
      process_epsilon(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_epsilon;
    } else if (strcmp(arg[iarg],"image") == 0) {
      process_image(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_image;
    } else if (strcmp(arg[iarg],"improper") == 0) {
      process_improper(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_improper;
    } else if (strcmp(arg[iarg],"length") == 0) {
      process_length(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_length;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      process_mass(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_mass;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      process_mol(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_mol;
    } else if (strcmp(arg[iarg],"omega") == 0) {
      process_omega(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_omega;
    } else if (strcmp(arg[iarg],"quat") == 0) {
      process_quat(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_quat;
    } else if (strcmp(arg[iarg],"quat/random") == 0) {
      process_quat_random(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_quat_random;
    } else if (strcmp(arg[iarg],"radius/electron") == 0) {
      process_radius_election(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_radius_election;
    } else if (strcmp(arg[iarg],"shape") == 0) {
      process_shape(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_shape;
    } else if (strcmp(arg[iarg],"smd/contact/radius") == 0) {
      process_smd_contact_radius(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_smd_contact_radius;
    } else if (strcmp(arg[iarg],"smd/mass/density") == 0) {
      process_smd_mass_density(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_density;
    } else if (strcmp(arg[iarg],"sph/cv") == 0) {
      process_sph_cv(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_sph_cv;
    } else if (strcmp(arg[iarg],"sph/e") == 0) {
      process_sph_e(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_sph_e;
    } else if (strcmp(arg[iarg],"sph/rho") == 0) {
      process_sph_rho(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_sph_rho;
    } else if ((strcmp(arg[iarg],"spin/atom") == 0) || (strcmp(arg[iarg],"spin") == 0)) {
      process_spin_atom(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_spin_atom;
    } else if ((strcmp(arg[iarg],"spin/atom/random") == 0) || (strcmp(arg[iarg],"spin/random") == 0)) {
      process_spin_atom_random(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_spin_atom_random;
    } else if (strcmp(arg[iarg],"spin/electron") == 0) {
      process_spin_electron(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_spin_electron;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      process_temperature(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_temperature;
    } else if (strcmp(arg[iarg],"theta") == 0) {
      process_theta(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_theta;
    } else if (strcmp(arg[iarg],"theta/random") == 0) {
      process_theta_random(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_theta_random;
    } else if (strcmp(arg[iarg],"tri") == 0) {
      process_tri(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_tri;
    } else if (strcmp(arg[iarg],"type") == 0) {
      process_type(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_type;
    } else if (strcmp(arg[iarg],"type/fraction") == 0) {
      process_type_fraction(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_type_fraction;
    } else if (strcmp(arg[iarg],"type/ratio") == 0) {
      process_type_ratio(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_type_ratio;
    } else if (strcmp(arg[iarg],"type/subset") == 0) {
      process_type_subset(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_type_subset;
    } else if (strcmp(arg[iarg],"volume") == 0) {
      process_volume(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_volume;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      process_vx(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_vx;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      process_vy(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_vy;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      process_vz(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_vz;
    } else if (strcmp(arg[iarg],"x") == 0) {
      process_x(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_x;
    } else if (strcmp(arg[iarg],"y") == 0) {
      process_y(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_y;
    } else if (strcmp(arg[iarg],"z") == 0) {
      process_z(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_z;
      
    } else if (utils::strmatch(arg[iarg],"^[id]2?_")) {
      process_custom(iarg,narg,arg);
      invoke_choice[naction++] = &Set2::invoke_custom;

    } else {
      error->all(FLERR,"Illegal set command keyword");
    }
  }
}

/* ----------------------------------------------------------------------
   loop over list of actions
   perform each on all selected atoms via call to invoke_choice() method
------------------------------------------------------------------------- */

void Set2::invoke_actions()
{
  for (int i = 0; i < naction; i++) {

    // overwrite dvalue, ivalue, xyzw value if variables defined
    // else the input script scalar value remains in place

    if (varflag) {
      if (varflag1) {
        dvalue = xvalue = vec1[i];
        ivalue = static_cast<int> (dvalue);
      }
      if (varflag2) yvalue = vec2[i];
      if (varflag3) zvalue = vec3[i];
      if (varflag4) wvalue = vec4[i];
    }
    
    (this->*invoke_choice[i])();
  }
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, MOLECULE, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set2::selection(int n)
{
  memory->destroy(select);
  memory->create(select,n,"set:select");

  if (style == ATOM_SELECT) {
    tagint *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == MOL_SELECT) {
    tagint *molecule = atom->molecule;
    for (int i = 0; i < n; i++)
      if (molecule[i] >= nlobig && molecule[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    region->prematch();
    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (region->match(x[i][0],x[i][1],x[i][2])) select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set owned atom properties directly
   either scalar or per-atom values from atom-style variable(s)
------------------------------------------------------------------------- */

void Set2::set(int keyword)
{
  // evaluate atom-style variable(s) if necessary

  vec1 = vec2 = vec3 = vec4 = nullptr;

  if (varflag) {
    int nlocal = atom->nlocal;
    if (varflag1) {
      memory->create(vec1,nlocal,"set:vec1");
      input->variable->compute_atom(ivar1,0,vec1,1,0);
    }
    if (varflag2) {
      memory->create(vec2,nlocal,"set:vec2");
      input->variable->compute_atom(ivar2,0,vec2,1,0);
    }
    if (varflag3) {
      memory->create(vec3,nlocal,"set:vec3");
      input->variable->compute_atom(ivar3,0,vec3,1,0);
    }
    if (varflag4) {
      memory->create(vec4,nlocal,"set:vec4");
      input->variable->compute_atom(ivar4,0,vec4,1,0);
    }
  }

  // check if properties of atoms in rigid bodies are updated
  // that are cached as per-body data
  
  switch (keyword) {
  case X:
  case Y:
  case Z:
  case MOLECULE:
  case MASS:
  case ANGMOM:
  case SHAPE:
  case DIAMETER:
  case DENSITY:
  case TEMPERATURE:
  case QUAT:
  case IMAGE:
    if (modify->check_rigid_list_overlap(select))
      error->warning(FLERR,"Changing a property of atoms in rigid bodies "
                     "that has no effect unless rigid bodies are re-initialized");
    break;
  default: // assume no conflict for all other properties
    break;
  }

  // clear up per-atom memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
  memory->destroy(vec4);
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom coordinates
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set2::setrandom(int keyword)
{
  int i;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));

  double **x = atom->x;
  int seed = ivalue;

  auto ranpark = new RanPark(lmp,1);
  auto ranmars = new RanMars(lmp,seed + comm->me);

  // set approx fraction of atom types to newtype

  if (keyword == TYPE_FRACTION) {
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (select[i]) {
        ranpark->reset(seed,x[i]);
        if (ranpark->uniform() > fraction) continue;
        atom->type[i] = newtype;
        count++;
      }

  // set exact count of atom types to newtype
  // for TYPE_RATIO, exact = fraction out of total eligible
  // for TYPE_SUBSET, exact = nsubset out of total eligible

  } else if (keyword == TYPE_RATIO || keyword == TYPE_SUBSET) {
    int nlocal = atom->nlocal;

    // count = number of eligible atoms I own

    count = 0;
    for (i = 0; i < nlocal; i++)
      if (select[i]) count++;

    // convert specified fraction to nsubset

    bigint bcount = count;
    bigint allcount;
    MPI_Allreduce(&bcount,&allcount,1,MPI_LMP_BIGINT,MPI_SUM,world);

    if (keyword == TYPE_RATIO) {
      nsubset = static_cast<bigint> (fraction * allcount);
    } else if (keyword == TYPE_SUBSET) {
      if (nsubset > allcount)
        error->all(FLERR,"Set type/subset value exceeds eligible atoms");
    }

    // make selection

    int *flag = memory->create(flag,count,"set:flag");
    int *work = memory->create(work,count,"set:work");

    ranmars->select_subset(nsubset,count,flag,work);

    // change types of selected atoms
    // flag vector from select_subset() is only for eligible atoms

    count = 0;
    int eligible = 0;
    for (i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (flag[eligible]) {
        atom->type[i] = newtype;
        count++;
      }
      eligible++;
    }

    // clean up

    memory->destroy(flag);
    memory->destroy(work);

  // set dipole moments to random orientations in 3d or 2d
  // dipole length is determined by dipole type array

  } else if (keyword == DIPOLE_RANDOM) {
    double **mu = atom->mu;
    int nlocal = atom->nlocal;

    double msq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          mu[i][0] = ranpark->uniform() - 0.5;
          mu[i][1] = ranpark->uniform() - 0.5;
          mu[i][2] = ranpark->uniform() - 0.5;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][2] *= scale;
          mu[i][3] = dvalue;
          count++;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          mu[i][0] = ranpark->uniform() - 0.5;
          mu[i][1] = ranpark->uniform() - 0.5;
          mu[i][2] = 0.0;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][3] = dvalue;
          count++;
        }
    }


  // set spin moments to random orientations in 3d or 2d
  // spin length is fixed to unity

  } else if (keyword == SPIN_RANDOM) {
    double **sp = atom->sp;
    int nlocal = atom->nlocal;

    double sp_sq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          sp[i][0] = ranpark->uniform() - 0.5;
          sp[i][1] = ranpark->uniform() - 0.5;
          sp[i][2] = ranpark->uniform() - 0.5;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1] + sp[i][2]*sp[i][2];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][2] *= scale;
          sp[i][3] = dvalue;
          count++;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          ranpark->reset(seed,x[i]);
          sp[i][0] = ranpark->uniform() - 0.5;
          sp[i][1] = ranpark->uniform() - 0.5;
          sp[i][2] = 0.0;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][3] = dvalue;
          count++;
        }
    }

  // set quaternions to random orientations in 3d and 2d

  } else if (keyword == QUAT_RANDOM) {
    int nlocal = atom->nlocal;
    double *quat;
    double **quat2;

    if (domain->dimension == 3) {
      double s,t1,t2,theta1,theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          else if (avec_tri && atom->tri[i] >= 0)
            quat = avec_tri->bonus[atom->tri[i]].quat;
          else if (avec_body && atom->body[i] >= 0)
            quat = avec_body->bonus[atom->body[i]].quat;
          else if (atom->quat_flag)
            quat2 = atom->quat;
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          ranpark->reset(seed,x[i]);
          s = ranpark->uniform();
          t1 = sqrt(1.0-s);
          t2 = sqrt(s);
          theta1 = 2.0*MY_PI*ranpark->uniform();
          theta2 = 2.0*MY_PI*ranpark->uniform();
          if (atom->quat_flag) {
            quat2[i][0] = cos(theta2)*t2;
            quat2[i][1] = sin(theta1)*t1;
            quat2[i][2] = cos(theta1)*t1;
            quat2[i][3] = sin(theta2)*t2;
          } else {
            quat[0] = cos(theta2)*t2;
            quat[1] = sin(theta1)*t1;
            quat[2] = cos(theta1)*t1;
            quat[3] = sin(theta2)*t2;
          }
          count++;
        }

    } else {
      double theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          else if (avec_body && atom->body[i] >= 0)
            quat = avec_body->bonus[atom->body[i]].quat;
          else if (atom->quat_flag)
            quat2 = atom->quat;
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          ranpark->reset(seed,x[i]);
          theta2 = MY_PI*ranpark->uniform();
          if (atom->quat_flag) {
            quat2[i][0] = cos(theta2);
            quat2[i][1] = 0.0;
            quat2[i][2] = 0.0;
            quat2[i][3] = sin(theta2);
          } else {
            quat[0] = cos(theta2);
            quat[1] = 0.0;
            quat[2] = 0.0;
            quat[3] = sin(theta2);
          }
          count++;
        }
    }

  // set theta to random orientation in 2d

  } else if (keyword == THETA_RANDOM) {
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (select[i]) {
        if (atom->line[i] < 0)
          error->one(FLERR,"Cannot set theta for atom that is not a line");
        ranpark->reset(seed,x[i]);
        avec_line->bonus[atom->line[i]].theta = MY_2PI*ranpark->uniform();
        count++;
      }
    }
  }

  delete ranpark;
  delete ranmars;
}

/* ---------------------------------------------------------------------- */

void Set2::topology(int keyword)
{
  int m,atom1,atom2,atom3,atom4;

  // error check

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR,"Cannot set bond topology types for atom style template");

  // border swap to acquire ghost atom info
  // enforce PBC before in case atoms are outside box
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0) utils::logmesg(lmp,"  system init for set ...\n");
  lmp->init();

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // select both owned and ghost atoms

  selection(atom->nlocal + atom->nghost);

  // for BOND, each of 2 atoms must be in group

  if (keyword == BOND) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++) {
        atom1 = atom->map(atom->bond_atom[i][m]);
        if (atom1 == -1) error->one(FLERR,"Bond atom missing in set command");
        if (select[i] && select[atom1]) {
          atom->bond_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for ANGLE, each of 3 atoms must be in group

  if (keyword == ANGLE) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++) {
        atom1 = atom->map(atom->angle_atom1[i][m]);
        atom2 = atom->map(atom->angle_atom2[i][m]);
        atom3 = atom->map(atom->angle_atom3[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Angle atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3]) {
          atom->angle_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for DIHEDRAL, each of 4 atoms must be in group

  if (keyword == DIHEDRAL) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++) {
        atom1 = atom->map(atom->dihedral_atom1[i][m]);
        atom2 = atom->map(atom->dihedral_atom2[i][m]);
        atom3 = atom->map(atom->dihedral_atom3[i][m]);
        atom4 = atom->map(atom->dihedral_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Dihedral atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->dihedral_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for IMPROPER, each of 4 atoms must be in group

  if (keyword == IMPROPER) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++) {
        atom1 = atom->map(atom->improper_atom1[i][m]);
        atom2 = atom->map(atom->improper_atom2[i][m]);
        atom3 = atom->map(atom->improper_atom3[i][m]);
        atom4 = atom->map(atom->improper_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Improper atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->improper_type[i][m] = ivalue;
          count++;
        }
      }
  }
}

/* ---------------------------------------------------------------------- */

void Set2::varparse(const char *name, int m)
{
  varflag = 1;
  int ivar = input->variable->find(name+2);

  if (ivar < 0)
    error->all(FLERR,"Variable name {} for set command does not exist", name);
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable {} for set command is invalid style", name);

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  } else if (m == 4) {
    varflag4 = 1; ivar4 = ivar;
  }
}

// ----------------------------------------------------------------------
// pairs of process/invoke methods for each keyword
// process method reads args, stores parameters in Action instance
// invoke method resets atoms properties using Action instance
// separate two operations so can be called by either set or fix set command
// ----------------------------------------------------------------------

void Set2::process_type(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set type", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
    ivalue = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
    delete[] typestr;
  }
  iarg += 2;
}

void Set2::invoke_type()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (varflag)
      if (ivalue <= 0 || ivalue > atom->ntypes)
        error->one(FLERR,"Invalid value in set command");
    type[i] = ivalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_type_fraction(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/fraction", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  newtype = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  fraction = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  ivalue = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (newtype <= 0 || newtype > atom->ntypes)
    error->all(FLERR,"Invalid type value {} in set type/fraction command", newtype);
  if (fraction < 0.0 || fraction > 1.0)
    error->all(FLERR,"Invalid fraction value {} in set type/fraction command", fraction);
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed {} in set type/fraction command", ivalue);
  iarg += 4;
}

void Set2::invoke_type_fraction()
{
  setrandom(TYPE_FRACTION);
}

/* ---------------------------------------------------------------------- */

void Set2::process_type_ratio(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/ratio", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  newtype = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  fraction = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  ivalue = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (newtype <= 0 || newtype > atom->ntypes)
    error->all(FLERR,"Invalid type value {} in set type/ratio command", newtype);
  if (fraction < 0.0 || fraction > 1.0)
    error->all(FLERR,"Invalid fraction value {} in set type/ratio command", fraction);
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed {} in set type/ratio command", ivalue);
  iarg += 4;
}

void Set2::invoke_type_ratio()
{
  setrandom(TYPE_RATIO);
}

/* ---------------------------------------------------------------------- */

void Set2::process_type_subset(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set type/subset", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ATOM,lmp);
  newtype = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  nsubset = utils::bnumeric(FLERR,arg[iarg+2],false,lmp);
  ivalue = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  if (newtype <= 0 || newtype > atom->ntypes)
    error->all(FLERR,"Invalid type value {} in set type/subset command", newtype);
  if (nsubset < 0)
    error->all(FLERR,"Invalid subset size {} in set type/subset command", nsubset);
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed {} in set type/subset command", ivalue);
  iarg += 4;
}

void Set2::invoke_type_subset()
{
  setrandom(TYPE_SUBSET);
}

/* ---------------------------------------------------------------------- */

void Set2::process_mol(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set mol", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->molecule_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}
 
void Set2::invoke_mol()
{
  int nlocal = atom->nlocal;
  int *molecule = atom->molecule;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    molecule[i] = ivalue;
  }
}
 
/* ---------------------------------------------------------------------- */

void Set2::process_x(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set x", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_x()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    x[i][0] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_y(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set y", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_y()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    x[i][1] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_z(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set z", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_z()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    x[i][2] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_vx(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vx", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_vx()
{
  int nlocal = atom->nlocal;
  double **v = atom->v;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    v[i][0] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_vy(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vy", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_vy()
{
  int nlocal = atom->nlocal;
  double **v = atom->v;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    v[i][1] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_vz(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set vz", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  iarg += 2;
}

void Set2::invoke_vz()
{
  int nlocal = atom->nlocal;
  double **v = atom->v;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    v[i][2] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_charge(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set charge", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->q_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_charge()
{
  int nlocal = atom->nlocal;
  double *q = atom->q;
  double *q_scaled = atom->q_scaled;
  double *epsilon = atom->epsilon;

  // ensure that scaled charges are consistent the new charge value

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    q[i] = dvalue;
    if (epsilon) q_scaled[i] = dvalue / epsilon[i];
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_mass(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set mass", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->rmass_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_mass()
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue <= 0.0) error->one(FLERR,"Invalid mass in set command");
    rmass[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_shape(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set shape", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (!atom->ellipsoid_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 4;
}

void Set2::invoke_shape()
{
  int nlocal = atom->nlocal;
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
      if (xvalue < 0.0 || yvalue < 0.0 || zvalue < 0.0)
        error->one(FLERR,"Invalid shape in set command");
      if (xvalue > 0.0 || yvalue > 0.0 || zvalue > 0.0) {
        if (xvalue == 0.0 || yvalue == 0.0 || zvalue == 0.0)
          error->one(FLERR,"Invalid shape in set command");
      }
      avec_ellipsoid->set_shape(i,0.5*xvalue,0.5*yvalue,0.5*zvalue);
  }
  
  // update global ellipsoid count

  bigint nlocal_bonus = avec_ellipsoid->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->nellipsoids,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set2::process_length(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set length", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_length()
{
  int nlocal = atom->nlocal;
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0) error->one(FLERR,"Invalid length in set command");
    avec_line->set_length(i,dvalue);
  }

  // update global line count

  bigint nlocal_bonus = avec_line->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->nlines,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set2::process_tri(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set tri", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->tri_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_tri()
{
  int nlocal = atom->nlocal;
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0) error->one(FLERR,"Invalid tri size in set command");
    avec_tri->set_equilateral(i,dvalue);
  }

  // update bonus tri count

  bigint nlocal_bonus = avec_tri->nlocal_bonus;
  MPI_Allreduce(&nlocal_bonus,&atom->ntris,1,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Set2::process_dipole(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set dipole", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (!atom->mu_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 4;
}

void Set2::invoke_dipole()
{
  int nlocal = atom->nlocal;
  double **mu = atom->mu;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    mu[i][0] = xvalue;
    mu[i][1] = yvalue;
    mu[i][2] = zvalue;
    mu[i][3] = sqrt(mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2]);
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_dipole_random(int &iarg, int narg, char **arg)
{
  if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set dipole/random", error);
  ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  dvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (!atom->mu_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed in set command");
  if (dvalue <= 0.0)
    error->all(FLERR,"Invalid dipole length in set command");
  iarg += 3;
}

void Set2::invoke_dipole_random()
{
  setrandom(DIPOLE_RANDOM);
}

/* ---------------------------------------------------------------------- */

void Set2::process_spin_atom(int &iarg, int narg, char **arg)
{
  if ((strcmp(arg[iarg],"spin") == 0) && (comm->me == 0))
    error->warning(FLERR, "Set attribute spin is deprecated. Please use spin/atom instead.");
  if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "set spin/atom", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else xvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else yvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (utils::strmatch(arg[iarg+4],"^v_")) varparse(arg[iarg+4],4);
  else zvalue = utils::numeric(FLERR,arg[iarg+4],false,lmp);
  if ((xvalue == 0.0) && (yvalue == 0.0) && (zvalue == 0.0))
    error->all(FLERR,"At least one spin vector component must be non-zero");
  if (!atom->sp_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (dvalue <= 0.0)
    error->all(FLERR,"Invalid spin magnitude {} in set {} command", dvalue, arg[iarg]);
  iarg += 5;
}

void Set2::invoke_spin_atom()
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;

  double norm;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0)
      error->one(FLERR,"Incorrect value for atomic spin magnitude: {}", dvalue);
    norm = 1.0/sqrt(xvalue*xvalue+yvalue*yvalue+zvalue*zvalue);
    sp[i][0] = norm*xvalue;
    sp[i][1] = norm*yvalue;
    sp[i][2] = norm*zvalue;
    sp[i][3] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_spin_atom_random(int &iarg, int narg, char **arg)
{
  if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set spin/atom/random", error);
  ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  dvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if ((strcmp(arg[iarg],"spin/random") == 0) && (comm->me == 0))
    error->warning(FLERR, "Set attribute spin/random is deprecated. "
                   "Please use spin/atom/random instead.");
  if (!atom->sp_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed {} in set {} command", ivalue, arg[iarg]);
  if (dvalue <= 0.0)
    error->all(FLERR,"Invalid spin magnitude {} in set {} command", dvalue, arg[iarg]);
  iarg += 3;
}

void Set2::invoke_spin_atom_random()
{
  setrandom(SPIN_RANDOM);
}

/* ---------------------------------------------------------------------- */

void Set2::process_radius_election(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set radius/electron", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->eradius_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_radius_election()
{
  int nlocal = atom->nlocal;
  double *eradius = atom->eradius;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0)
      error->one(FLERR,"Incorrect value for electron radius: {}", dvalue);
    eradius[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_spin_electron(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set spin/electron", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->spin_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_spin_electron()
{
  int nlocal = atom->nlocal;
  int *spin = atom->spin;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (ivalue < -1 || ivalue > 3)
      error->one(FLERR,"Incorrect value for electron spin: {}", ivalue);
    atom->spin[i] = ivalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_quat(int &iarg, int narg, char **arg)
{
  if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "set quat", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (utils::strmatch(arg[iarg+4],"^v_")) varparse(arg[iarg+4],4);
  else wvalue = utils::numeric(FLERR,arg[iarg+4],false,lmp);
  if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag && !atom->quat_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 5;
}

void Set2::invoke_quat()
{
  int nlocal = atom->nlocal;
  int *ellipsoid = atom->ellipsoid;
  int *tri = atom->tri;
  int *body = atom->body;
  double **quat = atom->quat;
  int quat_flag = atom->quat_flag;
  
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));

  int dimension = domain->dimension;
  double theta2,sintheta2;
  double *quat_one;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

      if (avec_ellipsoid && ellipsoid[i] >= 0)
        quat_one = avec_ellipsoid->bonus[ellipsoid[i]].quat;
      else if (avec_tri && tri[i] >= 0)
        quat_one = avec_tri->bonus[tri[i]].quat;
      else if (avec_body && body[i] >= 0)
        quat_one = avec_body->bonus[body[i]].quat;
      else if (quat_flag)
        quat_one = quat[i];
      else
        error->one(FLERR,"Cannot set quaternion for atom that has none");

      // quat rotation vector must be only in z dir for 2d systems

      if (dimension == 2 && (xvalue != 0.0 || yvalue != 0.0))
        error->one(FLERR,"Cannot set quaternion with xy components for 2d system");

      theta2 = MY_PI2 * wvalue/180.0;
      sintheta2 = sin(theta2);
      quat_one[0] = cos(theta2);
      quat_one[1] = xvalue * sintheta2;
      quat_one[2] = yvalue * sintheta2;
      quat_one[3] = zvalue * sintheta2;
      MathExtra::qnormalize(quat_one);
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_quat_random(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set quat/random", error);
  ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag && !atom->quat_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed in set command");
  iarg += 2;
}

void Set2::invoke_quat_random()
{
  setrandom(QUAT_RANDOM);
}

/* ---------------------------------------------------------------------- */

void Set2::process_theta(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set theta", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = DEG2RAD * utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_theta()
{
  int nlocal = atom->nlocal;
  int *line = atom->line;

  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (line[i] < 0)
      error->one(FLERR,"Cannot set theta for atom that is not a line");
    avec_line->bonus[atom->line[i]].theta = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_theta_random(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set theta/random", error);
  ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->line_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0)
    error->all(FLERR,"Invalid random number seed in set command");
  iarg += 2;
}

void Set2::invoke_theta_random()
{
  set(THETA_RANDOM);
}

/* ---------------------------------------------------------------------- */

void Set2::process_angmom(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set angmom", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (!atom->angmom_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 4;
}

void Set2::invoke_angmom()
{
  int nlocal = atom->nlocal;
  double **angmom = atom->angmom;

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    angmom[i][0] = xvalue;
    angmom[i][1] = yvalue;
    angmom[i][2] = zvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_omega(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set omega", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
  else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
  else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
  if (!atom->omega_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 4;
}

void Set2::invoke_omega()
{
  int nlocal = atom->nlocal;
  double **omega = atom->angmom;

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    omega[i][0] = xvalue;
    omega[i][1] = yvalue;
    omega[i][2] = zvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_diameter(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set diameter", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->radius_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_diameter()
{
  int nlocal = atom->nlocal;
  double *radius = atom->radius;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0) error->one(FLERR,"Invalid diameter in set command");
    radius[i] = 0.5 * dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_density(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set density", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->rmass_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (dvalue <= 0.0) error->all(FLERR,"Invalid density in set command");
  discflag = 0;
  if (strcmp(arg[iarg],"density/disc") == 0) {
    discflag = 1;
    if (domain->dimension != 2) error->all(FLERR,"Set density/disc requires 2d simulation");
  }
  iarg += 2;
}

void Set2::invoke_density()
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;

  int radius_flag = atom->radius_flag;
  int ellipsoid_flag = atom->ellipsoid_flag;
  int line_flag = atom->line_flag;
  int tri_flag = atom->tri_flag;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));

  // set rmass via density
  // if radius > 0.0, treat as sphere or disc
  // if shape > 0.0, treat as ellipsoid (or ellipse, when uncomment below)
  // if length > 0.0, treat as line
  // if area > 0.0, treat as tri
  // else set rmass to density directly

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue <= 0.0) error->one(FLERR,"Invalid density in set command");
    if (radius_flag && radius[i] > 0.0)
      if (discflag) rmass[i] = MY_PI*radius[i]*radius[i] * dvalue;
      else rmass[i] = 4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i] * dvalue;
    else if (ellipsoid_flag && ellipsoid[i] >= 0) {
      double *shape = avec_ellipsoid->bonus[ellipsoid[i]].shape;
      // enable 2d ellipse (versus 3d ellipsoid) when time integration
      //   options (fix nve/asphere, fix nh/asphere) are also implemented
      // if (discflag)
      // atom->rmass[i] = MY_PI*shape[0]*shape[1] * dvalue;
      // else
      rmass[i] = 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2] * dvalue;
    } else if (line_flag && line[i] >= 0) {
      double length = avec_line->bonus[line[i]].length;
      rmass[i] = length * dvalue;
    } else if (tri_flag && tri[i] >= 0) {
      double *c1 = avec_tri->bonus[tri[i]].c1;
      double *c2 = avec_tri->bonus[tri[i]].c2;
      double *c3 = avec_tri->bonus[tri[i]].c3;
      double c2mc1[3],c3mc1[3];
      MathExtra::sub3(c2,c1,c2mc1);
      MathExtra::sub3(c3,c1,c3mc1);
      double norm[3];
      MathExtra::cross3(c2mc1,c3mc1,norm);
      double area = 0.5 * MathExtra::len3(norm);
      rmass[i] = area * dvalue;
    } else rmass[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_temperature(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->temperature_flag)
    error->all(FLERR,"Cannot set this attribute for this atom style");
  iarg += 2;
}

void Set2::invoke_temperature()
{
  int nlocal = atom->nlocal;
  double *temperature = atom->temperature;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0) error->one(FLERR,"Invalid temperature in set command");
    temperature[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_volume(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set volume", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->vfrac_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (dvalue <= 0.0) error->all(FLERR,"Invalid volume in set command");
  iarg += 2;
}

void Set2::invoke_volume()
{
  int nlocal = atom->nlocal;
  double *vfrac = atom->vfrac;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    if (dvalue < 0.0) error->one(FLERR,"Invalid diameter in set command");
    vfrac[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_image(int &iarg, int narg, char **arg)
{
  if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set image", error);
  ximageflag = yimageflag = zimageflag = 0;
  if (strcmp(arg[iarg+1],"NULL") != 0) {
    ximageflag = 1;
    if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
    else ximage = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
  }
  if (strcmp(arg[iarg+2],"NULL") != 0) {
    yimageflag = 1;
    if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
    else yimage = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
  }
  if (strcmp(arg[iarg+3],"NULL") != 0) {
    zimageflag = 1;
    if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
    else zimage = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
  }
  if (ximageflag && ximage && !domain->xperiodic)
    error->all(FLERR,
               "Cannot set non-zero image flag for non-periodic dimension");
  if (yimageflag && yimage && !domain->yperiodic)
    error->all(FLERR,
               "Cannot set non-zero image flag for non-periodic dimension");
  if (zimageflag && zimage && !domain->zperiodic)
    error->all(FLERR,
               "Cannot set non-zero image flag for non-periodic dimension");
  iarg += 4;
}

void Set2::invoke_image()
{
  int nlocal = atom->nlocal;
  imageint *image = atom->image;

  int xbox,ybox,zbox;

  // reset any or all of 3 image flags

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;
    if (varflag1) ximage = static_cast<int>(xvalue);
    if (varflag2) yimage = static_cast<int>(yvalue);
    if (varflag3) zimage = static_cast<int>(zvalue);
    if (ximageflag) xbox = ximage;
    if (yimageflag) ybox = yimage;
    if (zimageflag) zbox = zimage;
    image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
      (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
      (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_bond(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set bond", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::BOND,lmp);
  ivalue = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0 || ivalue > atom->nbondtypes)
    error->all(FLERR,"Invalid value in set command");
  iarg += 2;
}

void Set2::invoke_bond()
{
  topology(BOND);
}

/* ---------------------------------------------------------------------- */

void Set2::process_angle(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set angle", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::ANGLE,lmp);
  ivalue = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0 || ivalue > atom->nangletypes)
    error->all(FLERR,"Invalid value in set command");
  iarg += 2;
}

void Set2::invoke_angle()
{
  topology(ANGLE);
}

/* ---------------------------------------------------------------------- */

void Set2::process_dihedral(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set dihedral", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::DIHEDRAL,lmp);
  ivalue = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0 || ivalue > atom->ndihedraltypes)
    error->all(FLERR,"Invalid value in set command");
  iarg += 2;
}

void Set2::invoke_dihedral()
{
  topology(DIHEDRAL);
}

/* ---------------------------------------------------------------------- */

void Set2::process_improper(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set improper", error);
  char *typestr = utils::expand_type(FLERR,arg[iarg+1],Atom::IMPROPER,lmp);
  ivalue = utils::inumeric(FLERR,typestr?typestr:arg[iarg+1],false,lmp);
  delete[] typestr;
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  if (ivalue <= 0 || ivalue > atom->nimpropertypes)
    error->all(FLERR,"Invalid value in set command");
  iarg += 2;
}

void Set2::invoke_improper()
{
  topology(IMPROPER);
}

/* ---------------------------------------------------------------------- */

void Set2::process_sph_e(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/e", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->esph_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_sph_e()
{
  int nlocal = atom->nlocal;
  double *esph = atom->esph;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    esph[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_sph_cv(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/cv", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->cv_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_sph_cv()
{
  int nlocal = atom->nlocal;
  double *cv = atom->cv;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    cv[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_sph_rho(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/rho", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->rho_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_sph_rho()
{
  int nlocal = atom->nlocal;
  double *rho = atom->rho;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    rho[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_edpd_temp(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/temp", error);
  if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
  }
  if (!atom->edpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_edpd_temp()
{
  int nlocal = atom->nlocal;
  double *edpd_temp = atom->edpd_temp;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    edpd_temp[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_edpd_cv(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/cv", error);
  if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
  }
  if (!atom->edpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_edpd_cv()
{
  int nlocal = atom->nlocal;
  double *edpd_cv = atom->edpd_cv;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    edpd_cv[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_cc(int &iarg, int narg, char **arg)
{
  if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set cc", error);
  // NOTE: seems like error to not set cc_index for all 3 cases
  //       syntax is iarg+1 is index, iarg+2 is value
  //       doc page does not talk about NULL as valid value
  //       check package DPD-MESO package examples for tDPD for use of this command
  if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    cc_index = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    dvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    if (cc_index < 1) error->all(FLERR,"Illegal set command");
  }
  if (!atom->tdpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 3;
}

void Set2::invoke_cc()
{
  int nlocal = atom->nlocal;
  double **cc = atom->cc;

  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    cc[i][cc_index-1] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_smd_mass_density(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set smd/mass/density", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->smd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_smd_mass_density()
{
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *vfrac = atom->vfrac;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    // set mass from volume and supplied mass density
    rmass[i] = vfrac[i] * dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_smd_contact_radius(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set smd/contact/radius", error);
  if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  if (!atom->smd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_smd_contact_radius()
{
  int nlocal = atom->nlocal;
  double *contact_radius = atom->contact_radius;
  
  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    contact_radius[i] = dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_dpd_theta(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set dpd/theta", error);
  if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
  }
  if (!atom->dpd_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_dpd_theta()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *dpdTheta = atom->dpdTheta;
  
  double tfactor = force->mvv2e / (domain->dimension * force->boltz);
  double onemass;
  double vx,vy,vz;
  
  // VARIABLE option or NULL to set temp to KE of particle
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    if (dvalue >= 0.0) dpdTheta[i] = dvalue;
    else {
      if (rmass) onemass = rmass[i];
      else onemass = mass[type[i]];
      vx = v[i][0];
      vy = v[i][1];
      vz = v[i][2];
      dpdTheta[i] = tfactor * onemass * (vx*vx + vy*vy + vz*vz);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_epsilon(int &iarg, int narg, char **arg)
{
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set epsilon", error);
  if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
  else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
  else {
    dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
  }
  if (!atom->dielectric_flag)
    error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
  iarg += 2;
}

void Set2::invoke_epsilon()
{
  int nlocal = atom->nlocal;
  double *epsilon = atom->epsilon;
  double *q = atom->q;
  double *q_scaled = atom->q_scaled;

  // assign local dielectric constant
  // also update scaled charge value

  // VARIABLE option
  
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;
    // NOTE: should it be error if dvalue < 0.0 ?
    if (dvalue >= 0.0) {
      epsilon[i] = dvalue;
      q_scaled[i] = q[i] / dvalue;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Set2::process_custom(int &iarg, int narg, char **arg)
{
  int flag,cols;
  ArgInfo argi(arg[iarg],ArgInfo::DNAME|ArgInfo::INAME);
  const char *pname = argi.get_name();
  if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set", error);
  index_custom = atom->find_custom(argi.get_name(),flag,cols);
  if (index_custom < 0)
    error->all(FLERR,"Set keyword or custom property {} does not exist",pname);

  switch (argi.get_type()) {

  case ArgInfo::INAME:
    if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
    else ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    if (flag != 0) error->all(FLERR,"Set command custom property {} is not integer",pname);
    
    if (argi.get_dim() == 0) {
      if (cols > 0)
        error->all(FLERR,"Set command custom integer property {} is not a vector",pname);
      custom_flag = IVEC;
    } else if (argi.get_dim() == 1) {
      if (cols == 0)
        error->all(FLERR,"Set command custom integer property {} is not an array",pname);
      icol_custom = argi.get_index1();
      if (icol_custom <= 0 || icol_custom > cols)
        error->all(FLERR,"Set command per-atom custom integer array {} is accessed "
                   "out-of-range",pname);
      custom_flag = IARRAY;
    } else error->all(FLERR,"Illegal set command");
    break;

  case ArgInfo::DNAME:
    if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
    else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (flag != 1) error->all(FLERR,"Custom property {} is not floating-point",argi.get_name());
    
    if (argi.get_dim() == 0) {
      if (cols > 0)
        error->all(FLERR,"Set command custom double property {} is not a vector",pname);
      custom_flag = DVEC;
    } else if (argi.get_dim() == 1) {
      if (cols == 0)
        error->all(FLERR,"Set command custom double property {} is not an array",pname);
      icol_custom = argi.get_index1();
      if (icol_custom <= 0 || icol_custom > cols)
        error->all(FLERR,"Set command per-atom custom double array {} is "
                   "accessed out-of-range",pname);
      custom_flag = DARRAY;
    } else error->all(FLERR,"Illegal set command");
    break;

  default:
    error->all(FLERR,"Illegal set command");
    break;
  }
      
  iarg += 2;
}

void Set2::invoke_custom()
{
  int nlocal = atom->nlocal;

  // VARIABLE option

  if (custom_flag == IVEC) {
    int *ivector = atom->ivector[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      ivector[i] = ivalue;
    }
  } else if (custom_flag == DVEC) {
    double *dvector = atom->dvector[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      dvector[i] = dvalue;
    }
  } else if (custom_flag == IARRAY) {
    int **iarray = atom->iarray[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      iarray[i][icol_custom-1] = ivalue;
    }
  } else if (custom_flag == DARRAY) {
    double **darray = atom->darray[index_custom];
    for (int i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      darray[i][icol_custom-1] = dvalue;
    }
  }
}
