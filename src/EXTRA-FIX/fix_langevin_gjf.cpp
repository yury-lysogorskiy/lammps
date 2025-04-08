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

/* ----------------------------------------------------------------------
   Contributing authors: Tim Linke & Niels Gronbech-Jensen (UC Davis)
------------------------------------------------------------------------- */

#include "fix_langevin.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NOBIAS, BIAS };
enum { CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixLangevin::FixLangevin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), gjfflag(0), gfactor1(nullptr), gfactor2(nullptr), ratio(nullptr),
    tstr(nullptr), flangevin(nullptr), tforce(nullptr), lv(nullptr), id_temp(nullptr), random(nullptr)
{
  if (narg < 8) error->all(FLERR, "Illegal fix langevin/gjf command");

  time_integrate = 1;
  restart_peratom = 1;
  // dynamic_group_allow = 1;
  // scalar_flag = 1;
  // global_freq = 1;
  // extscalar = 1;
  // ecouple_flag = 1;
  nevery = 1;

  if (utils::strmatch(arg[3], "^v_")) {
    tstr = utils::strdup(arg[3] + 2);
  } else {
    t_start = utils::numeric(FLERR, arg[3], false, lmp);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = utils::numeric(FLERR, arg[4], false, lmp);
  t_period = utils::numeric(FLERR, arg[5], false, lmp);
  seed = utils::inumeric(FLERR, arg[6], false, lmp);

  if (t_period <= 0.0) error->all(FLERR, "Fix langevin/gjf period must be > 0.0");
  if (seed <= 0) error->all(FLERR, "Illegal fix langevin/gjf command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp, seed + comm->me);

  // allocate per-type arrays for force prefactors

  // gfactor1 = new double[atom->ntypes + 1];
  // gfactor2 = new double[atom->ntypes + 1];
  // ratio = new double[atom->ntypes + 1];
  int GJmethods = 8 // number of currently implemented GJ methods

  // optional args

  for (int i = 1; i <= atom->ntypes; i++) ratio[i] = 1.0;
  osflag = 0;
  GJmethod = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "vel") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix langevin/gjf command");
      if (strcmp(arg[iarg + 1], "vfull") == 0) {
        osflag = 1;
      } else if (strcmp(arg[iarg + 1], "vhalf") == 0) {
        osflag = 0;
      } else
        error->all(FLERR, "Illegal fix langevin/gjf command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "method") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix langevin/gjf command");
      GJmethod = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (GJmethod <= 0 || GJmethod > GJmethods) error->all(FLERR, "Invalid GJ method choice in langevin/gjf command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix langevin/gjf command");
  }

  // set temperature = nullptr, user can override via fix_modify if wants bias

  id_temp = nullptr;
  temperature = nullptr;

  energy = 0.0;

  // flangevin is unallocated until first call to setup()
  // compute_scalar checks for this and returns 0.0
  // if flangevin_allocated is not set

  flangevin = nullptr;
  flangevin_allocated = 0;
  lv = nullptr;
  tforce = nullptr;
  maxatom1 = maxatom2 = 0;

  // setup atom-based array for lv
  // register with Atom class
  // no need to set peratom_flag, b/c data is for internal use only

  
  FixLangevin::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // initialize lv to zero

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    lv[i][0] = 0.0;
    lv[i][1] = 0.0;
    lv[i][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixLangevin::~FixLangevin()
{
  if (copymode) return;

  delete random;
  delete[] tstr;
  delete[] gfactor1;
  delete[] gfactor2;
  delete[] ratio;
  delete[] id_temp;
  memory->destroy(flangevin);
  memory->destroy(tforce);

  memory->destroy(lv);
  if (modify->get_fix_by_id(id)) atom->delete_callback(id, Atom::GROW);
}

/* ---------------------------------------------------------------------- */

int FixLangevin::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::init()
{
  if (id_temp) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature) {
      error->all(FLERR, "Temperature compute ID {} for fix {} does not exist", id_temp, style);
    } else {
      if (temperature->tempflag == 0)
        error->all(FLERR, "Compute ID {} for fix {} does not compute temperature", id_temp, style);
    }
  }
  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR, "Variable name {} for fix langevin does not exist", tstr);
    if (input->variable->equalstyle(tvar))
      tstyle = EQUAL;
    else if (input->variable->atomstyle(tvar))
      tstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix langevin is invalid style", tstr);
  }

  // set force prefactors

  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
      gfactor2[i] = sqrt(atom->mass[i]) / force->ftm2v;
      gfactor2[i] *= sqrt(2.0 * update->dt * force->boltz / t_period / force->mvv2e); // gjfflag
    }
  }

  if (temperature && temperature->tempbias)
    tbiasflag = BIAS;
  else
    tbiasflag = NOBIAS;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    nlevels_respa = (static_cast<Respa *>(update->integrate))->nlevels;
    if (gjfflag) error->all(FLERR, "Fix langevin gjf and run style respa are not compatible");
  }

  if (gjfflag) {
    gjfc2 = (1.0 - update->dt / 2.0 / t_period) / (1.0 + update->dt / 2.0 / t_period);
    gjfc1 = 1.0 / (1.0 + update->dt / 2.0 / t_period);
  }

  switch (GJmethod) {
    case 1:
      gjfc2 = (1.0 - update->dt / 2.0 / t_period) / (1.0 + update->dt / 2.0 / t_period);
      gjfc1 = 1.0 / (1.0 + update->dt / 2.0 / t_period);
      break;
    case 2:
      // Insert logic for method 2
      break;
    case 3:
      // Insert logic for method 3
      break;
    case 4:
      // Insert logic for method 4
      break;
    case 5:
      // Insert logic for method 5
      break;
    case 6:
      // Insert logic for method 6
      break;
    case 7:
      // Insert logic for method 7
      break;
    case 8:
      // Insert logic for method 8
      break;
    default:
      error->all(FLERR, "Fix langevin/gjf method not found");
      break;
}
}

/* ----------------------------------------------------------------------
  integrate position and velocity according to the GJF method
  in Grønbech-Jensen, J Stat Phys 191, 137 (2024). The general workflow is
    1. Langevin GJF Initial Integration
    2. Force Update
    3. Langevin GJF Final Integration
    4. Velocity Choice in end_of_step()
------------------------------------------------------------------------- */

void FixLangevin::initial_integrate(int /* vflag */)
{
  double gamma1,gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fran[3];

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double dtf = 0.5 * dt * ftm2v;
  double dtfm;
  double c1sqrt = sqrt(gjfc1);
  
  // NVE integrates position and velocity according to Eq. 8a, 8b
  // This function embeds the GJF formulation into the NVE framework, which corresponds to the GJF case c1=c3.

  //NVE
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
      }
  }

  // The initial NVE integration should always use the on-site velocity. Therefore, a velocity correction
  // must be done when using the half-step option.
  //----------
  if (!osflag) {
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          // Undo NVE integration
          x[i][0] -= dt * v[i][0];
          x[i][1] -= dt * v[i][1];
          x[i][2] -= dt * v[i][2];
          // Obtain Eq. 24a. lv[][] stores on-site velocity from previous timestep
          v[i][0] = lv[i][0] + dtfm * f[i][0];
          v[i][1] = lv[i][1] + dtfm * f[i][1];
          v[i][2] = lv[i][2] + dtfm * f[i][2];
          // Redo NVE integration with correct velocity
          x[i][0] += dt * v[i][0];
          x[i][1] += dt * v[i][1];
          x[i][2] += dt * v[i][2];
        }
  
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / mass[type[i]];
          // Undo NVE integration
          x[i][0] -= dt * v[i][0];
          x[i][1] -= dt * v[i][1];
          x[i][2] -= dt * v[i][2];
          // Obtain Eq. 24a. lv[][] stores on-site velocity from previous timestep
          v[i][0] = lv[i][0] + dtfm * f[i][0];
          v[i][1] = lv[i][1] + dtfm * f[i][1];
          v[i][2] = lv[i][2] + dtfm * f[i][2];
          // Redo NVE integration with correct velocity
          x[i][0] += dt * v[i][0];
          x[i][1] += dt * v[i][1];
          x[i][2] += dt * v[i][2];
        }
    }
  }
  //----------

  compute_target();

  if (tbiasflag == BIAS) temperature->compute_scalar();
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
      if (rmass) {
        gamma2 = sqrt(rmass[i]) * sqrt(2.0*dt*boltz/t_period/mvv2e) / ftm2v;
        gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
      } else {
        gamma2 = gfactor2[type[i]] * tsqrt;
      }
      fran[0] = gamma2*random->gaussian();
      fran[1] = gamma2*random->gaussian();
      fran[2] = gamma2*random->gaussian();
      
      // NVE integrator delivers Eq. 24a, but also overshoots position integration. Calculate Eq. 24b:
      x[i][0] -= 0.5 * dt * v[i][0];
      x[i][1] -= 0.5 * dt * v[i][1];
      x[i][2] -= 0.5 * dt * v[i][2];
      // Calculate Eq. 24c:
      if (tbiasflag == BIAS)
        temperature->remove_bias(i,v[i]);
      if (rmass) {
        lv[i][0] = c1sqrt*v[i][0] + ftm2v * (c1sqrt / (2.0 * rmass[i])) * fran[0];
        lv[i][1] = c1sqrt*v[i][1] + ftm2v * (c1sqrt / (2.0 * rmass[i])) * fran[1];
        lv[i][2] = c1sqrt*v[i][2] + ftm2v * (c1sqrt / (2.0 * rmass[i])) * fran[2];
      } else {
        lv[i][0] = c1sqrt*v[i][0] + ftm2v * (c1sqrt / (2.0 * mass[type[i]])) * fran[0];
        lv[i][1] = c1sqrt*v[i][1] + ftm2v * (c1sqrt / (2.0 * mass[type[i]])) * fran[1];
        lv[i][2] = c1sqrt*v[i][2] + ftm2v * (c1sqrt / (2.0 * mass[type[i]])) * fran[2];
      }
      if (tbiasflag == BIAS)
        temperature->restore_bias(i,v[i]);
        if (tbiasflag == BIAS)
        temperature->restore_bias(i,lv[i]);
      
      // Calculate Eq. 24d
      if (tbiasflag == BIAS) temperature->remove_bias(i, lv[i]);
      if (atom->rmass) {
        v[i][0] = (gjfc2 / c1sqrt) * lv[i][0] + ftm2v * (0.5 / rmass[i]) * fran[0];
        v[i][1] = (gjfc2 / c1sqrt) * lv[i][1] + ftm2v * (0.5 / rmass[i]) * fran[1];
        v[i][2] = (gjfc2 / c1sqrt) * lv[i][2] + ftm2v * (0.5 / rmass[i]) * fran[2];
      } else {
        v[i][0] = (gjfc2 / c1sqrt) * lv[i][0] + ftm2v * (0.5 / mass[type[i]]) * fran[0];
        v[i][1] = (gjfc2 / c1sqrt) * lv[i][1] + ftm2v * (0.5 / mass[type[i]]) * fran[1];
        v[i][2] = (gjfc2 / c1sqrt) * lv[i][2] + ftm2v * (0.5 / mass[type[i]]) * fran[2];
      }
      if (tbiasflag == BIAS) temperature->restore_bias(i, lv[i]);
      // Calculate Eq. 24e. NVE integrator then calculates Eq. 24f.
      x[i][0] += 0.5 * dt * v[i][0];
      x[i][1] += 0.5 * dt * v[i][1];
      x[i][2] += 0.5 * dt * v[i][2];
    }
  }
}

void FixLangevin::final_integrate()
{
  double dtfm;
  double dt = update->dt;
  double ftm2v = force->ftm2v;
  double dtf = 0.5 * dt * ftm2v;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }
}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixLangevin::compute_target()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // if variable temp, evaluate variable, wrap with clear/add
  // reallocate tforce array if necessary

  if (tstyle == CONSTANT) {
    t_target = t_start + delta * (t_stop-t_start);
    tsqrt = sqrt(t_target);
  } else {
    modify->clearstep_compute();
    if (tstyle == EQUAL) {
      t_target = input->variable->compute_equal(tvar);
      if (t_target < 0.0)
        error->one(FLERR, "Fix langevin variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (atom->nmax > maxatom2) {
        maxatom2 = atom->nmax;
        memory->destroy(tforce);
        memory->create(tforce,maxatom2,"langevin:tforce");
      }
      input->variable->compute_atom(tvar,igroup,tforce,1,0);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
            if (tforce[i] < 0.0)
              error->one(FLERR, "Fix langevin variable returned negative temperature");
    }
    modify->addstep_compute(update->ntimestep + 1);
  }
}

/* ----------------------------------------------------------------------
   tally energy transfer to thermal reservoir, select velocity for GJF
------------------------------------------------------------------------- */

void FixLangevin::end_of_step()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  energy_onestep = 0.0;

  if (tallyflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        energy_onestep += flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
                          flangevin[i][2]*v[i][2];
  }

  energy += energy_onestep*update->dt;

  // After the NVE integrator delivers 24f, either the on-site or half-step
  // velocity is used in remaining simulation tasks, depending on user input
  if (gjfflag && !osflag) {
    double tmp[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        // v is Eq. 24f
        tmp[0] = v[i][0];
        tmp[1] = v[i][1];
        tmp[2] = v[i][2];
        // Move on with half-step velocity
        v[i][0] = lv[i][0];
        v[i][1] = lv[i][1];
        v[i][2] = lv[i][2];
        // store Eq. 24f in lv for next timestep
        lv[i][0] = tmp[0];
        lv[i][1] = tmp[1];
        lv[i][2] = tmp[2];
      }
  }
}

// clang-format on
/* ---------------------------------------------------------------------- */

void FixLangevin::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::reset_dt()
{
  if (atom->mass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor2[i] = sqrt(atom->mass[i]) / force->ftm2v;
      if (gjfflag)
        gfactor2[i] *= sqrt(2.0 * update->dt * force->boltz / t_period / force->mvv2e); // sqrt(2*alpha*kT*dt)
      else
        gfactor2[i] *= sqrt(24.0 * force->boltz / t_period / update->dt / force->mvv2e);
      gfactor2[i] *= 1.0 / sqrt(ratio[i]);
    }
  }
  if (gjfflag) {
    gjfc2 = (1.0 - update->dt / 2.0 / t_period) / (1.0 + update->dt / 2.0 / t_period);
    gjfc1 = 1.0 / (1.0 + update->dt / 2.0 / t_period);
  }
}

/* ---------------------------------------------------------------------- */

int FixLangevin::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR, "Could not find fix_modify temperature compute ID: {}", id_temp);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature compute {} does not compute temperature", id_temp);
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group: {} vs {}",
                     group->names[igroup], group->names[temperature->igroup]);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixLangevin::compute_scalar()
{
  if (!tallyflag || !flangevin_allocated) return 0.0;

  // capture the very first energy transfer to thermal reservoir

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (update->ntimestep == update->beginstep) {
    energy_onestep = 0.0;
    if (!gjfflag) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          energy_onestep +=
              flangevin[i][0] * v[i][0] + flangevin[i][1] * v[i][1] + flangevin[i][2] * v[i][2];
      energy = 0.5 * energy_onestep * update->dt;
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (tbiasflag) temperature->remove_bias(i, lv[i]);
          energy_onestep +=
              flangevin[i][0] * lv[i][0] + flangevin[i][1] * lv[i][1] + flangevin[i][2] * lv[i][2];
          if (tbiasflag) temperature->restore_bias(i, lv[i]);
        }
     energy = -0.5 * energy_onestep * update->dt;
    }
  }

  // convert midstep energy back to previous fullstep energy

  double energy_me = energy - 0.5 * energy_onestep * update->dt;

  double energy_all;
  MPI_Allreduce(&energy_me, &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return -energy_all;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixLangevin::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "t_target") == 0) { return &t_target; }
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of tally array
------------------------------------------------------------------------- */

double FixLangevin::memory_usage()
{
  double bytes = 0.0;
  if (gjfflag) bytes += (double) atom->nmax * 3 * sizeof(double);
  if (tallyflag || osflag) bytes += (double) atom->nmax * 3 * sizeof(double);
  if (tforce) bytes += (double) atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for lv
------------------------------------------------------------------------- */

void FixLangevin::grow_arrays(int nmax)
{
  memory->grow(lv, nmax, 3, "fix_langevin:lv");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixLangevin::copy_arrays(int i, int j, int /*delflag*/)
{
  lv[j][0] = lv[i][0];
  lv[j][1] = lv[i][1];
  lv[j][2] = lv[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixLangevin::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = lv[i][0];
  buf[n++] = lv[i][1];
  buf[n++] = lv[i][2];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixLangevin::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  lv[nlocal][0] = buf[n++];
  lv[nlocal][1] = buf[n++];
  lv[nlocal][2] = buf[n++];
  return n;
}
