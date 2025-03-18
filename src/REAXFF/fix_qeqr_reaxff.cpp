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

/* ----------------------------------------------------------------------
   Contributing authors:
      Navraj S Lalli, Imperial College London (navrajsinghlalli@gmail.com)

------------------------------------------------------------------------- */

#include "fix_qeqr_reaxff.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_efield.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "region.h"
#include "respa.h"
#include "text_file_reader.h"
#include "update.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double ANGSTROM_TO_BOHRRADIUS = 1.8897261259;

/* ---------------------------------------------------------------------- */

FixQEqrReaxFF::FixQEqrReaxFF(LAMMPS *lmp, int narg, char **arg) :
  FixQtpieReaxFF(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixQEqrReaxFF::calc_chi_eff()
{
  memset(&chi_eff[0],0,atom->nmax*sizeof(double));

  const auto x = (const double * const *)atom->x;
  const int *type = atom->type;

  double dist,overlap,sum_n,sum_d,expa,expb,chia,phia,phib,p,m;
  int i,j;

  // check ghost atoms are stored up to the distance cutoff for overlap integrals
  const double comm_cutoff = MAX(neighbor->cutneighmax,comm->cutghostuser);
  if(comm_cutoff < dist_cutoff/ANGSTROM_TO_BOHRRADIUS) {
    error->all(FLERR,"comm cutoff = {} Angstrom is smaller than distance cutoff = {} Angstrom "
               "for overlap integrals in {}. Increase comm cutoff with comm_modify",
               comm_cutoff, dist_cutoff/ANGSTROM_TO_BOHRRADIUS, style);
  }

  // efield energy is in real units of kcal/mol, factor needed for conversion to eV
  const double qe2f = force->qe2f;
  const double factor = 1.0/qe2f;

  if (efield) {
    if (efield->varflag != FixEfield::CONSTANT)
      efield->update_efield_variables();

    // compute chi_eff for each local atom
    for (i = 0; i < nn; i++) {
      expa = gauss_exp[type[i]];
      chia = chi[type[i]];
      if (efield->varflag != FixEfield::ATOM) {
        phia = -factor*(x[i][0]*efield->ex  + x[i][1]*efield->ey + x[i][2]*efield->ez);
      } else { // atom-style potential from FixEfield
        phia = efield->efield[i][3];
      }

      sum_n = 0.0;
      sum_d = 0.0;

      for (j = 0; j < nt; j++) {
        dist = distance(x[i],x[j])*ANGSTROM_TO_BOHRRADIUS; // in atomic units

        if (dist < dist_cutoff) {
          expb = gauss_exp[type[j]];

          // overlap integral of two normalised 1s Gaussian type orbitals
          p = expa + expb;
          m = expa * expb / p;
          overlap = pow((4.0*m/p),0.75) * exp(-m*dist*dist);

          if (efield->varflag != FixEfield::ATOM) {
            phib = -factor*(x[j][0]*efield->ex  + x[j][1]*efield->ey + x[j][2]*efield->ez);
          } else { // atom-style potential from FixEfield
            phib = efield->efield[j][3];
          }
          sum_n += (chia + scale * (phia - phib)) * overlap;
          sum_d += overlap;
        }
      }
      chi_eff[i] = sum_n / sum_d;
    }
  } else {
    for (i = 0; i < nn; i++) {
      chi_eff[i] = chi[type[i]];
    }
  }
}
