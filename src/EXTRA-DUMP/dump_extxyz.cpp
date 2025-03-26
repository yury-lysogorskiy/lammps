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

#include "dump_extxyz.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int ONELINE = 128;
static constexpr int DELTA = 1048576;

/* ---------------------------------------------------------------------- */

DumpExtXYZ::DumpExtXYZ(LAMMPS *lmp, int narg, char **arg) : DumpXYZ(lmp, narg, arg)
{

  size_one = 11;
  delete[] format_default;
  format_default = utils::strdup("%s %g %g %g %g %g %g %g %g %g");

  // use type labels by default if present
  if (atom->labelmapflag) {
    typenames = new char *[ntypes + 1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = utils::strdup(atom->lmap->typelabel[itype - 1]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::init_style()
{
  if (!typenames)
    error->all(FLERR, Error::NOLASTLINE,
               "Must use either typelables or dump_modify element with dump style extxyz");

  DumpXYZ::init_style();
}

/* ---------------------------------------------------------------------- */

int DumpExtXYZ::modify_param(int narg, char **arg)
{
  int rv = DumpXYZ::modify_param(narg, arg);
  if (rv > 0) return rv;
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::write_header(bigint n)
{
  if (me == 0) {
    if (!fp) error->one(FLERR, Error::NOLASTLINE, "Must not use 'run pre no' after creating a new dump");

    std::string header;

    header = fmt::format("{}\nTimestep={}", n, update->ntimestep);
    if (time_flag) header += fmt::format(" Time={:.6f}", compute_time());
    header += fmt::format(" pbc=\"{} {} {}\"", domain->xperiodic ? "T" : "F",
                          domain->yperiodic ? "T" : "F", domain->zperiodic ? "T" : "F");
    header +=
        fmt::format(" Lattice=\"{:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\"", domain->xprd, 0.,
                    0., domain->xy, domain->yprd, 0., domain->xz, domain->yz, domain->zprd);
    header += " Properties=species:S:1:pos:R:3:vel:R:3:forces:R:3";
    utils::print(fp, header + "\n");
  }
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::pack(tagint *ids)
{
  int m, n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
      buf[m++] = f[i][0];
      buf[m++] = f[i][1];
      buf[m++] = f[i][2];

      if (ids) ids[n++] = tag[i];
    }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpExtXYZ::convert_string(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }

    offset +=
        snprintf(&sbuf[offset], maxsbuf - offset, format, typenames[static_cast<int>(mybuf[m + 1])],
                 mybuf[m + 2], mybuf[m + 3], mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7],
                 mybuf[m + 8], mybuf[m + 9], mybuf[m + 10]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
            mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8], mybuf[m + 9],
            mybuf[m + 10]);
    m += size_one;
  }
}
