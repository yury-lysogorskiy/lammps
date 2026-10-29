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
   grace/runtime and pace/runtime (CPU and KOKKOS): pair_style keywords, the model file of pair_coeff, the type map,
   the libgrace log callback and the rank count per node / GPU (grace_runtime_settings.h)
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#include "grace_runtime_settings.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "utils.h"

#include "grace/error.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void GraceRuntimeSettings::parse(LAMMPS *lmp, const std::string &style, unsigned backend, int narg, char **arg)
{
  *this = GraceRuntimeSettings();
  Error *error = lmp->error;
  for (int iarg = 0; iarg < narg;) {
    const Keyword *kw = nullptr;
    for (const auto &k : keywords)
      if (strcmp(arg[iarg], k.name) == 0) kw = &k;
    if (!kw) error->all(FLERR, "pair_style {}: unknown keyword {}", style, arg[iarg]);
    if (iarg + 1 + kw->nvalues > narg)
      utils::missing_cmd_args(FLERR, "pair_style " + style + " " + kw->name, error);
    char **v = arg + iarg + 1;
    iarg += 1 + kw->nvalues;
    if (!(kw->backends & backend)) {
      const char *other = backend == KK ? "the CPU" : "KOKKOS";
      // performance-only keywords of the other backend: ignored, so an input runs unchanged on either one
      if (!kw->perf_only) error->all(FLERR, "pair_style {}: keyword {} is only for {}", style, kw->name, other);
      if (lmp->comm->me == 0)
        error->warning(FLERR, "pair_style {}: keyword {} is only for {}, ignored", style, kw->name, other);
      continue;
    }
    given.push_back(kw);
    const bool aut = kw->nvalues > 0 && strcmp(v[0], "auto") == 0;
    switch (kw->id) {
      case SPLINE:
        spline = utils::inumeric(FLERR, v[0], false, lmp);
        if (spline != 0 && spline < 16)
          error->all(FLERR, "pair_style {}: spline N must be >= 16 (PACE: or 0)", style);
        break;
      case SPLINE_LAYERS:
        for (int k = 0; k < 3; k++) {
          spline_layer[k] = utils::inumeric(FLERR, v[k], false, lmp);
          if (spline_layer[k] < 16) error->all(FLERR, "pair_style {}: spline_layers values must be >= 16", style);
        }
        break;
      case SPLINE_MB:
        spline_mb = aut ? -1.0 : utils::numeric(FLERR, v[0], false, lmp);
        if (!aut && spline_mb < 0.0) error->all(FLERR, "pair_style {}: spline_mb must be >= 0 or auto", style);
        break;
      case STORE_FWD: store_fwd = aut ? -1 : utils::logical(FLERR, v[0], false, lmp); break;
      case FC_FOLD: fc_fold = utils::logical(FLERR, v[0], false, lmp); break;
      case STORE_L1: store_l1 = utils::logical(FLERR, v[0], false, lmp); break;
      case BLOCK:
        block = utils::inumeric(FLERR, v[0], false, lmp);
        if (block < 1 || block > 64) error->all(FLERR, "pair_style {}: block must be in 1..64", style);
        break;
      case FUSE_RADIAL: fuse_radial = utils::logical(FLERR, v[0], false, lmp); break;
      case NOSHAPE: noshape = true; break;
      case CHUNKSIZE:
        chunksize = aut ? -1 : utils::inumeric(FLERR, v[0], false, lmp);
        if (!aut && chunksize < 0)
          error->all(FLERR, "pair_style {}: chunksize must be auto or >= 0 (0: one chunk)", style);
        break;
      case DEBUG_NO_ENERGY_ONLY_CALC: debug_no_energy_only_calc = true; break;
    }
  }
}

void GraceRuntimeSettings::check_kind(LAMMPS *lmp, const std::string &style, grace::ModelKind kind,
                                      const std::string &path) const
{
  for (const Keyword *kw : given)
    if (!(kw->kinds & grace_kind_bit(kind))) {
      std::string applies;
      for (grace::ModelKind k : kinds_all)
        if (kw->kinds & grace_kind_bit(k))
          applies += (applies.empty() ? "" : ", ") + std::string(grace::model_kind_name(k));
      lmp->error->all(FLERR, "pair_style {}: keyword {} does not apply to the {} model {} (it applies to {})", style,
                      kw->name, grace::model_kind_name(kind), path, applies);
    }
  if (spline == 0 && kind != grace::ModelKind::Pace)
    lmp->error->all(FLERR, "pair_style {}: spline 0 (no radial tables) is for PACE models only", style);
}

namespace LAMMPS_NS {
namespace GraceRuntime {

ModelFile model_file(LAMMPS *lmp, const std::string &style, bool pace, const GraceRuntimeSettings &set,
                            int narg, char **arg)
{
  Error *error = lmp->error;
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients");
  ModelFile m{grace::ModelKind::Grace1L, arg[2], "", 0};
  if (pace && narg > 3 && utils::strmatch(arg[3], "\\.asi$")) {    // as pair_style pace/extrapolation
    m.asi = arg[3];
    m.nskip = 1;
  }
  const std::string sfx = style.substr(style.find("/runtime") + 8);    // "" or "/kk"
  if (pace && utils::strmatch(m.path, "\\.(yace|ace)$"))
    error->all(FLERR, "pair_style {} reads the B-basis .yaml of pacemaker, not {} (or use pair_style pace{})", style,
               m.path, sfx);
  try {
    m.kind = grace::model_kind(m.path);
  } catch (const grace::Error &e) {
    error->all(FLERR, "pair_style {}: {}", style, e.what());
  }
  if ((m.kind == grace::ModelKind::Pace) != pace)
    error->all(FLERR, "pair_style {}: {} is a {} model; use pair_style {}/runtime{}", style, m.path,
               grace::model_kind_name(m.kind), pace ? "grace" : "pace", sfx);
  set.check_kind(lmp, style, m.kind, m.path);
  if (lmp->comm->me == 0)
    utils::logmesg(lmp, "[{}] loading the {} model {}{}\n", style, grace::model_kind_name(m.kind), m.path,
                   m.asi.empty() ? "" : " with the active set " + m.asi);
  return m;
}

std::vector<int> type2mu(LAMMPS *lmp, const std::string &style, const int *map, char *const *elements,
                                const std::vector<std::string> &names)
{
  const int ntypes = lmp->atom->ntypes;
  std::vector<int> t2m(ntypes + 1, -1);
  for (int t = 1; t <= ntypes; t++) {
    if (map[t] < 0) continue;
    const auto it = std::find(names.begin(), names.end(), elements[map[t]]);
    if (it == names.end())
      lmp->error->all(FLERR, "pair_style {}: element {} is not in the model", style, elements[map[t]]);
    t2m[t] = (int) (it - names.begin());
  }
  return t2m;
}

void log(void *user, int level, const char *msg)
{
  auto *lmp = static_cast<LAMMPS *>(user);
  if (lmp->comm->me != 0) return;
  std::string m(msg);
  while (!m.empty() && m.back() == '\n') m.pop_back();
  if (level) lmp->error->warning(FLERR, "{}", m);
  else utils::logmesg(lmp, m + "\n");
}

int ranks_sharing(MPI_Comm world, const std::string &key)
{
#ifdef MPI_STUBS
  (void) world;
  (void) key;
  return 1;
#else
  MPI_Comm node;
  MPI_Comm_split_type(world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
  int nnode = 1, len = (int) key.size(), maxlen = 0;
  MPI_Comm_size(node, &nnode);
  if (key.empty()) {
    MPI_Comm_free(&node);
    return nnode;
  }
  MPI_Allreduce(&len, &maxlen, 1, MPI_INT, MPI_MAX, node);
  std::vector<char> mine(maxlen + 1, '\0'), all((size_t) nnode * (maxlen + 1), '\0');
  std::copy(key.begin(), key.end(), mine.begin());
  MPI_Allgather(mine.data(), maxlen + 1, MPI_CHAR, all.data(), maxlen + 1, MPI_CHAR, node);
  MPI_Comm_free(&node);
  int n = 0;
  for (int r = 0; r < nnode; r++)
    if (key == std::string(all.data() + (size_t) r * (maxlen + 1))) n++;
  return n;
#endif
}

}    // namespace GraceRuntime
}    // namespace LAMMPS_NS
#endif    // NO_GRACE_RUNTIME
