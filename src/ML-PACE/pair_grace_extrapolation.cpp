// pair_style grace/extrapolation
//
// Thin subclass of PairGRACE that activates the saved-model UQ head
// (TF signature `compute_uq`) and exposes per-atom gamma plus an optional
// HAL-style kappa-rescaled mixing of -dsigma/dr into real forces. All UQ logic
// lives in the parent class, gated by `request_extrapolation`.

#ifndef NO_GRACE_TF

#include "pair_grace_extrapolation.h"

#include "comm.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

PairGRACEExtrapolation::PairGRACEExtrapolation(LAMMPS *lmp) : PairGRACE(lmp)
{
  request_extrapolation = true;
  // flag_compute_gamma / flag_compute_uncertainty_force stay at parent default 0;
  // `fix pair grace/extrapolation gamma|uncertainty_force ...` toggles them via
  // extract("<name>_flag") every Nevery steps. kappa != 0 also forces UQ every
  // step (handled in compute()). Pair-force mode is required so we have per-bond
  // sigma-gradients for the kappa-rescale and for the per-atom uncertainty_force tally.
  pair_forces = true;
}

void PairGRACEExtrapolation::settings(int narg, char **arg)
{
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE/extrapolation potentials require 'metal' units");

  // Pull out grace/extrapolation-specific keywords first, then forward the rest
  // to the parent settings() parser. We rebuild a filtered argv so the parent
  // sees only the keywords it knows.
  std::vector<char *> forwarded;
  forwarded.reserve(narg);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kappa") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR,
                   "[GRACE/extrapolation] kappa requires a numeric argument (relative-force "
                   "uncertainty bias coefficient; e.g. 0.1 = 10 % of ||F^phys||)");
      kappa = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE/extrapolation] kappa = {} (relative-force HAL bias)\n",
                       kappa);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bias_virial") == 0) {
      bias_virial = true;
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE/extrapolation] bias_virial: kappa contribution also tallied "
                            "into the global virial (per-atom stress unchanged)\n");
      iarg += 1;
    } else if (strcmp(arg[iarg], "kappa_norm") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "[GRACE/extrapolation] kappa_norm requires 'max' or 'mean'");
      parse_kappa_norm(arg[iarg + 1], "[GRACE/extrapolation] kappa_norm");
      iarg += 2;
    } else if (strcmp(arg[iarg], "kappa_group") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "[GRACE/extrapolation] kappa_group requires a LAMMPS group name");
      parse_kappa_group(arg[iarg + 1], "[GRACE/extrapolation] kappa_group");
      iarg += 2;
    } else if (strcmp(arg[iarg], "no_pair_forces") == 0) {
      error->all(FLERR,
                 "[GRACE/extrapolation] 'no_pair_forces' is incompatible: kappa-rescaling of "
                 "sigma-forces requires per-bond gradients (pair_forces mode).");
    } else {
      forwarded.push_back(arg[iarg]);
      iarg += 1;
    }
  }

  PairGRACE::settings(forwarded.size(), forwarded.data());

  // Re-assert the invariants after parent settings() runs
  pair_forces = true;
  if (comm->me == 0) {
    utils::logmesg(lmp, "[GRACE/extrapolation] UQ head enabled (compute_uq); gamma/sigma exposed via "
                        "extract_peratom\n");
    utils::logmesg(lmp, "[GRACE/extrapolation] kappa_norm = {} (HAL force-ratio normalization)\n",
                   kappa_norm_str());
    utils::logmesg(lmp, "[GRACE/extrapolation] kappa_group = {} (atoms eligible for HAL bias)\n",
                   kappa_group_str());
  }
}

#endif    //#ifndef NO_GRACE_TF
