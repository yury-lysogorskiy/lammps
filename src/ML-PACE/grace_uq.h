/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Yury Lysogorskiy (ICAMS)
   Per-atom UQ outputs (gamma, atomic_sigma, gmm_cluster) of the native GRACE styles for fix pair. The models' UQ
   heads are in libgrace (grace/uq_head.h, model_1l.h, model_2l.h).
------------------------------------------------------------------------- */

#ifndef LMP_GRACE_UQ_H
#define LMP_GRACE_UQ_H

#include <cstring>

namespace LAMMPS_NS {

// per-atom UQ outputs exposed to fix pair (extract / extract_peratom), shared by grace/runtime, pace/runtime and their
// KOKKOS variants; grow() keeps at least n entries (n = nlocal on the steps fix pair requests UQ)
struct GraceUQOut {
  int flag_gamma = 0, flag_sigma = 0, flag_cluster = 0;
  bool gamma_only = false;    // the extrapolation grade only (PACE): no atomic_sigma / gmm_cluster
  int nmax = 0;
  double *gamma = nullptr, *sigma = nullptr, *cluster = nullptr;
  GraceUQOut() = default;
  GraceUQOut(const GraceUQOut &) = delete;
  GraceUQOut &operator=(const GraceUQOut &) = delete;
  bool any(bool has_uq) const { return has_uq && (flag_gamma || flag_sigma || flag_cluster); }
  void *flag(const char *str)
  {
    if (strcmp(str, "gamma_flag") == 0) return &flag_gamma;
    if (gamma_only) return nullptr;
    if (strcmp(str, "atomic_sigma_flag") == 0) return &flag_sigma;
    if (strcmp(str, "gmm_cluster_flag") == 0) return &flag_cluster;
    return nullptr;
  }
  void *peratom(const char *str)
  {
    if (strcmp(str, "gamma") == 0) return gamma;
    if (gamma_only) return nullptr;
    if (strcmp(str, "atomic_sigma") == 0) return sigma;
    if (strcmp(str, "gmm_cluster") == 0) return cluster;
    return nullptr;
  }
  void grow(int n)
  {
    if (n <= nmax) return;
    release();
    nmax = n;
    gamma = new double[n];
    if (gamma_only) return;
    sigma = new double[n];
    cluster = new double[n];
  }
  void release()
  {
    delete[] gamma;
    delete[] sigma;
    delete[] cluster;
    gamma = sigma = cluster = nullptr;
    nmax = 0;
  }
  ~GraceUQOut() { release(); }
};

}    // namespace LAMMPS_NS

#endif
