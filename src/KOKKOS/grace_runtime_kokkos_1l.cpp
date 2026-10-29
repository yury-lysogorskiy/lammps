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
   Contributing author: Yury Lysogorskiy (ICAMS)
   The GRACE-1L Kokkos engine of grace/runtime/kk (pair_grace_runtime_kokkos.h), instantiated here so it compiles in
   parallel with the style
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#include "pair_grace_runtime_kokkos.h"

#include "grace/kk/grace_1l_fused_kk_impl.h"

GRACE_KK_INSTANTIATE_1L(LMPDeviceType, LAMMPS_NS::GraceLammpsTraits<LMPDeviceType>);
#endif    // NO_GRACE_RUNTIME
