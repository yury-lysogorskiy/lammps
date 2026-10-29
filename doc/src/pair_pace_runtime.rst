.. index:: pair_style pace/runtime
.. index:: pair_style pace/runtime/kk

pair_style pace/runtime command
===============================

Accelerator Variants: *pace/runtime/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style pace/runtime keyword value ...
   pair_style pace/runtime/kk keyword value ...

* zero or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *spline* or *chunksize* or *debug_no_energy_only_calc*
       *spline* value = N
         N = intervals of the cubic radial tables per element pair (default 1024, at least 16),
             0 = radial functions evaluated directly
       *chunksize* value = *auto* or N
         auto = atoms per pass from the free GPU memory (default)
         N = atoms per pass, 0 = all atoms in one pass
         (GPU; ignored with a warning on the CPU)
       *debug_no_energy_only_calc* (GPU)
         run the force kernels also in energy-only requests (testing)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pace/runtime
   pair_coeff * * output_potential.yaml C

   pair_style pace/runtime/kk chunksize 20000
   pair_coeff * * output_potential.yaml H O

   pair_style pace/runtime/kk
   pair_coeff * * output_potential.yaml output_potential.asi C

Description
"""""""""""

Pair styles *pace/runtime* and *pace/runtime/kk* evaluate the same Atomic
Cluster Expansion (ACE) models as :doc:`pair_style pace <pair_pace>`
(C-tilde basis, Finnis-Sinclair embedding, inner cutoff by density or
ZBL core repulsion), with the kernels of the grace-runtime library
(libgrace).  The model is the B-basis ``.yaml`` written by pacemaker
(``output_potential.yaml``), converted to its C-tilde basis when it is
read.
The ``.yace`` / ``.ace`` formats are not read.  A GRACE model is an error that
names :doc:`pair_style grace/runtime <pair_grace_runtime>`.

*pace/runtime* runs on the CPU with OpenMP threads over groups of eight
atoms of one element, vectorized across neighbors and across the atoms
of a group; *pace/runtime/kk* runs on a GPU (Kokkos CUDA or HIP).  Both
compute in double precision.  The radial functions are tabulated as
cubic Hermite splines on two segments (from 0.1 to 0.25 of the pair
cutoff with half of the *spline* intervals, above it with *spline*
intervals) and evaluated directly below a tenth of the cutoff; with
*spline* 0 they are evaluated directly everywhere.  Forces then agree with
those of *pace* to about 1e-5 eV/Angstrom (1e-7 with *spline* 0).

Extrapolation grade
"""""""""""""""""""

With an active set, the ``.asi`` file of pacemaker after the ``.yaml``
as for :doc:`pair_style pace/extrapolation <pair_pace>`, both styles
compute the extrapolation grade *gamma* of *pace/extrapolation* (linear
or non-linear, from the shape of the active set) on the steps :doc:`fix pair <fix_pair>` requests it:

.. code-block:: LAMMPS

   pair_style pace/runtime
   pair_coeff * * output_potential.yaml output_potential.asi C
   fix pace_gamma all pair 10 pace/runtime gamma 1
   compute max_pace_gamma all reduce max f_pace_gamma

With *-sf kk* the unsuffixed name in the fix refers to *pace/runtime/kk*.  Only
*gamma* is available (no *atomic_sigma* or *gmm_cluster*).

Only a single pair_coeff command is used, with the model file (and the
active set, if any) followed by the element of each LAMMPS atom type;
NULL leaves a type out (for :doc:`pair_style hybrid <pair_hybrid>`).  A *pair_style* command with the same style after
*pair_coeff* rebuilds the evaluator from that model file with the new
keywords; *pair_coeff* need not be repeated (under
:doc:`pair_style hybrid <pair_hybrid>` it must, as the sub-styles are
created anew).

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These pair styles do not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

These pair styles do not write their information to :doc:`binary restart
files <restart>`, since it is stored in potential files.  Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

----------

Restrictions
""""""""""""

These pair styles are part of the ML-PACE package, built with the
grace-runtime styles (CMake option ``GRACE_RUNTIME``, on by default).
They require metal units, atom IDs and newton pair on.
*pace/runtime/kk* needs a GPU Kokkos backend.

Related commands
""""""""""""""""

:doc:`pair_style pace <pair_pace>`,
:doc:`pair_style grace/runtime <pair_grace_runtime>`

Default
"""""""

spline = 1024, chunksize = auto
