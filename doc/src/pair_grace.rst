.. index:: pair_style grace
.. index:: pair_style grace/1layer/chunk
.. index:: pair_style grace/2layer/chunk
.. index:: pair_style grace/2layer/parallel
.. index:: pair_style grace/fs
.. index:: pair_style grace/fs/kk

pair_style grace command
========================

.. versionadded:: 03Dec2025

pair_style grace/1layer/chunk command
======================================

.. versionadded:: 03Dec2025

pair_style grace/2layer/chunk command
======================================

.. versionadded:: 03Dec2025

pair_style grace/2layer/parallel command
=========================================

.. versionadded:: 03Dec2025

pair_style grace/fs command
============================

.. versionadded:: 03Dec2025

pair_style grace/fs/kk command
================================

.. versionadded:: 03Dec2025

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grace keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *padding* or *pad_verbose* or *pair_forces* or *max_number_of_reduction* or *reduce_padding* or *debug_no_energy_only_calc*
       *padding* value = fraction of neighbors to pad (default 0.01)
       *pad_verbose* = print messages when padding triggers recompilation
       *pair_forces* = compute pairwise forces (required for virials and MPI > 1)
       *max_number_of_reduction* value = maximum number of recompilations during padding reduction
       *reduce_padding* value = fraction to reduce padding by
       *debug_no_energy_only_calc* = disable energy-only optimization for debugging

.. code-block:: LAMMPS

   pair_style grace/1layer/chunk keyword values ...
   pair_style grace/2layer/chunk keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *padding* or *pad_verbose* or *max_number_of_reduction* or *reduce_padding* or *chunksize* or *debug_no_energy_only_calc*
       *padding* value = fraction of neighbors to pad (default 0.01)
       *pad_verbose* = print messages when padding triggers recompilation
       *max_number_of_reduction* value = maximum number of recompilations during padding reduction
       *reduce_padding* value = fraction to reduce padding by
       *chunksize* value = size of atom blocks processed by TensorFlow (default 4096)
       *debug_no_energy_only_calc* = disable energy-only optimization for debugging

.. code-block:: LAMMPS

   pair_style grace/2layer/parallel keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *padding* or *pad_verbose* or *max_number_of_reduction* or *reduce_padding* or *debug_no_energy_only_calc*
       *padding* value = fraction of neighbors to pad (default 0.01)
       *pad_verbose* = print messages when padding triggers recompilation
       *max_number_of_reduction* value = maximum number of recompilations during padding reduction
       *reduce_padding* value = fraction to reduce padding by
       *debug_no_energy_only_calc* = disable energy-only optimization for debugging

.. code-block:: LAMMPS

   pair_style grace/fs keyword values ...
   pair_style grace/fs/kk keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *extrapolation* or *chunksize* or *debug_no_energy_only_calc*
       *extrapolation* = compute extrapolation grade (requires .asi file in pair_coeff)
       *chunksize* value = size of atom blocks processed by the evaluator (default 4096)
       *debug_no_energy_only_calc* = disable energy-only optimization for debugging

Examples
""""""""

.. code-block:: LAMMPS

   pair_style grace
   pair_coeff * * /path/to/saved_model Al Li

   pair_style grace padding 0.05 pad_verbose pair_forces
   pair_coeff * * /path/to/grace_cache/AlLi_model Al Li

   pair_style grace/1layer/chunk
   pair_coeff * * /path/to/saved_model Al Li

   pair_style grace/1layer/chunk chunksize 2048 padding 0.05
   pair_coeff * * /path/to/saved_model Al Li

   pair_style grace/2layer/chunk
   pair_coeff * * /path/to/2layer_model Al Li

   pair_style grace/2layer/parallel
   pair_coeff * * /path/to/2layer_model Al Li

   pair_style grace/fs
   pair_coeff * * FS_model.yaml Mo Nb Ta W

   pair_style grace/fs extrapolation
   pair_coeff * * FS_model.yaml FS_model.asi Mo Nb Ta W

   pair_style grace/fs/kk
   pair_coeff * * FS_model.yaml Mo Nb Ta W

Description
"""""""""""

The *grace*, *grace/fs*, and their variants compute interactions using the
Graph Atomic Cluster Expansions (GRACE) framework :ref:`(Bochkarev24) <Bochkarev20241>`,
:ref:`(Lysogorskiy25) <Lysogorskiy20251>`.

**pair_style grace**

The *grace* style utilizes *libtensorflow* to load and execute GRACE models
saved in TensorFlow ``saved_model`` format. It is designed for single-layer
models and processes all local atoms in a single TensorFlow call.

Only a single pair_coeff command is used with the *grace* style which
specifies the directory containing the saved model followed by N additional
arguments specifying the mapping of GRACE model elements to LAMMPS atom types,
where N is the number of LAMMPS atom types:

* path to ``saved_model`` directory
* N element names = mapping of model elements to atom types

GRACE models in TensorFlow are Just-In-Time (JIT) compiled. This means the
first evaluation will be slower than subsequent steps. To maintain performance
if the number of neighbors changes, the style uses a padding strategy.

* **padding**: Sets the fraction of neighbors to pad (default is 0.01 or 1%). Increasing this can reduce the frequency of recompilations but increases the time and memory overhead.
* **pad_verbose**: If specified, LAMMPS will output messages whenever new padding levels trigger a recompilation. By default this is false.
* **pair_forces**: By default, the GRACE model provides total atomic forces. If *pair_forces* is enabled, the model calculates pairwise forces. This is **required** for calculating atomic virials (stress) and is automatically enforced if running on more than one MPI processor.
* **max_number_of_reduction** and **reduce_padding**: Control the heuristics for reducing the padding buffer size dynamically during the simulation.
* **chunksize**: Controls the size of atom blocks processed by TensorFlow. This helps manage peak memory usage for large systems or models.

**pair_style grace/1layer/chunk**

The *grace/1layer/chunk* style is a chunked variant of *grace* for single-layer
TensorFlow models. Instead of processing all atoms at once, atoms are processed
in blocks of *chunksize*. This reduces peak memory usage and is safe for MPI
parallelization without requiring the *pair_forces* keyword (virial/stress
computation is always available).

Accepted keywords are the same as *grace* except *pair_forces* is not needed
and is not accepted.

**pair_style grace/2layer/chunk**

The *grace/2layer/chunk* style supports two-layer GRACE TensorFlow models
using a chunked processing strategy. Two-layer models split computation into
a forward pass (layer 1: descriptors) and a backward pass (layer 2: energy/forces).
Between the two passes, inter-processor communication of per-atom features is
performed using MPI forward and reverse communication. The *chunksize* keyword
controls the block size for each layer independently.

pair_coeff syntax:

* path to ``saved_model`` directory (must contain ``forward_layer_1`` and ``backward_layer_2`` signatures)
* N element names = mapping of model elements to atom types

**pair_style grace/2layer/parallel**

The *grace/2layer/parallel* style is an alternative MPI-parallel implementation
for two-layer GRACE TensorFlow models. Unlike *grace/2layer/chunk*, this style
does not use chunked processing (no *chunksize* keyword) and instead performs
the full layer evaluation per MPI rank with explicit forward and reverse
communication of per-atom features between layers.

pair_coeff syntax is the same as *grace/2layer/chunk*.

**pair_style grace/fs**

The *grace/fs* style provides a native C++ implementation (product evaluator)
for the GRACE/FS family of models. It is lightweight and does not require
TensorFlow or GPUs. MPI parallelization is natively supported.

Only a single pair_coeff command is used with the *grace/fs* style:

* GRACE/FS coefficient file (.yaml format)
* (Optional) Active Set Inverted file (.asi format) if *extrapolation* keyword is used
* N element names = mapping of elements to atom types

**pair_style grace/fs/kk**

The *grace/fs/kk* style is the Kokkos-accelerated version of *grace/fs*,
supporting execution on GPUs (via CUDA/HIP) and multi-core CPUs (via OpenMP)
while maintaining efficient MPI parallelization. It accepts the same keywords
as *grace/fs*. Newton must be on and only a single CPU thread per MPI rank
is supported. The device and host variants are also available as
*grace/fs/kk/device* and *grace/fs/kk/host*.

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade is implemented in `pair_style grace/fs`
via the *extrapolation* keyword. It is based on the MaxVol algorithm.
In order to compute the extrapolation grade one needs to provide:

#. GRACE/FS potential in `.yaml` format
#. Active Set Inverted (ASI) file for the corresponding potential (`.asi` format)

Calculation of extrapolation grades requires matrix-vector multiplication
for each atom and is slower than the standard evaluation. The extrapolation
grade is accessed via `fix pair`, which requests to compute `gamma`.

Example of monitoring extrapolation warnings:

.. code-block:: LAMMPS

    pair_style  grace/fs extrapolation
    pair_coeff  * * FS_model.yaml FS_model.asi Mo Nb Ta W

    # Compute gamma every 100 steps, store in f_grace_gamma
    fix grace_gamma all pair 100 grace/fs gamma 1

    compute max_grace_gamma all reduce max f_grace_gamma
    variable dump_skip equal "c_max_grace_gamma < 5"

    dump grace_dump all custom 20 extrapolative_structures.dump id type x y z f_grace_gamma
    dump_modify grace_dump skip v_dump_skip

    variable max_grace_gamma equal c_max_grace_gamma
    fix extreme_extrapolation all halt 10 v_max_grace_gamma > 25

Here extrapolation grade gamma is computed every 100 steps and is stored
in the `f_grace_gamma` per-atom variable. The largest value of extrapolation
grade among all atoms in a structure is reduced to the `c_max_grace_gamma`
variable. Only if this value exceeds extrapolation threshold 5 will the
structure be dumped.

Energy-only calculation
"""""""""""""""""""""""

All styles support an optimized energy-only calculation mode. In this mode,
the styles skip the computation of atomic forces and virials, which significantly
reduces the computational cost. For TensorFlow-based GRACE models, this also
allows skipping the backward (gradient) pass.

The energy-only mode is automatically triggered when LAMMPS requests only
potential energy from the pair style, which is frequently used by Monte Carlo
algorithms implemented in the MC package:

* :doc:`fix atom/swap <fix_atom_swap>`
* :doc:`fix neighbor/swap <fix_neighbor_swap>`
* :doc:`fix widom <fix_widom>`
* :doc:`fix gcmc <fix_gcmc>`
* :doc:`fix sgcmc <fix_sgcmc>`

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential files/directories. Thus,
you need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

----------

Restrictions
""""""""""""

These pair styles are part of the ML-PACE package. They are only enabled if
LAMMPS was built with that package.

All styles require `metal` units.

*pair_style grace*, *grace/1layer/chunk*, *grace/2layer/chunk*, and
*grace/2layer/parallel* rely on the TensorFlow library (via *libtensorflow*).
While GPU usage is optional, TensorFlow is significantly less efficient when
running solely on the CPU.

*pair_style grace* requires *pair_forces* (enabled automatically for MPI > 1)
to compute virials and stress. The *grace/1layer/chunk*, *grace/2layer/chunk*,
and *grace/2layer/parallel* styles always support virial/stress computation.

*pair_style grace/fs* does not require TensorFlow. Kokkos support is available
via *pair_style grace/fs/kk*. The *grace/fs/kk* style requires Newton's third
law to be on (``newton on``) and supports only a single CPU thread per MPI rank.

Further read
""""""""""""

See `gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_ for more details.

Related commands
""""""""""""""""

:doc:`pair_style pace  <pair_pace>`,
:doc:`fix pair  <fix_pair>`

Default
"""""""

For *grace*: padding = 0.01, pair_forces is OFF (auto-enabled for MPI > 1), pad_verbose is OFF, max_number_of_reduction = 10, reduce_padding = 0.2.

For *grace/1layer/chunk* and *grace/2layer/chunk*: padding = 0.01, pad_verbose is OFF, max_number_of_reduction = 10, reduce_padding = 0.2, chunksize = 4096.

For *grace/2layer/parallel*: padding = 0.01, pad_verbose is OFF, max_number_of_reduction = 10, reduce_padding = 0.2.

For *grace/fs* and *grace/fs/kk*: extrapolation is OFF, chunksize = 4096.

.. _Bochkarev20241:

**(Bochkarev24)** Bochkarev, Lysogorskiy, Drautz, Phys Rev X, 14, 021036 (2024).

.. _Lysogorskiy20251:

**(Lysogorskiy25)** Lysogorskiy, Bochkarev, Drautz, arXiv:2508.17936 (2025).
