.. index:: pair_style grace
.. index:: pair_style grace/1layer/chunk
.. index:: pair_style grace/2layer/chunk
.. index:: pair_style grace/2layer/parallel
.. index:: pair_style grace/fs
.. index:: pair_style grace/fs/kk
.. index:: pair_style grace/1l
.. index:: pair_style grace/1l/kk
.. index:: pair_style grace/1l/cpu
.. index:: pair_style grace/1l/cpu/kk
.. index:: pair_style grace/2l
.. index:: pair_style grace/2l/kk
.. index:: pair_style grace/2l/cpu
.. index:: pair_style grace/2l/cpu/kk
.. index:: pair_style grace/3l/kk

pair_style grace command
========================

.. versionadded:: 03Dec2025

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grace keyword value ...
   pair_style grace/1layer/chunk keyword value ...
   pair_style grace/2layer/chunk keyword value ...
   pair_style grace/2layer/parallel keyword value ...
   pair_style grace/fs keyword value ...
   pair_style grace/fs/kk keyword value ...
   pair_style grace/1l/kk keyword value ...
   pair_style grace/1l/cpu/kk keyword value ...
   pair_style grace/2l/kk keyword value ...
   pair_style grace/2l/cpu/kk keyword value ...
   pair_style grace/3l/kk keyword value ...

Zero or more keyword/value pairs may be appended to the ``pair_style``
command.  The accepted keywords depend on the selected style.

.. list-table:: Pair style keywords
   :header-rows: 1
   :widths: 28 52 20

   * - Style
     - Accepted keywords
     - Notes
   * - ``grace``
     - ``padding``, ``pad_verbose``, ``pair_forces``,
       ``max_number_of_reduction``, ``reduce_padding``,
       ``debug_no_energy_only_calc``, ``kappa``, ``bias_virial``,
       ``kappa_norm``, ``kappa_group``
     - TensorFlow saved-model evaluator; uncertainty quantification (UQ)
       activates automatically when the model exports a UQ head
   * - ``grace/1layer/chunk``
     - ``padding``, ``pad_verbose``, ``max_number_of_reduction``,
       ``reduce_padding``, ``chunksize``, ``debug_no_energy_only_calc``
     - Chunked single-layer TensorFlow evaluator
   * - ``grace/2layer/chunk``
     - ``padding``, ``pad_verbose``, ``max_number_of_reduction``,
       ``reduce_padding``, ``chunksize``, ``debug_no_energy_only_calc``
     - Chunked two-layer TensorFlow evaluator
   * - ``grace/2layer/parallel``
     - ``padding``, ``pad_verbose``, ``max_number_of_reduction``,
       ``reduce_padding``, ``debug_no_energy_only_calc``
     - Two-layer MPI TensorFlow evaluator without chunking
   * - ``grace/fs`` and ``grace/fs/kk``
     - ``extrapolation``, ``chunksize``, ``debug_no_energy_only_calc``
     - Native GRACE/FS evaluator
   * - ``grace/{1,2}l/kk``, ``grace/{1,2}l/cpu/kk``, and ``grace/3l/kk``
     - ``chunksize``, ``debug_no_energy_only_calc``
     - Native Kokkos GRACE 1L/2L/3L evaluator (gracemaker ``.npz`` export)

The common TensorFlow keywords are:

* ``padding`` value = fraction of neighbors to pad.  The default is 0.01.
* ``pad_verbose`` = print messages when a new padding level triggers
  recompilation.
* ``max_number_of_reduction`` value = maximum number of recompilations
  during padding reduction.
* ``reduce_padding`` value = fraction by which to reduce the padding buffer.
* ``debug_no_energy_only_calc`` = disable the optimized energy-only code path.
  The optimized path is activated automatically whenever the caller of
  ``Pair::compute`` passes ``ENERGY_ONLY`` in ``eflag`` (the standard LAMMPS
  signal for "compute energy only, no forces or virials").  Callers that do
  this include the Monte Carlo fixes (:doc:`fix sgcmc <fix_sgcmc>`,
  :doc:`fix atom/swap <fix_atom_swap>`,
  :doc:`fix neighbor/swap <fix_neighbor_swap>`, :doc:`fix gcmc <fix_gcmc>`,
  :doc:`fix widom <fix_widom>`), :doc:`compute fep <compute_fep>`,
  :doc:`fix numdiff <fix_numdiff>`, and the DIELECTRIC polarize fixes.
  Setting ``debug_no_energy_only_calc`` forces the full force/virial kernel
  even when ``ENERGY_ONLY`` is set; it is intended for testing and debugging.

The ``grace`` style additionally accepts:

* ``pair_forces`` = compute pairwise forces.  This is required for virials
  and stress with ``grace`` and is enabled automatically when running on more
  than one MPI rank.

The chunked TensorFlow styles additionally accept:

* ``chunksize`` value = number of atoms processed in each TensorFlow block.
  The default is 4096.

For uncertainty quantification and biased dynamics, the ``grace`` style
additionally accepts the following keywords.  They are only meaningful when the
saved model exports a UQ head (see :ref:`pair_grace_uq`); on a non-UQ model they
configure nothing and UQ output is never produced.

* ``kappa`` value = relative-force bias coefficient :math:`\kappa`.  The
  default is 0, which disables biased dynamics.
* ``bias_virial`` = include the :math:`\kappa` bias contribution in the
  global virial.  Per-atom stress is unaffected.
* ``kappa_norm`` value = ``max`` or ``mean``.  This selects the norm used to
  rescale the uncertainty force in biased dynamics.  The default is ``max``.
* ``kappa_group`` value = LAMMPS atom group name restricting the
  :math:`\kappa` bias.  The default is ``all``.  See
  :ref:`pair_grace_bias_dynamics` for the precise semantics.

The ``grace/fs`` and ``grace/fs/kk`` styles additionally accept:

* ``extrapolation`` = compute the MaxVol extrapolation grade.  This requires
  an Active Set Inverted file in the ``pair_coeff`` command.
* ``chunksize`` value = number of atoms processed in each evaluator block.
  The default is 4096.

Pair coefficients
"""""""""""""""""

Each style documented on this page uses a single ``pair_coeff`` command with
``* *`` followed by a model file or directory and one element name per LAMMPS
atom type.

For TensorFlow GRACE styles:

.. code-block:: LAMMPS

   pair_coeff * * saved_model_dir elem1 elem2 ...

For ``grace/fs`` and ``grace/fs/kk`` without extrapolation:

.. code-block:: LAMMPS

   pair_coeff * * model.yaml elem1 elem2 ...

For ``grace/fs`` and ``grace/fs/kk`` with extrapolation:

.. code-block:: LAMMPS

   pair_coeff * * model.yaml model.asi elem1 elem2 ...

The number of element names must equal the number of LAMMPS atom types.  The
element names define the mapping from LAMMPS atom types to model elements.

Examples
""""""""

Basic TensorFlow GRACE model:

.. code-block:: LAMMPS

   pair_style grace
   pair_coeff * * /path/to/saved_model Al Li

TensorFlow GRACE model with explicit padding and pairwise forces:

.. code-block:: LAMMPS

   pair_style grace padding 0.05 pad_verbose pair_forces
   pair_coeff * * /path/to/grace_cache/AlLi_model Al Li

Chunked single-layer TensorFlow model:

.. code-block:: LAMMPS

   pair_style grace/1layer/chunk chunksize 2048 padding 0.05
   pair_coeff * * /path/to/saved_model Al Li

Chunked two-layer TensorFlow model:

.. code-block:: LAMMPS

   pair_style grace/2layer/chunk
   pair_coeff * * /path/to/2layer_model Al Li

Two-layer TensorFlow model with the non-chunked MPI-parallel implementation:

.. code-block:: LAMMPS

   pair_style grace/2layer/parallel
   pair_coeff * * /path/to/2layer_model Al Li

UQ-enabled TensorFlow model (UQ output activates via :doc:`fix pair <fix_pair>`):

.. code-block:: LAMMPS

   pair_style grace
   pair_coeff * * /path/to/uq_saved_model W

UQ-enabled TensorFlow model with biased dynamics:

.. code-block:: LAMMPS

   pair_style grace kappa 0.1
   pair_coeff * * /path/to/uq_saved_model W

Native GRACE/FS model:

.. code-block:: LAMMPS

   pair_style grace/fs
   pair_coeff * * FS_model.yaml Mo Nb Ta W

Native GRACE/FS model with MaxVol extrapolation grade:

.. code-block:: LAMMPS

   pair_style grace/fs extrapolation
   pair_coeff * * FS_model.yaml FS_model.asi Mo Nb Ta W

Kokkos-accelerated GRACE/FS model:

.. code-block:: LAMMPS

   pair_style grace/fs/kk
   pair_coeff * * FS_model.yaml Mo Nb Ta W

Native Kokkos GRACE 2L model on GPU (FP32, lower memory and faster):

.. code-block:: LAMMPS

   pair_style grace/2l/kk/fp32
   pair_coeff * * grace_2l_weights.npz Mo Nb Ta W

Native Kokkos GRACE 3L model on GPU (Mixed precision, the default):

.. code-block:: LAMMPS

   pair_style grace/3l/kk
   pair_coeff * * grace_3l_weights.npz Cu

Native Kokkos GRACE 1L model on CPU with custom chunk size:

.. code-block:: LAMMPS

   pair_style grace/1l/cpu/kk chunksize 2048
   pair_coeff * * grace_1l_weights.npz Cu

Description
"""""""""""

The ``grace`` pair styles compute interactions using the Graph Atomic Cluster
Expansion (GRACE) framework :ref:`(Bochkarev24) <Bochkarev20241>`,
:ref:`(Lysogorskiy25) <Lysogorskiy20251>`.

This page documents the following related pair styles:

* ``grace``
* ``grace/1layer/chunk``
* ``grace/2layer/chunk``
* ``grace/2layer/parallel``
* ``grace/fs``
* ``grace/fs/kk``
* ``grace/1l/kk`` and ``grace/1l/cpu/kk``
* ``grace/2l/kk`` and ``grace/2l/cpu/kk``
* ``grace/3l/kk``

Choosing a GRACE pair style
"""""""""""""""""""""""""""

.. list-table:: GRACE pair style variants
   :header-rows: 1
   :widths: 20 22 26 12 10 10

   * - Style
     - Model format
     - Main use case
     - MPI
     - Chunking
     - TensorFlow
   * - ``grace``
     - TensorFlow saved model
     - General model with simple setup; uncertainty quantification and biased
       dynamics when the model has a UQ head (biased dynamics: one rank only)
     - Yes for 1L, no for 2L
     - No
     - Required
   * - ``grace/1layer/chunk``
     - TensorFlow saved model
     - Single-layer model with lower peak memory use
     - Yes
     - Yes
     - Required
   * - ``grace/2layer/chunk``
     - TensorFlow saved model
     - Two-layer model with chunked evaluation
     - Yes
     - Yes
     - Required
   * - ``grace/2layer/parallel``
     - TensorFlow saved model
     - Two-layer model without chunking
     - Yes
     - No
     - Required
   * - ``grace/fs``
     - GRACE/FS YAML model
     - Native CPU evaluator
     - Yes
     - Yes
     - Not required
   * - ``grace/fs/kk``
     - GRACE/FS YAML model
     - Kokkos-accelerated native evaluator
     - Yes
     - Yes
     - Not required
   * - ``grace/1l/kk``
     - Native GRACE 1L ``.npz`` (gracemaker export)
     - Kokkos GPU evaluator for single-layer models
     - Yes
     - Yes
     - Not required
   * - ``grace/1l/cpu/kk``
     - Native GRACE 1L ``.npz`` (gracemaker export)
     - Kokkos CPU evaluator for single-layer models
     - Yes
     - Yes
     - Not required
   * - ``grace/2l/kk``
     - Native GRACE 2L ``.npz`` (gracemaker export)
     - Kokkos GPU evaluator for two-layer models
     - Yes
     - Yes
     - Not required
   * - ``grace/2l/cpu/kk``
     - Native GRACE 2L ``.npz`` (gracemaker export)
     - Kokkos CPU evaluator for two-layer models
     - Yes
     - Yes
     - Not required
   * - ``grace/3l/kk``
     - Native GRACE 3L ``.npz`` (gracemaker export)
     - Kokkos GPU evaluator for three-layer models
     - Yes
     - Yes
     - Not required

TensorFlow GRACE models
"""""""""""""""""""""""

The ``grace`` style uses ``libtensorflow`` to load and execute GRACE models
stored in TensorFlow ``saved_model`` format.  


GRACE TensorFlow models are just-in-time (JIT) compiled.  The first evaluation
is therefore slower than subsequent evaluations.  To maintain performance when
the number of neighbors changes, the TensorFlow-based styles use padding.  A
larger ``padding`` value can reduce the frequency of recompilations, but
increases memory use and the cost of each TensorFlow call.

The plain ``grace`` style does not chunk atoms.  To control peak memory through
chunked evaluation, use ``grace/1layer/chunk`` or ``grace/2layer/chunk``.
Both styles accept the ``chunksize`` keyword.

Chunking and MPI parallelization
""""""""""""""""""""""""""""""""

The ``grace/1layer/chunk`` style is a chunked variant of ``grace`` for
single-layer TensorFlow models.  Instead of processing all atoms at once, atoms
are processed in blocks of ``chunksize``.  This reduces peak memory use and is
safe for MPI parallelization without the ``pair_forces`` keyword.  Virial and
stress computation are always available.

The ``grace/2layer/chunk`` style supports two-layer TensorFlow GRACE models
using a chunked processing strategy.  Two-layer models split the computation
into a forward pass, which evaluates descriptors, and a backward pass, which
evaluates energies and forces.  Between the two passes, MPI forward and reverse
communication exchange per-atom features between ranks.  The ``chunksize``
keyword controls the block size used for each layer.

The ``grace/2layer/parallel`` style is an alternative MPI-parallel
implementation for two-layer TensorFlow GRACE models.  Each MPI rank performs
the full forward and backward layer evaluation locally, with explicit forward
and reverse communication of per-atom features between layers.  Unlike
``grace/2layer/chunk``, this style does not chunk TensorFlow calls and does not
accept the ``chunksize`` keyword.

For ``grace/2layer/chunk`` and ``grace/2layer/parallel``, the saved-model
directory must contain the TensorFlow signatures required by the two-layer
model, including ``forward_layer_1`` and ``backward_layer_2``.

GRACE/FS models
"""""""""""""""

The ``grace/fs`` style provides a native C++ implementation, the product
evaluator, for the GRACE/FS family of models.  It does not require TensorFlow
or GPUs and supports MPI parallelization.

The ``grace/fs/kk`` style is the Kokkos-accelerated version of ``grace/fs``.
It supports GPU execution through Kokkos backends such as CUDA or HIP and
multicore CPU execution through OpenMP.  Newton's third law must be enabled
with ``newton on`` and only one CPU thread per MPI rank is supported.  The
device and host variants are also available as ``grace/fs/kk/device`` and
``grace/fs/kk/host``.

When the ``extrapolation`` keyword is used, ``grace/fs`` and ``grace/fs/kk``
compute the MaxVol extrapolation grade :math:`\gamma` using the Active Set
Inverted (ASI) file supplied in the ``pair_coeff`` command.  This requires a
matrix-vector multiplication per atom and is slower than plain evaluation.

Native Kokkos GRACE 1L/2L models
""""""""""""""""""""""""""""""""

The ``grace/1l/kk`` and ``grace/2l/kk`` styles (and their ``/cpu/kk``
companions) are native Kokkos implementations of single-layer and two-layer
GRACE models.  Unlike the TensorFlow-based ``grace`` family, these styles
read weights from a LAMMPS-specific ``.npz`` archive produced by the
gracemaker export utilities; arbitrary TensorFlow saved models are not
accepted.  See `gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_
for the export commands.

The base styles (``grace/1l/kk``, ``grace/2l/kk``) are optimized for GPU
execution; the ``/cpu/kk`` variants (``grace/1l/cpu/kk``,
``grace/2l/cpu/kk``) are optimized for multicore CPU execution.  All four
share the same ``.npz`` model format.  Newton's third law must be enabled
with ``newton on``.

All ``grace/{1,2}l{,/cpu}/kk`` styles support both atom chunking via the
``chunksize`` keyword (default 4096) and MPI spatial decomposition.  The
``grace/3l/kk`` styles also support ``chunksize``, with a default cap selected
from the CUDA device memory: 4096 atoms below 40 GiB, 8192 atoms from 40 GiB,
and 16384 atoms from 70 GiB.  This lets larger devices use the single-chunk
backward fast path while retaining the explicit override for
memory-constrained workloads.
Chunking caps the per-rank atom batch processed in a single kernel call to
bound peak memory and limit register pressure on GPUs; it is independent of
MPI decomposition and works with any rank count.

For each base style, ``/device`` and ``/host`` suffixes pin execution to the
GPU or CPU regardless of the build's default Kokkos device, and ``/mixed``
and ``/fp32`` select the floating-point precision used inside the kernels:

* (default) -- 64-bit forward + backward (full accuracy).
* ``/mixed`` -- 64-bit forward geometrical features, 32-bit trainable features.
* ``/fp32`` -- 32-bit forward + backward (lowest memory and fastest,
  reduced accuracy; suitable for screening or preview runs).

For example, ``grace/2l/kk/fp32/device`` runs the two-layer GPU evaluator in
fp32, while ``grace/2l/cpu/kk/mixed/host`` runs the CPU evaluator in mixed
precision.

Native Kokkos GRACE 3L model
""""""""""""""""""""""""""""

The ``grace/3l/kk`` style is a native Kokkos implementation of three-layer
GRACE models, following the same ``.npz`` model format and gracemaker export
workflow as ``grace/{1,2}l/kk``.  Unlike the 1L/2L styles, ``grace/3l/kk`` has
no ``/cpu`` or ``/host`` variant: it is a GPU-only evaluator.  It supports the
``chunksize`` keyword and MPI spatial decomposition in the same way as
``grace/{1,2}l/kk``.

``grace/3l/kk`` (alias ``grace/3l/kk/device``) is Mixed precision by
default and is the validated target precision: NN weights and activations
are 32-bit, while geometry, spherical harmonics, and accumulation remain
64-bit.  ``grace/3l/kk/fp32`` (alias ``grace/3l/kk/fp32/device``) runs
fully in 32-bit for lower memory and faster, lower-accuracy screening runs.
There is no separate full 64-bit or ``/mixed``-suffixed style for 3L.

As with the other native Kokkos GRACE styles, ``grace/3l/kk`` requires
Newton's third law and a half neighbor list, e.g. ``-pk kokkos newton on
neigh half`` on the command line (or ``package kokkos newton on neigh
half`` in the input script).

.. _pair_grace_uq:

Uncertainty quantification
""""""""""""""""""""""""""

.. versionadded:: 07May2026

The ``grace`` style computes uncertainty quantification (UQ) whenever the loaded
GRACE saved model exports a ``compute_uq`` TensorFlow signature.  The UQ head is
detected automatically at ``pair_coeff`` time -- no separate pair style and no
opt-in keyword are required.  If the model additionally exports a
``compute_uq_gamma_only`` signature, the style uses it automatically on the
gamma-only fast path (see *Signature dispatch* below).

UQ is computed only when a per-atom UQ field is requested at run time through
:doc:`fix pair <fix_pair>`, or when ``kappa != 0``.  Otherwise -- including every
plain ``grace`` run, even on a UQ-capable model -- the regular non-UQ code path
runs at no extra per-step cost.  If a UQ field (or ``kappa``) is requested but the
model has no ``compute_uq`` head, LAMMPS stops with an error.

The style can expose four per-atom fields through :doc:`fix pair <fix_pair>`:

* ``gamma`` = scalar extrapolation grade :math:`\gamma` (Mahalanobis-based;
  see below).
* ``gmm_cluster`` = scalar GMM cluster index :math:`k^*` for the atom
  (integer, returned as a per-atom double because *fix pair* accepts only
  doubles).
* ``atomic_sigma`` = scalar raw per-atom uncertainty :math:`\sigma_i` as
  produced by the UQ head of the saved model (eV).  This is the underlying
  ensemble/dropout disagreement signal that ``uncertainty_force`` and the
  kappa rescale are derived from; exposing it directly is useful for
  calibration plots, threshold tuning, and diagnostics.
* ``uncertainty_force`` = three-vector uncertainty force
  :math:`\mathbf{F}^{\sigma}_i = \partial \sigma_{tot} / \partial \mathbf{r}_i`.
  This is the raw uncertainty force and is not scaled by ``kappa``.


The ``uncertainty_force`` output is accumulated consistently with Newton's
third law when Newton pair communication is enabled.

**Signature dispatch.** The pair style selects the TensorFlow signature based
on which UQ outputs are requested and on the value of ``kappa``, in order of
increasing cost:

* Regular ``compute`` -- when no UQ field is requested and ``kappa = 0``.
* ``compute_uq_gamma_only`` -- when only ``gamma``, ``gmm_cluster``, and/or
  ``atomic_sigma`` are requested and ``kappa = 0``.  Skips the
  :math:`\sigma`-gradient backward pass.  Used automatically when the saved
  model exports this signature; otherwise falls back to ``compute_uq``.
* ``compute_uq`` -- when ``uncertainty_force`` is requested or ``kappa`` is
  nonzero.  Includes the :math:`\sigma`-gradient backward pass.

Example with the ``grace`` UQ head, biased dynamics, and thermo
output containing the maximum extrapolation grade, maximum total force, and
maximum raw uncertainty-force magnitude:

.. code-block:: LAMMPS

   pair_style  grace kappa 0.1
   pair_coeff  * * /path/to/uq_saved_model Mo Nb Ta W

   # Per-atom UQ outputs.
   # uncertainty_force is exposed as f_uf[1], f_uf[2], and f_uf[3].
   fix gp  all pair 1 grace gamma             1
   fix gc  all pair 1 grace gmm_cluster       1
   fix as  all pair 1 grace atomic_sigma      1
   fix uf  all pair 1 grace uncertainty_force 1

   # Magnitude of the raw uncertainty force.
   variable ufmag atom sqrt(f_uf[1]^2 + f_uf[2]^2 + f_uf[3]^2)

   compute max_gamma all reduce max f_gp
   compute fmax_uq   all reduce max v_ufmag

   thermo_style custom step pe c_max_gamma fmax c_fmax_uq
   thermo 50

   # Dump per-atom UQ fields. The cluster index is stored as a per-atom
   # double; "%.0f" prints it as a whole number.
   dump uq all custom 100 uq.dump id type x y z f_gp f_gc f_as
   dump_modify uq format line "%d %d %g %g %g %g %.0f %g"

   # Stop the run when the extrapolation grade exceeds its threshold.
   fix stop all halt 10 c_max_gamma > 1.0

In this example, ``fmax`` is the maximum magnitude of the LAMMPS-stored force,
which includes the :math:`\kappa`-scaled bias when ``kappa`` is nonzero.
``c_fmax_uq`` is the maximum raw uncertainty-force magnitude.  The
``fix halt`` command stops the run when ``c_max_gamma`` exceeds 1.0.

For UQ-enabled TensorFlow GRACE models, the extrapolation grade is computed
from a per-element Gaussian mixture model fitted to latent feature vectors
during model export.  The latent feature vector :math:`\mathbf{z}` is a fixed
random projection of the rotationally-invariant (:math:`l=0`) energy-path
basis functions.  The exported saved model contains the projection matrix, the
fitted GMM parameters, and the per-cluster calibration thresholds.

For each local atom :math:`i` of element :math:`e` with latent vector
:math:`\mathbf{z}_i`, the UQ head assigns the atom to the nearest cluster by
Euclidean distance to the per-element centroids
:math:`\boldsymbol{\mu}_{e,k}`:

.. math::

   k^{\!*} = \arg\min_k \| \mathbf{z}_i - \boldsymbol{\mu}_{e,k} \|^2 .

It then computes ``atomic_sigma`` as the Mahalanobis distance of
:math:`\mathbf{z}_i` to the assigned centroid in the metric induced by the
cluster covariance :math:`\Sigma_{e,k^{\!*}}`:

.. math::

   \sigma_i =
   \sqrt{
   (\mathbf{z}_i - \boldsymbol{\mu}_{e,k^{\!*}})^\top
   \Sigma_{e,k^{\!*}}^{-1}
   (\mathbf{z}_i - \boldsymbol{\mu}_{e,k^{\!*}})
   } .

The extrapolation grade :math:`\gamma_i` is ``atomic_sigma`` normalized by the
calibrated 99th-percentile threshold :math:`\theta_{e,k^{\!*}}` of the assigned
cluster:

.. math::

   \gamma_i = \frac{\sigma_i}{\theta_{e,k^{\!*}}} .

.. _pair_grace_uq_kokkos:

Native Kokkos GRACE UQ
""""""""""""""""""""""

.. versionadded:: TBD

The native Kokkos GRACE styles ``grace/1l/kk``, ``grace/2l/kk``,
``grace/1l/cpu/kk``, and ``grace/2l/cpu/kk`` (including their ``/device``,
``/host``, ``/mixed``, and ``/fp32`` precision variants) compute the same
per-atom uncertainty quantities on the GPU, directly from a UQ-enabled Kokkos
``.npz`` model file.  The pair style detects UQ artifacts automatically.

These styles expose three per-atom fields through :doc:`fix pair <fix_pair>`,
with the same names and semantics as the TensorFlow ``grace`` style above:

* ``gamma`` = extrapolation grade :math:`\gamma_i`,
* ``atomic_sigma`` = raw Mahalanobis uncertainty :math:`\sigma_i`,
* ``gmm_cluster`` = GMM cluster index :math:`k^*` for the atom (an integer
  stored as a double).

The latent feature :math:`\mathbf{z}_i` feeding the GMM is built identically to
the TensorFlow ``grace`` UQ head, so ``gamma`` matches it to floating-point
precision.  It is the concatenation of two parts:

* a **random projection** of the rotationally-invariant (:math:`l=0`)
  energy-path B-basis :math:`\mathbf{B}_i`, after L2-normalization, through the
  stored projection matrix :math:`R`:
  :math:`\mathbf{p}_i = (\mathbf{B}_i / \|\mathbf{B}_i\|)\,R`; and
* one or more **log-norm density channels**
  :math:`\log(\|\cdot\| + 10^{-12})\cdot s_\rho`, one for the full basis and one
  per invariant block, where :math:`s_\rho` is the stored ``density_scale``.

For the two-layer styles the basis has two invariant blocks; the density
channels are ordered ``[full, block0, block1, ...]`` following the block order
baked into the artifacts (the projection-matrix row layout follows the same
order).

The Kokkos styles accept **only** UQ artifacts written with schema version 6
(the L2-normalized projection + density-channel feature described above; the
``.npz`` carries ``uq_schema_version = 6`` with ``uq_rp_normalize = 1``,
``uq_rp_add_density_channel = 1``, ``uq_feature_transform = 0``).  Older
artifacts are rejected at load time with a re-export message.  The artifacts
must also provide GMM clusters for **every** model element; an incomplete export
(an element with no clusters) is a hard load-time error rather than a silent
fallback.

Unlike the TensorFlow ``grace`` UQ path, the native Kokkos styles
provide **no** ``uncertainty_force``,  ``kappa``/HAL biasing, or
:math:`\partial\gamma/\partial\mathbf{r}` gradient: only the forward per-atom
scalars listed above are available.

.. code-block:: LAMMPS

   pair_style grace/1l/kk          # or grace/2l/kk, grace/1l/cpu/kk, grace/2l/cpu/kk
   pair_coeff * * grace_kokkos_uq.npz Mo Nb Ta W

   fix gp  all pair 1 grace/1l/kk gamma          1
   fix as  all pair 1 grace/1l/kk atomic_sigma   1
   fix gm  all pair 1 grace/1l/kk gmm_cluster    1

   compute max_gamma all reduce max f_gp
   thermo_style custom step pe c_max_gamma
   dump uq all custom 100 uq.dump id type x y z f_gp f_as f_gm

.. _pair_grace_bias_dynamics:

Bias-driven dynamics (experimental)
"""""""""""""""""""""""""""""""""""

.. warning::

   Bias-driven dynamics is experimental features and may be changed in the future.

.. warning::

   Bias-driven dynamics with ``kappa != 0`` is supported only on a single MPI
   rank.  LAMMPS stops with an error during initialization if ``kappa`` is
   nonzero and more than one MPI rank is used.

The ``kappa`` keyword sets :math:`\kappa \geq 0`.  When nonzero, it adds a
relative-force uncertainty bias on top of the physical force:

.. math::

   \mathbf{f}_i = \mathbf{F}^{phys}_i + s\,\mathbf{F}^{\sigma}_i,
   \qquad
   s = \kappa \cdot
       \frac{N(\|\mathbf{F}^{phys}_j\|) + \varepsilon}
            {N(\|\mathbf{F}^{\sigma}_j\|) + \varepsilon}

with :math:`\varepsilon = 10^{-8}`.  The scale factor :math:`s` is one
system-wide scalar applied uniformly to every atom.  The aggregator
:math:`N(\cdot)` over all atoms :math:`j` is selected by ``kappa_norm``:

* ``max`` = use the maximum force norm over atoms.  This is the default.
* ``mean`` = use the sum of force norms over atoms.

By default, the :math:`\kappa` bias is applied to every atom and the
aggregator :math:`N(\cdot)` runs over every atom.  When ``kappa_group`` names a
LAMMPS atom group, both the reduction and the force application are restricted
to atoms in that group: :math:`s` is built from the in-group atoms only and
:math:`s\,\mathbf{F}^{\sigma}_i` is added only to in-group atoms.  Atoms
outside the group experience the physical force only.  ``gamma``,
``gmm_cluster``, and ``uncertainty_force`` are still exposed for every local
atom regardless of ``kappa_group``.

By default, the global virial and per-atom stress reflect only the physical
force.  With ``bias_virial``, the uncertainty contribution is added to the
global virial.  Per-atom stress is built from the physical force only in either
case.

.. _pair_grace_runtime_tuning:

Runtime tuning
""""""""""""""

The following ``grace`` UQ / biased-dynamics parameters can be changed between
``run`` commands with :doc:`pair_modify <pair_modify>`:

.. code-block:: LAMMPS

   pair_modify kappa 0.5
   pair_modify bias_virial yes
   pair_modify kappa_norm mean
   pair_modify kappa_group boundary

The current value of ``kappa`` can be accessed from the pair style through the
LAMMPS ``extract`` interface using the key ``"kappa"``.

All other ``pair_modify`` keywords, such as ``compute`` and ``special``, are
forwarded to the standard pair-style parser unchanged.

GRACE/FS extrapolation example
""""""""""""""""""""""""""""""

The following input fragment computes the GRACE/FS MaxVol extrapolation grade
every 100 steps, stores it in a per-atom ``fix pair`` output, dumps only
structures whose maximum extrapolation grade exceeds 5, and stops when the
maximum grade exceeds 25.

.. code-block:: LAMMPS

   pair_style  grace/fs extrapolation
   pair_coeff  * * FS_model.yaml FS_model.asi Mo Nb Ta W

   fix grace_gamma all pair 100 grace/fs gamma 1

   compute max_grace_gamma all reduce max f_grace_gamma
   variable dump_skip equal "c_max_grace_gamma < 5"

   dump grace_dump all custom 20 extrapolative_structures.dump id type x y z f_grace_gamma
   dump_modify grace_dump skip v_dump_skip

   variable max_grace_gamma equal c_max_grace_gamma
   fix extreme_extrapolation all halt 10 v_max_grace_gamma > 25

Energy-only calculation
"""""""""""""""""""""""

All styles support an optimized energy-only calculation mode.  In this mode,
the styles skip atomic forces and virials, which reduces computational cost.
For TensorFlow-based GRACE models, this also avoids the backward gradient pass.

The energy-only mode is selected automatically when LAMMPS requests only the
potential energy from the pair style.  This is commonly used by Monte Carlo
algorithms implemented in the MC package.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Mixing is not used by these pair styles.  The element mapping and all model
parameters are specified by the model file or directory and the ``pair_coeff``
command.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
``shift``, ``table``, and ``tail`` options.

The ``grace`` style additionally accepts four
biased-dynamics-specific :doc:`pair_modify <pair_modify>` keywords:

* ``pair_modify kappa value`` = reset the relative-force biased-dynamics
  coefficient :math:`\kappa`.  A nonzero value is rejected with an error when
  more than one MPI rank is used.
* ``pair_modify bias_virial yes|no`` = toggle whether the :math:`\kappa`
  contribution is included in the global virial.  Per-atom stress is
  unaffected.
* ``pair_modify kappa_norm max|mean`` = switch the norm used in the
  :math:`\kappa` rescaling.
* ``pair_modify kappa_group name`` = restrict the :math:`\kappa` bias to
  atoms in the named LAMMPS group.  Pass ``all`` to remove the restriction.

These settings can be changed between ``run`` commands.

These pair styles do not write their information to :doc:`binary restart files
<restart>`, since the information is stored in model files or directories.  The
``pair_style`` and ``pair_coeff`` commands must therefore be specified again in
an input script that reads a restart file.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

* These styles are part of the ML-PACE package.  They are enabled only if
  LAMMPS is built with that package.
* All styles require :doc:`units metal <units>`.
* The TensorFlow-based styles require LAMMPS to be linked with
  ``libtensorflow``: ``grace``, ``grace/1layer/chunk``,
  ``grace/2layer/chunk``, and ``grace/2layer/parallel``.
* Uncertainty quantification output and biased dynamics in ``grace`` require a
  UQ-enabled saved model that exports the ``compute_uq`` signature.
* Biased dynamics with ``kappa != 0`` in ``grace`` is supported
  only on a single MPI rank.
* ``grace/fs`` does not require TensorFlow.
* ``grace/*/kk`` requires ``newton on`` and supports only one CPU thread per
  MPI rank.

Further reading
"""""""""""""""

See `gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_ for more
details.

Related commands
""""""""""""""""

:doc:`pair_style pace <pair_pace>`, :doc:`fix pair <fix_pair>`

Default
"""""""

.. list-table:: Default settings
   :header-rows: 1
   :widths: 32 68

   * - Style
     - Defaults
   * - ``grace``
     - ``padding = 0.01``, ``pair_forces = on`` (use ``no_pair_forces`` to
       disable; forces are always on under MPI with more than one rank),
       ``pad_verbose = off``, ``max_number_of_reduction = 10``,
       ``reduce_padding = 0.2``, ``kappa = 0``, ``kappa_norm = max``,
       ``kappa_group = all``, ``bias_virial = off``
   * - ``grace/1layer/chunk`` and ``grace/2layer/chunk``
     - ``padding = 0.01``, ``chunksize = 4096``, ``pad_verbose = off``,
       ``max_number_of_reduction = 10``, ``reduce_padding = 0.2``
   * - ``grace/2layer/parallel``
     - ``padding = 0.01``, ``pad_verbose = off``,
       ``max_number_of_reduction = 10``, ``reduce_padding = 0.2``
   * - ``grace/fs`` and ``grace/fs/kk``
     - ``extrapolation = off``, ``chunksize = 4096``

----------

.. _Bochkarev20241:

**(Bochkarev24)** Bochkarev, Lysogorskiy, Drautz, Phys. Rev. X, 14, 021036
(2024).

.. _Lysogorskiy20251:

**(Lysogorskiy25)** Lysogorskiy, Bochkarev, Drautz, arXiv:2508.17936 (2025).
