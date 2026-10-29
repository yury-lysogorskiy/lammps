.. index:: pair_style grace
.. index:: pair_style grace/fs
.. index:: pair_style grace/fs/kk

pair_style grace command
========================

.. versionadded:: 03Dec2025

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grace keyword value ...
   pair_style grace/fs keyword value ...
   pair_style grace/fs/kk keyword value ...

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
       ``debug_no_energy_only_calc``, ``tau``, ``tau_norm``,
       ``tau_group``, ``stress_tau``, ``bulk_modulus``
     - TensorFlow saved-model evaluator; uncertainty quantification (UQ)
       activates automatically when the model exports a UQ head
   * - ``grace/fs`` and ``grace/fs/kk``
     - ``extrapolation``, ``chunksize``, ``debug_no_energy_only_calc``
     - Native GRACE/FS evaluator

The keywords of the TensorFlow ``grace`` style are
(``debug_no_energy_only_calc`` also applies to ``grace/fs`` and
``grace/fs/kk``):

* ``padding`` value = fraction of neighbors to pad.  The default is 0.01.
* ``pad_verbose`` = print messages when a new padding level triggers
  recompilation.
* ``pair_forces`` or ``no_pair_forces`` = compute pairwise forces, or not.
  They are needed for virials and stress, are on by default, and are always
  on with more than one MPI rank.
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

For uncertainty quantification and biased dynamics, the ``grace`` style
additionally accepts the following keywords.  They are only meaningful when the
saved model exports a UQ head (see :ref:`pair_grace_uq`); on a non-UQ model they
configure nothing and UQ output is never produced.

* ``tau`` value = relative-force bias coefficient :math:`\tau`.  The
  default is 0, which disables biased dynamics.
* ``tau_norm`` value = ``sum`` or ``max``.  This selects how the per-atom
  force norms are aggregated to rescale the uncertainty force in biased
  dynamics.  The default is ``sum``.
* ``tau_group`` value = LAMMPS atom group name restricting the
  :math:`\tau` bias.  The default is ``all``.  See
  :ref:`pair_grace_bias_dynamics` for the precise semantics.
* ``stress_tau`` value = stress-bias coefficient :math:`\tau_s`.  The default
  is 0, which disables the stress bias.
* ``bulk_modulus`` value = bulk modulus :math:`B` (GPa) that sets the scale of
  the stress bias.  The default is 100.

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

For ``grace``:

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

UQ-enabled TensorFlow model (UQ output activates via :doc:`fix pair <fix_pair>`):

.. code-block:: LAMMPS

   pair_style grace
   pair_coeff * * /path/to/uq_saved_model W

UQ-enabled TensorFlow model with biased dynamics:

.. code-block:: LAMMPS

   pair_style grace tau 0.1
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

Description
"""""""""""

The ``grace`` pair styles compute interactions using the Graph Atomic Cluster
Expansion (GRACE) framework :ref:`(Bochkarev24) <Bochkarev20241>`,
:ref:`(Lysogorskiy25) <Lysogorskiy20251>`.

This page documents the following related pair styles:

* ``grace``
* ``grace/fs``
* ``grace/fs/kk``

The native GRACE-1L, 2L and 3L ``.npz`` models (gracemaker export) are
evaluated by :doc:`pair_style grace/runtime <pair_grace_runtime>`.

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
     - 1L: yes; 2L: one rank only (see the warning below)
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
   * - ``grace/runtime`` and ``grace/runtime/kk``
     - Native GRACE 1L / 2L / 3L ``.npz`` (gracemaker export)
     - CPU (OpenMP) and GPU (Kokkos) production runs, see
       :doc:`pair_style grace/runtime <pair_grace_runtime>`
     - Yes
     - GPU
     - Not required

TensorFlow GRACE models
"""""""""""""""""""""""

The ``grace`` style uses ``libtensorflow`` to load and execute GRACE models
stored in TensorFlow ``saved_model`` format.

.. warning::

   With a two-layer (GRACE-2L) saved model, ``pair_style grace`` gives
   correct results only on a single MPI rank.  On several ranks the energies,
   forces and stresses are wrong without any error or warning: for a
   216-atom W-Mo cell, energies differed by about 8 to 20 meV/atom and forces by
   up to 1.1 eV/Angstrom from the single-rank result on 2 and 4 ranks.  For
   GRACE-2L models on several MPI ranks or GPUs, export the model to the
   native ``kokkos.npz`` format and use
   :doc:`pair_style grace/runtime <pair_grace_runtime>`.  Single-layer
   (GRACE-1L) saved models are not affected.

GRACE TensorFlow models are just-in-time (JIT) compiled.  The first evaluation
is therefore slower than subsequent evaluations.  To maintain performance when
the number of neighbors changes, the ``grace`` style uses padding.  A
larger ``padding`` value can reduce the frequency of recompilations, but
increases memory use and the cost of each TensorFlow call.

The ``grace`` style does not chunk atoms.  The native model export with
:doc:`pair_style grace/runtime <pair_grace_runtime>` needs less memory (it
stores nothing per bond, and chunks atoms on the GPU) and runs two- and
three-layer models on several MPI ranks.

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
:doc:`fix pair <fix_pair>`, or when ``tau != 0``.  Otherwise -- including every
plain ``grace`` run, even on a UQ-capable model -- the regular non-UQ code path
runs at no extra per-step cost.  If a UQ field (or ``tau``) is requested but the
model has no ``compute_uq`` head, LAMMPS stops with an error.

The style can expose four per-atom fields through :doc:`fix pair <fix_pair>`:

* ``gamma`` = scalar extrapolation grade :math:`\gamma` (Mahalanobis-based;
  see below).
* ``gmm_cluster`` = scalar GMM cluster index :math:`k^*` for the atom
  (integer, returned as a per-atom double because *fix pair* accepts only
  doubles).
* ``atomic_sigma`` = scalar raw per-atom uncertainty :math:`\sigma_i`, the
  Mahalanobis distance of the atom's latent feature vector to its assigned
  GMM cluster (dimensionless; see below).  This is :math:`\gamma_i` before
  normalization by the per-cluster threshold; exposing it directly is useful
  for calibration plots, threshold tuning, and diagnostics.
* ``uncertainty_force`` = three-vector uncertainty force
  :math:`\mathbf{F}^{\sigma}_i = \partial \left( \sum_j \gamma_j \right) /
  \partial \mathbf{r}_i`, the gradient of the summed extrapolation grade with
  respect to the position of atom :math:`i`.  It points towards increasing
  uncertainty.  For a model exported without calibration thresholds,
  :math:`\sum_j \sigma_j` is differentiated instead.  This is the raw
  uncertainty force and is not scaled by ``tau``.


The ``uncertainty_force`` output is accumulated consistently with Newton's
third law when Newton pair communication is enabled.

**Signature dispatch.** The pair style selects the TensorFlow signature based
on which UQ outputs are requested and on the value of ``tau``, in order of
increasing cost:

* Regular ``compute`` -- when no UQ field is requested and ``tau = 0``.
* ``compute_uq_gamma_only`` -- when only ``gamma``, ``gmm_cluster``, and/or
  ``atomic_sigma`` are requested and ``tau = 0``.  Skips the
  :math:`\sigma`-gradient backward pass.  Used automatically when the saved
  model exports this signature; otherwise falls back to ``compute_uq``.
* ``compute_uq`` -- when ``uncertainty_force`` is requested or ``tau`` is
  nonzero.  Includes the :math:`\sigma`-gradient backward pass.

Example with the ``grace`` UQ head, biased dynamics, and thermo
output containing the maximum extrapolation grade, maximum total force, and
maximum raw uncertainty-force magnitude:

.. code-block:: LAMMPS

   pair_style  grace tau 0.1
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
which includes the :math:`\tau`-scaled bias when ``tau`` is nonzero.
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
calibrated threshold :math:`\theta_{e,k^{\!*}}` of the assigned cluster:

.. math::

   \gamma_i = \frac{\sigma_i}{\theta_{e,k^{\!*}}} .

The threshold is a robust outlier fence on the :math:`\sigma` distribution of
the training atoms assigned to that cluster,

.. math::

   \theta_{e,k} = \operatorname{median}(\sigma) + 3 \cdot 1.4826 \cdot
   \operatorname{MAD}(\sigma) ,

where MAD is the median absolute deviation (:math:`1.4826\,\mathrm{MAD}`
estimates one standard deviation for a normal distribution).  Unlike a tail
percentile, this threshold is not inflated by a few outlier environments in
the training set.  :math:`\gamma_i > 1` therefore marks an atom outside the
training distribution of its cluster.  This is the default; a model exported
with ``grace_uq build --threshold-percentile Q`` uses the Q-th percentile of the
training :math:`\sigma` instead, and its :math:`\gamma` values are not
comparable with those of a default model.

.. _pair_grace_bias_dynamics:

Bias-driven dynamics (experimental)
"""""""""""""""""""""""""""""""""""

.. warning::

   Bias-driven dynamics is experimental features and may be changed in the future.

.. warning::

   Bias-driven dynamics with ``tau != 0`` or ``stress_tau != 0`` is supported
   only on a single MPI rank.  LAMMPS stops with an error during
   initialization if either is nonzero and more than one MPI rank is used.

The ``tau`` keyword sets :math:`\tau \geq 0`.  When nonzero, it adds a
relative-force uncertainty bias on top of the physical force:

.. math::

   \mathbf{f}_i = \mathbf{F}^{phys}_i + s\,\mathbf{F}^{\sigma}_i,
   \qquad
   s = \tau \cdot
       \frac{N(\|\mathbf{F}^{phys}_j\|) + \varepsilon}
            {N(\|\mathbf{F}^{\sigma}_j\|) + \varepsilon}

with :math:`\varepsilon = 10^{-8}`.  The scale factor :math:`s` is one
system-wide scalar applied uniformly to every atom.  The aggregator
:math:`N(\cdot)` over all atoms :math:`j` is selected by ``tau_norm``:

* ``sum`` = use the sum of force norms over atoms.  This is the default.
* ``max`` = use the maximum force norm over atoms.

:math:`\mathbf{F}^{phys}_i` is the force of the ``grace`` model alone, so
with :doc:`pair_style hybrid/overlay <pair_hybrid>` the forces of other
sub-styles (for example ``zbl``) do not enter :math:`s`.

By default, the :math:`\tau` bias is applied to every atom and the
aggregator :math:`N(\cdot)` runs over every atom.  When ``tau_group`` names a
LAMMPS atom group, both the reduction and the force application are restricted
to atoms in that group: :math:`s` is built from the in-group atoms only and
:math:`s\,\mathbf{F}^{\sigma}_i` is added only to in-group atoms.  Atoms
outside the group experience the physical force only.  ``gamma``,
``gmm_cluster``, and ``uncertainty_force`` are still exposed for every local
atom regardless of ``tau_group``.

The bias force is not conservative and has no pairwise decomposition.  It
adds nothing to the energy, the global virial, or the per-atom stress.

The ``stress_tau`` keyword adds a separate bias on the cell, for use with a
barostat such as :doc:`fix npt <fix_nh>`.  It adds a pressure tensor of fixed
magnitude :math:`\tau_s B` along the unit direction of the virial of the
uncertainty force:

.. math::

   \mathbf{P} = \mathbf{P}^{phys} + \tau_s B\,
   \frac{\mathbf{W}^{\sigma}}{\|\mathbf{W}^{\sigma}\|},
   \qquad
   W^{\sigma}_{ab} = \sum_{\text{bonds } ij} r_{ij,a}\,F^{\sigma}_{ij,b} ,

where :math:`\mathbf{F}^{\sigma}_{ij}` is the per-bond uncertainty force and
the norm is taken over the six Voigt components.  The diagonal components are
:math:`W^{\sigma}_{aa} = \partial(\sum_j \gamma_j)/\partial\varepsilon_{aa}`,
so the bias expands the cell along directions in which tensile strain raises
the extrapolation grade and compresses it otherwise.  :math:`B` only sets the
scale of the bias and does not have to be the exact bulk modulus of the
system; the default of 100 GPa is typical for metals.  The stress bias is
applied to the global virial only, so per-atom stress is unaffected, and it is
independent of ``tau`` and ``tau_group``.

.. _pair_grace_runtime_tuning:

Runtime tuning
""""""""""""""

The following ``grace`` UQ / biased-dynamics parameters can be changed between
``run`` commands with :doc:`pair_modify <pair_modify>`:

.. code-block:: LAMMPS

   pair_modify tau 0.5
   pair_modify tau_norm max
   pair_modify tau_group boundary
   pair_modify stress_tau 0.02 bulk_modulus 150

The current value of ``tau`` can be accessed from the pair style through the
LAMMPS ``extract`` interface using the key ``"tau"``.

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

The ``grace`` style additionally accepts five
biased-dynamics-specific :doc:`pair_modify <pair_modify>` keywords:

* ``pair_modify tau value`` = reset the relative-force biased-dynamics
  coefficient :math:`\tau`.  A nonzero value is rejected with an error when
  more than one MPI rank is used.
* ``pair_modify tau_norm sum|max`` = switch the norm used in the
  :math:`\tau` rescaling.
* ``pair_modify tau_group name`` = restrict the :math:`\tau` bias to
  atoms in the named LAMMPS group.  Pass ``all`` to remove the restriction.
* ``pair_modify stress_tau value`` = reset the stress-bias coefficient
  :math:`\tau_s`.  A nonzero value is rejected with an error when more than
  one MPI rank is used.
* ``pair_modify bulk_modulus value`` = reset the bulk modulus :math:`B` (GPa)
  used by the stress bias.

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
* The TensorFlow-based style ``grace`` is built only when LAMMPS is
  configured with ``-D ENABLE_GRACE_TENSORFLOW=yes`` and is then linked with
  ``libtensorflow`` (see :ref:`ML-PACE <ml-pace>`).
* Uncertainty quantification output and biased dynamics in ``grace`` require a
  UQ-enabled saved model that exports the ``compute_uq`` signature.
* Biased dynamics with ``tau != 0`` or ``stress_tau != 0`` in ``grace`` is
  supported only on a single MPI rank.
* ``grace`` with a two-layer (GRACE-2L) saved model gives correct results only
  on a single MPI rank; LAMMPS does not check this.
* ``grace/fs`` does not require TensorFlow.
* ``grace/*/kk`` requires ``newton on`` and supports only one CPU thread per
  MPI rank.

Further reading
"""""""""""""""

See `gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_ for more
details.

Related commands
""""""""""""""""

:doc:`pair_style grace/runtime <pair_grace_runtime>`, :doc:`pair_style pace <pair_pace>`,
:doc:`fix pair <fix_pair>`

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
       ``reduce_padding = 0.2``, ``tau = 0``, ``tau_norm = sum``,
       ``tau_group = all``, ``stress_tau = 0``, ``bulk_modulus = 100``
   * - ``grace/fs`` and ``grace/fs/kk``
     - ``extrapolation = off``, ``chunksize = 4096``

----------

.. _Bochkarev20241:

**(Bochkarev24)** Bochkarev, Lysogorskiy, Drautz, Phys. Rev. X, 14, 021036
(2024).

.. _Lysogorskiy20251:

**(Lysogorskiy25)** Lysogorskiy, Bochkarev, Drautz, arXiv:2508.17936 (2025).
