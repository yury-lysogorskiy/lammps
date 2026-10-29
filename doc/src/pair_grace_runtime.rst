.. index:: pair_style grace/runtime
.. index:: pair_style grace/runtime/kk

pair_style grace/runtime command
================================

Accelerator Variants: *grace/runtime/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grace/runtime keyword value ...
   pair_style grace/runtime/kk keyword value ...

* zero or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *spline* or *spline_layers* or *spline_mb* or *store_fwd* or *fc_fold* or
               *store_l1* or *block* or *fuse_radial* or *noshape* or *chunksize* or
               *debug_no_energy_only_calc*
       *spline* value = N
         N = intervals of the cubic radial tables (at least 16); GRACE-3L: of layer 3,
             layers 1 and 2 use N/2
       *spline_layers* values = N1 N2 N3 (GRACE-3L)
         N1, N2, N3 = intervals of the layer 1, 2 and 3 tables (each at least 16)
       *spline_mb* value = *auto* or MB (GRACE-3L)
         MB = memory budget of the radial tables, 0 = always evaluate the radial MLP
       *store_fwd* value = *auto* or *yes* or *no* (GRACE-3L)
         keep the layer-1 and layer-2 forward state for the backward pass
       *fc_fold* value = *yes* or *no* (GRACE-3L)
         fold the linear layers that feed a tensor product into its projections
       *store_l1* value = *yes* or *no* (GRACE-2L)
         keep the layer-1 forward state for the backward pass
       *block* value = N (GRACE-3L, CPU)
         N = atoms of one element a thread processes together (1 to 64)
       *fuse_radial* value = *yes* or *no* (GRACE-3L, GPU)
         evaluate the radial tables inside the kernels that use them
       *noshape* (GRACE-1L and 2L, GPU)
         use the generic CG-product kernels (diagnostics)
       *chunksize* value = *auto* or N (GPU)
         auto = atoms per pass from the free GPU memory
         N = atoms per pass, 0 = all atoms in one pass
       *debug_no_energy_only_calc* (GPU)
         run the force kernels also in energy-only requests (testing)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style grace/runtime
   pair_coeff * * GRACE-2L-OMAT-medium/kokkos.npz Mo Nb Ta W

   pair_style grace/runtime
   pair_coeff * * GRACE-3L-OMAT-large/kokkos.npz Sr Ti O

   # GPU: lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.lmp
   pair_style grace/runtime chunksize 20000
   pair_coeff * * GRACE-1L-OMAT-large/kokkos.npz Al Li

Description
"""""""""""

Pair style *grace/runtime* evaluates the Graph Atomic Cluster Expansion
(GRACE) models :ref:`(Bochkarev24) <Bochkarev20242>`,
:ref:`(Lysogorskiy25) <Lysogorskiy20252>` with the kernels of the
grace-runtime library (libgrace).  The model is the ``.npz`` written by the
gracemaker export utilities (``kokkos.npz``, see
`gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_); TensorFlow
saved models are read by :doc:`pair_style grace <pair_grace>`.  The
``pair_coeff`` command reads the model file and picks the evaluator from it:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Model
     - *grace/runtime* (CPU)
     - *grace/runtime/kk* (GPU)
   * - GRACE-1L
     - fused single-layer evaluator, no communication between the layers
     - fused single-layer evaluator
   * - GRACE-2L
     - fused two-layer evaluator; the layer-2 messages of ghost atoms come
       from their owners by forward communication, their adjoints go back by
       reverse communication
     - the same, on the device
   * - GRACE-3L (GRACE-3L-OMAT-large architecture)
     - three-layer evaluator, two message layers
     - the same, on the device

ACE models (pacemaker) are evaluated by :doc:`pair_style pace/runtime
<pair_pace_runtime>`; a PACE model in *grace/runtime* is an error that names
that style.

All evaluators use one precision scheme: the bond vectors, their lengths and
unit vectors are computed in 64-bit (GRACE-3L on the GPU: also the cutoff
envelope and the Chebyshev radial basis), the network and the radial tables
(built in 64-bit) are evaluated in 32-bit, and energies, forces and virials are
accumulated in 64-bit.  Results agree with 64-bit evaluations of the same
models to 32-bit precision (about :math:`10^{-6}` eV/atom and
:math:`5 \cdot 10^{-6}` eV/Angstrom for the published foundation models).

The radial MLPs are replaced by per-interval cubic tables in
:math:`x = r / r_c` built at initialization (*spline*); the bond geometry,
spherical harmonics and radial functions are fused into per-atom kernels, so
nothing is stored per bond.  GRACE-3L models with many elements, whose tables
would exceed *spline_mb*, evaluate the radial MLP directly.

Keywords that do not apply to the model in ``pair_coeff`` stop with an error
naming the keyword and the model (for example *store_l1* with a GRACE-3L
model).  Keywords of the other backend that only tune speed or memory
(*block* on the GPU; *chunksize*, *noshape* and *fuse_radial* on the CPU) are
ignored with a warning, so an input runs unchanged with and without
``-sf kk``; *debug_no_energy_only_calc* on the CPU is an error.

CPU
"""

*grace/runtime* is a plain CPU style of the ML-PACE package; KOKKOS is not
needed.  Within each MPI rank, OpenMP threads (``OMP_NUM_THREADS``, one per
physical core recommended) process atoms of one element in small batches
(GRACE-3L: blocks of up to *block* atoms that run a whole layer, with the dense
contractions as matrix products), with SIMD over channels.  The style
vectorizes with the widest instruction set enabled at compile time, so LAMMPS
should be compiled for the CPU it runs on, e.g. with ``-march=native`` in
``CMAKE_CXX_FLAGS``; on x86-64 the style prints a warning when the CPU
supports AVX2 or AVX-512 but the binary was compiled without it.  GCC is the
recommended compiler.

For GRACE-3L, when LAMMPS is configured with CMake, ``-D GRACE_BLAS=AUTO``
(the default) links a BLAS library that exports ``cblas_sgemm`` (LP64
interface) if one is found, ``ON`` requires one and ``OFF`` uses the built-in
kernel.  Sequential MKL (``-D BLA_VENDOR=Intel10_64lp_seq``) or an OpenMP
build of OpenBLAS is recommended; an OpenBLAS built with pthreads cannot be
called from several OpenMP threads and is replaced by the built-in kernel at
run time, with a warning.  The log reports the GEMM path in use.  The default
*spline_mb auto* is a quarter of the host memory divided by the MPI ranks on
the host, at least 2048 MB (about 14 elements) and at most the available
memory per rank.

GRACE-2L with *store_l1 yes* and GRACE-3L with *store_fwd yes* keep the
forward state of all atoms (2L: layer 1, 3L: layers 1-2) while it fits half
the available host memory divided by the MPI ranks on the host (8 GB per rank
if the host memory cannot be queried), else recompute it in the backward pass,
with a warning.

GPU
"""

*grace/runtime/kk* (alias *grace/runtime/kk/device*) runs on GPU Kokkos
backends (CUDA, HIP); with a host-only Kokkos build it stops with an error.
It requires Newton's third law and a half neighbor list, e.g.
``-pk kokkos newton on neigh half``.  *chunksize* sets the number of local
atoms per engine pass, ``0`` for all in one pass, or ``auto`` (the default):
the style measures the work buffer memory per atom and takes all local atoms
in one pass while that fits 70% of the free GPU memory (shared evenly by the
MPI ranks on the same GPU), else the largest power of two that fits, at least
1024; the choice is printed to the log and revised only when the number of
atoms or neighbors grows by more than 10%.  Results depend on the chunk size
only through floating-point rounding.  For GRACE-1L and 2L the style uses
generated kernels for the CG products of the GRACE-2L-(SMAX-)OMAT medium and
large and GRACE-1L large model shapes (other shapes, and *noshape*, use generic
kernels).  For GRACE-3L the default *spline_mb auto* is a quarter of the GPU
memory per process, at least 2048 MB, and *store_fwd auto* keeps the forward
state when it fits.

With several GPUs, run one MPI rank per GPU, each seeing only its own GPU
(``-k on g 1``).  GRACE-2L and 3L exchange per-atom messages after each layer
between the ranks, through GPU buffers when LAMMPS uses a GPU-aware MPI
(``-pk kokkos gpu/aware on``, the default; ``off`` stages them through host
buffers).  For GPU-aware MPI on HPE Cray systems
(``cray-mpich``, ``MPICH_GPU_SUPPORT_ENABLED=1``), give the job step all
GPUs of a node and select one per rank with a wrapper script:
``srun --gpus-per-task`` or ``--gpu-bind`` puts each rank in its own GPU
cgroup, where ``cray-mpich`` cannot map the buffers of the other GPUs, and
every KOKKOS style aborts.  For example, on nodes with 8 AMD GPUs:

.. code-block:: bash

   #!/bin/bash
   # select_gpu: the rank with node-local index i sees only GPU i
   # (on NVIDIA GPUs: CUDA_VISIBLE_DEVICES)
   export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID
   exec "$@"

Make it executable (``chmod +x select_gpu``) and start LAMMPS through it:

.. code-block:: bash

   # job script
   #SBATCH --nodes=1
   #SBATCH --ntasks-per-node=8
   #SBATCH --gpus-per-node=8
   export MPICH_GPU_SUPPORT_ENABLED=1
   srun /path/to/select_gpu lmp -in in.lmp -k on g 1 -sf kk -pk kokkos newton on neigh half

The LAMMPS binary must be linked with the GPU transfer library of
``cray-mpich``; the Cray compiler wrappers do this when the GPU target module
is loaded at build time.

.. _pair_grace_uq_kokkos:

Uncertainty quantification
""""""""""""""""""""""""""

For a model file with UQ artifacts, *grace/runtime* and *grace/runtime/kk*
compute the per-atom uncertainty quantities of the TensorFlow
:ref:`grace UQ head <pair_grace_uq>` on the steps where
:doc:`fix pair <fix_pair>` requests them, in double precision; the CPU style
computes the random projection of the basis in single precision (its cost
dominates there) and the norms and the GMM distance in double precision,
which changes :math:`\gamma_i` by about :math:`10^{-5}` relative.

The style exposes three per-atom fields through :doc:`fix pair <fix_pair>`,
with the same names and semantics as the TensorFlow ``grace`` style:

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

For two- and three-layer models the basis has two or three invariant blocks;
the density channels are ordered ``[full, block0, block1, ...]`` following the
block order baked into the artifacts (the projection-matrix row layout follows
the same order).

Only UQ artifacts written with schema version 6 are accepted (the
L2-normalized projection + density-channel feature described above; the
``.npz`` carries ``uq_schema_version = 6`` with ``uq_rp_normalize = 1``,
``uq_rp_add_density_channel = 1``, ``uq_feature_transform = 0``); with older
artifacts the style runs without UQ and logs why.  The artifacts must also
provide GMM clusters for **every** model element.  Unlike the TensorFlow
``grace`` UQ path, there is **no** ``uncertainty_force``, ``tau``/HAL biasing,
or :math:`\partial\gamma/\partial\mathbf{r}` gradient.

.. code-block:: LAMMPS

   pair_style grace/runtime
   pair_coeff * * grace_kokkos_uq.npz Mo Nb Ta W

   fix gp  all pair 1 grace/runtime gamma          1
   fix as  all pair 1 grace/runtime atomic_sigma   1
   fix gm  all pair 1 grace/runtime gmm_cluster    1

   compute max_gamma all reduce max f_gp
   thermo_style custom step pe c_max_gamma
   dump uq all custom 100 uq.dump id type x y z f_gp f_as f_gm

With ``-sf kk`` the unsuffixed style name in ``fix pair`` refers to
*grace/runtime/kk*.

.. _pair_grace_decouple:

4D decoupling (alchemical free energies)
""""""""""""""""""""""""""""""""""""""""

*grace/runtime* and *grace/runtime/kk* can decouple two groups of atoms from
each other ("4D decoupling", P. Peña Cano), e.g. for solvation free energies
by thermodynamic integration (TI) or the Bennett acceptance ratio (BAR):

.. code-block:: LAMMPS

   pair_modify decouple/4d group1 group2 [eq]
   pair_modify decouple/4d off

Every bond between an atom of *group1* and an atom of *group2* is stretched
along its own direction, :math:`r \to r'(r, \lambda)`; bonds within a group and
bonds of atoms in neither group are unchanged.  The optional integer *eq*
selects the shift (:math:`r_c`: cutoff of the bond):

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - eq
     - :math:`r'`
   * - 1
     - :math:`r + \lambda r_c`
   * - 2
     - :math:`r + \lambda \max(0, r_c - r)`
   * - 3
     - :math:`\sqrt{r^2 + (\lambda r_c)^2}`
   * - 4 (default)
     - :math:`\sqrt{r^2 + \lambda^2 \max(0, r_c^2 - r^2)}`

At :math:`\lambda = 0` the bonds are unchanged.  A bond with :math:`r' \ge r_c`
no longer contributes (at :math:`\lambda = 1` this holds for every cross-group
bond with eq 1, 2 and 4).  The forces are the exact gradient of the stretched
energy with respect to the true positions, and the virial uses the true bond
vectors.  The UQ outputs (``gamma``, ``atomic_sigma``, ``gmm_cluster``) are
those of the stretched environments, i.e. of what the model is evaluated on.

:math:`\lambda` is a global parameter of the style named ``lambda`` (default 0),
set with :doc:`fix adapt <fix_adapt>`; :math:`\partial E / \partial \lambda` is
the first extra quantity of the style, available as ``c_ID[1]`` of
:doc:`compute pair <compute_pair>` (summed over the MPI ranks).  With a
suffix (``-sf kk``) the unsuffixed style name in ``fix adapt`` and
``compute pair`` refers to the accelerated style.  Thermodynamic integration
at one :math:`\lambda`:

.. code-block:: LAMMPS

   group solute molecule 1
   group solvent subtract all solute
   pair_style grace/runtime
   pair_coeff * * GRACE-2L-model/kokkos.npz H C N O
   pair_modify decouple/4d solute solvent
   variable lam equal 0.3
   fix lam all adapt 0 pair grace/runtime lambda * * v_lam
   compute dedl all pair grace/runtime
   thermo_style custom step temp pe c_dedl[1]

For BAR, the trajectory of one window can be re-evaluated at a neighbouring
:math:`\lambda` with :doc:`rerun <rerun>`; ``fix adapt`` sets
:math:`\lambda` before every snapshot:

.. code-block:: LAMMPS

   variable lam equal 0.4
   fix lam all adapt 0 pair grace/runtime lambda * * v_lam
   rerun window_0.3.dump dump x y z

No atom may belong to both groups (checked at every run).  In the energy-only
mode (Monte Carlo trial moves) :math:`\partial E / \partial \lambda` is not
computed and reads 0.

Energy-only calculation
"""""""""""""""""""""""

When the caller of the pair style asks for the energy only (``ENERGY_ONLY``,
e.g. the Monte Carlo fixes :doc:`fix atom/swap <fix_atom_swap>`,
:doc:`fix sgcmc <fix_sgcmc>` or :doc:`fix gcmc <fix_gcmc>`), the style skips
the backward pass: no forces and no virials.

Pair coefficients
"""""""""""""""""

Only a single pair_coeff command is used, with ``* *``, the model file and the
element of each LAMMPS atom type; with :doc:`pair_style hybrid <pair_hybrid>`,
``NULL`` leaves a type to the other sub-styles.  A *pair_style* command with the same style
after *pair_coeff* rebuilds the evaluator from the model file of that
*pair_coeff* with the new keywords (checked against the model as at
*pair_coeff*); *pair_coeff* need not be repeated.  Under
:doc:`pair_style hybrid <pair_hybrid>` a new *pair_style* creates the
sub-styles anew, so their *pair_coeff* commands must be given again.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Mixing is not used by these pair styles.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options; they accept the keyword ``decouple/4d``
(:ref:`4D decoupling <pair_grace_decouple>`).

These pair styles do not write their information to :doc:`binary restart
files <restart>`, since it is stored in model files.  The ``pair_style`` and
``pair_coeff`` commands must therefore be specified again in an input script
that reads a restart file.

The styles provide energies, forces, global and per-atom energy and virial
(``compute centroid/stress/atom`` uses the per-atom virial), and work with any
number of MPI ranks.

----------

Restrictions
""""""""""""

These pair styles are part of the ML-PACE package, built with the
grace-runtime styles (CMake option ``GRACE_RUNTIME``, on by default).  They
require metal units, atom IDs and newton pair on.  *grace/runtime/kk* needs a
GPU Kokkos backend.  GRACE-3L models must be of the GRACE-3L-OMAT-large
architecture (the foundation model and its fine-tunes); others are rejected
when ``pair_coeff`` reads them.

Related commands
""""""""""""""""

:doc:`pair_style grace <pair_grace>`,
:doc:`pair_style pace/runtime <pair_pace_runtime>`,
:doc:`fix pair <fix_pair>`, :doc:`fix adapt <fix_adapt>`,
:doc:`compute pair <compute_pair>`, :doc:`rerun <rerun>`

Default
"""""""

spline = 1024, spline_layers = from spline, spline_mb = auto, store_fwd = auto
(CPU: yes), fc_fold = yes, store_l1 = yes, block = 16, fuse_radial = yes,
chunksize = auto, decouple/4d off

----------

.. _Bochkarev20242:

**(Bochkarev24)** Bochkarev, Lysogorskiy, Drautz, Phys. Rev. X, 14, 021036
(2024).

.. _Lysogorskiy20252:

**(Lysogorskiy25)** Lysogorskiy, Bochkarev, Drautz, arXiv:2508.17936 (2025).
