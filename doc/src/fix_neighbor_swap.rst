.. index:: fix neighbor/swap

fix neighbor/swap command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID neighbor/swap N X seed T R0 voro keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* neighbor/swap = style name of this fix command
* N = invoke this fix every N steps
* X = number of swaps to attempt every N steps
* seed = random # seed (positive integer)
* T = scaling temperature of the MC swaps (temperature units)
* R0 = scaling swap probability of the MC swaps (distance units)
* voro = valid voronoi compute id (compute voronoi/atom)
* two or more keyword/value pairs may be appended to args
* keywords *types* and *diff* are mutually exclusive, but one must be specified
* keyword = *types* or *diff* or *ke* or *region* or *rates*

  .. parsed-literal::

       *types* values = two or more atom types (Integers in range [1,Ntypes] or type labels)
<<<<<<< HEAD
       *diff* values = one atom type
       *ke* value = *yes* or *no*
=======
       *ke* value = *no* or *yes*
         *no* = no conservation of kinetic energy after atom swaps
>>>>>>> a3bc1a6c0b09ee7a84167ab8aa626e4a93443846
         *yes* = kinetic energy is conserved after atom swaps
         *no* = no conservation of kinetic energy after atom swaps
       *region* value = region-ID
         region-ID = ID of region to use as an exchange/move volume
       *rates* values = V1 V2 . . . Vntypes values to conduct variable diffusion for different atom types (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   compute voroN all voronoi/atom neighbors yes
   fix mc all neighbor/swap 10 160 15238 1000.0 3.0 diff 2 voro voroN
   fix myFix all neighbor/swap 100 1 12345 298.0 3.0 region my_swap_region types 5 6 voro voroN
   fix kmc all neighbor/swap 1 100 345 1.0 3.0 diff 3 rates 3 1 6 voro voroN

Description
"""""""""""

.. versionadded:: TBD

<<<<<<< HEAD
This fix performs Monte-Carlo (MC) evaluations to enable kinetic 
Monte Carlo (kMC)-type behavior during MD simulation by allowing 
neighboring atoms to swap their positions. In constrast to the :doc:`fix
atom/swap <fix_atom_swap>` command which swaps pairs of atoms anywhere
in the simulation domain, the restriction of the MC swapping to 
neighbors enables a hybrid MD/kMC-like simulation.

Neighboring atoms are defined by using a Voronoi tesselation performed
by the :doc:`compute voronoi/atom <compute_voronoi_atom>` command.
Two atoms are neighbors if their Voronoi cells share a common face
(3d) or edge (2d).

The selection of a swap neighbor is made using a distance-based
criterion for weighting the selection probability of each swap, in the
same manner as kMC selects a next event using relative probabilities.
The acceptance or rejection of each swap is determined via the
Metropolis criterion after evaluating the change in system energy due
to the swap.

A detailed explanation of the original implementation of this
algorithm can be found in :ref:`(Tavenner 2023) <_TavennerMDkMC>`
where it was used to simulated accelerated diffusion in an MD context.
=======
This fix computes Monte Carlo (MC) evaluations to enable kinetic 
Monte Carlo (kMC) behavior during an MD simulation by allowing 
neighboring atoms to swap their positions. This creates a hybrid MD/kMC 
simulation
where atoms are swapped only with their neighbors, and the swapping acceptance is
perfomed by evaluating the change in system energy using the Metropolis criterion.
Neighboring atoms are selected using a Voronoi tesselation approach. A detailed
explination of the original implementation of this procedure can be found in
:ref:`(Tavenner 2023) <_TavennerMDkMC>` as originally intended for simulating
accelerated diffusion in an MD context.

Simulating inherently kinetically-limited behaviors which rely on rare events
(such as atomic diffusion) is challenging for traditional Molecular Dynamics
approaches, since thermal vibration motion occurs on a timescale much 
shorter than local structural changes such as vacancy hopping.
To address
this challenge, an approach from kMC simulations is adopted where rare events can
be sampled at selected rates. By selecting such swap behaviors, the process
of atomic diffusion can be approximated during an MD simulation, effectively
decoupling the MD atomic vibrational time and the timescale of atomic hopping.

To achieve such simulations, this algorithm takes the following approach. First,
the MD simulation is stopped *N* steps to attempt *X* atom swaps.
Voronoi neighbors are computed for all valid swap atoms. 
For each swap attempt, one atom
I is selected at random from the global list of valid atoms. One if its Voronoi neighbors
that is a valid atom to swap is then selected. The neighbor atom ID is communicated to all 
processors, as it may be owned by a different processor. 
The two atom types are then swapped, and the change in system energy from before 
the swap is used to evaluate the Metropolis acceptance criterion. This evaluation of the energy 
change is a global calculation, such that it has a computational cost similar to
that of an MD timestep. If the swap is accepted from the Metropolis criterion, the
atoms remain swapped. If the swap is rejected, the atoms are reverted to their original types. This
process of MC evaluation is repeated for a given number of iterations until the
original MD simulation is resumed from the new state, where any successfully
swapped atoms have changed type, though the global system balance is preserved.

A few key notes regarding this implementation are as follows. The parallel
efficiency of this algorithm is similar to that of other MC approaches, i.e,
the global potential energy must be calculated after each attempted swap.
Efficiency is sensitive to the maximum cutoff distance for the pair style, 
since the neighbor list will need to be rebuilt between swap events. 
Limitations are imposed on the Voronoi neighbors
to restrict swapping of atoms that are outside of a reasonable cutoff.
>>>>>>> a3bc1a6c0b09ee7a84167ab8aa626e4a93443846

Simulating inherently kinetically-limited behaviors which rely on rare
events (such as atomic diffusion in a solid) is challenging for
traditional MD since its relatively short timescale will not naturally
sample many events. This fix addresses this challenge by allowing rare
neighbor hopping events to be sampled in a kMC-like fashion at a much
faster rate (set by the specified *N* and *X* parameters).  This enables
the processes of atomic diffusion to be approximated during an MD
simulation, effectively decoupling the MD atomic vibrational timescale
and the atomic hopping (kMC event) timescale.

The algorithm implemented by this fix is as follows. The MD
simulation is paused every *N* stepsA Voronoi tesselation is
performed for the current atom configuration.  Then *X* atom swaps are
attempted, one after the other.  For each swap, an atom *I* is
selected randomly from the list of atom types specified by either the
*types* or *diff* keywords.  One of *I*'s Voronoi neighbors *J* is
selected using the distance-weighted probability for each neighbor
detailed below.  The *I,J* atom IDs are communicated to all processors
so that a global energy evaluation can be performed for the post-swap
state of the system.  The swap is accepted or rejected based on the
Metropolis criterion using the energy change of the system and the
specified temperature *T*.

Here are a few comments on the computational cost of the swapping
algorithm.

(1) The cost of a global energy evaluation is similar to that of an MD
timestep.

(2) Simliar to other MC algorithms in LAMMPS, an optimized parallel efficiency
is achieved with a smaller number of atoms per processor than would typically
be used in an standard MD simulation. This is because the per-energy evaluation
cost increases relative to the balance of MD/MC steps as indicated by (1), but
the communication cost remains relatively constant for a given number of MD steps.

(3) The MC portion of the simulation will run dramatically slower if
the pair style uses different cutoffs for different atom types (or
type pairs).  This is because each atom swap then requires a rebuild
of the neighbor list to ensure the post-swap global energy can be
computed correctly.

Limitations are imposed on selection of *I,J* atom pairs to avoid
swapping of atoms which are outside of a reasonable cutoff (e.g. due
to a Voronoi tesselation near free surfaces) though the
use of a distance-weighted probabiltiy scaling. 

----------

This section gives more details on other arguments and keywords.

The random number generator (RNG) used by all the processors for MC
operations is initialized with the specified *seed*.

The distance-based probability is weighted by the specified *R0* which
sets the radius :math:`r_0` in this formula

.. math::

    p_{ij} = e^{(\frac{r_{ij}}{r_0})^2}

where :math:`p_{ij}` is the probability of selecting atom :math:`j` to
swap with atom :math:`i`.  Typically, a value for *R0* around the
average nearest-neighbor spacing is appropriate.  Since this is simply
a proability weighting, the swapping behavior is not very sensitive to
the exact value of *R0*.

<<<<<<< HEAD
The keyword *types* takes two or more atom types as its values.  Only
atoms *I* of the first atom type will be selected.  Only atoms *J* of the
remaining atom types will be considered as potential swap partners.
=======
Typically, a value around the average nearest-neighbor spacing is appropriate
for *R0*. Since this is simply a probability weighting, behavior is not
particularly sensitive to the exact value of *R0*.
>>>>>>> a3bc1a6c0b09ee7a84167ab8aa626e4a93443846

The keyword *diff* take a single atom type as its value.  Only atoms
*I* of the that atom type will be selected.  Atoms *J* of all
remaining atom types will be considered as potential swap partners.
This includes the atom type specified with the *diff* keyword to
account for self-diffusive hops between two atoms of the same type.

<<<<<<< HEAD
The keyword *voro* is required.  Its *vID* value is the compute-ID of
a :doc:`compute voronoi/atom <compute_voronoi_atom>` command like
this:
=======
The keyword *diff* is used for implementation of simulated diffusion of
a given atom type as given by *diff type*. This command selects all atom
types as acceptable swap types to a centrally selected atom of type
*type*. This includes the atom type specified by the *diff* keyword to
account for self-diffusion hops of an atom type with itself.

Keyword *voro* is currently required, and is implemented as

.. code-block:: LAMMPS

   voro compute-ID

where *compute-ID* is the ID of a corresponding Voronoi computation with
neighbor list, i.e.
>>>>>>> a3bc1a6c0b09ee7a84167ab8aa626e4a93443846

.. code-block:: LAMMPS

    compute compute-ID group-ID voronoi/atom neighbors yes

Note that the *neighbors yes* option must be enabled for use with this
fix. The group-ID should include all the atoms which this fix will
potentialy select. I.e. the group-ID used in the voronoi compute should
include the same atoms as that indicated by the *types* keyword. If the
*diff* keyword is used, the group-ID should include atoms of all types
in the simulation.

The keyword *ke* takes *yes* (default) or *no* as its value.  It two
atoms are swapped with different masses, then a value of *yes* will
rescale their respective velocities to conserve the kinetic energy of
the system.  A value of *no* will perform no rescaling, so that
kinetic energy is not conserved.  See the restriction on this keyword
below.

The *region* keyword takes a *region-ID* as its value.  If specified,
then only atoms *I* and *J* within the geometric region will be
considered as swap partners.  See the :doc:`region <region>` command
for details.  This means the group-ID for the :doc:`compute
voronoi/atom <compute_voronoi_atom>` command also need only contain
atoms within the region.

<<<<<<< HEAD
The keyword *rates* can modify the swap rate based on the type of atom
*J*.  Ntype values must be specified, where Ntype = the number of atom
types in the system.  Each value is used to scale the probability
weighting given by the equation above.  In the third example command
above, a simulation has 3 atoms types.  Atom *I*s of type 1 are
eligible for swapping.  Swaps may occur with atom *J*s of all 3 types.
Assuming all *J* atoms are equidistant from an atom *I*, *J* atoms of
type 1 will be 3x more likely to be selected as a swap partner than
atoms of type 2.  And *J* atoms of type 3 will be 6.5x more likely to
be selected than atoms of type 2.  If the *rates* keyword is not used,
=======
for using *fix neighbor/swap* with *diff 2*.

If atoms in the specified group are not in the voro calculated group
they will not be considered for swapping.

The keyword *rates* can modify the swap rate for each swapped type by
values where the adjusted rates values are given in order of increasing
atom type.  The number of rates provided must equal the number of atom
types in the simulation.  In the third provided example above, a
simulation is conducted with three atom types where the third atom type
is the one sampled for attempted swaps. All three atom types are
considered valid swaps, but atoms of type 1 will be selected three times
as often as atoms of type 2. Conversely, atoms of type 3 are six times
more likely to be selected than atoms of type two and twice as likely as
atoms of type 1.

Finally, the *region* keyword is implemented as in other atomic fixes,
where the *region region-ID* command indicates that atom swaps only be
considered in the area given by *region-ID*. If only atoms of certain
groups are expected to be in this region, the corresponding compute
voronoi command can be adjusted accordingly.

Either the *types* or *diff* keyword must be specified to select atom
types for swapping
>>>>>>> a3bc1a6c0b09ee7a84167ab8aa626e4a93443846


Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the random number generator
seed, the next timestep for MC exchanges, the number of exchange
attempts and successes, etc.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a global vector of length 2, which can be accessed
by various :doc:`output commands <Howto_output>`.  The vector values are
the following global cumulative quantities:

  #. swap attempts
  #. swap accepts

The vector values calculated by this fix are "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.  Also this fix requires that the :ref:`VORONOI
package <PKG-VORONOI>` be installed, otherwise the fix will not be
compiled.

The :doc:`compute voronoi/atom <compute_voronoi_atom>` command
referenced by keyword *voro* must return neighboring atoms as
illustrated in the examples above.

If this fix is used with systems that do not have per-type masses
(e.g. atom style sphere), the *ke* keyword must be set to *off* since
the implemented algorithm will not be able to re-scale velocities
properly.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`compute voronoi/atom <compute_voronoi_atom>`
:doc:`delete_atoms <delete_atoms>`, :doc:`fix gcmc <fix_gcmc>`,
:doc:`fix atom/swap <fix_atom_swap>`, :doc:`fix mol/swap <fix_mol_swap>`,
:doc:`fix sgcmc <fix_sgcmc>`

Default
"""""""

The option defaults are *ke* = yes and *rates* = 1 for all atom types.

----------

.. _TavennerMDkMC:

**(Tavenner 2023)** J Tavenner, M Mendelev, J Lawson, Computational
 Materials Science, 218, 111929 (2023).
