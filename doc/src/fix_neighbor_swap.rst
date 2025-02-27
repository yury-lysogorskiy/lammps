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
* one or more keyword/value pairs may be appended to args
* keywords *types* and *diff* are mutually exclusive, but one must be specified
* keyword = *types* or *diff* or *rates* or *ke* or *region*

  .. parsed-literal::

       *types* values = two or more atom types (1-Ntypes or type label)
       *ke* value = *no* or *yes*
         *no* = no conservation of kinetic energy after atom swaps
         *yes* = kinetic energy is conserved after atom swaps
       *region* value = region-ID
         region-ID = ID of region to use as an exchange/move volume
       *diff* values = one atom type
       *rates* values = V1 V2 . . . Vntypes values to conduct variable diffusion for different atom types (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   fix mc all neighbor/swap 10 160 15238 1000.0 3.0 diff 2 voro voroN
   fix myFix all neighbor/swap 100 1 12345 298.0 3.0 region my_swap_region types 5 6 voro voroN
   fix kmc all neighbor/swap 1 100 345 1.0 3.0 diff 3 rates 3 1 6 voro voroN

Description
"""""""""""

.. versionadded:: TBD

This fix computes Monte-Carlo (MC) evaluations to enable kinetic 
Monte Carlo (kMC)-type behavior during MD simulation through only allowing 
neighboring atom swaps. This creates a hybrid type simulation of MDkMC simulation
where atoms are only swapped with their neighbors, but the swapping acceptance is
perfomed by evaluating the change in system energy using the Metropolis Criterion.
Neighboring atoms are selected using a Voronoi tesselation approach. A detailed
explination of the original implementation of this procedure can be found in
:ref:`(Tavenner 2023) <_TavennerMDkMC>` as originally intended for simulating
accelerated diffusion in an MD context.

Simulating inherently kineticly driven behaviors which rely on rare events
(such as atomic diffusion) is challenging for traditional Molecular Dynamics
approaches since simulations are restricted in their time-scale of events.
Since thermal vibration motion occurs on a timescale much shorter than the movement
of vacancies, such behaviors are challenging to model simultaneously. To address
this challenge, an approach from kMC simulations is adpoted where rare events can
be sampled at selected rates. By selecting such swap behaviors, the process
of atomic diffusion can be approximated during an MD simulation, effectively
decoupling the MD atomic vibrational time and the timescale of atomic hopping.

To achieve such simulations, this algorithm takes the following approach. First,
the MD simulation is stopped after a given number of steps to perform atom swaps.
Given this instantaneous configuration from the MD simulation, Voronoi neighbors
are computed for all valid swap atoms. From the list of valid swap atoms, one atom
I is selected at random across the entire simulation. One if its Voronoi neighbors
that is a valid atom to swap is then selected. The atom ID is communicated to all 
processors, such that if the neighbors are on different processors the swap still 
occurs. The two atom types are swapped, and the change in system energy from before 
the swap is compared using the Metropolis Criterion. This evaluation of the energy 
change is a global calculation, such that it has a computational cost similar to
that of an MD timestep. If the swap is accepted from the Metropolis Criterion, the
atoms remain swapped. Else, the atoms are returned to their original types. This
process of MC evaluation is repeated for a given number of iterations until the
original MD simulation is resumed from the new state, where any successfully
swapped atoms have changed type, though the global system balance is preserved.

A few key notes regarding this implementation are as follows. The parallel
efficiency of this algorithm is similar to that of other MC approaches. I.e,
due to the additional energy calculations for the MC steps, efficiency is
improved with a smaller number of atoms per processor than standalone MD simulation
since there is more weighting on the calculation of a given atomic domain and
minor additonal communication load. Communication of the atom ids to be swapped
between processors is negligible. Efficiency will additionally be much worse for
pair styles with different per-atom cutoffs, since the neighbor list will need to
be rebuilt between swap events. Limitations are imposed on the Voronoi neighbors
to restrict swapping of atoms which are outside of a reasonable cutoff.

Input Parameters Usage
"""""""""""

The fix is called every *N* timesteps and attempts *X* swaps. The system
is initialized with a random seed, using a temperature *T* for
evaluating the MC energy swaps. The distance-based probability is
weighted according to *R0* which sets the radius :math:`r_0` for the
weighting

.. math::

    p_{ij} = e^{(\frac{r_{ij}}{r_0})^2}

where :math:`p_{ij}` is the probability of selecting atoms :math:`i` and
:math:`j` for an evaluated swap.

Typically, a value around the average nearest-neighbor spacing is appropriate
for *R0*. Since this is simply a proability weighting, behavior is not
particularly sensitive to the exact value of *R0*.

The keyword *types* is submitted with two or more atom types as the
value.  Atoms of the first atom type are swapped with valid neighbors of
all the remaining atom types.

The keyword *diff* is used for implementation of simulated diffusion of
a given atom type as given by *diff type*. This command selects all atom
types as acceptable swap types to a centrally selected atom of type
*type*. This includes the atom type specified by the diff keyword to
account for self-diffusion hops of an atom type with itself.

Keyword *voro* is currently required, and is implemented as

.. code-block:: LAMMPS

   voro compute-ID

where *compute-ID* is the ID of a corresponding Voronoi computation with
neighbor list, i.e.

.. code-block:: LAMMPS

    compute compute-ID group-ID voronoi/atom neighbors yes

The group selected for computing *voro* should correspond to all the
potential atoms to be swapped at the initial step, i.e.

.. code-block:: LAMMPS

   group group-ID type 2

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

Keyword Summary
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

types = Select random atom matching first type as type I, remaining
atom types are valid for selecting atom J.
diff = Select random atom of this type as atom I, all atoms are valid
for type J.
ke = re-scale velocities when atoms are swapped based on difference in
mass
region = select only atoms I and J from region
rates = pre-factor modification to the J atom selection probability
based on atom type.


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
doc page for more info.  Also this fix requires that the
:ref:`VORONOI package <PKG-VORONOI>` is installed, otherwise the fix
will not be compiled.

The voronoi command specified by *voro* must return neighboring atoms.

When this fix is used with a :doc:`hybrid pair style <pair_hybrid>`
system, only swaps between atom types of the same sub-style (or
combination of sub-styles) are permitted.

If this fix is used with systems that do not have per-type masses
(e.g. atom style sphere), the ke flag must be set to off since the
implemented algorithm will not be able to re-scale velocity properly.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`compute voronoi/atom <compute_voronoi_atom>`
:doc:`delete_atoms <delete_atoms>`, :doc:`fix gcmc <fix_gcmc>`,
:doc:`fix atom/swap <fix_atom_swap>`, :doc:`fix mol/swap <fix_mol_swap>`,
:doc:`fix sgcmc <fix_sgcmc>`

Default
"""""""

The option defaults are *ke* = yes, *rates* = 1 for all
atom types.

----------

.. _TavennerMDkMC:

**(Tavenner 2023)** J Tavenner, M Mendelev, J Lawson, Computational Materials Science, 218, 111929 (2023).
