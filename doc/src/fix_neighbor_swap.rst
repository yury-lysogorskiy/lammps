.. index:: fix neighbor/swap

fix neighbor/swap command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID neighbor/swap N X seed T R keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* neighbor/swap = style name of this fix command
* N = invoke this fix every N steps
* X = number of swaps to attempt every N steps
* seed = random # seed (positive integer)
* T = scaling temperature of the MC swaps (temperature units)
* R = scaling swap probability of the MC swaps (distance units)
* one or more keyword/value pairs may be appended to args
* keyword = *types* or *mu* or *ke* or *semi-grand* or *region*

  .. parsed-literal::

       *types* values = two or more atom types (1-Ntypes or type label)
       *ke* value = *no* or *yes*
         *no* = no conservation of kinetic energy after atom swaps
         *yes* = kinetic energy is conserved after atom swaps
       *region* value = region-ID
         region-ID = ID of region to use as an exchange/move volume
       *diff* values = one atom type
       *voro* values = valid voronoi compute id (compute voronoi/atom)
       *rates* values = Ntype values to conduct variable diffusion for different atom types (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   fix mc all neighbor/swap 10 160 15238 1000.0 diff 2 voro voroN
   fix myFix all neighbor/swap 100 1 12345 298.0 region my_swap_region types 5 6 voro voroN
   fix kmc all neighbor/swap 1 100 345 1.0 diff 3 rates 3 1 6 voro voroN

Description
"""""""""""

Computes MC evaluations to enable kinetic Monte Carlo (kMC) behavior
during MD simulation through only allowing neighboring atom swaps.
Neighboring atoms are selected using a Voronoi tesselation approach. This
implementation is as described in :ref:`(Tavenner) <Tavenner>`.

The fix is called every *N* timesteps and attempts *X* swaps. The system
is initialized with a random seed, using a temperature *T* for evaluating
the MC energy swaps. The distance-based probability is weighted according
to *R* which sets the radius :math:`r_0` for the weighting

.. math::

    p_{ij} = e^{(\frac{r_{ij}}{r_0})^2}

where :math:`p_{ij}` is the probability of selecting atoms :math:`i` and :math:`j` for an
evaluated swap.

The keyword *types* is submitted with two or more atom types as the value.
Atoms of the first atom type are swapped with valid neighbors of all the remaining
atom types.

The keyword *diff* is used for implementation of simulated diffusion of a given atom type
as given by *diff type*. This command selects all atom types as acceptable swap types to a
centrally selected atom of type *type*. This includes the atom type specified by the diff
keyword to account for self-diffusion hops of an atom type with itself.

Keyword *voro* is currently required, and is implemented as 

.. code-block:: LAMMPS
    voro compute-ID
    
where *compute-ID* is the ID of a corresponding Voronoi computation with neighbor list, i.e.

.. code-block:: LAMMPS

    compute compute-ID group-ID voronoi/atom neighbors yes

The group selected for computing *voro* should correspond to all the potential atoms to
be swapped at the initial step, i.e.

.. code-block:: LAMMPS
    group group-ID type 2

for using *fix neighbor/swap* with *diff 2*.

The keyword *rates* can modify the swap rate for each swapped type by values 
where the adjusted rates values are given in order of increasing atom type. 
The number of rates provided must equal the number of atom types in the simulaton.
In the third provided example above, a simulation is conducted with three atom types
where the third atom type is the one sampled for attempted swaps. All three atom
types are considered valid swaps, but atoms of type 1 will be selected three times
as often as atoms of type 2. Conversely, atoms of type 3 are six times more likely to
be selected than atoms of type two and twice as likely as atoms of type 1. 

Finally, the *region* keyword is implemented as in other atomic fixes, where
the *region region-ID* command indicates that atom swaps only be considered in the area 
given by *region-ID*. If only atoms of certain groups are expected to be in this region,
the corresponding compute voronoi command can be adjusted accordingly.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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
doc page for more info.

Voronoi compute must be enabled on build. See :doc:`compute voronoi/atom <compute_voronoi_atom>`.
A vaild voronoi command which returns neighboring atoms must be used
and referenced with the *voro* keyword.

When this fix is used with a :doc:`hybrid pair style <pair_hybrid>`
system, only swaps between atom types of the same sub-style (or
combination of sub-styles) are permitted.

If this fix is used with systems that do not have per-type masses
(e.g. atom style sphere), the ke flag must be set to off since the implemented
algorithm will not be able to re-scale velocity properly.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`compute voronoi/atom <compute_voronoi_atom>`
:doc:`delete_atoms <delete_atoms>`, :doc:`fix gcmc <fix_gcmc>`,
:doc:`fix atom/swap <fix_atom_swap>`, :doc:`fix mol/swap <fix_mol_swap>`, :doc:`fix sgcmc <fix_sgcmc>`

Default
"""""""

The option defaults are *ke* = yes, *diff* = no, *rates* = 1 for
all atom types.

----------

.. Tavenner:

**(Tavenner)** J Tavenner, M Mendelev, J Lawson, Computational Materials Science, 218, 111929 (2023).
