.. index:: pair_style lambda_input
.. index:: pair_style lambda_input/csp

pair_style lambda_input command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lambda_input cutoff

* lambda_input = style name of this pair style
* cutoff = global cutoff (distance units)

pair_style lambda_input/csp command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lambda_input/csp lattice keyword args

* lambda_input/csp = style name of this pair style
* lattice = *fcc* or *bcc* or integer

  .. parsed-literal::

       *fcc* = use 12 nearest neighbors to calculate the CSP like in a perfect fcc lattice
       *bcc* = use 8 nearest neighbors to calculate the CSP like in a perfect bcc lattice
       integer = use N nearest neighbors to calculate the CSP

* zero or more keyword/args pairs may be appended
* keyword = *cutoff* or *N_buffer*

  .. parsed-literal::

       *cutoff* args = cutoff
         cutoff = distance in which neighbouring atoms are considered (> 0)
       *N_buffer* args = N_buffer
         N_buffer = numer of additional neighbours, which are included in the j-j+N/2 calculation

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lambda_input/csp fcc
   pair_style lambda_input/csp fcc cutoff 5.0
   pair_style lambda_input/csp bcc cutoff 5.0 N_buffer 2
   pair_style lambda_input/csp 14

Description
"""""""""""

This pair_styles calculates :math:`\lambda_i^\text{input}(t)`, which
is required for :doc:`fix lambda <fix_lambda>`.

The pair_style lambda_input sets :math:`\lambda_i^\text{input}(t) = 0`.

The pair_style lambda_input/csp calculates
:math:`\lambda_i^\text{input}(t) = \text{CSP}_i(t)`.
The centro-symmetry parameter (CSP) :ref:`(Kelchner) <Kelchner_2>` is described
in :doc:`compute centro/atom <compute_centro_atom>`.

The lattice argument is described in
:doc:`compute centro/atom <compute_centro_atom>` and determines
the number of neighboring atoms that are used to compute the CSP.
The *N_buffer* argument allows to include more neighboring atoms in
the calculation of the contributions from the pair j,j+N/2 to the CSP as
discussed in :ref:`(Immel) <Immel2025_6>`.

The computation of :math:`\lambda_i^\text{input}(t)` is done by this
pair_style instead of by :doc:`fix lambda <fix_lambda>`, as this computation
takes time and this pair_style can be included in the load-balancing via
:doc:`fix apip_atom_weight <fix_apip_atom_weight>`.

A code example for the calculation of the switching parameter for an adaptive-
precision potential is given in the following:
The adaptive-precision potential is created
by combining :doc:`pair_style eam/fs/apip <pair_eam_apip>`
and :doc:`pair_style pace/apip/precise <pair_pace_apip>`.
The input, from which the switching parameter is calculated, is provided
by this pair_style.
The switching parameter is calculated by :doc:`fix lambda <fix_lambda>`,
whereas the spatial
transition zone of the switching parameter is calculated by
:doc:`pair_style lambda/zone <pair_lambda_zone>`.

.. code-block:: LAMMPS

   pair_style hybrid/overlay eam/fs/apip pace/apip/precise lambda_input/csp fcc cutoff 5.0 lambda 12.0
   pair_coeff * * eam/fs/apip Cu.eam.fs Cu
   pair_coeff * * pace/apip Cu_precise.yace Cu
   pair_coeff * * lambda_input/csp
   pair_coeff * * lambda
   fix 2 all lambda 3.0 3.5 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style writes no information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands need
to be specified in an input script that reads a restart file.

This pair style does not support the use of the *inner*, *middle*,
and *outer* keywords of the :doc:`run_style respa <run_style>` command.

----------

Restrictions
""""""""""""
This fix is part of the APIP package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute centro/atom <compute_centro_atom>`,
:doc:`fix lambda <fix_lambda>`,
:doc:`fix lambda_thermostat <fix_lambda_thermostat>`,
:doc:`pair_style lambda/zone <pair_lambda_zone>`,
:doc:`pair_style eam/apip <pair_eam_apip>`,
:doc:`pair_style pace/apip  <pair_pace_apip>`,
:doc:`fix apip_atom_weight <fix_apip_atom_weight>`

Default
"""""""

N_buffer=0, cutoff=5.0

----------

.. _Kelchner_2:

**(Kelchner)** Kelchner, Plimpton, Hamilton, Phys Rev B, 58, 11085 (1998).

.. _Immel2025_6:

**(Immel)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)
