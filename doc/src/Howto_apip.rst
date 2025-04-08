Adaptive-precision interatomic potentials (APIP)
================================================

The :ref:`APIP package <PKG-APIP>` allows to use adaptive-precision potentials
according to :ref:`(Immel) <Immel2025_1>`.
The potential energy :math:`E_i` of an atom :math:`i` of an adaptive-precision
interatomic potential is given by :ref:`(Immel) <Immel2025_1>`

.. math::

   E_i = \lambda_i E_i^\text{(fast)} + (1-\lambda_i) E_i^\text{(precise)},

whereas :math:`E_i^\text{(fast)}` is the potential energy of atom :math:`i`
according to a fast interatomic potential,
:math:`E_i^\text{(precise)}` is the potential energy according to a precise
interatomic potential and :math:`\lambda_i\in[0,1]` is the
switching parameter that decides which potential energy is used.

The currently implemented potentials are:

.. list-table::
   :header-rows: 1

   * - Fast potential
     - Precise potential
   * - :doc:`ACE <pair_pace_apip>`
     - :doc:`ACE <pair_pace_apip>`
   * - :doc:`EAM <pair_eam_apip>`
     -

In theory, any short-range potential can be used for an adaptive-precision
interatomic potential. How to implement a new (fast or precise)
adaptive-precision
potential is explained in :ref:`here <implementing_new_apip_styles>`.

To run a simulation with an adaptive-precision potential, one needs the
following components:

  #. :doc:`atom_style apip <atom_style>` so that the switching parameter :math:`\lambda_i` can be stored.
  #. A fast potential: :doc:`eam/apip <pair_eam_apip>` or :doc:`pace/apip/fast <pair_pace_apip>`.
  #. A precise potential: :doc:`pace/apip/precise <pair_pace_apip>`.
  #. :doc:`pair_style lambda_input  <pair_lambda_input>` to calculate :math:`\lambda_i^\text{input}`, from which :math:`\lambda_i` is calculated.
  #. :doc:`fix lambda <fix_lambda>` to calculate the switching parameter.
  #. :doc:`pair_style lambda/zone <pair_lambda_zone>` to calculate the spatial transition zone of the switching parameter.
  #. :doc:`pair_style hybrid/overlay <pair_hybrid>` to combine the previously mentioned pair_styles.
  #. :doc:`fix lambda_thermostat <fix_lambda_thermostat>` to conserve the energy when switching parameters change.
  #. :doc:`fix apip_atom_weight <fix_apip_atom_weight>` to approximate the load caused by every atom, as the computations of the pair_styles are only required for a subset of atoms.
  #. :doc:`fix balance <fix_balance>` to perform dynamic load balancing with the calculated load.


Example
"""""""
.. note::

   How to select the values of the parameters of an adaptive-precision
   interatomic potential is discussed in detail in :ref:`(Immel) <Immel2025_1>`.

The affected parts of a LAMMPS script can look as follows:

.. code-block:: LAMMPS

   atom_style apip
   comm_style tiled

   pair_style hybrid/overlay eam/fs/apip pace/apip/precise lambda_input/csp fcc cutoff 5.0 lambda/zone 12.0
   pair_coeff * * eam/fs/apip Cu.eam.fs Cu
   pair_coeff * * pace/apip Cu.yace Cu
   pair_coeff * * lambda_input/csp
   pair_coeff * * lambda/zone

   fix 2 all lambda 2.5 3.0 time_averaged_zone 4.0 12.0 110 110 min_delta_lambda 0.01
   fix 3 all lambda_thermostat N_rescaling 200
   fix 4 all apip_atom_weight 100 eam ace lambda_input lambda all

   variable myweight atom f_4

   fix 5 all balance 100 1.1 rcb weight var myweight

First, the :doc:`atom_style <atom_style>` and the communication style are set.

.. note::
   Note, that :doc:`comm_style <comm_style>` *tiled* is required for the style *rcb* of
   :doc:`fix balance <fix_balance>`, but not for APIP.
   However, the flexibility offered by the balancing style *rcb*, compared to the
   balancing style *shift*, is advantageous for APIP.

An adaptive-precision EAM-ACE potential, for which the switching parameter
:math:`\lambda` is calculated from the CSP is defined via
:doc:`pair_style hybrid/overlay <pair_hybrid>`.
The fixes ensure that the switching parameter is calculated, the energy conserved,
the weight for the load balancing calculated and the load-balancing itself is done.

----------

.. _implementing_new_apip_styles:

Implementing new APIP pair styles
"""""""""""""""""""""""""""""""""

One can introduce adaptive-precision to an existing pair style by modifying
the original pair style.
One should calculate the force
:math:`F_i = - \nabla_i \sum_j E_j^\text{original}` for a fast potential or
:math:`F_i = - (1-\nabla_i) \sum_j E_j^\text{original}` for a precise
potential from the original potential
energy :math:`E_j^\text{original}` to see where the switching parameter
:math:`\lambda_i` needs to be introduced in the force calculation.
The switching parameter :math:`\lambda_i` is known for all atoms :math:`i`
in force calculation routine.
One needs to introduce an abortion criterion based on :math:`\lambda_i` to
ensure that all not required calculations are skipped and compute time can
be saved.
Furthermore, one needs to provide the number of calculations and measure the
computation time.
Communication within the force calculation needs to be prevented to allow
effective load-balancing.
With communication, the load balancer cannot balance few calculations of the
precise potential on one processor with many computations of the fast
potential on another processor.

All changes in the pair_style pace/apip compared to the pair_style pace
are annotated and commented.
Thus, the pair_style pace/apip can serve as an example for the implementation
of new adaptive-precision potentials.

----------

.. _Immel2025_1:

**(Immel)** Immel, Drautz and Sutmann, J Chem Phys, 162, 114119 (2025)
