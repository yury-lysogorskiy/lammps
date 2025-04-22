.. index:: pair_style lj/pirani
.. index:: pair_style lj/pirani/omp

pair_style lj/pirani command
============================

Accelerator Variants: *lj/pirani/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/pirani cutoff

* lj/pirani = name of the pair style
* cutoff = global cutoff (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/pirani 10.0
   pair_coeff 1 1 4.0 7.0 6.0 3.5 0.0045

Description
"""""""""""

.. versionadded:: TBD

Pair style *lj/pirani* computes pairwise interactions from an Improved
Lennard-Jones (ILJ) potential according to :ref:`(Pirani) <Pirani>`:


.. math::

   x = r/R_m   \\
   n_x = \alpha*x^2 + \beta   \\
   \gamma \equiv m  \\

  V(x) = \varepsilon \cdot \left( \frac{\gamma}{ n_x - \gamma}  \left(\frac{1}{x} \right)^{n_x}
          -  \frac{n_x}{n_x - \gamma}  \left(\frac{1}{x} \right)^{\gamma} \right) \qquad r < r_c

:math:`r_c` is the cutoff.


An additional parameter, :math:`\alpha`, has been introduced in order
to be able to recover the traditional Lennard-Jones (LJ) 12-6 with an adequate
choice of parameters. With :math:`R_m \equiv r_0 = \sigma \cdot 2^{1 / 6}`,
:math:`\alpha = 0`, :math:`\beta = 12` and :math:`\gamma = 6`
it is straightforward to prove that LJ 12-6 is obtained.


This potential provides some advantages with respect to the LJ
potential and can be really useful for molecular dynamics simulations,
as one can see from :ref:`(Pirani) <Pirani>`.
It can be used for neutral-neutral (:math:`\gamma = 6`),
ion-neutral (:math:`\gamma = 4`) or ion-ion systems (:math:`\gamma = 1`).
It removes most of the issues at short- and long-range of the LJ model.


It is possible to verify that using (:math:`\alpha= 4`), (:math:`\beta= 6`)
and (:math:`\gamma = 6`), at the equilibrium distance,
the first and second derivatives of ILJ coincide with those of LJ 12-6
( and the reduced force constant amounts to the typical 72).
In this case, LJ provides a long-range coefficient with a factor of 2 compared
with the ILJ. Also, the short-range interaction is overestimated by LJ.
The ILJ potential solves both problems.


The analysis of a diverse amount of systems verified that (:math:`\alpha= 4`)
works very well. In some special cases (e.g. those involving very small
multiple charged ions) this factor may take a slightly different value.
The parameter (:math:`\beta`) codifies the hardness (polarizability) of the
interacting partners, and for neutral-neutral systems it ranges from 6 to 11.
Moreover, the modulation of (:math:`\beta`) permits to indirectly consider the
role of further interaction components (such as the charge transfer in the
perturbative limit) and mitigates the effect of some uncertainty in the data.


The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\alpha` (dimensionless)
* :math:`\beta` (dimensionless)
* :math:`\gamma` (dimensionless)
* :math:`R_m` (distance units)
* :math:`\epsilon` (energy units)
* cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Thus, coefficients for all I,J
pairs must be specified explicitly.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options are not relevant for
this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style supports the use of the *inner*, *middle*, and
*outer* keywords of the :doc:`run_style respa <run_style>` command,
meaning the pairwise forces can be partitioned by distance at different
levels of the rRESPA hierarchy. See the :doc:`run_style <run_style>`
command for details.


----------

Restrictions
""""""""""""

This pair style is only enabled if LAMMPS was built with the EXTRA-PAIR
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut <pair_lj>`

Default
"""""""

none

--------------

.. _Pirani:

**(Pirani)** F. Pirani, S. Brizi, L. Roncaratti, P. Casavecchia, D. Cappelletti and F. Vecchiocattivi,
Phys. Chem. Chem. Phys., 2008, 10, 5489-5503.
