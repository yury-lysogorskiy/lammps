.. index:: fix ave/moments

fix ave/moments command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ave/moments Nevery Nrepeat Nfreq value1 value2 ... moment1 moment2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/moments = style name of this fix command
* Nevery = use input values every this many time steps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many time steps
* one or more input variables can be listed
* value = v_name

  .. parsed-literal::

       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID, I can include wildcard (see below)
       v_name = value calculated by an equal-style variable with name
       v_name[I] = value calculated by a vector-style variable with name, I can include wildcard (see below)

* one or more moments to compute can be listed
* moment = *mean* or *stddev* or *variance* or *skew* or *kurtosis*, see exact definitions below.
* zero or more keyword/arg pairs may be appended
* keyword = *start* or *delay*

  .. parsed-literal::

       *start* args = Nstart
         Nstart = invoke first after this time step
       *history* args = Ndelay
         Ndelay = keep a history of up to Ndelay invocations

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/moments 1 1000 100 v_volume mean sttdev
   fix 1 all ave/moments 1 200 1000 v_volume variance kurtosis history 10

Description
"""""""""""

Using one or more variables as input every few time steps, calculate the
moments of the underlying distribution based on the samples collected over
a time step window. The definitions of the moments calculated are given below.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an equal-style or vector-style
:doc:`variable <variable>`.  In each case, the compute, fix, or variable must
produce a global quantity, quantity, not a per-atom or local quantity.
If you wish to spatial- or time-average or histogram per-atom quantities
from a compute, fix, or variable, then see the :doc:`fix ave/chunk <fix_ave_chunk>`,
:doc:`fix ave/atom <fix_ave_atom>`, or :doc:`fix ave/histo <fix_ave_histo>` commands.
If you wish to sum a per-atom quantity into a single global quantity, see the
:doc:`compute reduce <compute_reduce>` command.

:doc:`Computes <compute>` that produce global quantities are those which
do not have the word *atom* in their style name.  Only a few
:doc:`fixes <fix>` produce global quantities.  See the doc pages for
individual fixes for info on which ones produce such values.
:doc:`Variables <variable>` of style *equal* and *vector* are the only
ones that can be used with this fix.  Variables of style *atom* cannot
be used, since they produce per-atom values.

The input values must all be scalars or vectors with a bracketed term appended,
indicating the :math:`I^\text{th}` value of the vector is used.

The result of this fix can be accessed as a vector, containing the interleaved
moments of each input in order.  The first requested moment of input 1
has index 1, the second index 2, the first of input 2 has index 3
and so on.

----------

For input values from a compute or fix or variable, the bracketed
index I can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the size of the vector,
then an asterisk with no numeric values means all indices from 1 to :math:`N`.
A leading asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to :math:`N` (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual elements of the
vector or cells of the array had been listed one by one.  For examples, see the
description of this capability in :doc:`fix ave/time <fix_ave_time>`.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and :math:`N_\text{freq}`
arguments specify on what time steps the input values will be used in order to
contribute to the average.  The final statistics are generated on
time steps that are a multiple of :math:`N_\text{freq}`\ .  The average is over
a window of up to :math:`N_\text{repeat}` quantities, computed in the preceding
portion of the simulation every :math:`N_\text{every}` time steps.

.. note::

    Contrary to some fix ave/* commands, the values of this fix are not restricted by any special relation:
    it is valid to have a window larger than :math:`N_\text{freq}` as well as the other way around.

For example, if :math:`N_\text{freq}=100` and :math:`N_\text{repeat}=5` (and
:math:`N_\text{every}=1`), then values from time steps 96, 97, 98, 99, and 100
will be used. This means some intervening time steps do not contribute to the result.
If :math:`N_\text{freq}=5` and :math:`N_\text{repeat}=10`, then values will
first be calculated on step 5 from steps 1-5, on step 10 from 1-10, on
step 15 from 5-15 and so on, forming a rolling average.

----------

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If no bracketed term is appended,
the global scalar calculated by the compute is used.  If a bracketed term is
appended, the Ith element of the global vector calculated by the compute is used.
See the discussion above for how I can be specified with a wildcard
asterisk to effectively specify multiple values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If no bracketed term is appended,
the global scalar calculated by the fix is used.  If a bracketed term is
appended, the Ith element of the global vector calculated by the fix is used.
See the discussion above for how I can be specified with a wildcard asterisk to
effectively specify multiple values.

Note that some fixes only produce their values on certain time steps,
which must be compatible with *Nevery*, else an error will result.
Users can also write code for their own fix styles and :doc:`add them to LAMMPS <Modify>`.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script. Only equal-style or vector-style
variables can be used, which both produce global values.  Vector-style
variables require a bracketed term to specify the Ith element of the
vector calculated by the variable.

Note that variables of style *equal* and *vector* define a formula
which can reference individual atom properties or thermodynamic
keywords, or they can invoke other computes, fixes, or variables when
they are evaluated, so this is a very general means of specifying
quantities to time average.

----------

The moments are output in the order requested in the arguments following the last
input.  Any number and order of moments can be specified, although it does not make
much sense to specify the same moment multiple times.  All moments are computed in
terms of corrected sample (not population) cumulants :math:`k_{1..4}` (see
:ref:`(Cramer)<Cramer1>`), the standardized moments follow :ref:`(Joanes)<Joanes1>`.

For *mean*, the arithmetic mean :math:`\bar{x} = \frac{1}{n} \sum_{i=1}^{n} x_i` is calculated.

For *variance*, the Bessel-corrected sample variance
:math:`var = k_2 = \frac{1}{n - 1} \sum_{i=1}^{n} (x_i - \bar{x})^2` is calculated.

For *stddev*, the Bessel-corrected sample standard deviation
:math:`stddev = \sqrt{k_2}` is calculated.

For *skew*, the adjusted Fisher--Pearson standardized moment
:math:`G_1 = \frac{k_3}{k_2^{3/2}} = \frac{k_3}{stddev^3}` is calculated.

For *kurtosis*, the adjusted Fisher--Pearson standardized moment
:math:`G_2 = \frac{k_4}{k_2^2}` is calculated.

----------

Fix invocation and output can be modified by optional keywords.

The *start* keyword specifies that the first invocation should be no earlier than
the step number given (but will still occur on a multiple of *Nfreq*).
The default is step 0.  Often input values can be 0.0 at time 0, so setting
*start* to a larger value can avoid including a 0.0 in a longer series.

The *history* allows keeping a record of previous results.  By default, only
the most recent invocation is accessible.

For example, this will output values which are delayed by 10 invocations,
meaning 10000 time steps:

.. code-block:: LAMMPS

   fix 1 all ave/moments 1 200 1000 v_volume mean history 10

The previous results can be accessed by additional rows on the fix output
array, containing the N-th last evaluation result.  For example, the most recent
result of the first input value would be accessed as "f_name[1][1]",
"f_name[1][4]" is the 4th most recent and so on.  Vector access is always the
same as the first array row, corresponding to the most recent result.

This fix can be used in conjunction with :doc:`fix halt <fix_halt>` to stop
a run automatically if a quantity is converged to within some limit:

.. code-block:: LAMMPS

   variable target equal etot
   fix aveg all ave/moments 1 200 1000 v_target mean stddev history 10
   variable stopcond equal "abs(f_aveg[1]-f_aveg[1][10])<f_aveg[2]"
   fix fhalt all halt 1000 v_stopcond == 1

In this example, every 1000 time steps, the average and standard deviation
of the total energy over the previous 200 time steps are calculated.  If the
difference between the most recent and 10-th most recent average is lower than
the most recent standard deviation, the run is stopped.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

This fix produces a global vector and global array which can be accessed
by various :doc:`output commands <Howto_output>`.
The values can be accessed on any time step, but may not be current.

A vector is produced with # of elements = number of moments * number of inputs.
The moments are output in the order given on fix definition.  An array is
produced having # of rows = value of *history* and # of columns = same as
vector output, using the same ordering.

Each element can be either "intensive" or "extensive", depending on whether
the values contributing to the element are "intensive" or "extensive". If a
compute or fix provides the value being time averaged, then the compute or
fix determines whether the value is intensive or extensive; see the page
for that compute or fix for further info.  Values produced by a variable
are treated as intensive.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This compute is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix ave/time <fix_ave_time>`,

Default
"""""""

The option defaults are history = 1, start = 0.

----------

.. _Cramer1:

**(Cramer)** Cramer, Mathematical Methods of Statistics, Princeton University Press (1946).

.. _Joanes1:

**(Joanes)** Joanes, Gill, The Statistician, 47, 183--189 (1998).
