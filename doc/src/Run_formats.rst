File formats used by LAMMPS
===========================

This page provides a general overview of the kinds of files and file
formats that LAMMPS is reading and writing.

Character Encoding
^^^^^^^^^^^^^^^^^^

For text files, LAMMPS uses `ASCII character encoding
<https://en.wikipedia.org/wiki/ASCII>`_ which represents the digits 0 to
9, the lower and upper case letters a to z, some common punctuation and
other symbols and a few whitespace characters including a regular "space
character", "line feed", "carriage return", "tabulator". These are all
represented by bytes with a value smaller than 128 and only 95 of those
128 values represent printable characters.  This is sufficient to represent
most English text, but misses accented characters or umlauts or Greek
symbols and more.

Modern text often uses `UTF-8 character encoding
<https://en.wikipedia.org/wiki/UTF-8>`_ instead. This is a way to
represent many more different characters as defined by the Unicode
standard.  This is compatible with ASCII, since the first 128 values are
identical with the ASCIII encoding.  It is important to note, however,
that there are Unicode characters that look similar or even identical to
ASCII characters, but have a different representation.  As a general
rule, these characters are not correctly recognized by LAMMPS.  For some
parts of LAMMPS' text processing, translation tables with known
"lookalike" characters are used that transparently substitute non-ASCII
characters with their ASCII equivalents.  Non-ASCII lookalike characters
are often used by web browsers or PDF viewers to improve the readability
of text. Thus, when using copy-n-paste to transfer text from such an
application to your input file, you may unintentionally create text that
is not exclusively using ASCII encoding and may cause errors when LAMMPS
is trying to read it.

Lines with non-printable and non-ASCII characters in text files can be
detected for example with a (Linux) command like the following:

.. code-block:: bash

   env LC_ALL=C grep -n '[^ -~]' some_file.txt

Number Formatting
^^^^^^^^^^^^^^^^^

Different countries and languages have different conventions to format
numbers.  While in some regions commas are used for fractions and points
to indicate thousand, million and so on, this is reversed in other
regions.  Modern operating systems have facilities to adjust input and
output accordingly.  The exact rules are often applied according to the
value of the ``$LANG`` environment variable (e.g. "en_US.utf8").

For the sake of simplicity of the implementation and transferability of
results, LAMMPS does not support this and instead expects numbers being
formatted in the generic or "C" locale.  The "C" locale has no
punctuation for thousand, million and so on and uses a decimal point for
fractions.  One thousand would be represented as "1000.0" and not as
"1,000.0" nor as "1.000,0".

LAMMPS also only accepts integer numbers when an integer is required,
so using "1.0" is not accepted; you have to use "1" instead.

For floating point numbers in scientific notation, the Fortran double
precision notation "1.1d3" is not accepted either; you have to use
"1100", "1100.0" or "1.1e3".

Input file
^^^^^^^^^^

A LAMMPS input file is a text file with commands. It is read
line-by-line and each line is processed *immediately*.  Before looking
for commands and executing them, there is a pre-processing step where
comments (text starting with a pound sign '#') are removed,
`${variable}` and `$(expression)` constructs are expanded or evaluated,
and lines that end in the ampersand character '&' are combined with the
next line (similar to Fortran 90 free format source code).  After the
pre-processing, lines are split into "words" and the first word must be a
:doc:`command <Commands_all>` and everything .  Below are some example lines:

.. code-block:: LAMMPS

   # full line comment

   # some global settings
   units           lj
   atom_style      atomic
   # ^^ command    ^^ argument(s)

   variable        x index 1       # may be overridden from command line with -var x <value>
   variable        xx equal 20*$x  # variable "xx" is always 20 time "x"

   lattice         fcc 0.8442

   # multi-line command, uses spacing from "lattice" command
   region          box block 0.0 ${xx} &
                             0.0 40.0  &
                             0.0 30.0
   # create simulation box and fillwith atoms according to lattice setting
   create_box      1 box
   create_atoms    1 box

   # set force field and parameters
   mass            1 1.0
   pair_style      lj/cut 2.5
   pair_coeff      1 1 1.0 1.0 2.5

   # run simulation
   fix             1 all nve
   run             1000

The pivotal command in this example input is the :doc:`create_box
command <create_box>`.  It defines the simulation system and many
parameters that go with it: units, atom style, number of atom types (and
other types) and more.  Those settings are *locked in* after the box is
created.  Commands that change these kind of settings are only allowed
**before** a simulation box is created and many other commands are only
allowed **after** the simulation box is defined (e.g. :doc:`pair_coeff
<pair_coeff>`).  Very few commands (e.g. :doc:`pair_style <pair_style>`)
may be used in either part of the input.  The :doc:`read_data
<read_data>` and :doc:`read_restart <read_restart>` commands also create
the system box and thus have a similar pivotal function.

The LAMMPS input syntax has minimal support for conditionals and loops,
but if more complex operations are required, it is recommended to use
the library interface, e.g. :doc:`from Python using the LAMMPS Python
module <Python_run>`.

There is a frequent misconception about the :doc:`if command <if>`:
this is a command for conditional execution **outside** a run or
minimization.  To trigger actions on specific conditions **during**
a run is a non-trivial operation that usually requires adopting one
of the available fix commands or creating a new one.

LAMMPS commands change the internal state and thus the order of commands
matters and reordering them can produce different results.  For example,
the region defined by the :doc:`region command <region>` in the example
above depends on the :doc:`lattice setting <lattice>` and thus its
dimensions will be different depending on the order of the two commands.

Each line must have an "end-of-line" character (line feed or carriage
return plus line feed).  Some text editors do not automatically insert
one which may cause LAMMPS to ignore the last command.  It is thus
recommended, to always have an empty line at the end of an input file.

The specific details describing how LAMMPS input is processed and parsed
are explained in :doc:`Commands_parse`.

Data file
^^^^^^^^^

A LAMMPS data file contains a description of a system suitable for
reading with the :doc:`read_data command <read_data>`.  This is commonly
used for setting up more complex and particularly molecular systems
which can be difficult to achieve with the commands :doc:`create_box
<create_box>` and :doc:`create_atoms <create_atoms>` alone.  Also, data
files can be used as a portable alternatives to a :doc:`binary restart
file <restart>`.  A restart file can be converted into a data file
from the :doc:`command line <Run_options>`.

The file is generally structured into a header section at the very
beginning of the file and multiple titled sections like "Atoms",
Masses", "Pair Coeffs", and so on.  The data file **always** starts
with a "title" line, which will be **ignored** by LAMMPS.  Omitting
the title line can lead to unexpected behavior as then a line of
the header with an actual setting may be ignored.  This is often a
line with the "atoms" keyword, which results in LAMMPS assuming that
there are no atoms in the data file and thus throwing an error on the
contents of the "Atoms" section.  The title line may contain some
keywords that can be used by external programs to convey information
about the system, that is not required and not read by LAMMPS.

Data files may contain comments, which start with the pound sign '#'.
There must be at least one blank between a valid keyword and the pound
sign.

.. code-block:: bash

   LAMMPS Title line (ignored)
   # full line comment

           10  atoms # comment
            4  atom types

    -36.840194 64.211560 xlo xhi
    -41.013691 68.385058 ylo yhi
    -29.768095 57.139462 zlo zhi

   Masses

     1 12.0110
     2 12.0110
     3 15.9990
     4  1.0080

   Pair Coeffs

     1    0.110000    3.563595    0.110000    3.563595
     2    0.080000    3.670503    0.010000    3.385415
     3    0.120000    3.029056    0.120000    2.494516
     4    0.022000    2.351973    0.022000    2.351973

   Atoms # full

         1      1       1       0.560   43.99993  58.52678  36.78550   0   0   0
         2      1       2      -0.270   45.10395  58.23499  35.86693   0   0   0
         3      1       3      -0.510   43.81519  59.54928  37.43995   0   0   0
         4      1       4       0.090   45.71714  57.34797  36.13434   0   0   0
         5      1       4       0.090   45.72261  59.13657  35.67007   0   0   0
         6      1       4       0.090   44.66624  58.09539  34.85538   0   0   0
         7      1       3      -0.470   43.28193  57.47427  36.91953   0   0   0
         8      1       4       0.070   42.07157  57.45486  37.62418   0   0   0
         9      1       1       0.510   42.19985  57.57789  39.12163   0   0   0
        10      1       1       0.510   41.88641  58.62251  39.70398   0   0   0
   #  ^^atomID ^^molID ^^type  ^^charge ^^xcoord  ^^ycoord  ^^ycoord  ^^image^^flags

Molecule file
^^^^^^^^^^^^^


Potential file
^^^^^^^^^^^^^^


Restart file
^^^^^^^^^^^^


Dump file
^^^^^^^^^

