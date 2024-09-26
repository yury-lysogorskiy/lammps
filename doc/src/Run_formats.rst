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
rule, these characters are not be correctly recognized by LAMMPS.  For
some parts of LAMMPS' text processing, there are translation tables with
known "lookalike" characters that will transparently substitute them
with their ASCII equivalents.  Non-ASCII lookalike characters are often
used by web browsers or PDF viewers to improve the readability of
text. Thus, when using copy-n-paste to transfer text from such an
application to your input file, you may unintentionally create text that
is not fully in ASCII encoding and may cause errors when LAMMPS is
trying to read it.

Lines with non-printable and non-ASCII characters in text files can be
detected for example with:

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
`${variable}` and `$(expression)` constructs are expanded or evaluated
and lines that end in the ampersand character '&' are combined with the
next line (similar to Fortran 90 free format source code).

The LAMMPS input syntax has minimal support for conditionals and loops,
but if more complex operations are required, it is recommended to use
the library interface, e.g. :doc:`from Python using the LAMMPS Python
module <Python_run>`.

There is a frequent misconception about the :doc:`if command <if>`:
this is a command for conditional execution **outside** a run or
minimization.  To trigger actions on specific conditions **during**
a run is a non-trivial operation that usually requires adopting one
of the available fix commands or creating a new one.

LAMMPS commands can change the internal state and thus the order of
commands matters and reordering them can produce different results.

Each line must have an "end-of-line" character (line feed or carriage
return plus line feed).  Some text editors do not automatically insert
one which may have the result that LAMMPS ignores the last command.
It is thus recommended, to always have an empty line at the end of an
input file.

The specific details describing how LAMMPS input is processed and parsed
are explained in :doc:`Commands_parse`.

Data file
^^^^^^^^^


Molecule file
^^^^^^^^^^^^^


Potential file
^^^^^^^^^^^^^^


Restart file
^^^^^^^^^^^^


Dump file
^^^^^^^^^

