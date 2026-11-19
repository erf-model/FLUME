Inputs
======

FLUME reads an AMReX inputs file, given as the first argument. Keys can also be
given or overridden on the command line:

.. code-block:: bash

   mpiexec -n 4 ./flume inputs flume.max_grid_size=32

All of FLUME's keys are in the ``flume`` namespace. Every key has a range check
that stops the run with a message naming the key; no value is clamped. A key in
the ``flume`` namespace that FLUME does not know also stops the run at start-up,
so a misspelled key is not ignored.

A key in any other namespace that nothing reads (``flum.n_cell``, or an AMReX
key FLUME does not use, such as ``geometry.prob_lo``) is listed and stops the
run when it ends, after the output is written. FLUME turns on AMReX's
``amrex.abort_on_unused_inputs`` for this; ``amrex.abort_on_unused_inputs = 0``
turns it off.

A value of the wrong type (text for a number, ``1.5`` for an integer) is
stopped by AMReX's parser. Its message shows the value and may not name the
key. The parser accepts expressions, so ``1e3`` is a valid integer and ``16/2``
a valid number.

.. list-table::
   :header-rows: 1
   :widths: 24 12 14 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.n_cell``
     - 2 integers
     - required
     - Number of cells in x and y. Each must be from 1 to 1073741824
       (2\ :sup:`30`).
   * - ``flume.prob_lo``
     - 2 reals
     - ``0 0``
     - Lower corner of the domain in x and y, in metres. Must be finite.
   * - ``flume.prob_hi``
     - 2 reals
     - required
     - Upper corner of the domain in x and y, in metres. Must be finite and
       greater than ``flume.prob_lo`` in both directions. See the two checks
       below the table.
   * - ``flume.max_grid_size``
     - integer
     - ``64``
     - Largest box edge, in cells, in x and y. Must be from 1 to 2147483647.
   * - ``flume.plot_file``
     - string
     - ``plt``
     - Plotfile name prefix. Must not be empty. The plotfile of step 0 is
       ``<prefix>00000``.

The keys that take two values stop the run when given any other number of
values, and so do the keys that take one.

Two more checks use ``flume.prob_lo``, ``flume.prob_hi`` and ``flume.n_cell``
together, in each direction:

* The corners must be at most one eighth of the largest floating-point number
  in magnitude (2.2e307 in double precision, 4.2e37 in single).
* The cell size, ``(prob_hi - prob_lo) / n_cell``, must be a normal
  floating-point number and at least machine epsilon times the larger corner
  in magnitude. Below that, neighbouring cell edges are the same number. In
  double precision the limit is about 1e-9 m at coordinates of 4.4e6 m. In
  single precision it is about 0.5 m there.

The integer keys have no joint limit on the number of boxes. A very large
``flume.n_cell`` with a small ``flume.max_grid_size`` can ask for more boxes
than memory holds.

Single precision and large coordinates
--------------------------------------

The corners are stored in the precision of the build. In single precision a
corner at UTM-size coordinates is rounded to the nearest representable number,
which is 0.5 m apart at 4.4e6 m: ``4400123.37`` becomes ``4400123.5``. The cell
size then differs from the one requested. Use double precision, or corners near
the origin, when this matters.

Output
------

The run writes one AMReX plotfile, ``<prefix>00000``, with one variable. If a
plotfile of that name exists, AMReX first renames it to
``<prefix>00000.old.<number>``.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Variable
     - Unit
     - Meaning
   * - ``depth``
     - m
     - Water depth. Zero everywhere until the solver exists.

The plotfile is three-dimensional with one cell in z. Its z extent, 0 to 1 m,
carries no meaning.
