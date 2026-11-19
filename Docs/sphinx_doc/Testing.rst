Testing
=======

Running the tests
-----------------

.. code-block:: bash

   ctest --test-dir build-rel --output-on-failure

Every test runs in its own directory, ``<build>/Tests/run/<name>``, which is
removed and made again at the start of the test. Nothing is written into the
source tree. ``Tests/run_case.py`` runs the executable and then the check; the
test fails if either fails.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Label
     - Tests
   * - ``regression``
     - Runs of the executable whose output is checked against pinned numbers.
   * - ``abort``
     - Runs that must stop with a non-zero exit code and a message naming the
       input.

Select a label with ``ctest -L <label>``.

Current tests
-------------

``skeleton_np1``, ``skeleton_np2``, ``skeleton_np4``
   ``Tests/Skeleton/inputs`` on 1, 2 and 4 MPI ranks: 48 by 32 cells of 2 m on
   [-10, 86] by [5, 69] m, in 6 boxes of at most 16 cells, written to
   ``skel00000``. Every grid and output key has a value that is not its default.
   The run has no water and no stop, so it takes no step.
   ``Tests/check_plotfile.py`` reads the plotfile's headers and its binary
   cell data and checks the variable names, the step and time, the domain, the
   cell size, that the boxes cover every cell once, that the boxes were
   written to as many data files as the run has ranks, and every cell of every
   variable against a closed form (here: bed at 0, dry). Without MPI only the
   1-rank test is registered.

   The check compares the cell size with a relative tolerance of 1e-12 and the
   corners with 1e-12 of the domain extent. Cell values are compared within 4
   machine epsilons of the precision the data are written in, times the larger
   of the bed and surface elevations. The pinned numbers are exact in single
   and double precision. The check also fails if the run directory holds a
   plotfile of a step the test does not list.

``skeleton_defaults``
   ``Tests/Skeleton/inputs_defaults`` gives only the required keys: 128 by 65
   cells. The check pins the defaults: lower corner at the origin, 4 boxes
   (``flume.max_grid_size = 64``), the plotfile name ``plt00000``, a bed at
   0 and a dry domain.

``uniform_np1``, ``uniform_np2``, ``uniform_np4``
   ``Tests/Uniform/inputs`` on 1, 2 and 4 ranks: water 2 m deep on a flat bed
   moving at (6, -2) m/s, with :math:`g` = 8 m/s\ :sup:`2` on 2 m cells. The
   wave speed is 4 m/s, the rate (6 + 4) / 2 + (2 + 4) / 2 = 8 1/s and the
   step 0.5 / 8 = 0.0625 s. Sixteen steps reach the stop time of 1 s, and
   ``flume.plot_int = 5`` writes steps 0, 5, 10, 15 and 16. Every cell of all
   five plotfiles must hold the initial state, and the step numbers and times
   must be these.

``uniform_defaults``, ``uniform_gravity_default``, ``uniform_last_step_cut``, ``uniform_stop_time_roundoff``, ``uniform_max_step_first``, ``uniform_negative_u``, ``uniform_negative_v``, ``uniform_unequal_cells``, ``uniform_stop_time_zero``, ``uniform_max_step_zero``
   The uniform state with other controls: the defaults of ``flume.cfl``,
   ``flume.plot_int`` and ``flume.gravity`` (one step of
   0.5 / sqrt(9.81 x 2.5) s); CFL 1 with a stop time of 0.3125 s, where the
   third step is shortened from 0.125 s to 0.0625 s (read from the log,
   because the state does not use the step yet); ten steps of 0.1 s whose sum
   in double precision is 0.9999999999999999 s, where the tenth must end on
   the stop time of 1 s with no eleventh step; a step limit that ends the run
   before the stop time; velocities whose larger component is negative, in x
   and in y; cells of 2 m by 1 m, where each direction's rate uses its own
   cell size; and a stop time or a step limit of zero, which give one
   plotfile.

``dry_stop_time``
   A dry domain over a sloping bed, on cells of 2 m by 1 m, with a stop time
   of 5 s: one step, to 5 s.

``sloping_bed_np1``, ``sloping_bed_np2``, ``sloping_bed_np4``
   ``Tests/SlopingBed/inputs`` on 1, 2 and 4 ranks: a plane bed at 100 m at
   the lower corner with slopes -1/32 and 1/64, under a level surface at
   98.921875 m, on 50 by 30 cells of 2 m in 8 boxes of uneven size. The bed
   rises through the surface: cells with 4 i - 2 j <= 68 are dry, and the bed
   of cell (17, 0) is exactly at the surface. Every cell is checked against
   the closed form for the bed, the depth, the discharges and the surface
   elevation. The deepest cell, with 2 m of water, sets the step (0.0625 s)
   and lies in one box, and some boxes are dry. On 2 and 4 ranks the ranks
   agree on the step only through the reduction of the largest rate over the
   ranks: without it a rank that holds only dry boxes finds no wave and stops
   the run.

``abort_*``
   At least one test per check in :doc:`Inputs`: missing keys, wrong numbers of values,
   values outside the range (including integers that would wrap), corners that
   are not finite, not ordered, too large, or too close for the cell size, an
   unknown key in the ``flume`` namespace and an unused key in another
   namespace; the ranges of the time-loop, bed and initial-state keys; a dry
   domain with a step limit and no stop time; finite inputs that give a
   bed, an initial state or a wave speed that is not finite; and a time step
   that is zero or not finite. The wrong-type values that AMReX's parser
   stops are not tested, and neither is a negative depth, which no input can
   give until the flux kernel exists.

Adding a test
-------------

``Tests/CMakeLists.txt`` has two helpers.

.. code-block:: cmake

   flume_add_regtest(<name> NPROCS <n> INPUTS <file> [ARGS <key=value>...]
                     [OUTPUT <regular expression>]
                     CHECK <command and its arguments>)

   flume_add_abort_test(<name> EXPECT <regular expression> ARGS <key=value>...)

The check command runs in the test's directory. With ``OUTPUT``, the output of
the run must also match the regular expression. A new test is first shown to
fail against the defect it guards or against a one-line change to the code,
and only then committed.

Continuous integration
----------------------

Three workflows run on every pull request to ``development``.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Workflow
     - Jobs
   * - ``ci.yml``
     - GCC and Clang with MPI, GCC with OpenMP, and GCC Debug with AMReX's
       testing mode: build with warnings as errors, then every test. No MPI
       with single precision: build and the 1-rank tests. CUDA, HIP and SYCL:
       compile only, because the hosted runners have no GPU.
   * - ``style.yml``
     - ``codespell``, no tabs, no trailing whitespace, and ``actionlint`` on
       the workflow files.
   * - ``docs.yml``
     - This documentation, with warnings as errors.
