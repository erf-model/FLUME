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
   ``flume.plot_int = 5`` writes steps 0, 5, 10, 15 and 16. The domain is
   periodic in x and y, so the flow stays uniform: the two faces of a cell
   carry the same flux and their difference is exactly zero. Every cell of
   all five plotfiles must hold the initial state, and the step numbers and
   times must be these. With a wall in place of a periodic side the flow
   reflects and the test fails.

``uniform_channel_x``, ``uniform_channel_y``
   The uniform state at 8 m/s along a channel: periodic along the flow, walls
   on the two sides. The flow stays uniform. With the two values of
   ``flume.periodic`` exchanged it meets a wall.

``uniform_defaults``, ``uniform_gravity_default``, ``uniform_last_step_cut``, ``uniform_stop_time_roundoff``, ``uniform_max_step_first``, ``uniform_negative_u``, ``uniform_negative_v``, ``uniform_unequal_cells``, ``uniform_stop_time_zero``, ``uniform_max_step_zero``
   The uniform state with other controls: the defaults of ``flume.cfl``,
   ``flume.plot_int`` and ``flume.gravity`` (one step of
   0.5 / sqrt(9.81 x 2.5) s); CFL 1 with a stop time of 0.3125 s, where the
   third step is shortened from 0.125 s to 0.0625 s (read from the log,
   because a uniform state does not show the size of a step); ten steps of 0.1 s whose sum
   in double precision is 0.9999999999999999 s, where the tenth must end on
   the stop time of 1 s with no eleventh step; a step limit that ends the run
   before the stop time; velocities whose larger component is negative, in x
   and in y; cells of 2 m by 1 m, where each direction's rate uses its own
   cell size; and a stop time or a step limit of zero, which give one
   plotfile.

``dry_stop_time``
   A dry domain over a sloping bed, on cells of 2 m by 1 m, with a stop time
   of 5 s: one step, to 5 s.

``dry_stop_time_huge``
   A dry domain with a stop time of 10\ :sup:`308` s, on 2 ranks: one step
   to that time. No kernel runs on a dry domain, so nothing is multiplied by
   the step.

``sloping_bed_initial``
   The initial state of ``Tests/SlopingBed/inputs`` with a velocity of
   (6, -2) m/s: the discharges of plotfile 0 are the depth times that
   velocity in every wet cell and zero in every dry one.

``dam_initial_on_centre``, ``dam_initial_on_centre_y``
   The initial state of ``Tests/Basin/inputs`` with the dam on a line of
   cell centres and another velocity on its high side. Every cell of
   plotfile 0 is checked: a cell whose centre is at the dam belongs to the
   high side, in depth and in velocity, and the last place allowed for a dam
   (the last cell centre) runs.

``sloping_bed_np1``, ``sloping_bed_np2``, ``sloping_bed_np4``
   ``Tests/SlopingBed/inputs`` on 1, 2 and 4 ranks: a lake at rest that does
   not cover its bed. A plane bed at 100 m at the lower corner with slopes
   -1/32 and 1/64, under a level surface at 98.921875 m, on 50 by 30 cells of
   2 m in 8 boxes of uneven size, between walls. The bed rises through the
   surface: cells with 4 i - 2 j <= 68 are dry, and the bed of cell (17, 0)
   is exactly at the surface. After three steps every cell is checked against
   the closed form for the bed, the depth, zero discharges and the surface
   elevation. The deepest cell, with 2 m of water, sets the step (0.125 s)
   and lies in one box, and some boxes are dry. On 2 and 4 ranks the ranks
   agree on the step only through the reduction of the largest rate over the
   ranks: without it a rank that holds only dry boxes finds no wave and stops
   the run.

``lake_at_rest_np1``, ``lake_at_rest_np2``, ``lake_at_rest_np4``, ``lake_at_rest_islands``
   ``Tests/LakeAtRest/inputs``: a level surface at 13 m over a rough bed (a
   plane moved by up to 0.5 m in every cell), between walls, for 50 steps.
   Depth and discharges must be the initial ones within 16 machine epsilons
   of 13 m. ``Tests/check_plotfile.py`` computes the same rough bed as the
   solver. Without the bed-slope term of the hydrostatic reconstruction the
   discharges reach 1 m\ :sup:`2`/s. ``lake_at_rest_islands`` lowers the
   surface to 10.2 m, so that the high parts of the bed are dry and wet and
   dry cells lie side by side over the whole domain.

``dam_break_x_np1``, ``dam_break_x_np2``, ``dam_break_x_np4``, ``dam_break_y``
   ``Tests/DamBreak/inputs``: still water 2 m and 1 m deep on a flat bed on
   the two sides of a dam, 400 cells of 0.25 m along the flow, compared at
   5 s with the solution of Stoker (1957) by ``Tests/check_dam_break.py``.
   The constant state between the rarefaction and the shock is 1.4538 m deep
   and moves at 1.3058 m/s, and the shock moves at 4.1831 m/s. The check
   bounds the L1 error of the depth and of the discharge (2 %; the run gives
   0.8 % and 1.6 %), the constant state (0.5 %; 0.2 %) and the place of the
   shock (1 cell; 0.2), requires that nothing varies along the dam, and
   requires that the sum of the depth has not changed by more than 16
   machine epsilons of itself. ``dam_break_y`` runs the dam along y with the
   deep water on the high side, so the shock runs towards low y.

``dam_break_dry_x``, ``dam_break_dry_y``
   The same deck with no water on the low side of the dam, compared at 5 s
   with the solution of Ritter (1892), whose front has run 44.3 m and has not
   reached the wall. These runs use the dry-bed wave speeds of the flux. The
   check bounds the L1 error of the depth, as a fraction of 2 m, and of the
   discharge, as a fraction of 2 m times 4.43 m/s (1 %; the run gives 0.46 %
   and 0.39 %), requires that no depth is negative and that mass is kept as
   above, and requires that the depth beyond the exact front stays below
   10\ :sup:`-11` of 2 m (the run gives 10\ :sup:`-13` m: the scheme carries a
   trace of water one cell per step). The thin front of the run lags the exact
   one, as is usual at first order: the depth passes 1 mm about 7 m behind it.

``dam_break_dry_first_step_x``, ``dam_break_dry_first_step_y``
   The first step of the two runs above, which has a closed form. With still
   water of depth :math:`h` and celerity :math:`c = \sqrt{g h}` on one side
   of a face and a dry bed on the other, the dry-bed wave speeds are
   :math:`-c` and :math:`2c`, and the flux is :math:`2 c h / 3` of mass and
   :math:`g h^2 / 3` of momentum. The step is :math:`\mathrm{cfl} / (c /
   \Delta x + c / \Delta y)`, so with ``flume.cfl`` = 0.5, square cells and
   2 m of water the first dry cell ends 1/3 m deep and the last wet cell
   5/3 m deep, with discharges :math:`c h / 12` and :math:`c h / 24`; every
   other cell is unchanged. The
   check allows 64 machine epsilons. After the first step the bed ahead of the
   front is thinly wet and these wave speeds have little effect, so the 5 s
   runs do not test them; this one does.

``shear_x_np1``, ``shear_x_np4``, ``shear_x_reversed``, ``shear_y``
   ``Tests/Shear/inputs``: a uniform flow at 2 m/s (Froude number 0.5, so
   each face blends the states of its two sides) in a periodic domain
   carries a jump of the velocity across the flow, checked at 5 s by
   ``Tests/check_shear.py``. Depth and discharge along the flow must stay
   uniform within 4 machine epsilons. The velocity across the flow must stay
   between its two initial values, its sum must not change, and its two
   jumps must be within one cell of where the flow has carried them. One of
   the two jumps starts at the periodic side. ``shear_x_reversed`` and
   ``shear_y`` run the flow towards the low coordinate, so the upwind side of
   a face is its high side. The tests fail if the momentum along a face is
   not taken from the upwind side, whichever side that is.

``basin_mass_np1``, ``basin_mass_np2``, ``basin_mass_np4``
   ``Tests/Basin/inputs``: a dam break with moving water over a rough,
   sloping bed between four walls, for 20 s, in which the waves cross the
   basin and come back. The sum of the depth must not change by more than 16
   machine epsilons of itself, and the depth must have changed by at least
   0.1 m somewhere. A wall that does not reverse the normal discharge lets
   water through and fails the test.

``abort_*``
   At least one test per check in :doc:`Inputs`: missing keys, wrong numbers
   of values, values outside the range (including integers that would wrap),
   corners that are not finite, not ordered, too large, or too close for the
   cell size, an unknown key in the ``flume`` namespace and an unused key in
   another namespace; the ranges of the time-loop, boundary, bed, dam and
   initial-state keys, the keys a dam needs, and a misspelled dam key, which
   is reported as unknown before any key that needs it; a dry domain with a
   step limit and no stop time; finite inputs that give a bed, an initial
   state or a wave speed that is not finite; a depth that goes negative
   with ``flume.cfl`` = 0.8 (a column of wet cells between dry ones); and a
   time step that is zero or not finite. The wrong-type values that AMReX's parser
   stops are not tested.

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
