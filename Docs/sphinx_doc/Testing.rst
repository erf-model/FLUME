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
   A dry domain with a stop time of 10\ :sup:`308` s (2\ :sup:`127` s in
   single precision), on 2 ranks: one step
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

``lake_at_rest_islands_dry_depth``
   The islands with ``flume.dry_depth`` = 0.05 m, which many cells along the
   shores are below. The same bounds hold: the dry depth acts on the depth a
   cell presents to a face, which is the same on the two sides of a face
   under a level surface.

``lake_at_rest_dry_depth_at_a_face``
   The same lake near elevation zero (bed 0.3 m at the corner, surface
   0.7 m). There the depths that the two cells of a face present differ in
   the last digit on 468 of the 3000 faces. ``flume.dry_depth`` is set to one
   of those depths (0.94279296875 m, in double precision), so one side of
   that face is above it and the other is not. The cut is decided once for
   the face, and the lake stays at rest (discharges within 64 machine
   epsilons of the scale; the depths that differ in the last digit leave
   16). Decided for each side alone, the lake moves by 0.34 m in 50 steps.
   In single precision the value is not at a face and the test only repeats
   the lake at rest.

``draining_column_cfl_0.8``, ``draining_column_cfl_1.0``
   ``Tests/DamBreak/inputs`` with the dam 0.3 m from a periodic side: one
   column of cells holds 2 m of water, with a dry bed on both sides. Without
   the outflow limit the column would lose 4/3 of ``flume.cfl`` times its
   depth in the first step, more than it holds for both values. The runs
   reach 5 s (about 130 steps); the run stops on a negative depth, the sum of
   the depth must not change by more than 16 machine epsilons of itself, and
   no cell at or below ``flume.dry_depth`` may hold a discharge.

``draining_column_first_step``
   The first step of the column at ``flume.cfl`` = 1, which has a closed
   form, checked by ``Tests/check_cells.py`` within 64 machine epsilons. The
   column would lose :math:`4 n h / 3` with :math:`n` = 0.999, so the fluxes
   on its two faces are scaled by :math:`3 / (4 n)`: the column ends empty,
   and each neighbour holds half its water (1 m) and a discharge of
   :math:`c h / 4` = 2.2147 m\ :sup:`2`/s away from it. Every other cell is
   unchanged. The test fails if the momentum flux is not scaled with the
   mass flux, or if the factor is taken from the wrong cell.

``film_passes_no_water``
   ``Tests/Film/inputs``: a tilted surface over a flat bed on 8 x 4 cells.
   Along x there are two dry cells, a
   cell with a film of 5 x 10\ :sup:`-7` m (below the default
   ``flume.dry_depth``), then water of 0.5 m and more. After one step the two
   dry cells hold exactly nothing: the film presents no depth to its faces.
   With the cut at zero the second cell would receive
   4 x 10\ :sup:`-11` m.

``film_behind_a_bed_step``
   The same with a bed that falls by 0.0009995 m per cell along x. The third
   cell then holds 1 mm of water, above the tolerance, of which only
   5 x 10\ :sup:`-7` m stands above the bed of the second cell. The two dry
   cells hold exactly nothing after the step: the tolerance acts on the
   depth at the face, not on the depth of the cell.

``draining_column_round_off``
   A column 9.987529 m deep on a slightly rough bed at ``flume.cfl`` =
   0.9332, to 1 s. In double precision the first step leaves two cells that
   the outflow limit empties a rounding error below zero. The step sets such
   a depth to zero; without that the run stops at step 1. The deck was found
   by a search over 60 random decks, 7 of which need this.

``thin_film_time_step``
   Still water 1 m deep beside a film of 10\ :sup:`-8` m that starts at
   10\ :sup:`4` m/s, to 0.03125 s. The film is below the default
   ``flume.dry_depth``: the deep water sets the step (0.019955 s, so 2
   steps; 1130 with ``flume.dry_depth`` = 0), and after the run the film cells hold no
   discharge and the mass is unchanged. This pins the default of
   ``flume.dry_depth``.

``bowl_planar_half``, ``bowl_planar_half_fine``, ``bowl_planar``
   ``Tests/Bowl/inputs``: water in a paraboloid bowl
   (``flume.bed.bowl_depth`` :math:`h_0` = 0.1 m, ``flume.bed.bowl_radius``
   :math:`a` = 1 m) in a 4 m square, with a shoreline that moves over the dry
   bed, checked by ``Tests/check_bowl.py`` against Thacker's solutions as
   written by Delestre et al. (2013, section 4.2.2). With :math:`X`,
   :math:`Y` measured from the centre, :math:`r^2 = X^2 + Y^2` and the bed
   :math:`z = -h_0 (1 - r^2 / a^2)`, the planar solution is a plane surface
   that turns about the axis of the bowl:

   .. math::

      h = \max\left(\frac{\eta h_0}{a^2}\left(2 X \cos\omega t
          + 2 Y \sin\omega t - \eta\right) - z,\; 0\right) , \qquad
      (u, v) = \eta\,\omega\,(-\sin\omega t,\; \cos\omega t) , \qquad
      \omega = \frac{\sqrt{2 g h_0}}{a} ,

   with :math:`\eta` = 0.5 m and a period of 4.4857 s. The check compares the
   bed and the first plotfile with the closed form (8 machine epsilons),
   requires that no depth is negative, that the mass is kept (64 epsilons)
   and that dry cells hold no discharge, and bounds the L1 errors of the
   depth and of the discharge over the cells that hold water at the start or
   at the end, in the run or in the solution,
   divided by the mean initial depth (times :math:`\sqrt{g h_0}` for the
   discharge). The first-order scheme damps this motion (the second-order
   runs are below):

   .. list-table::
      :header-rows: 1
      :widths: 30 20 15 15 20

      * - Run
        - Cells
        - Depth
        - Discharge
        - Bounds
      * - half a period
        - 100 x 100
        - 13.5 %
        - 19.7 %
        - 15 %, 22 %
      * - half a period
        - 200 x 200
        - 7.3 %
        - 10.8 %
        - 8.5 %, 12.5 %
      * - three periods
        - 100 x 100
        - 57 %
        - 77 %
        - 62 %, 85 %

   The error halves when the cell size halves (6.4 % of the depth after one
   period on 400 x 400 cells; 12.3 % on 200 x 200), and the bounds of the
   fine run are below the errors of the coarse one. The three-period run is
   a test of the shoreline over 1798 steps (at most 2000): depth
   non-negative and mass kept. Its error bound and its step limit only guard
   against a change (``thin_film_time_step`` is the test of the time step:
   in the bowl the films are slow, and the step count is the same with
   ``flume.dry_depth = 0``).

``bowl_radial_half``, ``bowl_radial_half_fine``, ``bowl_radial``
   The same bowl with Thacker's radially symmetric solution: a paraboloid
   surface that rises and falls, starting at rest
   (``flume.init.surface_curvature``). With
   :math:`A = (a^2 - r_0^2) / (a^2 + r_0^2)`, :math:`r_0` = 0.8 m and
   :math:`D = 1 - A \cos\omega t`,

   .. math::

      h = \max\left(h_0\left[\frac{\sqrt{1 - A^2}}{D} - 1
          - \frac{r^2}{a^2}\left(\frac{1 - A^2}{D^2} - 1\right)\right] - z,\; 0\right) ,
      \qquad
      (u, v) = \frac{\omega A \sin\omega t}{2 D}\,(X, Y) , \qquad
      \omega = \frac{\sqrt{8 g h_0}}{a} ,

   with a period of 2.2429 s. Errors as above: 4.2 % (depth) and 2.8 %
   (discharge) after half a period on 100 x 100 cells (bounds 4.8 % and
   3.3 %), 2.2 % and 1.8 % on 200 x 200 (2.6 % and 2.1 %), and 20.7 % and
   4.5 % after three periods on 100 x 100 cells in 668 steps (23 % and
   5.5 %; at most 750 steps).

``bowl_mound_initial``
   A mound (a negative ``flume.bed.bowl_depth``) under a tilted surface on
   cells of 0.04 m by 0.05 m away from the origin, not advanced: the bed and
   the initial state against the closed form a quarter of a period on, when
   the surface is tilted along y.

``dam_break_x_np1``, ``dam_break_x_np2``, ``dam_break_x_np4``, ``dam_break_y``
   ``Tests/DamBreak/inputs``: still water 2 m and 1 m deep on a flat bed on
   the two sides of a dam, 400 cells of 0.25 m along the flow, compared at
   5 s with the exact solution by ``Tests/check_dam_break.py`` (the dam break
   on a wet domain of Delestre et al., 2013, section 4.1.1). The solution
   is a rarefaction into the deep water, a constant state and a shock into
   the shallow water. With :math:`h_1` and :math:`h_0` the deep and shallow
   depths, :math:`c = \sqrt{g h}`, and :math:`h_m`, :math:`u_m` the constant
   state, :math:`u + 2c` keeps its value across the rarefaction, and mass and
   momentum are kept across the shock, which moves at :math:`S`:

   .. math::

      u_m = 2\,(c_1 - c_m) = (h_m - h_0)\sqrt{\frac{g\,(h_m + h_0)}{2\,h_m h_0}} ,
      \qquad S = \frac{h_m u_m}{h_m - h_0} .

   The script solves the first relation for :math:`h_m` by bisection. Here
   the constant state is 1.4538 m deep and moves at 1.3058 m/s, and the
   shock moves at 4.1831 m/s. The check
   bounds the L1 error of the depth and of the discharge (2 %; the run gives
   0.8 % and 1.6 %), the constant state (0.5 %; 0.2 %) and the place of the
   shock (1 cell; 0.2), requires that nothing varies along the dam, and
   requires that the sum of the depth has not changed by more than 16
   machine epsilons of itself. ``dam_break_y`` runs the dam along y with the
   deep water on the high side, so the shock runs towards low y.

``dam_break_dry_x``, ``dam_break_dry_y``
   The same deck with no water on the low side of the dam, compared at 5 s
   with Ritter's solution (the dam break on a dry domain of Delestre et al.,
   2013, section 4.1.2), whose front has run 44.3 m and has not
   reached the wall. These runs use the dry-bed wave speeds of the flux. The
   check bounds the L1 error of the depth, as a fraction of 2 m, and of the
   discharge, as a fraction of 2 m times 4.43 m/s (1 %; the run gives 0.46 %
   and 0.39 %), requires that no depth is negative and that mass is kept as
   above, and requires that the depth beyond the exact front is zero: a cell
   no deeper than ``flume.dry_depth`` passes no water on. The thin front of the run lags the exact
   one, as is usual at first order: the depth passes 1 mm about 7 m behind it.

``dam_break_dry_x_no_tolerance``
   ``dam_break_dry_x`` with ``flume.dry_depth = 0``, the scheme without a
   tolerance. The errors are the same (0.46 %, 0.39 %). The scheme then
   carries a trace of water one cell per step: the depth beyond the exact
   front must be above zero and below 10\ :sup:`-11` of 2 m (the run gives
   10\ :sup:`-13` m). With ``dam_break_dry_x`` this shows that the value of
   the input reaches the fluxes.

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
   check allows 64 machine epsilons. The 5 s runs bound errors only; this one pins the
   wave speeds.

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

``wave_x_first_order``, ``wave_x_second_order``, ``wave_x_second_order_fine``, ``wave_x_second_order_finest``, ``wave_y_second_order``, ``wave_y_second_order_fine``, ``wave_x_second_order_wide_cells``, ``wave_x_second_order_one_cell_wide``
   ``Tests/Wave/inputs``: the rate of convergence on a smooth solution. One
   sine wave of 0.1 m on water 1 m deep over a flat bed, 10 m long and
   periodic, in water that moves at (0.5, -0.3) m/s
   (``flume.init.wave_amplitude``, :doc:`Model`), compared at 1 s with the
   exact solution by ``Tests/check_wave.py``.

   Over a flat bed and along one direction the shallow-water equations can
   be written, with :math:`c = \sqrt{g h}`, as

   .. math::

      \left(\frac{\partial}{\partial t} + (u \pm c)\frac{\partial}{\partial x}\right)
      (u \pm 2c) = 0 .

   The initial state has :math:`u - 2c = u_0 - 2 c_H` everywhere
   (:math:`c_H = \sqrt{g H}`), and the equation with the lower signs keeps it
   so. The other equation then says that :math:`u + 2c`, and so :math:`u`
   and :math:`c` each, do not change along a path that moves at
   :math:`u + c = u_0 + 3c - 2 c_H`. This speed is constant on the path, so
   the path is a straight line: the depth at :math:`x` at the time :math:`t`
   is the initial depth at :math:`x_0`, where

   .. math::

      x_0 + \bigl(u_0 + 3 c(x_0) - 2 c_H\bigr)\, t = x .

   The script solves this for :math:`x_0` by bisection. Deeper water moves
   faster, so the front of the wave steepens, and the solution is smooth
   until two paths meet, at :math:`t = -1 / \min\,(3\, \partial c / \partial x_0)`.
   This is no earlier than
   :math:`L \sqrt{H - |A|} / (3 \pi |A| \sqrt{g})` = 3.2 s here, with
   :math:`L` the length of the wave. The check stops if the time of the run
   is not below it.

   The check compares the first plotfile with the wave (64 machine
   epsilons), requires that the state does not vary across the wave and that
   the mass is kept (16 epsilons), and bounds the mean over the cells of the
   error of the depth, divided by :math:`|A|`, and of the discharge along the
   wave, divided by :math:`|A|\, c_H`. The cells are square. The runs give:

   .. list-table::
      :header-rows: 1
      :widths: 28 14 14 14 14 16

      * - Run
        - Cells along
        - Steps
        - Depth
        - Discharge
        - Bounds
      * - order 1, x
        - 100
        - 154
        - 0.0353
        - 0.0408
        - 0.045, 0.052
      * - order 2, x
        - 100
        - 154
        - 0.00350
        - 0.00413
        - 0.0044, 0.0052
      * - order 2, x
        - 200
        - 307
        - 0.00094
        - 0.00111
        - 0.0012, 0.0014
      * - order 2, x
        - 400
        - 614
        - 0.00025
        - 0.00030
        - 0.00032, 0.00038
      * - order 2, y
        - 100
        - 142
        - 0.00271
        - 0.00252
        - 0.0034, 0.0032
      * - order 2, y
        - 200
        - 283
        - 0.00073
        - 0.00068
        - 0.00092, 0.00086
      * - order 2, x, cells 250 m wide, :math:`A` = -0.1 m
        - 100
        - 82
        - 0.00364
        - 0.00428
        - 0.0046, 0.0054
      * - order 2, x, one cell wide
        - 100
        - 154
        - 0.00350
        - 0.00413
        - 0.0044, 0.0052

   The run one cell wide has a periodic direction narrower than its two
   ghost layers. The error of the second-order runs falls by 3.7 each time the cell size
   is halved, a rate of 1.9. It is below 2 because the minmod limiter sets
   the differences to zero at the highest and at the lowest point of the wave (:doc:`Model`).
   The first-order run is ten times further from the solution on the same
   cells; its test also has a lower bound (0.025), so it fails if the
   default scheme changes. The bounds of each second-order run are below
   half the errors of the run on cells twice as large, and below the
   first-order errors, so a scheme that is first order in space or in time
   fails them. On the cells 250 m wide the direction of the wave alone sets
   the time step, at ``flume.cfl = 0.5``, the largest value of the second
   order: at ``flume.cfl = 1``, which the reader stops, the errors on these
   cells are 0.0153 and 0.0181 (:doc:`Model`).

``*_second_order*``
   The tests above repeated with ``flume.order = 2``, on the same decks.

   * ``uniform_second_order_np1``, ``_np4``: the uniform state is unchanged
     in every plotfile.
   * ``sloping_bed_second_order_np1``, ``_np2``, ``_np4``,
     ``lake_at_rest_second_order_np1``, ``_np2``, ``_np4``,
     ``lake_at_rest_islands_second_order``,
     ``lake_at_rest_islands_dry_depth_second_order``,
     ``lake_at_rest_dry_depth_at_a_face_second_order``: water at rest stays
     at rest, with the tolerances of the first-order tests (16 machine
     epsilons of the surface elevation; 64 with the dry depth at a face).
   * ``dam_break_x_second_order_np1``, ``_np4``, ``dam_break_y_second_order``:
     L1 errors of 0.24 % (depth) and 0.50 % (discharge), against 0.80 % and
     1.59 % at first order; constant state within 0.03 % (0.19 %); shock
     0.1 cell from its place. Bounds: 0.65 %, 0.1 % and 1 cell.
   * ``dam_break_dry_x_second_order``, ``dam_break_dry_y_second_order``:
     0.13 % and 0.12 % against Ritter's solution (0.46 % and 0.39 %), bound
     0.2 %; no negative depth and no water beyond the exact front.
   * ``draining_column_second_order``: no negative
     depth and the same mass at 5 s.
   * ``film_passes_no_water_second_order``,
     ``film_behind_a_bed_step_second_order``: after one step (two stages)
     the first cell holds exactly nothing. The film takes water in the first
     stage, so the second cell is reached in the second.
   * ``shear_x_second_order_np1``, ``_np4``,
     ``shear_x_reversed_second_order``, ``shear_y_second_order``: the bounds
     of the first-order tests; the jumps are 0.04 cell from their place
     (0.13 at first order). In addition at most 10 cells may hold a velocity
     across the flow that is more than a tenth of the jump away from both of
     its two values (``--width-max``): the runs give 8 cells for the two
     layers together, and 16 at first order or when the velocity across the
     direction of a face is not reconstructed.
   * ``dam_break_x_second_order_one_cell_wide``: the dam break between
     walls one cell apart, where the second ghost layer outside each wall
     mirrors the first layer outside the other; every value starts as a
     signalling NaN.
   * ``thin_film_time_step_second_order``: the checks of the first-order
     test; the dry rule acts in each of the two stages.
   * ``basin_mass_second_order_np1``, ``_np2``, ``_np4``: the mass is kept
     to 16 machine epsilons between four walls.
   * ``basin_mass_second_order_thin_boxes``,
     ``basin_mass_second_order_one_cell_boxes``: the same check with boxes
     two cells wide on 41 x 25 cells, and one cell wide on 20 x 12. A box
     that ends one cell short of a wall holds the first layer outside the
     wall in its second ghost layer, and its last face reads it. The runs
     start with every value set to a signalling NaN
     (``amrex.init_snan = 1``); they fail (the mass changes by 0.9 m\ :sup:`3`
     of 1267) if only the boxes that touch a wall are filled.
   * ``bowl_planar_half_second_order``, ``bowl_planar_half_fine_second_order``,
     ``bowl_planar_second_order``, ``bowl_radial_half_second_order``,
     ``bowl_radial_half_fine_second_order``, ``bowl_radial_second_order``:
     Thacker's solutions, with the errors below.

   .. list-table::
      :header-rows: 1
      :widths: 30 14 14 14 14 14

      * - Run
        - Cells
        - Depth
        - Discharge
        - Depth, order 1
        - Discharge, order 1
      * - planar, half a period
        - 100 x 100
        - 1.71 %
        - 2.12 %
        - 13.5 %
        - 19.7 %
      * - planar, half a period
        - 200 x 200
        - 0.61 %
        - 0.71 %
        - 7.3 %
        - 10.8 %
      * - planar, three periods
        - 100 x 100
        - 3.67 %
        - 5.56 %
        - 56.5 %
        - 76.8 %
      * - radial, half a period
        - 100 x 100
        - 1.02 %
        - 0.71 %
        - 4.2 %
        - 2.8 %
      * - radial, half a period
        - 200 x 200
        - 0.29 %
        - 0.20 %
        - 2.2 %
        - 1.8 %
      * - radial, three periods
        - 100 x 100
        - 3.80 %
        - 1.54 %
        - 20.7 %
        - 4.5 %

   With a moving shoreline the error falls by 2.8 (planar) and 3.5 (radial)
   when the cell size is halved, rates of 1.5 and 1.8. The bounds are about
   1.3 times the errors; those of the fine runs are below half the errors of
   the coarse ones.

   Not repeated at second order: the closed forms of one first-order step
   (``dam_break_dry_first_step_*``, ``draining_column_first_step``), which
   do not hold for two stages, ``dam_break_dry_x_no_tolerance``, whose
   lower bound on the trace of water ahead of the front is that of the
   first-order scheme (at second order the trace is 10\ :sup:`-44` m), and
   ``draining_column_round_off`` and the draining columns at ``flume.cfl``
   0.8 and 1, which run above the 0.5 that ``flume.order = 2`` allows.

``friction_channel*``, ``friction_thin_film_second_order``
   ``Tests/Friction/inputs``: water 1 m deep starts at rest on 8 by 4 cells
   of 1 m in 8 boxes, periodic both ways, on a bed that falls by 1 mm per
   metre along x, with ``flume.bed.continue_slope = 1`` and
   :math:`n` = 0.03 s/m\ :sup:`1/3`. The fluxes are the same on every face,
   so the depth stays 1 m and the discharge obeys
   :math:`dq/dt = g h S - g n^2 q^2 / h^{7/3}`. From
   :math:`q(0) = q_0` this gives

   .. math::

      q(t) = q_n \tanh\left(\frac{t}{T_r} + \operatorname{artanh}\frac{q_0}{q_n}\right)
      \ (q_0 < q_n), \qquad
      q(t) = q_n \coth\left(\frac{t}{T_r} + \operatorname{arcoth}\frac{q_0}{q_n}\right)
      \ (q_0 > q_n),

   with Manning's normal discharge :math:`q_n = h^{5/3} S^{1/2} / n`
   (1.05409 m\ :sup:`2`/s) and :math:`T_r = q_n / (g h S)` (107.451 s).
   ``Tests/check_friction.py`` checks the depth in every cell within 64
   epsilons of the depth plus the bed, and the discharge vector against
   :math:`q(t)` down the slope. Largest relative error of the discharge in
   a cell, at 100 s unless noted:

   .. list-table::
      :header-rows: 1

      * - Test
        - Run
        - Error
        - Bound
      * - ``friction_channel``
        - order 1
        - 5.64e-4
        - 6e-4
      * - ``friction_channel_steady``
        - order 1, 1500 s, 2 ranks
        - 2.500313e-4
        - 2.5002e-4 to 2.5004e-4 (64 epsilons wider)
      * - ``friction_channel_second_order``
        - order 2
        - 1.66e-4
        - 2e-4
      * - ``friction_channel_second_order_half_step``
        - order 2, ``flume.cfl = 0.25``, 2 ranks
        - 8.31e-5
        - 1e-4
      * - ``friction_channel_second_order_steady``
        - order 2, 1500 s, 3 ranks
        - 1.2e-13
        - 1024 epsilons
      * - ``friction_channel_second_order_diagonal``
        - order 2, slope (-0.0006, 0.0008), 3 ranks
        - 1.61e-4
        - 2e-4
      * - ``friction_channel_second_order_y``
        - order 2, slope along y, 2 ranks
        - 1.66e-4
        - 2e-4
      * - ``friction_channel_second_order_slowing``
        - order 2, starts at 2 m/s
        - 1.46e-4
        - 2e-4
      * - ``friction_channel_second_order_walls``
        - order 2, walls along y, SNaN fills and trap, 2 ranks
        - 1.66e-4
        - 2e-4
      * - ``friction_channel_second_order_gravity``
        - order 2, ``flume.gravity = 5`` (:math:`T_r` = 210.819 s)
        - 9.80e-5
        - 1.2e-4
      * - ``friction_thin_film_second_order``
        - order 2, 1 mm on a slope of 0.1, :math:`n` = 0.05, 10 s
        - 1.1e-9
        - 1e-8 + 4096 epsilons

   The steady error at first order is the bias of the first-order bed-slope
   term, :math:`1 - \sqrt{1 - S \Delta x / (2 h)}` = 2.50031e-4
   (:doc:`Model`), against 2.50031e-4 measured; the test pins it from both
   sides. The transient bounds at second order add 1024 epsilons, for
   rounding in single precision. At second order the transient error
   halves with the step, as friction is first order in time. In the thin
   film :math:`T_r` = 0.0645 s and the step is about 2 s (2.5 s for the
   first step); the run takes 5 steps (checked), and an explicit friction
   term would reverse the discharge (:math:`\Delta t\, g n^2 |u| / h^{4/3}`
   = 31). The film is not run at first order, where cells of 1 m drive it
   with 0.5 % of :math:`g h S` (:doc:`Model`). The walls test runs with
   ``amrex.init_snan=1 amrex.fpe_trap_invalid=1``: it stops on a corner
   ghost cell beyond a periodic side and a wall if the continued slope is
   added there.

``friction_film_uneven_second_order``
   As the thin film, but with the surface falling by 0.0999 per metre, so
   the film is 1.05 to 1.75 mm deep, on 2 ranks; 6 steps to 10 s.
   ``Tests/check_film_local_manning.py`` checks that every depth stays
   above 0, the mass is kept to 64 epsilons, and the discharge in every
   cell is the local Manning discharge :math:`h^{5/3} S^{1/2} / n` of the
   surface slope :math:`S` (centred differences) within 3e-3; it is
   2.36e-3. Friction balances the surface slope when :math:`T_r` is much
   shorter than the step; the rest is the inertia of the film and the
   discretisation. The test guards the second stage, which starts from an
   undamped discharge (:doc:`Model`).

``friction_lake_at_rest_islands_order1``, ``_order2``, ``friction_basin_mass_second_order``
   Guards, not checks of the friction formula: they pass with friction off.
   ``lake_at_rest_islands`` with :math:`n` = 0.05 at both orders: friction
   does not move still water, nor put a NaN in a dry cell. The closed basin
   at second order with :math:`n` = 0.05 on 3 ranks: the mass is kept to 16
   epsilons.

``dam_break_dry_x_outflow``, ``dam_break_dry_y_outflow``, ``_second_order``
   ``dam_break_dry_x`` and ``dam_break_dry_y`` with an outflow side ahead of
   the front (``xhi``, and ``ylo`` for the dam along y), run to 9 s. The
   front leaves the domain at 5.64 s; the rarefaction reaches the wall behind
   it at 11.3 s. The outflow is supercritical, so Ritter's solution of an
   unbounded domain still holds. ``check_dam_break.py --open-side`` compares
   with it and requires the sum of the depth to drop. L1 errors at first and
   second order: 0.45 % and 0.10 % of the depth, 0.30 % and 0.07 % of the
   discharge; bounds 0.6 % and 0.2 %. With a wall in place of the outflow
   the errors are 0.99 % and 0.69 %, above the bounds.

``dam_break_x_outflow``, ``dam_break_x_outflow_second_order``
   The wet-bed dam break with the dam at 70 m and an outflow side at
   ``xhi``, run to 10 s. The shock leaves at 7.17 s and the flow behind it
   is subcritical. The weak wave it sends back (:doc:`Model`) leaves the
   constant state 0.17 % low in depth and 0.28 % high in velocity next to
   the side at first order, 0.68 % and 1.97 % at second order. Bounds: L1
   error 2 % and 0.5 % (the runs give 0.67 % and 0.20 % of the depth step,
   1.21 % and 0.32 % of the discharge), constant state 0.5 % and 2.5 %. With
   a wall in place of the outflow the reflected shock gives L1 errors of 6 %
   (depth) and 11 % (discharge) at first order, 5.5 % and 10.1 % at second
   order. The shock is not checked once it has left.

``lake_at_rest_segments``, ``lake_at_rest_segments_second_order``, ``lake_at_rest_islands_segments_second_order``
   The lake at rest over the rough bed with outflow segments on parts of
   ``xlo`` and ``yhi`` and on all of ``xhi``, and an explicit wall segment
   on part of ``ylo``, on 3 (and 2) ranks, every value a signalling NaN at
   the start. Still water stays still to 16 epsilons of the surface, with
   islands too.

``outflow_whole_side``, ``outflow_split_side``
   A dam break over a dry bed on 400 by 12 cells in 100 boxes on 3 ranks, at
   second order, run to 7 s with an outflow side at ``xhi``: once as one
   segment, once as two outflow segments that meet inside a box. The
   plotfiles of the two runs are identical to the last bit
   (``Tests/compare_plotfiles.py``).

``partial_outflow``, ``partial_outflow_second_order``
   The dry-bed dam break on 200 by 40 cells (100 m by 20 m) in boxes of up
   to 16 cells, on 2 and 3 ranks, run to 9 s. The side at ``xhi`` is an
   outflow from y = 10 m to 20 m and a wall below: at first order an
   explicit wall segment listed after the outflow segment, at second order
   the gap no segment covers. ``Tests/check_partial_outflow.py`` checks the
   rows at least 5 m from y = 10 m. In the open rows the L1 error against
   Ritter's solution is at most 1 % (first order; the run gives 0.75 % of
   the depth and 0.49 % of the discharge) and 0.3 % (second order; 0.20 %
   and 0.13 %). In the wall rows the water arriving stops in a bore; the
   depth at the wall is within 6 % of the bore depth from the jump
   conditions with the arriving Ritter state (the runs give 3.9 % and
   4.3 %). The sum of the depth must drop.

``outflow_boxes_one_box``, ``outflow_boxes_one_cell``
   The basin dam break (rough bed sloping in x, water moving at an angle)
   on 20 by 12 cells, periodic in y, at second order, with an outflow
   segment on part of ``xhi``, every value a signalling NaN at the start,
   run to 8 s: in one box on one rank, where the sum of the depth must drop
   by the outflow, and in boxes of one cell on 2 ranks. Boxes of one cell
   are thinner than the two ghost layers, and some end one cell short of the
   side. Every value of every cell is the same in the two runs
   (``Tests/compare_plotfiles.py --cells``).

``boundary_segment_cells``
   The printed cells of two segments whose ends fall on cell centres: a
   centre on the first coordinate is in the segment, a centre on the second
   is not.

``raster_bed_plane``, ``raster_bed_coarse``
   ``Tests/Raster/inputs_bed``: 48 by 32 cells of 1 m in boxes of 8, no
   step. ``raster_bed_plane`` reads a plane
   :math:`z = 2 + 0.03 x - 0.02 y + 0.001 x y` given at the centres of 5 m
   cells (``xllcenter``), with columns beyond the domain and a
   ``NODATA_value`` cell that no cell reads; bilinear interpolation
   reproduces the plane (observed largest error 8.9e-16 m).
   ``raster_bed_coarse`` reads 6 by 8 cells of 8 by 4 m (``dx`` and ``dy``)
   that hold the domain exactly, with values that hold no plane; the cells
   within half a raster cell of the sides take the edge values (observed
   largest error 0). ``raster_bed_relative_path`` reads the same file by a
   name relative to the directory the run starts in.
   ``Tests/Raster/check_raster.py`` computes the bed from the file, in
   double precision, and allows 64 epsilons of the largest bed plus 1 m.
   The numbers in brackets here are what the runs give, not bounds.

``raster_streets``, ``_second_order``, ``_mask``, ``_small_boxes``
   ``Tests/Raster/inputs_streets``: a dam break on 40 by 40 cells of 1 m
   with walls, water 1 m deep west of x = 8 m, run to 8 s. A bed falling to
   the east (0.075 to 0.025 m between raster centres) from a raster of 20 m
   cells, and blocks of about 4 by 4 m with 4 m streets east of x = 12 m,
   3, 5 or 8 m high, from a raster of 2 m cells whose corner is at
   (-0.75, -0.75): every grid centre lies a quarter or three quarters of a
   raster cell from a raster edge, so a shift of half a raster cell in the
   sampling moves the footprints. 256 building cells (checked in the
   output). In every plotfile the bed is the raster plus the heights, every
   building cell holds no water (depth exactly 0), and the sum of the depth
   is the same within 64 epsilons; at 8 s a street cell beyond x = 14 m
   holds more than 1 cm (0.381 m at order 1). ``_mask`` gives the same
   footprints as a 0/1 mask with ``flume.buildings.height = 5``.
   ``_small_boxes`` runs in 100 boxes of 4 by 4 cells and gives the same
   cells as ``raster_streets`` in 9 boxes (``compare_plotfiles.py --cells``).

``raster_lake_at_rest_order_1``, ``_order_2``, ``raster_streets_overtop``
   ``Tests/Raster/inputs_lake``: the streets with the surface at 4 m and no
   dam, run to 4 s. The blocks of 3 m are under 1 m of water and those of 5
   and 8 m stand dry (160 cells). In every plotfile every cell holds the
   depth :math:`\max(4 - z_b, 0)` and no discharge, within 64 epsilons
   (observed: depth error 8.9e-16 m, discharge 1.6e-14 m\ :sup:`2`/s).
   ``raster_streets_overtop`` runs the dam break over the mask raised by
   0.3 m: water goes over the buildings (0.189 m deep at most on a building
   cell) and the sum of the depth does not change.

``raster_manning_rows``, ``_second_order``, ``raster_manning_uniform_*``
   The friction channel (``friction_channel*``) with the Manning coefficient
   from a map, one value per row of 1 m: 0.02, 0.03, 0.045 and 0.06 from the
   south. The water in a row does not exchange momentum with its neighbours
   (no velocity across the rows), so each row follows :math:`q(t)` with its
   own :math:`n` (``check_friction.py --manning-rows``). Largest relative
   error 5.63e-4 at order 1 (bound 6e-4) and 1.91e-4 at order 2 (bound
   2e-4), as for :math:`n` = 0.03 alone. A uniform map of 0.03 gives the
   plotfiles of ``flume.friction.manning = 0.03`` to the last bit.
   ``raster_manning_ties`` reads a map of 0.1 m rows starting at y = -2 m:
   every grid centre lies on the edge between two rows, and its position
   computed from the part of the raster kept is a whole number less a few
   epsilons. The row above holds the coefficient of the grid row and the
   row below 0.09; each grid row follows :math:`q(t)` with its own
   :math:`n`, so every centre took the row above.

``abort_raster_*``
   Each check of the reader, on files in ``Tests/Raster/bad``: an unknown
   or repeated header keyword, no ``ncols`` or ``nrows``, a fraction for
   ``ncols``, a header value that is not a number, both a corner and a
   centre in x and in y, both ``cellsize`` and ``dx``, ``dx`` without
   ``dy``, a cell size of 0, an extent that is not finite, a header that
   asks for more values than the file can hold, a value that is not a
   number, too few and too many values, a raster that does not hold the
   domain in x and in y, a cell that reads ``NODATA_value`` (the count: 36
   cells for the bed, which reads it with a weight, and 16 for the Manning
   map and the buildings), a file that cannot be opened, a directory, and an
   empty name for each of the three keys. Then the keys ``flume.bed.file``
   cannot be used with, ``flume.friction.manning`` with
   ``flume.friction.manning_file``, a negative Manning coefficient or
   building height, and ``flume.buildings.height`` of 0 or without
   ``flume.buildings.file``. In a single-precision build,
   ``abort_raster_beyond_float`` reads a value of 1e39.

``rain_flat_constant_np1``, ``_np2``, ``rain_flat_hyetograph_*``, ``rain_flat_no_stop_time``, ``rain_dry_*``
   Rain on a closed, flat, dry domain (``Tests/Rain/inputs_flat``: 24 by 16
   cells of 2 m in 6 boxes, walls on all sides). The rain is the same in
   every cell, so the water stays level and still: in every plotfile every
   cell must hold the rain fallen by the plotfile's time, within 4 epsilons
   of that depth per step, and both discharges must be exactly 0
   (``Tests/Rain/check_rain.py --uniform``, which computes the rain from
   the rate or the file itself). Observed: equal to the last bit in double
   precision.

   * ``rain_flat_constant_np1``, ``_np2``: 36 mm/h to 600 s, 6 mm, on 1 and
     2 ranks. The domain starts dry: the first step computes no fluxes (the
     output says so) and the rain still falls.
   * ``rain_flat_hyetograph_step``: ``hyetograph_steps.txt``, 36 mm/h for
     600 s, 72 mm/h for 600 s, then none, held between entries: 18 mm at
     1800 s. ``_linear``: ``hyetograph_triangle.txt``, from 0 to 72 mm/h and
     back in 1200 s, linear: 12 mm at 1800 s. ``_step_file_linear`` reads the
     step file with linear interpolation. Every plotfile is checked at its
     own time, on the rising and falling sides and after the rain.
   * ``rain_flat_hyetograph_late_start``: ``hyetograph_late.txt``, whose
     first time, 600 s, is after the start, with 36 mm/h: that rate falls
     from :math:`t = 0` (12 mm at 1800 s), and the run prints a note.
   * ``rain_flat_deep_water``: 0.36 mm/h on 1 m of still water
     (``inputs_flat_deep``), on 2 ranks. The rain of a step, about
     1.6e-8 m, is below half the spacing of single-precision numbers near
     1 m; the compensated sum keeps it. Every cell holds 1 m plus the rain
     within 4 epsilons of the depth, not scaled by the number of steps.
   * ``rain_flat_no_stop_time``: the dry domain of ``Tests/Skeleton/inputs``
     with three steps and no stop time: the rain sets the step.
   * ``rain_dry_burst_no_stop_time``: no stop time and a rate of 0 at the
     start and at the end, with 1 mm/s from 100 s to 101 s
     (``hyetograph_burst.txt``). Before 100 s only the burst, an entry of
     the table inside the interval the step looks ahead over, gives a
     finite step; the domain gets 1 mm.
   * ``rain_dry_start_discharge``: the Skeleton domain starting with
     5e-7 m of water (below ``flume.dry_depth``) moving at (1, -2) m/s,
     under 36 mm/h for 5 steps. The domain is dry at first, so no fluxes
     are computed; after every step both discharges are exactly 0 and the
     depth is 5e-7 m plus the rain (the initial state, plotfile 0, is not
     checked).
   * ``rain_dry_tiny_rate``: 1e-300 mm/h on the dry Skeleton domain with no
     stop time, for 2 steps. The Newton iteration of the rain-limited step
     converges from the smaller one-term root; from the larger it would
     need 171 iterations, more than the 100 allowed, and the run would
     stop. In single precision the depth, about 4e-305 m, rounds to 0, and
     the check expects 0.

``rain_basin_volume``
   The closed basin of ``basin_mass_*`` (rough sloping bed, moving water,
   waves off the walls) under the triangular hyetograph, on 2 ranks: the sum
   of the depth over the cells grows by the number of cells times the rain
   fallen, within 4 epsilons of the sum per step (observed: 8.9e-13 of an
   added 6.4e-3). In single precision that tolerance is larger than the
   rain on this deep water; there ``rain_flat_deep_water`` checks the rain
   added to deep water.

``rain_plane_kinematic_wave``, ``_order2``
   Rain on an impermeable plane (``Tests/Rain/inputs_plane``): 100 m long,
   falling at :math:`S_0 = 0.01` towards an outflow side, a wall upstream,
   Manning :math:`n = 0.01` s/m\ :sup:`1/3`, rain :math:`i = 50` mm/h on a dry
   bed from :math:`t = 0`, cells of 0.1 m along the slope, on 2 ranks. The
   kinematic-wave approximation keeps only the bed slope and friction in the
   momentum equation, :math:`q = \alpha h^{5/3}` with
   :math:`\alpha = S_0^{1/2} / n`, and

   .. math::

      \frac{\partial h}{\partial t} + \frac{\partial q}{\partial x} = i .

   Along the characteristics :math:`dx/dt = \tfrac53 \alpha h^{2/3}` the depth
   grows as :math:`dh/dt = i`. A characteristic that starts at the wall
   (:math:`x = 0`, :math:`h = 0`) reaches :math:`x` at
   :math:`t_x = (x / (\alpha\, i^{2/3}))^{3/5}`; before that the depth at
   :math:`x` is that of a plane with no upstream inflow, and after it the
   steady depth:

   .. math::

      h(x, t) = \begin{cases} i\, t, & t < t_x \\
                (i\, x / \alpha)^{3/5}, & t \ge t_x \end{cases}
      \qquad
      q(x, t) = \alpha\, h(x, t)^{5/3} .

   Here :math:`t_{100} = 349` s; the outlet depth at equilibrium is 4.8 mm,
   with a Froude number of 1.3, so the flow leaves supercritical and the
   outflow side sets nothing the flow inside depends on. The kinematic
   wave neglects the pressure term, whose ratio to the bed slope,
   :math:`\partial h / \partial x / S_0 \approx 0.6\, h / (x S_0)`, is
   0.003 at the outlet. In plotfiles to 900 s the check computes the L1
   error of the depth over the cells, relative to the L1 norm of the exact
   depth, and the error of the outlet discharge (mean of the last column)
   relative to :math:`i L`:

   .. list-table::
      :header-rows: 1

      * - Order, cell, cfl
        - depth L1 (largest)
        - outlet (largest)
      * - 1, 0.2 m, 0.5
        - 10.2 %
        - 11.1 %
      * - **1, 0.1 m, 0.5 (test)**
        - **5.5 %**
        - **5.7 %**
      * - 2, 0.2 m, 0.5
        - 2.8 %
        - 3.6 %
      * - 2, 0.1 m, 0.5
        - 1.2 %
        - 2.6 %
      * - **2, 0.1 m, 0.125 (test)**
        - **0.22 %**
        - **1.3 %**

   At order 1 the errors fall with the cell size at about first order; the
   largest depth error is at the upstream end, where the film is thinner
   than the bed step of a cell (Delestre et al., 2012), and the outlet error
   is largest on the rising limb. At order 2 and cfl 0.5 the error is mostly
   that of friction split from the stages (below), which halves with the
   step, not with the cell: so the order-2 test runs at cfl 0.125, where a
   rain rate off by 1 % or a Manning coefficient off by 2 % raises the depth
   L1 error above 1 %. Bounds: 6 % and 6.5 % at order 1, 0.5 % and 1.6 % at
   order 2. At order 1 the outlet discharge at equilibrium is
   :math:`i L` to round-off (observed 1.2e-13; bound 1e-9, 1e-4 in single
   precision): the outlet cell's discharge is the flux through the side,
   and at equilibrium the side passes all the rain on the plane. At order 2
   the cell value is a cell average, 1.2 % below the flux; the steady
   discharge in the cells is about 3 % below :math:`i x` at 0.2 m and
   1.3 % at 0.1 m. That matches friction acting once per step on the state
   after the two stages, with a step of about 4 % of the film's friction
   time :math:`h^{4/3} / (g n^2 |u|)`: at 0.1 m it is 2.75 %, 1.38 %,
   0.69 % and 0.35 % at a quarter of the plane for cfl 0.5, 0.25, 0.125 and
   0.0625, and the same at equal steps on 0.2 m cells. At order 1 the
   implicit friction balances the bed slope exactly at steady state.

   On a plane at 1 cm per metre with 1 m cells and :math:`n = 0.03`, the bed
   falls by more than the film depth from one cell to the next over most of
   the plane; the hydrostatic reconstruction
   then slows the film (depth error 27 % at equilibrium; Delestre et al.,
   2012), and the outflow is subcritical, which a zero-gradient side does
   not hold (water ponds in the last cells). The test avoids both.

``abort_rain_*``
   Each check of the rain inputs: the three reserved and an unknown
   ``flume.rain.type``; a negative or infinite ``flume.rain.rate``; a missing
   rate or file for its type; a rate, file or interpolation given with the
   other type or with none; an unknown interpolation; a misspelled key, which is reported as unknown
   before the missing rate; an empty file name; a
   file that cannot be opened (with an absolute and with a relative path,
   whose message names the directory it is taken from) and a directory; and on files in
   ``Tests/Rain/bad`` a line of one or three items, a time that is not a
   number or beyond double precision, an infinite rate, a negative rate, a
   time below or equal to the one before, and a file with only comments and
   blank lines. ``abort_rain_zero_dry_without_stop_time``: rain of 0 mm/h on
   a dry domain with a step limit and no stop time stops with the message
   that asks for ``flume.stop_time``.

``abort_*``
   At least one test per check in :doc:`Inputs`: missing keys, wrong numbers
   of values, values outside the range (including integers that would wrap),
   corners that are not finite, not ordered, too large, or too close for the
   cell size, an unknown key in the ``flume`` namespace and an unused key in
   another namespace; the ranges of the time-loop, boundary, bed, dam and
   initial-state keys, the keys a dam needs, and a misspelled dam key, which
   is reported as unknown before any key that needs it; a dry domain with a
   step limit and no stop time; finite inputs that give a bed, an initial
   state or a wave speed that is not finite; the ranges of
   ``flume.dry_depth``, of the bowl and of the tilted and curved surface, the
   keys these need and the keys they cannot be used with; water that is all
   below ``flume.dry_depth`` with a step limit and no stop time; a
   time step that is zero or not finite; and the range of ``flume.order``,
   the keys the wave needs and the keys it cannot be used with; the range of
   ``flume.friction.manning`` and of ``flume.bed.continue_slope``, and
   ``flume.bed.continue_slope = 1`` without a periodic side across which
   the bed slopes; each check of ``flume.bc.*`` (names, sides, types, the
   two reserved types, a periodic side, ranges, a range with no cell centre,
   overlapping segments, the keys of a segment that is not listed, a
   misspelled key of a segment). The
   wrong-type values that AMReX's parser
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
     - GCC and Clang with MPI, and GCC Debug with AMReX's
       testing mode: build with warnings as errors, then every test. No MPI
       with single precision: build and the 1-rank tests. CUDA, HIP and SYCL:
       compile only, because the hosted runners have no GPU.
   * - ``style.yml``
     - ``codespell``, no tabs, no trailing whitespace, and ``actionlint`` on
       the workflow files.
   * - ``docs.yml``
     - This documentation, with warnings as errors.
