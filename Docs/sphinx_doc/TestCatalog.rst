Test catalogue
==============

Tests named ``*_np1``, ``*_np2``, ``*_np4`` run on 1, 2 and 4 MPI ranks.
Tests named ``*_second_order`` repeat a test with ``flume.order = 2``.
"eps" is the machine epsilon of the build.

Grid, output and time loop
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``skeleton_np*``, ``skeleton_defaults``
     - Plotfile headers and every cell against a closed form: names, step,
       time, domain, cell size, box coverage, one data file per rank. The
       defaults of the grid and output keys.
   * - ``skeleton_parmparse_*``
     - A run given ``amrex.parmparse.verbose`` (``= 1``, ``= 0``) or
       ``amrex.parmparse.v`` ends with exit code 0 and reports no unused
       key.
   * - ``uniform_np*``, ``uniform_channel_x``, ``_y``
     - A uniform flow in a periodic domain, or along a channel between
       walls, is unchanged in every plotfile; step count and times as
       computed by hand (step 0.0625 s, 16 steps).
   * - ``uniform_defaults``, ``uniform_gravity_default``,
       ``uniform_last_step_cut``, ``uniform_stop_time_roundoff``,
       ``uniform_max_step_first``, ``uniform_negative_u``, ``_v``,
       ``uniform_unequal_cells``, ``uniform_stop_time_zero``,
       ``uniform_max_step_zero``
     - Defaults of ``flume.cfl``, ``flume.plot_int`` and ``flume.gravity``;
       the last step shortened or lengthened onto the stop time; a step
       limit before the stop time; negative velocities; cells of 2 m by
       1 m; a stop time or step limit of zero.
   * - ``dry_stop_time``, ``dry_stop_time_huge``
     - A dry domain takes one step to the stop time, also for 10\ :sup:`308`
       s.
   * - ``max_dt_plane``
     - Rain of 50 mm/h on a sloping plane, 64 x 8 cells of 1 m, 600 steps
       of ``flume.max_dt = 0.25`` s on 2 ranks (``Tests/ERF``). Each
       plotfile is at step N and time 0.25 N exactly, so every step was
       ``flume.max_dt``. The stand-alone run of ``erf_flume_plane``.
   * - ``max_dt_dry_no_stop_time``
     - A dry domain with no rain and no stop time takes steps of
       ``flume.max_dt``: 3 steps of 2.5 s.
   * - ``abort_max_dt_zero``, ``abort_max_dt_negative``
     - ``flume.max_dt`` of 0 or below stops the run.

Initial state and water at rest
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``sloping_bed_initial``, ``dam_initial_on_centre``, ``_y``,
       ``bowl_mound_initial``
     - Every cell of the initial state against its closed form; a cell
       centre on the dam belongs to the high side.
   * - ``sloping_bed_np*``
     - A lake at rest that does not cover its sloping bed, in 8 uneven
       boxes: every cell after three steps; the step is reduced over ranks
       that hold only dry boxes.
   * - ``lake_at_rest_np*``, ``lake_at_rest_islands``,
       ``lake_at_rest_islands_dry_depth``
     - Level surface over a rough bed, 50 steps: unchanged within 16 eps of
       the surface elevation; with dry islands; with a dry depth of 0.05 m.
   * - ``lake_at_rest_dry_depth_at_a_face``
     - The dry depth set to the depth one side of a face presents: the cut
       is decided once per face and the lake stays at rest (64 eps).
   * - ``lake_at_rest_segments``, ``lake_at_rest_islands_segments``
     - The lake with outflow and wall segments, every value a signalling
       NaN at the start.

Wetting and drying
------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``draining_column_cfl_0.8``, ``_cfl_1.0``
     - One column of water between dry cells, to 5 s: no negative depth,
       mass within 16 eps, no discharge in dry cells.
   * - ``draining_column_first_step``
     - The closed form of the first step at CFL 1: fluxes scaled by
       :math:`3 / (4 n)`, the column ends empty (64 eps).
   * - ``draining_column_round_off``
     - A depth the outflow limit leaves a rounding error below zero is set
       to zero.
   * - ``film_passes_no_water``, ``film_behind_a_bed_step``
     - A film below the dry depth passes no water on; the tolerance acts on
       the depth at the face, not on the depth of the cell.
   * - ``thin_film_time_step``
     - A fast film below the dry depth does not set the time step (2 steps;
       1130 with a dry depth of 0).

Dam break and shear
-------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``dam_break_x_np*``, ``dam_break_y``
     - Wet-bed exact solution at 5 s: L1 errors, constant state, shock
       position, mass (:doc:`Verification`).
   * - ``dam_break_dry_x``, ``_y``, ``dam_break_dry_x_no_tolerance``
     - Ritter's solution: L1 errors, no negative depth, no water beyond the
       exact front; with a dry depth of 0, a trace below 10\ :sup:`-11` of
       the depth.
   * - ``dam_break_dry_first_step_x``, ``_y``
     - The closed form of the first step: pins the dry-bed wave speeds.
   * - ``dam_break_x_second_order_one_cell_wide``
     - Walls one cell apart, narrower than the two ghost layers.
   * - ``shear_x_np*``, ``shear_x_reversed``, ``shear_y``
     - A jump of the tangential velocity carried by a uniform flow: bounded,
       conserved, within one cell of its place; at most 10 cells wide at
       order 2.
   * - ``basin_mass_np*``, ``basin_mass_second_order_thin_boxes``,
       ``_one_cell_boxes``
     - A dam break over a rough bed between four walls, 20 s: mass within
       16 eps, also in boxes two cells and one cell wide.

Moving shoreline and convergence
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``bowl_planar_half``, ``_half_fine``, ``bowl_planar``,
       ``bowl_radial_half``, ``_half_fine``, ``bowl_radial``
     - Thacker's planar and radial solutions on 100 and 200 cells, half a
       period and three periods: L1 errors, no negative depth, mass.
   * - ``wave_x_first_order``, ``wave_x_second_order``, ``_fine``,
       ``_finest``, ``wave_y_second_order``, ``_fine``,
       ``wave_x_second_order_wide_cells``, ``_one_cell_wide``
     - Simple wave at 1 s on 100, 200 and 400 cells: each bound is below
       half the error of the run on cells twice as large, so a first-order
       scheme fails.

Friction
--------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``friction_channel``, ``_steady``, ``_second_order``,
       ``_second_order_half_step``, ``_steady``, ``_diagonal``, ``_y``,
       ``_slowing``, ``_walls``, ``_gravity``
     - Discharge against :math:`q(t)` in an inclined periodic channel:
       transient and steady, slope along x, y and diagonal, starting above
       the normal discharge, between walls, with another gravity.
   * - ``friction_thin_film_second_order``,
       ``friction_film_uneven_second_order``
     - A stiff film: 5 steps of 2 s against a friction time of 0.0645 s; a
       film of uneven depth follows the local Manning discharge within
       3e-3.
   * - ``friction_lake_at_rest_islands_order1``, ``_order2``,
       ``friction_basin_mass_second_order``
     - Friction does not move still water or change the mass.

Boundaries
----------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``dam_break_dry_x_outflow``, ``_y_outflow``
     - A dry-bed front leaves through an outflow side; Ritter's solution
       holds inside.
   * - ``dam_break_x_outflow``
     - A shock leaves through an outflow side; the reflected wave is
       bounded.
   * - ``outflow_whole_side``, ``outflow_split_side``
     - One segment and two that meet inside a box give identical plotfiles.
   * - ``partial_outflow``
     - A side that is half outflow, half wall: Ritter's solution in the
       open rows, the bore depth in the wall rows.
   * - ``outflow_boxes_one_box``, ``outflow_boxes_one_cell``
     - One box and boxes of one cell give the same cells.
   * - ``boundary_segment_cells``
     - A centre on the first coordinate of a range is in the segment; a
       centre on the second is not.

Inflow hydrographs
------------------

Decks in ``Tests/Inflow``, checked by ``check_inflow.py`` against closed
forms and by ``check_balance.py`` against the integral of the hydrograph.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``inflow_uniform_subcritical_order1``, ``_order2``
     - Uniform subcritical flow fed with its own discharge stays uniform for
       3419 steps on 4 ranks (4096 eps).
   * - ``inflow_rarefaction_{xlo,xhi,ylo,yhi}_order1``, ``_order2``
     - The discharge drops: the depth of the plateau next to each of the
       four sides follows from the invariant kept across the rarefaction.
   * - ``inflow_dry_fan_order1``, ``_order2``
     - Supercritical inflow onto a dry bed: the critical depth at the face
       and the rarefaction to the dry front.
   * - ``inflow_normal_depth_steep_order1``, ``_order2``
     - Steady inflow down a steep slope with friction reaches the normal
       depth.
   * - ``inflow_basin_volume_np4``, ``_swapped``, ``inflow_basin_late_start``
     - Two segments on a closed basin: each inflow column, and the stored
       water, equal the integral of the hydrograph (``linear``, ``step``,
       and a file that starts after 0 s) to round-off.
   * - ``inflow_basin_lake_at_rest``
     - Segments with no discharge keep still water still.
   * - ``inflow_dry_no_stop_time``
     - The hydrograph sets the step of a dry domain with no stop time.
   * - ``inflow_restart_reference``, ``inflow_restart_np2``,
       ``_np2_np1``
     - A restart while hydrograph ``a`` rises gives every plotfile to the
       last bit, and the log byte for byte on the same ranks.

Rasters
-------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``raster_bed_plane``, ``raster_bed_coarse``,
       ``raster_bed_relative_path``
     - Bilinear interpolation reproduces a plane; edge values are held
       within half a raster cell of the sides.
   * - ``raster_streets``, ``_mask``, ``_small_boxes``,
       ``raster_streets_overtop``
     - A dam break through streets: 256 building cells stay dry, mass
       within 64 eps; heights and a mask give the same footprints; water
       goes over low buildings.
   * - ``raster_lake_at_rest_order_1``, ``_order_2``
     - Water at rest among and over buildings.
   * - ``raster_manning_rows``, ``raster_manning_uniform_*``,
       ``raster_manning_ties``
     - Each row follows :math:`q(t)` with its own coefficient; a uniform
       map equals the scalar key to the last bit; a centre on a raster edge
       takes the cell above.

Rain
----

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``rain_flat_constant_np*``, ``rain_flat_hyetograph_step``,
       ``_linear``, ``_step_file_linear``, ``_late_start``
     - Every cell of a flat, closed domain holds the rain fallen, at every
       plotfile time (4 eps per step); discharges exactly 0.
   * - ``rain_flat_deep_water``
     - Rain below the rounding of the depth is kept by the compensated sum.
   * - ``rain_flat_no_stop_time``, ``rain_dry_burst_no_stop_time``,
       ``rain_dry_tiny_rate``, ``rain_dry_start_discharge``
     - A dry domain under rain gets a finite step, also from a burst ahead
       and from a rate of 1e-300 mm/h; discharges of a dry domain are 0.
   * - ``rain_basin_volume``
     - The volume of a closed basin with moving water grows by the rain
       fallen.
   * - ``rain_plane_kinematic_wave``, ``_order2``
     - Depth and outlet discharge against the kinematic-wave solution to
       900 s.

Rain that varies in space
-------------------------

Decks in ``Tests/RainField``, checked by ``check_rain_field.py``, which
computes the rain from the grid files or the cell's numbers independently of
the solver.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``rain_gridded_cells_linear_np1``, ``_np2``,
       ``rain_gridded_cells_fine``
     - One step on a closed, flat, dry domain: every cell holds the rain of
       its area, on grids that share no edge with the solver's (16 eps).
   * - ``rain_gridded_cells_aligned_nodata_ring``
     - Grid edges on solver edges, NODATA just outside the domain: no abort,
       every cell within 32 eps.
   * - ``rain_gridded_volume_step``, ``_volume_linear_shifted``
     - Water in the domain to 1200 s equals the integral of three grids,
       held and linear (4 eps per step).
   * - ``rain_gridded_dry_burst_no_stop_time``
     - A grid ahead of the run's time sets the step of a dry domain.
   * - ``rain_cell_at_rest_np*``, ``rain_cell_fast_step``,
       ``rain_cell_volume``
     - The moving cell against the closed form: at rest (16 eps), 29 m in
       one step (1e-11), leaving the domain over 1200 s (volume).
   * - ``rain_cell_zero_peak_long_step``
     - A cell of peak 0 takes one step of 1e7 s and adds no water.
   * - ``rain_replay_unit``
     - The writer of ``flume.rain.replay_file`` on its own
       (``Tests/RainReplay``, a driver linked to the FLUME library, 2
       ranks): four host steps with negative depths, a step of no rain and a
       step that starts 1e-8 s before the end of the last. The index to the
       byte and the grids to the bit against the rates computed by the
       check; restarts from the checkpoints after steps 2 and 3 write them
       again to the byte; a restart with an index that is not the run's, and
       a step that starts before the end of the last, stop.

Infiltration
------------

Decks in ``Tests/Infiltration``, checked by ``check_infiltration.py``, which
computes the closed forms independently of the solver (bisection for the
Green-Ampt equation).

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``infiltration_green_ampt_ponded``
     - 0.5 m of still water on silt loam for 6 h: every cell against the
       implicit solution :math:`K_s t = F - M \ln(1 + F/M)`.
   * - ``infiltration_green_ampt_ponded_moving``
     - The same water moving on a periodic domain: the velocity is kept.
   * - ``infiltration_green_ampt_rain``
     - 36 mm/h on dry sandy loam: all rain infiltrates up to
       :math:`F_s` at 1437 s, then ponding (round-off).
   * - ``infiltration_curve_number_rain``, ``_ia_ratio``
     - CN 75 (:math:`\lambda = 0.2`) and CN 65.31 (:math:`\lambda = 0.05`)
       under constant rain: the cumulative loss against the formula.
   * - ``infiltration_curve_number_many_steps``
     - CN 75 on 2 cm cells, 279217 steps: the loss within a few epsilons,
       with a bound that does not grow with the steps (compensated sums).
   * - ``infiltration_curve_number_map``
     - A map of curve numbers with flow and buildings: each cell loses the
       loss of its own rain; buildings and CN 100 lose nothing.
   * - ``infiltration_green_ampt_balance``
     - A sloped, walled domain, second order, hyetograph, conductivity map
       and buildings: surface water plus infiltrated water equals the rain,
       to round-off.
   * - ``infiltration_green_ampt_cap``, ``_cap_moving``
     - A 1 mm film on a soil that could take 28 m in a step: the cell loses
       exactly its water and holds no discharge.
   * - ``infiltration_green_ampt_burst``
     - Rain in the first step only, on cells that infiltration empties: the
       rain's negative carry does not take a depth below 0 in the next,
       dry step, and the water balances.

Storm-drain inlets
------------------

Decks in ``Tests/Drains``, checked by ``check_drains.py``, which computes the
exact solution through its inverse (the time to reach a depth),
independently of the solver.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``drains_uniform``, ``_weir``, ``_coefficients``
     - Level, still water over an inlet in each of 4 cells: every cell
       against the exact solution, from above the switch depth (orifice,
       then weir), from below it (weir only), and with other coefficients
       (16 eps of the starting depth).
   * - ``balance_drains_uniform``
     - The ``drained`` and ``stored`` columns against the closed form.
   * - ``drains_basin``, ``_order2``
     - One inlet in a closed basin of 256 cells: the drained volume within
       0.95 and 1 times that of one cell of the basin's area, at both orders.
   * - ``drains_basin_boxes_4_np3``
     - The basin on 3 ranks in 16 boxes: the plotfiles to the last bit.
   * - ``drains_restart_np4_np2``
     - The basin continued from a checkpoint on 2 ranks: the plotfiles to
       the last bit, ``drained_depth`` among them.
   * - ``drains_film``
     - A 1 mm film over an inlet that an explicit step would let take 13 mm:
       no depth below 0, and the depth drained in the first step that of
       the exact solution (16 eps of 1 mm).
   * - ``abort_restart_drains_moved``, ``_coefficient``
     - A restart with the inlet in another cell, or with another weir
       coefficient, stops with a message naming the inlets.

Water balance, flood maps and gauges
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - every regression test with a log
     - The residual of ``water_balance.txt`` within the bound of
       :doc:`theory/WaterBalance`, at every line, and equal to the residual
       recomputed from the columns.
   * - ``balance_rain_columns_np*``
     - Rain on a closed basin: the rain column and the stored water against
       the rain rate times the time and the area.
   * - ``balance_rain_gauges_np*``
     - The same basin on a raised bed, with maps and gauges: maximum depth,
       zero speeds and the arrival time in every cell; three gauges (on a
       face, on the high corner, on the low edge) read their cells, the
       rain fallen and the surface elevation at every line.
   * - ``gauge_face_rounding``
     - A gauge whose coordinates round just below two faces reads the cell
       above both.
   * - ``gauge_cells_sloped_bed``
     - Three gauges in three boxes on 2 ranks, on a sloped bed with water
       in motion, hold exactly the values of their cells in the plotfile.
   * - ``gauge_large_coordinates``
     - A gauge 500 km from the origin, inside a cell, stays in that cell in
       single precision.
   * - ``balance_green_ampt_ponded``
     - The infiltrated column of ponded Green-Ampt infiltration against the
       closed form.
   * - ``balance_outflow_side``, ``_two_segments``, ``_low_side``
     - Uniform supercritical flow: the water out of a whole high side, out
       of two segments at second order, and out of a low side, against the
       discharge times the length and the time.
   * - ``balance_periodic``
     - A periodic domain keeps its water exactly; no outflow column.
   * - ``balance_log_interval``
     - ``flume.log_int = 5``: lines at the multiples of 5 and at the last
       step.
   * - ``balance_log_on``, ``balance_log_off``
     - ``flume.log_int = 0`` writes no log, no gauge file and no summary
       line, and the plotfiles are byte for byte those of the run with a
       line every step.
   * - ``balance_every_term``, ``_log_interval``
     - Rain, infiltration and outflow at once, with a line every step and
       every 25 steps: the log without any one term fails the check on at
       least 90 % of the pairs of lines where that term changes.
   * - ``flood_maps_dam_break``, ``_second_order``
     - The flood maps and a gauge of the wet-bed dam break against the exact
       Riemann solution.

Checkpoint and restart
----------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``restart_reference_np1``, ``_np2``, ``restart_boxes_8_np1``
     - A moving storm cell over a sloped bed with a bowl, second order,
       curve-number infiltration, an outflow segment, flood maps and two
       gauges: the run on 2 ranks, and on 1 rank in boxes of 8 cells, gives
       every value of every plotfile of the 1-rank run to the last bit.
   * - ``restart_first_order_reference``, ``restart_first_order_no_log``
     - A dam break with rain and Green-Ampt infiltration, first order:
       the plotfiles of the 2-rank run with the log, and of the 1-rank run
       in boxes of 8 cells without it, agree to the last bit.
   * - ``restart_np1_np1``, ``restart_np2_np2``
     - A run continued from a checkpoint on the same ranks and boxes gives
       the plotfiles to the last bit and the log and gauge files byte for
       byte; the first part logs past the checkpoint (the logs are cut
       back, and a line cut short at their end is removed), or stops at it
       (the checkpoint of the last step).
   * - ``restart_np1_np1_off_interval``
     - The same from the checkpoint of the first part's last step, 405, off
       the log interval: its log line is removed.
   * - ``restart_np1_np2``, ``restart_np2_np1_boxes_8``,
       ``restart_np2_np2_boxes_8``
     - The same on another rank count or other boxes: plotfiles and gauge
       files to the last bit, the log within :math:`10^{-14}` of the largest
       value of each column.
   * - ``restart_first_order_np2``
     - The same for Green-Ampt at first order.
   * - ``abort_restart_*``, ``abort_check_*``
     - Each check of the checkpoint keys, and a checkpoint that does not
       match the inputs: a plotfile, another grid, domain, bed, source
       list, set of periodic sides, water-balance columns, outflow or
       discharge segment or gauge cell, maps asked of a checkpoint without them, and a log
       asked of a checkpoint without a balance. Copies of a checkpoint
       altered by ``Tests/Restart/make_bad_checkpoints.py``: another
       version, the other precision, a key given twice, a negative step, and a NaN in the state.

Input checks
------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``abort_*``
     - At least one test per check in :doc:`Inputs`: missing keys, wrong
       numbers of values, values out of range, unknown and unused keys,
       keys that need or exclude each other, fields that are not finite.
   * - ``abort_raster_*``
     - Each check of the raster reader, on files in ``Tests/Raster/bad``.
   * - ``abort_rain_*``
     - Each check of the rain inputs and of the hyetograph reader, on files
       in ``Tests/Rain/bad``.
   * - ``abort_infiltration_*``
     - Each infiltration key: ranges, missing and excluded keys, a value
       with a file, the curve number without rain, and files with values
       out of range or NODATA cells.
   * - ``abort_drains_*``
     - Each drain key: names, missing keys, ranges, a location outside the
       domain, an unknown key, a key of a name not listed, two inlets in one
       cell and an inlet in a building cell.
   * - ``abort_log_*``, ``abort_flood_maps_*``, ``abort_arrival_depth_*``,
       ``abort_gauge_*``, ``abort_gauges_*``
     - Each check of the log, map and gauge keys: ranges, a non-finite
       arrival depth or one below the dry depth, keys that need another,
       gauge names (also two that differ only in case), missing and
       unlisted locations, and points outside the domain.
   * - ``abort_inflow_*``
     - Each key of a discharge segment (a missing or empty file, two
       values, a file or an interpolation with a wall or outflow segment,
       an unknown interpolation), each
       check of the hydrograph reader on the files of ``Tests/Rain/bad``, a
       dry domain without a stop time and with no more inflow, and a
       restart whose segment covers other faces.
   * - ``abort_rain_gridded_*``, ``abort_rain_cell_*``
     - Each check of the index reader and the grids (files in
       ``Tests/RainField/bad``), each cell key, cell keys with another type,
       and a step that needs more than :math:`10^4` sub-intervals.

Not repeated at second order: the closed forms of one first-order step, the
lower bound of ``dam_break_dry_x_no_tolerance``, and the draining columns at
CFL 0.8 and 1, which are above the 0.5 that ``flume.order = 2`` allows.

ERF built with FLUME
--------------------

Registered only with ``FLUME_WITH_ERF=ON`` (``Tests/ERF``, ``Tests/ERFRain``,
``Tests/ERFReplay``, ``Tests/ERFStorm``). FLUME runs as
ERF's flood engine (:ref:`sec:InsideERF`).

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Tests
     - What is checked
   * - ``erf_flume_plane``
     - ``erf_exec`` with FLUME as its flood engine on the flood grid of
       ``max_dt_plane``
       (ERF's 32 x 4 columns refined by 2), ``erf.fixed_dt = 0.25``, 2 ranks.
       Every flood plotfile (every cell, time and step), the water balance
       and the two gauge files are those of the stand-alone run, to the last
       bit (``max_dt_plane``). Then ``erf.restart`` from ERF's checkpoint of step 300 writes
       the same files again, and cuts the logs back to that step. The
       ``flood_*`` fields of ERF's plotfile are the means of FLUME's cells
       over each column, within 4 machine epsilons.
   * - ``erf_flume_plane_substeps``
     - The same flood with ``erf.fixed_dt = 2`` and no cap on FLUME's step:
       FLUME takes its CFL steps (more than one in each ERF step), lands on
       the end of each ERF step, and the water balance closes to round-off.
   * - ``erf_flume_plane_short_steps``
     - No cap on FLUME's step and ``erf.fixed_dt = 0.25``, shorter than
       FLUME's own step (about 1 s): each ERF step is one FLUME step, and the
       plotfiles, water balance and gauges are those of ``max_dt_plane``, to
       the last bit. 2 ranks.
   * - ``erf_flume_plane_step_1.5``, ``_2.0``, ``_3.75``, ``_15``
     - No cap on FLUME's step, ERF's step 1.5, 2, 3.75 and 15 s to 150 s, 4
       ranks: each flood plotfile is at ERF's time, after more FLUME steps
       than ERF steps, and the water balance closes to round-off.
   * - ``flood_plane_dt_05``, ``flood_plane_dt_0125``, ``flood_plane_own_steps``
     - The stand-alone plane with every step 0.5 s, every step 0.125 s, and
       FLUME's own steps (99 to 150 s), 4 ranks: references of
       ``erf_flume_plane_step_lengths``.
   * - ``erf_flume_plane_step_lengths``
     - At 150 s, the depth and discharge of each ``erf_flume_plane_step_*``
       run against the time-converged flood :math:`2F_{0.125} - F_{0.25}`
       (Richardson extrapolation of ``flood_plane_dt_0125`` and
       ``max_dt_plane``). The time error is first order:
       :math:`|F_{0.5} - F_{0.25}| / |F_{0.25} - F_{0.125}|` is 2
       :math:`\pm` 0.1 (L1; measured 2.0003 for depth, 2.0088 for
       discharge). No run inside ERF is further from the time-converged
       flood than the run with FLUME's own steps, within 2 %, and each
       differs from that run by at most half of its time error (L1 and
       largest value). The 2 % margin is observed, not proven: the
       stand-alone run with every step capped at 2 s is 1.0103 times
       further than FLUME's own steps in the largest depth difference.
       Measured, depth L1 relative
       to the largest depth: own steps 7.9 :math:`\times` 10\ :sup:`-4`; ERF
       steps of 1.5, 2, 3.75 and 15 s 6.1, 6.7, 7.2 and 7.8 :math:`\times`
       10\ :sup:`-4`.
   * - ``erf_flume_plane_substeps_restart``
     - ``erf_flume_plane_step_2.0`` with a checkpoint at ERF step 37, after
       51 FLUME steps (the step in the checkpoint, and the step FLUME reports
       on restart): the plotfiles and logs are those of the
       run without it, to the last bit, and so are they after ``erf.restart``
       from ``chk00037``. 4 ranks.
   * - ``erf_flume_plane_ranks_restart``
     - ``erf_flume_plane_substeps_restart``, but ``erf.restart`` from
       ``chk00037`` on 3 ranks instead of 4 (the same 4 boxes, one rank
       with two). The flood plotfiles and gauge files are those of
       ``erf_flume_plane_step_2.0`` to the last bit. The sums of
       ``water_balance.txt`` are taken over the ranks in another order
       after the restart: the same steps and times, each value within
       10\ :sup:`-14` of its column's largest, the residual within
       10\ :sup:`-14` of the largest volume (measured 2.3 :math:`\times`
       10\ :sup:`-16`). The restart's output must report 3 MPI processes.
   * - ``erf_flume_start_without_flood``
     - ERF runs 37 steps of 2 s with ``erf.flood.enable = false`` and
       writes ``chk00037`` at 74 s, without FLUME state, then continues
       from it with the flood on to ``max_step = 75``, FLUME starting from
       a flat surface at 0.5 m. FLUME writes its summary at ERF's step 75;
       the water balance and gauges start at 74 s, at step 0; and the 112
       cells wet at the start (centres beyond x = 50 m on the bed falling
       1 cm per metre from 1 m) have arrival time 74 s, no other cell an
       earlier one. 4 ranks.
   * - ``abort_erf_flume_checkpoint_title``
     - FLUME reads the step of a restart from ERF's checkpoint ``Header``.
       ERF skips its title line, so a copy of that ``chk00037`` with
       another title (``erf_flume_retitle_checkpoint``) restarts in ERF,
       and FLUME stops with the message that names the header.
   * - ``erf_flume_plane_stop_time``
     - ERF ends on ``stop_time = 52.2`` after 174 steps of 0.3 s, which sum
       to 1.4 :math:`\times` 10\ :sup:`-13` s less, with ``flume.max_dt``
       10\ :sup:`-8` s below ERF's step. FLUME takes one step in each ERF
       step and writes its summary after the 174th.
   * - ``abort_erf_flume_host_*``
     - Each of the 12 keys that ERF sets (:ref:`sec:InsideERF`) stops
       ``erf_exec`` with the message that names ERF's key.
   * - ``erf_flume_rain_cells``
     - ``flume.rain.type = microphysics`` (``Tests/ERFRain``): Kessler rain
       from ERF on 16 x 16 columns, each over 2 x 2 flood cells, 2 ranks.
       No water moves (``flume.dry_depth = 1``) and FLUME takes two steps
       in each ERF step. After 60 and 120 ERF steps every flood cell holds
       the accumulated precipitation of its column (ERF's
       ``precip_total_accum``) divided by 1000, within
       10\ :sup:`-12` of the largest; the precipitation varies 3.3-fold
       over the columns.
   * - ``erf_flume_rain_volume``
     - The same rain on a dam break of 0.1 m over a flat periodic bed: the
       water moves. The rain of the water balance is the precipitation
       ERF accumulated over the grid, within the rain each cell may still
       hold in its rounding carry (3.7 :math:`\times` 10\ :sup:`-13`
       m\ :sup:`3` of 0.63 m\ :sup:`3`). Then ``erf.restart`` from the
       checkpoint of step 60 writes the same flood plotfile of step 120 and
       water balance, to the last bit.
   * - ``erf_flume_rain_frozen``
     - WDM6 snow and graupel from ERF over ground below freezing, 2 ranks:
       every flood cell holds its column's precipitation, liquid and frozen,
       divided by 1000, and the frozen part is at least a fifth of the total
       (all of it here).
   * - ``abort_erf_flume_rain_no_moisture``
     - ``flume.rain.type = microphysics`` with no moisture model: ERF stops
       at the start with its message.
   * - ``abort_erf_flume_rain_rate``
     - ``flume.rain.rate`` with ``microphysics`` inside ERF stops the run
       with the message that names the key.
   * - ``erf_flume_rain_replay``
     - ``flume.rain.replay_file`` (``Tests/ERFReplay``): the dam break of
       ``erf_flume_rain_volume`` on 4 ranks, ``erf.fixed_dt = 0.25``, one
       FLUME step per ERF step, 120 steps. The index holds 119 grids of
       rain, closes at 30 s, and its grids hold the rain of the water
       balance. The stand-alone run that reads it back (16 boxes on 4
       ranks) gives every depth of the coupled flood at 15 s and 30 s to
       64 eps (h + R), with R the largest rain of a cell (largest
       difference 0.008 of the bound), and the same rain to
       1.6 :math:`\times` 10\ :sup:`-16` of it. Run on to 60 s, it adds no
       rain after 30 s while the water moves. ``erf.restart`` from step 60
       writes the same index and grids, to the byte; the restart without
       the index, with the index cut back to 7.5 s, or from a checkpoint
       without what the writer holds back, stops.
   * - ``abort_erf_flume_rain_replay_directory``,
       ``abort_erf_flume_rain_replay_hash``
     - ``flume.rain.replay_file`` that does not end in a file name, or
       whose name holds a ``#``, stops the run.
   * - ``erf_flume_storm_catchment``, ``erf_flume_storm``
     - The coupled example (``Tests/ERFStorm``, :ref:`sec:CoupledExample`):
       a convective cell from ERF's Kessler microphysics over a catchment
       of 320 km\ :sup:`2`, 960 ERF steps of 5 s on 4 ranks, one FLUME step
       each. ERF's rain on the catchment, at each of the 17 plt2d
       plotfiles, is that of the hyetograph made from the rain FLUME wrote,
       to 10\ :sup:`-12` of the total. The stand-alone replay (8 h) gives
       the coupled flood at 40 and 80 min to 64 eps (h + R). The hyetograph
       run (8 h) closes its water balance (factor 256: every cell wet at
       once, all rain infiltrated from about 18 to 36 min) and rains the
       hyetograph's depth at every step. The replay's outlet peak is at
       least 3 times the hyetograph run's and at least 1 h earlier (86.0
       against 16.4 m\ :sup:`3`/s, 2.9 h, in one run). The fixture writes
       the rasters.
