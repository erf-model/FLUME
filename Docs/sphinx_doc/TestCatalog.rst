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
   * - ``abort_log_*``, ``abort_flood_maps_*``, ``abort_arrival_depth_*``,
       ``abort_gauge_*``, ``abort_gauges_*``
     - Each check of the log, map and gauge keys: ranges, a non-finite
       arrival depth or one below the dry depth, keys that need another,
       gauge names (also two that differ only in case), missing and
       unlisted locations, and points outside the domain.
   * - ``abort_rain_gridded_*``, ``abort_rain_cell_*``
     - Each check of the index reader and the grids (files in
       ``Tests/RainField/bad``), each cell key, cell keys with another type,
       and a step that needs more than :math:`10^4` sub-intervals.

Not repeated at second order: the closed forms of one first-order step, the
lower bound of ``dam_break_dry_x_no_tolerance``, and the draining columns at
CFL 0.8 and 1, which are above the 0.5 that ``flume.order = 2`` allows.
