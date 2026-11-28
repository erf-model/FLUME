# ---------------------------------------------------------------------------
# Inflows inside the domain: points and lines (flume.inflows.*;
# Tests/InflowSource/README.md). Included from Tests/CMakeLists.txt, which
# defines the functions and ${rok}, after Tests/Restart/RestartTests.cmake,
# which defines flume_add_restart_test, ${check_rs} and rs_abort.
#
# Tests/InflowSource/check_inflow_source.py checks the first step from a dry
# domain against an exact clipping of the line and the closed-form step;
# Tests/check_balance.py pins the inflow columns of the water balance against
# the integrals of the hydrographs; Tests/Inflow/check_inflow.py compares a
# channel with the normal depth.
# ---------------------------------------------------------------------------
set(is_dir ${CMAKE_CURRENT_SOURCE_DIR}/InflowSource)
set(is_basin ${is_dir}/inputs_basin)
set(is_files flume.inflows.a.file=${is_dir}/hydrograph_triangle.txt
             flume.inflows.b.file=${is_dir}/hydrograph_steps.txt)
set(check_is ${Python3_EXECUTABLE} ${is_dir}/check_inflow_source.py)
set(check_is_bal ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py)
set(check_is_chan ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Inflow/check_inflow.py)
set(is_dry flume.inflows.names=l flume.inflows.l.file=${is_dir}/hydrograph_half.txt)

# --- The first step: the parts of a line, and the time step ------------------
# One step of 0.5 m3/s from a dry, closed domain of 2 x 1 m cells
# (inputs_dry) along a line across 16 cells and 4 of the 6 boxes, on 4
# ranks: the step is the root of the rain-limited cubic with the largest
# part of a cell (0.122) times Q over the cell area as the rain rate (0.716 s;
# Flood::compute_dt), and each cell holds its part of the line times Q dt
# over its area, all within 1e-12 (the run: 1.6e-16 of the largest depth).
flume_add_regtest(inflow_source_line_first_step NPROCS 4
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} "flume.inflows.l.line=1.3 0.7 14.2 8.9"
    OUTPUT "inflow l: hydrograph from [^\n]*hydrograph_half\\.txt, 1 entries, linear in time, up to 0\\.5 m3/s into 16 cell\\(s\\) along the line from \\(1\\.3, 0\\.7\\) to \\(14\\.2, 8\\.9\\) m"
    CHECK ${check_is} dry --q 0.5 --line 1.3 0.7 14.2 8.9)
# The line the other way round, from (14.2, 8.9) to (1.3, 0.7), on 1 rank:
# the same parts.
flume_add_regtest(inflow_source_line_reversed NPROCS 1
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} "flume.inflows.l.line=14.2 8.9 1.3 0.7"
    CHECK ${check_is} dry --q 0.5 --line 14.2 8.9 1.3 0.7)
# A line along x from the low to the high edge of the domain, y = 3.5 m,
# through the middle of row 3: eight cells, each 1/8.
flume_add_regtest(inflow_source_line_edge_to_edge NPROCS 2
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} "flume.inflows.l.line=0 3.5 16 3.5"
    OUTPUT "into 8 cell\\(s\\)"
    CHECK ${check_is} dry --q 0.5 --line 0 3.5 16 3.5)
# A point at (9.1, 6.4) m: cell (4, 6) takes it all.
flume_add_regtest(inflow_source_point_first_step NPROCS 2
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} "flume.inflows.l.location=9.1 6.4"
    OUTPUT "up to 0\\.5 m3/s into 1 cell\\(s\\), \\(4, 6\\), at \\(9\\.1, 6\\.4\\) m"
    CHECK ${check_is} dry --q 0.5 --point 9.1 6.4)
# A line along the face x = 0.6 m of cells of 0.3 m (10 x 10 cells on 3 m,
# where 0.6 / 3 * 10 rounds to 1.9999999999999998): the cells to the right
# of the face, (2, 1) and (2, 2), take it, as the point (0.6, 1.5) takes
# cell (2, 5).
flume_add_regtest(inflow_source_line_on_face NPROCS 1
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} flume.n_cell=10\ 10 flume.prob_hi=3\ 3 "flume.inflows.l.line=0.6 0.45 0.6 0.75"
    OUTPUT "into 2 cell\\(s\\)"
    CHECK ${check_is} dry --q 0.5 --line 0.6 0.45 0.6 0.75 --cells 2)
# A diagonal of 30 x 45 cells of 0.1 m, from (3, 0) to (0, 4.5), through a
# corner of the grid at every third cell: 60 cells, none of those it only
# touches at a corner (without dropping the slivers of rounding there, 70).
flume_add_regtest(inflow_source_line_through_corners NPROCS 2
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} flume.n_cell=30\ 45 flume.prob_hi=3.0\ 4.5 "flume.inflows.l.line=3.0 0 0 4.5"
    OUTPUT "into 60 cell\\(s\\)"
    CHECK ${check_is} dry --q 0.5 --line 3.0 0 0 4.5 --cells 60)
# Rain of 360000 mm/h (0.1 m/s) with the point of
# inflow_source_point_first_step (0.25 m/s into its cell): the step is the
# root of the cubic with the sum of the two rates, 0.35 m/s, not the larger.
flume_add_regtest(inflow_source_point_with_rain NPROCS 2
    INPUTS ${is_dir}/inputs_dry
    ARGS ${is_dry} "flume.inflows.l.location=9.1 6.4" flume.rain.type=constant
         flume.rain.rate=360000
    CHECK ${check_is} dry --q 0.5 --point 9.1 6.4 --rain 360000)

# --- Mass in = the integral of the hydrograph --------------------------------
# The closed basin of inputs_basin: point a takes hydrograph_triangle.txt
# (linear: 200 m3), line b hydrograph_steps.txt (step: 260 m3). At 300 s the
# inflow columns hold the integrals and the basin 800 + 460 m3, to round-off,
# on 4 ranks in 9 boxes. The steps do not land on the times of the files.
flume_add_regtest(inflow_source_basin_volume_np4 NPROCS 4
    INPUTS ${is_basin}
    ARGS ${is_files}
    OUTPUT "inflow a: hydrograph from [^\n]*hydrograph_triangle\\.txt, 3 entries, linear in time, up to 2 m3/s into 1 cell\\(s\\), \\(5, 10\\), at \\(11, 21\\) m"
    CHECK ${check_is_bal} --expect-columns stored inflow_a inflow_b
          --expect-final inflow_a=200 inflow_b=260 stored=1260)
# The other interpolation for each file: 280 m3 and 205 m3, at first order.
flume_add_regtest(inflow_source_basin_volume_swapped NPROCS 1
    INPUTS ${is_basin}
    ARGS ${is_files} flume.inflows.a.interpolation=step flume.inflows.b.interpolation=linear
         flume.order=1
    OUTPUT "inflow b: hydrograph from [^\n]*hydrograph_steps\\.txt, 3 entries, linear in time, up to 3 m3/s into 25 cell\\(s\\)"
    CHECK ${check_is_bal} --expect-final inflow_a=280 inflow_b=205 stored=1285)
# A file whose first time, 50 s, is after the start of the run: its
# discharge, 2 m3/s, enters from t = 0 (600 m3 by 300 s), and a note says so.
# Line b takes no discharge (hydrograph_zero.txt) and no note.
flume_add_regtest(inflow_source_basin_late_start NPROCS 2
    INPUTS ${is_basin}
    ARGS flume.inflows.a.file=${is_dir}/hydrograph_late.txt
         flume.inflows.b.file=${is_dir}/hydrograph_zero.txt
    OUTPUT "inflow a: note: the first time of [^\n]*hydrograph_late\\.txt is 50 s; its discharge, 2 m3/s, enters from the start of the run"
    NO_OUTPUT "inflow b: note"
    CHECK ${check_is_bal} --expect-final inflow_a=600 inflow_b=0 stored=1400)
# With rain on the basin as well (60 mm/h for the whole run: 5 mm, 8 m3 over
# the 1600 m2 of the basin by 300 s): the rain and the inflows are counted
# apart, and the step bound adds their rates.
flume_add_regtest(inflow_source_basin_with_rain NPROCS 3
    INPUTS ${is_basin}
    ARGS ${is_files} flume.rain.type=constant flume.rain.rate=60
    CHECK ${check_is_bal} --expect-columns stored rain inflow_a inflow_b
          --expect-final rain=8 inflow_a=200 inflow_b=260 stored=1268)

# A tiny inflow into deep water (inputs_deep): each step's depth is below
# the rounding of the 1 m depth, and the compensated sum keeps it. After
# 20000 steps the cell holds 1 m plus Q t / A (1.6e-12 m, about 7300
# epsilons) within 4 epsilons of 1 m; without the carry it would hold 1 m.
flume_add_regtest(inflow_source_deep_water NPROCS 1
    INPUTS ${is_dir}/inputs_deep
    ARGS flume.inflows.p.file=${is_dir}/hydrograph_tiny.txt
    CHECK ${check_is} deep --q 1e-15 --point 0.5 0.5 --deep 1.0)

# The basin on 3 ranks in boxes of 4 (25 boxes in place of 9): every value of
# every plotfile is that of a run on 2 ranks in boxes of 8 to the last bit
# (each cell of an inflow takes its depth on its own; one maximum sets the
# step).
flume_add_regtest(inflow_source_basin_reference NPROCS 2
    INPUTS ${is_basin}
    ARGS ${is_files} flume.check_int=100 flume.plot_file=plt
    CHECK ${check_rs} --expect chk00100)
if(TEST inflow_source_basin_reference)
    set_tests_properties(inflow_source_basin_reference PROPERTIES
                         FIXTURES_SETUP inflow_source_basin_reference)
endif()
flume_add_regtest(inflow_source_basin_boxes_4_np3 NPROCS 3
    INPUTS ${is_basin}
    ARGS ${is_files} flume.max_grid_size=4 flume.plot_file=plt
    CHECK ${check_rs} --compare ${FLUME_TEST_RUN_DIR}/inflow_source_basin_reference)
if(TEST inflow_source_basin_boxes_4_np3)
    set_tests_properties(inflow_source_basin_boxes_4_np3 PROPERTIES
                         FIXTURES_REQUIRED inflow_source_basin_reference)
endif()

# --- Restart ----------------------------------------------------------------
# A restart from the checkpoint at step 100 (hydrograph a rising, b at its
# first 1 m3/s) gives every plotfile to the last bit and, on the same ranks,
# the log byte for byte; on another rank count the log agrees to 1e-14 of
# each column. The carry of each inflow's depth and its volume come from the
# checkpoint.
flume_add_restart_test(inflow_source_restart_np2 REFERENCE inflow_source_basin_reference
    INPUTS ${is_basin}
    FIRST_NPROCS 2 RESTART_NPROCS 2
    FIRST_ARGS ${is_files} flume.check_int=100 flume.max_step=150 flume.plot_file=plt
    RESTART_ARGS ${is_files} flume.check_int=100 flume.restart=chk00100 flume.plot_file=plt
    CHECK_ARGS --gauges 0 --rewound)
flume_add_restart_test(inflow_source_restart_np2_np3 REFERENCE inflow_source_basin_reference
    INPUTS ${is_basin}
    FIRST_NPROCS 2 RESTART_NPROCS 3
    FIRST_ARGS ${is_files} flume.check_int=100 flume.max_step=150 flume.plot_file=plt
    RESTART_ARGS ${is_files} flume.check_int=100 flume.restart=chk00100 flume.plot_file=plt
    CHECK_ARGS --gauges 0 --rewound --balance-rtol 1e-14)
# A restart with the point moved to another cell, or with the line's end
# moved, stops: the carry and the volume belong to the cells of the inflow.
set(is_chk flume.restart=${FLUME_TEST_RUN_DIR}/inflow_source_basin_reference/chk00100)
rs_abort(abort_inflow_source_restart_moved inflow_source_basin_reference
    EXPECT "its inflow a:point:5:10 differs from that of the inputs [(]a:point:6:10[)]"
    ARGS ${is_basin} ${is_files} ${is_chk} "flume.inflows.a.location=13 21")
rs_abort(abort_inflow_source_restart_line inflow_source_basin_reference
    EXPECT "its inflow b:line:3\\.[0-9]+:5\\.[0-9]+:25\\.[0-9]+:31\\.[0-9]+ differs from that of the inputs [(]b:line:3\\.[0-9]+:5\\.[0-9]+:25\\.[0-9]+:32[)]"
    ARGS ${is_basin} ${is_files} ${is_chk} "flume.inflows.b.line=3.3 5.1 25.7 32")

# --- Normal depth -----------------------------------------------------------
# A steep channel with Manning friction, fed by 4 m3/s (1 m2/s) along a line
# across it through the first column of cells, the upstream side a wall
# (inputs_channel): the normal depth, 0.394424 m, from 40 m on, as with the
# discharge segment of inflow_normal_depth_steep_order*. The water enters
# with no momentum and the bed accelerates it. At first order the depth is
# 0.77 % high, at second order 0.17 % (0.29 % in the discharge);
# tolerances 1 % and 0.5 %. Cells from 40 m to 190 m.
foreach(case "1|0.01" "2|0.005")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 o)
    list(GET pair 1 tol)
    flume_add_regtest(inflow_source_normal_depth_order${o} NPROCS 2
        INPUTS ${is_dir}/inputs_channel
        ARGS flume.order=${o} flume.inflows.up.file=${is_dir}/hydrograph_steady.txt
        CHECK ${check_is_chan} chan --side xlo --q 1 --mode normal --manning 0.03 --slope 0.02
              --range 40 190 --rtol ${tol})
endforeach()

# --- Aborts -----------------------------------------------------------------
# The keys of an inflow, on the 8 x 8 cells of [0, 40]^2 m, ${rok}.
set(isk "flume\\.inflows\\.f")
set(is_file flume.inflows.f.file=${is_dir}/hydrograph_half.txt)
set(is_pt flume.inflows.names=f flume.inflows.f.location=12\ 12 ${is_file})
flume_add_abort_test(abort_inflow_source_name
    EXPECT "flume\\.inflows\\.names: a name may hold only letters, digits and _; it is \"f-1\""
    ARGS ${rok} flume.inflows.names=f-1)
flume_add_abort_test(abort_inflow_source_name_twice
    EXPECT "flume\\.inflows\\.names lists f twice"
    ARGS ${rok} ${is_pt} flume.inflows.names=f\ f)
flume_add_abort_test(abort_inflow_source_name_case
    EXPECT "flume\\.inflows\\.names lists f and F, which differ only in case"
    ARGS ${rok} ${is_pt} flume.inflows.names=f\ F flume.inflows.F.location=13\ 13
         flume.inflows.F.file=${is_dir}/hydrograph_half.txt)
flume_add_abort_test(abort_inflow_source_location_count
    EXPECT "${isk}\\.location needs exactly 2 values \\(x and y\\); it has 3"
    ARGS ${rok} ${is_pt} flume.inflows.f.location=12\ 12\ 1)
flume_add_abort_test(abort_inflow_source_missing_location
    EXPECT "flume\\.inflows\\.names lists f, which needs ${isk}\\.location \\(x y of a point, m\\) or ${isk}\\.line \\(x0 y0 x1 y1 of a line, m\\)"
    ARGS ${rok} flume.inflows.names=f ${is_file})
flume_add_abort_test(abort_inflow_source_missing_file
    EXPECT "flume\\.inflows\\.names lists f, which needs ${isk}\\.file: lines of a time in s and a discharge in m3/s"
    ARGS ${rok} flume.inflows.names=f flume.inflows.f.location=12\ 12)
flume_add_abort_test(abort_inflow_source_point_and_line
    EXPECT "${isk}\\.location and ${isk}\\.line cannot both be given: an inflow is a point or a line"
    ARGS ${rok} ${is_pt} flume.inflows.f.line=1\ 1\ 2\ 2)
flume_add_abort_test(abort_inflow_source_location_outside
    EXPECT "${isk}\\.location must be inside the domain: from 0 to 40 in y \\(flume\\.prob_lo and flume\\.prob_hi\\); it is 40\\.5"
    ARGS ${rok} ${is_pt} flume.inflows.f.location=12\ 40.5)
flume_add_abort_test(abort_inflow_source_location_not_finite
    EXPECT "${isk}\\.location must be finite in x and y"
    ARGS ${rok} ${is_pt} flume.inflows.f.location=nan\ 12)
flume_add_abort_test(abort_inflow_source_line_three_values
    EXPECT "${isk}\\.line needs exactly 4 values \\(x0 y0 x1 y1\\); it has 3"
    ARGS ${rok} flume.inflows.names=f ${is_file} flume.inflows.f.line=1\ 1\ 2)
flume_add_abort_test(abort_inflow_source_line_not_finite
    EXPECT "${isk}\\.line must be finite: x0 y0 x1 y1, m"
    ARGS ${rok} flume.inflows.names=f ${is_file} flume.inflows.f.line=1\ 1\ inf\ 2)
flume_add_abort_test(abort_inflow_source_line_outside
    EXPECT "${isk}\\.line must be inside the domain: from 0 to 40 in x \\(flume\\.prob_lo and flume\\.prob_hi\\); it is -1"
    ARGS ${rok} flume.inflows.names=f ${is_file} flume.inflows.f.line=1\ 1\ -1\ 2)
flume_add_abort_test(abort_inflow_source_line_one_point
    EXPECT "${isk}\\.line has both ends at the same point: use ${isk}\\.location for a point"
    ARGS ${rok} flume.inflows.names=f ${is_file} flume.inflows.f.line=3\ 4\ 3\ 4)
flume_add_abort_test(abort_inflow_source_interpolation_unknown
    EXPECT "${isk}\\.interpolation must be step or linear; it is cubic"
    ARGS ${rok} ${is_pt} flume.inflows.f.interpolation=cubic)
flume_add_abort_test(abort_inflow_source_interpolation_two_values
    EXPECT "${isk}\\.interpolation needs exactly 1 value"
    ARGS ${rok} ${is_pt} flume.inflows.f.interpolation=step\ linear)
flume_add_abort_test(abort_inflow_source_file_not_found
    EXPECT "${isk}\\.file \\([^)]*no_such_file\\.txt\\): the file cannot be opened"
    ARGS ${rok} flume.inflows.names=f flume.inflows.f.location=12\ 12
         flume.inflows.f.file=${is_dir}/no_such_file.txt)
flume_add_abort_test(abort_inflow_source_file_negative
    EXPECT "${isk}\\.file \\([^)]*\\), line 2: the discharge must be >= 0 m3/s; it is -1"
    ARGS ${rok} flume.inflows.names=f flume.inflows.f.location=12\ 12
         flume.inflows.f.file=${CMAKE_CURRENT_SOURCE_DIR}/Rain/bad/negative.txt)
flume_add_abort_test(abort_inflow_source_unknown_key
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.inflows\\.f\\.locaton"
    ARGS ${rok} ${is_pt} flume.inflows.f.locaton=1\ 1)
flume_add_abort_test(abort_inflow_source_unlisted_name
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.inflows\\.other\\.location"
    ARGS ${rok} ${is_pt} flume.inflows.other.location=1\ 1)
# A point in a building cell, and a line that crosses one: cell (12, 30)
# of the streets of Tests/Raster/streets_heights.asc.
set(is_town flume.n_cell=40\ 40 flume.prob_hi=40\ 40
            flume.buildings.file=${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc)
flume_add_abort_test(abort_inflow_source_point_building
    EXPECT "${isk}\\.location puts inflow f in a building cell \\(12, 30\\) \\(flume\\.buildings\\.file\\)"
    ARGS ${is_town} ${is_pt} flume.inflows.f.location=12.5\ 30.5)
flume_add_abort_test(abort_inflow_source_line_building
    EXPECT "${isk}\\.line puts inflow f in a building cell \\(11, 30\\) \\(flume\\.buildings\\.file\\)"
    ARGS ${is_town} flume.inflows.names=f ${is_file} flume.inflows.f.line=0.5\ 30.5\ 39.5\ 30.5)
# A dry domain, no stop time, and an inflow of no discharge: nothing sets
# the step.
flume_add_abort_test(abort_inflow_source_zero_dry_without_stop_time
    EXPECT "and no rain falls and no water enters through a discharge segment or an inflow [(]flume\\.inflows\\.names[)] from this time on, so nothing sets the time step: set flume\\.stop_time"
    ARGS ${rok} flume.inflows.names=f flume.inflows.f.location=12\ 12
         flume.inflows.f.file=${is_dir}/hydrograph_zero.txt flume.max_step=1)
