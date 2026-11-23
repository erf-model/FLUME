# Tests of gridded rain and of the moving rain cell, included from
# Tests/CMakeLists.txt, whose functions (flume_add_regtest,
# flume_add_abort_test) and ${ok} deck they use. CMAKE_CURRENT_SOURCE_DIR is
# Tests/. See Tests/RainField/README.md.

# Tests/RainField/check_rain_field.py computes the rain from the grids or the
# cell's numbers, independently of the solver: the grids by the area each
# grid cell shares with a solver cell or with the domain, the time by
# splitting [0, t] at the index times, the cell by error functions and
# adaptive quadrature. The decks are Tests/RainField/inputs_flat: a closed,
# flat, dry domain of 24 x 16 cells of 2 m in 6 boxes.
set(check_rain_field ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/RainField/check_rain_field.py)
set(rf_dir ${CMAKE_CURRENT_SOURCE_DIR}/RainField)

# Every cell after one step from a dry start (no fluxes in that step): the
# rain that fell on it. index_offset.txt holds grids at -300 s and 600 s, so
# with linear interpolation the step uses both; grid_a's cells (3.7 x 2.9 m
# from (-5.3, -4.1)) share no edge with the solver's, and its NODATA cell
# lies outside the domain. index_fine.txt holds grid_c alone, 1.3 m cells
# given by their centre: several grid cells per solver cell. Within 16
# epsilons of the largest depth.
foreach(np 1 2)
    flume_add_regtest(rain_gridded_cells_linear_np${np} NPROCS ${np}
        INPUTS ${rf_dir}/inputs_flat
        ARGS flume.rain.file=${rf_dir}/index_offset.txt flume.rain.interpolation=linear
             flume.max_step=1 flume.plot_int=1
        OUTPUT "rain: gridded from [^\n]*index_offset\\.txt, 2 grid\\(s\\), linear in time"
        CHECK ${check_rain_field} flat --index ${rf_dir}/index_offset.txt --linear --cells 1)
endforeach()
flume_add_regtest(rain_gridded_cells_fine NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS flume.rain.file=${rf_dir}/index_fine.txt flume.max_step=1 flume.plot_int=1
    CHECK ${check_rain_field} flat --index ${rf_dir}/index_fine.txt --cells 1)
# A domain of 14 m on 25 cells whose edges meet the 1 m cells of
# grid_aligned.asc, ringed by NODATA cells just outside it: formed in raster
# cells, the last solver edge rounds to 15.000000000000002, a sliver into the
# ring that must not count as read. Every cell within 32 epsilons of the
# largest depth (seen: 17; the 0.56 m edges are not exact in binary).
flume_add_regtest(rain_gridded_cells_aligned_nodata_ring NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS flume.rain.file=${rf_dir}/index_aligned.txt "flume.n_cell=25 25" "flume.prob_hi=14 14"
         flume.max_step=1 flume.plot_int=1
    CHECK ${check_rain_field} flat --index ${rf_dir}/index_aligned.txt --cells 1 --tol-eps 32)
# The water in the closed domain, in every plotfile to 1200 s, against the
# integral of the input over the domain and over the time: three grids on
# two grids of their own (index_steps.txt), held between 0, 300 and 900 s,
# then linear. The second run moves the domain to (10, 5)-(58, 37), so it
# covers another part of the grids. Within 4 epsilons times the water times
# the number of the step (seen: 7e-16 of the water, 3 epsilons, at step 679).
flume_add_regtest(rain_gridded_volume_step NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS flume.rain.file=${rf_dir}/index_steps.txt
    OUTPUT "rain: gridded from [^\n]*index_steps\\.txt, 3 grid\\(s\\), step in time"
    CHECK ${check_rain_field} flat --index ${rf_dir}/index_steps.txt --volume --tol-eps 4
          --n-plotfiles 10)
flume_add_regtest(rain_gridded_volume_linear_shifted NPROCS 1
    INPUTS ${rf_dir}/inputs_flat
    ARGS flume.rain.file=${rf_dir}/index_steps.txt flume.rain.interpolation=linear
         "flume.prob_lo=10 5" "flume.prob_hi=58 37"
    CHECK ${check_rain_field} flat --index ${rf_dir}/index_steps.txt --linear --volume
          --tol-eps 4 --n-plotfiles 10)
# A dry domain with no stop time (Tests/Skeleton/inputs) and grids of no
# rain at 0 s and 101 s around a burst of about 1 mm/s at 100 s, held
# (index_burst.txt): only the burst, a grid ahead of the run's time, sets a
# finite step before 100 s (without it the run stops: nothing sets the
# step). The water in the closed domain against the rain.
flume_add_regtest(rain_gridded_dry_burst_no_stop_time NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.rain.type=gridded flume.rain.file=${rf_dir}/index_burst.txt
         flume.max_step=60 flume.plot_int=10
    CHECK ${check_rain_field} skel --index ${rf_dir}/index_burst.txt --volume --tol-eps 4
          --n-plotfiles 7)

# The moving cell: 120 mm/h at the centre, sigma 6 m. At rest at (40, 26),
# near the north-east corner, so the domain holds part of it: every cell
# after one step within 16 epsilons of the largest depth (the time integral
# is exact for a cell at rest). Moving at (3, 1.5) m/s from (10, 8): the
# centre moves 29 m in the first step, over 20 sub-intervals of the
# Gauss-Legendre rule; every cell within 1e-11 of the largest depth (seen:
# 1.1e-13; with one sub-interval, 1.3e-2). Moving at (0.04, 0.015) m/s from
# (10, 8) to (58, 26) at 1200 s, so it leaves the domain: the water in the
# domain in every plotfile against the integral of the rate over the domain
# and the time, within 4 epsilons times the step (seen: 1.4e-16).
set(rf_cell flume.rain.type=moving_cell flume.rain.cell.peak_rate=120 flume.rain.cell.sigma=6)
foreach(np 1 2)
    flume_add_regtest(rain_cell_at_rest_np${np} NPROCS ${np}
        INPUTS ${rf_dir}/inputs_flat
        ARGS ${rf_cell} "flume.rain.cell.center=40 26" flume.max_step=1 flume.plot_int=1
        OUTPUT "rain: moving cell, peak 120 mm/h, sigma 6 m, centre \\(40, 26\\) m at t = 0, velocity \\(0, 0\\) m/s"
        CHECK ${check_rain_field} flat --cell 120 6 40 26 0 0 --cells 1)
endforeach()
flume_add_regtest(rain_cell_fast_step NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS ${rf_cell} "flume.rain.cell.center=10 8" "flume.rain.cell.velocity=3 1.5"
         flume.max_step=1 flume.plot_int=1
    CHECK ${check_rain_field} flat --cell 120 6 10 8 3 1.5 --cells 1 --rtol 1e-11)
# A cell of no rain moving at 1000 m/s over a dry domain: nothing bounds the
# step, which runs to the stop time, 1e7 s, in one step (1.6e10 sub-intervals
# if they were counted). No water.
flume_add_regtest(rain_cell_zero_peak_long_step NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS flume.rain.type=moving_cell flume.rain.cell.peak_rate=0 flume.rain.cell.sigma=1
         "flume.rain.cell.center=10 8" "flume.rain.cell.velocity=1000 0" flume.stop_time=1e7
    OUTPUT "FLUME: 1 step\\(s\\), time 10000000 s"
    CHECK ${check_rain_field} flat --cell 0 1 10 8 1000 0 --volume)
flume_add_regtest(rain_cell_volume NPROCS 2
    INPUTS ${rf_dir}/inputs_flat
    ARGS ${rf_cell} "flume.rain.cell.center=10 8" "flume.rain.cell.velocity=0.04 0.015"
         flume.plot_int=100
    CHECK ${check_rain_field} flat --cell 120 6 10 8 0.04 0.015 --volume --tol-eps 4
          --n-plotfiles 6)

# Inputs. ${ok} is the 8 x 8 deck of the input checks in Tests/CMakeLists.txt;
# the grids are checked on the 48 x 32 m domain of the decks above (${rf_ok}).
set(rf_bad ${rf_dir}/bad)
set(rf_ok flume.n_cell=24\ 16 flume.prob_hi=48\ 32)
set(gbad "flume\\.rain\\.file \\([^)]*\\)")
flume_add_abort_test(abort_rain_gridded_file_missing_key
    EXPECT "flume\\.rain\\.type = gridded needs flume\\.rain\\.file: an index file, lines of a time in s and an ESRI ASCII grid of the rain rate in mm/h"
    ARGS ${ok} flume.rain.type=gridded)
flume_add_abort_test(abort_rain_gridded_file_empty_name
    EXPECT "flume\\.rain\\.file must not be empty: the name of an index file"
    ARGS ${ok} flume.rain.type=gridded "flume.rain.file=\"\"")
flume_add_abort_test(abort_rain_gridded_index_not_found
    EXPECT "flume\\.rain\\.file \\([^)]*no_such_index\\.txt\\): the file cannot be opened"
    ARGS ${ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/no_such_index.txt)
flume_add_abort_test(abort_rain_gridded_index_directory
    EXPECT "flume\\.rain\\.file \\([^)]*bad\\): the path is not a regular file"
    ARGS ${ok} flume.rain.type=gridded flume.rain.file=${rf_bad})
foreach(case
        "one_item|, line 2: a line holds a time in s and the path of an ESRI ASCII grid of the rain rate in mm/h; it holds 1 item\\(s\\)"
        "three_items|, line 2: a line holds a time in s and the path of an ESRI ASCII grid of the rain rate in mm/h; it holds 3 item\\(s\\)"
        "time_not_a_number|, line 2: the time must be a finite number; it is \"later\""
        "time_infinite|, line 2: the time must be a finite number; it is \"inf\""
        "not_increasing|, line 3: the time must be above the time of the entry before \\(600 s\\); it is 300"
        "empty|: the file holds no entry: it needs at least one line with a time in s and the path of an ESRI ASCII grid of the rain rate in mm/h")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 file)
    list(GET pair 1 expect)
    flume_add_abort_test(abort_rain_gridded_index_${file} EXPECT "${gbad}${expect}"
        ARGS ${ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/index_${file}.txt)
endforeach()
# The grids: a path relative to the index is taken from the index's
# directory; each grid is checked at start-up.
flume_add_abort_test(abort_rain_gridded_grid_missing
    EXPECT "${gbad}, entry 2 at 300 s: grid \\([^)]*bad/no_such_grid\\.asc\\): the file cannot be opened"
    ARGS ${rf_ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/index_missing_grid.txt)
# A rain grid that does not hold the domain (it stops at x = 40 m; the
# domain of ${rf_ok} runs to 48 m in x).
flume_add_abort_test(abort_rain_gridded_not_covering
    EXPECT "${gbad}, entry 2 at 300 s: grid \\([^)]*grid_short\\.asc\\): the raster covers 0 to 40 in x, which does not hold the domain, 0 to 48"
    ARGS ${rf_ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/index_short.txt)
flume_add_abort_test(abort_rain_gridded_negative
    EXPECT "${gbad}, entry 1 at 0 s: grid \\([^)]*grid_negative\\.asc\\): the rain rate must be >= 0 mm/h; it is -2\\.5 in the cell whose west edge is at x = 20 and south edge at y = 12"
    ARGS ${rf_ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/index_negative.txt)
flume_add_abort_test(abort_rain_gridded_nodata
    EXPECT "${gbad}, entry 1 at 0 s: grid \\([^)]*grid_nodata\\.asc\\): [0-9]+ cell\\(s\\) of the domain read a NODATA_value cell of the raster"
    ARGS ${rf_ok} flume.rain.type=gridded flume.rain.file=${rf_bad}/index_nodata.txt)

flume_add_abort_test(abort_rain_cell_peak_missing
    EXPECT "flume\\.rain\\.type = moving_cell needs flume\\.rain\\.cell\\.peak_rate: the rain rate at the centre of the cell in mm/h"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.sigma=5 "flume.rain.cell.center=1 1")
flume_add_abort_test(abort_rain_cell_sigma_missing
    EXPECT "flume\\.rain\\.type = moving_cell needs flume\\.rain\\.cell\\.sigma: the standard deviation of the cell's Gaussian rate in m"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 "flume.rain.cell.center=1 1")
flume_add_abort_test(abort_rain_cell_center_missing
    EXPECT "flume\\.rain\\.type = moving_cell needs flume\\.rain\\.cell\\.center: x and y of the centre of the cell at t = 0 in m"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=5)
flume_add_abort_test(abort_rain_cell_peak_negative
    EXPECT "flume\\.rain\\.cell\\.peak_rate must be >= 0 \\(mm/h\\); it is -1"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=-1 flume.rain.cell.sigma=5
         "flume.rain.cell.center=1 1")
flume_add_abort_test(abort_rain_cell_sigma_zero
    EXPECT "flume\\.rain\\.cell\\.sigma must be > 0 \\(m\\); it is 0"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=0
         "flume.rain.cell.center=1 1")
flume_add_abort_test(abort_rain_cell_center_one_value
    EXPECT "flume\\.rain\\.cell\\.center needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=5
         flume.rain.cell.center=1)
flume_add_abort_test(abort_rain_cell_velocity_inf
    EXPECT "flume\\.rain\\.cell\\.velocity must be finite in x and y"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=5
         "flume.rain.cell.center=1 1" "flume.rain.cell.velocity=0 inf")
foreach(key "peak_rate=10" "sigma=5" "center=1 1" "velocity=1 1")
    string(REGEX REPLACE "=.*" "" name "${key}")
    flume_add_abort_test(abort_rain_cell_${name}_without_type
        EXPECT "flume\\.rain\\.cell\\.${name} is read only with flume\\.rain\\.type = moving_cell"
        ARGS ${ok} flume.rain.type=constant flume.rain.rate=1 "flume.rain.cell.${key}")
endforeach()
# Without any rain type, a cell key is still refused.
flume_add_abort_test(abort_rain_cell_peak_rate_no_type
    EXPECT "flume\\.rain\\.cell\\.peak_rate is read only with flume\\.rain\\.type = moving_cell"
    ARGS ${ok} flume.rain.cell.peak_rate=10)
flume_add_abort_test(abort_rain_cell_file
    EXPECT "flume\\.rain\\.file is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=5
         "flume.rain.cell.center=1 1" flume.rain.file=${rf_dir}/index_fine.txt)
# A faint cell moving at 1000 m/s with sigma 1 m over a dry domain: the
# rain-limited step (160 s) moves it 1.6e5 sigma.
flume_add_abort_test(abort_rain_cell_step_too_long
    EXPECT "flume\\.rain\\.cell: the cell moves [0-9.e+]+ m in one step of [0-9.e+]+ s, more than 2500 times flume\\.rain\\.cell\\.sigma \\(1 m\\)"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=1e-9 flume.rain.cell.sigma=1
         "flume.rain.cell.center=10 8" "flume.rain.cell.velocity=1000 0" flume.stop_time=1e7)
# A gridded run takes flume.rain.interpolation; moving_cell does not.
flume_add_abort_test(abort_rain_cell_interpolation
    EXPECT "flume\\.rain\\.interpolation is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.type=moving_cell flume.rain.cell.peak_rate=10 flume.rain.cell.sigma=5
         "flume.rain.cell.center=1 1" flume.rain.interpolation=linear)
