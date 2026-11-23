# Tests of prescribed rain, included from Tests/CMakeLists.txt, whose
# functions (flume_add_regtest, flume_add_abort_test) and ${ok} deck they use.
# CMAKE_CURRENT_SOURCE_DIR is Tests/.

# Tests/Rain/check_rain.py computes the rain fallen since t = 0 from the rate
# or the hyetograph file itself. On the closed, flat, dry domain of
# Tests/Rain/inputs_flat (6 boxes) the water stays level and still: every cell
# of every plotfile holds the rain fallen by its time, and both discharges
# are exactly 0. 36 mm/h for 600 s is 6 mm. The domain starts dry, so the
# first step computes no fluxes and the rain still falls. On 1 and 2 ranks.
set(check_rain ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Rain/check_rain.py)
set(rain_dir ${CMAKE_CURRENT_SOURCE_DIR}/Rain)
foreach(np 1 2)
    flume_add_regtest(rain_flat_constant_np${np} NPROCS ${np}
        INPUTS ${rain_dir}/inputs_flat
        OUTPUT "No depth is above flume\\.dry_depth [^\n]*: the fluxes are not computed"
        CHECK ${check_rain} flat --rate 36 --uniform --n-plotfiles 6)
endforeach()
# Hyetographs on the same domain: 36 mm/h for 600 s, 72 mm/h for 600 s, then
# none (held between entries): 18 mm at 1800 s, 6 mm at 600 s; and a
# triangle from 0 to 72 mm/h and back in 1200 s (linear between entries):
# 12 mm at 1800 s, 1.5 mm at 300 s. Every plotfile is checked at its own time,
# on the rising and the falling side and after the rain.
flume_add_regtest(rain_flat_hyetograph_step NPROCS 2
    INPUTS ${rain_dir}/inputs_flat_steps
    ARGS flume.rain.file=${rain_dir}/hyetograph_steps.txt
    OUTPUT "rain: hyetograph from [^\n]*hyetograph_steps\\.txt, step in time, up to 72 mm/h"
    CHECK ${check_rain} flat_steps --hyetograph ${rain_dir}/hyetograph_steps.txt --uniform
          --n-plotfiles 10)
flume_add_regtest(rain_flat_hyetograph_linear NPROCS 1
    INPUTS ${rain_dir}/inputs_flat_triangle
    ARGS flume.rain.file=${rain_dir}/hyetograph_triangle.txt
    OUTPUT "rain: hyetograph from [^\n]*hyetograph_triangle\\.txt, linear in time, up to 72 mm/h"
    CHECK ${check_rain} flat_triangle --hyetograph ${rain_dir}/hyetograph_triangle.txt --linear
          --uniform --n-plotfiles 9)
# The step file read with linear interpolation: another integral (the
# check reads the file the same way).
flume_add_regtest(rain_flat_hyetograph_step_file_linear NPROCS 1
    INPUTS ${rain_dir}/inputs_flat_steps
    ARGS flume.rain.file=${rain_dir}/hyetograph_steps.txt flume.rain.interpolation=linear
    CHECK ${check_rain} flat_steps --hyetograph ${rain_dir}/hyetograph_steps.txt --linear --uniform)
# A file whose first time, 600 s, is after the start of the run, with a rate
# above 0 (hyetograph_late.txt): that rate falls from t = 0, and the run
# says so. 36 mm/h for 1200 s is 12 mm at 1800 s.
flume_add_regtest(rain_flat_hyetograph_late_start NPROCS 1
    INPUTS ${rain_dir}/inputs_flat_steps
    ARGS flume.rain.file=${rain_dir}/hyetograph_late.txt
    OUTPUT "rain: note: the first time of [^\n]*hyetograph_late\\.txt is 600 s; its rate, 36 mm/h, falls from the start of the run"
    CHECK ${check_rain} flat_steps --hyetograph ${rain_dir}/hyetograph_late.txt --uniform)
# Light rain on 1 m of still water (Tests/Rain/inputs_flat_deep): each
# step's rain is below the rounding of the depth in single precision, and
# the compensated sum keeps it. Every cell holds 1 m plus the rain within 4
# epsilons of the depth, whatever the number of steps.
flume_add_regtest(rain_flat_deep_water NPROCS 2
    INPUTS ${rain_dir}/inputs_flat_deep
    CHECK ${check_rain} flat_deep --rate 0.36 --uniform --initial-depth 1.0)
# Rain with no stop time on the dry, closed domain of Tests/Skeleton/inputs:
# the rain sets the step of a dry domain, and three steps are taken.
flume_add_regtest(rain_flat_no_stop_time NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.rain.type=constant flume.rain.rate=36 flume.max_step=3 flume.plot_int=1
    CHECK ${check_rain} skel --rate 36 --uniform --n-plotfiles 4)
# A dry domain, no stop time, and a hyetograph whose rate is 0 at the start
# and at the end, with a burst of 1 mm/s from 100 s to 101 s
# (hyetograph_burst.txt). Only the burst, an entry inside the interval the
# step looks ahead over, sets a finite step before 100 s; the domain gets
# 1 mm and stays still.
flume_add_regtest(rain_dry_burst_no_stop_time NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.rain.type=hyetograph flume.rain.file=${rain_dir}/hyetograph_burst.txt
         flume.max_step=60 flume.plot_int=10
    CHECK ${check_rain} skel --hyetograph ${rain_dir}/hyetograph_burst.txt --uniform --n-plotfiles 7)
# A dry domain whose initial state holds a discharge in its dry cells
# (5e-7 m of water moving at (1, -2) m/s, below flume.dry_depth): the first
# step computes no fluxes and sets the discharges to 0, so the water the
# rain adds starts still. Every plotfile after the first holds the initial
# depth plus the rain and no discharge.
flume_add_regtest(rain_dry_start_discharge NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.init.surface_elevation=5e-7 "flume.init.velocity=1 -2"
         flume.rain.type=constant flume.rain.rate=36 flume.max_step=5 flume.plot_int=1
    CHECK ${check_rain} skel --rate 36 --uniform --initial-depth 5e-7 --skip-initial
          --n-plotfiles 6)
# A rain rate of 1e-300 mm/h on a dry domain with no stop time: the
# rain-limited step is set by the dry-depth term, and its Newton iteration
# converges from the smaller of the two one-term roots (from the larger it
# would need 171, more than the 100 allowed: the run would stop). The step
# is computed in double; in single precision the depth, about 4e-305 m,
# rounds to 0, so the check expects 0 there.
if(FLUME_PRECISION STREQUAL "SINGLE")
    set(rain_tiny_check 0)
else()
    set(rain_tiny_check 1e-300)
endif()
flume_add_regtest(rain_dry_tiny_rate NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.rain.type=constant flume.rain.rate=1e-300 flume.max_step=2 flume.plot_int=1
    CHECK ${check_rain} skel --rate ${rain_tiny_check} --uniform --n-plotfiles 3)
# The closed basin of Tests/Basin/inputs (rough sloping bed, moving water,
# waves off the walls) under the triangle: the water in the basin grows by
# the rain over its area, to round-off.
flume_add_regtest(rain_basin_volume NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.rain.type=hyetograph flume.rain.file=${rain_dir}/hyetograph_triangle.txt
         flume.rain.interpolation=linear
    CHECK ${check_rain} basin --hyetograph ${rain_dir}/hyetograph_triangle.txt --linear
          --added-volume)
# Rain on an impermeable plane (Tests/Rain/inputs_plane) against the
# kinematic-wave solution (Docs/sphinx_doc/Testing.rst), to 900 s, past
# equilibrium at 349 s. First order, every 400 steps: depth L1 error at most
# 5.5 %, outlet discharge within 5.7 % of i L; at equilibrium the outlet
# passes i L to round-off. Bounds 6 % and 6.5 %. Second order at cfl 0.125,
# every 1600 steps: 0.22 % and 1.3 %, bounds 0.5 % and 1.6 % (the outlet
# value is a cell average, not the flux through the side, so no equilibrium
# check). At cfl 0.5 the second-order error is mostly the error of friction
# split from the stages, which halves with the step; at 0.125 a rain rate
# off by 1 % or a Manning coefficient off by 2 % gives a depth L1 error
# above 1 %. The first-order error halves with the cell size.
if(FLUME_PRECISION STREQUAL "SINGLE")
    set(rain_eq_rtol 1e-4)
else()
    set(rain_eq_rtol 1e-9)
endif()
set(rain_kw --rate 50 --kinematic-wave 0.01 0.01 100 --n-plotfiles 4)
flume_add_regtest(rain_plane_kinematic_wave NPROCS 2
    INPUTS ${rain_dir}/inputs_plane
    CHECK ${check_rain} plane ${rain_kw} --depth-l1-max 0.06 --q-out-max 0.065
          --equilibrium-rtol ${rain_eq_rtol})
flume_add_regtest(rain_plane_kinematic_wave_order2 NPROCS 2
    INPUTS ${rain_dir}/inputs_plane
    ARGS flume.order=2 flume.cfl=0.125 flume.plot_int=1600
    CHECK ${check_rain} plane ${rain_kw} --depth-l1-max 0.005 --q-out-max 0.016)

# Rain inputs. ${ok} is the 8 x 8 deck of the input checks above.
set(rain_bad ${rain_dir}/bad)
set(rbad "flume\\.rain\\.file \\([^)]*\\)")
flume_add_abort_test(abort_rain_type_microphysics
    EXPECT "flume\\.rain\\.type = microphysics is reserved for a later version and is not available yet: use none, constant, hyetograph, gridded or moving_cell"
    ARGS ${ok} flume.rain.type=microphysics)
flume_add_abort_test(abort_rain_type_unknown
    EXPECT "flume\\.rain\\.type must be none, constant, hyetograph, gridded or moving_cell; it is drizzle"
    ARGS ${ok} flume.rain.type=drizzle)
flume_add_abort_test(abort_rain_rate_negative
    EXPECT "flume\\.rain\\.rate must be >= 0 \\(mm/h\\); it is -1"
    ARGS ${ok} flume.rain.type=constant flume.rain.rate=-1)
flume_add_abort_test(abort_rain_rate_inf
    EXPECT "flume\\.rain\\.rate must be finite"
    ARGS ${ok} flume.rain.type=constant flume.rain.rate=inf)
flume_add_abort_test(abort_rain_rate_missing
    EXPECT "flume\\.rain\\.type = constant needs flume\\.rain\\.rate: the rain rate in mm/h"
    ARGS ${ok} flume.rain.type=constant)
flume_add_abort_test(abort_rain_rate_without_type
    EXPECT "flume\\.rain\\.rate is read only with flume\\.rain\\.type = constant"
    ARGS ${ok} flume.rain.rate=10)
flume_add_abort_test(abort_rain_rate_with_hyetograph
    EXPECT "flume\\.rain\\.rate is read only with flume\\.rain\\.type = constant"
    ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=${rain_dir}/hyetograph_steps.txt
         flume.rain.rate=10)
flume_add_abort_test(abort_rain_file_missing_key
    EXPECT "flume\\.rain\\.type = hyetograph needs flume\\.rain\\.file: lines of a time in s and a rain rate in mm/h"
    ARGS ${ok} flume.rain.type=hyetograph)
flume_add_abort_test(abort_rain_file_with_constant
    EXPECT "flume\\.rain\\.file is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.type=constant flume.rain.rate=10
         flume.rain.file=${rain_dir}/hyetograph_steps.txt)
flume_add_abort_test(abort_rain_interpolation_with_constant
    EXPECT "flume\\.rain\\.interpolation is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.type=constant flume.rain.rate=10 flume.rain.interpolation=linear)
flume_add_abort_test(abort_rain_file_without_type
    EXPECT "flume\\.rain\\.file is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.file=${rain_dir}/hyetograph_steps.txt)
flume_add_abort_test(abort_rain_interpolation_without_type
    EXPECT "flume\\.rain\\.interpolation is read only with flume\\.rain\\.type = hyetograph or gridded"
    ARGS ${ok} flume.rain.interpolation=linear)
flume_add_abort_test(abort_rain_interpolation_unknown
    EXPECT "flume\\.rain\\.interpolation must be step or linear; it is cubic"
    ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=${rain_dir}/hyetograph_steps.txt
         flume.rain.interpolation=cubic)
flume_add_abort_test(abort_rain_file_empty_name
    EXPECT "flume\\.rain\\.file must not be empty: the name of a hyetograph file"
    ARGS ${ok} flume.rain.type=hyetograph "flume.rain.file=\"\"")
flume_add_abort_test(abort_rain_file_not_found
    EXPECT "flume\\.rain\\.file \\([^)]*no_such_file\\.txt\\): the file cannot be opened"
    ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=${rain_bad}/no_such_file.txt)
# A relative path is taken from the directory the run starts in, which the
# message names.
flume_add_abort_test(abort_rain_file_not_found_relative
    EXPECT "flume\\.rain\\.file \\(no_such_file\\.txt\\): the file cannot be opened \\(a relative path, taken from /"
    ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=no_such_file.txt)
flume_add_abort_test(abort_rain_file_directory
    EXPECT "flume\\.rain\\.file \\([^)]*bad\\): the path is not a regular file"
    ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=${rain_bad})
foreach(case
        "one_item|, line 2: a line holds two numbers, a time in s and a rain rate in mm/h; it holds 1 item\\(s\\)"
        "three_items|, line 2: a line holds two numbers, a time in s and a rain rate in mm/h; it holds 3 item\\(s\\)"
        "time_not_a_number|, line 2: the time must be a finite number; it is \"sixty\""
        "time_beyond_double|, line 2: the time must be a finite number; it is \"1e400\""
        "rate_inf|, line 2: the rain rate must be a finite number; it is \"inf\""
        "negative|, line 2: the rain rate must be >= 0 mm/h; it is -1"
        "not_increasing|, line 3: the time must be above the time of the entry before \\(60 s\\); it is 30"
        "equal_times|, line 3: the time must be above the time of the entry before \\(60 s\\); it is 60"
        "empty|: the file holds no entry: it needs at least one line with a time in s and a rain rate in mm/h")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 file)
    list(GET pair 1 expect)
    flume_add_abort_test(abort_rain_file_${file} EXPECT "${rbad}${expect}"
        ARGS ${ok} flume.rain.type=hyetograph flume.rain.file=${rain_bad}/${file}.txt)
endforeach()
# A misspelled key is reported as unknown, not as the rate that constant
# rain needs.
flume_add_abort_test(abort_rain_key_misspelled
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.rain\\.rat "
    ARGS ${ok} flume.rain.type=constant flume.rain.rat=10)
# A dry domain with no rain and no stop time: nothing sets the step. Rain of
# 0 mm/h is no rain.
flume_add_abort_test(abort_rain_zero_dry_without_stop_time
    EXPECT "and no rain falls from this time on, so nothing sets the time step: set flume\\.stop_time"
    ARGS ${ok} flume.rain.type=constant flume.rain.rate=0 flume.max_step=1)
