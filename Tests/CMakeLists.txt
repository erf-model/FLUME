# ---------------------------------------------------------------------------
# Tests. Labels:
#   regression  whole-application runs whose output is checked
#   abort       runs that must stop with a message naming the input
# Every test runs in its own directory under <build>/Tests/run/<name>, so
# tests can run in parallel and nothing is written into the source tree.
# Tests/run_case.py makes the directory, runs the executable and then the
# check; the test fails if either fails.
# ---------------------------------------------------------------------------

set(FLUME_TEST_RUN_DIR ${CMAKE_BINARY_DIR}/Tests/run)
set(FLUME_RUN_CASE ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py)

# Passed to every run. Linux traps invalid floating-point operations; the
# arena starts empty, as in the GPU CI of the sibling codes.
set(FLUME_TEST_ARGS amrex.fpe_trap_invalid=1 amrex.the_arena_init_size=0)

# The launcher, as one ;-separated list: empty without MPI.
function(flume_launcher out nprocs)
    if(FLUME_MPI)
        separate_arguments(preflags NATIVE_COMMAND "${FLUME_MPIEXEC_PREFLAGS}")
        set(${out} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nprocs}
                   ${preflags} ${MPIEXEC_PREFLAGS} PARENT_SCOPE)
        # MPIEXEC_POSTFLAGS is not used.
    else()
        set(${out} "" PARENT_SCOPE)
    endif()
endfunction()

# flume_add_regtest(<name> NPROCS <n> INPUTS <file> [ARGS <key=value>...]
#                   [OUTPUT <regular expression the run's output must match>]
#                   CHECK <command and arguments, run in the test directory>)
# Without MPI only the NPROCS 1 tests are registered.
function(flume_add_regtest name)
    cmake_parse_arguments(T "" "NPROCS;INPUTS;OUTPUT" "ARGS;CHECK" ${ARGN})
    if(NOT T_NPROCS OR NOT T_INPUTS OR NOT T_CHECK)
        message(FATAL_ERROR "flume_add_regtest(${name}): NPROCS, INPUTS and CHECK are required")
    endif()
    if(NOT FLUME_MPI AND NOT T_NPROCS EQUAL 1)
        return()
    endif()
    flume_launcher(launcher ${T_NPROCS})
    set(output "")
    if(T_OUTPUT)
        set(output --expect-output ${T_OUTPUT})
    endif()
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name} ${output}
                --run ${launcher} $<TARGET_FILE:flume> ${T_INPUTS} ${T_ARGS} ${FLUME_TEST_ARGS}
                --check ${T_CHECK})
    set_tests_properties(${name} PROPERTIES LABELS "regression" PROCESSORS ${T_NPROCS} TIMEOUT 120)
endfunction()

# flume_add_abort_test(<name> EXPECT <regular expression> ARGS <key=value>...)
# The run must end with a non-zero exit code and print the expression. The
# keys come from the command line only, so each test states its whole deck.
# The backtrace is not resolved to source lines (amrex.call_addr2line=0),
# which takes seconds per abort.
function(flume_add_abort_test name)
    cmake_parse_arguments(T "" "EXPECT" "ARGS" ${ARGN})
    if(NOT T_EXPECT OR NOT T_ARGS)
        message(FATAL_ERROR "flume_add_abort_test(${name}): EXPECT and ARGS are required")
    endif()
    flume_launcher(launcher 1)
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name}
                --expect-abort ${T_EXPECT}
                --run ${launcher} $<TARGET_FILE:flume> ${T_ARGS} ${FLUME_TEST_ARGS}
                      amrex.call_addr2line=0)
    set_tests_properties(${name} PROPERTIES LABELS "abort" TIMEOUT 120)
endfunction()

if(FLUME_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    if(NOT MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "The tests need an MPI launcher: set MPIEXEC_EXECUTABLE, "
                            "or configure with -DFLUME_ENABLE_TESTS=OFF")
    endif()
endif()

set(check_plotfile ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_plotfile.py)

# --- Skeleton: the grid and the plotfile ------------------------------------
# 48 x 32 cells on [-10, 86] x [5, 69] m, boxes of at most 16 cells: 6 boxes of
# 2 m cells, in a plotfile named skel00000. No water and no stop: the run
# takes no step. The check pins the grid and reads every cell: bed 0, dry.
# Run on 1, 2 and 4 ranks: each rank that owns a box writes one data file, so
# the number of data files shows that the run used the ranks it was given.
set(skel_grid --n-cell 48 32 --prob-lo -10 5 --prob-hi 86 69 --dx 2 2
              --max-grid-size 16 --n-boxes 6)
foreach(np 1 2 4)
    flume_add_regtest(skeleton_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
        CHECK ${check_plotfile} skel --steps 0 --times 0 ${skel_grid} --n-files ${np})
endforeach()

# Defaults: prob_lo = 0 0, max_grid_size = 64 (128 x 65 cells: 4 boxes; any
# other default gives another number), plot_file = plt, bed at 0, dry.
flume_add_regtest(skeleton_defaults NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs_defaults
    CHECK ${check_plotfile} plt --steps 0 --times 0
          --n-cell 128 65 --prob-lo 0 0 --prob-hi 320 162.5 --dx 2.5 2.5
          --max-grid-size 64 --n-boxes 4 --n-files 1)

# --- Uniform state ------------------------------------------------------------
# Tests/Uniform/inputs: depth 2 m, velocity (6, -2) m/s, g = 8 m/s2, 2 m cells,
# CFL 0.5: the rate is (6 + 4) / 2 + (2 + 4) / 2 = 8 1/s, the step is 0.0625 s,
# 16 steps reach the stop time of 1 s, and plot_int = 5 writes steps 0, 5, 10,
# 15 and 16. Every cell of every plotfile must hold the initial state on 1, 2
# and 4 ranks.
set(uni_state --surface 2 --velocity 6 -2)
foreach(np 1 2 4)
    flume_add_regtest(uniform_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
        CHECK ${check_plotfile} uni --steps 0 5 10 15 16 --times 0 0.3125 0.625 0.9375 1
              ${skel_grid} --n-files ${np} ${uni_state})
endforeach()
# The same state with other controls. Defaults: cfl = 0.5 and plot_int = 0
# (first and last plotfile only); g = 9.81 m/s2 would give another step, so g
# stays 8. Six steps of 0.0625 s.
flume_add_regtest(uniform_defaults NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.init.surface_elevation=2 flume.init.velocity=6\ -2 flume.gravity=8 flume.max_step=6
    CHECK ${check_plotfile} skel --steps 0 6 --times 0 0.375 ${skel_grid} --n-files 1 ${uni_state})
# The default gravity, 9.81 m/s2: still water 2.5 m deep (g h = 24.525) on
# 2 m cells, CFL 0.5, one step of 0.5 / (2 sqrt(24.525) / 2) s. The time alone
# has a tolerance that covers single precision; a gravity of 9.8 or 9.82 m/s2
# moves it by 5e-4.
flume_add_regtest(uniform_gravity_default NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.init.surface_elevation=2.5 flume.max_step=1
    CHECK ${check_plotfile} skel --steps 0 1 --times 0 0.10096375546923044 --time-rtol 1e-6
          ${skel_grid} --n-files 1 --surface 2.5)
# CFL 1 and a stop time that is not a whole number of steps: steps of 0.125 s
# to 0.25 s, then one of 0.0625 s that lands on 0.3125 s.
# plot_int = 1 writes every step. The state does not use the step yet, so the
# shortened step is read from the log.
flume_add_regtest(uniform_last_step_cut NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.cfl=1 flume.stop_time=0.3125 flume.plot_int=1 flume.max_step=7
    OUTPUT "Step 3: time 0\\.3125 s, dt 0\\.0625 s"
    CHECK ${check_plotfile} uni --steps 0 1 2 3 --times 0 0.125 0.25 0.3125
          ${skel_grid} --n-files 1 ${uni_state})
# Steps whose sum falls short of the stop time by round-off: still water,
# rate 4 1/s, CFL 0.4, steps of 0.1 s. In double precision ten of them sum to
# 0.9999999999999999 s; the tenth step must land on 1 s, with no eleventh
# step of 1e-16 s. (In single precision the sum passes 1 s.)
flume_add_regtest(uniform_stop_time_roundoff NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.init.velocity=0\ 0 flume.cfl=0.4 flume.plot_int=0
    OUTPUT "FLUME: 10 step\\(s\\), time 1 s"
    CHECK ${check_plotfile} uni --steps 0 10 --times 0 1 ${skel_grid} --n-files 1 --surface 2)
flume_add_regtest(uniform_max_step_first NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.max_step=2
    CHECK ${check_plotfile} uni --steps 0 2 --times 0 0.125 ${skel_grid} --n-files 1 ${uni_state})
# Velocities whose larger component is negative, in x and then in y: the rate
# uses |u| and |v|. Both give (6 + 4) / 2 + (2 + 4) / 2 = 8 1/s and steps of
# 0.0625 s.
flume_add_regtest(uniform_negative_u NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.init.velocity=-6\ 2 flume.max_step=2
    CHECK ${check_plotfile} uni --steps 0 2 --times 0 0.125 ${skel_grid} --n-files 1
          --surface 2 --velocity -6 2)
flume_add_regtest(uniform_negative_v NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.init.velocity=2\ -6 flume.max_step=2
    CHECK ${check_plotfile} uni --steps 0 2 --times 0 0.125 ${skel_grid} --n-files 1
          --surface 2 --velocity 2 -6)
# Cells of 2 m by 1 m and a velocity of (4, 0) m/s: the rate is
# (4 + 4) / 2 = 4 1/s in x plus (0 + 4) / 1 = 4 1/s in y, and the step is
# 0.5 / 8 = 0.0625 s. With the two cell sizes exchanged it would be 0.05 s.
flume_add_regtest(uniform_unequal_cells NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.n_cell=48\ 64 flume.init.velocity=4\ 0 flume.max_step=2
    CHECK ${check_plotfile} uni --steps 0 2 --times 0 0.125
          --n-cell 48 64 --prob-lo -10 5 --prob-hi 86 69 --dx 2 1
          --max-grid-size 16 --n-boxes 12 --n-files 1 --surface 2 --velocity 4 0)
# A stop time or a step limit of zero: no step, one plotfile.
flume_add_regtest(uniform_stop_time_zero NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.stop_time=0
    CHECK ${check_plotfile} uni --steps 0 --times 0 ${skel_grid} --n-files 1 ${uni_state})
flume_add_regtest(uniform_max_step_zero NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.max_step=0
    CHECK ${check_plotfile} uni --steps 0 --times 0 ${skel_grid} --n-files 1 ${uni_state})
# A dry domain has no wave to set a step: one step to the stop time. Cells of
# 2 m by 1 m over a sloping bed, so the bed uses each direction's cell size.
flume_add_regtest(dry_stop_time NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.n_cell=48\ 64 flume.stop_time=5 flume.bed.elevation=-3 flume.bed.slope=0.5\ -0.25
    CHECK ${check_plotfile} skel --steps 0 1 --times 0 5
          --n-cell 48 64 --prob-lo -10 5 --prob-hi 86 69 --dx 2 1
          --max-grid-size 16 --n-boxes 12 --n-files 1
          --bed-elevation -3 --bed-slope 0.5 -0.25)

# --- Sloping bed, partly dry --------------------------------------------------
# Tests/SlopingBed/inputs: a plane bed that rises through a level surface, on
# 50 x 30 cells in 8 boxes of uneven size. The deepest cell sets the step and
# is in one box only, and some boxes are dry, so on 2 and 4 ranks the ranks
# agree on the step only through the reduction of the largest rate over the
# ranks. Every cell is checked against the
# closed form: bed, depth (zero where the bed is above the surface),
# discharges and surface elevation on the absolute datum.
foreach(np 1 2 4)
    flume_add_regtest(sloping_bed_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/SlopingBed/inputs
        CHECK ${check_plotfile} slope --steps 0 3 --times 0 0.1875
              --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
              --max-grid-size 16 --n-boxes 8 --n-files ${np}
              --bed-elevation 100 --bed-slope -0.03125 0.015625
              --surface 98.921875 --velocity 6 -2)
endforeach()

# --- Inputs that must abort ---------------------------------------------------
# Every check in Source/FLUME_Inputs.cpp has at least one test.
set(ok flume.n_cell=8\ 8 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_missing   EXPECT "flume\\.n_cell is required"
    ARGS flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_zero      EXPECT "flume\\.n_cell must be >= 1 and <= 1073741824 in x and y; it is 0"
    ARGS flume.n_cell=8\ 0 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_negative  EXPECT "flume\\.n_cell must be >= 1 and <= 1073741824 in x and y; it is -8"
    ARGS flume.n_cell=-8\ 8 flume.prob_hi=8\ 8)
# 2^30 + 1, and 2^32 + 8, which an int would wrap to 8.
flume_add_abort_test(abort_n_cell_above_max EXPECT "flume\\.n_cell must be .* it is 1073741825"
    ARGS flume.n_cell=8\ 1073741825 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_wraps     EXPECT "flume\\.n_cell must be .* it is 4294967304"
    ARGS flume.n_cell=4294967304\ 8 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_three     EXPECT "flume\\.n_cell needs exactly 2 values \\(x and y\\); it has 3"
    ARGS flume.n_cell=8\ 8\ 1 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_prob_lo_one      EXPECT "flume\\.prob_lo needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${ok} flume.prob_lo=0)
flume_add_abort_test(abort_prob_lo_nan      EXPECT "flume\\.prob_lo must be finite"
    ARGS ${ok} flume.prob_lo=nan\ 0)
flume_add_abort_test(abort_prob_hi_missing  EXPECT "flume\\.prob_hi is required"
    ARGS flume.n_cell=8\ 8)
flume_add_abort_test(abort_prob_hi_one      EXPECT "flume\\.prob_hi needs exactly 2 values \\(x and y\\); it has 1"
    ARGS flume.n_cell=8\ 8 flume.prob_hi=8)
flume_add_abort_test(abort_prob_hi_inf      EXPECT "flume\\.prob_hi must be finite"
    ARGS flume.n_cell=8\ 8 flume.prob_hi=inf\ 8)
flume_add_abort_test(abort_prob_hi_equal    EXPECT "flume\\.prob_hi must be greater than flume\\.prob_lo"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=0\ 8 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_prob_hi_below    EXPECT "flume\\.prob_hi must be greater than flume\\.prob_lo"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=9\ 0 flume.prob_hi=8\ 8)
# A domain one floating-point step wide in y, and a coordinate above the
# largest the build allows: the numbers depend on the precision.
if(FLUME_PRECISION STREQUAL "SINGLE")
    set(one_step 1.0000001)
    set(too_big 1e38)
    set(tiny_extent 4e-30)
    set(huge_extent 4e18)
    set(tiny_depth 1e-44)
else()
    set(one_step 1.0000000000000002)
    set(too_big 1e308)
    set(tiny_extent 4e-300)
    set(huge_extent 4e155)
    set(tiny_depth 1e-320)
endif()
flume_add_abort_test(abort_cell_size_unresolved
    EXPECT "flume\\.prob_lo, flume\\.prob_hi and flume\\.n_cell give a cell size in y that .* cannot resolve"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=0\ 1 flume.prob_hi=8\ ${one_step})
flume_add_abort_test(abort_prob_too_large
    EXPECT "flume\\.prob_lo and flume\\.prob_hi must be at most .* in magnitude"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=-${too_big}\ 0 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_max_grid_size_zero  EXPECT "flume\\.max_grid_size must be >= 1 and <= 2147483647; it is 0"
    ARGS ${ok} flume.max_grid_size=0)
# 2^32 + 1, which an int would wrap to 1.
flume_add_abort_test(abort_max_grid_size_wraps EXPECT "flume\\.max_grid_size must be .* it is 4294967297"
    ARGS ${ok} flume.max_grid_size=4294967297)
flume_add_abort_test(abort_max_grid_size_two   EXPECT "flume\\.max_grid_size needs exactly 1 value; it has 2"
    ARGS ${ok} flume.max_grid_size=4\ 2)
flume_add_abort_test(abort_plot_file_empty  EXPECT "flume\\.plot_file must not be empty"
    ARGS ${ok} flume.plot_file=\"\")
flume_add_abort_test(abort_plot_file_two    EXPECT "flume\\.plot_file needs exactly 1 value; it has 2"
    ARGS ${ok} flume.plot_file=a\ b)
flume_add_abort_test(abort_unknown_key      EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.n_cel "
    ARGS ${ok} flume.n_cel=8\ 8)
# A key outside the flume namespace that nothing reads: reported when the run
# ends (amrex.abort_on_unused_inputs, which main() turns on).
flume_add_abort_test(abort_unused_key_other_namespace
    EXPECT "flum\\.n_cell[\\s\\S]*unused ParmParse variables"
    ARGS ${ok} flum.n_cell=8\ 8)

# --- Time loop, bed and initial state ------------------------------------------
flume_add_abort_test(abort_plot_int_negative EXPECT "flume\\.plot_int must be >= 0 and <= 2147483647; it is -1"
    ARGS ${ok} flume.plot_int=-1)
flume_add_abort_test(abort_stop_time_negative EXPECT "flume\\.stop_time must be >= 0; it is -1"
    ARGS ${ok} flume.stop_time=-1)
flume_add_abort_test(abort_stop_time_inf    EXPECT "flume\\.stop_time must be finite"
    ARGS ${ok} flume.stop_time=inf)
flume_add_abort_test(abort_stop_time_two    EXPECT "flume\\.stop_time needs exactly 1 value; it has 2"
    ARGS ${ok} flume.stop_time=1\ 2)
flume_add_abort_test(abort_max_step_negative EXPECT "flume\\.max_step must be >= 0 and <= 2147483647; it is -1"
    ARGS ${ok} flume.max_step=-1)
# 2^32 + 2, which an int would wrap to 2.
flume_add_abort_test(abort_max_step_wraps   EXPECT "flume\\.max_step must be .* it is 4294967298"
    ARGS ${ok} flume.max_step=4294967298)
flume_add_abort_test(abort_cfl_zero         EXPECT "flume\\.cfl must be > 0 and <= 1; it is 0"
    ARGS ${ok} flume.cfl=0)
flume_add_abort_test(abort_cfl_above_one    EXPECT "flume\\.cfl must be > 0 and <= 1; it is 1\\.5"
    ARGS ${ok} flume.cfl=1.5)
flume_add_abort_test(abort_cfl_negative     EXPECT "flume\\.cfl must be > 0 and <= 1; it is -0\\.5"
    ARGS ${ok} flume.cfl=-0.5)
flume_add_abort_test(abort_cfl_nan          EXPECT "flume\\.cfl must be finite"
    ARGS ${ok} flume.cfl=nan)
flume_add_abort_test(abort_gravity_zero     EXPECT "flume\\.gravity must be > 0; it is 0"
    ARGS ${ok} flume.gravity=0)
flume_add_abort_test(abort_gravity_negative EXPECT "flume\\.gravity must be > 0; it is -9\\.8"
    ARGS ${ok} flume.gravity=-9.81)
flume_add_abort_test(abort_gravity_inf      EXPECT "flume\\.gravity must be finite"
    ARGS ${ok} flume.gravity=inf)
flume_add_abort_test(abort_gravity_nan      EXPECT "flume\\.gravity must be finite"
    ARGS ${ok} flume.gravity=nan)
flume_add_abort_test(abort_bed_elevation_nan EXPECT "flume\\.bed\\.elevation must be finite"
    ARGS ${ok} flume.bed.elevation=nan)
flume_add_abort_test(abort_bed_slope_one    EXPECT "flume\\.bed\\.slope needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${ok} flume.bed.slope=0.1)
flume_add_abort_test(abort_bed_slope_inf    EXPECT "flume\\.bed\\.slope must be finite in x and y"
    ARGS ${ok} flume.bed.slope=0\ inf)
flume_add_abort_test(abort_bed_slope_x_nan  EXPECT "flume\\.bed\\.slope must be finite in x and y"
    ARGS ${ok} flume.bed.slope=nan\ 0)
flume_add_abort_test(abort_surface_nan      EXPECT "flume\\.init\\.surface_elevation must be finite"
    ARGS ${ok} flume.init.surface_elevation=nan)
flume_add_abort_test(abort_surface_two      EXPECT "flume\\.init\\.surface_elevation needs exactly 1 value; it has 2"
    ARGS ${ok} flume.init.surface_elevation=1\ 2)
flume_add_abort_test(abort_velocity_one     EXPECT "flume\\.init\\.velocity needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.velocity=1)
flume_add_abort_test(abort_velocity_nan     EXPECT "flume\\.init\\.velocity must be finite in x and y"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.velocity=nan\ 0)
flume_add_abort_test(abort_velocity_y_inf   EXPECT "flume\\.init\\.velocity must be finite in x and y"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.velocity=0\ inf)
flume_add_abort_test(abort_velocity_without_surface
    EXPECT "flume\\.init\\.velocity needs flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.init.velocity=1\ 0)
# A dry domain with a step limit and no stop time: nothing sets the step.
flume_add_abort_test(abort_dry_without_stop_time EXPECT "set flume\\.stop_time"
    ARGS ${ok} flume.max_step=1)
# Finite inputs whose product is not: a slope that takes the bed past the
# largest number, a depth times a velocity that does, and a depth whose wave
# speed does. ${too_big} is 1e308 in double precision and 1e38 in single.
flume_add_abort_test(abort_bed_overflows
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation and flume\\.bed\\.slope"
    ARGS ${ok} flume.bed.slope=${too_big}\ 0)
# A bed far below the surface: the bed is reported, not the depth built on it.
flume_add_abort_test(abort_bed_overflows_below_surface
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation and flume\\.bed\\.slope"
    ARGS ${ok} flume.bed.slope=-${too_big}\ 0 flume.init.surface_elevation=1)
flume_add_abort_test(abort_initial_state_overflows
    EXPECT "The initial depth or discharge is not finite in 64 cell\\(s\\): check flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.init.surface_elevation=${too_big} flume.init.velocity=1e10\ 0)
flume_add_abort_test(abort_wave_speed_overflows
    EXPECT "the wave speed is not finite, in 64 cell\\(s\\) at step 0"
    ARGS ${ok} flume.init.surface_elevation=${too_big} flume.gravity=100 flume.max_step=1)
# A time step that cannot advance the run. Cells so small that cfl / rate
# underflows to zero, and water so shallow on cells so large that it
# overflows. The numbers depend on the precision.
flume_add_abort_test(abort_time_step_zero
    EXPECT "The time step, 0 s, is zero, not finite or too small to advance the time, 0 s, at step 0: check flume\\.cfl"
    ARGS flume.n_cell=4\ 4 flume.prob_hi=${tiny_extent}\ ${tiny_extent}
         flume.init.surface_elevation=1 flume.cfl=1e-30 flume.stop_time=1)
flume_add_abort_test(abort_time_step_not_finite
    EXPECT "The time step, inf s, is zero, not finite or too small to advance the time, 0 s, at step 0:"
    ARGS flume.n_cell=4\ 4 flume.prob_hi=${huge_extent}\ ${huge_extent}
         flume.init.surface_elevation=${tiny_depth} flume.max_step=3)
