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
# 15 and 16. The domain is periodic in x and y, so every cell of every plotfile
# must hold the initial state on 1, 2 and 4 ranks: the fluxes on the two faces
# of a cell are the same number, and their difference is exactly zero.
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
         flume.periodic=1\ 1
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
# plot_int = 1 writes every step. A uniform state does not show the size of
# a step, so the shortened step is read from the log.
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
# A channel: periodic along the flow and walls on its two sides. The flow is
# along the walls and stays uniform. With the two periodic flags exchanged the
# flow would meet a wall. Once along x and once along y, at 8 m/s: the rate is
# (8 + 4) / 2 + 4 / 2 = 8 1/s and four steps of 0.0625 s end at 0.25 s.
flume_add_regtest(uniform_channel_x NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.periodic=1\ 0 flume.init.velocity=8\ 0 flume.max_step=4
    CHECK ${check_plotfile} uni --steps 0 4 --times 0 0.25 ${skel_grid} --n-files 1
          --surface 2 --velocity 8 0)
flume_add_regtest(uniform_channel_y NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.periodic=0\ 1 flume.init.velocity=0\ -8 flume.max_step=4
    CHECK ${check_plotfile} uni --steps 0 4 --times 0 0.25 ${skel_grid} --n-files 2
          --surface 2 --velocity 0 -8)
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
# The same with a stop time near the largest number: the step of a dry domain
# is not bounded, and a state update with it would give 0 times infinity. The
# run only advances the time, and the plotfile holds zeros, not NaN.
# The plotfile holds the time in the precision of the build: 2^127 is near
# the largest single-precision number and is held exactly.
if(FLUME_PRECISION STREQUAL "SINGLE")
    set(huge_time 170141183460469231731687303715884105728)
else()
    set(huge_time 1e308)
endif()
flume_add_regtest(dry_stop_time_huge NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.stop_time=${huge_time}
    CHECK ${check_plotfile} skel --steps 0 1 --times 0 ${huge_time} ${skel_grid} --n-files 2)

# --- Lake at rest --------------------------------------------------------------
# Tests/SlopingBed/inputs: a plane bed that rises through a level surface, on
# 50 x 30 cells in 8 boxes of uneven size, between walls; the water is at
# rest. The deepest cell sets the step and is in one box only, and some boxes
# are dry, so on 2 and 4 ranks the ranks agree on the step only through the
# reduction of the largest rate over the ranks. After three steps every cell
# is checked against the closed form: bed, depth (zero where the bed is above
# the surface), zero discharges and the surface elevation on the absolute
# datum. Every number of the deck is exact in the arithmetic, and the
# tolerance is 4 epsilons of 100 m.
foreach(np 1 2 4)
    flume_add_regtest(sloping_bed_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/SlopingBed/inputs
        CHECK ${check_plotfile} slope --steps 0 3 --times 0 0.375
              --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
              --max-grid-size 16 --n-boxes 8 --n-files ${np}
              --bed-elevation 100 --bed-slope -0.03125 0.015625
              --surface 98.921875)
endforeach()
# Tests/LakeAtRest/inputs: a level surface at 13 m over a rough bed (a plane
# moved by up to 0.5 m in every cell), between walls. Fifty steps; the deepest
# cell has 4.52404296875 m of water, so with g = 9.81 m/s2 and 2 m cells a step
# is 0.5 / sqrt(g h) = 0.07505370572066064 s. Depth and discharges must be the
# initial ones within 16 epsilons of 13 m (4.6e-14 in double precision): the
# bed-slope term balances the difference of the fluxes to round-off. Without
# that term the discharges reach 1 m2/s.
set(lake_grid --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
              --max-grid-size 16 --n-boxes 8
              --bed-elevation 10 --bed-slope 0.01 -0.02 --bed-noise 0.5 --tol-eps 16)
foreach(np 1 2 4)
    flume_add_regtest(lake_at_rest_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
        CHECK ${check_plotfile} lake --steps 0 50 --times 0 3.752685286033032 --time-rtol 1e-6
              ${lake_grid} --n-files ${np} --surface 13)
endforeach()
# The same bed with the surface at 10.2 m: the high parts of the rough bed are
# islands, and wet and dry cells lie side by side over the whole domain. The
# deepest cell has 1.72404296875 m of water: steps of 0.12157980929135076 s.
flume_add_regtest(lake_at_rest_islands NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.init.surface_elevation=10.2
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
          ${lake_grid} --n-files 2 --surface 10.2)
# The islands with a dry tolerance of 0.05 m, which many cells at the shores
# are below. The tolerance acts on the depth a cell presents to a face, which
# is the same on both sides of a face under a level surface, so the two sides
# are cut together and the water stays at rest. The deepest cell still sets
# the step.
flume_add_regtest(lake_at_rest_islands_dry_depth NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.init.surface_elevation=10.2 flume.dry_depth=0.05
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
          ${lake_grid} --n-files 2 --surface 10.2)
# The same lake near elevation zero (bed 0.3 m at the corner, surface 0.7 m),
# where the depths the two cells of a face present differ in the last digit
# on 468 of 3000 faces, with a tolerance equal to one of those depths in
# double precision. The cut is decided once for a face; decided for each side
# alone, one side of that face is cut and the lake moves by 0.34 m in 50 steps.
# The depths that differ in the last digit leave discharges of 16 epsilons of
# the scale; the check allows 64.
flume_add_regtest(lake_at_rest_dry_depth_at_a_face NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.bed.elevation=0.3 flume.init.surface_elevation=0.7
         flume.dry_depth=0.94279296875
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 5.754374485 --time-rtol 1e-6
          --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
          --max-grid-size 16 --n-boxes 8 --bed-elevation 0.3 --bed-slope 0.01 -0.02
          --bed-noise 0.5 --tol-eps 64 --n-files 2 --surface 0.7)

# --- Dam break -------------------------------------------------------------------
# Tests/DamBreak/inputs: still water 2 m and 1 m deep on a flat bed, 400 cells
# of 0.25 m along the flow, compared at 5 s with the exact solution (Docs/sphinx_doc/Testing.rst):
# hm = 1.4538 m, um = 1.3058 m/s, shock speed 4.1831 m/s. The scheme is first
# order; at this cell size the L1 error is 0.8 % of the depth step and 1.6 % of
# the discharge, the constant state is within 0.2 %, and the shock is 0.2 cell
# from its place. The bounds are 2 %, 0.5 % and 1 cell. The check also sums
# the depth: mass is conserved to round-off.
set(check_dam_break ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_dam_break.py)
set(wet_bed_bounds --l1-max 0.02 --plateau-rtol 0.005 --shock-cells 1)
foreach(np 1 2 4)
    flume_add_regtest(dam_break_x_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
        CHECK ${check_dam_break} dam --time 5 --direction x --position 50
              --surface-lo 2 --surface-hi 1 ${wet_bed_bounds})
endforeach()
# The same along y, with the deep water on the high side: the shock runs
# towards low y.
flume_add_regtest(dam_break_y NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=1 flume.init.dam_surface_elevation=2
    CHECK ${check_dam_break} dam --time 5 --direction y --position 50
          --surface-lo 1 --surface-hi 2 ${wet_bed_bounds})
# Over a dry bed (the low side has no water): Ritter's solution (Delestre et
# al., 2013), whose front is at 94.3 m (x) or 5.7 m (y) at 5 s, short of the wall. This is
# the test of the dry-bed wave speeds and of the depth staying non-negative.
# The L1 error is 0.46 % of the depth and 0.39 % of the discharge scale
# (bound 1 %). Beyond the exact front the depth is zero: a cell no deeper
# than flume.dry_depth passes no water on. With flume.dry_depth = 0
# (dam_break_dry_x_no_tolerance) the scheme carries a trace one cell per
# step, 1.3e-13 m there: above zero and below 1e-11 of 2 m.
set(ritter_bounds --l1-max 0.01 --front-max 0)
flume_add_regtest(dam_break_dry_x NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_surface_elevation=0
    CHECK ${check_dam_break} dam --time 5 --direction x --position 50
          --surface-lo 2 --surface-hi 0 ${ritter_bounds})
flume_add_regtest(dam_break_dry_x_no_tolerance NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_surface_elevation=0 flume.dry_depth=0
    CHECK ${check_dam_break} dam --time 5 --direction x --position 50
          --surface-lo 2 --surface-hi 0 --l1-max 0.01 --front-max 1e-11 --front-min 1e-30)
flume_add_regtest(dam_break_dry_y NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=0 flume.init.dam_surface_elevation=2
    CHECK ${check_dam_break} dam --time 5 --direction y --position 50
          --surface-lo 0 --surface-hi 2 ${ritter_bounds})
# The first step of the two runs above has a closed form that holds the
# dry-bed wave speeds (-c and 2 c, c = sqrt(g h)): the face at the dam passes
# 2 c h / 3 of mass and g h^2 / 3 of momentum, and the step is
# cfl / (c / dx + c / dy). With cfl = 0.5, square cells and h = 2 m the first
# dry cell ends 1/3 m deep and the last wet cell 5/3 m deep. The 5 s runs
# bound errors only; this one pins the speeds.
flume_add_regtest(dam_break_dry_first_step_x NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_surface_elevation=0 flume.max_step=1
    CHECK ${check_dam_break} dam --first-step-cfl 0.5 --direction x --position 50
          --surface-lo 2 --surface-hi 0)
flume_add_regtest(dam_break_dry_first_step_y NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=0 flume.init.dam_surface_elevation=2 flume.max_step=1
    CHECK ${check_dam_break} dam --first-step-cfl 0.5 --direction y --position 50
          --surface-lo 0 --surface-hi 2)

# --- Wetting and drying ------------------------------------------------------------
# Water 2 m deep in the one column of cells below the dam, with a dry bed
# beyond it and across the periodic side. The column loses water through both
# faces: 4 cfl / 3 of its depth in the first step without the limit on the
# outflow, more than it holds for flume.cfl above 0.75
# (Docs/sphinx_doc/Model.rst). The run reaches 5 s with no negative depth (the
# run stops on one) and the same mass.
foreach(cfl 0.8 1.0)
    flume_add_regtest(draining_column_cfl_${cfl} NPROCS 2
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
        ARGS flume.periodic=1\ 0 flume.prob_hi=100\ 1000 flume.init.dam_position=0.3
             flume.init.dam_surface_elevation=0 flume.cfl=${cfl}
        CHECK ${check_dam_break} dam --time 5 --direction x --position 0.3
              --surface-lo 2 --surface-hi 0 --mass-only --min-change 1.9 --dry-depth 1e-6)
endforeach()
# The first step of the column at flume.cfl = 1 has a closed form. With
# c = sqrt(g h) and n = cfl dy / (dx + dy) = 0.999 the column would lose
# 4 n h / 3, so its two faces are scaled by 3 / (4 n): the column ends empty
# and each neighbour holds h / 2 = 1 m and, the momentum flux g h^2 / 3 being
# scaled with the mass flux, a discharge of c h / 4 = 2.2147 m2/s away from
# the column. Cells of 0.25 m by 250 m; the column is cell 0, its neighbours
# cells 1 and 399 (across the periodic side).
set(check_cells ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_cells.py)
set(column_cells "")
foreach(j 0 1 2 3)
    list(APPEND column_cells --cell 0 ${j} 0 0 0 --cell 1 ${j} 1 2.2147234590350102 0
                             --cell 399 ${j} 1 -2.2147234590350102 0)
endforeach()
flume_add_regtest(draining_column_first_step NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.periodic=1\ 0 flume.prob_hi=100\ 1000 flume.init.dam_position=0.3
         flume.init.dam_surface_elevation=0 flume.cfl=1 flume.max_step=1
    CHECK ${check_cells} dam00001 ${column_cells} --others-unchanged dam00000)
# Tests/Film/inputs: a tilted surface over a flat bed, on 8 x 4 cells of 1 m: along x the first
# two cells are dry, the third holds a film of 5e-7 m (below the default
# flume.dry_depth) and the others 0.5 m and more. After one step the two dry
# cells hold exactly nothing: the film passes no water on. With the face
# depth cut at zero the second cell would receive 4e-11 m.
set(film_cells "")
foreach(j 0 1 2 3)
    list(APPEND film_cells --cell 0 ${j} 0 0 0 --cell 1 ${j} 0 0 0)
endforeach()
flume_add_regtest(film_passes_no_water NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Film/inputs
    CHECK ${check_cells} film00001 ${film_cells} --tol-eps 0)
# The same with a bed that falls by 0.0009995 m per cell along x: the third
# cell holds 1 mm of water, above the tolerance, but only 5e-7 m of it stands
# above the bed of the second cell. The tolerance acts on that depth, and the
# two dry cells again hold exactly nothing after the step. A cut on the depth
# of the cell would let water through.
flume_add_regtest(film_behind_a_bed_step NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Film/inputs
    ARGS flume.bed.slope=-0.0009995\ 0 flume.init.surface_elevation=0.74850125
    CHECK ${check_cells} film00001 ${film_cells} --tol-eps 0)
# A column 9.987529 m deep on a slightly rough bed, periodic both ways, at
# flume.cfl = 0.9332. In double precision the first step leaves two cells
# that the limit empties a rounding error below zero (found by a search over
# random decks); the step sets such a depth to zero and the run goes on.
flume_add_regtest(draining_column_round_off NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.periodic=1\ 1 flume.n_cell=400\ 7 flume.prob_hi=100\ 1000
         flume.init.dam_position=0.3 flume.init.surface_elevation=9.987529
         flume.init.dam_surface_elevation=0 flume.cfl=0.9332 flume.bed.noise=0.001
         flume.stop_time=1
    CHECK ${check_dam_break} dam --time 1 --direction x --position 0.3
          --surface-lo 9.987529 --surface-hi 0 --mass-only --min-change 5 --dry-depth 1e-6)
# A film of 1e-8 m that starts at 1e4 m/s beside still water 1 m deep. The
# film is below the default flume.dry_depth (1e-6 m): it does not enter the
# time step, which the deep water sets at 0.019955 s (2 steps to 0.03125 s,
# a time single precision holds exactly; 1130 with flume.dry_depth = 0), and after the run its cells hold no discharge.
flume_add_regtest(thin_film_time_step NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.surface_elevation=1 flume.init.dam_surface_elevation=1e-8
         flume.init.dam_velocity=1e4\ 0 flume.stop_time=0.03125
    OUTPUT "FLUME: 2 step\\(s\\)"
    CHECK ${check_dam_break} dam --time 0.03125 --direction x --position 50
          --surface-lo 1 --surface-hi 1e-8 --velocity-hi 1e4 0
          --mass-only --min-change 0.05 --dry-depth 1e-6)

# Tests/Bowl/inputs: water in a paraboloid bowl with a moving shoreline,
# against Thacker's periodic solutions (Delestre et al., 2013, section 4.2.2;
# Docs/sphinx_doc/Testing.rst). Planar: a tilted plane surface that turns
# about the axis. The scheme is first order and damps the motion: after half
# a period the L1 error of the depth is 13.5 % of the mean depth on 100 x 100
# cells and 7.3 % on 200 x 200 (discharge 19.7 % and 10.8 %); the bounds of
# the fine run are below the errors of the coarse one. After three periods on
# 100 x 100 cells it is 57 % (discharge 77 %), in 1798 steps.
set(check_bowl ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_bowl.py)
set(planar bowl --shape planar --bowl-depth 0.1 --bowl-radius 1 --eta 0.5
           --time-rtol 1e-6 --dry-depth 1e-6)
flume_add_regtest(bowl_planar_half NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.stop_time=2.242850732733187
    CHECK ${check_bowl} ${planar} --periods 0.5 --l1-max 0.15 0.22 --l1-min 0.05)
flume_add_regtest(bowl_planar_half_fine NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.stop_time=2.242850732733187 flume.n_cell=200\ 200 flume.max_grid_size=64
    CHECK ${check_bowl} ${planar} --periods 0.5 --l1-max 0.085 0.125 --l1-min 0.03)
flume_add_regtest(bowl_planar NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    CHECK ${check_bowl} ${planar} --periods 3 --l1-max 0.62 0.85 --max-steps 2000)
# Radial: a paraboloid surface that rises and falls, period
# 2 pi / sqrt(8 g h0) = 2.2429 s, with r0 = 0.8 m: surface 0.025 m at the
# centre, curvature -0.05625 1/m, at rest. Errors after half a period: 4.2 %
# and 2.8 % on 100 x 100 cells, 2.2 % and 1.8 % on 200 x 200; after three
# periods on 100 x 100 cells 20.7 % and 4.5 %, in 668 steps.
set(radial_deck flume.init.surface_elevation=0.025 flume.init.surface_slope=0\ 0
                flume.init.surface_curvature=-0.05625 flume.init.velocity=0\ 0)
set(radial bowl --shape radial --bowl-depth 0.1 --bowl-radius 1 --r0 0.8
           --time-rtol 1e-6 --dry-depth 1e-6)
flume_add_regtest(bowl_radial_half NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS ${radial_deck} flume.stop_time=1.1214253663665935
    CHECK ${check_bowl} ${radial} --periods 0.5 --l1-max 0.048 0.033 --l1-min 0.015)
flume_add_regtest(bowl_radial_half_fine NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS ${radial_deck} flume.stop_time=1.1214253663665935 flume.n_cell=200\ 200
         flume.max_grid_size=64
    CHECK ${check_bowl} ${radial} --periods 0.5 --l1-max 0.026 0.021 --l1-min 0.008)
flume_add_regtest(bowl_radial NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS ${radial_deck} flume.stop_time=6.72855219819956
    CHECK ${check_bowl} ${radial} --periods 3 --l1-max 0.23 0.055 --max-steps 750)
# A mound in place of the bowl (a negative depth) on cells of 0.04 by 0.05 m
# away from the origin, not advanced: the initial state against the closed
# form a quarter of a period on, when the surface is tilted along y.
flume_add_regtest(bowl_mound_initial NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.bed.bowl_depth=-0.1 flume.bed.bowl_radius=2 flume.n_cell=50\ 80
         flume.prob_lo=10\ -3 flume.prob_hi=12\ 1
         flume.init.surface_elevation=0.00625 flume.init.surface_slope=0\ -0.025
         flume.init.velocity=-0.3501785258978626\ 0 flume.max_step=0
    CHECK ${check_bowl} bowl --shape planar --bowl-depth -0.1 --bowl-radius 2 --eta 0.5
          --initial-only --initial-phase 0.25)

# --- Shear layer -------------------------------------------------------------------
# Tests/Shear/inputs: a uniform, subcritical flow at 2 m/s along x carries a jump of the
# velocity along y from 1 to -1 m/s, in a periodic domain. Depth and discharge
# along the flow must stay uniform to 4 epsilons. The velocity across the flow
# must stay between -1 and 1 m/s, its sum must not change, and its two jumps
# must be within one cell of where the flow has carried them, one of them
# through the periodic side. The momentum along a face is carried by the mass
# flux from its upwind side; taken from the other side the jump grows and the
# bounds fail.
set(check_shear ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_shear.py)
foreach(np 1 4)
    flume_add_regtest(shear_x_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
        CHECK ${check_shear} shear --time 5 --direction x --position 50 --depth 2
              --speed 2 --along-lo 1 --along-hi -1)
endforeach()
# The flow reversed: the upwind side of every face is now its high side.
flume_add_regtest(shear_x_reversed NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
    ARGS flume.init.velocity=-2\ 1 flume.init.dam_velocity=-2\ -1
    CHECK ${check_shear} shear --time 5 --direction x --position 50 --depth 2
          --speed -2 --along-lo 1 --along-hi -1)
# Along y, with the dam across y, jumps of 3 to 0.5 m/s, and the flow towards
# low y.
flume_add_regtest(shear_y NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
    ARGS flume.n_cell=8\ 100 flume.prob_hi=8\ 100 flume.init.dam_direction=y
         flume.init.velocity=3\ -2 flume.init.dam_velocity=0.5\ -2
    CHECK ${check_shear} shear --time 5 --direction y --position 50 --depth 2
          --speed -2 --along-lo 3 --along-hi 0.5)

# --- Mass in a closed basin ------------------------------------------------------
# Tests/Basin/inputs: a dam break with moving water over a rough, sloping bed
# between four walls, run until the waves have come back from the walls. The
# sum of the depth over the cells must not change by more than 16 epsilons of
# itself, and the depth must have changed by at least 0.1 m somewhere.
foreach(np 1 2 4)
    flume_add_regtest(basin_mass_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
        CHECK ${check_dam_break} basin --time 20 --direction x --position 15
              --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
              --mass-only --min-change 0.1)
endforeach()

# --- Initial states, no step ---------------------------------------------------
# The dam of the basin moved to 15.5 m, a cell centre, with another velocity
# beyond it: the cell whose centre is at the dam is on the high side. The
# check compares every cell of the first plotfile with the dam initial state
# (it does so in every dam test; here the dam is not on a face).
flume_add_regtest(dam_initial_on_centre NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.init.dam_position=15.5 flume.init.dam_velocity=-1\ 2 flume.max_step=0
    CHECK ${check_dam_break} basin --direction x --position 15.5 --surface-lo 2 --surface-hi 1.5
          --velocity-lo 0.5 -0.3 --velocity-hi -1 2 --initial-only)
flume_add_regtest(dam_initial_on_centre_y NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.init.dam_position=7.5 flume.init.dam_direction=y flume.init.dam_velocity=-1\ 2
         flume.max_step=0
    CHECK ${check_dam_break} basin --direction y --position 7.5 --surface-lo 2 --surface-hi 1.5
          --velocity-lo 0.5 -0.3 --velocity-hi -1 2 --initial-only)
# Tests/SlopingBed/inputs with moving water, not advanced: the discharges are
# the depth times the velocity in a domain that is partly dry.
flume_add_regtest(sloping_bed_initial NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/SlopingBed/inputs
    ARGS flume.init.velocity=6\ -2 flume.max_step=0
    CHECK ${check_plotfile} slope --steps 0 --times 0
          --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
          --max-grid-size 16 --n-boxes 8 --n-files 2
          --bed-elevation 100 --bed-slope -0.03125 0.015625
          --surface 98.921875 --velocity 6 -2)

# --- Second order ------------------------------------------------------------------
# flume.order = 2: limited linear reconstruction and two stages per step
# (Docs/sphinx_doc/Model.rst). First the rate of convergence on a smooth
# solution, then the tests above repeated at second order.
#
# Tests/Wave/inputs: a simple wave of 0.1 m on water 1 m deep, periodic,
# compared at 1 s with its exact solution (the depth of each point travels
# at u + c; Docs/sphinx_doc/Testing.rst). L1 errors of the depth and of the
# discharge, over the amplitude:
#                100 cells            200 cells            400 cells
#   order 1    0.0353  0.0408
#   order 2    0.00350 0.00413      0.00094 0.00111      0.00025 0.00030
# on square cells (along y: 0.00271 0.00252 and 0.00073 0.00068). The error of
# the second-order runs falls by 3.7 when the cell size is halved, a rate of
# 1.9: the minmod limiter is first order at the highest and at the lowest
# point of the wave.
# The bounds of each second-order run are below half the errors of the run
# with cells twice as large, and below the errors of the first-order run on
# the same cells: a scheme that falls back to first order fails them.
set(check_wave ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_wave.py)
set(wave_x wave --direction x --surface 1 --amplitude 0.1 --time 1 --velocity 0.5 -0.3)
set(wave_y wave --direction y --surface 1 --amplitude 0.1 --time 1 --velocity 0.5 -0.3)
set(wave_deck_y flume.prob_hi=0.4\ 10 flume.init.wave_direction=y)
flume_add_regtest(wave_x_first_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    CHECK ${check_wave} ${wave_x} --l1-max 0.045 0.052 --l1-min 0.025)
flume_add_regtest(wave_x_second_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2
    CHECK ${check_wave} ${wave_x} --l1-max 0.0044 0.0052)
flume_add_regtest(wave_x_second_order_fine NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.n_cell=200\ 8
    CHECK ${check_wave} ${wave_x} --l1-max 0.0012 0.0014)
flume_add_regtest(wave_x_second_order_finest NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.n_cell=400\ 16
    CHECK ${check_wave} ${wave_x} --l1-max 0.00032 0.00038)
# Along y, where the water moves at -0.3 m/s along the wave.
flume_add_regtest(wave_y_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.n_cell=4\ 100 ${wave_deck_y}
    CHECK ${check_wave} ${wave_y} --l1-max 0.0034 0.0032)
flume_add_regtest(wave_y_second_order_fine NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.n_cell=8\ 200 ${wave_deck_y}
    CHECK ${check_wave} ${wave_y} --l1-max 0.00092 0.00086)
# A wave of negative amplitude on cells 250 m wide: the time step is set by
# the direction of the wave alone, at the largest flume.cfl of the second
# order (0.5), in 82 steps: 0.00364 and 0.00428. At flume.cfl = 1, which the reader
# stops, the errors on these cells are 0.0153 and 0.0181.
flume_add_regtest(wave_x_second_order_wide_cells NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.prob_hi=10\ 1000 flume.init.wave_amplitude=-0.1
    CHECK ${check_wave} wave --direction x --surface 1 --amplitude -0.1 --time 1
          --velocity 0.5 -0.3 --l1-max 0.0046 0.0054)
# A periodic domain one cell wide, narrower than its two ghost layers: the
# errors of the run four cells wide.
flume_add_regtest(wave_x_second_order_one_cell_wide NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Wave/inputs
    ARGS flume.order=2 flume.n_cell=100\ 1 flume.prob_hi=10\ 0.1
    CHECK ${check_wave} ${wave_x} --l1-max 0.0044 0.0052)

# The uniform state stays as it is: every reconstructed difference is zero.
# Two ghost layers are exchanged across the periodic sides and between boxes.
foreach(np 1 4)
    flume_add_regtest(uniform_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
        ARGS flume.order=2
        CHECK ${check_plotfile} uni --steps 0 5 10 15 16 --times 0 0.3125 0.625 0.9375 1
              ${skel_grid} --n-files ${np} ${uni_state})
endforeach()
# Water at rest stays at rest to round-off: the difference of the surface is
# zero in every cell, whatever the bed. The decks, steps and tolerances of
# the first-order tests; boxes of 16 cells, so the second ghost layer of a
# box lies in its neighbour, and at a wall mirrors the second cell inside.
foreach(np 1 2 4)
    flume_add_regtest(sloping_bed_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/SlopingBed/inputs
        ARGS flume.order=2
        CHECK ${check_plotfile} slope --steps 0 3 --times 0 0.375
              --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
              --max-grid-size 16 --n-boxes 8 --n-files ${np}
              --bed-elevation 100 --bed-slope -0.03125 0.015625
              --surface 98.921875)
    flume_add_regtest(lake_at_rest_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
        ARGS flume.order=2
        CHECK ${check_plotfile} lake --steps 0 50 --times 0 3.752685286033032 --time-rtol 1e-6
              ${lake_grid} --n-files ${np} --surface 13)
endforeach()
flume_add_regtest(lake_at_rest_islands_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.order=2 flume.init.surface_elevation=10.2
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
          ${lake_grid} --n-files 2 --surface 10.2)
flume_add_regtest(lake_at_rest_islands_dry_depth_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.order=2 flume.init.surface_elevation=10.2 flume.dry_depth=0.05
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
          ${lake_grid} --n-files 2 --surface 10.2)
flume_add_regtest(lake_at_rest_dry_depth_at_a_face_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.order=2 flume.bed.elevation=0.3 flume.init.surface_elevation=0.7
         flume.dry_depth=0.94279296875
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 5.754374485 --time-rtol 1e-6
          --n-cell 50 30 --prob-lo 1000 -60 --prob-hi 1100 0 --dx 2 2
          --max-grid-size 16 --n-boxes 8 --bed-elevation 0.3 --bed-slope 0.01 -0.02
          --bed-noise 0.5 --tol-eps 64 --n-files 2 --surface 0.7)

# The dam breaks of Tests/DamBreak/inputs. Over a wet bed the L1 error is
# 0.24 % of the depth step and 0.50 % of the discharge (first order: 0.80 %
# and 1.59 %), and the constant state is within 0.03 % (0.19 %). The bounds
# are 0.65 % and 0.1 %: the first-order run is above them.
set(wet_bed_bounds_2 --l1-max 0.0065 --plateau-rtol 0.001 --shock-cells 1)
foreach(np 1 4)
    flume_add_regtest(dam_break_x_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
        ARGS flume.order=2
        CHECK ${check_dam_break} dam --time 5 --direction x --position 50
              --surface-lo 2 --surface-hi 1 ${wet_bed_bounds_2})
endforeach()
# Between walls one cell apart: the second ghost layer outside each wall
# mirrors the first layer outside the other. Every value starts as a
# signalling NaN.
flume_add_regtest(dam_break_x_second_order_one_cell_wide NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.n_cell=400\ 1 flume.prob_hi=100\ 0.25 amrex.init_snan=1
    CHECK ${check_dam_break} dam --time 5 --direction x --position 50
          --surface-lo 2 --surface-hi 1 ${wet_bed_bounds_2})
flume_add_regtest(dam_break_y_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=1 flume.init.dam_surface_elevation=2
    CHECK ${check_dam_break} dam --time 5 --direction y --position 50
          --surface-lo 1 --surface-hi 2 ${wet_bed_bounds_2})
# Over a dry bed: 0.13 % and 0.12 % (first order: 0.46 % and 0.39 %), bound
# 0.2 %. No water beyond the exact front.
set(ritter_bounds_2 --l1-max 0.002 --front-max 0)
flume_add_regtest(dam_break_dry_x_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.init.dam_surface_elevation=0
    CHECK ${check_dam_break} dam --time 5 --direction x --position 50
          --surface-lo 2 --surface-hi 0 ${ritter_bounds_2})
flume_add_regtest(dam_break_dry_y_second_order NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=0 flume.init.dam_surface_elevation=2
    CHECK ${check_dam_break} dam --time 5 --direction y --position 50
          --surface-lo 0 --surface-hi 2 ${ritter_bounds_2})
# The draining column: each of the two stages keeps the depth at zero or
# above, and the mass is the same. The cells are 250 m wide, so the step is
# that of the direction of the flow at flume.cfl = 0.5, the largest value of
# the second order.
flume_add_regtest(draining_column_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.periodic=1\ 0 flume.prob_hi=100\ 1000
         flume.init.dam_position=0.3 flume.init.dam_surface_elevation=0 flume.cfl=0.5
    CHECK ${check_dam_break} dam --time 5 --direction x --position 0.3
          --surface-lo 2 --surface-hi 0 --mass-only --min-change 1.9 --dry-depth 1e-6)
# The film of Tests/Film/inputs. In the first stage the film (third cell)
# passes no water on and takes water from the fourth cell; in the second it
# is wet and passes water to the second cell. The first cell holds exactly
# nothing after the step. With the face depth cut at zero the second cell
# would pass on in the second stage the trace it took in the first.
set(film_cells_2 "")
foreach(j 0 1 2 3)
    list(APPEND film_cells_2 --cell 0 ${j} 0 0 0)
endforeach()
flume_add_regtest(film_passes_no_water_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Film/inputs
    ARGS flume.order=2
    CHECK ${check_cells} film00001 ${film_cells_2} --tol-eps 0)
flume_add_regtest(film_behind_a_bed_step_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Film/inputs
    ARGS flume.order=2 flume.bed.slope=-0.0009995\ 0 flume.init.surface_elevation=0.74850125
    CHECK ${check_cells} film00001 ${film_cells_2} --tol-eps 0)

# Thacker's solutions in the bowl (Tests/Bowl/inputs). L1 errors of the depth
# and of the discharge (first order in brackets):
#   planar, half a period, 100 x 100   1.71 %  2.12 %   (13.5 %  19.7 %)
#   planar, half a period, 200 x 200   0.61 %  0.71 %   ( 7.3 %  10.8 %)
#   planar, three periods, 100 x 100   3.67 %  5.56 %   (56.5 %  76.8 %)
#   radial, half a period, 100 x 100   1.02 %  0.71 %   ( 4.2 %   2.8 %)
#   radial, half a period, 200 x 200   0.29 %  0.20 %   ( 2.2 %   1.8 %)
#   radial, three periods, 100 x 100   3.80 %  1.54 %   (20.7 %   4.5 %)
# The bounds of the fine runs are below half the errors of the coarse ones,
# and every bound is below the first-order error of the same run.
flume_add_regtest(bowl_planar_half_second_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2 flume.stop_time=2.242850732733187
    CHECK ${check_bowl} ${planar} --periods 0.5 --l1-max 0.022 0.027 --l1-min 0.005)
flume_add_regtest(bowl_planar_half_fine_second_order NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2 flume.stop_time=2.242850732733187 flume.n_cell=200\ 200
         flume.max_grid_size=64
    CHECK ${check_bowl} ${planar} --periods 0.5 --l1-max 0.008 0.0092 --l1-min 0.002)
flume_add_regtest(bowl_planar_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2
    CHECK ${check_bowl} ${planar} --periods 3 --l1-max 0.047 0.07 --max-steps 2100)
flume_add_regtest(bowl_radial_half_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2 ${radial_deck} flume.stop_time=1.1214253663665935
    CHECK ${check_bowl} ${radial} --periods 0.5 --l1-max 0.013 0.0092 --l1-min 0.003)
flume_add_regtest(bowl_radial_half_fine_second_order NPROCS 4
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2 ${radial_deck} flume.stop_time=1.1214253663665935 flume.n_cell=200\ 200
         flume.max_grid_size=64
    CHECK ${check_bowl} ${radial} --periods 0.5 --l1-max 0.0038 0.0026 --l1-min 0.001)
flume_add_regtest(bowl_radial_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Bowl/inputs
    ARGS flume.order=2 ${radial_deck} flume.stop_time=6.72855219819956
    CHECK ${check_bowl} ${radial} --periods 3 --l1-max 0.048 0.02 --max-steps 750)

# The shear layer: the velocity across the flow gets a limited difference of
# its own, and must stay between its two values.
foreach(np 1 4)
    flume_add_regtest(shear_x_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
        ARGS flume.order=2
        CHECK ${check_shear} shear --time 5 --direction x --position 50 --depth 2
              --speed 2 --along-lo 1 --along-hi -1 --width-max 10)
endforeach()
flume_add_regtest(shear_x_reversed_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
    ARGS flume.order=2 flume.init.velocity=-2\ 1 flume.init.dam_velocity=-2\ -1
    CHECK ${check_shear} shear --time 5 --direction x --position 50 --depth 2
          --speed -2 --along-lo 1 --along-hi -1 --width-max 10)
flume_add_regtest(shear_y_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Shear/inputs
    ARGS flume.order=2 flume.n_cell=8\ 100 flume.prob_hi=8\ 100 flume.init.dam_direction=y
         flume.init.velocity=3\ -2 flume.init.dam_velocity=0.5\ -2
    CHECK ${check_shear} shear --time 5 --direction y --position 50 --depth 2
          --speed -2 --along-lo 3 --along-hi 0.5 --width-max 10)
# The film below flume.dry_depth beside still water: the dry rule acts in
# each of the two stages.
flume_add_regtest(thin_film_time_step_second_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.init.surface_elevation=1 flume.init.dam_surface_elevation=1e-8
         flume.init.dam_velocity=1e4\ 0 flume.stop_time=0.03125
    OUTPUT "FLUME: 2 step\\(s\\)"
    CHECK ${check_dam_break} dam --time 0.03125 --direction x --position 50
          --surface-lo 1 --surface-hi 1e-8 --velocity-hi 1e4 0
          --mass-only --min-change 0.05 --dry-depth 1e-6)
# Mass in the closed basin: two ghost layers at the four walls, rough bed.
foreach(np 1 2 4)
    flume_add_regtest(basin_mass_second_order_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
        ARGS flume.order=2
        CHECK ${check_dam_break} basin --time 20 --direction x --position 15
              --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
              --mass-only --min-change 0.1)
endforeach()
# Boxes two cells wide on an odd number of cells: the last box along each
# direction is one cell wide, so the box before it holds the first layer
# outside the wall in its second ghost layer. That layer is read by the last
# face of the box; left as it was allocated, the two boxes of a face disagree
# on its flux and the mass changes (by 0.9 m3 of 1267 here).
flume_add_regtest(basin_mass_second_order_thin_boxes NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.order=2 flume.n_cell=41\ 25 flume.prob_hi=41\ 25 flume.max_grid_size=2
         amrex.init_snan=1
    CHECK ${check_dam_break} basin --time 20 --direction x --position 15
          --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
          --mass-only --min-change 0.1)
flume_add_regtest(basin_mass_second_order_one_cell_boxes NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.order=2 flume.n_cell=20\ 12 flume.prob_hi=40\ 24 flume.max_grid_size=1
         amrex.init_snan=1
    CHECK ${check_dam_break} basin --time 20 --direction x --position 15
          --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
          --mass-only --min-change 0.1)

# --- Manning friction ---------------------------------------------------------------
# Tests/Friction/inputs: water 1 m deep starts at rest in an endless channel
# whose bed falls by 1 mm per metre (flume.bed.continue_slope = 1), with
# n = 0.03 s/m^(1/3). The depth stays 1 m, and the discharge follows
# q_n tanh(t / T_r), q_n = 1.05409 m2/s from Manning's formula and
# T_r = 107.451 s (Tests/check_friction.py; Docs/sphinx_doc/Testing.rst).
# Largest relative error of the discharge in a cell, at 100 s unless noted:
#   order 1                          5.64e-4
#   order 1, steady (1500 s)         2.50e-4  (the bed-slope term of the
#                                             first-order hydrostatic
#                                             reconstruction drives with
#                                             g S (h - S dx / 2))
#   order 2                          1.66e-4
#   order 2, half the step           8.31e-5  (friction is first order in time)
#   order 2, steady (1500 s)         1.2e-13
#   order 2, slope along x and y     1.61e-4
#   order 2, starting at 2 m/s       1.46e-4  (slows down to q_n)
set(check_friction ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_friction.py)
set(fric_x --depth 1 --slope -0.001 0 --manning 0.03)
flume_add_regtest(friction_channel NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    OUTPUT "FLUME: sources: manning_friction"
    CHECK ${check_friction} fric --time 100 ${fric_x} --q-rtol 6e-4)
flume_add_regtest(friction_channel_steady NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.stop_time=1500
    CHECK ${check_friction} fric --time 1500 ${fric_x} --q-rtol 2.5004e-4 --q-rtol-min 2.5002e-4
          --q-rtol-eps 64)
flume_add_regtest(friction_channel_second_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2
    CHECK ${check_friction} fric --time 100 ${fric_x} --q-rtol 2e-4 --q-rtol-eps 1024)
flume_add_regtest(friction_channel_second_order_half_step NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.cfl=0.25
    CHECK ${check_friction} fric --time 100 ${fric_x} --q-rtol 1e-4 --q-rtol-eps 1024)
flume_add_regtest(friction_channel_second_order_steady NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.stop_time=1500
    CHECK ${check_friction} fric --time 1500 ${fric_x} --q-rtol 0 --q-rtol-eps 1024)
flume_add_regtest(friction_channel_second_order_diagonal NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.bed.slope=-0.0006\ 0.0008 flume.init.surface_slope=-0.0006\ 0.0008
         flume.init.surface_elevation=10.9992
    CHECK ${check_friction} fric --time 100 --depth 1 --slope -0.0006 0.0008 --manning 0.03
          --q-rtol 2e-4 --q-rtol-eps 1024)
flume_add_regtest(friction_channel_second_order_y NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.bed.slope=0\ 0.001 flume.init.surface_slope=0\ 0.001
         flume.init.surface_elevation=11.002
    CHECK ${check_friction} fric --time 100 --depth 1 --slope 0 0.001 --manning 0.03
          --q-rtol 2e-4 --q-rtol-eps 1024)
# Side walls along y (periodic in x only), with SNaN fills and the trap on
# invalid operations: the corner ghost cells beyond the periodic side and a
# wall are not touched by the continued slope. Slip walls do not change
# uniform flow along them.
flume_add_regtest(friction_channel_second_order_walls NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.periodic=1\ 0 amrex.init_snan=1 amrex.fpe_trap_invalid=1
    CHECK ${check_friction} fric --time 100 ${fric_x} --q-rtol 2e-4 --q-rtol-eps 1024)
# Gravity of 5 m/s2: q_n does not change, T_r = 210.819 s.
flume_add_regtest(friction_channel_second_order_gravity NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.gravity=5
    CHECK ${check_friction} fric --time 100 ${fric_x} --gravity 5 --q-rtol 1.2e-4
          --q-rtol-eps 1024)
flume_add_regtest(friction_channel_second_order_slowing NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.init.velocity=2\ 0
    CHECK ${check_friction} fric --time 100 ${fric_x} --speed0 2 --q-rtol 2e-4 --q-rtol-eps 1024)
# A film 1 mm deep on a slope of 1 in 10 with n = 0.05: q_n = 6.32e-5 m2/s
# and T_r = 0.0645 s, while the CFL step is about 2 s (2.5 s for the first
# step). Five steps reach the normal discharge (error 1.1e-9); an explicit
# friction term, at dt g n^2 |u| / h^(4/3) = 31, would reverse and grow it.
# At first order the bed-slope term on cells of 1 m (S dx = 100 h) drives
# this film with g h^2 / (2 dx), 0.5 % of g h S: the test is second order
# only.
flume_add_regtest(friction_thin_film_second_order NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.stop_time=10 flume.friction.manning=0.05
         flume.bed.slope=-0.1\ 0 flume.init.surface_slope=-0.1\ 0
         flume.init.surface_elevation=9.601
    OUTPUT "FLUME: 5 step\\(s\\)"
    CHECK ${check_friction} fric --time 10 --depth 0.001 --slope -0.1 0 --manning 0.05
          --q-rtol 1e-8 --q-rtol-eps 4096)
# A film of uneven depth (1.05 to 1.75 mm; the surface falls by 0.0999 per
# metre over a bed falling by 0.1) on the same slope at second order: the
# second stage starts from a discharge friction has not damped yet. The run
# stays positive, keeps its mass, and follows the local Manning discharge of
# the surface slope within 2.4e-3 at 10 s (6 steps).
set(check_film ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_film_local_manning.py)
flume_add_regtest(friction_film_uneven_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Friction/inputs
    ARGS flume.order=2 flume.stop_time=10 flume.friction.manning=0.05
         flume.bed.slope=-0.1\ 0 flume.init.surface_slope=-0.0999\ 0
         flume.init.surface_elevation=9.601
    CHECK ${check_film} fric --time 10 --bed-slope -0.1 --manning 0.05 --rtol 3e-3)
# Guards, not checks of the friction formula (they pass with friction off):
# friction does not move still water, nor put a NaN in a dry cell: the lake
# with islands (wet and dry cells side by side) with friction, at both
# orders, as lake_at_rest_islands.
foreach(order 1 2)
    flume_add_regtest(friction_lake_at_rest_islands_order${order} NPROCS 2
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
        ARGS flume.init.surface_elevation=10.2 flume.friction.manning=0.05 flume.order=${order}
        CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
              ${lake_grid} --n-files 2 --surface 10.2)
endforeach()
# Guard: friction does not change the mass, at second order on 3 ranks.
flume_add_regtest(friction_basin_mass_second_order NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS flume.order=2 flume.friction.manning=0.05
    CHECK ${check_dam_break} basin --time 20 --direction x --position 15
          --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
          --mass-only --min-change 0.1)

# --- Boundary segments ----------------------------------------------------------
# flume.bc.*: a side that is not periodic is made of segments, wall or
# outflow; the cells no segment covers are walls (Docs/sphinx_doc/Model.rst).
#
# The dam breaks of Tests/DamBreak/inputs with an outflow side, run past the
# time the waves reach it. The exact solution of an unbounded domain still
# holds: check_dam_break.py --open-side also requires the mass to drop.
# Over a dry bed the front leaves the domain at 5.64 s and the outflow is
# supercritical; at 9 s the rarefaction has not reached the wall behind it
# (11.3 s). L1 errors (first order, second order): 0.45 % and 0.10 % of the
# depth, 0.30 % and 0.07 % of the discharge. With a wall in place of the
# outflow they are 0.99 % and 0.69 %, so the bounds are 0.6 % and 0.2 %.
flume_add_regtest(dam_break_dry_x_outflow NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_surface_elevation=0 flume.stop_time=9
         flume.bc.segments=out flume.bc.out.side=xhi flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 9 --direction x --position 50
          --surface-lo 2 --surface-hi 0 --l1-max 0.006 --front-max 0 --open-side)
flume_add_regtest(dam_break_dry_x_outflow_second_order NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.init.dam_surface_elevation=0 flume.stop_time=9
         flume.bc.segments=out flume.bc.out.side=xhi flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 9 --direction x --position 50
          --surface-lo 2 --surface-hi 0 --l1-max 0.002 --front-max 0 --open-side)
# The same along y, with the water running out through the low side.
flume_add_regtest(dam_break_dry_y_outflow NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=0 flume.init.dam_surface_elevation=2 flume.stop_time=9
         flume.bc.segments=out flume.bc.out.side=ylo flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 9 --direction y --position 50
          --surface-lo 0 --surface-hi 2 --l1-max 0.006 --front-max 0 --open-side)
flume_add_regtest(dam_break_dry_y_outflow_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.n_cell=4\ 400 flume.prob_hi=1\ 100 flume.init.dam_direction=y
         flume.init.surface_elevation=0 flume.init.dam_surface_elevation=2 flume.stop_time=9
         flume.bc.segments=out flume.bc.out.side=ylo flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 9 --direction y --position 50
          --surface-lo 0 --surface-hi 2 --l1-max 0.002 --front-max 0 --open-side)
# Over a wet bed with the dam at 70 m: the shock leaves at 7.17 s and the flow
# behind it is subcritical (Froude number 0.35). As the shock crosses the
# outflow side a weak wave runs back into the constant state: at 10 s the
# depth there is 0.17 % low and the velocity 0.28 % high at first order,
# 0.68 % and 1.97 % at second order. L1 errors: 0.67 % and 0.20 % of the
# depth step, 1.21 % and 0.32 % of the discharge. With a wall in place of the
# outflow the reflected shock gives 6.0 % and 11.0 % (first order). Bounds:
# L1 2 % and 0.5 %, constant state 0.5 % and 2.5 %.
flume_add_regtest(dam_break_x_outflow NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_position=70 flume.stop_time=10
         flume.bc.segments=out flume.bc.out.side=xhi flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 10 --direction x --position 70
          --surface-lo 2 --surface-hi 1 --l1-max 0.02 --plateau-rtol 0.005 --shock-cells 1
          --open-side)
flume_add_regtest(dam_break_x_outflow_second_order NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.init.dam_position=70 flume.stop_time=10
         flume.bc.segments=out flume.bc.out.side=xhi flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 10 --direction x --position 70
          --surface-lo 2 --surface-hi 1 --l1-max 0.005 --plateau-rtol 0.025 --shock-cells 1
          --open-side)
# Water at rest on the rough bed of Tests/LakeAtRest/inputs stays at rest
# with outflow on parts of three sides and an explicit wall segment on the
# fourth: an outflow ghost cell copies the depth and the bed of the boundary
# cell, so the surface stays level across the side. Three ranks, eight boxes,
# every value a signalling NaN at the start; the times are those of the runs
# between walls.
set(lake_segments flume.bc.segments=w\ n\ e\ s
    flume.bc.w.side=xlo flume.bc.w.type=outflow flume.bc.w.range=-60\ -30
    flume.bc.n.side=yhi flume.bc.n.type=outflow flume.bc.n.range=1020\ 1070
    flume.bc.e.side=xhi flume.bc.e.type=outflow
    flume.bc.s.side=ylo flume.bc.s.type=wall flume.bc.s.range=1000\ 1050 amrex.init_snan=1)
flume_add_regtest(lake_at_rest_segments NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS ${lake_segments}
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 3.752685286033032 --time-rtol 1e-6
          ${lake_grid} --n-files 3 --surface 13)
flume_add_regtest(lake_at_rest_segments_second_order NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.order=2 ${lake_segments}
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 3.752685286033032 --time-rtol 1e-6
          ${lake_grid} --n-files 3 --surface 13)
flume_add_regtest(lake_at_rest_islands_segments_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/LakeAtRest/inputs
    ARGS flume.order=2 flume.init.surface_elevation=10.2 ${lake_segments}
    CHECK ${check_plotfile} lake --steps 0 50 --times 0 6.078990464567538 --time-rtol 1e-6
          ${lake_grid} --n-files 2 --surface 10.2)
# One outflow segment over a whole side, and the same cells as two outflow
# segments that meet inside a box: the plotfiles are identical to the last
# bit. 400 x 12 cells in 100 boxes of up to 8 x 8 on three ranks, second
# order, every value a signalling NaN at the start: a boundary cell that no
# piece fills stops the run.
set(split_deck flume.order=2 flume.n_cell=400\ 12 flume.prob_hi=100\ 3 flume.max_grid_size=8
    flume.init.dam_surface_elevation=0 flume.stop_time=7 amrex.init_snan=1)
flume_add_regtest(outflow_whole_side NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS ${split_deck} flume.bc.segments=out flume.bc.out.side=xhi flume.bc.out.type=outflow
    CHECK ${check_dam_break} dam --time 7 --direction x --position 50
          --surface-lo 2 --surface-hi 0 --mass-only --min-change 1 --open-side)
flume_add_regtest(outflow_split_side NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS ${split_deck} flume.bc.segments=b\ a
         flume.bc.a.side=xhi flume.bc.a.type=outflow flume.bc.a.range=0\ 1.1
         flume.bc.b.side=xhi flume.bc.b.type=outflow flume.bc.b.range=1.1\ 3
    OUTPUT "cells 0 to 3 in y"
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_plotfiles.py dam
          ${FLUME_TEST_RUN_DIR}/outflow_whole_side)
# A build without MPI registers no test on more than one rank.
if(TEST outflow_whole_side AND TEST outflow_split_side)
    set_tests_properties(outflow_whole_side PROPERTIES FIXTURES_SETUP outflow_whole)
    set_tests_properties(outflow_split_side PROPERTIES FIXTURES_REQUIRED outflow_whole)
endif()
# A dry-bed dam break, 2 m deep at x < 50 m, on 200 x 40 cells of 0.5 m. The
# side at x = 100 m is an outflow for y >= 10 m and a wall below: away from
# y = 10 m the open rows follow Ritter's solution, and the wall rows hold the
# bore that stops the arriving water (Tests/check_partial_outflow.py). First
# order: the wall is a segment of its own, listed after the outflow, so the
# pieces are sorted along the side. Second order: the wall is the cells that
# no segment covers.
set(partial_deck flume.n_cell=200\ 40 flume.prob_hi=100\ 20 flume.max_grid_size=16 flume.stop_time=9
    flume.init.dam_surface_elevation=0 flume.bc.up.side=xhi flume.bc.up.type=outflow
    flume.bc.up.range=10\ 20)
set(check_partial ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_partial_outflow.py dam
    --time 9 --position 50 --depth 2 --split 10 --margin 5)
flume_add_regtest(partial_outflow NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS ${partial_deck} flume.bc.segments=up\ down flume.bc.down.side=xhi flume.bc.down.type=wall
         flume.bc.down.range=0\ 10
    CHECK ${check_partial} --l1-max 0.01 --bore-rtol 0.06)
flume_add_regtest(partial_outflow_second_order NPROCS 3
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS ${partial_deck} flume.order=2 flume.bc.segments=up
    CHECK ${check_partial} --l1-max 0.003 --bore-rtol 0.06)
# The basin dam break (rough bed sloping in x, water moving at an angle) on
# 20 x 12 cells, periodic in y, with an outflow on part of the side at
# x = 40 m, second order, signalling NaN start: run in one box, where water
# must leave, and in boxes of one cell, where every cell value must be the
# same. Boxes of one cell are thinner than the two ghost layers, and some end
# one cell short of the side.
set(boxes_deck flume.n_cell=20\ 12 flume.periodic=0\ 1 flume.bed.slope=0.01\ 0 flume.order=2
    flume.stop_time=8 flume.bc.segments=s flume.bc.s.side=xhi flume.bc.s.type=outflow
    flume.bc.s.range=4\ 16 amrex.init_snan=1)
flume_add_regtest(outflow_boxes_one_box NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS ${boxes_deck} flume.max_grid_size=64
    CHECK ${check_dam_break} basin --time 8 --direction x --position 15
          --surface-lo 2 --surface-hi 1.5 --velocity-lo 0.5 -0.3 --velocity-hi 0.5 -0.3
          --mass-only --min-change 0.1 --open-side)
flume_add_regtest(outflow_boxes_one_cell NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Basin/inputs
    ARGS ${boxes_deck} flume.max_grid_size=1
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_plotfiles.py basin
          ${FLUME_TEST_RUN_DIR}/outflow_boxes_one_box --cells)
if(TEST outflow_boxes_one_box AND TEST outflow_boxes_one_cell)
    set_tests_properties(outflow_boxes_one_box PROPERTIES FIXTURES_SETUP outflow_one_box)
    set_tests_properties(outflow_boxes_one_cell PROPERTIES FIXTURES_REQUIRED outflow_one_box)
endif()
# The cells of a segment: those whose centre is at the first coordinate of
# flume.bc.<name>.range or above, and below the second. The cell centres are
# at 0.125 + 0.25 j; both ends of segment lo and the start of hi fall on a
# centre.
flume_add_regtest(boundary_segment_cells NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.max_step=0 flume.bc.segments=lo\ hi
         flume.bc.lo.side=ylo flume.bc.lo.type=outflow flume.bc.lo.range=0.125\ 50.125
         flume.bc.hi.side=ylo flume.bc.hi.type=wall flume.bc.hi.range=50.125\ 100
    OUTPUT "(?s)segment lo: outflow on ylo, cells 0 to 199 in x.*segment hi: wall on ylo, cells 200 to 399 in x.*sources: none"
    CHECK ${check_dam_break} dam --initial-only --direction x --position 50
          --surface-lo 2 --surface-hi 1)

# --- Raster input: bed, buildings and Manning map ----------------------------
# ESRI ASCII grids written by Tests/Raster/make_rasters.py (synthetic, no
# downloaded data; Tests/Raster/README.md). Every file key is passed on the
# command line, as an absolute path. Tests/Raster/check_raster.py computes the
# expected bed from the files in double precision. All on one rank.
set(raster ${CMAKE_CURRENT_SOURCE_DIR}/Raster)
set(check_raster ${Python3_EXECUTABLE} ${raster}/check_raster.py)
# A plane z = 2 + 0.03 x - 0.02 y + 0.001 x y given at the centres of 5 m
# cells (xllcenter): bilinear interpolation reproduces it at every cell
# centre. Largest error 8.9e-16 m. The raster reaches beyond the domain, and
# its NODATA cell lies outside the cells the domain reads.
flume_add_regtest(raster_bed_plane NPROCS 1
    INPUTS ${raster}/inputs_bed
    ARGS flume.bed.file=${raster}/bed_plane.asc
    OUTPUT "FLUME: bed from .*bed_plane\\.asc"
    CHECK ${check_raster} bed --plane 2 0.03 -0.02 0.001)
# A coarse raster of 8 x 4 m cells (dx and dy) that holds the domain exactly:
# the cells within half a raster cell of the sides take the edge values.
# Largest error 0.
flume_add_regtest(raster_bed_coarse NPROCS 1
    INPUTS ${raster}/inputs_bed
    ARGS flume.bed.file=${raster}/bed_coarse.asc
    CHECK ${check_raster} bed --bed ${raster}/bed_coarse.asc)
# A dam break into streets between blocks of 3, 5 and 8 m (a height raster
# with NODATA -1 in no cell, its corner a quarter of a raster cell off the
# grid's edges): the bed is the raster falling to the east plus the
# heights, the 256 building cells stay dry (depth exactly 0) in every
# plotfile, the total depth does not change (walls), and water reaches the
# streets beyond x = 14 m (observed: 0.381 m deep at most at 8 s at order 1,
# 0.371 m at order 2; largest bed error 1.4e-17 m).
set(streets flume.bed.file=${raster}/streets_bed.asc
            flume.buildings.file=${raster}/streets_heights.asc)
set(streets_check --bed ${raster}/streets_bed.asc --buildings ${raster}/streets_heights.asc
                  --dry-buildings --wet-beyond 14)
flume_add_regtest(raster_streets NPROCS 1
    INPUTS ${raster}/inputs_streets
    ARGS ${streets}
    OUTPUT "FLUME: buildings from .*streets_heights\\.asc: 256 cell"
    CHECK ${check_raster} str ${streets_check})
flume_add_regtest(raster_streets_second_order NPROCS 1
    INPUTS ${raster}/inputs_streets
    ARGS ${streets} flume.order=2
    CHECK ${check_raster} str ${streets_check})
# The same footprints as a 0/1 mask, raised by buildings.height.
flume_add_regtest(raster_streets_mask NPROCS 1
    INPUTS ${raster}/inputs_streets
    ARGS flume.bed.file=${raster}/streets_bed.asc
         flume.buildings.file=${raster}/streets_mask.asc flume.buildings.height=5
    OUTPUT "FLUME: buildings from .*streets_mask\\.asc: 256 cell"
    CHECK ${check_raster} str --bed ${raster}/streets_bed.asc
          --buildings ${raster}/streets_mask.asc --building-height 5 --dry-buildings --wet-beyond 14)
# The lake at rest among the same buildings, the surface at 4 m: the blocks
# of 3 m are under 1 m of water, those of 5 and 8 m stand dry. Every cell
# keeps the depth max(4 - bed, 0) and no discharge, to round-off (observed:
# 160 dry cells, depth error 8.9e-16 m and discharge 1.6e-14 m^2/s at most).
foreach(order 1 2)
    flume_add_regtest(raster_lake_at_rest_order_${order} NPROCS 1
        INPUTS ${raster}/inputs_lake
        ARGS ${streets} flume.order=${order}
        CHECK ${check_raster} lake --bed ${raster}/streets_bed.asc
              --buildings ${raster}/streets_heights.asc --surface 4 --conserved)
endforeach()
# Buildings of 0.3 m (the mask, raised by buildings.height) under the 1 m
# dam break: water goes over them (observed: 0.189 m deep at most on a
# building cell), and the total depth does not change.
flume_add_regtest(raster_streets_overtop NPROCS 1
    INPUTS ${raster}/inputs_streets
    ARGS flume.bed.file=${raster}/streets_bed.asc
         flume.buildings.file=${raster}/streets_mask.asc flume.buildings.height=0.3
    CHECK ${check_raster} str --bed ${raster}/streets_bed.asc
          --buildings ${raster}/streets_mask.asc --building-height 0.3 --wet-buildings --conserved)
# A file name relative to the directory the run starts in.
file(RELATIVE_PATH bed_relative ${FLUME_TEST_RUN_DIR}/raster_bed_relative_path ${raster}/bed_coarse.asc)
flume_add_regtest(raster_bed_relative_path NPROCS 1
    INPUTS ${raster}/inputs_bed
    ARGS flume.bed.file=${bed_relative}
    CHECK ${check_raster} bed --bed ${raster}/bed_coarse.asc)
# 100 boxes of 4 x 4 cells give the same cells as 9 boxes of up to 16 x 16:
# the rasters are sampled per cell, not per box.
flume_add_regtest(raster_streets_small_boxes NPROCS 1
    INPUTS ${raster}/inputs_streets
    ARGS ${streets} flume.max_grid_size=4
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_plotfiles.py str
          ${FLUME_TEST_RUN_DIR}/raster_streets --cells)
if(TEST raster_streets AND TEST raster_streets_small_boxes)
    set_tests_properties(raster_streets PROPERTIES FIXTURES_SETUP raster_streets)
    set_tests_properties(raster_streets_small_boxes PROPERTIES FIXTURES_REQUIRED raster_streets)
endif()
# The friction channel with a Manning map: one coefficient per row of 1 m
# (0.02, 0.03, 0.045, 0.06 from the south). Each row follows the closed form
# with its own coefficient. Largest relative error 5.63e-4 at order 1 (as
# friction_channel) and 1.91e-4 at order 2.
# The smallest and largest coefficient as the run prints them: in single
# precision 0.02 is 0.01999999955.
set(n_lo "0\\.0(2|1999999[0-9]*)")
set(n_hi "0\\.0(6|5999999[0-9]*)")
set(manning_rows --depth 1 --slope -0.001 0 --manning-rows 0.02 0.03 0.045 0.06)
flume_add_regtest(raster_manning_rows NPROCS 1
    INPUTS ${raster}/inputs_manning
    ARGS flume.friction.manning_file=${raster}/manning_rows.asc
    OUTPUT "FLUME: Manning friction from .*manning_rows\\.asc, n from ${n_lo} to ${n_hi}"
    CHECK ${check_friction} fric --time 100 ${manning_rows} --q-rtol 6e-4)
flume_add_regtest(raster_manning_rows_second_order NPROCS 1
    INPUTS ${raster}/inputs_manning
    ARGS flume.order=2 flume.friction.manning_file=${raster}/manning_rows.asc
    CHECK ${check_friction} fric --time 100 ${manning_rows} --q-rtol 2e-4 --q-rtol-eps 1024)
# A uniform map of 0.03 gives the plotfiles of flume.friction.manning = 0.03,
# bit for bit.
flume_add_regtest(raster_manning_uniform_value NPROCS 1
    INPUTS ${raster}/inputs_manning
    ARGS flume.friction.manning=0.03
    CHECK ${check_friction} fric --time 100 ${fric_x} --q-rtol 6e-4)
flume_add_regtest(raster_manning_uniform_file NPROCS 1
    INPUTS ${raster}/inputs_manning
    ARGS flume.friction.manning_file=${raster}/manning_uniform.asc
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_plotfiles.py fric
          ${FLUME_TEST_RUN_DIR}/raster_manning_uniform_value)
# Ties: every grid centre lies on the edge between two rows of a raster of
# 0.1 m rows, and computed from the part kept its position is a whole
# number less a few epsilons. The row above holds the coefficient of the
# grid row; the row below holds 0.09. Each row follows the closed form with
# its own coefficient, as raster_manning_rows.
flume_add_regtest(raster_manning_ties NPROCS 1
    INPUTS ${raster}/inputs_manning
    ARGS flume.friction.manning_file=${raster}/manning_ties.asc
    OUTPUT "FLUME: Manning friction from .*manning_ties\\.asc, n from ${n_lo} to ${n_hi}"
    CHECK ${check_friction} fric --time 100 ${manning_rows} --q-rtol 6e-4)
if(TEST raster_manning_uniform_value AND TEST raster_manning_uniform_file)
    set_tests_properties(raster_manning_uniform_value PROPERTIES FIXTURES_SETUP raster_manning_value)
    set_tests_properties(raster_manning_uniform_file PROPERTIES FIXTURES_REQUIRED raster_manning_value)
endif()

# --- Inputs that must abort ---------------------------------------------------
# Every check in Source/FLUME_Inputs.cpp has at least one test.
set(ok flume.n_cell=8\ 8 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_missing   EXPECT "flume\\.n_cell is required"
    ARGS flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_zero      EXPECT "flume\\.n_cell must be >= 1 and <= 1073741824 in x and y. it is 0"
    ARGS flume.n_cell=8\ 0 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_negative  EXPECT "flume\\.n_cell must be >= 1 and <= 1073741824 in x and y. it is -8"
    ARGS flume.n_cell=-8\ 8 flume.prob_hi=8\ 8)
# 2^30 + 1, and 2^32 + 8, which an int would wrap to 8.
flume_add_abort_test(abort_n_cell_above_max EXPECT "flume\\.n_cell must be .* it is 1073741825"
    ARGS flume.n_cell=8\ 1073741825 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_wraps     EXPECT "flume\\.n_cell must be .* it is 4294967304"
    ARGS flume.n_cell=4294967304\ 8 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_n_cell_three     EXPECT "flume\\.n_cell needs exactly 2 values \\(x and y\\). it has 3"
    ARGS flume.n_cell=8\ 8\ 1 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_prob_lo_one      EXPECT "flume\\.prob_lo needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${ok} flume.prob_lo=0)
flume_add_abort_test(abort_prob_lo_nan      EXPECT "flume\\.prob_lo must be finite"
    ARGS ${ok} flume.prob_lo=nan\ 0)
flume_add_abort_test(abort_prob_hi_missing  EXPECT "flume\\.prob_hi is required"
    ARGS flume.n_cell=8\ 8)
flume_add_abort_test(abort_prob_hi_one      EXPECT "flume\\.prob_hi needs exactly 2 values \\(x and y\\). it has 1"
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
    set(tiny_radius 1e-20)
    set(near_max 3e38)
else()
    set(one_step 1.0000000000000002)
    set(too_big 1e308)
    set(tiny_extent 4e-300)
    set(huge_extent 4e155)
    set(tiny_depth 1e-320)
    set(tiny_radius 1e-160)
    set(near_max 1.7e308)
endif()
flume_add_abort_test(abort_cell_size_unresolved
    EXPECT "flume\\.prob_lo, flume\\.prob_hi and flume\\.n_cell give a cell size in y that .* cannot resolve"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=0\ 1 flume.prob_hi=8\ ${one_step})
flume_add_abort_test(abort_prob_too_large
    EXPECT "flume\\.prob_lo and flume\\.prob_hi must be at most .* in magnitude"
    ARGS flume.n_cell=8\ 8 flume.prob_lo=-${too_big}\ 0 flume.prob_hi=8\ 8)
flume_add_abort_test(abort_max_grid_size_zero  EXPECT "flume\\.max_grid_size must be >= 1 and <= 2147483647. it is 0"
    ARGS ${ok} flume.max_grid_size=0)
# 2^32 + 1, which an int would wrap to 1.
flume_add_abort_test(abort_max_grid_size_wraps EXPECT "flume\\.max_grid_size must be .* it is 4294967297"
    ARGS ${ok} flume.max_grid_size=4294967297)
flume_add_abort_test(abort_max_grid_size_two   EXPECT "flume\\.max_grid_size needs exactly 1 value. it has 2"
    ARGS ${ok} flume.max_grid_size=4\ 2)
flume_add_abort_test(abort_plot_file_empty  EXPECT "flume\\.plot_file must not be empty"
    ARGS ${ok} flume.plot_file=\"\")
flume_add_abort_test(abort_plot_file_two    EXPECT "flume\\.plot_file needs exactly 1 value. it has 2"
    ARGS ${ok} flume.plot_file=a\ b)
flume_add_abort_test(abort_unknown_key      EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.n_cel "
    ARGS ${ok} flume.n_cel=8\ 8)
# A key outside the flume namespace that nothing reads: reported when the run
# ends (amrex.abort_on_unused_inputs, which main() turns on). The list of keys
# goes to stdout and the abort message to stderr, which reach the pipe in
# either order, so the pattern asks for both in any order.
flume_add_abort_test(abort_unused_key_other_namespace
    EXPECT "^(?=[\\s\\S]*flum\\.n_cell)[\\s\\S]*unused ParmParse variables"
    ARGS ${ok} flum.n_cell=8\ 8)

# --- Time loop, bed and initial state ------------------------------------------
flume_add_abort_test(abort_plot_int_negative EXPECT "flume\\.plot_int must be >= 0 and <= 2147483647. it is -1"
    ARGS ${ok} flume.plot_int=-1)
flume_add_abort_test(abort_stop_time_negative EXPECT "flume\\.stop_time must be >= 0. it is -1"
    ARGS ${ok} flume.stop_time=-1)
flume_add_abort_test(abort_stop_time_inf    EXPECT "flume\\.stop_time must be finite"
    ARGS ${ok} flume.stop_time=inf)
flume_add_abort_test(abort_stop_time_two    EXPECT "flume\\.stop_time needs exactly 1 value. it has 2"
    ARGS ${ok} flume.stop_time=1\ 2)
flume_add_abort_test(abort_max_step_negative EXPECT "flume\\.max_step must be >= 0 and <= 2147483647. it is -1"
    ARGS ${ok} flume.max_step=-1)
# 2^32 + 2, which an int would wrap to 2.
flume_add_abort_test(abort_max_step_wraps   EXPECT "flume\\.max_step must be .* it is 4294967298"
    ARGS ${ok} flume.max_step=4294967298)
flume_add_abort_test(abort_cfl_zero         EXPECT "flume\\.cfl must be > 0 and <= 1. it is 0"
    ARGS ${ok} flume.cfl=0)
flume_add_abort_test(abort_cfl_above_one    EXPECT "flume\\.cfl must be > 0 and <= 1. it is 1\\.5"
    ARGS ${ok} flume.cfl=1.5)
flume_add_abort_test(abort_cfl_negative     EXPECT "flume\\.cfl must be > 0 and <= 1. it is -0\\.5"
    ARGS ${ok} flume.cfl=-0.5)
flume_add_abort_test(abort_cfl_nan          EXPECT "flume\\.cfl must be finite"
    ARGS ${ok} flume.cfl=nan)
flume_add_abort_test(abort_gravity_zero     EXPECT "flume\\.gravity must be > 0. it is 0"
    ARGS ${ok} flume.gravity=0)
flume_add_abort_test(abort_gravity_negative EXPECT "flume\\.gravity must be > 0. it is -9\\.8"
    ARGS ${ok} flume.gravity=-9.81)
flume_add_abort_test(abort_gravity_inf      EXPECT "flume\\.gravity must be finite"
    ARGS ${ok} flume.gravity=inf)
flume_add_abort_test(abort_gravity_nan      EXPECT "flume\\.gravity must be finite"
    ARGS ${ok} flume.gravity=nan)
flume_add_abort_test(abort_bed_elevation_nan EXPECT "flume\\.bed\\.elevation must be finite"
    ARGS ${ok} flume.bed.elevation=nan)
flume_add_abort_test(abort_bed_slope_one    EXPECT "flume\\.bed\\.slope needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${ok} flume.bed.slope=0.1)
flume_add_abort_test(abort_bed_slope_inf    EXPECT "flume\\.bed\\.slope must be finite in x and y"
    ARGS ${ok} flume.bed.slope=0\ inf)
flume_add_abort_test(abort_bed_slope_x_nan  EXPECT "flume\\.bed\\.slope must be finite in x and y"
    ARGS ${ok} flume.bed.slope=nan\ 0)
flume_add_abort_test(abort_surface_nan      EXPECT "flume\\.init\\.surface_elevation must be finite"
    ARGS ${ok} flume.init.surface_elevation=nan)
flume_add_abort_test(abort_surface_two      EXPECT "flume\\.init\\.surface_elevation needs exactly 1 value. it has 2"
    ARGS ${ok} flume.init.surface_elevation=1\ 2)
flume_add_abort_test(abort_velocity_one     EXPECT "flume\\.init\\.velocity needs exactly 2 values \\(x and y\\). it has 1"
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
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation, flume\\.bed\\.slope"
    ARGS ${ok} flume.bed.slope=${too_big}\ 0)
# A bed far below the surface: the bed is reported, not the depth built on it.
flume_add_abort_test(abort_bed_overflows_below_surface
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation, flume\\.bed\\.slope"
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
         flume.init.surface_elevation=${tiny_depth} flume.dry_depth=0 flume.max_step=3)

# --- Wetting and drying, bowl and curved surface ---------------------------------
flume_add_abort_test(abort_dry_depth_negative EXPECT "flume\\.dry_depth must be >= 0. it is -0\\.5"
    ARGS ${ok} flume.dry_depth=-0.5)
flume_add_abort_test(abort_dry_depth_nan    EXPECT "flume\\.dry_depth must be finite"
    ARGS ${ok} flume.dry_depth=nan)
flume_add_abort_test(abort_dry_depth_two    EXPECT "flume\\.dry_depth needs exactly 1 value. it has 2"
    ARGS ${ok} flume.dry_depth=1e-6\ 1e-6)
# Water that is all below the tolerance is a dry domain.
flume_add_abort_test(abort_all_below_dry_depth
    EXPECT "no depth is above flume\\.dry_depth.*set flume\\.stop_time"
    ARGS ${ok} flume.init.surface_elevation=0.5 flume.dry_depth=0.5 flume.max_step=3)
flume_add_abort_test(abort_bowl_depth_inf   EXPECT "flume\\.bed\\.bowl_depth must be finite"
    ARGS ${ok} flume.bed.bowl_depth=inf flume.bed.bowl_radius=1)
flume_add_abort_test(abort_bowl_radius_zero EXPECT "flume\\.bed\\.bowl_radius must be > 0. it is 0"
    ARGS ${ok} flume.bed.bowl_depth=1 flume.bed.bowl_radius=0)
flume_add_abort_test(abort_bowl_radius_negative EXPECT "flume\\.bed\\.bowl_radius must be > 0. it is -2"
    ARGS ${ok} flume.bed.bowl_depth=1 flume.bed.bowl_radius=-2)
flume_add_abort_test(abort_bowl_radius_tiny
    EXPECT "flume\\.bed\\.bowl_radius must have a square that is a normal number. it is "
    ARGS ${ok} flume.bed.bowl_depth=1 flume.bed.bowl_radius=${tiny_radius})
flume_add_abort_test(abort_bowl_radius_nan  EXPECT "flume\\.bed\\.bowl_radius must be finite"
    ARGS ${ok} flume.bed.bowl_depth=1 flume.bed.bowl_radius=nan)
flume_add_abort_test(abort_bowl_depth_without_radius
    EXPECT "flume\\.bed\\.bowl_depth needs flume\\.bed\\.bowl_radius"
    ARGS ${ok} flume.bed.bowl_depth=1)
flume_add_abort_test(abort_bowl_radius_without_depth
    EXPECT "flume\\.bed\\.bowl_radius needs flume\\.bed\\.bowl_depth"
    ARGS ${ok} flume.bed.bowl_radius=1)
flume_add_abort_test(abort_bowl_radius_with_depth_zero
    EXPECT "flume\\.bed\\.bowl_radius needs flume\\.bed\\.bowl_depth"
    ARGS ${ok} flume.bed.bowl_radius=1 flume.bed.bowl_depth=0)
# A radius so small that the bowl overflows at the cell centres.
flume_add_abort_test(abort_bowl_overflows
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check .*flume\\.bed\\.bowl_radius"
    ARGS ${ok} flume.bed.bowl_depth=${near_max} flume.bed.bowl_radius=1e-3)
flume_add_abort_test(abort_surface_slope_one
    EXPECT "flume\\.init\\.surface_slope needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.surface_slope=0.1)
flume_add_abort_test(abort_surface_slope_nan
    EXPECT "flume\\.init\\.surface_slope must be finite in x and y"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.surface_slope=0\ nan)
flume_add_abort_test(abort_surface_curvature_inf
    EXPECT "flume\\.init\\.surface_curvature must be finite"
    ARGS ${ok} flume.init.surface_elevation=1 flume.init.surface_curvature=inf)
flume_add_abort_test(abort_surface_slope_without_surface
    EXPECT "flume\\.init\\.surface_slope needs flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.init.surface_slope=0.1\ 0)
flume_add_abort_test(abort_surface_curvature_without_surface
    EXPECT "flume\\.init\\.surface_curvature needs flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.init.surface_curvature=0.1)
flume_add_abort_test(abort_surface_slope_with_dam
    EXPECT "flume\\.init\\.surface_slope cannot be used with flume\\.init\\.dam_position"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_position=4
         flume.init.dam_surface_elevation=1 flume.init.surface_slope=0.1\ 0)
flume_add_abort_test(abort_surface_curvature_with_dam
    EXPECT "flume\\.init\\.surface_curvature cannot be used with flume\\.init\\.dam_position"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_position=4
         flume.init.dam_surface_elevation=1 flume.init.surface_curvature=0.1)

# --- Order and wave ----------------------------------------------------------------
flume_add_abort_test(abort_cfl_above_half_second_order
    EXPECT "flume\\.cfl must be <= 0\\.5 with flume\\.order = 2. it is 0\\.75"
    ARGS ${ok} flume.order=2 flume.cfl=0.75)
flume_add_abort_test(abort_order_zero       EXPECT "flume\\.order must be >= 1 and <= 2. it is 0"
    ARGS ${ok} flume.order=0)
flume_add_abort_test(abort_order_three      EXPECT "flume\\.order must be >= 1 and <= 2. it is 3"
    ARGS ${ok} flume.order=3)
flume_add_abort_test(abort_order_two_values EXPECT "flume\\.order needs exactly 1 value. it has 2"
    ARGS ${ok} flume.order=1\ 2)
set(wave_ok ${ok} flume.periodic=1\ 1 flume.init.surface_elevation=1)
flume_add_abort_test(abort_wave_amplitude_nan
    EXPECT "flume\\.init\\.wave_amplitude must be finite"
    ARGS ${wave_ok} flume.init.wave_amplitude=nan)
flume_add_abort_test(abort_wave_amplitude_two
    EXPECT "flume\\.init\\.wave_amplitude needs exactly 1 value. it has 2"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1\ 0.1)
flume_add_abort_test(abort_wave_amplitude_without_surface
    EXPECT "flume\\.init\\.wave_amplitude needs flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.periodic=1\ 1 flume.init.wave_amplitude=0.1)
flume_add_abort_test(abort_wave_amplitude_not_periodic
    EXPECT "flume\\.init\\.wave_amplitude needs flume\\.periodic = 1 along flume\\.init\\.wave_direction"
    ARGS ${ok} flume.periodic=0\ 1 flume.init.surface_elevation=1 flume.init.wave_amplitude=0.1)
flume_add_abort_test(abort_wave_amplitude_not_periodic_y
    EXPECT "flume\\.init\\.wave_amplitude needs flume\\.periodic = 1 along flume\\.init\\.wave_direction"
    ARGS ${ok} flume.periodic=1\ 0 flume.init.surface_elevation=1 flume.init.wave_amplitude=0.1
         flume.init.wave_direction=y)
flume_add_abort_test(abort_wave_amplitude_with_slope
    EXPECT "flume\\.init\\.wave_amplitude cannot be used with flume\\.init\\.surface_slope or flume\\.init\\.surface_curvature"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1 flume.init.surface_slope=0.1\ 0)
flume_add_abort_test(abort_wave_amplitude_with_curvature
    EXPECT "flume\\.init\\.wave_amplitude cannot be used with flume\\.init\\.surface_slope or flume\\.init\\.surface_curvature"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1 flume.init.surface_curvature=0.1)
flume_add_abort_test(abort_wave_amplitude_with_dam
    EXPECT "flume\\.init\\.wave_amplitude cannot be used with flume\\.init\\.dam_position"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1 flume.init.dam_position=4
         flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_wave_direction_z
    EXPECT "flume\\.init\\.wave_direction must be x or y. it is z"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1 flume.init.wave_direction=z)
flume_add_abort_test(abort_wave_direction_two
    EXPECT "flume\\.init\\.wave_direction needs exactly 1 value. it has 2"
    ARGS ${wave_ok} flume.init.wave_amplitude=0.1 flume.init.wave_direction=x\ y)
flume_add_abort_test(abort_wave_direction_without_amplitude
    EXPECT "flume\\.init\\.wave_direction needs flume\\.init\\.wave_amplitude"
    ARGS ${wave_ok} flume.init.wave_direction=y)
flume_add_abort_test(abort_wave_direction_with_amplitude_zero
    EXPECT "flume\\.init\\.wave_direction needs flume\\.init\\.wave_amplitude"
    ARGS ${wave_ok} flume.init.wave_amplitude=0 flume.init.wave_direction=y)
# An amplitude near the largest number: gravity times the depth overflows in
# the velocity of the wave, and the discharge is not finite.
flume_add_abort_test(abort_wave_overflows
    EXPECT "The initial depth or discharge is not finite in [0-9]+ cell\\(s\\): check .*flume\\.init\\.wave_amplitude"
    ARGS ${wave_ok} flume.init.wave_amplitude=${near_max})

# --- Boundaries, rough bed and dam -----------------------------------------------
flume_add_abort_test(abort_periodic_one     EXPECT "flume\\.periodic needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${ok} flume.periodic=1)
flume_add_abort_test(abort_periodic_two     EXPECT "flume\\.periodic must be >= 0 and <= 1 in x and y. it is 2"
    ARGS ${ok} flume.periodic=0\ 2)
flume_add_abort_test(abort_periodic_negative EXPECT "flume\\.periodic must be >= 0 and <= 1 in x and y. it is -1"
    ARGS ${ok} flume.periodic=-1\ 0)
flume_add_abort_test(abort_bed_noise_negative EXPECT "flume\\.bed\\.noise must be >= 0. it is -0\\.5"
    ARGS ${ok} flume.bed.noise=-0.5)
flume_add_abort_test(abort_bed_noise_nan    EXPECT "flume\\.bed\\.noise must be finite"
    ARGS ${ok} flume.bed.noise=nan)
flume_add_abort_test(abort_bed_noise_two    EXPECT "flume\\.bed\\.noise needs exactly 1 value. it has 2"
    ARGS ${ok} flume.bed.noise=0.1\ 0.2)
# A bed and a roughness that are each finite and whose sum is not, in the
# cells that the roughness raises. ${near_max} is 1.7e308 or 3e38.
flume_add_abort_test(abort_bed_noise_overflows
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation, flume\\.bed\\.slope, flume\\.bed\\.bowl_depth, flume\\.bed\\.bowl_radius and flume\\.bed\\.noise"
    ARGS ${ok} flume.bed.elevation=${near_max} flume.bed.noise=${near_max})
set(dam ${ok} flume.init.surface_elevation=2 flume.init.dam_position=4 flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_dam_without_surface
    EXPECT "flume\\.init\\.dam_position needs flume\\.init\\.surface_elevation"
    ARGS ${ok} flume.init.dam_position=4 flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_dam_without_dam_surface
    EXPECT "flume\\.init\\.dam_position needs flume\\.init\\.dam_surface_elevation"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_position=4)
flume_add_abort_test(abort_dam_surface_without_dam
    EXPECT "flume\\.init\\.dam_surface_elevation needs flume\\.init\\.dam_position"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_dam_direction_without_dam
    EXPECT "flume\\.init\\.dam_direction needs flume\\.init\\.dam_position"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_direction=y)
flume_add_abort_test(abort_dam_velocity_without_dam
    EXPECT "flume\\.init\\.dam_velocity needs flume\\.init\\.dam_position"
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_velocity=1\ 0)
flume_add_abort_test(abort_dam_velocity_one
    EXPECT "flume\\.init\\.dam_velocity needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${dam} flume.init.dam_velocity=1)
flume_add_abort_test(abort_dam_velocity_nan
    EXPECT "flume\\.init\\.dam_velocity must be finite in x and y"
    ARGS ${dam} flume.init.dam_velocity=0\ nan)
flume_add_abort_test(abort_dam_direction_unknown
    EXPECT "flume\\.init\\.dam_direction must be x or y. it is z"
    ARGS ${dam} flume.init.dam_direction=z)
flume_add_abort_test(abort_dam_direction_two
    EXPECT "flume\\.init\\.dam_direction needs exactly 1 value. it has 2"
    ARGS ${dam} flume.init.dam_direction=x\ y)
flume_add_abort_test(abort_dam_position_nan
    EXPECT "flume\\.init\\.dam_position must be finite"
    ARGS ${dam} flume.init.dam_position=nan)
flume_add_abort_test(abort_dam_surface_inf
    EXPECT "flume\\.init\\.dam_surface_elevation must be finite"
    ARGS ${dam} flume.init.dam_surface_elevation=inf)
# The dam must leave a cell centre on each side, along its own direction. With
# 8 cells of 1 m the centres are at 0.5 to 7.5 m: a dam at 0.5 m has no cell
# below it (a centre at the dam is on the high side), one at 7.75 m none
# above, and 7.5 m is the last place allowed (abort tests cannot show that it
# runs; dam_initial_on_centre below puts a dam on a centre). 20 m is inside a
# domain of 8 m by 32 m only along y.
set(dam_range "flume\\.init\\.dam_position must leave a cell centre on each side of the dam: it must be > 0\\.5 and <= 7\\.5 in")
flume_add_abort_test(abort_dam_position_at_corner
    EXPECT "${dam_range} x .*. it is 8"
    ARGS ${dam} flume.init.dam_position=8)
flume_add_abort_test(abort_dam_position_first_centre
    EXPECT "${dam_range} x .*. it is 0\\.5"
    ARGS ${dam} flume.init.dam_position=0.5)
flume_add_abort_test(abort_dam_position_past_last_centre
    EXPECT "${dam_range} x .*. it is 7\\.75"
    ARGS ${dam} flume.init.dam_position=7.75)
flume_add_abort_test(abort_dam_position_below
    EXPECT "${dam_range} y .*. it is -1"
    ARGS ${dam} flume.init.dam_position=-1 flume.init.dam_direction=y)
# A misspelled key is reported as unknown, not as the key another one needs.
flume_add_abort_test(abort_dam_key_misspelled
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.init\\.dam_positon "
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_positon=4 flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_dam_position_other_direction
    EXPECT "${dam_range} x .*. it is 20"
    ARGS flume.n_cell=8\ 8 flume.prob_hi=8\ 32 flume.init.surface_elevation=2
         flume.init.dam_position=20 flume.init.dam_surface_elevation=1)

# --- Friction and continued slope ---------------------------------------------------
flume_add_abort_test(abort_manning_negative
    EXPECT "flume\\.friction\\.manning must be >= 0. it is -0\\.5"
    ARGS ${ok} flume.friction.manning=-0.5)
flume_add_abort_test(abort_manning_nan      EXPECT "flume\\.friction\\.manning must be finite"
    ARGS ${ok} flume.friction.manning=nan)
flume_add_abort_test(abort_manning_inf      EXPECT "flume\\.friction\\.manning must be finite"
    ARGS ${ok} flume.friction.manning=inf)
flume_add_abort_test(abort_manning_two
    EXPECT "flume\\.friction\\.manning needs exactly 1 value. it has 2"
    ARGS ${ok} flume.friction.manning=0.03\ 0.04)
flume_add_abort_test(abort_continue_slope_above_one
    EXPECT "flume\\.bed\\.continue_slope must be >= 0 and <= 1. it is 2"
    ARGS ${ok} flume.bed.continue_slope=2)
flume_add_abort_test(abort_continue_slope_minus_one
    EXPECT "flume\\.bed\\.continue_slope must be >= 0 and <= 1. it is -1"
    ARGS ${ok} flume.bed.continue_slope=-1)
flume_add_abort_test(abort_continue_slope_two_values
    EXPECT "flume\\.bed\\.continue_slope needs exactly 1 value. it has 2"
    ARGS ${ok} flume.bed.continue_slope=1\ 1)
set(continue_needs "flume\\.bed\\.continue_slope = 1 needs flume\\.periodic = 1 in a direction where flume\\.bed\\.slope is not 0")
# Walls on all sides; periodic in x with the slope along y; a slope along a
# periodic x that is zero.
flume_add_abort_test(abort_continue_slope_walls
    EXPECT "${continue_needs}"
    ARGS ${ok} flume.bed.slope=0.01\ 0.01 flume.bed.continue_slope=1)
flume_add_abort_test(abort_continue_slope_other_direction
    EXPECT "${continue_needs}"
    ARGS ${ok} flume.periodic=1\ 0 flume.bed.slope=0\ 0.01 flume.bed.continue_slope=1)
flume_add_abort_test(abort_continue_slope_flat
    EXPECT "${continue_needs}"
    ARGS ${ok} flume.periodic=1\ 1 flume.bed.continue_slope=1)

# --- Boundary segments ----------------------------------------------------------
set(seg flume.bc.segments=a)
flume_add_abort_test(abort_segment_name
    EXPECT "flume\\.bc\\.segments: a name may hold only letters, digits and _. it is .a-b."
    ARGS ${ok} flume.bc.segments=a-b)
flume_add_abort_test(abort_segment_twice   EXPECT "flume\\.bc\\.segments lists a twice"
    ARGS ${ok} flume.bc.segments=a\ a flume.bc.a.side=xlo flume.bc.a.type=wall)
flume_add_abort_test(abort_segment_side_missing
    EXPECT "flume\\.bc\\.a\\.side is required: xlo, xhi, ylo or yhi"
    ARGS ${ok} ${seg} flume.bc.a.type=wall)
# A misspelled side is reported as an unknown key, not as a missing side.
flume_add_abort_test(abort_segment_side_misspelled
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.bc\\.a\\.sid "
    ARGS ${ok} ${seg} flume.bc.a.sid=xlo flume.bc.a.type=wall flume.bc.a.range=1\ 2)
flume_add_abort_test(abort_segment_side_unknown
    EXPECT "flume\\.bc\\.a\\.side must be xlo, xhi, ylo or yhi. it is top"
    ARGS ${ok} ${seg} flume.bc.a.side=top flume.bc.a.type=wall)
flume_add_abort_test(abort_segment_side_two
    EXPECT "flume\\.bc\\.a\\.side needs exactly 1 value. it has 2"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo\ xhi flume.bc.a.type=wall)
flume_add_abort_test(abort_segment_type_missing
    EXPECT "flume\\.bc\\.a\\.type is required: wall or outflow"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo)
flume_add_abort_test(abort_segment_type_unknown
    EXPECT "flume\\.bc\\.a\\.type must be wall or outflow. it is inflow"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo flume.bc.a.type=inflow)
# The names of the two types that later versions add.
flume_add_abort_test(abort_segment_type_discharge
    EXPECT "flume\\.bc\\.a\\.type = discharge is reserved for a later version \\(a prescribed inflow discharge from a time series\\)"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo flume.bc.a.type=discharge)
flume_add_abort_test(abort_segment_type_stage
    EXPECT "flume\\.bc\\.a\\.type = stage is reserved for a later version \\(a prescribed water level from a time series\\)"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo flume.bc.a.type=stage)
flume_add_abort_test(abort_segment_periodic_side
    EXPECT "flume\\.bc\\.a\\.side is yhi, but flume\\.periodic is 1 in y: a periodic side has no boundary to set"
    ARGS ${ok} flume.periodic=0\ 1 ${seg} flume.bc.a.side=yhi flume.bc.a.type=outflow)
flume_add_abort_test(abort_segment_range_one
    EXPECT "flume\\.bc\\.a\\.range needs exactly 2 values \\(from and to, in y\\). it has 1"
    ARGS ${ok} ${seg} flume.bc.a.side=xlo flume.bc.a.type=outflow flume.bc.a.range=4)
set(range_bad "flume\\.bc\\.a\\.range must be two coordinates along the side, in x, the first below the second, both from 0 to 8 \\(flume\\.prob_lo and flume\\.prob_hi\\). it is")
flume_add_abort_test(abort_segment_range_outside  EXPECT "${range_bad} -1 4"
    ARGS ${ok} ${seg} flume.bc.a.side=ylo flume.bc.a.type=outflow flume.bc.a.range=-1\ 4)
flume_add_abort_test(abort_segment_range_beyond   EXPECT "${range_bad} 4 8\\.5"
    ARGS ${ok} ${seg} flume.bc.a.side=yhi flume.bc.a.type=outflow flume.bc.a.range=4\ 8.5)
flume_add_abort_test(abort_segment_range_reversed EXPECT "${range_bad} 6 2"
    ARGS ${ok} ${seg} flume.bc.a.side=ylo flume.bc.a.type=outflow flume.bc.a.range=6\ 2)
flume_add_abort_test(abort_segment_range_equal    EXPECT "${range_bad} 3 3"
    ARGS ${ok} ${seg} flume.bc.a.side=ylo flume.bc.a.type=outflow flume.bc.a.range=3\ 3)
flume_add_abort_test(abort_segment_range_nan
    EXPECT "flume\\.bc\\.a\\.range must be finite. it is nan 3"
    ARGS ${ok} ${seg} flume.bc.a.side=ylo flume.bc.a.type=outflow flume.bc.a.range=nan\ 3)
# Cells of 1 m: centres at 0.5, 1.5, ...; none in [2.75, 3.25).
flume_add_abort_test(abort_segment_no_cell
    EXPECT "flume\\.bc\\.a\\.range holds no cell centre: a segment takes the cells whose centre is at the first coordinate or above and below the second. it is 2\\.75 3\\.25"
    ARGS ${ok} ${seg} flume.bc.a.side=ylo flume.bc.a.type=outflow flume.bc.a.range=2.75\ 3.25)
flume_add_abort_test(abort_segment_overlap
    EXPECT "flume\\.bc\\.a and flume\\.bc\\.b overlap on side xlo: a boundary cell can belong to one segment only"
    ARGS ${ok} flume.bc.segments=a\ b flume.bc.a.side=xlo flume.bc.a.type=outflow flume.bc.a.range=0\ 4
         flume.bc.b.side=xlo flume.bc.b.type=wall flume.bc.b.range=3.5\ 8)
# The keys of a segment that flume.bc.segments does not list are unknown.
flume_add_abort_test(abort_segment_not_listed
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.bc\\.b\\.side "
    ARGS ${ok} ${seg} flume.bc.a.side=xlo flume.bc.a.type=wall flume.bc.b.side=xhi)
flume_add_abort_test(abort_segment_no_list
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.bc\\.a\\.side "
    ARGS ${ok} flume.bc.a.side=xlo)

# --- Raster files ---------------------------------------------------------------
# The bad files are in Tests/Raster/bad, each on the domain [0, 40]^2 of
# 8 x 8 cells of 5 m. nodata.asc has 2 x 2 cells of 20 m, the south-west one
# NODATA: bilinear interpolation gives it a weight in every cell whose centre
# is below 30 m in x and in y (6 x 6 cells); the Manning map takes the cell
# holding the centre (4 x 4 cells). ParmParse takes an empty value only in
# quotes.
set(rok flume.n_cell=8\ 8 flume.prob_hi=40\ 40)
set(bad ${CMAKE_CURRENT_SOURCE_DIR}/Raster/bad)
set(rfile "flume\\.bed\\.file \\([^)]*\\): ")
foreach(case
        "unknown_key|unknown header keyword cellsise"
        "twice|header keyword ncols is given twice"
        "no_ncols|the header has no ncols"
        "ncols_fraction|ncols must be a whole number from 1 to 2\\^30; it is 2\\.5"
        "corner_and_center|the header needs exactly one of xllcorner and xllcenter"
        "cellsize_and_dx|the header needs cellsize, or dx and dy, but not both"
        "cellsize_zero|the cell size \\(cellsize, or dx and dy\\) must be > 0"
        "not_a_number|value 4 is not a finite number: \"x\""
        "too_few|the file holds 3 values; the header says"
        "too_many|the file holds 5 values; the header says"
        "small|the raster covers 0 to 38 in x, which does not hold the domain, 0 to 40"
        "short_y|the raster covers 0 to 38 in y, which does not hold the domain, 0 to 40"
        "header_not_number|header keyword cellsize needs a finite number; it has \"twenty\""
        "dx_no_dy|the header needs cellsize, or dx and dy, but not both"
        "y_corner_and_center|the header needs exactly one of yllcorner and yllcenter"
        "no_nrows|the header has no nrows"
        "extent_inf|the raster's extent is not finite"
        "huge_header|the header says 1073741824 x 1073741824 = 1152921504606846976 values, more than a file of 78 bytes can hold"
        "nodata|36 cell\\(s\\) of the domain read a NODATA_value cell of the raster")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 file)
    list(GET pair 1 expect)
    flume_add_abort_test(abort_raster_${file} EXPECT "${rfile}${expect}"
        ARGS ${rok} flume.bed.file=${bad}/${file}.asc)
endforeach()
flume_add_abort_test(abort_raster_missing
    EXPECT "flume\\.bed\\.file \\([^)]*no_such_file\\.asc\\): the file cannot be opened"
    ARGS ${rok} flume.bed.file=${bad}/no_such_file.asc)
flume_add_abort_test(abort_raster_directory
    EXPECT "flume\\.bed\\.file \\([^)]*bad\\): the path is not a regular file"
    ARGS ${rok} flume.bed.file=${bad})
# A value beyond the range of single precision (1e39) is read as it is in a
# double build; in a single build the run aborts. Columns and rows count
# from 0.
if(FLUME_PRECISION STREQUAL "SINGLE")
    flume_add_abort_test(abort_raster_beyond_float
        EXPECT "${rfile}value 3 \\(column 0 from the west, row 0 from the south\\) is (1e\\+39|9\\.99+[0-9]*e\\+38), beyond the range of the build's precision"
        ARGS ${rok} flume.bed.file=${bad}/beyond_float.asc)
endif()
foreach(key buildings.file friction.manning_file)
    string(REPLACE "." "\\." key_re "${key}")
    string(REPLACE "." "_" key_name "${key}")
    flume_add_abort_test(abort_raster_${key_name}_empty
        EXPECT "flume\\.${key_re} must not be empty: the name of an ESRI ASCII grid file"
        ARGS ${rok} "flume.${key}=\"\"")
endforeach()
# Each key of the analytic bed with a bed file. continue_slope = 0 is given
# but does nothing, and still aborts. bowl_radius alone aborts first because
# it needs bowl_depth; with bowl_depth the message names bowl_depth.
foreach(case "elevation|1" "continue_slope|0" "bowl_depth|1 flume.bed.bowl_radius=10")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 name)
    list(GET pair 1 value)
    separate_arguments(value)
    flume_add_abort_test(abort_raster_bed_file_with_${name}
        EXPECT "flume\\.bed\\.file cannot be used with flume\\.bed\\.${name}: the file gives the bed of every cell"
        ARGS ${rok} flume.bed.file=${raster}/streets_bed.asc flume.bed.${name}=${value})
endforeach()
flume_add_abort_test(abort_raster_buildings_nodata
    EXPECT "flume\\.buildings\\.file \\([^)]*\\): 16 cell\\(s\\) of the domain read a NODATA_value cell"
    ARGS ${rok} flume.buildings.file=${bad}/nodata.asc)
flume_add_abort_test(abort_raster_bed_file_empty
    EXPECT "flume\\.bed\\.file must not be empty: the name of an ESRI ASCII grid file"
    ARGS ${rok} "flume.bed.file=\"\"")
flume_add_abort_test(abort_raster_bed_file_with_slope
    EXPECT "flume\\.bed\\.file cannot be used with flume\\.bed\\.slope: the file gives the bed of every cell"
    ARGS ${rok} flume.bed.file=${raster}/streets_bed.asc flume.bed.slope=0.001\ 0)
flume_add_abort_test(abort_raster_bed_file_with_noise
    EXPECT "flume\\.bed\\.file cannot be used with flume\\.bed\\.noise"
    ARGS ${rok} flume.bed.file=${raster}/streets_bed.asc flume.bed.noise=0.01)
flume_add_abort_test(abort_raster_manning_file_with_manning
    EXPECT "flume\\.friction\\.manning_file cannot be used with flume\\.friction\\.manning"
    ARGS ${rok} flume.friction.manning_file=${raster}/manning_uniform.asc flume.friction.manning=0.03)
flume_add_abort_test(abort_raster_manning_negative
    EXPECT "flume\\.friction\\.manning_file \\([^)]*\\): 16 cell\\(s\\) of the domain take a Manning coefficient below 0"
    ARGS ${rok} flume.friction.manning_file=${bad}/negative.asc)
flume_add_abort_test(abort_raster_manning_nodata
    EXPECT "flume\\.friction\\.manning_file \\([^)]*\\): 16 cell\\(s\\) of the domain read a NODATA_value cell"
    ARGS ${rok} flume.friction.manning_file=${bad}/nodata.asc)
flume_add_abort_test(abort_raster_buildings_height_no_file
    EXPECT "flume\\.buildings\\.height needs flume\\.buildings\\.file"
    ARGS ${rok} flume.buildings.height=5)
flume_add_abort_test(abort_raster_buildings_height_zero
    EXPECT "flume\\.buildings\\.height must be > 0; it is 0"
    ARGS ${rok} flume.buildings.file=${raster}/streets_mask.asc flume.buildings.height=0)
flume_add_abort_test(abort_raster_buildings_negative
    EXPECT "flume\\.buildings\\.file \\([^)]*\\): 16 cell\\(s\\) of the domain take a building height below 0"
    ARGS ${rok} flume.buildings.file=${bad}/negative.asc)


# --- Prescribed rain ---------------------------------------------------------------
# In their own file, which uses the functions and ${ok} above.
include(${CMAKE_CURRENT_SOURCE_DIR}/Rain/RainTests.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/RainField/RainFieldTests.cmake)
