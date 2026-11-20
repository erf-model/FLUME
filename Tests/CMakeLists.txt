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
flume_add_regtest(dry_stop_time_huge NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.stop_time=1e308
    CHECK ${check_plotfile} skel --steps 0 1 --times 0 1e308 ${skel_grid} --n-files 2)

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
# (bound 1 %); beyond the exact front the depth stays below 1e-11 of 2 m (the
# run gives 6.5e-14 of 2 m).
set(ritter_bounds --l1-max 0.01 --front-max 1e-11)
flume_add_regtest(dam_break_dry_x NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.init.dam_surface_elevation=0
    CHECK ${check_dam_break} dam --time 5 --direction x --position 50
          --surface-lo 2 --surface-hi 0 ${ritter_bounds})
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
# dry cell ends 1/3 m deep and the last wet cell 5/3 m deep. Later steps hardly use these speeds (the bed ahead of the
# front is then thinly wet), so the 5 s runs do not test them.
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
    set(near_max 3e38)
else()
    set(one_step 1.0000000000000002)
    set(too_big 1e308)
    set(tiny_extent 4e-300)
    set(huge_extent 4e155)
    set(tiny_depth 1e-320)
    set(near_max 1.7e308)
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
# A step too long for the flow. Water 2 m deep in the one column of cells
# below the dam, with a dry bed beyond it and across the periodic side: the
# column loses water through both faces and ends near h (1 - 4 cfl / 3),
# which is negative for flume.cfl = 0.8 (Docs/sphinx_doc/Model.rst).
flume_add_abort_test(abort_depth_negative
    EXPECT "The depth is negative or not finite, or the wave speed is not finite, in 4 cell\\(s\\) at step 1: .*lower flume\\.cfl"
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs flume.periodic=1\ 0
         flume.prob_hi=100\ 1000 flume.init.dam_position=0.3
         flume.init.dam_surface_elevation=0 flume.cfl=0.8 flume.max_step=3)
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

# --- Boundaries, rough bed and dam -----------------------------------------------
flume_add_abort_test(abort_periodic_one     EXPECT "flume\\.periodic needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${ok} flume.periodic=1)
flume_add_abort_test(abort_periodic_two     EXPECT "flume\\.periodic must be >= 0 and <= 1 in x and y; it is 2"
    ARGS ${ok} flume.periodic=0\ 2)
flume_add_abort_test(abort_periodic_negative EXPECT "flume\\.periodic must be >= 0 and <= 1 in x and y; it is -1"
    ARGS ${ok} flume.periodic=-1\ 0)
flume_add_abort_test(abort_bed_noise_negative EXPECT "flume\\.bed\\.noise must be >= 0; it is -0\\.5"
    ARGS ${ok} flume.bed.noise=-0.5)
flume_add_abort_test(abort_bed_noise_nan    EXPECT "flume\\.bed\\.noise must be finite"
    ARGS ${ok} flume.bed.noise=nan)
flume_add_abort_test(abort_bed_noise_two    EXPECT "flume\\.bed\\.noise needs exactly 1 value; it has 2"
    ARGS ${ok} flume.bed.noise=0.1\ 0.2)
# A bed and a roughness that are each finite and whose sum is not, in the
# cells that the roughness raises. ${near_max} is 1.7e308 or 3e38.
flume_add_abort_test(abort_bed_noise_overflows
    EXPECT "The bed elevation is not finite in [0-9]+ cell\\(s\\): check flume\\.bed\\.elevation, flume\\.bed\\.slope and flume\\.bed\\.noise"
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
    EXPECT "flume\\.init\\.dam_velocity needs exactly 2 values \\(x and y\\); it has 1"
    ARGS ${dam} flume.init.dam_velocity=1)
flume_add_abort_test(abort_dam_velocity_nan
    EXPECT "flume\\.init\\.dam_velocity must be finite in x and y"
    ARGS ${dam} flume.init.dam_velocity=0\ nan)
flume_add_abort_test(abort_dam_direction_unknown
    EXPECT "flume\\.init\\.dam_direction must be x or y; it is z"
    ARGS ${dam} flume.init.dam_direction=z)
flume_add_abort_test(abort_dam_direction_two
    EXPECT "flume\\.init\\.dam_direction needs exactly 1 value; it has 2"
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
    EXPECT "${dam_range} x .*; it is 8"
    ARGS ${dam} flume.init.dam_position=8)
flume_add_abort_test(abort_dam_position_first_centre
    EXPECT "${dam_range} x .*; it is 0\\.5"
    ARGS ${dam} flume.init.dam_position=0.5)
flume_add_abort_test(abort_dam_position_past_last_centre
    EXPECT "${dam_range} x .*; it is 7\\.75"
    ARGS ${dam} flume.init.dam_position=7.75)
flume_add_abort_test(abort_dam_position_below
    EXPECT "${dam_range} y .*; it is -1"
    ARGS ${dam} flume.init.dam_position=-1 flume.init.dam_direction=y)
# A misspelled key is reported as unknown, not as the key another one needs.
flume_add_abort_test(abort_dam_key_misspelled
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.init\\.dam_positon "
    ARGS ${ok} flume.init.surface_elevation=2 flume.init.dam_positon=4 flume.init.dam_surface_elevation=1)
flume_add_abort_test(abort_dam_position_other_direction
    EXPECT "${dam_range} x .*; it is 20"
    ARGS flume.n_cell=8\ 8 flume.prob_hi=8\ 32 flume.init.surface_elevation=2
         flume.init.dam_position=20 flume.init.dam_surface_elevation=1)
