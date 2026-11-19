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
#                   CHECK <command and arguments, run in the test directory>)
# Without MPI only the NPROCS 1 tests are registered.
function(flume_add_regtest name)
    cmake_parse_arguments(T "" "NPROCS;INPUTS" "ARGS;CHECK" ${ARGN})
    if(NOT T_NPROCS OR NOT T_INPUTS OR NOT T_CHECK)
        message(FATAL_ERROR "flume_add_regtest(${name}): NPROCS, INPUTS and CHECK are required")
    endif()
    if(NOT FLUME_MPI AND NOT T_NPROCS EQUAL 1)
        return()
    endif()
    flume_launcher(launcher ${T_NPROCS})
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name}
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

# --- Skeleton: the grid and the plotfile ------------------------------------
# 48 x 32 cells on [-10, 86] x [5, 69] m, boxes of at most 16 cells: 6 boxes of
# 2 m cells, in a plotfile named skel00000. The check reads the plotfile
# headers and pins these numbers; the depth must be zero in every box. Run on
# 1, 2 and 4 ranks: each rank that owns a box writes one data file, so the
# number of data files shows that the run used the ranks it was given.
set(check_plotfile ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_plotfile.py)
foreach(np 1 2 4)
    flume_add_regtest(skeleton_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
        CHECK ${check_plotfile} skel00000
              --n-cell 48 32 --prob-lo -10 5 --prob-hi 86 69 --dx 2 2
              --max-grid-size 16 --n-boxes 6 --n-files ${np} --var depth --value 0)
endforeach()

# Defaults: prob_lo = 0 0, max_grid_size = 64 (128 x 65 cells: 4 boxes; any
# other default gives another number), plot_file = plt.
flume_add_regtest(skeleton_defaults NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs_defaults
    CHECK ${check_plotfile} plt00000
          --n-cell 128 65 --prob-lo 0 0 --prob-hi 320 162.5 --dx 2.5 2.5
          --max-grid-size 64 --n-boxes 4 --n-files 1 --var depth --value 0)

# --- Inputs that must abort ---------------------------------------------------
# One test per check in Source/FLUME_Inputs.cpp.
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
else()
    set(one_step 1.0000000000000002)
    set(too_big 1e308)
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
