# FLUME as ERF's flood engine (erf.flood.engine = flume), with FLUME_WITH_ERF.
# Included from Tests/CMakeLists.txt after flume_add_regtest and
# flume_add_abort_test are defined. See Tests/ERF/README.md.
#
# The ERF runs go through Tests/run_case.py, as FLUME's own tests do: each in
# its own directory under Tests/run, the water balance of the run checked to
# round-off, then the check command.

set(erf_case_dir ${CMAKE_CURRENT_SOURCE_DIR}/ERF)
set(check_coupled ${Python3_EXECUTABLE} ${erf_case_dir}/check_coupled.py)
# FLUME's half of each deck, read by both the stand-alone run and ERF.
set(erf_flume_plane FILE=${erf_case_dir}/flume_plane)
set(erf_logs water_balance.txt gauge_mid.txt gauge_outlet.txt)
flume_launcher(launcher_np2 2)

# erf_flume_add_regtest(<name> NPROCS <n> INPUTS <file> [ARGS <key=value>...]
#                       [EXPECT_OUTPUT <regular expression>]
#                       CHECK <command and arguments, run in the test directory>)
# flume_add_regtest with erf_exec in place of flume; with EXPECT_OUTPUT, the
# output of the run must match the expression.
function(erf_flume_add_regtest name)
    cmake_parse_arguments(T "" "NPROCS;INPUTS;EXPECT_OUTPUT" "ARGS;CHECK" ${ARGN})
    if(NOT FLUME_MPI AND NOT T_NPROCS EQUAL 1)
        return()
    endif()
    flume_launcher(launcher ${T_NPROCS})
    set(expect)
    if(T_EXPECT_OUTPUT)
        set(expect --expect-output ${T_EXPECT_OUTPUT})
    endif()
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name} ${expect}
                --run ${launcher} $<TARGET_FILE:erf_exec> ${T_INPUTS} ${T_ARGS} ${FLUME_TEST_ARGS}
                --check ${T_CHECK})
    set_tests_properties(${name} PROPERTIES LABELS "erf;regression" PROCESSORS ${T_NPROCS}
                         TIMEOUT 600)
endfunction()

# erf_flume_add_abort_test(<name> EXPECT <regular expression> ARGS <key=value>...)
# erf_exec on the coupled deck, one rank, with the keys given: the run must
# stop and print the expression.
function(erf_flume_add_abort_test name)
    cmake_parse_arguments(T "" "EXPECT" "ARGS" ${ARGN})
    flume_launcher(launcher 1)
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name}
                --expect-abort ${T_EXPECT}
                --run ${launcher} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                      ${erf_flume_plane} max_step=1 ${T_ARGS} ${FLUME_TEST_ARGS}
                      amrex.call_addr2line=0)
    set_tests_properties(${name} PROPERTIES LABELS "erf;abort" TIMEOUT 120)
endfunction()

# --- Rain on a plane inside ERF ------------------------------------------------
# The stand-alone run, max_dt_plane, is registered in Tests/CMakeLists.txt, in
# every build.

# erf_flume_plane: ERF with FLUME as its flood engine on the same flood grid
# (ERF's 32 x 4 columns refined by 2), ERF's fixed step 0.25 s. Every flood
# plotfile, the water balance and the gauge files are those of the stand-alone
# run, to the last bit. The flood fields of ERF's plotfile are the means of
# FLUME's cells over each column. Then ERF continues from its checkpoint of step 300,
# with FLUME's state in it, and gives the same files again.
erf_flume_add_regtest(erf_flume_plane NPROCS 2
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_flume_plane}
    CHECK ${check_coupled} flood --dt 0.25
          --reference ${FLUME_TEST_RUN_DIR}/max_dt_plane --log ${erf_logs}
          --column-mean plt00600 flood00600 2
          --restart chk00300 ${launcher_np2} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                    ${erf_flume_plane} erf.restart=chk00300 ${FLUME_TEST_ARGS})
if(TEST erf_flume_plane)
    set_tests_properties(erf_flume_plane PROPERTIES FIXTURES_REQUIRED max_dt_plane)
endif()

# erf_flume_plane_substeps: the same flood with ERF's step 2 s and no cap on
# FLUME's (flume.max_dt = 1000 s): FLUME takes its CFL steps, about 1 s, and
# lands on the end of every ERF step. Each flood plotfile is at ERF's time
# after more FLUME steps than ERF steps, and the water balance closes to
# round-off (run_case.py). The flood is not that of the stand-alone run.
erf_flume_add_regtest(erf_flume_plane_substeps NPROCS 2
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_flume_plane} flume.max_dt=1000 erf.fixed_dt=2.0 max_step=75
         erf.check_int=-1 erf.plot_int_1=75 erf.flood.plot_int=25
    CHECK ${check_coupled} flood --dt 2.0 --substeps)

# erf_flume_plane_stop_time: ERF ends on its stop time, 52.2 s, after 174
# steps of 0.3 s, which sum to 52.19999999999986 s: short of 52.2 s by
# 1.4e-13 s, three times 4 machine epsilons of the stop time and less than
# 1e-6 of the step (ERF stops within 1e-6 of a step of its stop time).
# FLUME's cap, flume.max_dt = 0.29999999 s, falls short of ERF's step by
# 1e-8 s, under 1e-6 of the step, so FLUME lands on the end of each ERF step
# in one step. FLUME must see that the 174th step is the last one and write
# its summary after 174 steps; the output must say so.
erf_flume_add_regtest(erf_flume_plane_stop_time NPROCS 2
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_flume_plane} flume.max_dt=0.29999999 erf.fixed_dt=0.3 stop_time=52.2
         max_step=1000 erf.check_int=-1 erf.plot_int_1=-1 erf.flood.plot_int=58
    EXPECT_OUTPUT "FLUME: 174 step\\(s\\), time 52\\.2 s"
    CHECK ${CMAKE_COMMAND} -E true)

# --- What FLUME does not take inside ERF ---------------------------------------
# Each key ERF sets stops the run, with a message naming ERF's key.
foreach(case "n_cell=64 8|amr\\.n_cell" "prob_lo=0 0|geometry\\.prob_lo"
             "prob_hi=64 8|geometry\\.prob_hi" "periodic=0 1|geometry\\.is_periodic"
             "max_grid_size=32|amr\\.max_grid_size" "plot_file=f|use erf\\.flood\\.plot_file"
             "plot_int=10|use erf\\.flood\\.plot_int" "check_file=fchk|use erf\\.check_file"
             "check_int=10|use erf\\.check_int" "restart=fchk00010|use erf\\.restart"
             "stop_time=10|ERF's stop_time ends" "max_step=10|ERF's max_step ends")
    string(REGEX REPLACE "\\|.*" "" key_value "${case}")
    string(REGEX REPLACE ".*\\|" "" erf_key "${case}")
    string(REGEX REPLACE "=.*" "" key "${key_value}")
    erf_flume_add_abort_test(abort_erf_flume_host_${key}
        EXPECT "flume\\.${key} cannot be set when FLUME runs inside ERF \\(erf\\.flood\\.engine = flume\\):.*${erf_key}"
        ARGS "flume.${key_value}")
endforeach()
