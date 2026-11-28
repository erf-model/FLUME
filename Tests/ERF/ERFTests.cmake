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

# --- ERF's step and FLUME's substeps ------------------------------------------
# Mostly on 4 ranks: ERF's 32 x 4 columns in boxes of 8 x 4 (flood boxes of
# 16 x 8), the new stand-alone runs in boxes of 16. max_dt_plane and the runs
# above use 2 ranks and boxes of 32. The plotfiles of runs on these different
# boxes are compared cell by cell to the last bit here, so FLUME's cells did
# not depend on the boxes in these runs. The sums in water_balance.txt do
# depend on them, in the last digits, so logs are compared only between runs
# on the same boxes.
flume_launcher(launcher_np4 4)
set(erf_np4 "amr.max_grid_size=8 4 8")
set(erf_uncapped ${erf_flume_plane} flume.max_dt=1000 erf.check_int=-1)

# The stand-alone runs of check_step_lengths.py: every step 0.5 s and 0.125 s
# (with max_dt_plane, 0.25 s) for the time-converged flood and the order of
# the time error, and FLUME's own steps.
flume_add_regtest(flood_plane_dt_05 NPROCS 4
    INPUTS ${erf_case_dir}/inputs_standalone
    ARGS ${erf_flume_plane} flume.max_dt=0.5 flume.max_grid_size=16 flume.plot_int=100
    CHECK ${check_coupled} flood --dt 0.5)
flume_add_regtest(flood_plane_dt_0125 NPROCS 4
    INPUTS ${erf_case_dir}/inputs_standalone
    ARGS ${erf_flume_plane} flume.max_dt=0.125 flume.max_grid_size=16 flume.plot_int=400
    CHECK ${check_coupled} flood --dt 0.125)
flume_add_regtest(flood_plane_own_steps NPROCS 4
    INPUTS ${erf_case_dir}/inputs_standalone
    ARGS ${erf_flume_plane} flume.max_dt=1000 flume.max_grid_size=16
    OUTPUT "FLUME: 99 step\\(s\\), time 150 s"
    CHECK ${CMAKE_COMMAND} -E true)

# erf_flume_plane_short_steps: no cap on FLUME's step, ERF's step 0.25 s, four
# times shorter than FLUME's own (about 1 s): each ERF step is one FLUME step
# that ends on it, and the flood is that of max_dt_plane (every step capped at
# 0.25 s) to the last bit: plotfiles, water balance and gauges. On the boxes
# and ranks of max_dt_plane (2), so that the logs can be compared byte for byte.
erf_flume_add_regtest(erf_flume_plane_short_steps NPROCS 2
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_uncapped} erf.fixed_dt=0.25
    CHECK ${check_coupled} flood --dt 0.25
          --reference ${FLUME_TEST_RUN_DIR}/max_dt_plane --log ${erf_logs})
if(TEST erf_flume_plane_short_steps)
    set_tests_properties(erf_flume_plane_short_steps PROPERTIES FIXTURES_REQUIRED max_dt_plane)
endif()

# erf_flume_plane_step_<dt>: ERF's steps of 1.5 s (about 1.5 FLUME steps),
# 2 s (erf_flume_plane_substeps on 4 ranks), 3.75 s (not a multiple of
# FLUME's) and 15 s (about 10 FLUME steps) to 150 s, FLUME landing on the end
# of each: the runs of erf_flume_plane_step_lengths.
set(step_length_runs)
set(step_length_fixtures max_dt_plane flood_plane_dt_05 flood_plane_dt_0125 flood_plane_own_steps)
foreach(case "1.5|100|50" "2.0|75|25" "3.75|40|20" "15|10|5")
    string(REPLACE "|" ";" case "${case}")
    list(GET case 0 dt)
    list(GET case 1 nsteps)
    list(GET case 2 plot)
    erf_flume_add_regtest(erf_flume_plane_step_${dt} NPROCS 4
        INPUTS ${erf_case_dir}/inputs_erf
        ARGS ${erf_uncapped} ${erf_np4} erf.fixed_dt=${dt} max_step=${nsteps}
             erf.plot_int_1=${nsteps} erf.flood.plot_int=${plot}
        CHECK ${check_coupled} flood --dt ${dt} --substeps)
    list(APPEND step_length_runs ${FLUME_TEST_RUN_DIR}/erf_flume_plane_step_${dt})
    list(APPEND step_length_fixtures erf_flume_plane_step_${dt})
endforeach()

foreach(t IN LISTS step_length_fixtures)
    if(TEST ${t})
        set_tests_properties(${t} PROPERTIES FIXTURES_SETUP ${t})
    endif()
endforeach()

# erf_flume_plane_step_lengths: at 150 s, no run inside ERF is further from the
# time-converged flood than the stand-alone run with FLUME's own steps, and
# each is within half of that run's time error of it (check_step_lengths.py;
# the first-order time error is checked first, from three stand-alone runs).
if(FLUME_MPI)
    add_test(NAME erf_flume_plane_step_lengths
        COMMAND ${Python3_EXECUTABLE} ${erf_case_dir}/check_step_lengths.py flood
                --coarse ${FLUME_TEST_RUN_DIR}/flood_plane_dt_05
                --fine ${FLUME_TEST_RUN_DIR}/max_dt_plane
                --finer ${FLUME_TEST_RUN_DIR}/flood_plane_dt_0125
                --own ${FLUME_TEST_RUN_DIR}/flood_plane_own_steps
                --runs ${step_length_runs})
    set_tests_properties(erf_flume_plane_step_lengths PROPERTIES LABELS "erf;regression"
                         FIXTURES_REQUIRED "${step_length_fixtures}" TIMEOUT 120)
endif()

# erf_flume_plane_substeps_restart: erf_flume_plane_step_2.0 (no checkpoint)
# with a checkpoint at ERF step 37, after 51 FLUME substeps. The flood
# plotfiles and logs are those of the run without a checkpoint, to the
# last bit; then ERF continues from chk00037 and gives them again.
set(substeps_args flume.max_dt=1000 erf.fixed_dt=2.0 max_step=75
                  erf.plot_int_1=75 erf.flood.plot_int=25 ${erf_np4})
erf_flume_add_regtest(erf_flume_plane_substeps_restart NPROCS 4
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_flume_plane} ${substeps_args} erf.check_int=37
    CHECK ${check_coupled} flood --dt 2.0 --substeps
          --reference ${FLUME_TEST_RUN_DIR}/erf_flume_plane_step_2.0 --log ${erf_logs}
          --restart chk00037 ${launcher_np4} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                    ${erf_flume_plane} ${substeps_args} erf.check_int=-1 erf.restart=chk00037
                    ${FLUME_TEST_ARGS})
if(TEST erf_flume_plane_substeps_restart)
    set_tests_properties(erf_flume_plane_substeps_restart PROPERTIES
                         FIXTURES_REQUIRED erf_flume_plane_step_2.0)
endif()

# erf_flume_plane_ranks_restart: erf_flume_plane_substeps_restart, but ERF
# continues from chk00037 on 3 ranks instead of 4: the same 4 boxes, one rank
# with two of them. The flood plotfiles and gauge files are those of
# erf_flume_plane_step_2.0 to the last bit. The sums of water_balance.txt are
# taken over the ranks in another order after the restart, so they agree to
# 1e-14 of each column's largest value, the residual to 1e-14 of the largest
# volume (2.3e-16 in this run), not byte for byte. The restart's output must
# report 3 MPI processes.
flume_launcher(launcher_np3 3)
erf_flume_add_regtest(erf_flume_plane_ranks_restart NPROCS 4
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${erf_flume_plane} ${substeps_args} erf.check_int=37
    CHECK ${check_coupled} flood --dt 2.0 --substeps
          --reference ${FLUME_TEST_RUN_DIR}/erf_flume_plane_step_2.0
          --restart-balance-rtol 1e-14 --restart-ranks 3 --log ${erf_logs}
          --restart chk00037 ${launcher_np3} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                    ${erf_flume_plane} ${substeps_args} erf.check_int=-1 erf.restart=chk00037
                    ${FLUME_TEST_ARGS})
if(TEST erf_flume_plane_ranks_restart)
    set_tests_properties(erf_flume_plane_ranks_restart PROPERTIES
                         FIXTURES_REQUIRED erf_flume_plane_step_2.0)
endif()

# erf_flume_start_without_flood: ERF runs 37 steps of 2 s without the flood
# (erf.flood.enable = false) and writes chk00037 at 74 s, with no FLUME state
# in it. It then continues from chk00037 with the flood on, to max_step = 75
# (150 s), FLUME starting from a flat water surface at 0.5 m. The bed falls
# from 1 m by 1 cm per metre, so the cells with centres beyond x = 50 m, 14
# columns of 8, are wet at the start: 112 cells, 0.005 to 0.135 m deep.
# FLUME's summary must be written at ERF's step 75; its logs and gauges
# start at 74 s; the arrival time of those 112 cells is 74 s and that of no
# other cell is earlier (check_start_without_flood.py).
set(start_args ${erf_flume_plane} ${substeps_args} erf.check_int=37)
erf_flume_add_regtest(erf_flume_start_without_flood NPROCS 4
    INPUTS ${erf_case_dir}/inputs_erf
    ARGS ${start_args} erf.flood.enable=false max_step=37
    CHECK ${Python3_EXECUTABLE} ${erf_case_dir}/check_start_without_flood.py flood
          --start-time 74 --end-time 150 --surface 0.5 --wet 112
          --restart chk00037 ${launcher_np4} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                    ${start_args} erf.check_int=-1 erf.restart=chk00037
                    flume.flood_maps=1 flume.init.surface_elevation=0.5 ${FLUME_TEST_ARGS})
if(TEST erf_flume_start_without_flood)
    set_tests_properties(erf_flume_start_without_flood PROPERTIES
                         FIXTURES_SETUP erf_flume_start_without_flood)
endif()

# abort_erf_flume_checkpoint_title: FLUME reads the step of a restart from
# the Header of ERF's checkpoint, which ERF writes as "Checkpoint file for
# ERF", the finest level, four lines of component counts and the step of
# each level. ERF itself skips the title line, so a copy of chk00037 of
# erf_flume_start_without_flood with another title restarts in ERF, and
# FLUME must stop on it rather than read the lines that follow.
if(FLUME_MPI)
    set(retitled_chk ${FLUME_TEST_RUN_DIR}/erf_flume_retitled_chk00037)
    add_test(NAME erf_flume_retitle_checkpoint
        COMMAND ${Python3_EXECUTABLE} ${erf_case_dir}/retitle_checkpoint.py
                ${FLUME_TEST_RUN_DIR}/erf_flume_start_without_flood/chk00037 ${retitled_chk}
                "Checkpoint file for something else")
    set_tests_properties(erf_flume_retitle_checkpoint PROPERTIES LABELS "erf;abort" TIMEOUT 60
                         FIXTURES_REQUIRED erf_flume_start_without_flood
                         FIXTURES_SETUP erf_flume_retitled_chk)
    erf_flume_add_abort_test(abort_erf_flume_checkpoint_title
        EXPECT "flume: ERF's checkpoint header .*erf_flume_retitled_chk00037/Header does not start with the line \"Checkpoint file for ERF\": FLUME reads the step of the restart from it"
        ARGS erf.restart=${retitled_chk})
    set_tests_properties(abort_erf_flume_checkpoint_title PROPERTIES
                         FIXTURES_REQUIRED erf_flume_retitled_chk)
endif()

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
