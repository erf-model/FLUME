# ---------------------------------------------------------------------------
# Checkpoint and restart (Tests/Restart/README.md). Included from
# Tests/CMakeLists.txt, which defines the functions and ${ok}.
#
# A reference run goes to the end without a break and writes checkpoints. A
# restart test runs the first part of the same run, which writes a
# checkpoint and may go past it, then continues from the checkpoint
# (Tests/Restart/check_restart.py): every plotfile must be that of the
# reference to the last bit, the gauge files byte for byte, and
# water_balance.txt byte for byte on the same rank count and boxes, or within
# round-off of the sums over the ranks on others.
# ---------------------------------------------------------------------------
set(rs_dir ${CMAKE_CURRENT_SOURCE_DIR}/Restart)
set(rs_main ${rs_dir}/inputs)
set(rs_first ${rs_dir}/inputs_first_order)
set(check_rs ${Python3_EXECUTABLE} ${rs_dir}/check_restart.py)

# flume_add_restart_test(<name> REFERENCE <test> INPUTS <file>
#                        FIRST_NPROCS <n> RESTART_NPROCS <m>
#                        [FIRST_ARGS <key=value>...] [RESTART_ARGS <key=value>...]
#                        [CHECK_ARGS <options of check_restart.py>...])
# The first part runs on FIRST_NPROCS ranks with FIRST_ARGS, the restart on
# RESTART_NPROCS ranks with RESTART_ARGS (which name the checkpoint with
# flume.restart). Registered only when the reference test is.
function(flume_add_restart_test name)
    cmake_parse_arguments(T "" "REFERENCE;INPUTS;FIRST_NPROCS;RESTART_NPROCS"
                          "FIRST_ARGS;RESTART_ARGS;CHECK_ARGS" ${ARGN})
    if(NOT TEST ${T_REFERENCE})
        return()
    endif()
    if(NOT FLUME_MPI AND NOT T_RESTART_NPROCS EQUAL 1)
        return()
    endif()
    flume_launcher(restart_launcher ${T_RESTART_NPROCS})
    flume_add_regtest(${name} NPROCS ${T_FIRST_NPROCS}
        INPUTS ${T_INPUTS}
        ARGS ${T_FIRST_ARGS}
        CHECK ${check_rs} --reference ${FLUME_TEST_RUN_DIR}/${T_REFERENCE} ${T_CHECK_ARGS}
              --restart ${restart_launcher} $<TARGET_FILE:flume> ${T_INPUTS} ${T_RESTART_ARGS}
                        ${FLUME_TEST_ARGS})
    if(TEST ${name})
        set(np ${T_FIRST_NPROCS})
        if(T_RESTART_NPROCS GREATER np)
            set(np ${T_RESTART_NPROCS})
        endif()
        set_tests_properties(${name} PROPERTIES PROCESSORS ${np}
                             FIXTURES_REQUIRED ${T_REFERENCE})
    endif()
endfunction()

# --- References and the independence of the decomposition ------------------
# Tests/Restart/inputs: 856 steps to 240 s, checkpoints at steps 400, 800
# and 856 (the last). The same run on 2 ranks, and on 1 rank in boxes of 8
# cells (24 boxes in place of 6), gives every value of every plotfile to the
# last bit (one reduction, a maximum, sets the step).
flume_add_regtest(restart_reference_np1 NPROCS 1
    INPUTS ${rs_main}
    ARGS flume.check_int=400
    CHECK ${check_rs} --expect chk00400 chk00800 chk00856)
set_tests_properties(restart_reference_np1 PROPERTIES FIXTURES_SETUP restart_reference_np1)
flume_add_regtest(restart_reference_np2 NPROCS 2
    INPUTS ${rs_main}
    ARGS flume.check_int=400
    CHECK ${check_rs} --compare ${FLUME_TEST_RUN_DIR}/restart_reference_np1)
flume_add_regtest(restart_boxes_8_np1 NPROCS 1
    INPUTS ${rs_main}
    ARGS flume.max_grid_size=8
    CHECK ${check_rs} --compare ${FLUME_TEST_RUN_DIR}/restart_reference_np1)
foreach(t restart_reference_np2 restart_boxes_8_np1)
    if(TEST ${t})
        set_tests_properties(${t} PROPERTIES FIXTURES_REQUIRED restart_reference_np1)
    endif()
endforeach()
if(TEST restart_reference_np2)
    set_tests_properties(restart_reference_np2 PROPERTIES FIXTURES_SETUP restart_reference_np2)
endif()

# Tests/Restart/inputs_first_order: 785 steps to 60 s, a checkpoint every 300
# steps. The run without a log, on 1 rank in boxes of 8, gives the plotfiles
# of the 2-rank run; its checkpoint, which keeps no balance, serves an abort
# test below.
flume_add_regtest(restart_first_order_reference NPROCS 2
    INPUTS ${rs_first}
    ARGS flume.check_int=300
    CHECK ${check_rs} --expect chk00300 chk00600 chk00785)
if(TEST restart_first_order_reference)
    set_tests_properties(restart_first_order_reference PROPERTIES
                         FIXTURES_SETUP restart_first_order_reference)
    flume_add_regtest(restart_first_order_no_log NPROCS 1
        INPUTS ${rs_first}
        ARGS flume.max_grid_size=8 flume.log_int=0 flume.check_int=300
        CHECK ${check_rs} --compare ${FLUME_TEST_RUN_DIR}/restart_first_order_reference)
    set_tests_properties(restart_first_order_no_log PROPERTIES
                         FIXTURES_REQUIRED restart_first_order_reference
                         FIXTURES_SETUP restart_first_order_no_log)
endif()

# --- Restarts -----------------------------------------------------------------
# The same rank count and boxes: the first part goes on to step 550, past the
# checkpoint of step 400, so the restart cuts the logs back to step 400. Each
# log also ends in a line cut short, which the restart removes.
flume_add_restart_test(restart_np1_np1 REFERENCE restart_reference_np1 INPUTS ${rs_main}
    FIRST_NPROCS 1 RESTART_NPROCS 1
    FIRST_ARGS flume.check_int=400 flume.max_step=550
    RESTART_ARGS flume.check_int=400 flume.restart=chk00400
    CHECK_ARGS --gauges 2 --rewound --cut-line)
# The checkpoint of the last step: the first part stops at step 400 with
# flume.check_int = 1000.
flume_add_restart_test(restart_np2_np2 REFERENCE restart_reference_np2 INPUTS ${rs_main}
    FIRST_NPROCS 2 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=1000 flume.max_step=400
    RESTART_ARGS flume.check_int=400 flume.restart=chk00400
    CHECK_ARGS --gauges 2)
# The checkpoint of a last step that is not a multiple of flume.log_int (10)
# or flume.plot_int (100): the first part wrote log lines and a plotfile at
# step 405 that the run without a break does not have; the restart removes
# the log lines (plt00405 stays).
flume_add_restart_test(restart_np1_np1_off_interval REFERENCE restart_reference_np1
    INPUTS ${rs_main}
    FIRST_NPROCS 1 RESTART_NPROCS 1
    FIRST_ARGS flume.check_int=1000 flume.max_step=405
    RESTART_ARGS flume.check_int=400 flume.restart=chk00405
    CHECK_ARGS --gauges 2)
# Another rank count, other boxes, or both: the plotfiles and gauges to the
# last bit, the balance to round-off (relative to the largest value of a
# column).
flume_add_restart_test(restart_np1_np2 REFERENCE restart_reference_np1 INPUTS ${rs_main}
    FIRST_NPROCS 1 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=400 flume.max_step=550
    RESTART_ARGS flume.restart=chk00400
    CHECK_ARGS --gauges 2 --rewound --balance-rtol 1e-14)
flume_add_restart_test(restart_np2_np1_boxes_8 REFERENCE restart_reference_np1 INPUTS ${rs_main}
    FIRST_NPROCS 2 RESTART_NPROCS 1
    FIRST_ARGS flume.check_int=400 flume.max_step=550
    RESTART_ARGS flume.restart=chk00400 flume.max_grid_size=8
    CHECK_ARGS --gauges 2 --rewound --balance-rtol 1e-14)
flume_add_restart_test(restart_np2_np2_boxes_8 REFERENCE restart_reference_np1 INPUTS ${rs_main}
    FIRST_NPROCS 2 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=400 flume.max_step=550
    RESTART_ARGS flume.restart=chk00400 flume.max_grid_size=8
    CHECK_ARGS --gauges 2 --rewound --balance-rtol 1e-14)
# First order, Green-Ampt.
flume_add_restart_test(restart_first_order_np2 REFERENCE restart_first_order_reference
    INPUTS ${rs_first}
    FIRST_NPROCS 2 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=300 flume.max_step=450
    RESTART_ARGS flume.check_int=300 flume.restart=chk00300
    CHECK_ARGS --gauges 0 --rewound)

# --- Inputs and checkpoints that do not match ---------------------------------
flume_add_abort_test(abort_check_int_negative
    EXPECT "flume\\.check_int must be >= 0 and <= 2147483647. it is -1"
    ARGS ${ok} flume.check_int=-1)
flume_add_abort_test(abort_check_int_two
    EXPECT "flume\\.check_int needs exactly 1 value. it has 2"
    ARGS ${ok} flume.check_int=1\ 2)
flume_add_abort_test(abort_check_file_empty EXPECT "flume\\.check_file must not be empty"
    ARGS ${ok} flume.check_file=\"\" flume.check_int=1)
flume_add_abort_test(abort_check_file_without_int
    EXPECT "flume\\.check_file needs flume\\.check_int > 0"
    ARGS ${ok} flume.check_file=mychk)
flume_add_abort_test(abort_check_file_plot_file
    EXPECT "flume\\.check_file and flume\\.plot_file are both 'out': they must differ"
    ARGS ${ok} flume.check_int=10 flume.check_file=out flume.plot_file=out)
flume_add_abort_test(abort_restart_empty EXPECT "flume\\.restart must not be empty"
    ARGS ${ok} flume.restart=\"\")
flume_add_abort_test(abort_restart_two EXPECT "flume\\.restart needs exactly 1 value. it has 2"
    ARGS ${ok} flume.restart=a\ b)
flume_add_abort_test(abort_restart_missing
    EXPECT "flume\\.restart: no_such_chk: no file no_such_chk/Header: not a FLUME checkpoint"
    ARGS ${ok} flume.restart=no_such_chk)

set(rs_ref1 ${FLUME_TEST_RUN_DIR}/restart_reference_np1)
set(rs_ref_first ${FLUME_TEST_RUN_DIR}/restart_first_order_reference)
# A function, not a macro: a macro would remove one more level of the
# backslashes in EXPECT.
function(rs_abort name fixture)
    flume_add_abort_test(${name} ${ARGN})
    if(TEST ${name} AND TEST ${fixture})
        set_tests_properties(${name} PROPERTIES FIXTURES_REQUIRED ${fixture})
    elseif(TEST ${name})
        set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    endif()
endfunction()
rs_abort(abort_restart_plotfile restart_reference_np1
    EXPECT "plt00000/Header is not 'FLUME_checkpoint': not a FLUME checkpoint"
    ARGS ${rs_main} flume.restart=${rs_ref1}/plt00000)
rs_abort(abort_restart_grid restart_reference_np1
    EXPECT "its grid is 48 x 32 cells, not the 24 x 32 of flume\\.n_cell"
    ARGS ${rs_main} flume.restart=${rs_ref1}/chk00400 flume.n_cell=24\ 32)
rs_abort(abort_restart_domain restart_reference_np1
    EXPECT "its domain differs from that of flume\\.prob_lo and flume\\.prob_hi"
    ARGS ${rs_main} flume.restart=${rs_ref1}/chk00400 flume.prob_hi=96\ 64.5)
rs_abort(abort_restart_bed restart_reference_np1
    EXPECT "its bed differs from the bed of the inputs"
    ARGS ${rs_main} flume.restart=${rs_ref1}/chk00400 flume.bed.noise=0.03)
rs_abort(abort_restart_sources restart_reference_np1
    EXPECT "its sources \\(rain curve_number_infiltration manning_friction\\) differ from those of the inputs \\(rain green_ampt_infiltration manning_friction\\)"
    ARGS ${rs_first} flume.restart=${rs_ref1}/chk00400)
rs_abort(abort_restart_periodic restart_first_order_reference
    EXPECT "its periodic sides differ from those of flume\\.periodic"
    ARGS ${rs_first} flume.restart=${rs_ref_first}/chk00300 flume.periodic=0\ 1)
rs_abort(abort_restart_columns restart_first_order_reference
    EXPECT "the water-balance columns of the checkpoint \\( rain infiltrated \\) differ"
    ARGS ${rs_first} flume.restart=${rs_ref_first}/chk00300 flume.bc.segments=out
         flume.bc.out.side=xhi flume.bc.out.type=outflow)
rs_abort(abort_restart_maps restart_first_order_reference
    EXPECT "written with flume\\.flood_maps = 0: it holds no maps to continue"
    ARGS ${rs_first} flume.restart=${rs_ref_first}/chk00300 flume.flood_maps=1)
rs_abort(abort_restart_segments restart_reference_np1
    EXPECT "its outflow segments [(] out:xhi:4:11 [)] cover other faces than those of the inputs [(] out:xhi:4:10 [)]"
    ARGS ${rs_main} flume.restart=${rs_ref1}/chk00400 flume.bc.out.range=8\ 22)
rs_abort(abort_restart_gauges restart_reference_np1
    EXPECT "its gauges [(]mid:23:15 low:46:7[)] differ from those of the inputs [(]mid:24:15 low:46:7[)]"
    ARGS ${rs_main} flume.restart=${rs_ref1}/chk00400 flume.gauges.mid.location=49\ 31)
rs_abort(abort_restart_no_balance restart_first_order_no_log
    EXPECT "written with flume\\.log_int = 0, which keeps no water balance"
    ARGS ${rs_first}
         flume.restart=${FLUME_TEST_RUN_DIR}/restart_first_order_no_log/chk00300)

# Checkpoints with one change each (Tests/Restart/make_bad_checkpoints.py),
# copied from chk00300 of restart_first_order_reference.
if(TEST restart_first_order_reference)
    set(rs_bad ${FLUME_TEST_RUN_DIR}/restart_bad_checkpoints)
    add_test(NAME restart_bad_checkpoints
             COMMAND ${CMAKE_COMMAND} -E make_directory ${rs_bad})
    add_test(NAME restart_bad_checkpoints_write
             COMMAND ${Python3_EXECUTABLE} ${rs_dir}/make_bad_checkpoints.py
                     ${rs_ref_first}/chk00300
             WORKING_DIRECTORY ${rs_bad})
    set_tests_properties(restart_bad_checkpoints PROPERTIES
                         FIXTURES_REQUIRED restart_first_order_reference
                         FIXTURES_SETUP restart_bad_checkpoints)
    set_tests_properties(restart_bad_checkpoints_write PROPERTIES
                         FIXTURES_REQUIRED "restart_first_order_reference;restart_bad_checkpoints"
                         FIXTURES_SETUP restart_bad_checkpoints_write
                         LABELS regression)
endif()
rs_abort(abort_restart_version restart_bad_checkpoints_write
    EXPECT "written by checkpoint version 2. this build reads 1"
    ARGS ${rs_first} flume.restart=${rs_bad}/bad_version)
rs_abort(abort_restart_precision restart_bad_checkpoints_write
    EXPECT "written by a (single|double)-precision build. this build is (double|single) precision"
    ARGS ${rs_first} flume.restart=${rs_bad}/bad_precision)
rs_abort(abort_restart_duplicate restart_bad_checkpoints_write
    EXPECT "the key time appears twice in"
    ARGS ${rs_first} flume.restart=${rs_bad}/bad_duplicate)
rs_abort(abort_restart_step restart_bad_checkpoints_write
    EXPECT "bad_step: step is not valid"
    ARGS ${rs_first} flume.restart=${rs_bad}/bad_step)
rs_abort(abort_restart_nan restart_bad_checkpoints_write
    EXPECT "its depth or discharge is not finite in 1 cell"
    ARGS ${rs_first} flume.restart=${rs_bad}/bad_nan)
