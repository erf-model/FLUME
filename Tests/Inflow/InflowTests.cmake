# ---------------------------------------------------------------------------
# Inflow hydrographs on boundary segments (flume.bc.<name>.type = discharge;
# Tests/Inflow/README.md). Included from Tests/CMakeLists.txt, which defines
# the functions and ${ok}, after Tests/Restart/RestartTests.cmake, which
# defines flume_add_restart_test.
#
# Tests/Inflow/check_inflow.py compares the last plotfile of a channel run with a
# closed form; Tests/check_balance.py pins the inflow column of the water
# balance against the integral of the hydrograph.
# ---------------------------------------------------------------------------
set(in_dir ${CMAKE_CURRENT_SOURCE_DIR}/Inflow)
set(in_channel ${in_dir}/inputs_channel)
set(in_basin ${in_dir}/inputs_basin)
set(check_in ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Inflow/check_inflow.py)
set(check_bal ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py)

# --- Subcritical inflow -----------------------------------------------------
# Uniform flow, 1 m deep at 1 m/s (Froude number 0.32), on the flat channel
# without friction, fed with 1 m2/s: an exact steady state of the scheme. The
# inflow face keeps it: every cell holds the initial state to round-off
# after 300 s (3419 steps at first order), on 4 ranks in 4 boxes. A copy
# (zero_gradient) side would also keep this state; the test guards the rounding of
# the face depth and the momentum flux, and the other tests the inflow. The
# rounding of the face depth gathers over the steps: about 250 epsilons in
# the depth and 1030 in the discharge at both orders; the bound is 4096.
foreach(o 1 2)
    flume_add_regtest(inflow_uniform_subcritical_order${o} NPROCS 4
        INPUTS ${in_channel}
        ARGS flume.order=${o} "flume.n_cell=200 2" flume.max_grid_size=50 flume.stop_time=300
             flume.bc.in.file=${in_dir}/hydrograph_steady.txt
             flume.init.surface_elevation=11 "flume.init.velocity=1 0"
        CHECK ${check_in} chan --side xlo --q 1 --mode uniform --depth 1 --tol-eps 4096)
endforeach()
# The same flow when the discharge drops to 0.5 m2/s at t = 0
# (hydrograph_half.txt): a rarefaction runs downstream, and next to the side
# the depth is h_b = 0.869019 m, from u - 2 c kept across it, with the new
# discharge. At 20 s the plateau reaches 70 m; cells up to 50 m are compared,
# within 0.1 % (measured: 0.022 % in the depth and 0.068 % in the discharge
# at first order, 0.006 % and 0.028 % at second). On each
# of the four sides (the flow turned to enter across it), at both orders.
set(in_rare ARGS flume.stop_time=20 flume.bc.in.file=${in_dir}/hydrograph_half.txt
             flume.init.surface_elevation=11)
set(in_rare_check --q 0.5 --mode plateau --h0 1 --u0 1 --range 0 50 --rtol 1e-3)
foreach(o 1 2)
    flume_add_regtest(inflow_rarefaction_xlo_order${o} NPROCS 2
        INPUTS ${in_channel}
        ${in_rare} flume.order=${o} "flume.init.velocity=1 0"
        CHECK ${check_in} chan --side xlo ${in_rare_check})
    flume_add_regtest(inflow_rarefaction_xhi_order${o} NPROCS 2
        INPUTS ${in_channel}
        ${in_rare} flume.order=${o} "flume.init.velocity=-1 0"
             flume.bc.in.side=xhi flume.bc.out.side=xlo
        CHECK ${check_in} chan --side xhi ${in_rare_check})
    flume_add_regtest(inflow_rarefaction_ylo_order${o} NPROCS 2
        INPUTS ${in_channel}
        ${in_rare} flume.order=${o} "flume.init.velocity=0 1"
             "flume.n_cell=2 400" "flume.prob_hi=4 200" "flume.periodic=1 0"
             flume.bc.in.side=ylo flume.bc.out.side=yhi
        CHECK ${check_in} chan --side ylo ${in_rare_check})
    flume_add_regtest(inflow_rarefaction_yhi_order${o} NPROCS 2
        INPUTS ${in_channel}
        ${in_rare} flume.order=${o} "flume.init.velocity=0 -1"
             "flume.n_cell=2 400" "flume.prob_hi=4 200" "flume.periodic=1 0"
             flume.bc.in.side=yhi flume.bc.out.side=ylo
        CHECK ${check_in} chan --side yhi ${in_rare_check})
endforeach()

# --- Supercritical inflow ---------------------------------------------------
# 1 m2/s onto the dry, flat channel without friction: the face is at the
# critical depth (0.4671 m) and a rarefaction runs from it to the dry front
# (3 c_c t = 193 m at 30 s), with u + 2 c kept across it. Cells up to 60 m,
# within 2 % at first order and 1 % at second (measured: 1.1 % and 0.75 %,
# largest next to the side, where the fan starts at the critical speed).
foreach(case "1|0.02" "2|0.01")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 o)
    list(GET pair 1 tol)
    flume_add_regtest(inflow_dry_fan_order${o} NPROCS 2
        INPUTS ${in_channel}
        ARGS flume.order=${o} flume.stop_time=30
             flume.bc.in.file=${in_dir}/hydrograph_steady.txt
        CHECK ${check_in} chan --side xlo --q 1 --mode fan --time 30 --range 0 60 --rtol ${tol})
endforeach()
# Normal depth: 1 m2/s down a slope of 0.02 with n = 0.03 s/m^(1/3), from
# dry, 600 s, 200 cells of 1 m. The normal depth, 0.394424 m (Froude number
# 1.29), is reached from 40 m on. At first order the bed source drives the
# flow with g S (h - S dx / 2), which raises the depth by 0.77 %; at second
# order the depth is 0.21 % low (0.11 % and 0.06 % at 0.5 m and 0.25 m
# cells) and the discharge within 0.35 %. Tolerances 1 % and 0.5 %. The outflow side bends the last cells;
# cells from 40 m to 190 m are compared.
foreach(case "1|0.01" "2|0.005")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 o)
    list(GET pair 1 tol)
    flume_add_regtest(inflow_normal_depth_steep_order${o} NPROCS 2
        INPUTS ${in_channel}
        ARGS flume.order=${o} "flume.n_cell=200 2" flume.stop_time=600
             flume.friction.manning=0.03 "flume.bed.slope=-0.02 0"
             flume.bc.in.file=${in_dir}/hydrograph_steady.txt
        CHECK ${check_in} chan --side xlo --q 1 --mode normal --manning 0.03 --slope 0.02
              --range 40 190 --rtol ${tol})
endforeach()

# --- Mass in = the integral of the hydrograph --------------------------------
# The closed basin of inputs_basin: segment a takes hydrograph_triangle.txt
# (linear: 200 m3), segment b hydrograph_steps.txt (step: 260 m3). At 300 s
# the inflow columns hold the integrals and the basin 800 + 460 m3, to
# round-off, on 4 ranks in 9 boxes. The steps do not land on the times of the
# files.
flume_add_regtest(inflow_basin_volume_np4 NPROCS 4
    INPUTS ${in_basin}
    ARGS flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
         flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    OUTPUT "inflow a: hydrograph from [^\n]*hydrograph_triangle\\.txt, 3 entries, linear in time, up to 2 m3/s over 20 m of xlo"
    CHECK ${check_bal} --expect-columns stored in_a in_b
          --expect-final in_a=200 in_b=260 stored=1260)
# The other interpolation for each file: 280 m3 and 205 m3.
flume_add_regtest(inflow_basin_volume_swapped NPROCS 1
    INPUTS ${in_basin}
    ARGS flume.bc.a.file=${in_dir}/hydrograph_triangle.txt flume.bc.a.interpolation=step
         flume.bc.b.file=${in_dir}/hydrograph_steps.txt flume.bc.b.interpolation=linear
    OUTPUT "inflow b: hydrograph from [^\n]*hydrograph_steps\\.txt, 3 entries, linear in time, up to 3 m3/s over 10 m of yhi"
    CHECK ${check_bal} --expect-final in_a=280 in_b=205 stored=1285)
# A file whose first time, 50 s, is after the start of the run, with a
# discharge above 0: it enters from t = 0 (600 m3 by 300 s), and the run says
# so.
flume_add_regtest(inflow_basin_late_start NPROCS 2
    INPUTS ${in_basin}
    ARGS flume.bc.a.file=${in_dir}/hydrograph_late.txt
         flume.bc.b.file=${in_dir}/hydrograph_zero.txt
    OUTPUT "inflow a: note: the first time of [^\n]*hydrograph_late\\.txt is 50 s; its discharge, 2 m3/s, enters from the start of the run"
    CHECK ${check_bal} --expect-final in_a=600 in_b=0 stored=1400)
# Lake at rest: no discharge on either segment. The face of a segment then
# keeps the depth of the cell inside it with no velocity, so the still water
# stays still: every cell holds 0.5 m and no discharge to round-off. A wall
# would pass too: this guards the q = 0 branch of the face depth (no 0 / 0,
# no spurious momentum flux), not the inflow itself.
flume_add_regtest(inflow_basin_lake_at_rest NPROCS 2
    INPUTS ${in_basin}
    ARGS flume.bc.a.file=${in_dir}/hydrograph_zero.txt
         flume.bc.b.file=${in_dir}/hydrograph_zero.txt
    CHECK ${check_in} basin --side xlo --q 0 --mode uniform --depth 0.5)
# The dry channel fed by a hydrograph, with no stop time (inputs_channel has
# none): the discharge sets the step of the dry domain, and three steps are
# taken.
flume_add_regtest(inflow_dry_no_stop_time NPROCS 1
    INPUTS ${in_channel}
    ARGS flume.max_step=3 flume.log_int=1 flume.bc.in.file=${in_dir}/hydrograph_steady.txt
    CHECK ${check_bal} --expect-last-step 3)

# --- Restart ----------------------------------------------------------------
# inputs_basin to 300 s with checkpoints every 100 steps; a restart from the
# checkpoint at step 100 (20.9 s: hydrograph a rising, b at its first
# 1 m3/s) gives every plotfile to the last bit and, on the same ranks, the
# log byte for byte; on another rank count the log agrees to 1e-14 of each
# column (Tests/Restart/check_restart.py, which reads plotfiles named plt).
flume_add_regtest(inflow_restart_reference NPROCS 2
    INPUTS ${in_basin}
    ARGS flume.check_int=100 flume.plot_file=plt flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
         flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    CHECK ${check_rs} --expect chk00100)
if(TEST inflow_restart_reference)
    set_tests_properties(inflow_restart_reference PROPERTIES
                         FIXTURES_SETUP inflow_restart_reference)
endif()
flume_add_restart_test(inflow_restart_np2 REFERENCE inflow_restart_reference INPUTS ${in_basin}
    FIRST_NPROCS 2 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=100 flume.max_step=150 flume.plot_file=plt
               flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
               flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    RESTART_ARGS flume.check_int=100 flume.restart=chk00100 flume.plot_file=plt
                 flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
                 flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    CHECK_ARGS --gauges 0 --rewound)
flume_add_restart_test(inflow_restart_np2_np1 REFERENCE inflow_restart_reference INPUTS ${in_basin}
    FIRST_NPROCS 2 RESTART_NPROCS 1
    FIRST_ARGS flume.check_int=100 flume.max_step=150 flume.plot_file=plt
               flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
               flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    RESTART_ARGS flume.check_int=100 flume.restart=chk00100 flume.plot_file=plt
                 flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
                 flume.bc.b.file=${in_dir}/hydrograph_steps.txt
    CHECK_ARGS --gauges 0 --rewound --balance-rtol 1e-14)

# --- Aborts -----------------------------------------------------------------
# The keys of a discharge segment, on the 8 x 8 m deck ${ok}.
set(in_seg flume.bc.segments=a flume.bc.a.side=xlo)
flume_add_abort_test(abort_inflow_file_missing_key
    EXPECT "flume\\.bc\\.a\\.type = discharge needs flume\\.bc\\.a\\.file: lines of a time in s and a discharge in m3/s"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge)
flume_add_abort_test(abort_inflow_file_with_outflow
    EXPECT "flume\\.bc\\.a\\.file is read only with flume\\.bc\\.a\\.type = discharge"
    ARGS ${ok} ${in_seg} flume.bc.a.type=outflow flume.bc.a.file=${in_dir}/hydrograph_steady.txt)
flume_add_abort_test(abort_inflow_interpolation_with_wall
    EXPECT "flume\\.bc\\.a\\.interpolation is read only with flume\\.bc\\.a\\.type = discharge"
    ARGS ${ok} ${in_seg} flume.bc.a.type=wall flume.bc.a.interpolation=step)
# The type is checked before the value of the interpolation.
flume_add_abort_test(abort_inflow_interpolation_unknown_with_wall
    EXPECT "flume\\.bc\\.a\\.interpolation is read only with flume\\.bc\\.a\\.type = discharge"
    ARGS ${ok} ${in_seg} flume.bc.a.type=wall flume.bc.a.interpolation=cubic)
flume_add_abort_test(abort_inflow_file_two_values
    EXPECT "flume\\.bc\\.a\\.file needs exactly 1 value"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge "flume.bc.a.file=x y")
flume_add_abort_test(abort_inflow_interpolation_two_values
    EXPECT "flume\\.bc\\.a\\.interpolation needs exactly 1 value"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_dir}/hydrograph_steady.txt
         "flume.bc.a.interpolation=step linear")
flume_add_abort_test(abort_inflow_interpolation_unknown
    EXPECT "flume\\.bc\\.a\\.interpolation must be step or linear. it is cubic"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_dir}/hydrograph_steady.txt
         flume.bc.a.interpolation=cubic)
flume_add_abort_test(abort_inflow_file_empty_name
    EXPECT "flume\\.bc\\.a\\.file must not be empty: the name of a hydrograph file"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge "flume.bc.a.file=\"\"")
# The file itself: the bad files of the hyetograph tests (Tests/Rain/bad),
# which the same reader takes.
set(in_bad ${CMAKE_CURRENT_SOURCE_DIR}/Rain/bad)
set(ibad "flume\\.bc\\.a\\.file \\([^)]*\\)")
flume_add_abort_test(abort_inflow_file_not_found
    EXPECT "flume\\.bc\\.a\\.file \\([^)]*no_such_file\\.txt\\): the file cannot be opened"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_bad}/no_such_file.txt)
flume_add_abort_test(abort_inflow_file_directory
    EXPECT "flume\\.bc\\.a\\.file \\([^)]*bad\\): the path is not a regular file"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_bad})
foreach(case
        "one_item|, line 2: a line holds two numbers, a time in s and a discharge in m3/s; it holds 1 item\\(s\\)"
        "three_items|, line 2: a line holds two numbers, a time in s and a discharge in m3/s; it holds 3 item\\(s\\)"
        "time_not_a_number|, line 2: the time must be a finite number; it is \"sixty\""
        "time_beyond_double|, line 2: the time must be a finite number; it is \"1e400\""
        "rate_inf|, line 2: the discharge must be a finite number; it is \"inf\""
        "negative|, line 2: the discharge must be >= 0 m3/s; it is -1"
        "not_increasing|, line 3: the time must be above the time of the entry before \\(60 s\\); it is 30"
        "equal_times|, line 3: the time must be above the time of the entry before \\(60 s\\); it is 60"
        "empty|: the file holds no entry: it needs at least one line with a time in s and a discharge in m3/s")
    string(REPLACE "|" ";" pair "${case}")
    list(GET pair 0 file)
    list(GET pair 1 expect)
    flume_add_abort_test(abort_inflow_file_${file} EXPECT "${ibad}${expect}"
        ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_bad}/${file}.txt)
endforeach()
# A dry domain, no stop time, and no discharge from now on: nothing sets the
# step.
flume_add_abort_test(abort_inflow_zero_dry_without_stop_time
    EXPECT "and no rain falls and no water enters through a discharge segment or an inflow [(]flume\\.inflows\\.names[)] from this time on, so nothing sets the time step: set flume\\.stop_time"
    ARGS ${ok} ${in_seg} flume.bc.a.type=discharge flume.bc.a.file=${in_dir}/hydrograph_zero.txt
         flume.max_step=1)
# A restart whose segment b covers more faces than in the run that wrote
# the checkpoint: the balance of its inflow cannot go on.
rs_abort(abort_inflow_restart_segment_range inflow_restart_reference
    EXPECT "its outflow and discharge segments [(] a:xlo:5:14 b:yhi:0:4 [)] cover other faces than those of the inputs [(] a:xlo:5:14 b:yhi:0:9 [)]"
    ARGS ${in_basin} flume.restart=${FLUME_TEST_RUN_DIR}/inflow_restart_reference/chk00100
         flume.bc.a.file=${in_dir}/hydrograph_triangle.txt
         flume.bc.b.file=${in_dir}/hydrograph_steps.txt "flume.bc.b.range=0 20")
