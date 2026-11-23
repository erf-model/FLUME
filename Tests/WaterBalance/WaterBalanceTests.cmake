# ---------------------------------------------------------------------------
# Water-balance log, flood maps and gauges (Tests/WaterBalance/README.md).
# Included from Tests/CMakeLists.txt, which defines the functions and ${ok}.
# Every regression run that leaves water_balance.txt also has its residual
# checked by Tests/run_case.py; these tests pin the columns themselves
# against closed forms.
# ---------------------------------------------------------------------------
set(wb_dir ${CMAKE_CURRENT_SOURCE_DIR}/WaterBalance)
set(check_wb ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py)
set(check_flood ${Python3_EXECUTABLE} ${wb_dir}/check_flood_outputs.py)

# Rain on a closed basin (Tests/Rain/inputs_flat, and the same on a bed at
# 5 m with maps and gauges in inputs_rain_gauges): the rain column is
# 36 mm/h x 600 s x 1536 m2 = 9.216 m3, all of it stored. Every cell's
# maximum depth is the 6 mm fallen, its speeds 0, its arrival time the first
# step end with more than 3 mm (within one step after 300 s), and the three
# gauges read the rain fallen at every step, and 5 m plus that as the
# surface elevation. The face gauge (x = 24 m, between
# cells 11 and 12) reads cell 12; the corner gauge (48, 32) m, on the high
# edges, the last cell (23, 15); the gauge at (3.5, 0) m, on the low edge,
# cell (1, 0). On 1 and 2 ranks.
foreach(np 1 2)
    flume_add_regtest(balance_rain_gauges_np${np} NPROCS ${np}
        INPUTS ${wb_dir}/inputs_rain_gauges
        OUTPUT "gauge face at \\(24, 7\\) m reads cell \\(12, 3\\)"
        CHECK ${check_flood} rg --rain-basin 36 0.003 --bed 5
              --gauge face 12 3 --gauge corner 23 15 --gauge inside 1 0)
    flume_add_regtest(balance_rain_columns_np${np} NPROCS ${np}
        INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Rain/inputs_flat
        ARGS flume.plot_int=100000
        CHECK ${check_wb} --expect-columns stored rain --expect-final rain=9.216 stored=9.216
              --rtol 1e-14)
endforeach()

# A gauge on two faces whose coordinates round below them: on 30 cells of
# 0.1 m, x = 0.3 m gives (x - 0) / 3 m x 30 = 2.9999999999999996 and
# y = 0.7 m gives 6.999999999999999 in double; within 4 epsilons of the
# largest coordinate of a face the point is on it and reads the cell above,
# (3, 7), not (2, 6).
flume_add_regtest(gauge_face_rounding NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.n_cell=30\ 30 flume.prob_lo=0\ 0 flume.prob_hi=3\ 3
         flume.gauges.names=g flume.gauges.g.location=0.3\ 0.7
    OUTPUT "gauge g at \\(0\\.[0-9]+, 0\\.[0-9]+\\) m reads cell \\(3, 7\\)"
    CHECK ${check_wb} --expect-columns stored --expect-lines 1)

# Gauges on a field that differs in every cell: the skeleton grid (6 boxes
# of 16 x 16 cells of 2 m on [-10, 86] x [5, 69] m) with a bed sloped in x
# and y under a level surface at 3 m, moving at (0.4, -0.3) m/s; no step is
# taken. Three gauges in three boxes, on 2 ranks, must hold exactly the
# depth, surface elevation and discharge over depth of their cells in the
# step-0 plotfile: (0.5, 10) m is cell (5, 2), (60.1, 50.3) m cell (35, 22),
# and (85.9, 68.9) m the last cell (47, 31).
flume_add_regtest(gauge_cells_sloped_bed NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.bed.slope=0.01\ 0.02 flume.init.surface_elevation=3 flume.init.velocity=0.4\ -0.3
         flume.gauges.names=a\ b\ c flume.gauges.a.location=0.5\ 10
         flume.gauges.b.location=60.1\ 50.3 flume.gauges.c.location=85.9\ 68.9
    CHECK ${Python3_EXECUTABLE} ${wb_dir}/check_gauges.py skel --gauge a 5 2 --gauge b 35 22
          --gauge c 47 31)

# A gauge far from the origin in single precision: on prob_lo = 500000 m and
# cells of 1 m, x = 500003.7 m is cell 3. Its rounding (1/32 m in single
# precision) is a large part of a cell, but the point is not near a face, so
# it stays in cell 3.
flume_add_regtest(gauge_large_coordinates NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Skeleton/inputs
    ARGS flume.n_cell=16\ 16 flume.prob_lo=500000\ 4100000 flume.prob_hi=500016\ 4100016
         flume.gauges.names=g flume.gauges.g.location=500003.7\ 4100008.5
    OUTPUT "gauge g at \\([0-9.e+]+, [0-9.e+]+\\) m reads cell \\(3, 8\\)"
    CHECK ${check_wb} --expect-columns stored --expect-lines 1)

# Ponded Green-Ampt infiltration (Tests/Infiltration/inputs_ponded): 6 h on
# 102400 m2 infiltrate F = 91.7517138 mm, 9395.3754951 m3, from the
# implicit equation K t = F - M ln(1 + F / M); the stored water is the rest
# of the 51200 m3. The start lags by 3.7e-6 s (the first step;
# Tests/Infiltration/InfiltrationTests.cmake), 1.1e-6 m3, inside 1e-10 of
# the stored volume.
flume_add_regtest(balance_green_ampt_ponded NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Infiltration/inputs_ponded
    ARGS flume.plot_int=100000
    CHECK ${check_wb} --expect-columns stored infiltrated
          --expect-final infiltrated=9395.375495144983 stored=41804.624504855017 --rtol 1e-10)

# Uniform supercritical flow out of a side (inputs_outflow): 20 m3 per
# metre of outflow side leave by 2 s, so 400 m3 of the 6000 m3. Out of the
# whole side at x = 300 m at first order; out of two segments of 5 m and
# 15 m of it at second order (100 and 300 m3; the two stages of a step each
# count half of it); out of the low side y = 0 at second order on the
# transposed domain, where the discharge, and the flux, is negative.
flume_add_regtest(balance_outflow_side NPROCS 1
    INPUTS ${wb_dir}/inputs_outflow
    ARGS flume.bc.segments=east flume.bc.east.side=xhi flume.bc.east.type=outflow
    CHECK ${check_wb} --expect-columns stored out_east --expect-final out_east=400 stored=5600
          --rtol 1e-13)
flume_add_regtest(balance_outflow_two_segments NPROCS 2
    INPUTS ${wb_dir}/inputs_outflow
    ARGS flume.order=2 flume.bc.segments=a\ b
         flume.bc.a.side=xhi flume.bc.a.type=outflow flume.bc.a.range=0\ 5
         flume.bc.b.side=xhi flume.bc.b.type=outflow flume.bc.b.range=5\ 20
    CHECK ${check_wb} --expect-columns stored out_a out_b
          --expect-final out_a=100 out_b=300 stored=5600 --rtol 1e-13)
flume_add_regtest(balance_outflow_low_side NPROCS 2
    INPUTS ${wb_dir}/inputs_outflow
    ARGS flume.order=2 flume.n_cell=20\ 300 flume.prob_hi=20\ 300 flume.init.velocity=0\ -10
         flume.bc.segments=south flume.bc.south.side=ylo flume.bc.south.type=outflow
    CHECK ${check_wb} --expect-columns stored out_south --expect-final out_south=400 stored=5600
          --rtol 1e-13)
# A periodic domain has no outflow column: the uniform flow of
# Tests/Uniform/inputs keeps its 12288 m3 (2 m on 6144 m2) to the last bit.
flume_add_regtest(balance_periodic NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.plot_int=100000
    CHECK ${check_wb} --expect-columns stored --expect-final stored=12288 --rtol 0 --factor 0)
# log_int = 5 on 16 steps: lines at steps 0, 5, 10 and 15, and the last
# step, 16.
flume_add_regtest(balance_log_interval NPROCS 1
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Uniform/inputs
    ARGS flume.plot_int=100000 flume.log_int=5
    CHECK ${check_wb} --expect-columns stored --expect-last-step 16 --expect-lines 5)

# The log off (flume.log_int = 0) on the every-term case for 400 steps:
# no water_balance.txt and no summary line, and plotfiles byte for byte
# those of the same run with a line every step.
set(wb_log_args flume.order=2 flume.plot_int=100 flume.max_step=400
    flume.infiltration.type=green_ampt flume.infiltration.green_ampt.conductivity=2
    flume.infiltration.green_ampt.suction=110 flume.infiltration.green_ampt.deficit=0.3)
flume_add_regtest(balance_log_on NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Rain/inputs_plane
    ARGS ${wb_log_args}
    OUTPUT "FLUME: water balance \\(m3\\)"
    CHECK ${check_wb} --expect-last-step 400 --expect-lines 401)
flume_add_regtest(balance_log_off NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Rain/inputs_plane
    ARGS ${wb_log_args} flume.log_int=0
    NO_OUTPUT "water balance"
    CHECK ${Python3_EXECUTABLE} ${wb_dir}/check_no_log.py plane ${FLUME_TEST_RUN_DIR}/balance_log_on)
if(TEST balance_log_on AND TEST balance_log_off)
    set_tests_properties(balance_log_on PROPERTIES FIXTURES_SETUP balance_log_on)
    set_tests_properties(balance_log_off PROPERTIES FIXTURES_REQUIRED balance_log_on)
endif()

# Every term at once: rain (50 mm/h for 900 s on 200 m2, 2.5 m3), Green-Ampt
# infiltration and an outflow side, on the plane of Tests/Rain/inputs_plane,
# second order. The residual stays at round-off, and the log without any one
# of its three terms fails the check (Tests/WaterBalance/check_dropped_terms.py).
flume_add_regtest(balance_every_term NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Rain/inputs_plane
    ARGS flume.order=2 flume.plot_int=100000 flume.infiltration.type=green_ampt
         flume.infiltration.green_ampt.conductivity=2 flume.infiltration.green_ampt.suction=110
         flume.infiltration.green_ampt.deficit=0.3
    CHECK ${Python3_EXECUTABLE} ${wb_dir}/check_dropped_terms.py --min-columns 3)
# The same with a line every 25 steps, where the bound grows with the steps
# between lines: each term left out still fails.
flume_add_regtest(balance_every_term_log_interval NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/Rain/inputs_plane
    ARGS flume.order=2 flume.plot_int=100000 flume.infiltration.type=green_ampt
         flume.infiltration.green_ampt.conductivity=2 flume.infiltration.green_ampt.suction=110
         flume.infiltration.green_ampt.deficit=0.3 flume.log_int=25
    CHECK ${Python3_EXECUTABLE} ${wb_dir}/check_dropped_terms.py --min-columns 3)

# Flood maps of the wet-bed dam break (Tests/DamBreak/inputs): 2 m against
# 1 m at x = 50 m. Behind the shock (S = 4.18 m/s) the state is hm =
# 1.454 m and um = 1.306 m/s; with an arrival depth of 1.25 m, between the
# two depths, the arrival time of each cell the shock has passed is
# (x - 50 m) / S within one cell over S, its maximum depth and speed are hm
# and um within 2 % (and their product, the maximum depth times speed), and
# the cells ahead of it hold -1. A gauge at x = 60 m (cell 240; on the face
# y = 0.5 m, so row 2) ends in that state. First and second order.
flume_add_regtest(flood_maps_dam_break NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.flood_maps=1 flume.flood_maps.arrival_depth=1.25
         flume.gauges.names=behind flume.gauges.behind.location=60\ 0.5
    CHECK ${check_flood} dam --dam-break 2 1 50 9.81 1.25 --gauge behind 240 2)
flume_add_regtest(flood_maps_dam_break_second_order NPROCS 2
    INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/DamBreak/inputs
    ARGS flume.order=2 flume.flood_maps=1 flume.flood_maps.arrival_depth=1.25
         flume.gauges.names=behind flume.gauges.behind.location=60\ 0.5
    CHECK ${check_flood} dam --dam-break 2 1 50 9.81 1.25 --gauge behind 240 2)

# The inputs.
flume_add_abort_test(abort_log_int_negative
    EXPECT "flume\\.log_int must be >= 0 and <= 2147483647. it is -1"
    ARGS ${ok} flume.log_int=-1)
flume_add_abort_test(abort_flood_maps_two
    EXPECT "flume\\.flood_maps must be >= 0 and <= 1. it is 2"
    ARGS ${ok} flume.flood_maps=2)
flume_add_abort_test(abort_arrival_depth_negative
    EXPECT "flume\\.flood_maps\\.arrival_depth must be >= flume\\.dry_depth \\([0-9.e-]+\\). it is -0\\.001"
    ARGS ${ok} flume.flood_maps=1 flume.flood_maps.arrival_depth=-0.001)
flume_add_abort_test(abort_arrival_depth_below_dry
    EXPECT "flume\\.flood_maps\\.arrival_depth must be >= flume\\.dry_depth \\(0\\.0(1|0999[0-9]*)\\). it is 0\\.00(5|499[0-9]*)"
    ARGS ${ok} flume.dry_depth=0.01 flume.flood_maps=1 flume.flood_maps.arrival_depth=0.005)
flume_add_abort_test(abort_gauge_case
    EXPECT "flume\\.gauges\\.names lists Ab and aB, which differ only in case"
    ARGS ${ok} flume.gauges.names=Ab\ aB flume.gauges.Ab.location=1\ 1 flume.gauges.aB.location=2\ 2)
flume_add_abort_test(abort_arrival_depth_nan
    EXPECT "flume\\.flood_maps\\.arrival_depth must be finite"
    ARGS ${ok} flume.flood_maps=1 flume.flood_maps.arrival_depth=nan)
flume_add_abort_test(abort_arrival_depth_without_maps
    EXPECT "flume\\.flood_maps\\.arrival_depth needs flume\\.flood_maps = 1"
    ARGS ${ok} flume.flood_maps.arrival_depth=0.01)
flume_add_abort_test(abort_gauge_name
    EXPECT "flume\\.gauges\\.names: a name may hold only letters, digits and _. it is \"a-b\""
    ARGS ${ok} flume.gauges.names=a-b)
flume_add_abort_test(abort_gauge_twice
    EXPECT "flume\\.gauges\\.names lists g twice"
    ARGS ${ok} flume.gauges.names=g\ g flume.gauges.g.location=1\ 1)
flume_add_abort_test(abort_gauge_no_location
    EXPECT "flume\\.gauges\\.names lists g, which needs flume\\.gauges\\.g\\.location"
    ARGS ${ok} flume.gauges.names=g)
flume_add_abort_test(abort_gauge_outside_x
    EXPECT "flume\\.gauges\\.g\\.location must be inside the domain: from 0 to 8 in x [^\n]* it is 8\\.5"
    ARGS ${ok} flume.gauges.names=g flume.gauges.g.location=8.5\ 1)
flume_add_abort_test(abort_gauge_outside_y
    EXPECT "flume\\.gauges\\.g\\.location must be inside the domain: from 0 to 8 in y [^\n]* it is -1"
    ARGS ${ok} flume.gauges.names=g flume.gauges.g.location=1\ -1)
flume_add_abort_test(abort_gauge_location_one_value
    EXPECT "flume\\.gauges\\.g\\.location needs exactly 2 values \\(x and y\\). it has 1"
    ARGS ${ok} flume.gauges.names=g flume.gauges.g.location=1)
flume_add_abort_test(abort_gauge_unlisted
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^\n]*\n *flume\\.gauges\\.h\\.location"
    ARGS ${ok} flume.gauges.names=g flume.gauges.g.location=1\ 1 flume.gauges.h.location=1\ 1)
flume_add_abort_test(abort_gauges_without_log
    EXPECT "flume\\.gauges\\.names needs flume\\.log_int > 0"
    ARGS ${ok} flume.log_int=0 flume.gauges.names=g flume.gauges.g.location=1\ 1)
