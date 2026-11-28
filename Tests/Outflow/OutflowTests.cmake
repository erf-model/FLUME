# ---------------------------------------------------------------------------
# Open sides that let no water in (flume.bc.<name>.type = outflow) and the
# plain copy that does (zero_gradient); Tests/Outflow/README.md. Included
# from Tests/CMakeLists.txt, which defines flume_add_regtest.
#
# Tests/Outflow/check_outflow.py reads the out_<name> columns of the water
# balance: with --no-inflow none may fall or go below zero and the domain may
# not store more than it held and was given; with --takes-in the comparison
# run must have taken water in.
# ---------------------------------------------------------------------------
set(out_dir ${CMAKE_CURRENT_SOURCE_DIR}/Outflow)
set(check_out ${Python3_EXECUTABLE} ${out_dir}/check_outflow.py)

# --- A bed that steps down to the side ---------------------------------------
# Still water, 0.6 m deep in the boundary cells of side xlo and 0.5 m inside,
# moving at 1e-4 m/s into the domain (Tests/Outflow/inputs_step). With the
# plain copy the boundary face passes the discharge of the whole depth and
# the face inside only that above the step, so the boundary cell fills, its
# ghost cells rise with it and the inflow grows: by 300 s the 200.8 m3 have
# become 8974 m3 (zero_gradient; at least 10 times is required). An outflow
# segment acts as a wall while the discharge of the boundary cell points
# inward: no water enters, at either order. The flow then turns and the
# pond drains through the side: 147 m3 (143 m3 at second order) have left
# by 300 s, and at least 100 m3 is required, which a wall would not pass.
foreach(o 1 2)
    flume_add_regtest(outflow_step_no_inflow_order${o} NPROCS 1
        INPUTS ${out_dir}/inputs_step
        ARGS flume.order=${o} flume.bed.file=${out_dir}/bed_step.asc flume.bc.w.type=outflow
        OUTPUT "boundary segment w: outflow on xlo, cells 0 to 1 in y"
        CHECK ${check_out} --no-inflow --min-out 100)
endforeach()
flume_add_regtest(outflow_step_zero_gradient_takes_in NPROCS 1
    INPUTS ${out_dir}/inputs_step
    ARGS flume.bed.file=${out_dir}/bed_step.asc flume.bc.w.type=zero_gradient
    OUTPUT "boundary segment w: zero_gradient on xlo, cells 0 to 1 in y"
    CHECK ${check_out} --takes-in 10)

# --- A bed that falls into the domain from both open sides -------------------
# Rain of 50 mm/h for 1200 s on the valley of Tests/Outflow/inputs_valley,
# 200 m2: 50 / 3.6e6 m/s * 200 m2 * 1200 s = 10/3 m3. The film runs away
# from sides xlo and xhi, so an outflow segment keeps all of it: the stored
# volume is the rain to 1e-12, and nothing crosses either side. With the
# plain copy each side lets in 2.0e-3 m3, so the domain stores 1.0012 times
# the rain (zero_gradient; at least 1.001 times is required). Two ranks,
# both orders.
foreach(o 1 2)
    flume_add_regtest(outflow_valley_keeps_rain_order${o} NPROCS 2
        INPUTS ${out_dir}/inputs_valley
        ARGS flume.order=${o} flume.bed.file=${out_dir}/bed_valley.asc
             flume.bc.west.type=outflow flume.bc.east.type=outflow
        CHECK ${check_out} --no-inflow --expect-stored 3.3333333333333335)
endforeach()
flume_add_regtest(outflow_valley_zero_gradient_takes_in NPROCS 2
    INPUTS ${out_dir}/inputs_valley
    ARGS flume.bed.file=${out_dir}/bed_valley.asc
         flume.bc.west.type=zero_gradient flume.bc.east.type=zero_gradient
    CHECK ${check_out} --takes-in 1.001)

# --- A pond beside two open sides ---------------------------------------------
# Rain on the rough plane of Tests/Outflow/inputs_pond collects against the
# wall at x = 200 m and reaches sides ylo and yhi. With the plain copy the
# pond takes water in through one side from about 890 s on: at 1200 s it
# stores 8069 m3 from 333 m3 of rain on the 3 ranks of the test (the growth
# starts from round-off, so the number depends on the boxes; zero_gradient;
# at least 5 times the rain is required). With outflow segments no water enters and the pond
# drains through both sides (42 m3 by 1200 s; at least 30 m3 is required).
# Three ranks and 8 boxes, so boxes of different
# ranks meet both sides; both orders.
foreach(o 1 2)
    flume_add_regtest(outflow_pond_no_inflow_order${o} NPROCS 3
        INPUTS ${out_dir}/inputs_pond
        ARGS flume.order=${o} flume.bc.s.type=outflow flume.bc.n.type=outflow
        CHECK ${check_out} --no-inflow --min-out 30)
endforeach()
flume_add_regtest(outflow_pond_zero_gradient_takes_in NPROCS 3
    INPUTS ${out_dir}/inputs_pond
    ARGS flume.bc.s.type=zero_gradient flume.bc.n.type=zero_gradient
    CHECK ${check_out} --takes-in 5)
