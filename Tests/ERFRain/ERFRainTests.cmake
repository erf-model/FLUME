# Rain from ERF's microphysics (flume.rain.type = microphysics), with
# FLUME_WITH_ERF. Included from Tests/CMakeLists.txt after ERF/ERFTests.cmake,
# whose erf_flume_add_regtest and erf_flume_add_abort_test it uses. See
# Tests/ERFRain/README.md.

set(erf_rain_dir ${CMAKE_CURRENT_SOURCE_DIR}/ERFRain)
set(check_erf_rain ${Python3_EXECUTABLE} ${erf_rain_dir}/check_rain.py)
# ERF reads the sounding from the directory of the run.
set(erf_rain_sounding erf.input_sounding_file=${erf_rain_dir}/sounding_moist)

# erf_flume_rain_cells: Kessler rain, 16 x 16 columns, each over 2 x 2 flood
# cells. No water moves (flume.dry_depth = 1 m). FLUME takes two steps in
# each ERF step of about 0.34 s (flume.max_dt = 0.2 s). After 60 and 120 ERF
# steps each flood cell holds the accumulated precipitation of its column /
# 1000 m, and the precipitation varies at least 1.5-fold over the columns.
erf_flume_add_regtest(erf_flume_rain_cells NPROCS 2
    INPUTS ${erf_rain_dir}/inputs_kessler
    ARGS FILE=${erf_rain_dir}/flume_cells ${erf_rain_sounding} erf.check_int=-1
    CHECK ${check_erf_rain} cells 2 60 120 --min-spread 1.5)

# erf_flume_rain_volume: the same rain on a flat periodic bed with friction,
# on a dam break of 0.1 m: the water moves. The rain of the water balance (closed to round-off by
# run_case.py) is the precipitation ERF accumulated over the grid. Then ERF
# continues from its checkpoint of step 60, with FLUME's state and ERF's
# precipitation reference in it, and gives the same flood plotfiles and
# water balance to the last bit.
erf_flume_add_regtest(erf_flume_rain_volume NPROCS 2
    INPUTS ${erf_rain_dir}/inputs_kessler
    ARGS FILE=${erf_rain_dir}/flume_basin ${erf_rain_sounding}
    CHECK ${check_erf_rain} volume 2 120
          --restart chk00060 ${launcher_np2} $<TARGET_FILE:erf_exec> ${erf_rain_dir}/inputs_kessler
                    FILE=${erf_rain_dir}/flume_basin ${erf_rain_sounding} erf.restart=chk00060
                    ${FLUME_TEST_ARGS})

# erf_flume_rain_frozen: WDM6 rain, snow and graupel over ground below
# freezing; no water moves. Each flood cell holds its column's precipitation,
# liquid and frozen as liquid water, / 1000 m, and the frozen part is at least
# a fifth of the total (it is all of it in this run; erf_flume_rain_cells has
# only liquid rain, so between them each part is checked).
erf_flume_add_regtest(erf_flume_rain_frozen NPROCS 2
    INPUTS ${erf_rain_dir}/inputs_wdm6
    ARGS FILE=${erf_rain_dir}/flume_cells
    CHECK ${check_erf_rain} cells 1 100 200 --min-frozen 0.2)

# Without moisture there is no precipitation to take: ERF stops the run at
# the start, with its own message. The atmosphere of Tests/ERF (no moisture
# model) with the FLUME deck of erf_flume_rain_cells.
flume_launcher(launcher_np1 1)
add_test(NAME abort_erf_flume_rain_no_moisture
    COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
            --dir ${FLUME_TEST_RUN_DIR}/abort_erf_flume_rain_no_moisture
            --expect-abort "takes the precipitation of every column from the microphysics surface precipitation accumulators, but erf\\.moisture_model = \"None\" provides none"
            --run ${launcher_np1} $<TARGET_FILE:erf_exec> ${erf_case_dir}/inputs_erf
                  FILE=${erf_rain_dir}/flume_cells max_step=1 ${FLUME_TEST_ARGS}
                  amrex.call_addr2line=0)
set_tests_properties(abort_erf_flume_rain_no_moisture PROPERTIES LABELS "erf;abort" TIMEOUT 120)

# Inside ERF, flume.rain.rate with microphysics stops the run, naming the key
# (the shared check of the flume.rain.* keys; a stand-alone run stops first,
# on the rain type).
add_test(NAME abort_erf_flume_rain_rate
    COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
            --dir ${FLUME_TEST_RUN_DIR}/abort_erf_flume_rain_rate
            --expect-abort "flume\\.rain\\.rate is read only with flume\\.rain\\.type = constant"
            --run ${launcher_np1} $<TARGET_FILE:erf_exec> ${erf_rain_dir}/inputs_kessler
                  FILE=${erf_rain_dir}/flume_cells ${erf_rain_sounding} flume.rain.rate=1
                  max_step=1 erf.check_int=-1 ${FLUME_TEST_ARGS} amrex.call_addr2line=0)
set_tests_properties(abort_erf_flume_rain_rate PROPERTIES LABELS "erf;abort" TIMEOUT 120)
