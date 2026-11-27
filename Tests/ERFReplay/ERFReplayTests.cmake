# Replay of the rain from ERF's microphysics (flume.rain.replay_file), with
# FLUME_WITH_ERF. Included from Tests/CMakeLists.txt after
# ERFRain/ERFRainTests.cmake, whose deck, sounding and helpers it uses. See
# Tests/ERFReplay/README.md.

set(erf_replay_dir ${CMAKE_CURRENT_SOURCE_DIR}/ERFReplay)
set(check_erf_replay ${Python3_EXECUTABLE} ${erf_replay_dir}/check_replay.py)
# The volume case of Tests/ERFRain on 4 ranks (ERF's boxes of 8 x 8 x 32
# cells), with ERF's step fixed at 0.25 s (below ERF's CFL limit, about
# 0.34 s) and FLUME's at most that: one FLUME step per ERF step. The rain is
# written to rain/replay.txt.
set(erf_replay_args FILE=${erf_rain_dir}/flume_basin ${erf_rain_sounding}
    amr.max_grid_size=8\ 8\ 32 erf.fixed_dt=0.25 flume.max_dt=0.25
    flume.rain.replay_file=rain/replay.txt)

# erf_flume_rain_replay: 120 ERF steps (30 s) of Kessler rain on the dam
# break, writing the rain. The index closes at 30 s and its grids hold the
# rain of the water balance. The stand-alone run (inputs_replay, 4 ranks, 16
# boxes) that reads it back gives the coupled flood at 15 s and 30 s to
# round-off, and runs on to 60 s with no more rain while the water moves.
# ERF then continues from its checkpoint of step 60 and writes the same index
# and grids, to the byte; the restart without the index, or from a
# checkpoint without what the writer holds back, stops.
erf_flume_add_regtest(erf_flume_rain_replay NPROCS 4
    INPUTS ${erf_rain_dir}/inputs_kessler
    ARGS ${erf_replay_args}
    CHECK ${check_erf_replay} rain/replay.txt 0.25 120
          --replay ${launcher_np4} $<TARGET_FILE:flume> ${erf_replay_dir}/inputs_replay
                   FILE=${erf_rain_dir}/flume_basin flume.rain.type=gridded ${FLUME_TEST_ARGS}
          --restart chk00060 ${launcher_np4} $<TARGET_FILE:erf_exec> ${erf_rain_dir}/inputs_kessler
                    ${erf_replay_args} erf.restart=chk00060 ${FLUME_TEST_ARGS})

# flume.rain.replay_file must name a file (the grids are named after it)
# whose name a line of the index can hold (no blank, no '#'). One rank, one
# step of the cell-by-cell deck of Tests/ERFRain.
function(erf_flume_replay_abort name value expect)
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
                --dir ${FLUME_TEST_RUN_DIR}/${name}
                --expect-abort ${expect}
                --run ${launcher_np1} $<TARGET_FILE:erf_exec> ${erf_rain_dir}/inputs_kessler
                      FILE=${erf_rain_dir}/flume_cells ${erf_rain_sounding}
                      flume.rain.replay_file=${value} max_step=1 erf.check_int=-1
                      ${FLUME_TEST_ARGS} amrex.call_addr2line=0)
    set_tests_properties(${name} PROPERTIES LABELS "erf;abort" TIMEOUT 120)
endfunction()
erf_flume_replay_abort(abort_erf_flume_rain_replay_directory rain/
    "flume\\.rain\\.replay_file: rain/: not the name of a file")
erf_flume_replay_abort(abort_erf_flume_rain_replay_hash "rain/a#b.txt"
    "flume\\.rain\\.replay_file: rain/a#b\\.txt: its name holds a blank or a '#'")
