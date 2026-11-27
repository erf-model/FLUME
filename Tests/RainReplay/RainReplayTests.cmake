# The writer of flume.rain.replay_file (flume::RainReplay) on its own: a
# small driver linked to flume_lib hands it negative depths, a step of no
# rain, a step that starts an instant before the end of the last one, and a
# restart while a negative depth is held. Runs in every build; the coupled
# test is Tests/ERFReplay. See Tests/RainReplay/README.md.

add_executable(replay_unit ${CMAKE_CURRENT_SOURCE_DIR}/RainReplay/replay_unit.cpp)
target_link_libraries(replay_unit PRIVATE flume_lib)

# rain_replay_unit: the index and grids of four steps, to the byte and the
# bit; the restart from the checkpoint after step 2 writes them again; a
# restart with an index that is not the run's, and a step that starts before
# the end of the last, stop. Two ranks (two boxes each) with MPI.
if(FLUME_MPI)
    set(replay_unit_np 2)
else()
    set(replay_unit_np 1)
endif()
flume_launcher(replay_unit_launcher ${replay_unit_np})
add_test(NAME rain_replay_unit
    COMMAND ${Python3_EXECUTABLE} ${FLUME_RUN_CASE}
            --dir ${FLUME_TEST_RUN_DIR}/rain_replay_unit
            --run ${replay_unit_launcher} $<TARGET_FILE:replay_unit> run
            --check ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/RainReplay/check_replay_unit.py
                    ${replay_unit_launcher} $<TARGET_FILE:replay_unit>)
set_tests_properties(rain_replay_unit PROPERTIES LABELS "regression"
                     PROCESSORS ${replay_unit_np} TIMEOUT 120)
