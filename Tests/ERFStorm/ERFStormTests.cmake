# The coupled example: a convective cell from ERF's Kessler microphysics over
# an urban catchment, with the replay of its rain and the same storm as a
# hyetograph, with FLUME_WITH_ERF. Included from Tests/CMakeLists.txt after
# ERF/ERFTests.cmake, which defines launcher_np4.
# See Tests/ERFStorm/README.md.

set(erf_storm_dir ${CMAKE_CURRENT_SOURCE_DIR}/ERFStorm)
set(erf_storm_rasters ${FLUME_TEST_RUN_DIR}/erf_flume_storm_catchment)
# FLUME's half of the case, with the rasters of prepare_catchment.py.
set(erf_storm_flume FILE=${erf_storm_dir}/flume_catchment
    flume.bed.file=${erf_storm_rasters}/bed.asc
    flume.friction.manning_file=${erf_storm_rasters}/manning.asc
    flume.infiltration.cn.number_file=${erf_storm_rasters}/cn.asc)

if(FLUME_MPI)
    # The bed, curve number and Manning n of the catchment (256 x 128 cells).
    add_test(NAME erf_flume_storm_catchment
        COMMAND ${Python3_EXECUTABLE} ${erf_storm_dir}/prepare_catchment.py ${erf_storm_rasters})
    set_tests_properties(erf_flume_storm_catchment PROPERTIES LABELS "erf;setup"
                         FIXTURES_SETUP erf_flume_storm_catchment)
endif()

# erf_flume_storm: 960 ERF steps of 5 s (80 min) on 4 ranks, one FLUME step
# each, the rain written to rain/storm.txt. ERF's rain on the catchment, at
# every plt2d plotfile, is that of the hyetograph made from the written rain,
# to 1e-12; the stand-alone replay (4 ranks, 8 boxes, 8 h) gives the coupled
# flood at 40 and 80 min to round-off; the hyetograph run closes its balance
# and rains the hyetograph's depth; the replay's outlet peak is at least 3
# times the hyetograph run's and at least 1 h earlier.
erf_flume_add_regtest(erf_flume_storm NPROCS 4
    INPUTS ${erf_storm_dir}/inputs_storm
    ARGS ${erf_storm_flume} erf.input_sounding_file=${erf_storm_dir}/sounding_storm
    CHECK ${Python3_EXECUTABLE} ${erf_storm_dir}/check_storm.py rain/storm.txt 5 960
          --replay ${launcher_np4} $<TARGET_FILE:flume> ${erf_storm_dir}/inputs_standalone
                   ${erf_storm_flume} flume.rain.type=gridded ${FLUME_TEST_ARGS}
          --hyetograph ${launcher_np4} $<TARGET_FILE:flume> ${erf_storm_dir}/inputs_standalone
                       ${erf_storm_flume} flume.rain.type=hyetograph ${FLUME_TEST_ARGS})
if(TEST erf_flume_storm)
    set_tests_properties(erf_flume_storm PROPERTIES FIXTURES_REQUIRED erf_flume_storm_catchment)
endif()
