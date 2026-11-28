# The town case: the rain of the storm case (Tests/ERFStorm) over a town of
# buildings, with ERF's immersed forcing for the air and FLUME's buildings
# for the flood, with FLUME_WITH_ERF. Included from Tests/CMakeLists.txt
# after ERFStorm/ERFStormTests.cmake, whose erf_flume_storm writes the rain.
# See Tests/ERFStormTown/README.md.

set(erf_town_dir ${CMAKE_CURRENT_SOURCE_DIR}/ERFStormTown)
set(erf_town_setup ${FLUME_TEST_RUN_DIR}/erf_flume_town_setup)

if(TEST erf_flume_storm)
    # The square and the window of most rain in the storm's town, the index of
    # their rain, and the buildings of ERF and of FLUME.
    add_test(NAME erf_flume_town_setup
        COMMAND ${Python3_EXECUTABLE} ${erf_town_dir}/prepare_town.py
                ${FLUME_TEST_RUN_DIR}/erf_flume_storm/rain/storm.txt ${erf_town_setup})
    set_tests_properties(erf_flume_town_setup PROPERTIES LABELS "erf;setup"
                         FIXTURES_REQUIRED erf_flume_storm_rain
                         FIXTURES_SETUP erf_flume_town_setup)

    # erf_flume_town: 1600 ERF steps of 0.75 s (1200 s) on 4 ranks, three FLUME
    # steps of 0.25 s each. The rain of the coupled run is that of the storm's
    # grids over the square; its water balance closes; ERF's solid cells are
    # FLUME's building cells, and the air inside the buildings is at rest;
    # FLUME alone on the same grid (4 ranks) gives the coupled flood to the
    # bit; FLUME without the buildings closes its balance on the same rain.
    # The water balances are checked with the factor 1024 (check_town.py
    # says why).
    erf_flume_add_regtest(erf_flume_town NPROCS 4 BALANCE_FACTOR 1024
        INPUTS ${erf_town_dir}/inputs_town
        ARGS FILE=${erf_town_setup}/town.inputs FILE=${erf_town_dir}/flume_town
             erf.input_sounding_file=${erf_town_dir}/sounding_town
        CHECK ${Python3_EXECUTABLE} ${erf_town_dir}/check_town.py ${erf_town_setup} 1600
              --standalone ${launcher_np4} $<TARGET_FILE:flume> ${erf_town_dir}/inputs_standalone
                           FILE=${erf_town_dir}/flume_town FILE=${erf_town_setup}/town_flume.inputs
                           ${FLUME_TEST_ARGS}
              --no-buildings ${launcher_np4} $<TARGET_FILE:flume> ${erf_town_dir}/inputs_standalone
                             FILE=${erf_town_dir}/flume_town
                             FILE=${erf_town_setup}/town_flume_open.inputs ${FLUME_TEST_ARGS})
    if(TEST erf_flume_town)
        set_tests_properties(erf_flume_town PROPERTIES FIXTURES_REQUIRED erf_flume_town_setup
                             TIMEOUT 1500)
    endif()
endif()
