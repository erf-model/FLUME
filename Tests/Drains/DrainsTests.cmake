# Tests of the storm-drain inlets (flume.drains.*), included from
# Tests/CMakeLists.txt, whose functions (flume_add_regtest,
# flume_add_abort_test; flume_add_restart_test of Tests/Restart) and ${rok}
# deck they use.
# CMAKE_CURRENT_SOURCE_DIR is Tests/. Tests/Drains/README.md describes the
# cases; Tests/Drains/check_drains.py computes the expected numbers from the
# inverse of the closed form, independently of the solver.
set(check_dr ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Drains/check_drains.py)
set(dr_dir ${CMAKE_CURRENT_SOURCE_DIR}/Drains)

# Four cells, each with an inlet, still water 0.8 m deep
# (Tests/Drains/inputs_uniform): orifice flow to 0.4469 m at 30.4 s, then
# weir flow, to 8.3 mm at 600 s. The water stays level and still, so every
# cell must hold the closed-form depth at the time of each plotfile, to 16
# epsilons of 0.8 m (the run is within 1 epsilon), with no discharge.
flume_add_regtest(drains_uniform NPROCS 4
    INPUTS ${dr_dir}/inputs_uniform
    OUTPUT "drain inlet d at \\(15, 15\\) m in cell \\(1, 1\\): perimeter 2 m, open area 0\\.5 m2, weir coefficient 1\\.66, orifice coefficient 0\\.67; orifice above a depth of 0\\.44694720[0-9]* m"
    CHECK ${check_dr} uni --uniform 0.8 --grate 2 0.5 --all-inlets --n-plotfiles 6)
# The same from 0.3 m, below the switch depth: weir flow only.
flume_add_regtest(drains_uniform_weir NPROCS 1
    INPUTS ${dr_dir}/inputs_uniform
    ARGS flume.init.surface_elevation=0.3
    CHECK ${check_dr} uni --uniform 0.3 --grate 2 0.5 --all-inlets --n-plotfiles 5)
# The same with the coefficients set for every inlet (1.5 and 0.6, the
# switch at 0.4429 m), on 1 rank.
set(dr_coef "")
foreach(n a b c d)
    list(APPEND dr_coef flume.drains.${n}.weir_coefficient=1.5
                        flume.drains.${n}.orifice_coefficient=0.6)
endforeach()
flume_add_regtest(drains_uniform_coefficients NPROCS 1
    INPUTS ${dr_dir}/inputs_uniform
    ARGS ${dr_coef}
    OUTPUT "weir coefficient 1\\.5, orifice coefficient 0\\.6; orifice above a depth of 0\\.44294469[0-9]* m"
    CHECK ${check_dr} uni --uniform 0.8 --grate 2 0.5 1.5 0.6 --all-inlets --n-plotfiles 6)
# The water balance of drains_uniform: 400 m2 x (0.8 m - 8.3395 mm) drained
# by 600 s (the closed form, Tests/Drains/check_drains.py), the rest stored.
flume_add_regtest(balance_drains_uniform NPROCS 4
    INPUTS ${dr_dir}/inputs_uniform
    ARGS flume.plot_int=100000
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py
          --expect-columns stored drained
          --expect-final drained=316.6641978030054 stored=3.335802196994654 --rtol 1e-12)

# One inlet in a closed basin of 16 x 16 cells (Tests/Drains/inputs_basin),
# still water 0.3 m deep, 20 min. The inlet cell is drawn below the mean
# depth, so the basin drains less than one cell of its area would: the
# drained volume must be 0.95 to 1 times that of the lumped closed form at
# every plotfile (the run gives 0.979 to 0.987). The upper bound holds while
# the inlet cell is the lowest point of the water (the rate grows with the
# depth); the lower bound is set from the run, with a margin. At first order on 4 ranks
# (it writes the reference checkpoints of the restart test below) and at
# second order on 2.
flume_add_regtest(drains_basin NPROCS 4
    INPUTS ${dr_dir}/inputs_basin
    ARGS flume.check_int=2000
    CHECK ${check_dr} plt --basin 0.3 0.95 1.0 --grate 0.5 0.05 --inlet 8 8 --n-plotfiles 14)
if(TEST drains_basin)
    set_tests_properties(drains_basin PROPERTIES FIXTURES_SETUP drains_basin)
endif()
flume_add_regtest(drains_basin_order2 NPROCS 2
    INPUTS ${dr_dir}/inputs_basin
    ARGS flume.order=2
    CHECK ${check_dr} plt --basin 0.3 0.95 1.0 --grate 0.5 0.05 --inlet 8 8 --n-plotfiles 14)
# The basin on 3 ranks in boxes of 4 (16 boxes in place of 4): every value of
# every plotfile is that of drains_basin to the last bit.
flume_add_regtest(drains_basin_boxes_4_np3 NPROCS 3
    INPUTS ${dr_dir}/inputs_basin
    ARGS flume.max_grid_size=4
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Restart/check_restart.py
          --compare ${FLUME_TEST_RUN_DIR}/drains_basin)
if(TEST drains_basin_boxes_4_np3)
    set_tests_properties(drains_basin_boxes_4_np3 PROPERTIES FIXTURES_REQUIRED drains_basin)
endif()
# Continued from the checkpoint of step 2000 on 2 ranks: the plotfiles to the
# last bit, drained_depth among them, and the balance to round-off. The
# drained depth and its carry come from the checkpoint.
flume_add_restart_test(drains_restart_np4_np2 REFERENCE drains_basin INPUTS ${dr_dir}/inputs_basin
    FIRST_NPROCS 4 RESTART_NPROCS 2
    FIRST_ARGS flume.check_int=2000 flume.max_step=2500
    RESTART_ARGS flume.restart=chk02000
    CHECK_ARGS --gauges 0 --rewound --balance-rtol 1e-14)

# A 1 mm film and an inlet that could take 13 mm in a step
# (Tests/Drains/inputs_film): no depth below 0 in any of the 21 plotfiles,
# every cell keeps its depth plus drained_depth, and the depth drained in the
# first step (still water, nothing flows in) is that of the closed form, not
# the whole film.
flume_add_regtest(drains_film NPROCS 1
    INPUTS ${dr_dir}/inputs_film
    CHECK ${check_dr} film --grate 100 1 --inlet 1 2 --first-step 0.001 --n-plotfiles 21)

# A restart with other inlets stops: the inlet of drains_basin moved to cell
# (3, 3), or given another weir coefficient.
function(dr_restart_abort name)
    flume_add_abort_test(${name} ${ARGN})
    if(TEST ${name} AND TEST drains_basin)
        set_tests_properties(${name} PROPERTIES FIXTURES_REQUIRED drains_basin)
    elseif(TEST ${name})
        set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    endif()
endfunction()
set(dr_chk flume.restart=${FLUME_TEST_RUN_DIR}/drains_basin/chk02000)
dr_restart_abort(abort_restart_drains_moved
    EXPECT "its drain inlets [(]grate:8:8:[^)]*[)] differ from those of the inputs [(]grate:3:3:"
    ARGS ${dr_dir}/inputs_basin ${dr_chk} flume.drains.grate.location=3.5\ 3.5)
dr_restart_abort(abort_restart_drains_coefficient
    EXPECT "its drain inlets [(]grate:8:8:0\\.5:0\\.05[0-9]*:1\\.6[0-9]*:0\\.67[0-9]*[)] differ from those of the inputs [(]grate:8:8:0\\.5:0\\.05[0-9]*:1\\.5:0\\.67[0-9]*[)]"
    ARGS ${dr_dir}/inputs_basin ${dr_chk} flume.drains.grate.weir_coefficient=1.5)

# Drain inputs, on the 8 x 8 cells of [0, 40]^2 m (cells of 25 m2).
set(dr flume.drains.names=g flume.drains.g.location=12\ 12 flume.drains.g.perimeter=1
       flume.drains.g.area=0.1)
set(drk "flume\\.drains\\.g")
flume_add_abort_test(abort_drains_name
    EXPECT "flume\\.drains\\.names: a name may hold only letters, digits and _; it is \"g-1\""
    ARGS ${rok} flume.drains.names=g-1)
flume_add_abort_test(abort_drains_name_twice
    EXPECT "flume\\.drains\\.names lists g twice"
    ARGS ${rok} ${dr} flume.drains.names=g\ g)
set(dr_loc flume.drains.g.location=12\ 12)
set(dr_per flume.drains.g.perimeter=1)
set(dr_area flume.drains.g.area=0.1)
flume_add_abort_test(abort_drains_missing_location
    EXPECT "flume\\.drains\\.names lists g, which needs ${drk}\\.location: the x and y of the inlet, m"
    ARGS ${rok} flume.drains.names=g ${dr_per} ${dr_area})
flume_add_abort_test(abort_drains_missing_perimeter
    EXPECT "flume\\.drains\\.names lists g, which needs ${drk}\\.perimeter: the weir length of the inlet, m"
    ARGS ${rok} flume.drains.names=g ${dr_loc} ${dr_area})
flume_add_abort_test(abort_drains_missing_area
    EXPECT "flume\\.drains\\.names lists g, which needs ${drk}\\.area: the open area of the inlet, m2"
    ARGS ${rok} flume.drains.names=g ${dr_loc} ${dr_per})
flume_add_abort_test(abort_drains_location_outside
    EXPECT "${drk}\\.location must be inside the domain: from 0 to 40 in y \\(flume\\.prob_lo and flume\\.prob_hi\\); it is 40\\.5"
    ARGS ${rok} ${dr} flume.drains.g.location=12\ 40.5)
flume_add_abort_test(abort_drains_perimeter_zero
    EXPECT "${drk}\\.perimeter must be > 0; it is 0"
    ARGS ${rok} ${dr} flume.drains.g.perimeter=0)
flume_add_abort_test(abort_drains_area_zero
    EXPECT "${drk}\\.area must be > 0 and <= the area of a cell \\(25 m2\\); it is 0"
    ARGS ${rok} ${dr} flume.drains.g.area=0)
flume_add_abort_test(abort_drains_area_cell
    EXPECT "${drk}\\.area must be > 0 and <= the area of a cell \\(25 m2\\); it is 25\\.5"
    ARGS ${rok} ${dr} flume.drains.g.area=25.5)
flume_add_abort_test(abort_drains_weir_coefficient
    EXPECT "${drk}\\.weir_coefficient must be > 0; it is -1"
    ARGS ${rok} ${dr} flume.drains.g.weir_coefficient=-1)
flume_add_abort_test(abort_drains_orifice_coefficient_zero
    EXPECT "${drk}\\.orifice_coefficient must be > 0 and <= 1; it is 0"
    ARGS ${rok} ${dr} flume.drains.g.orifice_coefficient=0)
flume_add_abort_test(abort_drains_orifice_coefficient_above_1
    EXPECT "${drk}\\.orifice_coefficient must be > 0 and <= 1; it is 1\\.2"
    ARGS ${rok} ${dr} flume.drains.g.orifice_coefficient=1.2)
flume_add_abort_test(abort_drains_unknown_key
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.drains\\.g\\.perimiter"
    ARGS ${rok} ${dr} flume.drains.g.perimiter=1)
# A key of a name that flume.drains.names does not list.
flume_add_abort_test(abort_drains_unlisted_name
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.drains\\.other\\.location"
    ARGS ${rok} ${dr} flume.drains.other.location=1\ 1)
# Two inlets in one cell (12, 12) and (14, 14) m: cell (2, 2) of 5 m.
flume_add_abort_test(abort_drains_same_cell
    EXPECT "flume\\.drains\\.h\\.location puts inlet h in cell \\(2, 2\\), which holds inlet g already"
    ARGS ${rok} ${dr} flume.drains.names=g\ h flume.drains.h.location=14\ 14
         flume.drains.h.perimeter=1 flume.drains.h.area=0.1)
# An inlet in a building cell: cell (12, 30) of the streets of
# Tests/Raster/inputs_streets.
flume_add_abort_test(abort_drains_building
    EXPECT "${drk}\\.location puts inlet g in a building cell \\(flume\\.buildings\\.file\\)"
    ARGS flume.n_cell=40\ 40 flume.prob_hi=40\ 40 ${dr} flume.drains.g.location=12.5\ 30.5
         flume.buildings.file=${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc)
