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

# --- The approach-flow law and inlets from a list file ----------------------
# Sixteen cells, each with an inlet of the approach-flow law, in water 0.5 m
# deep that moves at (1.2, 0.5) m/s (Tests/Drains/inputs_approach, the
# inlets from Tests/Drains/inlets_all.txt). The water stays level and keeps
# its velocity, so every cell must hold the closed-form depth at the time of
# each plotfile, to 64 epsilons of 0.5 m (the runs are within 15: the powers
# of the closed form round more than the square roots of the weir law), with a
# discharge of the depth times the velocity. It is the reference of the
# tests below.
set(check_same ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Restart/check_restart.py)
set(ap_all flume.drains.street.file=${dr_dir}/inlets_all.txt)
set(ap_check --uniform 0.5 --velocity 1.2 0.5 --all-inlets --tol-eps 64)
flume_add_regtest(drains_approach_uniform NPROCS 4
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all}
    OUTPUT "drain inlets street: 16 from [^ ]*inlets_all\\.txt in 16 cell\\(s\\), 0 with more than one inlet: approach flow, area 0\\.4 m2, coefficient 0\\.302, Froude exponent -0\\.816"
    CHECK ${check_dr} plt ${ap_check} --approach 0.4 --n-plotfiles 7)
if(TEST drains_approach_uniform)
    set_tests_properties(drains_approach_uniform PROPERTIES FIXTURES_SETUP drains_approach_uniform)
endif()
# Two inlets of half the area in every cell (Tests/Drains/inlets_twice.txt):
# a cell drains at the rate of its two inlets, so every value of every
# plotfile is that of drains_approach_uniform to the last bit.
flume_add_regtest(drains_approach_merge NPROCS 2
    INPUTS ${dr_dir}/inputs_approach
    ARGS flume.drains.street.file=${dr_dir}/inlets_twice.txt flume.drains.street.area=0.2
    OUTPUT "drain inlets street: 32 from [^ ]*inlets_twice\\.txt in 16 cell\\(s\\), 16 with more than one inlet"
    CHECK ${check_same} --compare ${FLUME_TEST_RUN_DIR}/drains_approach_uniform)
if(TEST drains_approach_merge)
    set_tests_properties(drains_approach_merge PROPERTIES FIXTURES_REQUIRED drains_approach_uniform)
endif()
# The same with another coefficient and exponent (0.226, the value Annis et
# al. calibrated, and -0.5).
flume_add_regtest(drains_approach_coefficients NPROCS 1
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.drains.street.coefficient=0.226 flume.drains.street.froude_exponent=-0.5
    CHECK ${check_dr} plt ${ap_check} --approach 0.4 0.226 -0.5 --n-plotfiles 7)
# To 120 s with a plotfile at every step: the cells empty at 87.05 s. No
# depth is below 0, no cell is empty in the plotfiles before that time,
# every cell is empty after it.
flume_add_regtest(drains_approach_empty NPROCS 2
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.stop_time=120 flume.plot_int=1
    CHECK ${check_dr} plt ${ap_check} --approach 0.4 --n-plotfiles 100)
# The same in 6492 steps (a CFL number of 0.011): a step lands in the last
# 0.04 s before 87.05 s, where the closed form is at or below the dry depth
# of 1e-6 m. The step takes that film, so the cells are empty afterwards
# (a film the law could not take would stay: it holds no discharge).
flume_add_regtest(drains_approach_empty_small_steps NPROCS 4
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.stop_time=120 flume.cfl=0.011 flume.plot_int=500
    CHECK ${check_dr} plt ${ap_check} --approach 0.4 --n-plotfiles 10)
# An exponent of 0, the upper end of its range: Q = k u A, a linear law
# (to 0.4195 m at 60 s).
flume_add_regtest(drains_approach_exponent_0 NPROCS 1
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.drains.street.froude_exponent=0
    CHECK ${check_dr} plt ${ap_check} --approach 0.4 0.302 0 --n-plotfiles 7)
# Two names with other coefficients and exponents in one run
# (Tests/Drains/inputs_approach_two), one step: each cell holds the closed
# form of its own inlet.
flume_add_regtest(drains_approach_two_names NPROCS 2
    INPUTS ${dr_dir}/inputs_approach_two
    ARGS flume.drains.west.file=${dr_dir}/inlets_west.txt
         flume.drains.east.file=${dr_dir}/inlets_east.txt
    CHECK ${check_dr} plt --approach-first-step 0.5 --velocity 1.2 0.5 --all-inlets
          --group 0.4 0.302 -0.816 0 1 --group 0.4 0.226 -0.5 2 3 --n-plotfiles 2)
# Published numbers: for the grate of their tests at full size (0.75 m by
# 0.45 m, A = 0.3375 m2) in water 0.4 m deep, Xia et al. (2022, section
# 3.4) give 0.1166 m3/s at 0.1 m/s and 0.1570 m3/s at 0.5 m/s. --published
# checks that the rate of the check script gives both to 0.5 % (it gives
# 0.1166 and 0.1567 with g = 9.81 m/s2), and the run checks one step of
# FLUME with that grate at 0.5 m/s against the closed form of that rate.
flume_add_regtest(drains_approach_published NPROCS 1
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.init.surface_elevation=0.4 "flume.init.velocity=0.5 0.0"
         flume.drains.street.area=0.3375 flume.max_step=1 flume.plot_int=1
    CHECK ${check_dr} plt --approach-first-step 0.4 --velocity 0.5 0.0 --all-inlets
          --group 0.3375 0.302 -0.816 0 3 --published --n-plotfiles 2)
# Still water: with an exponent above -1 the law takes nothing, and the
# depths stay those of the first plotfile to the last bit.
flume_add_regtest(drains_approach_still NPROCS 2
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.init.velocity=0\ 0
    CHECK ${check_dr} plt --approach 0.4 --still --all-inlets --n-plotfiles 2)
# Still water with an exponent of -1: the rate is k A sqrt(g h), without the
# speed, and the cells follow that closed form (to 0.0640 m at 60 s).
flume_add_regtest(drains_approach_exponent_minus_1 NPROCS 2
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.init.velocity=0\ 0 flume.drains.street.froude_exponent=-1
    CHECK ${check_dr} plt --uniform 0.5 --all-inlets --tol-eps 64 --approach 0.4 0.302 -1
          --n-plotfiles 5)
# The water balance of drains_approach_uniform: 400 m2 x (0.5 m - 69.4318 mm)
# drained by 60 s (the closed form), the rest stored.
flume_add_regtest(balance_drains_approach NPROCS 4
    INPUTS ${dr_dir}/inputs_approach
    ARGS ${ap_all} flume.plot_int=100000
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py
          --expect-columns stored drained
          --expect-final drained=172.22726599458983 stored=27.772734005410154 --rtol 1e-12)

# The weir and orifice law is as it was. drains_basin with its type given by
# name, with its inlet from a list file, and with the inlet as two halves in
# one cell of a list file: every value of every plotfile is that of
# drains_basin to the last bit.
flume_add_regtest(drains_basin_type NPROCS 2
    INPUTS ${dr_dir}/inputs_basin
    ARGS flume.drains.grate.type=weir_orifice
    CHECK ${check_same} --compare ${FLUME_TEST_RUN_DIR}/drains_basin)
flume_add_regtest(drains_basin_file NPROCS 2
    INPUTS ${dr_dir}/inputs_basin_list
    ARGS flume.drains.grate.file=${dr_dir}/inlets_basin.txt
    OUTPUT "drain inlets grate: 1 from [^ ]*inlets_basin\\.txt in 1 cell\\(s\\), 0 with more than one inlet: perimeter 0\\.5 m, open area 0\\.05 m2, weir coefficient 1\\.66, orifice coefficient 0\\.67; orifice above a depth of 0\\.1787"
    CHECK ${check_same} --compare ${FLUME_TEST_RUN_DIR}/drains_basin)
flume_add_regtest(drains_basin_file_halves NPROCS 2
    INPUTS ${dr_dir}/inputs_basin_list
    ARGS flume.drains.grate.file=${dr_dir}/inlets_basin_halves.txt
         flume.drains.grate.perimeter=0.25 flume.drains.grate.area=0.025
    OUTPUT "drain inlets grate: 2 from [^ ]*inlets_basin_halves\\.txt in 1 cell\\(s\\), 1 with more than one inlet"
    CHECK ${check_same} --compare ${FLUME_TEST_RUN_DIR}/drains_basin)
foreach(t drains_basin_type drains_basin_file drains_basin_file_halves)
    if(TEST ${t})
        set_tests_properties(${t} PROPERTIES FIXTURES_REQUIRED drains_basin)
    endif()
endforeach()

# One approach-flow inlet in water that sloshes between walls
# (Tests/Drains/inputs_approach_basin), 60 s: only cell (8, 8) drains, it
# drains more than 0, and the water is kept to round-off. On 4 ranks with
# the inlet from flume.drains.grate.location (it writes the reference
# checkpoints of the restart test below); on 3 ranks in boxes of 4 with the
# inlet from a list file, to the last bit of the first.
set(ap_loc flume.drains.grate.location=8.5\ 8.5)
set(ap_file flume.drains.grate.file=${dr_dir}/inlets_basin.txt)
flume_add_regtest(drains_approach_basin NPROCS 4
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS ${ap_loc} flume.check_int=200
    OUTPUT "drain inlet grate at \\(8\\.5, 8\\.5\\) m in cell \\(8, 8\\): approach flow, area 0\\.2 m2, coefficient 0\\.302, Froude exponent -0\\.816"
    CHECK ${check_dr} plt --approach 0.2 --inlet 8 8 --drained-positive --n-plotfiles 6)
if(TEST drains_approach_basin)
    set_tests_properties(drains_approach_basin PROPERTIES FIXTURES_SETUP drains_approach_basin)
endif()
flume_add_regtest(drains_approach_basin_file_np3 NPROCS 3
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS ${ap_file} flume.max_grid_size=4
    CHECK ${check_same} --compare ${FLUME_TEST_RUN_DIR}/drains_approach_basin)
if(TEST drains_approach_basin_file_np3)
    set_tests_properties(drains_approach_basin_file_np3 PROPERTIES
                         FIXTURES_REQUIRED drains_approach_basin)
endif()
# The drained volume of drains_approach_basin, 3.9504 m3 of 76.8: the value
# the run gives, kept as a regression value to 1e-6 (the case has no closed
# form; the volume changes by 4 % with the step, Docs/sphinx_doc/Verification.rst).
flume_add_regtest(balance_drains_approach_basin NPROCS 2
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS ${ap_loc} flume.plot_int=100000
    CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_balance.py
          --expect-columns stored drained
          --expect-final drained=3.95037629 stored=72.84962371 --rtol 1e-6)
# The inlet at second order.
flume_add_regtest(drains_approach_basin_order2 NPROCS 2
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS ${ap_loc} flume.order=2
    CHECK ${check_dr} plt --approach 0.2 --inlet 8 8 --drained-positive --n-plotfiles 2)
# Continued from the checkpoint of step 200 on 2 ranks, the inlet from the
# list file in both parts: the plotfiles to the last bit, the balance to
# round-off. The restart reads the same inlet from another file, saved with
# a byte-order mark, line ends of CR and LF and a tab: the name and the
# form of the file may change. The file is written here, in the build tree:
# the style checks of the repository allow no tab and no CR in its files.
string(ASCII 239 187 191 dr_bom)
set(dr_bom_list ${CMAKE_CURRENT_BINARY_DIR}/inlets_basin_bom.txt)
file(WRITE ${dr_bom_list}
     "${dr_bom}# The inlet of the basin, as an editor on Windows may save it.\r\n8.5\t8.5\r\n")
flume_add_restart_test(drains_approach_restart_np4_np2 REFERENCE drains_approach_basin
    INPUTS ${dr_dir}/inputs_approach_basin
    FIRST_NPROCS 4 RESTART_NPROCS 2
    FIRST_ARGS ${ap_file} flume.check_int=200 flume.max_step=300
    RESTART_ARGS flume.drains.grate.file=${dr_bom_list} flume.restart=chk00200
    CHECK_ARGS --gauges 0 --rewound --balance-rtol 1e-14)
# A file-list checkpoint for the restarts that must stop.
flume_add_regtest(drains_approach_basin_file NPROCS 1
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS ${ap_file} flume.check_int=200 flume.max_step=200
    CHECK ${check_dr} plt --approach 0.2 --inlet 8 8 --drained-positive --n-plotfiles 3)
if(TEST drains_approach_basin_file)
    set_tests_properties(drains_approach_basin_file PROPERTIES
                         FIXTURES_SETUP drains_approach_basin_file)
endif()
# A restart with another list (a second inlet), another law, another
# coefficient, or the inlet from a location in place of the list stops.
function(ap_restart_abort name)
    flume_add_abort_test(${name} ${ARGN})
    if(TEST ${name} AND TEST drains_approach_basin_file)
        set_tests_properties(${name} PROPERTIES FIXTURES_REQUIRED drains_approach_basin_file)
    elseif(TEST ${name})
        set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    endif()
endfunction()
set(ap_chk flume.restart=${FLUME_TEST_RUN_DIR}/drains_approach_basin_file/chk00200)
ap_restart_abort(abort_restart_drains_file_list
    EXPECT "its drain inlets [(]grate:file:1:1:[0-9a-f]+:approach_flow:[^)]*[)] differ from those of the inputs [(]grate:file:2:2:[0-9a-f]+:approach_flow:"
    ARGS ${dr_dir}/inputs_approach_basin ${ap_chk}
         flume.drains.grate.file=${dr_dir}/inlets_basin_two.txt)
# One inlet in one cell in both lists, but another cell: only the hash of
# the cells differs.
ap_restart_abort(abort_restart_drains_file_moved
    EXPECT "its drain inlets [(]grate:file:1:1:66b2e98c45d21c84:[^)]*[)] differ from those of the inputs [(]grate:file:1:1:(?!66b2e98c45d21c84)[0-9a-f]+:approach_flow:"
    ARGS ${dr_dir}/inputs_approach_basin ${ap_chk}
         flume.drains.grate.file=${dr_dir}/inlets_basin_moved.txt)
ap_restart_abort(abort_restart_drains_file_location
    EXPECT "its drain inlets [(]grate:file:1:1:[^)]*[)] differ from those of the inputs [(]grate:8:8:approach_flow:"
    ARGS ${dr_dir}/inputs_approach_basin ${ap_chk} ${ap_loc})
ap_restart_abort(abort_restart_drains_approach_coefficient
    EXPECT "its drain inlets [(]grate:file:[^)]*:0\\.30[0-9]*:-0\\.81[0-9]*[)] differ from those of the inputs [(]grate:file:[^)]*:0\\.22[0-9]*:-0\\.81[0-9]*[)]"
    ARGS ${dr_dir}/inputs_approach_basin ${ap_chk} ${ap_file} flume.drains.grate.coefficient=0.226)
ap_restart_abort(abort_restart_drains_approach_exponent
    EXPECT "its drain inlets [(]grate:file:[^)]*:-0\\.81[0-9]*[)] differ from those of the inputs [(]grate:file:[^)]*:-0\\.8[0-9]*[)]"
    ARGS ${dr_dir}/inputs_approach_basin ${ap_chk} ${ap_file}
         flume.drains.grate.froude_exponent=-0.8)
ap_restart_abort(abort_restart_drains_law
    EXPECT "its drain inlets [(]grate:file:[^)]*approach_flow[^)]*[)] differ from those of the inputs [(]grate:file:1:1:[0-9a-f]+:1:0\\.2[0-9]*:1\\.6"
    ARGS ${dr_dir}/inputs_basin_list ${ap_chk} ${ap_file} flume.plot_int=100
         flume.stop_time=60 flume.init.velocity=0.5\ 0.2
         flume.drains.grate.perimeter=1 flume.drains.grate.area=0.2)

# The same 3 inlets in the same 2 cells, but two in cell (8, 8) and one in
# cell (3, 3) in the checkpoint, and one and two in the inputs: only the
# hash differs.
flume_add_regtest(drains_approach_basin_counts NPROCS 1
    INPUTS ${dr_dir}/inputs_approach_basin
    ARGS flume.drains.grate.file=${dr_dir}/inlets_basin_2_1.txt flume.check_int=50
         flume.max_step=50 flume.plot_int=50
    OUTPUT "drain inlets grate: 3 from [^ ]*inlets_basin_2_1\\.txt in 2 cell\\(s\\), 1 with more than one inlet"
    CHECK ${check_dr} plt --approach 0.2 --inlet 8 8 --inlet 3 3 --drained-positive --n-plotfiles 2)
if(TEST drains_approach_basin_counts)
    set_tests_properties(drains_approach_basin_counts PROPERTIES
                         FIXTURES_SETUP drains_approach_basin_counts)
    flume_add_abort_test(abort_restart_drains_file_counts
        EXPECT "its drain inlets [(]grate:file:3:2:([0-9a-f]+):[^)]*[)] differ from those of the inputs [(]grate:file:3:2:(?!\\1)[0-9a-f]+:approach_flow:"
        ARGS ${dr_dir}/inputs_approach_basin
             flume.restart=${FLUME_TEST_RUN_DIR}/drains_approach_basin_counts/chk00050
             flume.drains.grate.file=${dr_dir}/inlets_basin_1_2.txt)
    set_tests_properties(abort_restart_drains_file_counts PROPERTIES
                         FIXTURES_REQUIRED drains_approach_basin_counts)
endif()

# Inputs of the approach-flow law and of the list file, on the 8 x 8 cells
# of [0, 40]^2 m (cells of 25 m2).
set(dra flume.drains.names=g flume.drains.g.type=approach_flow flume.drains.g.location=12\ 12
        flume.drains.g.area=0.4)
set(drf flume.drains.names=g flume.drains.g.type=approach_flow flume.drains.g.area=0.4)
set(bad_dr ${dr_dir}/bad)
flume_add_abort_test(abort_drains_type
    EXPECT "${drk}\\.type must be weir_orifice or approach_flow; it is grate"
    ARGS ${rok} ${dr} flume.drains.g.type=grate)
flume_add_abort_test(abort_drains_location_and_file
    EXPECT "${drk}\\.location and ${drk}\\.file are both given: an inlet takes one of the two"
    ARGS ${rok} ${dr} flume.drains.g.file=${dr_dir}/inlets_all.txt)
flume_add_abort_test(abort_drains_neither_location_nor_file
    EXPECT "flume\\.drains\\.names lists g, which needs ${drk}\\.location: the x and y of the inlet, m \\(or ${drk}\\.file: lines of an x and a y in m, one inlet each\\)"
    ARGS ${rok} ${drf})
flume_add_abort_test(abort_drains_approach_missing_area
    EXPECT "flume\\.drains\\.names lists g, which needs ${drk}\\.area: the area of the grate of the inlet, m2"
    ARGS ${rok} flume.drains.names=g flume.drains.g.type=approach_flow ${dr_loc})
flume_add_abort_test(abort_drains_approach_area_cell
    EXPECT "${drk}\\.area must be > 0 and <= the area of a cell \\(25 m2\\); it is 25\\.5"
    ARGS ${rok} ${dra} flume.drains.g.area=25.5)
flume_add_abort_test(abort_drains_coefficient_zero
    EXPECT "${drk}\\.coefficient must be > 0; it is 0"
    ARGS ${rok} ${dra} flume.drains.g.coefficient=0)
flume_add_abort_test(abort_drains_coefficient_nan
    EXPECT "${drk}\\.coefficient must be finite"
    ARGS ${rok} ${dra} flume.drains.g.coefficient=nan)
flume_add_abort_test(abort_drains_froude_exponent_below
    EXPECT "${drk}\\.froude_exponent must be >= -1 and <= 0; it is -1\\.5"
    ARGS ${rok} ${dra} flume.drains.g.froude_exponent=-1.5)
flume_add_abort_test(abort_drains_froude_exponent_above
    EXPECT "${drk}\\.froude_exponent must be >= -1 and <= 0; it is 0\\.5"
    ARGS ${rok} ${dra} flume.drains.g.froude_exponent=0.5)
# A key of the other law.
foreach(k perimeter weir_coefficient orifice_coefficient)
    flume_add_abort_test(abort_drains_approach_law_${k}
        EXPECT "${drk}\\.${k} is read only with ${drk}\\.type = weir_orifice; the type is approach_flow"
        ARGS ${rok} ${dra} flume.drains.g.${k}=0.5)
endforeach()
foreach(k coefficient froude_exponent)
    flume_add_abort_test(abort_drains_weir_law_${k}
        EXPECT "${drk}\\.${k} is read only with ${drk}\\.type = approach_flow; the type is weir_orifice"
        ARGS ${rok} ${dr} flume.drains.g.${k}=-0.5)
endforeach()
# The list file.
flume_add_abort_test(abort_drains_file_two_names
    EXPECT "${drk}\\.file needs exactly 1 value"
    ARGS ${rok} ${drf} flume.drains.g.file=a.txt\ b.txt)
flume_add_abort_test(abort_drains_file_missing
    EXPECT "${drk}\\.file \\(no_such_file\\.txt\\): the file cannot be opened \\(a relative path, taken from "
    ARGS ${rok} ${drf} flume.drains.g.file=no_such_file.txt)
flume_add_abort_test(abort_drains_file_directory
    EXPECT "${drk}\\.file \\([^)]*\\): the path is not a regular file"
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr})
flume_add_abort_test(abort_drains_file_items
    EXPECT "${drk}\\.file \\([^)]*inlets_items\\.txt\\), line 3: a line holds the x and the y of an inlet in m; it holds 3 item\\(s\\)"
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr}/inlets_items.txt)
flume_add_abort_test(abort_drains_file_number
    EXPECT "${drk}\\.file \\([^)]*inlets_number\\.txt\\), line 2: y must be a finite number; it is \"north\""
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr}/inlets_number.txt)
flume_add_abort_test(abort_drains_file_inf
    EXPECT "${drk}\\.file \\([^)]*inlets_inf\\.txt\\), line 2: y must be a finite number; it is \"inf\""
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr}/inlets_inf.txt)
flume_add_abort_test(abort_drains_file_hex
    EXPECT "${drk}\\.file \\([^)]*inlets_hex\\.txt\\), line 2: x must be a finite number; it is \"0x10\""
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr}/inlets_hex.txt)
# A coefficient so small that the rate of the cell is 0 (in single
# precision the coefficient itself reads as 0, and the range check stops).
if(NOT FLUME_PRECISION STREQUAL "SINGLE")
    flume_add_abort_test(abort_drains_rate_zero
        EXPECT "flume\\.drains\\.g: the rate of the inlet in cell \\(2, 2\\) is 0 in the precision of this build"
        ARGS ${rok} ${dra} flume.drains.g.coefficient=1e-323)
endif()
# The fourth line of the file, after an empty line and a comment; the
# domain is [0, 20]^2 m here.
flume_add_abort_test(abort_drains_file_outside
    EXPECT "${drk}\\.file \\([^)]*inlets_outside\\.txt\\), line 4: the inlet must be inside the domain: from 0 to 20 in x \\(flume\\.prob_lo and flume\\.prob_hi\\); it is 20\\.5"
    ARGS flume.n_cell=4\ 4 flume.prob_hi=20\ 20 ${drf}
         flume.drains.g.file=${bad_dr}/inlets_outside.txt)
flume_add_abort_test(abort_drains_file_empty
    EXPECT "${drk}\\.file \\([^)]*inlets_empty\\.txt\\): the file holds no inlet: it needs at least one line with an x and a y in m"
    ARGS ${rok} ${drf} flume.drains.g.file=${bad_dr}/inlets_empty.txt)
# The two inlets of a cell of 25 m2, 13 m2 each.
flume_add_abort_test(abort_drains_file_area_sum
    EXPECT "${drk}\\.file \\([^)]*inlets_twice\\.txt\\), line 4: cell \\(0, 0\\) holds 2 inlets, whose areas \\(${drk}\\.area\\) sum to 26 m2: more than the area of a cell \\(25 m2\\)"
    ARGS ${rok} ${drf} flume.drains.g.file=${dr_dir}/inlets_twice.txt flume.drains.g.area=13)
# An inlet of a list in the cell of another name, in either order.
flume_add_abort_test(abort_drains_file_same_cell
    EXPECT "flume\\.drains\\.h\\.file \\([^)]*inlets_basin\\.txt\\), line 2 puts an inlet in cell \\(1, 1\\), which holds inlet g already: only the inlets of one file share a cell"
    ARGS ${rok} ${dr} flume.drains.g.location=8\ 8 flume.drains.names=g\ h
         flume.drains.h.type=approach_flow flume.drains.h.area=0.4
         flume.drains.h.file=${dr_dir}/inlets_basin.txt)
flume_add_abort_test(abort_drains_location_in_file_cell
    EXPECT "${drk}\\.location puts inlet g in cell \\(1, 1\\), which holds inlet h already"
    ARGS ${rok} ${dr} flume.drains.g.location=8\ 8 flume.drains.names=h\ g
         flume.drains.h.type=approach_flow flume.drains.h.area=0.4
         flume.drains.h.file=${dr_dir}/inlets_basin.txt)
# A list inlet in a building cell: line 3 of the list, (12.5, 30.5) m, is
# cell (12, 30) of the streets of Tests/Raster/inputs_streets.
flume_add_abort_test(abort_drains_file_building
    EXPECT "${drk}\\.file \\([^)]*inlets_building\\.txt\\), line 3 puts an inlet in a building cell \\(flume\\.buildings\\.file\\)"
    ARGS flume.n_cell=40\ 40 flume.prob_hi=40\ 40 ${drf}
         flume.drains.g.file=${bad_dr}/inlets_building.txt
         flume.buildings.file=${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc)
