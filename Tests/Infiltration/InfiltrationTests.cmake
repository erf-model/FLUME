# Tests of infiltration (flume.infiltration.*), included from
# Tests/CMakeLists.txt, whose functions (flume_add_regtest,
# flume_add_abort_test) and ${rok} and ${bad} decks they use.
# CMAKE_CURRENT_SOURCE_DIR is Tests/. Tests/Infiltration/README.md describes
# the cases; Tests/Infiltration/check_infiltration.py computes the expected
# numbers from the closed forms, independently of the solver.
set(check_inf ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Infiltration/check_infiltration.py)
set(inf_dir ${CMAKE_CURRENT_SOURCE_DIR}/Infiltration)

# Ponded Green-Ampt infiltration (Tests/Infiltration/inputs_ponded): 0.5 m
# of still water on silt loam, 6 h, every 250 steps. Every cell holds the
# cumulative infiltration of the implicit equation K t = F - M ln(1 + F / M),
# solved by bisection, and the depth 0.5 m - F. The first step infiltrates
# the ponded water as rain of 0.5 m / dt and saturates the surface after
# K M dt^2 / H0^2 = 7.5e-6 s instead of at once; exact ponding takes half
# that, so F lags by K M dt^2 / (2 H0^2) = 3.7e-6 s (2.9e-11 m). The
# integrated equation is exact from there. The bound is 1e-5 s on top of
# round-off. In every cell the depth plus infiltrated_depth is 0.5 m to 4
# epsilons.
flume_add_regtest(infiltration_green_ampt_ponded NPROCS 2
    INPUTS ${inf_dir}/inputs_ponded
    OUTPUT "Green-Ampt infiltration: conductivity 6\\.604 mm/h, suction 166\\.878 mm, deficit 0\\.3"
    CHECK ${check_inf} ponded --green-ampt 6.604 166.878 0.3 --initial-depth 0.5
          --time-shift 1e-5 --ledger-eps 4 --n-plotfiles 10)
# The same water moving at (0.5, 0.25) m/s on a periodic domain, with no
# friction: the state stays uniform, and the water that infiltrates takes
# its share of the discharge, so the velocity stays (0.5, 0.25) m/s.
flume_add_regtest(infiltration_green_ampt_ponded_moving NPROCS 2
    INPUTS ${inf_dir}/inputs_ponded
    ARGS flume.periodic=1\ 1 flume.init.velocity=0.5\ 0.25
    CHECK ${check_inf} ponded --green-ampt 6.604 166.878 0.3 --initial-depth 0.5
          --time-shift 1e-5 --velocity 0.5 0.25 --ledger-eps 4 --n-plotfiles 10)
# Green-Ampt under constant rain (inputs_rain_green_ampt): 36 mm/h on dry
# sandy loam. All the rain infiltrates up to F_s = 14.37 mm at 1437 s; then
# the surface is saturated and water ponds. The step that crosses F_s is
# split there, so the run follows the closed form to round-off (8e-17 m).
flume_add_regtest(infiltration_green_ampt_rain NPROCS 4
    INPUTS ${inf_dir}/inputs_rain_green_ampt
    CHECK ${check_inf} rga --green-ampt 10.922 109.982 0.3 --rain 36 --ledger-eps 4 --n-plotfiles 8)
# Curve number 75 under constant rain (inputs_curve_number): no water stays
# on the surface until P = Ia = 16.93 mm at 1693 s; at 72 mm the cumulative
# loss is P - (P - Ia)^2 / (P - Ia + S) = 50.30 mm.
flume_add_regtest(infiltration_curve_number_rain NPROCS 2
    INPUTS ${inf_dir}/inputs_curve_number
    OUTPUT "curve-number infiltration: CN 75, Ia / S = 0\\.2"
    CHECK ${check_inf} cn --curve-number 75 --rain 36 --flat --ledger-eps 4 --n-plotfiles 7)
# The same with Ia = 0.05 S and CN 65.31, the lambda = 0.05 equivalent of
# CN 75 (Woodward et al., 2003, Table 4).
flume_add_regtest(infiltration_curve_number_ia_ratio NPROCS 1
    INPUTS ${inf_dir}/inputs_curve_number
    ARGS flume.infiltration.cn.number=65.31 flume.infiltration.cn.ia_ratio=0.05
    OUTPUT "curve-number infiltration: CN 65\\.31, Ia / S = 0\\.05"
    CHECK ${check_inf} cn --curve-number 65.31 --ia-ratio 0.05 --rain 36 --flat --ledger-eps 4
          --n-plotfiles 7)
# The same rain and soil on cells of 2 cm, so 279217 steps of 2.6e-7 m of rain
# on average. The rain total and the depth each cell lost stay within a few
# epsilons, with a tolerance that does not grow with the steps: the rain
# total is a compensated sum, and the depth lost is the depth removed.
flume_add_regtest(infiltration_curve_number_many_steps NPROCS 1
    INPUTS ${inf_dir}/inputs_curve_number
    ARGS flume.plot_file=cnmany flume.n_cell=4\ 4 flume.prob_hi=0.08\ 0.08
         flume.plot_int=100000
    CHECK ${check_inf} cnmany --curve-number 75 --rain 36 --flat --steady-tol --ledger-eps 4
          --n-plotfiles 4)
# Curve numbers from a map, with flow and buildings (inputs_cn_map): water
# runs between cells, yet each cell loses the curve-number loss of its own
# rain, and the 256 building cells and the CN = 100 column lose nothing.
flume_add_regtest(infiltration_curve_number_map NPROCS 4
    INPUTS ${inf_dir}/inputs_cn_map
    ARGS flume.infiltration.cn.number_file=${inf_dir}/curve_number_columns.asc
         flume.buildings.file=${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc
    OUTPUT "CN from 34\\.74 to 100 \\([^)]*curve_number_columns\\.asc\\), Ia / S = 0\\.2; 256 building cell\\(s\\) impervious"
    CHECK ${check_inf} cnm --curve-number 0 --cn-map ${inf_dir}/curve_number_columns.asc
          --rain 36 --impervious ${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc
          --n-plotfiles 5)
# The water balance (inputs_balance): a sloped, walled domain, second order,
# the hyetograph of Tests/Rain/hyetograph_steps.txt, a conductivity map and
# buildings, on 3 ranks (9 boxes). In every plotfile the water on the
# surface plus the water infiltrated is the rain that fell, to round-off;
# the building cells and the K = 0 columns infiltrate nothing.
flume_add_regtest(infiltration_green_ampt_balance NPROCS 3
    INPUTS ${inf_dir}/inputs_balance
    ARGS flume.rain.file=${CMAKE_CURRENT_SOURCE_DIR}/Rain/hyetograph_steps.txt
         flume.infiltration.green_ampt.conductivity_file=${inf_dir}/conductivity_columns.asc
         flume.buildings.file=${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc
    OUTPUT "conductivity from 0 to 120\\.396 \\([^)]*conductivity_columns\\.asc\\) mm/h"
    CHECK ${check_inf} bal --balance --hyetograph ${CMAKE_CURRENT_SOURCE_DIR}/Rain/hyetograph_steps.txt
          --impervious ${CMAKE_CURRENT_SOURCE_DIR}/Raster/streets_heights.asc
          --impervious-zero ${inf_dir}/conductivity_columns.asc --n-plotfiles 5)
# A cell never loses more water than it holds: 1 mm of water on a soil that
# could take 28 m in the first step. After it every cell is dry, holds no
# discharge, and has infiltrated exactly its 1 mm. Two lines keep the loss
# within the water: the branch for a rate of water below K, and the cap on
# the loss (a guard against rounding, since the ponded rate never exceeds
# the rate of water). The test fails when both are removed; each alone
# covers for the other.
flume_add_regtest(infiltration_green_ampt_cap NPROCS 1
    INPUTS ${inf_dir}/inputs_ponded
    ARGS flume.plot_file=cap flume.init.surface_elevation=0.001
         flume.infiltration.green_ampt.conductivity=1e6 flume.max_step=3 flume.plot_int=1
    CHECK ${check_inf} cap --cap --n-plotfiles 3)

# The same film moving at (1, -2) m/s on a periodic domain (so it stays
# uniform): a cell left dry holds no discharge.
flume_add_regtest(infiltration_green_ampt_cap_moving NPROCS 2
    INPUTS ${inf_dir}/inputs_ponded
    ARGS flume.plot_file=cap flume.init.surface_elevation=0.001 flume.init.velocity=1\ -2
         flume.periodic=1\ 1
         flume.infiltration.green_ampt.conductivity=1e6 flume.max_step=3 flume.plot_int=1
    CHECK ${check_inf} cap --cap --n-plotfiles 3)

# Rain that rounds, on a cell that infiltration then empties
# (Tests/Infiltration/hyetograph_burst.txt): 1.5 mm of water on a tilted
# surface (so each cell holds a different depth) on a soil that takes
# everything, and rain during the first 1 s of the first step only. Adding
# the rain to the water rounds, and the rain source keeps what was rounded
# away for the next step; in about half the cells that carry is negative.
# Infiltration then empties every cell. In the second step there is no rain:
# the carry waits instead of taking the depth below 0 (the run would abort in
# some 400 of the 1600 cells), and the water balances.
flume_add_regtest(infiltration_green_ampt_burst NPROCS 2
    INPUTS ${inf_dir}/inputs_ponded
    ARGS flume.plot_file=burst flume.plot_int=1 flume.max_step=3 flume.n_cell=40\ 40
         flume.prob_hi=40\ 40 flume.max_grid_size=14 flume.init.surface_elevation=0.0015
         flume.init.surface_slope=0.0000123\ 0.0000071 flume.rain.type=hyetograph
         flume.rain.file=${inf_dir}/hyetograph_burst.txt
         flume.infiltration.green_ampt.conductivity=1e6
    CHECK ${check_inf} burst --balance --hyetograph ${inf_dir}/hyetograph_burst.txt
          --n-plotfiles 3)

# Infiltration inputs, on the 8 x 8 cells of [0, 40]^2 m of the raster checks.
set(ga flume.infiltration.type=green_ampt flume.infiltration.green_ampt.conductivity=10
       flume.infiltration.green_ampt.suction=100 flume.infiltration.green_ampt.deficit=0.3)
set(cn flume.infiltration.type=curve_number flume.infiltration.cn.number=75
       flume.rain.type=constant flume.rain.rate=10)
set(gak "flume\\.infiltration\\.green_ampt")
set(cnk "flume\\.infiltration\\.cn")
flume_add_abort_test(abort_infiltration_type_unknown
    EXPECT "flume\\.infiltration\\.type must be none, green_ampt or curve_number; it is horton"
    ARGS ${rok} flume.infiltration.type=horton)
foreach(case
        "conductivity|-1|>= 0 \\(mm/h\\)"
        "suction|-1|>= 0 \\(mm\\)"
        "deficit|-0.1|>= 0 and <= 1"
        "deficit|1.5|>= 0 and <= 1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 value)
    list(GET c 2 range)
    flume_add_abort_test(abort_infiltration_green_ampt_${name}_${value}
        EXPECT "${gak}\\.${name} must be ${range}; it is ${value}"
        ARGS ${rok} ${ga} flume.infiltration.green_ampt.${name}=${value})
endforeach()
# Each parameter is required with its model.
foreach(name conductivity suction deficit)
    set(deck ${ga})
    list(FILTER deck EXCLUDE REGEX "green_ampt\\.${name}=")
    flume_add_abort_test(abort_infiltration_green_ampt_${name}_missing
        EXPECT "flume\\.infiltration\\.type = green_ampt needs ${gak}\\.${name} or ${gak}\\.${name}_file"
        ARGS ${rok} ${deck})
endforeach()
foreach(case "number|0" "number|101" "ia_ratio|1" "ia_ratio|-0.2")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 value)
    if(name STREQUAL "number")
        set(range "> 0 and <= 100")
    else()
        set(range ">= 0 and < 1")
    endif()
    flume_add_abort_test(abort_infiltration_cn_${name}_${value}
        EXPECT "${cnk}\\.${name} must be ${range}; it is ${value}"
        ARGS ${rok} ${cn} flume.infiltration.cn.${name}=${value})
endforeach()
flume_add_abort_test(abort_infiltration_conductivity_inf
    EXPECT "${gak}\\.conductivity must be finite"
    ARGS ${rok} ${ga} flume.infiltration.green_ampt.conductivity=inf)
flume_add_abort_test(abort_infiltration_cn_number_missing
    EXPECT "flume\\.infiltration\\.type = curve_number needs ${cnk}\\.number or ${cnk}\\.number_file"
    ARGS ${rok} flume.infiltration.type=curve_number flume.rain.type=constant flume.rain.rate=10)
flume_add_abort_test(abort_infiltration_cn_without_rain
    EXPECT "flume\\.infiltration\\.type = curve_number needs rain \\(flume\\.rain\\.type\\)"
    ARGS ${rok} flume.infiltration.type=curve_number flume.infiltration.cn.number=75)
flume_add_abort_test(abort_infiltration_value_with_file
    EXPECT "${gak}\\.conductivity_file cannot be used with ${gak}\\.conductivity: the file gives the value of every cell"
    ARGS ${rok} ${ga} flume.infiltration.green_ampt.conductivity_file=${inf_dir}/conductivity_columns.asc)
flume_add_abort_test(abort_infiltration_green_ampt_key_without_type
    EXPECT "${gak}\\.suction is read only with flume\\.infiltration\\.type = green_ampt"
    ARGS ${rok} flume.infiltration.green_ampt.suction=100)
flume_add_abort_test(abort_infiltration_green_ampt_key_with_cn
    EXPECT "${gak}\\.conductivity_file is read only with flume\\.infiltration\\.type = green_ampt"
    ARGS ${rok} ${cn} flume.infiltration.green_ampt.conductivity_file=${inf_dir}/conductivity_columns.asc)
flume_add_abort_test(abort_infiltration_cn_key_with_green_ampt
    EXPECT "${cnk}\\.number is read only with flume\\.infiltration\\.type = curve_number"
    ARGS ${rok} ${ga} flume.infiltration.cn.number=75)
flume_add_abort_test(abort_infiltration_ia_ratio_with_green_ampt
    EXPECT "${cnk}\\.ia_ratio is read only with flume\\.infiltration\\.type = curve_number"
    ARGS ${rok} ${ga} flume.infiltration.cn.ia_ratio=0.05)
flume_add_abort_test(abort_infiltration_file_empty
    EXPECT "${gak}\\.deficit_file must not be empty: the name of an ESRI ASCII grid file"
    ARGS ${rok} flume.infiltration.type=green_ampt flume.infiltration.green_ampt.conductivity=10
         flume.infiltration.green_ampt.suction=100 "flume.infiltration.green_ampt.deficit_file=\"\"")
# Files whose values are out of range: one quadrant of 4 x 4 cells.
set(inf_bad ${inf_dir}/bad)
foreach(case
        "green_ampt.deficit|${inf_bad}/deficit_above_one.asc|>= 0 and <= 1"
        "green_ampt.conductivity|${bad}/negative.asc|>= 0 \\(mm/h\\)"
        "cn.number|${inf_bad}/curve_number_zero.asc|> 0 and <= 100"
        "cn.number|${inf_bad}/curve_number_above_100.asc|> 0 and <= 100")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 key)
    list(GET c 1 file)
    list(GET c 2 range)
    string(REPLACE "." "\\." key_re "${key}")
    get_filename_component(stem "${file}" NAME_WE)
    string(REPLACE "." "_" key_name "${key}")
    if(key MATCHES "^cn")
        set(deck ${cn})
    else()
        set(deck ${ga})
    endif()
    list(FILTER deck EXCLUDE REGEX "^flume\\.infiltration\\.${key_re}=")
    flume_add_abort_test(abort_infiltration_${key_name}_file_${stem}
        EXPECT "flume\\.infiltration\\.${key_re}_file \\([^)]*${stem}\\.asc\\): 16 cell\\(s\\) of the domain take a value that is not ${range}"
        ARGS ${rok} ${deck} flume.infiltration.${key}_file=${file})
endforeach()
flume_add_abort_test(abort_infiltration_file_nodata
    EXPECT "flume\\.infiltration\\.green_ampt\\.suction_file \\([^)]*\\): 16 cell\\(s\\) of the domain read a NODATA_value cell"
    ARGS ${rok} flume.infiltration.type=green_ampt flume.infiltration.green_ampt.conductivity=10
         flume.infiltration.green_ampt.deficit=0.3 flume.infiltration.green_ampt.suction_file=${bad}/nodata.asc)
flume_add_abort_test(abort_infiltration_unknown_key
    EXPECT "Unknown key\\(s\\) in the flume namespace:[^a-z]*flume\\.infiltration\\.cn\\.value"
    ARGS ${rok} ${cn} flume.infiltration.cn.value=75)
