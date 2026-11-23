# Water-balance, flood-map and gauge tests

## Purpose

Check the water-balance log (`water_balance.txt`, written every
`flume.log_int` steps), the flood maps (`flume.flood_maps = 1`) and the
gauge time series (`flume.gauges.*`). The log holds the stored volume, the
water each source has added or removed, the water that has left through
each outflow segment, and the residual of the balance. The tests pin those
columns against closed forms: rain on a closed basin, ponded Green-Ampt
infiltration, and uniform flow out of a side. A further test runs rain,
infiltration and outflow together and checks that the log fails the
balance check as soon as any one term is left out. The flood maps and
gauges are checked on the same rain basin and on the wet-bed dam break; the
gauges also against the plotfile cells on a sloped bed, and with
coordinates far from the origin. A run with the log off must give the same
plotfiles, to the last byte, as one with the log on.
Every malformed input stops the run with a message naming the key.

Every other regression run that leaves a log also has its residual checked
by `Tests/run_case.py` (with `Tests/check_balance.py`). The log, the maps
and the gauge files are described in `Docs/sphinx_doc/Running.rst`, the
balance and its bound in `Docs/sphinx_doc/theory/WaterBalance.rst`, and the
inputs in `Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `inputs_rain_gauges` | The domain and rain of `Tests/Rain/inputs_flat` (24 x 16 cells of 2 m in 6 boxes, walls, 36 mm/h to 600 s) on a bed at 5 m, with the flood maps (arrival depth 3 mm) and three gauges: on a face between two cells, on the high corner and on the low edge. |
| `inputs_outflow` | Water 1 m deep at 10 m/s on a flat, frictionless bed, 300 x 20 cells of 1 m in boxes of at most 16 cells, to 2 s. The outflow segments are given on the command line. |
| `check_flood_outputs.py` | Reads the last plotfile of a run and the gauge files. `--rain-basin` checks the maps and gauges of the rain basin against the rain fallen; `--dam-break` checks the maps and a gauge against the exact Riemann solution of the dam break. |
| `check_dropped_terms.py` | Checks that a log passes `check_balance.py`, then writes the log again without each added or removed term in turn and checks that each of those fails on at least 90 % of the pairs of lines where that term changes. |
| `check_gauges.py` | Checks that each gauge read its cell and holds exactly that cell's values in the last plotfile. |
| `check_no_log.py` | Checks that a run with `flume.log_int = 0` left no log or gauge file, and that its plotfiles are those of a reference run with the log (`Tests/compare_plotfiles.py`). |
| `WaterBalanceTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The tests also run `Tests/Rain/inputs_flat`, `Tests/Rain/inputs_plane`,
`Tests/Infiltration/inputs_ponded`, `Tests/Uniform/inputs`,
`Tests/DamBreak/inputs` and `Tests/Skeleton/inputs`, with the changes the
CMake file gives on the command line.

## Data and acknowledgements

No data are downloaded, and nothing needs to be. The decks are synthetic,
written for these tests. The ponded Green-Ampt test uses the silt-loam
parameters of Table 4-7 of L. A. Rossman and W. C. Huber (2016), *Storm
Water Management Model Reference Manual, Volume I: Hydrology (Revised)*,
EPA/600/R-15/162A, U.S. Environmental Protection Agency, a public document,
as `Tests/Infiltration/README.md` sets out. We thank the U.S. Environmental
Protection Agency for making the manual freely available.

## How to run

From the build directory:

```bash
ctest -R "balance_|flood_maps|gauge_|abort_(log|flood|arrival|gauge)" --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/WaterBalance/inputs_rain_gauges
python3 <repo>/Tests/check_balance.py water_balance.txt
python3 <repo>/Tests/WaterBalance/check_flood_outputs.py rg --rain-basin 36 0.003 --bed 5 \
    --gauge face 12 3 --gauge corner 23 15 --gauge inside 1 0
```

```bash
mpiexec -n 1 <build>/flume <repo>/Tests/WaterBalance/inputs_outflow \
    flume.bc.segments=east flume.bc.east.side=xhi flume.bc.east.type=outflow
python3 <repo>/Tests/check_balance.py water_balance.txt \
    --expect-final out_east=400 stored=5600 --rtol 1e-13
```

The segment inputs above are those of `balance_outflow_side`; the CMake
file has the exact command line of each test.

## Expected Results

| Test | Result |
| --- | --- |
| `balance_rain_gauges_np1`, `_np2` | 384 cells hold a maximum depth of 6 mm and speeds of 0; arrival at 301.1 s, one step after the 300 s of the closed form; the three gauges read the cells (12, 3), (23, 15) and (1, 0) and the rain fallen at every line |
| `balance_rain_columns_np1`, `_np2` | rain 9.216 m3, all of it stored; difference 7.1e-15 m3 |
| `gauge_face_rounding` | the gauge at (0.3, 0.7) m on 0.1 m cells reads cell (3, 7), not (2, 6) |
| `gauge_cells_sloped_bed` | gauges a, b, c read cells (5, 2), (35, 22) and (47, 31) and hold their depth, surface elevation and velocity exactly |
| `gauge_large_coordinates` | the gauge at x = 500003.7 m on 1 m cells from 500000 m reads cell 3, also in single precision |
| `balance_green_ampt_ponded` | 9395.3755 m3 infiltrated in 6 h; 1.1e-6 m3 from the closed form (the 3.7e-6 s lag of the first step) |
| `balance_outflow_side` | 400 m3 out of the side at x = 300 m, 5600 m3 stored; difference 4.5e-13 m3 |
| `balance_outflow_two_segments` | 100 m3 and 300 m3 out of segments of 5 m and 15 m, second order |
| `balance_outflow_low_side` | 400 m3 out of the side at y = 0, second order |
| `balance_periodic` | 12288 m3 stored at every line, residual exactly 0 |
| `balance_log_interval` | lines at steps 0, 5, 10, 15 and 16 |
| `balance_log_on`, `balance_log_off` | 401 lines with the log; none, no summary line and the same plotfiles without it |
| `balance_every_term` | rain, infiltration and outflow together; residual -5.6e-16 m3 after 5315 steps; the log without any one term fails on every pair of lines where it changes |
| `balance_every_term_log_interval` | the same with a line every 25 steps (214 lines) |
| `flood_maps_dam_break` | 68 cells behind the shock: arrival within 0.024 s of (x - 50 m) / S (0.060 s allowed), maximum depth within 0.40 %, the speed maps within 1.5 %; 108 cells ahead hold -1 |
| `flood_maps_dam_break_second_order` | the same within 0.014 s, 0.027 % and 0.11 % |

The residual change between two lines stays within 64 machine epsilons of
the water the steps hold and move (`Tests/check_balance.py`). The largest
ratio in these tests is 12.3, in `balance_log_on`; over all regression
tests it is 14.5 in double and 15.3 in single precision
(`Docs/sphinx_doc/Verification.rst`). The numbers in the table are what the
runs give in double precision.

## CTest

`balance_rain_gauges_np1`, `balance_rain_gauges_np2`,
`balance_rain_columns_np1`, `balance_rain_columns_np2`,
`gauge_face_rounding`, `gauge_cells_sloped_bed`, `gauge_large_coordinates`,
`balance_green_ampt_ponded`, `balance_outflow_side`,
`balance_outflow_two_segments`, `balance_outflow_low_side`,
`balance_periodic`, `balance_log_interval`, `balance_log_on`,
`balance_log_off`, `balance_every_term`, `balance_every_term_log_interval`,
`flood_maps_dam_break`, `flood_maps_dam_break_second_order`, and the abort
tests `abort_log_int_negative`, `abort_flood_maps_two`,
`abort_arrival_depth_negative`, `abort_arrival_depth_below_dry`,
`abort_arrival_depth_nan`,
`abort_arrival_depth_without_maps`, `abort_gauge_name`,
`abort_gauge_twice`, `abort_gauge_case`, `abort_gauge_no_location`, `abort_gauge_outside_x`,
`abort_gauge_outside_y`, `abort_gauge_location_one_value`,
`abort_gauge_unlisted` and `abort_gauges_without_log`. The
`balance_rain_*` tests run on 1 and on 2 ranks, the others on 1 or 2 (the
CMake file gives each), and the abort tests on 1.
