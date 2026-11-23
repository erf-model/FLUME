# Infiltration tests

## Purpose

Check the two infiltration models of `flume.infiltration.type`: Green-Ampt
(`green_ampt`) and the runoff curve number as a loss on the rain
(`curve_number`). Ponded Green-Ampt infiltration and Green-Ampt under rain
are checked against the closed forms, the curve-number loss against its
formula, and the water balance of a domain with flow, a conductivity map
and buildings to round-off. Further tests check that a cell never loses more
water than it holds, that buildings and impervious cells lose nothing, that
water which infiltrates takes its share of the discharge, and that every
malformed input or file stops the run with a message naming the key and the
file. The models are in `Docs/sphinx_doc/theory/SourceTerms.rst`
(Infiltration) and the inputs in `Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `check_infiltration.py` | Reads every plotfile of a run and checks it: Green-Ampt against the closed forms solved by bisection (`--green-ampt`), the curve-number loss (`--curve-number`, with `--cn-map`), the water balance of a closed domain (`--balance`), and a film that infiltrates in one step (`--cap`); `--impervious` and `--impervious-zero` check that the cells a raster marks lose nothing; `--ledger-eps` checks each cell's depth plus infiltrated depth against the water it got, and `--steady-tol` keeps the bound from growing with the steps. |
| `make_infiltration_grids.py` | Writes the `.asc` files below; running it again writes them byte for byte. |
| `inputs_ponded` | 8 x 8 cells of 40 m in 4 boxes, closed and flat; 0.5 m of still water on silt loam, to 6 h. |
| `inputs_rain_green_ampt` | 8 x 8 cells of 10 m, closed, flat and dry; 36 mm/h on sandy loam, to 2 h. |
| `inputs_curve_number` | The same domain; 36 mm/h with CN 75, to 2 h. |
| `inputs_cn_map` | 40 x 40 cells of 1 m in 9 boxes, walls, bed sloping towards x = 0 and y = 0, n = 0.03; 36 mm/h with a curve-number map, to 30 min. |
| `inputs_balance` | The same domain at second order with a hyetograph and a conductivity map, to 30 min. |
| `hyetograph_burst.txt` | 3600 mm/h for the first 1 s, then dry. |
| `conductivity_columns.asc` | 8 x 8 cells of 5 m: columns of K = 0, 1.016, 6.604, 10.922, 29.972, 120.396, 0 and 3.302 mm/h. |
| `curve_number_columns.asc` | 8 x 8 cells of 5 m: columns of CN = 100, 98, 90, 75, 65.31, 50, 34.74 and 85. |
| `bad/*.asc` | 2 x 2 cells of 20 m, one quadrant out of range: a deficit of 1.5, a curve number of 0 and of 101. |
| `InfiltrationTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The map and balance tests take the buildings of
`Tests/Raster/streets_heights.asc` and the hyetograph
`Tests/Rain/hyetograph_steps.txt`; the abort tests also use
`Tests/Raster/bad/negative.asc` and `nodata.asc`.

## Data and acknowledgements

No data are downloaded, and nothing needs to be. The grids and decks are
synthetic, written for these tests. The soil parameters are the
representative values of Table 4-7 of L. A. Rossman and W. C. Huber (2016),
*Storm Water Management Model Reference Manual, Volume I: Hydrology
(Revised)*, EPA/600/R-15/162A, U.S. Environmental Protection Agency, a
public document. That table gives the values of W. J. Rawls,
D. L. Brakensiek and N. Miller (1983), Green-Ampt infiltration parameters
from soils data, *Journal of Hydraulic Engineering* 109(1), 62-70. They are
converted here from inches to millimetres:

| Soil | K, in/h | psi, in | K, mm/h | psi, mm |
| --- | --- | --- | --- | --- |
| Clay loam | 0.04 | 8.22 | 1.016 | 208.788 |
| Silt loam | 0.26 | 6.57 | 6.604 | 166.878 |
| Sandy loam | 0.43 | 4.33 | 10.922 | 109.982 |
| Loamy sand | 1.18 | 2.41 | 29.972 | 61.214 |
| Sand | 4.74 | 1.95 | 120.396 | 49.53 |
| Loam | 0.13 | 3.50 | 3.302 | 88.9 |

The tests use the conductivities of every class and the suction heads of
silt loam and sandy loam.

The curve-number method follows the NRCS *National Engineering Handbook,
Part 630 Hydrology*, chapter 10 (July 2004), a public document. The
curve numbers 65.31 and 34.74 are the lambda = 0.05 equivalents of 75 and
50 in Table 4 of D. E. Woodward, R. H. Hawkins, R. Jiang,
A. T. Hjelmfelt, J. A. Van Mullem and Q. D. Quan (2003), Runoff curve number
method: examination of the initial abstraction ratio, World Water &
Environmental Resources Congress 2003, ASCE,
https://doi.org/10.1061/40685(2003)308. We thank the U.S. Environmental
Protection Agency and the USDA Natural Resources Conservation Service for
making their manuals freely available.

## How to run

From the build directory:

```bash
ctest -R infiltration --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/Infiltration/inputs_ponded
python3 <repo>/Tests/Infiltration/check_infiltration.py ponded \
    --green-ampt 6.604 166.878 0.3 --initial-depth 0.5 --time-shift 1e-5
```

```bash
mpiexec -n 3 <build>/flume <repo>/Tests/Infiltration/inputs_balance \
    flume.rain.file=<repo>/Tests/Rain/hyetograph_steps.txt \
    flume.infiltration.green_ampt.conductivity_file=<repo>/Tests/Infiltration/conductivity_columns.asc \
    flume.buildings.file=<repo>/Tests/Raster/streets_heights.asc
python3 <repo>/Tests/Infiltration/check_infiltration.py bal --balance \
    --hyetograph <repo>/Tests/Rain/hyetograph_steps.txt \
    --impervious <repo>/Tests/Raster/streets_heights.asc \
    --impervious-zero <repo>/Tests/Infiltration/conductivity_columns.asc
```

## Expected Results

| Test | Result |
| --- | --- |
| `infiltration_green_ampt_ponded` | F = 91.75 mm at 6 h; largest difference 2.9e-11 m (a time shift of 3.7e-6 s) |
| `infiltration_green_ampt_ponded_moving` | the same F, velocity (0.5, 0.25) m/s kept; largest difference 2.2e-11 m |
| `infiltration_green_ampt_rain` | F = 50.60 mm and 21.40 mm on the surface at 2 h; largest difference 7.6e-17 m |
| `infiltration_curve_number_rain` | loss 50.30 mm of 72 mm at 2 h; largest difference 1.4e-17 m |
| `infiltration_curve_number_ia_ratio` | CN 65.31 with Ia = 0.05 S; largest difference 1.4e-17 m |
| `infiltration_curve_number_many_steps` | CN 75 on 2 cm cells, 279217 steps; largest difference 1.4e-17 m (4.5e-9 m in single precision, 1.2e-5 m without the compensated rain total) |
| `infiltration_curve_number_map` | every cell its own loss, buildings and CN 100 none; largest difference 3.5e-18 m |
| `infiltration_green_ampt_balance` | 18.04 m infiltrated and 10.76 m on the surface (sums over the cells) of 28.8 m of rain; difference 5.7e-14 |
| `infiltration_green_ampt_cap`, `_cap_moving` | 1 mm infiltrated exactly, no water or discharge left |
| `infiltration_green_ampt_burst` | every cell emptied in the first step stays at 0 in the second; without the fix to the rain carry, about 400 of the 1600 cells go below 0 |

In the cases with a closed form, each cell's depth plus infiltrated depth is
the water it got to 4 epsilons (`--ledger-eps 4`). The numbers are what the
runs give in double precision, not the bounds the
checks use (`InfiltrationTests.cmake`, `Docs/sphinx_doc/Verification.rst`).

## CTest

`infiltration_green_ampt_ponded`, `infiltration_green_ampt_ponded_moving`,
`infiltration_green_ampt_rain`, `infiltration_curve_number_rain`,
`infiltration_curve_number_ia_ratio`, `infiltration_curve_number_many_steps`,
`infiltration_curve_number_map`, `infiltration_green_ampt_balance`,
`infiltration_green_ampt_cap`, `infiltration_green_ampt_cap_moving`,
`infiltration_green_ampt_burst`, and the `abort_infiltration_*` tests.
They run on 1 to 4 ranks.
