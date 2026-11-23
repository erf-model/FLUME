# Storm-drain inlet tests

## Purpose

Check the storm-drain inlets of `flume.drains.*`: grates in a sag that take
the smaller of the weir and the orifice capacity at the depth of their cell,
and send the water out of the model. Level, still water over inlets in
every cell is checked against the exact solution of the depth equation, from
a start above the depth where the two rates are equal and from one below it.
One inlet in a closed basin is checked against the closed form for one cell
of the basin's area, which is an upper bound. A 1 mm film over an inlet far
too large for it checks that no depth goes below 0. Further tests pin the
water balance, a restart on another rank count, the independence of the
decomposition, and that every malformed input stops the run with a message
naming the key. The model is in `Docs/sphinx_doc/theory/SourceTerms.rst`
(Storm-drain inlets) and the inputs in `Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `check_drains.py` | Reads every plotfile of a run and checks that no depth is below 0, that `drained_depth` is 0 away from the inlets (`--inlet I J`, or `--all-inlets`), and that the depth plus the drained depth summed over the cells stays that of the first plotfile (`--ledger-eps`). `--uniform H0` checks every cell against the exact solution, through its inverse (the time to reach a depth), and `--basin H0 LO HI` the drained volume against the lumped closed form, found by bisection on the same inverse, and `--first-step H0` the depth each inlet cell drains in the first step. `--grate P A [C_w C_o]` gives the inlet. |
| `inputs_uniform` | 2 x 2 cells of 10 m in 4 boxes, closed and flat, still water 0.8 m deep, an inlet of P = 2 m and A = 0.5 m2 in every cell, to 600 s. The switch depth is 0.4469 m. |
| `inputs_basin` | 16 x 16 cells of 1 m in 4 boxes, closed and flat, still water 0.3 m deep, one inlet of P = 0.5 m and A = 0.05 m2 in cell (8, 8), to 1200 s. The switch depth is 0.1788 m. |
| `inputs_film` | 4 x 4 cells of 1 m in 4 boxes, closed and flat, a 1 mm film, one inlet of P = 100 m and A = 1 m2 in cell (1, 2), 20 steps. |
| `DrainsTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The abort test of an inlet in a building cell takes the buildings of
`Tests/Raster/streets_heights.asc`.

## Data and acknowledgements

No data are downloaded, and nothing needs to be. The decks are synthetic,
written for these tests. The capacities and the default coefficients
(C_w = 1.66 in SI units, C_o = 0.67) are those of grates in a sag in
S. A. Brown, J. D. Schall, J. L. Morris, C. L. Doherty, S. M. Stein and
J. C. Warner (2009), *Urban Drainage Design Manual*, Hydraulic Engineering
Circular No. 22, third edition, FHWA-NHI-10-009, U.S. Department of
Transportation, Federal Highway Administration, a public document. We thank
the Federal Highway Administration for making it freely available.

## How to run

From the build directory:

```bash
ctest -R drains --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/Drains/inputs_uniform
python3 <repo>/Tests/Drains/check_drains.py uni --uniform 0.8 --grate 2 0.5 --all-inlets
```

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/Drains/inputs_basin
python3 <repo>/Tests/Drains/check_drains.py plt --basin 0.3 0.95 1.0 --grate 0.5 0.05 --inlet 8 8
```

## Expected Results

| Test | Result |
| --- | --- |
| `drains_uniform` | orifice flow to 0.4469 m at 30.4 s, then weir flow, to 8.34 mm at 600 s; largest difference 0.8 epsilons of 0.8 m |
| `drains_uniform_weir` | weir flow only from 0.3 m, to 7.20 mm at 600 s; largest difference 1.0 epsilon of 0.3 m |
| `drains_uniform_coefficients` | C_w = 1.5 and C_o = 0.6 for every inlet, switch at 0.4429 m; largest difference 3.4 epsilons of 0.8 m |
| `balance_drains_uniform` | 316.6641978 m3 drained and 3.3358022 m3 stored, the closed form to 1e-15 |
| `drains_basin` | drained volume 0.979 to 0.987 of the lumped closed form (bounds 0.95 to 1: the upper one holds while the inlet cell is the lowest; the lower one is set from the run) |
| `drains_basin_order2` | 0.976 to 0.984 of the lumped closed form at second order |
| `drains_basin_boxes_4_np3` | every value of every plotfile that of `drains_basin` to the last bit |
| `drains_restart_np4_np2` | continued from step 2000 on 2 ranks: the plotfiles to the last bit, the balance within 1e-14 |
| `drains_film` | smallest depth 1.2e-5 m, none below 0; an explicit step at the weir rate would take 13 mm from the 1 mm film; the depth drained in the first step is the exact one to 0.2 epsilons of 1 mm |
| `abort_restart_drains_moved`, `_coefficient` | a restart from `drains_basin` with the inlet in cell (3, 3), or a weir coefficient of 1.5, stops |

In every plotfile, the depth plus the drained depth summed over the cells
is that of the first plotfile to 8 epsilons; in the basin runs, where the
water moves, to 8 + sqrt(n) epsilons after n steps (1 epsilon in double
precision and 17 in single precision over 6095 steps). The numbers are what the runs
give in double precision, not the bounds the checks use
(`DrainsTests.cmake`, `Docs/sphinx_doc/Verification.rst`). The basin drains
less than the lumped form because the inlet cell is drawn below the mean
depth of the basin.

## CTest

`drains_uniform`, `drains_uniform_weir`, `drains_uniform_coefficients`,
`balance_drains_uniform`, `drains_basin`, `drains_basin_order2`,
`drains_basin_boxes_4_np3`, `drains_restart_np4_np2`, `drains_film`,
`abort_restart_drains_moved`, `abort_restart_drains_coefficient`, and the
`abort_drains_*` tests. They run on 1 to 4 ranks.
