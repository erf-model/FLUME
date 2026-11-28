# Storm-drain inlet tests

## Purpose

Check the storm-drain inlets of `flume.drains.*`, which send water out of
the model, under their two capacity laws: `weir_orifice`, grates in a sag
that take the smaller of the weir and the orifice capacity at the depth of
their cell, and `approach_flow`, Q = k u A Fr^alpha from the speed u and
the Froude number of the water of the cell. Inlets are given one by one
(`location`) or as a list file (`file`). The weir and orifice tests come
first. Level, still water over inlets in
every cell is checked against the exact solution of the depth equation, from
a start above the depth where the two rates are equal and from one below it.
One inlet in a closed basin is checked against the closed form for one cell
of the basin's area, which is an upper bound. A 1 mm film over an inlet far
too large for it checks that no depth goes below 0. Further tests pin the
water balance, a restart on another rank count, the independence of the
decomposition, and that every malformed input stops the run with a message
naming the key. For the approach-flow law, moving water over an inlet in
every periodic cell is checked against the exact solution of its depth
equation, up to and past the time the cells empty, and still water checks
that nothing drains unless the exponent is -1. Old-against-new runs check
that the weir and orifice law gives the same bits with the type named and
with the inlet from a list file, and that inlets of a list that share a
cell are added. The model is in `Docs/sphinx_doc/theory/SourceTerms.rst`
(Storm-drain inlets) and the inputs in `Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `check_drains.py` | Reads every plotfile of a run and checks that no depth is below 0, that `drained_depth` is 0 away from the inlets (`--inlet I J`, or `--all-inlets`), and that the depth plus the drained depth summed over the cells stays that of the first plotfile (`--ledger-eps`). `--uniform H0` checks every cell against the exact solution, through its inverse (the time to reach a depth), and `--basin H0 LO HI` the drained volume against the lumped closed form, found by bisection on the same inverse, and `--first-step H0` the depth each inlet cell drains in the first step (`--approach-first-step H0` with one `--group A k alpha ILO IHI` per name for approach-flow inlets with other parameters). `--grate P A [C_w C_o]` gives a weir and orifice inlet; `--approach A [k alpha]` an approach-flow inlet in water of `--velocity UX UY`, for which `--uniform` also checks the discharge and the emptying time. `--still` checks that nothing drains or moves, `--drained-positive` that each inlet cell drained. |
| `inputs_uniform` | 2 x 2 cells of 10 m in 4 boxes, closed and flat, still water 0.8 m deep, an inlet of P = 2 m and A = 0.5 m2 in every cell, to 600 s. The switch depth is 0.4469 m. |
| `inputs_basin` | 16 x 16 cells of 1 m in 4 boxes, closed and flat, still water 0.3 m deep, one inlet of P = 0.5 m and A = 0.05 m2 in cell (8, 8), to 1200 s. The switch depth is 0.1788 m. |
| `inputs_film` | 4 x 4 cells of 1 m in 4 boxes, closed and flat, a 1 mm film, one inlet of P = 100 m and A = 1 m2 in cell (1, 2), 20 steps. |
| `inputs_basin_list` | The basin of `inputs_basin` without the location of its inlet: the tests give a list file. |
| `inputs_approach` | 4 x 4 cells of 5 m in 4 boxes, periodic and flat, water 0.5 m deep at (1.2, 0.5) m/s, approach-flow inlets of 0.4 m2 from a list file given by the tests, to 60 s. The cells empty at 87.05 s. |
| `inputs_approach_two` | The cells of `inputs_approach`, one step, two names: `west` with the default coefficient and exponent, `east` with 0.226 and -0.5; the tests give the two lists. |
| `inputs_approach_basin` | 16 x 16 cells of 1 m in 4 boxes, closed and flat, water 0.3 m deep that starts at (0.5, 0.2) m/s, one approach-flow inlet of 0.2 m2 placed by the tests in cell (8, 8), to 60 s. |
| `inlets_all.txt`, `inlets_twice.txt`, `inlets_west.txt`, `inlets_east.txt` | Lists for `inputs_approach` and `inputs_approach_two`: one inlet at the centre of each cell; two inlets in each cell; the western and the eastern 8 cells. |
| `inlets_basin.txt`, `inlets_basin_halves.txt`, `inlets_basin_two.txt`, `inlets_basin_moved.txt` | Lists for the basins: the inlet of cell (8, 8); two inlets in that cell; that inlet and one in cell (3, 3); one inlet in cell (3, 3). |
| `inlets_basin_2_1.txt`, `inlets_basin_1_2.txt` | Two inlets in cell (8, 8) and one in cell (3, 3), and the reverse. |
| `bad/` | List files that must stop the run: a line of three items, a word, a hexadecimal number, an infinite number, a point outside the domain, no inlet, an inlet in a building cell. |
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

The approach-flow law and its default coefficients (k = 0.302,
alpha = -0.816) are those of J. Xia, B. Dong, M. Zhou, R. Ahmadian,
R. A. Falconer, Q. Li and X. Zhang (2022), A unified formula for discharge
capacity of street inlets for urban flood management, *Journal of
Hydrology* 609, 127667, https://doi.org/10.1016/j.jhydrol.2022.127667, in
the form given by S. Annis, M. G. Badas and G. Mascaro (2026), Increasing
the fidelity of hyperlocal simulations of urban pluvial flooding through
street flooding observations, *Advances in Water Resources* 209, 105223,
https://doi.org/10.1016/j.advwatres.2026.105223 (equation 3). The
coefficient 0.226 of `drains_approach_coefficients` is the value Annis et
al. calibrated for Brooklyn. We thank both groups of authors. Only the
formula, these three numbers and the two discharges Xia et al. print for
their grate (0.1166 and 0.1570 m3/s, checked by `drains_approach_published`)
are used; no data, figure or text of either paper is in the repository.

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

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/Drains/inputs_approach flume.drains.street.file=<repo>/Tests/Drains/inlets_all.txt
python3 <repo>/Tests/Drains/check_drains.py plt --uniform 0.5 --approach 0.4 --velocity 1.2 0.5 --all-inlets --tol-eps 64
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
| `drains_approach_uniform` | to 69.43 mm at 60 s in every cell; largest difference 1.3 epsilons of 0.5 m; the discharge is the depth times (1.2, 0.5) m/s |
| `drains_approach_coefficients` | k = 0.226 and alpha = -0.5, to 171.3 mm at 60 s; largest difference 2.9 epsilons of 0.5 m |
| `drains_approach_merge` | two inlets of 0.2 m2 per cell: every value of every plotfile that of `drains_approach_uniform` to the last bit |
| `drains_approach_empty` | to 120 s in 145 steps, a plotfile at each: wet in every plotfile until the closed form reaches the dry depth, empty in every one after 87.05 s, 200 m3 drained, no depth below 0; largest difference 1.6 epsilons |
| `drains_approach_empty_small_steps` | the same in 6492 steps: a step lands where the closed form is below the dry depth, the step takes that film and the cells end empty; largest difference 11.6 epsilons of 0.5 m |
| `drains_approach_exponent_0` | alpha = 0, to 419.5 mm at 60 s; largest difference 0.6 epsilons of 0.5 m |
| `drains_approach_published` | Xia et al. (2022) give 0.1166 and 0.1570 m3/s for their grate (0.3375 m2) in water 0.4 m deep at 0.1 and 0.5 m/s; the rate of the check script gives 0.1166 and 0.1567, and one step of FLUME with that grate follows its closed form |
| `drains_approach_two_names` | after one step the western cells hold the depth of k = 0.302, alpha = -0.816 and the eastern cells that of k = 0.226, alpha = -0.5; largest difference 0.2 epsilons of 0.5 m |
| `drains_approach_still` | still water, default exponent: nothing drains, every depth that of the start to the last bit |
| `drains_approach_exponent_minus_1` | still water, alpha = -1: to 64.05 mm at 60 s; largest difference 0.7 epsilons of 0.5 m |
| `balance_drains_approach` | 172.2272660 m3 drained and 27.7727340 m3 stored, the closed form to 1e-12 |
| `drains_basin_type`, `_file`, `_file_halves` | every value of every plotfile that of `drains_basin` to the last bit |
| `drains_approach_basin` | 3.950 m3 of 76.8 drained in 60 s (456 steps), all from cell (8, 8); residual 9e-15 m3 |
| `balance_drains_approach_basin` | drained 3.95037629 m3 and stored 72.84962371 m3 to 1e-6: a regression value of this grid and step, not a closed form |
| `drains_approach_basin_order2` | 4.310 m3 drained at second order |
| `drains_approach_basin_file_np3` | inlet from a list file, 3 ranks, 16 boxes: every value that of `drains_approach_basin` to the last bit |
| `drains_approach_restart_np4_np2` | continued from step 200 on 2 ranks, the list read from a file CMake writes in the build tree, with a byte-order mark, line ends of CR and LF and a tab: the plotfiles to the last bit, the balance within 1e-14 |
| `abort_restart_drains_file_list`, `_file_moved`, `_file_location`, `_approach_coefficient`, `_approach_exponent`, `_law` | a restart from `drains_approach_basin_file` with a second inlet, the inlet in another cell, a location in place of the list, k = 0.226, alpha = -0.8, or the weir and orifice law, stops |
| `drains_approach_basin_counts`, `abort_restart_drains_file_counts` | 50 steps with two inlets in cell (8, 8) and one in cell (3, 3); a restart with one and two stops: only the hash of the list differs |

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
`abort_restart_drains_moved`, `abort_restart_drains_coefficient`,
`drains_approach_uniform`, `drains_approach_merge`,
`drains_approach_coefficients`, `drains_approach_empty`,
`drains_approach_still`, `drains_approach_exponent_minus_1`,
`balance_drains_approach`, `drains_basin_type`, `drains_basin_file`,
`drains_basin_file_halves`, `drains_approach_basin`,
`drains_approach_basin_file_np3`, `drains_approach_basin_order2`,
`drains_approach_restart_np4_np2`, `drains_approach_basin_file`, the five
`abort_restart_drains_*` tests of the approach-flow basin, and the
`abort_drains_*` tests. They run on 1 to 4 ranks.
