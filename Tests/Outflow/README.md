# Outflow side tests

## Purpose

Check that a boundary segment of type `outflow`
(`flume.bc.<name>.type = outflow`) lets water leave and lets none in, and
keep the plain copy it replaced reachable as type `zero_gradient`. An
`outflow` segment copies the boundary cell into its ghost cells where the
discharge of that cell is zero or points out of the domain, and mirrors the
cells inside, as a wall, where it points in. A `zero_gradient` segment
copies whatever the direction of the flow, and nothing bounds the water that
then enters. Each case runs with both types: the `outflow` run must take no
water in, and the `zero_gradient` run must take in the amount found before
the change, so the tests fail if either type turns into the other. The model
is in `Docs/sphinx_doc/theory/Boundaries.rst` (Outflow, Zero gradient) and
the inputs in `Docs/sphinx_doc/Inputs.rst` (Boundaries).

## Files

| File | What it is |
| --- | --- |
| `check_outflow.py` | Reads the water-balance log. `--no-inflow`: no `out_<name>` column goes below zero or falls from one line to the next, and the stored volume is at most the volume at step 0 plus the water added, all within 64 eps of the largest volume of the log. `--takes-in FACTOR`: an `out_<name>` column ends below zero and the stored volume is at least FACTOR times the volume held and added. `--expect-stored V`: the stored volume of the last line is V. `--min-out V`: the `out_<name>` columns of the last line add up to at least V, which a wall would not give. |
| `make_beds.py` | Writes the two bed rasters below. |
| `bed_step.asc` | 50 x 2 cells of 2 m: bed 0 m in the first column, 0.1 m in the others. |
| `bed_valley.asc` | 100 x 2 cells of 1 m: bed 0.01 \|x - 50 m\|, lowest in the middle. |
| `inputs_step` | The step raster, segment `w` on xlo, walls elsewhere, n = 0.015, still water with its surface at 0.6 m moving at 1e-4 m/s into the domain, to 300 s. The tests give the type of `w` and the raster. |
| `inputs_valley` | The valley raster, segments `west` (xlo) and `east` (xhi), dry, rain of 50 mm/h, to 1200 s. |
| `inputs_pond` | 100 x 50 cells of 2 m in boxes of at most 32, a plane at 2 m falling 1 cm per metre in x with 0.05 m of noise, n = 0.015, wall on xlo and xhi, segments `s` (ylo) and `n` (yhi), dry, rain of 50 mm/h, to 1200 s. |
| `OutflowTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

## Data and acknowledgements

Both rasters are synthetic, written by `make_beds.py` for these tests. No
external data are used or redistributed, and nothing needs to be downloaded.

## How to run

From the build directory:

```bash
ctest -R "^outflow_(step|valley|pond)" --output-on-failure
```

By hand, from an empty directory (`flume.bed.file` takes an absolute path,
or a path relative to the directory the run starts in):

```bash
mpiexec -n 1 <build>/flume <repo>/Tests/Outflow/inputs_step \
    flume.bed.file=<repo>/Tests/Outflow/bed_step.asc flume.bc.w.type=outflow
python3 <repo>/Tests/Outflow/check_outflow.py --no-inflow
```

```bash
mpiexec -n 3 <build>/flume <repo>/Tests/Outflow/inputs_pond \
    flume.bc.s.type=zero_gradient flume.bc.n.type=zero_gradient
python3 <repo>/Tests/Outflow/check_outflow.py --takes-in 5
```

To write the rasters again:

```bash
python3 <repo>/Tests/Outflow/make_beds.py
```

## Expected Results

Volumes of the last line of `water_balance.txt`, double precision. The
tests also pass in single precision; the seed of the step case, 1e-4 m/s,
is large enough for the growth to start there too (1e-6 m/s is not).

| Test | Result |
| --- | --- |
| `outflow_step_no_inflow_order1` | 200.8 m3 at step 0; at 300 s 53.76 m3 stored and 147.04 m3 out of `w` (at least 100 required); the column never falls |
| `outflow_step_no_inflow_order2` | 57.70 m3 stored, 143.10 m3 out of `w` |
| `outflow_step_zero_gradient_takes_in` | 8973.77 m3 stored, -8772.97 m3 out of `w` (44.7 times the 200.8 m3) |
| `outflow_valley_keeps_rain_order1`, `_order2` | 3.3333333333333321 m3 stored of 10/3 m3 of rain; out of `west` and `east` 0 to 2e-23 m3 in size |
| `outflow_valley_zero_gradient_takes_in` | 3.3373 m3 stored of 10/3 m3 of rain; -2.0e-3 m3 out of each side (1.0012 times the rain; 1.001 required) |
| `outflow_pond_no_inflow_order1` | 291.56 m3 stored of 333.33 m3 of rain; 19.85 m3 out of `s`, 21.92 m3 out of `n` (at least 30 in all required) |
| `outflow_pond_no_inflow_order2` | 290.39 m3 stored; 20.47 m3 out of `s`, 22.47 m3 out of `n` |
| `outflow_pond_zero_gradient_takes_in` | 8069 m3 stored of 333.33 m3 of rain; -12322 m3 out of `s`, 4586 m3 out of `n` (24 times the rain) |

In the step case the water moves inward at first, so the `outflow` side acts
as a wall; the disturbance later turns the flow outward and the pond drains
through the side (3 m3 remain at 3600 s). In the `zero_gradient` runs the growth starts
from the seed (step) or from round-off (pond), so the numbers of the pond
run depend on the boxes and the ranks; the tests ask for a factor, not a
value. Run to 3600 s, the `zero_gradient` pond stores 1.1e7 m3 from 1000 m3
of rain.

## CTest

`outflow_step_no_inflow_order1`, `outflow_step_no_inflow_order2`,
`outflow_step_zero_gradient_takes_in`, `outflow_valley_keeps_rain_order1`,
`outflow_valley_keeps_rain_order2`, `outflow_valley_zero_gradient_takes_in`,
`outflow_pond_no_inflow_order1`,
`outflow_pond_no_inflow_order2`, `outflow_pond_zero_gradient_takes_in`.

Each `no_inflow` and `keeps_rain` test fails with the ghost fill before the
change (the plain copy), with the direction test reversed, with the sign of
the high side dropped, with the test read from a tangential component, and
with a mirror that does not negate the normal discharge (at least one test
of the group fails for each). The step and pond tests also fail if an
outflow side is always a wall, since they require water to leave.
