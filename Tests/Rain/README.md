# Prescribed rain tests

## Purpose

Check the rain that FLUME takes from `flume.rain.type = constant` (a rate in
mm/h) and `flume.rain.type = hyetograph` (a rate in time from a file, held
or linear between entries): that a closed domain receives exactly the rain
that falls, at every plotfile time; that a dry domain under rain gets a
finite time step; that rain on an impermeable plane follows the
kinematic-wave solution; and that every malformed input or file stops the
run with a message naming the key, the file and the line. The model is in
`Docs/sphinx_doc/Model.rst` (Rain, Time step) and the inputs in
`Docs/sphinx_doc/Inputs.rst` (Hyetograph files).

## Files

| File | What it is |
| --- | --- |
| `check_rain.py` | Reads every plotfile of a run and checks it against the rain fallen by its time (computed here from the rate or the file): uniform depth and zero discharge (`--uniform`), the volume added to a closed domain (`--added-volume`), or the kinematic-wave solution on a plane (`--kinematic-wave`). |
| `inputs_flat` | 24 x 16 cells of 2 m in 6 boxes, walls, flat dry bed; 36 mm/h to 600 s. |
| `inputs_flat_steps` | The same domain with `hyetograph_steps.txt`, held between entries, to 1800 s. |
| `inputs_flat_triangle` | The same domain with `hyetograph_triangle.txt`, linear between entries, to 1800 s. |
| `inputs_flat_deep` | The same domain with 1 m of still water and 0.36 mm/h, to 600 s. |
| `inputs_plane` | 1000 x 2 cells of 0.1 m by 1 m, bed slope 0.01 towards an outflow side, n = 0.01, 50 mm/h to 900 s. |
| `hyetograph_steps.txt` | 36 mm/h from 0 s, 72 mm/h from 600 s, 0 from 1200 s; has comments. |
| `hyetograph_triangle.txt` | 0, 72 and 0 mm/h at 0, 600 and 1200 s; has a blank line. |
| `hyetograph_late.txt` | 36 mm/h from 600 s, 0 from 1200 s, with no line at 0 s. |
| `hyetograph_burst.txt` | 0 mm/h, 3600 mm/h from 100 s to 101 s, then 0. |
| `bad/*.txt` | One malformed hyetograph per check of the reader. |
| `RainTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The closed-basin test uses `Tests/Basin/inputs`, and the tests without a stop
time use `Tests/Skeleton/inputs`.

## Data and acknowledgements

Every rate and file is synthetic, written for these tests. No external data
are used or redistributed, and nothing needs to be downloaded.

## How to run

From the build directory:

```bash
ctest -R rain --output-on-failure
```

By hand, from an empty directory (`flume.rain.file` takes an absolute path,
or a path relative to the directory the run starts in):

```bash
<build>/flume <repo>/Tests/Rain/inputs_flat_steps \
    flume.rain.file=<repo>/Tests/Rain/hyetograph_steps.txt
python3 <repo>/Tests/Rain/check_rain.py flat_steps \
    --hyetograph <repo>/Tests/Rain/hyetograph_steps.txt --uniform
```

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/Rain/inputs_plane
python3 <repo>/Tests/Rain/check_rain.py plane --rate 50 \
    --kinematic-wave 0.01 0.01 100 --depth-l1-max 0.06 --q-out-max 0.065 \
    --equilibrium-rtol 1e-9
```

## Expected Results

| Test | Result |
| --- | --- |
| `rain_flat_constant_np1`, `_np2` | 6 mm in every cell at 600 s, discharges 0 |
| `rain_flat_hyetograph_step` | 18 mm at 1800 s; every plotfile equal to the rain by its time |
| `rain_flat_hyetograph_linear` | 12 mm at 1800 s; largest difference 1.7e-18 m |
| `rain_flat_hyetograph_step_file_linear` | 15 mm at 1800 s |
| `rain_flat_hyetograph_late_start` | 12 mm at 1800 s; the run prints a note on the first time |
| `rain_flat_deep_water` | 1 m plus 0.06 mm at 600 s, within 4 epsilons of the depth |
| `rain_flat_no_stop_time` | three steps on a dry domain without a stop time |
| `rain_dry_burst_no_stop_time` | 1 mm after the burst; finite steps before it |
| `rain_dry_start_discharge` | discharges exactly 0 after every step from a moving dry start |
| `rain_dry_tiny_rate` | two steps at 1e-300 mm/h; the step iteration converges |
| `rain_basin_volume` | added depth sum 6.4e-3, difference 8.9e-13 |
| `rain_plane_kinematic_wave` | depth L1 5.5 %, outlet 5.7 % of i L; equilibrium outflow i L to 1.2e-13 |
| `rain_plane_kinematic_wave_order2` | at cfl 0.125: depth L1 0.22 %, outlet 1.3 % of i L |

The numbers are what the runs give in double precision, not the bounds the
checks use (`RainTests.cmake`, `Docs/sphinx_doc/Testing.rst`).

## CTest

`rain_flat_constant_np1`, `rain_flat_constant_np2`,
`rain_flat_hyetograph_step`, `rain_flat_hyetograph_linear`,
`rain_flat_hyetograph_step_file_linear`, `rain_flat_hyetograph_late_start`,
`rain_flat_deep_water`, `rain_flat_no_stop_time`,
`rain_dry_burst_no_stop_time`, `rain_dry_start_discharge`,
`rain_dry_tiny_rate`, `rain_basin_volume`,
`rain_plane_kinematic_wave`, `rain_plane_kinematic_wave_order2`, and the
`abort_rain_*` tests. They run on 1 or 2 ranks. See
`Docs/sphinx_doc/Testing.rst`.
