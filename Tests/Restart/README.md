# Checkpoint and restart tests

## Purpose

Check that a run continued from a checkpoint (`flume.restart`) gives what
the same run gives without a break. A reference run goes to the end and
writes checkpoints (`flume.check_int`). A restart test runs the first part
of the run, which writes a checkpoint and may go on past it, then continues
from the checkpoint in the same directory. Every plotfile must be that of
the reference to the last bit, every variable of every cell, and the gauge
files byte for byte. `water_balance.txt` must be the same byte for byte on
the same rank count and boxes; on another rank count or other boxes its
values must agree to round-off, because the outflow and rain totals are
sums over the ranks. The restart is tested on the same and on another rank
count, with other boxes, from a checkpoint the first part went past (the
logs are cut back, and a line cut short at their end is removed), and from
the checkpoint of its last step, on and off the log interval.

The reference runs also check that the results do not depend on the
decomposition: the same run on 1 and 2 ranks, and in boxes of 16 and of 8
cells, gives the same plotfiles to the last bit.

Every checkpoint that does not match the inputs, and every damaged one the
tests make, stops the run with a message naming `flume.restart`. The checkpoint and the restart are
described in `Docs/sphinx_doc/Running.rst`, the inputs in
`Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `inputs` | A rain cell of 1000 mm/h moving over a tilted, rough bed with a hollow, 48 x 32 cells of 2 m in 6 boxes, second order, to 240 s: curve-number infiltration (CN 95), Manning friction, an outflow segment on the high-x side, the flood maps and two gauges. Every quantity a checkpoint holds changes during the run. |
| `inputs_first_order` | A dam break with constant rain and Green-Ampt infiltration on the same grid, first order, walls on every side, to 60 s; no maps and no gauges. |
| `check_restart.py` | Runs the restart in the directory of the first part and compares its plotfiles, gauge files and log with those of the reference (`--reference`, `--restart`); checks that checkpoints were written (`--expect`); or compares the plotfiles of two runs cell by cell (`--compare`). |
| `make_bad_checkpoints.py` | Copies a checkpoint and alters each copy in one way (version, precision, a key given twice, a negative step, a NaN in the state), for the abort tests. |
| `RestartTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

## Data and acknowledgements

No data are downloaded, and nothing needs to be. Both decks are synthetic,
written for these tests.

## How to run

From the build directory:

```bash
ctest -R "restart|check_int|check_file" --output-on-failure
```

By hand, from an empty directory (the reference in `ref/`):

```bash
mkdir ref && (cd ref && mpiexec -n 1 <build>/flume <repo>/Tests/Restart/inputs flume.check_int=400)
mpiexec -n 1 <build>/flume <repo>/Tests/Restart/inputs flume.check_int=400 flume.max_step=550
python3 <repo>/Tests/Restart/check_restart.py --reference ref --gauges 2 --balance-rtol 1e-14 --rewound \
    --restart mpiexec -n 2 <build>/flume <repo>/Tests/Restart/inputs flume.restart=chk00400
```

These are the runs of `restart_np1_np2`; the CMake file has the exact
command line of each test.

## Expected Results

| Test | Result |
| --- | --- |
| `restart_reference_np1` | 856 steps to 240 s; checkpoints `chk00400`, `chk00800` and `chk00856`; 10 plotfiles; at the end 72.11 m3 stored, 129.54 m3 of rain, 48.95 m3 infiltrated, 8.48 m3 out, residual 3.9e-14 m3 |
| `restart_reference_np2`, `restart_boxes_8_np1` | the 10 plotfiles of `restart_reference_np1`, 15360 cells, every value the same |
| `restart_first_order_reference` | 785 steps to 60 s on 2 ranks; checkpoints `chk00300`, `chk00600` and `chk00785`; at the end 7287.92 m3 stored, 10.24 m3 of rain, 19.75 m3 infiltrated, residual -3.6e-13 m3 |
| `restart_first_order_no_log` | its plotfiles on 1 rank in boxes of 8, without the log |
| `restart_np1_np1` | the first part logs to step 550, and a line cut short is added to each log; the restart from step 400 gives the plotfiles, the two gauge files and the log of the reference, byte for byte |
| `restart_np2_np2` | the same from the checkpoint the first part writes at its last step, 400 |
| `restart_np1_np1_off_interval` | the same from the checkpoint of the last step, 405, which is not a multiple of the log interval |
| `restart_np1_np2`, `restart_np2_np1_boxes_8`, `restart_np2_np2_boxes_8` | the plotfiles and gauge files to the bit; the log to 4.4e-16 of the largest value of a column (1e-14 allowed) |
| `restart_first_order_np2` | the first part logs to step 450; the restart from step 300 gives the plotfiles and the log byte for byte |

The numbers in the table are what the runs give in double precision.

## CTest

`restart_reference_np1`, `restart_reference_np2`, `restart_boxes_8_np1`,
`restart_first_order_reference`, `restart_first_order_no_log`,
`restart_np1_np1`, `restart_np2_np2`, `restart_np1_np2`,
`restart_np1_np1_off_interval`, `restart_np2_np1_boxes_8`,
`restart_np2_np2_boxes_8`, `restart_first_order_np2`, the fixtures
`restart_bad_checkpoints` and `restart_bad_checkpoints_write`, and the abort
tests `abort_check_int_negative`, `abort_check_int_two`,
`abort_check_file_empty`, `abort_check_file_without_int`,
`abort_check_file_plot_file`, `abort_restart_empty`, `abort_restart_two`,
`abort_restart_missing`, `abort_restart_plotfile`, `abort_restart_grid`,
`abort_restart_domain`, `abort_restart_bed`, `abort_restart_sources`,
`abort_restart_periodic`, `abort_restart_columns`, `abort_restart_maps`,
`abort_restart_segments`, `abort_restart_gauges`, `abort_restart_no_balance`,
`abort_restart_version`, `abort_restart_precision`, `abort_restart_duplicate`,
`abort_restart_step` and `abort_restart_nan`. The tests run on 1 or
2 ranks (the CMake file gives each); the abort tests run on 1. The abort
tests that read a checkpoint need the reference run that wrote it.
