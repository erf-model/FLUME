# Rain from ERF's microphysics

## Purpose

Check FLUME inside ERF with `flume.rain.type = microphysics`: the rain of
each flood cell is the surface precipitation ERF's microphysics gives its
atmospheric column, liquid and frozen, as liquid water (1 kg/m2 is 1 mm).
The checks are that each cell gets its own column's precipitation, over
FLUME steps shorter than ERF's; that the rain of the water balance is the
precipitation over the grid when the water moves; that the frozen part is
counted; that a restart from ERF's checkpoint gives the same flood to the
last bit; and that ERF stops a run whose microphysics has no precipitation.
The frozen case checks the bookkeeping only: FLUME has no snowpack, and its
ground below freezing is where adding snow as water at once holds least.
The checks are bounded for double precision; a single-precision build fails
them.
The coupling is described in `Docs/sphinx_doc/Running.rst` (Inside ERF).

## Files

| File | What it is |
| --- | --- |
| `inputs_kessler` | ERF: Kessler warm rain from a supersaturated layer below 150 m in a laterally periodic box, 1600 x 1600 x 1024 m in 16 x 16 x 32 cells (boxes of 8 x 16 x 32), a cool bubble off the centre that makes the rain vary over the columns, `erf.fixed_dt = 0.5` (cut to about 0.34 s by ERF's CFL limit), 120 steps, a checkpoint at step 60, 2D plotfiles with `precip_total_accum` and `precip_frozen_accum`, `erf.flood.grid_ratio = 2`. |
| `sounding_moist` | The sounding of `inputs_kessler`. |
| `inputs_wdm6` | ERF: a cold, shallow WDM6 moist bubble over ground below freezing, 20 x 0.4 x 5 km in 200 x 4 x 50 cells (two boxes), 200 steps, `erf.flood.grid_ratio = 1`. |
| `flume_cells` | FLUME's half of the cell-by-cell runs: `flume.rain.type = microphysics` on a flat bed, `flume.dry_depth = 1` so that no water moves, `flume.max_dt = 0.2` so that each ERF step is two FLUME steps. |
| `flume_basin` | FLUME's half of the volume run: `flume.rain.type = microphysics` on a dam break of 0.1 m over a flat periodic bed (surface 0.1 m for x < 800 m, dry beyond), Manning n = 0.03. |
| `check_rain.py` | `cells`: each flood cell holds its column's `precip_total_accum` / 1000, with the spread of the precipitation over the columns and its frozen part bounded below. `volume`: the rain of the water balance is the precipitation summed over the grid, and the water moved; with `--restart`, ERF continued from its checkpoint gives the same flood plotfiles and water balance. |
| `ERFRainTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt` with `FLUME_WITH_ERF=ON`. |

`check_rain.py` reads the plotfiles with `Tests/check_plotfile.py` and
`Tests/ERF/check_coupled.py`, and the water balance with
`Tests/check_balance.py`.

## Data and acknowledgements

No external data are used or redistributed, and nothing needs to be
downloaded. The atmospheres are those of two ERF regression tests, from the
ERF repository (`Submodules/ERF`, BSD-3-Clause-LBNL licence), with thanks to
the ERF developers:

- `inputs_kessler` and `sounding_moist` are ERF's
  `Tests/FloodEngineExample/inputs_flood_example` and `sounding_moist`, with
  FLUME as the engine and the plot and checkpoint keys these tests need.
- `inputs_wdm6` is ERF's
  `Tests/test_files/WDM6_ColdBubble_SnowAccum/WDM6_ColdBubble_SnowAccum.i`,
  with two boxes, 2D plotfiles and FLUME as the flood engine.

## How to run

From the build directory of a `-DFLUME_WITH_ERF=ON` build:

```bash
ctest -R "erf_flume_rain|abort_erf_flume_rain" --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 2 <build>/ERF/Exec/erf_exec <repo>/Tests/ERFRain/inputs_kessler \
    FILE=<repo>/Tests/ERFRain/flume_cells \
    erf.input_sounding_file=<repo>/Tests/ERFRain/sounding_moist erf.check_int=-1
python3 <repo>/Tests/ERFRain/check_rain.py cells 2 60 120 --min-spread 1.5
```

## Expected Results

| Test | Result |
| --- | --- |
| `erf_flume_rain_cells` | 120 ERF steps to 40.5 s, 240 FLUME steps; 1000 x depth within 1.1e-19 kg/m2 of its column's precipitation in every cell (largest 7.6e-4 kg/m2); spread 3.3 over the columns; water-balance residual 1.1e-16 m3 |
| `erf_flume_rain_volume` | rain 0.6346442873734 m3 against 0.6346442873737 m3 of precipitation over the grid (difference 3.7e-13 m3: the rain still held in the cells' rounding carry, bound 5.8e-11 m3); the water moved by up to 0.028 m in a cell; after the restart from `chk00060`, the flood plotfile of step 120 and the water balance identical |
| `erf_flume_rain_frozen` | 200 ERF steps to 72.8 s, 400 FLUME steps; 1000 x depth within 8.7e-19 kg/m2 of its column's precipitation in every cell (largest 6.3e-3 kg/m2); frozen part 1.0 of the total |
| `abort_erf_flume_rain_no_moisture` | ERF stops at the start: the flood engine takes the precipitation, but `erf.moisture_model = "None"` provides none |
| `abort_erf_flume_rain_rate` | `flume.rain.rate` with `microphysics` inside ERF stops the run, naming the key |

Each test fails against these one-line changes to the coupling:

| Change | Test that fails |
| --- | --- |
| Only the liquid part | `erf_flume_rain_frozen` |
| Only the frozen part | `erf_flume_rain_cells` |
| kg/m2 taken as m | `erf_flume_rain_cells` |
| Each cell given the column two cells over | `erf_flume_rain_cells` |
| The whole ERF step's precipitation on every FLUME step | `erf_flume_rain_cells` |
| The rate spread over twice the ERF step | `erf_flume_rain_cells`, `erf_flume_rain_volume` |
| The engine does not ask ERF for the precipitation | `erf_flume_rain_cells`, `abort_erf_flume_rain_no_moisture` |
| `flume.rain.rate` accepted with `microphysics` | `abort_erf_flume_rain_rate` |

## CTest

`erf_flume_rain_cells`, `erf_flume_rain_volume`, `erf_flume_rain_frozen`
(2 ranks each), `abort_erf_flume_rain_no_moisture` and
`abort_erf_flume_rain_rate` (1 rank each). The
stand-alone checks of the key are `abort_rain_type_microphysics` and
`abort_rain_type_unknown` (`Tests/Rain`). See
`Docs/sphinx_doc/TestCatalog.rst`.
