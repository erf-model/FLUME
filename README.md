# FLUME: FLood and Urban inundation ModulE

FLUME is a two-dimensional shallow-water flood solver for the
[Energy Research and Forecasting (ERF)](https://github.com/erf-model/ERF)
atmospheric code, built on [AMReX](https://github.com/AMReX-Codes/amrex).
It is planned as a fine 2D surface grid that shares the horizontal
decomposition of the atmosphere and is advanced inside ERF's time loop, the way
the fire front is in [FLARE](https://github.com/erf-model/FLARE). Rainfall-driven
flooding comes first. River inflow and coastal surge come later.

The expansion of the name is a draft.

## Status

The repository holds the build system, a test harness, the continuous
integration and the documentation. The executable, `flume`, reads an inputs
file, builds the grid, the bed and the initial state, and advances the
shallow-water equations with a first-order finite-volume scheme (HLLC fluxes
with hydrostatic reconstruction) between walls or periodic sides, writing
plotfiles. It is tested on a lake at rest over a rough bed, on the dam-break
solution of Stoker and on mass conservation in a closed basin.

## Planned scope

- Full shallow-water equations, finite volume, HLLC flux with hydrostatic
  reconstruction; wetting and drying; first and second order.
- Bed and free surface on an absolute vertical datum.
- Manning friction, infiltration (Green-Ampt, curve number), storm-drain inlets.
- Rain as a constant, a hyetograph, a gridded time-varying field, or the surface
  precipitation of ERF's microphysics.
- Terrain, roughness and building footprints from rasters.
- A water-balance log, flood outputs (maximum depth, speed, arrival time),
  checkpoint and restart.

Of these, the first-order scheme and the absolute datum exist; the rest does
not yet.

## Layout

| Path                | Contents                                                          |
| ------------------- | ----------------------------------------------------------------- |
| `Source`            | FLUME's sources.                                                  |
| `Tests`             | Test decks, the test runner and the check scripts.                |
| `Docs/sphinx_doc`   | The documentation (Sphinx): building, model, inputs, testing.     |
| `.github/workflows` | Continuous integration: builds and tests, style, documentation.   |
| `Submodules/ERF`    | ERF, pinned to the tagged release 26.11. AMReX comes through ERF. |

FLUME uses ERF the way [GALE](https://github.com/erf-model/GALE) and FLARE do:
ERF is a submodule, FLUME holds only its own files, and every change FLUME needs
in ERF goes to ERF as a pull request.

## Building and testing

FLUME needs CMake 3.25 or newer, a C++20 compiler and, for the tests, Python 3.
MPI is optional and on by default. Two submodules are needed: ERF, and the AMReX
that ERF pins.

```bash
git clone https://github.com/erf-model/FLUME.git
cd FLUME
git submodule update --init Submodules/ERF
git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
```

The stand-alone build compiles no ERF code. The CMake options, the inputs, the
tests and how to move the ERF submodule are in `Docs/sphinx_doc`:

```bash
python -m pip install -r Docs/sphinx_doc/requirements.txt
sphinx-build -W --keep-going -b html Docs/sphinx_doc Docs/sphinx_doc/_build/html
```

## License

FLUME is released under ERF's license. `LICENSE.md` and `copyright.txt` are the
same files as in ERF, GALE and FLARE.
