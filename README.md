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

This repository is a skeleton. It holds this README and ERF as a git submodule.
There is no solver code, build system, test or documentation yet.

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

None of these exist yet.

## Layout

| Path             | Contents                                                        |
| ---------------- | --------------------------------------------------------------- |
| `Submodules/ERF` | ERF, pinned to the tagged release 26.11. AMReX comes through ERF. |

FLUME uses ERF the way [GALE](https://github.com/erf-model/GALE) and FLARE do:
ERF is a submodule, FLUME holds only its own files, and every change FLUME needs
in ERF goes to ERF as a pull request.

## Getting the code

Clone with the submodules, which brings ERF and, through ERF, AMReX:

```bash
git clone --recursive https://github.com/erf-model/FLUME.git
```

In an existing clone:

```bash
git submodule update --init --recursive
```

To move the submodule to another ERF tag or commit:

```bash
git -C Submodules/ERF fetch --tags
git -C Submodules/ERF checkout <tag-or-commit>
git add Submodules/ERF
git commit -m "Move ERF to <tag-or-commit>"
```

## License

FLUME is released under ERF's license. `LICENSE.md` and `copyright.txt` are the
same files as in ERF, GALE and FLARE.
