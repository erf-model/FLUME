# Rain on the streets of Gowanus, Brooklyn

## Purpose

A first run of FLUME on a real place: a square of 1 km in Gowanus, Brooklyn,
New York, with its terrain and its 1861 building footprints, on cells of 2 m.
The square holds 4th Avenue around Carroll Street, where a FloodNet street
sensor measured 0.91 m (35.94 in) of water on 29 September 2023, the Gowanus
Canal in its north-west corner and the rise to Park Slope in the east.

The case shows how the pieces work together. It is not a hindcast:

- **The rain is simplified.** One rate over the whole square, 50 mm/h for
  one hour (50 mm, 50000 m3), then 30 min with no rain. It is not the rain of
  a storm.
- **No sewers, no infiltration.** All the rain runs over the surface.
- **A closed square.** Every side is a wall: no water comes in from the
  streets uphill and none leaves.

The same rain is run twice: with the buildings, and on the bare terrain.

![The Gowanus case](../../Docs/sphinx_doc/figures/gowanus.png)

The figure: (a) the bed and the buildings, with the three FloodNet sensors
inside the square; (b) the largest depth with the buildings; (c) without
them; (d) the depth at the three sensor positions.

## Files

| File | What it is |
| --- | --- |
| `prepare_gowanus.py` | Downloads the terrain and the footprints (`--download`) and writes the case: `bed.asc`, `buildings.asc`, `rain.txt`, `gowanus.inputs`, `gowanus_open.inputs` (no buildings) and `case.txt`. Needs numpy, matplotlib, Pillow and pyproj. |
| `inputs_gowanus` | The deck: 500 x 500 cells in boxes of 250 for four ranks, first order, Manning n = 0.015, walls, flood maps, a hyetograph. |
| `plot_gowanus.py` | Draws the figure and prints the numbers of the tables below. Needs numpy and matplotlib. |

### How the case is built

- **Square.** UTM zone 18N (EPSG:26918): x from 585400 to 586400 m, y from
  4502900 to 4503900 m. FLUME runs on these coordinates.
- **Bed.** The 1 m bare-earth elevation is averaged over 2 x 2 cells onto
  the flood cells, with 20 m around the square so the bilinear reader has
  the whole domain. The canal is the lidar's water surface (bed below 0 m);
  it is left out of the depth numbers.
- **Buildings.** A flood cell whose centre is inside a footprint is raised
  by the roof height (`flume.buildings.file`); a footprint with no height
  gets 3 m. 37.6 % of the square is building. Roofs shed their rain onto the
  cells beside them.
- **Gaps.** Footprints that almost touch leave slivers of a few cells that
  take roof water and have no way out. Groups of fewer than 25 open cells
  shut in by buildings are closed (546 cells). Without this the deepest cell
  was a sliver 5.1 m deep.
- **Gauges.** One at each FloodNet sensor inside the square: Carroll St /
  4th Av, 4th Av / Garfield Pl, 3rd Ave / Union St.

### Why the sides are walls

The first runs had `outflow` on every side and went wrong: after about
6 min water came in through the east and north sides, and by 20 min the
square held 2e8 m3 from 2e4 m3 of rain. An outflow side copies the boundary
cell to its ghost cell, so where a street falls into the square the film of
rain flows in from the ghost cell, with nothing to limit it. Opening only
the runs of cells where the bed falls outward on both sides of the boundary
worked with the buildings (3900 m3 of the 50000 m3 left) and not without
them (98000 m3 came in through one segment). With walls both runs keep
their rain to round-off. `Docs/sphinx_doc/theory/Boundaries.rst` lists the
limits of the outflow side; this case adds inflow on real terrain to them.

## Data and acknowledgements

Nothing is redistributed in the repository: `prepare_gowanus.py --download`
fetches the two files below (about 5 MB and 2 MB). The figure is drawn from
them.

| Data | Source | Terms |
| --- | --- | --- |
| Terrain, 1 m bare earth, NAVD88 | U.S. Geological Survey, 3D Elevation Program (3DEP), dynamic elevation service of The National Map, `elevation.nationalmap.gov` | U.S. public domain |
| Building footprints and roof heights | City of New York, Office of Technology and Innovation, *Building Footprints*, NYC Open Data, dataset `5zhs-2jue`, `data.cityofnewyork.us` | NYC Open Data terms of use; the dataset page names no licence, so the file is not committed |
| Sensor positions and the depth of 29 September 2023 | FloodNet (New York University, City University of New York and the City of New York), *FloodNet: Sensor Deployment Metadata* (`kb2e-tjy3`) and *Street Flooding Events Measured by FloodNet Sensors* (`aq7i-eu5q`), NYC Open Data | NYC Open Data terms of use; three positions and one depth are quoted here |

Thanks to the U.S. Geological Survey for the 3DEP elevation, to the City of
New York for the open building footprints, and to the FloodNet team for
their open street-flood measurements. The choice of this place follows
Annis, Badas and Mascaro, who model two New York basins against FloodNet
sensors with LISFLOOD-FP; none of their files is used here.

## How to run

From an empty directory, with a stand-alone MPI build:

```bash
python3 <repo>/Tests/Gowanus/prepare_gowanus.py $PWD/data $PWD/case --download
```

```bash
mkdir with && cd with && mpiexec -n 4 <build>/flume <repo>/Tests/Gowanus/inputs_gowanus FILE=../case/gowanus.inputs && cd ..
```

```bash
mkdir without && cd without && mpiexec -n 4 <build>/flume <repo>/Tests/Gowanus/inputs_gowanus FILE=../case/gowanus_open.inputs && cd ..
```

```bash
python3 <repo>/Tests/Gowanus/plot_gowanus.py case with without gowanus.png
```

Each run takes about 50000 steps and 6 min on 4 ranks (Apple clang, macOS
arm64, Release, the two runs side by side).

## Expected Results

From the runs of 28 November 2026 (FLUME at 185aa7a):

| Run | Steps | Rain, m3 | Stored at 5400 s, m3 | Residual, m3 |
| --- | --- | --- | --- | --- |
| With buildings | 50742 | 50000.0 | 50000.0 | -3.4e-10 |
| Without buildings | 49101 | 50000.0 | 50000.0 | -1.6e-10 |

Largest depth of the run, off the buildings and the canal:

| Run | Area above 0.15 m, ha | above 0.3 m | above 0.5 m | above 1 m |
| --- | --- | --- | --- | --- |
| **With buildings** | **8.43** | **4.58** | **2.19** | **0.31** |
| Without buildings | 5.49 | 2.69 | 1.11 | 0.13 |

Peak depth at the sensor positions (time of the peak):

| Run | Carroll St / 4th Av | 4th Av / Garfield Pl | 3rd Ave / Union St |
| --- | --- | --- | --- |
| **With buildings** | **0.83 m (63 min)** | **0.55 m (63 min)** | **0.22 m (61 min)** |
| Without buildings | 0.70 m (68 min) | 0.42 m (68 min) | 0.21 m (88 min) |

What the runs show:

- The water collects on 4th Avenue between Garfield Place and Carroll
  Street, at the two FloodNet sensors there. The terrain alone picks this
  place: the run without buildings pools there too.
- The buildings deepen it. Roofs shed their rain into the streets and the
  walls hold it there: 0.83 m against 0.70 m at Carroll Street, and 1.5 to
  2.4 times the street area above each depth.
- The pool is still 0.76 m deep 30 min after the rain: on the surface it
  has no way out. In the city it leaves through the sewer.
- The canal rises 0.6 m.

The 0.83 m at Carroll Street is close to the 0.91 m the sensor measured on
29 September 2023. This is not a validation: the rain here is not that
storm's, there is no sewer, and the square is closed. The deepest cells
(3.1 m with the buildings) are pits of the terrain, such as ramps and
excavations, not streets.

## CTest

None. The case needs the two downloads. A CTest would need a small synthetic
copy of the rasters or data that may be committed.
