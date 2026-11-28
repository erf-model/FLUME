#!/usr/bin/env python3
"""Prepare the Gowanus street-flood case: terrain, buildings, rain and inputs.

Usage:

    prepare_gowanus.py <data directory> <case directory> [--download]

The case is a square of 1 km in Gowanus, Brooklyn, New York, on cells of 2 m,
in UTM zone 18N (EPSG:26918, metres): x from 585400 to 586400, y from 4502900
to 4503900. It holds 4th Avenue from about 1st Street to Douglass Street, the
Gowanus Canal along its west side and the rise to Park Slope in the east.

<data directory> holds the two downloaded files:

- dem.tif: the bed, USGS 3D Elevation Program (bare earth, m, NAVD88), a
  GeoTIFF of 1040 x 1040 cells of 1 m that holds the square and 20 m around it;
- buildings.geojson: the footprints and roof heights (ft) of New York City's
  Building Footprints, over the square.

With --download the script fetches them first (about 5 MB and 2 MB); without
it they must be there. It writes into <case directory>:

- bed.asc: the bed on cells of 2 m, the mean of the 2 x 2 cells of 1 m;
- buildings.asc: the height above the bed of each cell of 2 m whose centre is
  inside a footprint, 0 elsewhere (flume.buildings.file); groups of fewer
  than GAP_CELLS open cells shut in by buildings are closed;
- rain.txt: the hyetograph, RAIN_RATE mm/h for RAIN_TIME s, then none;
- gowanus.inputs and gowanus_open.inputs: the file names and a gauge at each
  of the three FloodNet sensors inside the square (a sensor in a building
  cell reads the nearest open cell), with and without the buildings, read
  with FILE=;
- case.txt: the counts and the rain, for the README and the figure.

The rain is not that of a storm: it is one rate over the whole square, to
show the flood the terrain and the buildings give. There are no sewers.

Needs numpy, matplotlib (its Path, for the footprints), Pillow (the GeoTIFF)
and pyproj (longitude and latitude to UTM).
"""
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from matplotlib.path import Path as Outline
from PIL import Image
from pyproj import Transformer

X0, Y0 = 585400.0, 4502900.0   # lower corner of the square, UTM 18N, m
SIZE = 1000.0                  # m
DX = 2.0                       # flood cells, m
MARGIN = 20.0                  # of the bed raster around the square, m
FOOT = 0.3048                  # m
MIN_HEIGHT = 3.0               # the least height of a footprint (also one with none), m
GAP_CELLS = 25                 # open groups below this, shut in by buildings, are closed
WATER = 0.0                    # bed below this is the canal (the lidar's water surface), m
RAIN_RATE = 50.0               # mm/h
RAIN_TIME = 3600.0             # s
STOP_TIME = 5400.0             # s
# The FloodNet sensors inside the square: name, longitude, latitude (FloodNet:
# Sensor Deployment Metadata, NYC Open Data, kb2e-tjy3, read on 2026-11-28).
GAUGES = (("carroll_4th", -73.984089, 40.676281),
          ("garfield_4th", -73.983944, 40.675855),
          ("union_3rd", -73.985293, 40.678369))

DEM_URL = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"
BUILDINGS_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
LIMIT = 50000                  # footprints asked for; the square holds about 2000


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def download(data):
    data.mkdir(parents=True, exist_ok=True)
    lo = (X0 - MARGIN, Y0 - MARGIN)
    hi = (X0 + SIZE + MARGIN, Y0 + SIZE + MARGIN)
    n = int(SIZE + 2 * MARGIN)
    dem = {"bbox": f"{lo[0]:.0f},{lo[1]:.0f},{hi[0]:.0f},{hi[1]:.0f}",
           "bboxSR": 26918, "imageSR": 26918, "size": f"{n},{n}", "format": "tiff",
           "pixelType": "F32", "interpolation": "RSP_BilinearInterpolation", "f": "image"}
    to_ll = Transformer.from_crs(26918, 4326, always_xy=True)
    corners = [to_ll.transform(x, y) for x in (lo[0], hi[0]) for y in (lo[1], hi[1])]
    west, east = min(c[0] for c in corners), max(c[0] for c in corners)
    south, north = min(c[1] for c in corners), max(c[1] for c in corners)
    bld = {"$where": f"within_box(the_geom, {north:.5f}, {west:.5f}, {south:.5f}, {east:.5f})",
           "$limit": LIMIT}
    for url, query, name in ((DEM_URL, dem, "dem.tif"), (BUILDINGS_URL, bld, "buildings.geojson")):
        full = url + "?" + urllib.parse.urlencode(query)
        print(f"downloading {name} from {url}")
        with urllib.request.urlopen(full, timeout=120) as r:
            (data / name).write_bytes(r.read())


def read_bed(path):
    """The bed on the cells of 2 m of the square and MARGIN around it, north first."""
    try:
        im = Image.open(path)
        tags = im.tag_v2
    except (OSError, AttributeError) as e:
        fail(f"{path}: not a TIFF image ({e}); download it again")
    scale, tie, keys = tags.get(33550), tags.get(33922), tags.get(34735)
    if scale is None or tie is None or keys is None:
        fail(f"{path}: not a GeoTIFF (no pixel scale, tie point or keys)")
    # GeoKey 1025 (raster type) must be 1: a pixel is an area, not a point.
    kinds = [keys[k + 3] for k in range(4, len(keys) - 3, 4) if keys[k] == 1025]
    if kinds != [1] or tuple(tie[:3]) != (0.0, 0.0, 0.0):
        fail(f"{path}: raster type {kinds} tied at pixel {tuple(tie[:3])}, "
             f"not pixel-is-area tied at its corner")
    n = int(SIZE + 2 * MARGIN)
    want = (X0 - MARGIN, Y0 + SIZE + MARGIN)
    if im.size != (n, n) or scale[:2] != (1.0, 1.0) or (tie[3], tie[4]) != want:
        fail(f"{path}: {im.size} cells of {scale[:2]} m from {tie[3:5]}, "
             f"not {n} x {n} of 1 m from {want}")
    a = np.array(im, dtype=np.float64)
    if not np.isfinite(a).all() or a.min() < -50.0 or a.max() > 500.0:
        fail(f"{path}: elevations from {a.min():g} to {a.max():g} m; missing data?")
    # Rounded as bed.asc is written, so the counts here are those of the run.
    return np.round(a.reshape(n // 2, 2, n // 2, 2).mean(axis=(1, 3)), 3)


def read_buildings(path):
    """Heights above the bed on the flood cells, [j, i] from the south-west, and the count."""
    try:
        feats = json.loads(Path(path).read_text())["features"]
    except (ValueError, KeyError, TypeError) as e:
        fail(f"{path}: not a GeoJSON feature collection ({e!r}); download it again")
    if len(feats) >= LIMIT:
        fail(f"{path}: {len(feats)} footprints, the limit of the query; some are missing")
    # WGS 84 to NAD83 with no shift: about 1 m, below the cell.
    to_utm = Transformer.from_crs(4326, 26918, always_xy=True)
    n = int(SIZE / DX)
    height = np.zeros((n, n))
    used = 0
    for f in feats:
        g = f["geometry"]
        if g is None or g["type"] not in ("Polygon", "MultiPolygon"):
            continue
        raw = f["properties"].get("height_roof")
        h = max(float(raw) * FOOT if raw not in (None, "") else 0.0, MIN_HEIGHT)
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        hit = False
        for rings in polys:
            xy = [np.column_stack(to_utm.transform(*np.array(r).T[:2])) for r in rings]
            lo, hi = xy[0].min(axis=0), xy[0].max(axis=0)
            i0 = max(int(np.floor((lo[0] - X0) / DX)), 0)
            i1 = min(int(np.ceil((hi[0] - X0) / DX)), n)
            j0 = max(int(np.floor((lo[1] - Y0) / DX)), 0)
            j1 = min(int(np.ceil((hi[1] - Y0) / DX)), n)
            if i0 >= i1 or j0 >= j1:
                continue
            xc = X0 + (np.arange(i0, i1) + 0.5) * DX
            yc = Y0 + (np.arange(j0, j1) + 0.5) * DX
            pts = np.column_stack([a.ravel() for a in np.meshgrid(xc, yc)])
            # inside the outer ring and in none of the holes
            inside = Outline(xy[0]).contains_points(pts)
            for hole in xy[1:]:
                inside &= ~Outline(hole).contains_points(pts)
            inside = inside.reshape(len(yc), len(xc))
            if inside.any():
                hit = True
                part = height[j0:j1, i0:i1]
                part[inside] = np.maximum(part[inside], h)
        used += hit
    if not used:
        fail(f"{path}: no footprint of its {len(feats)} covers a cell of the square")
    return height, used


def fill_gaps(height):
    """Close the open cells that buildings shut in, in groups below GAP_CELLS.

    Footprints that almost touch leave slivers of one or two cells between
    them. A sliver takes the rain of the roofs around it and has no way out,
    so it fills to metres. Each group of open cells joined by their sides and
    smaller than GAP_CELLS gets the height of the lowest building beside it.
    Returns the number of cells closed.
    """
    nj, ni = height.shape
    seen = height > 0.0
    closed = 0
    for j0, i0 in zip(*np.nonzero(~seen)):
        if seen[j0, i0]:
            continue
        seen[j0, i0] = True
        group, stack, beside = [], [(j0, i0)], []
        while stack:
            j, i = stack.pop()
            group.append((j, i))
            for jj, ii in ((j + 1, i), (j - 1, i), (j, i + 1), (j, i - 1)):
                if 0 <= jj < nj and 0 <= ii < ni:
                    if height[jj, ii] > 0.0:
                        beside.append(height[jj, ii])
                    elif not seen[jj, ii]:
                        seen[jj, ii] = True
                        stack.append((jj, ii))
        if len(group) < GAP_CELLS and beside:
            for j, i in group:
                height[j, i] = -min(beside)  # below 0 until every group is found
            closed += len(group)
    np.abs(height, out=height)
    return closed


def write_asc(path, a_north_first, x0, y0, cell, fmt):
    rows, cols = a_north_first.shape
    with open(path, "w") as f:
        f.write(f"ncols {cols}\nnrows {rows}\nxllcorner {x0:.1f}\nyllcorner {y0:.1f}\n"
                f"cellsize {cell:g}\n")
        np.savetxt(f, a_north_first, fmt=fmt)


def main():
    args = [a for a in sys.argv[1:] if a != "--download"]
    if len(args) != 2 or any(a.startswith("--") for a in args):
        sys.exit("usage: prepare_gowanus.py <data directory> <case directory> [--download]")
    data, case = Path(args[0]).resolve(), Path(args[1]).resolve()
    if "--download" in sys.argv[1:]:
        download(data)
    for name in ("dem.tif", "buildings.geojson"):
        if not (data / name).is_file():
            fail(f"{data / name} is missing; run with --download")
    case.mkdir(parents=True, exist_ok=True)

    bed = read_bed(data / "dem.tif")
    write_asc(case / "bed.asc", bed, X0 - MARGIN, Y0 - MARGIN, DX, "%.3f")

    height, used = read_buildings(data / "buildings.geojson")
    closed = fill_gaps(height)
    write_asc(case / "buildings.asc", height[::-1], X0, Y0, DX, "%.2f")

    (case / "rain.txt").write_text(
        f"# time (s)  rate (mm/h): one rate over the square, not a storm's\n"
        f"0 {RAIN_RATE:g}\n{RAIN_TIME:g} 0\n")

    m = int(MARGIN / DX)
    south = bed[::-1][m:-m, m:-m]  # the bed of the flood cells, [j, i] from the south-west
    to_utm = Transformer.from_crs(4326, 26918, always_xy=True)
    jj, ii = np.mgrid[0:height.shape[0], 0:height.shape[1]]
    gauges = []
    for name, lon, lat in GAUGES:
        gx, gy = to_utm.transform(lon, lat)
        gi, gj = int(np.floor((gx - X0) / DX)), int(np.floor((gy - Y0) / DX))
        if not (0 <= gi < height.shape[1] and 0 <= gj < height.shape[0]):
            fail(f"the gauge {name} at ({gx:.1f}, {gy:.1f}) is outside the square")
        moved = 0.0
        if height[gj, gi] > 0.0:
            # a sensor on a pole at a wall: read the nearest open cell
            dist = np.where(height == 0.0, np.hypot(ii - gi, jj - gj), np.inf)
            gj, gi = np.unravel_index(np.argmin(dist), dist.shape)
            moved = dist[gj, gi] * DX
            gx, gy = X0 + (gi + 0.5) * DX, Y0 + (gj + 0.5) * DX
        gauges.append((name, gx, gy, south[gj, gi], moved))

    common = (f"flume.prob_lo = {X0:.1f} {Y0:.1f}\n"
              f"flume.prob_hi = {X0 + SIZE:.1f} {Y0 + SIZE:.1f}\n"
              f"flume.stop_time = {STOP_TIME:g}\n"
              f"flume.bed.file = \"{case / 'bed.asc'}\"\n"
              f"flume.rain.file = \"{case / 'rain.txt'}\"\n"
              f"flume.gauges.names = {' '.join(g[0] for g in gauges)}\n"
              + "".join(f"flume.gauges.{g[0]}.location = {g[1]:.2f} {g[2]:.2f}\n"
                        for g in gauges))
    (case / "gowanus.inputs").write_text(
        common + f"flume.buildings.file = \"{case / 'buildings.asc'}\"\n")
    (case / "gowanus_open.inputs").write_text(common)

    cells = int((height > 0).sum())
    rain_mm = RAIN_RATE * RAIN_TIME / 3600.0
    (case / "case.txt").write_text(
        f"x0 {X0:.1f}\ny0 {Y0:.1f}\nsize {SIZE:g}\ndx {DX:g}\n"
        f"bed_min {south.min():.3f}\nbed_max {south.max():.3f}\n"
        f"footprints {used}\nbuilding_cells {cells}\ngap_cells_closed {closed}\n"
        f"water_below {WATER:g}\nwater_cells {int((south < WATER).sum())}\n"
        f"building_fraction {cells / height.size:.4f}\n"
        f"rain_rate {RAIN_RATE:g}\nrain_time {RAIN_TIME:g}\nstop_time {STOP_TIME:g}\n"
        f"rain_mm {rain_mm:g}\nrain_m3 {rain_mm * 1e-3 * SIZE * SIZE:g}\n"
        + "".join(f"gauge {g[0]} {g[1]:.2f} {g[2]:.2f} bed {g[3]:.3f} moved {g[4]:g}\n"
                  for g in gauges))
    print((case / "case.txt").read_text(), end="")


if __name__ == "__main__":
    main()
