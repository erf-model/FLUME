#!/usr/bin/env python3
"""Prepare the town case from the rain of the storm case.

Usage:

    prepare_town.py <index> <directory>

<index> is the index of gridded rain that the coupled storm run writes
(Tests/ERFStorm, rain/storm.txt): lines of a time in s and an ESRI ASCII grid
of the rain rate in mm/h, 250 m cells. The script picks, inside the storm
case's town (x from -10 to -2 km, y from -12 to -4 km), the square of
SIZE = 1 km on the 250 m lines and the window of WINDOW = 1200 s, starting at
an entry of the index, with the most rain. It writes into <directory>:

- rain.txt: the index of that window, its times shifted to start at 0, the
  paths those of the storm's grids (absolute), and a grid of no rain at
  WINDOW;
- town.inputs: geometry.prob_lo and geometry.prob_hi of the square (ERF's
  domain, 240 m high), erf.buildings_file_name and flume.buildings.file,
  read with FILE=;
- buildings_erf.txt: ERF's building heights (erf.buildings_type =
  ImmersedForcing) on ERF's nodes of 10 m;
- buildings.asc: FLUME's building heights (flume.buildings.file) on the
  flood cells of 2.5 m;
- town_flume.inputs: flume.prob_lo, flume.prob_hi, flume.buildings.file
  and flume.rain.file, for the stand-alone run with the buildings;
- town_flume_open.inputs: the same without flume.buildings.file, for the
  stand-alone run without them;
- window.txt: the square, the window and its rain, for the check.

The town: 10 x 10 blocks of 100 m, each a building of 60 x 60 m in its middle
(x and y from 20 to 80 m in the block) with streets of 40 m between them. The
heights are 10 to 40 m in steps of 5 m, set by the block's place. ERF's nodes
inside or on a footprint get its height. FLUME's cells whose centre is inside
a footprint get its height: 24 x 24 cells of 2.5 m, the 6 x 6 ERF cells that
the footprint fills.

The choice is made from the rain itself, so the case holds where the storm's
cell moves on another compiler or processor.
"""
import bisect
import math
import sys
from pathlib import Path

TOWN = (-10000.0, -2000.0, -12000.0, -4000.0)  # the storm case's town, m
RAIN_CELL = 250.0
SIZE = 1000.0
WINDOW = 1200.0
HEIGHT = 240.0
ERF_DX = 10.0
FLOOD_DX = 2.5
BLOCK = 100.0
FOOT = (20.0, 80.0)  # footprint in its block, m


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def read_index(path):
    rows = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        s = line.split("#")[0].split()
        if not s:
            continue
        if len(s) != 2:
            fail(f"{path}:{n}: not a time and a path")
        rows.append((float(s[0]), s[1]))
    if len(rows) < 2:
        fail(f"{path}: fewer than 2 entries")
    for a, b in zip(rows, rows[1:]):
        if not b[0] > a[0]:
            fail(f"{path}: the times do not increase at {b[0]}")
    return rows


def read_asc(path):
    """An ESRI ASCII grid: the header and the rows, north first."""
    lines = Path(path).read_text().split("\n")
    head = {}
    n = 0
    while n < len(lines) and lines[n].split() and lines[n].split()[0][0].isalpha():
        k, v = lines[n].split()
        head[k.lower()] = float(v)
        n += 1
    if "nodata_value" in head:
        fail(f"{path}: a NODATA_value; the storm's grids have none")
    rows = [[float(v) for v in line.split()] for line in lines[n:] if line.strip()]
    if len(rows) != int(head["nrows"]) or any(len(r) != int(head["ncols"]) for r in rows):
        fail(f"{path}: the rows do not match the header")
    return head, rows


class Grid:
    """A rain grid: the rate of the 250 m cell at (x, y), mm/h."""

    def __init__(self, path):
        h, rows = read_asc(path)
        if h["cellsize"] != RAIN_CELL:
            fail(f"{path}: cells of {h['cellsize']:g} m, not {RAIN_CELL:g} m")
        self.x0, self.y0 = h["xllcorner"], h["yllcorner"]
        self.rows = rows
        self.nrows = len(rows)

    def rate(self, x, y):
        i = math.floor((x - self.x0) / RAIN_CELL)
        j = math.floor((y - self.y0) / RAIN_CELL)
        if not (0 <= i < len(self.rows[0]) and 0 <= j < self.nrows):
            fail(f"the rain grid does not hold ({x:g}, {y:g})")
        return self.rows[self.nrows - 1 - j][i]


def squares():
    """The lower corners of the squares of SIZE inside TOWN on the 250 m lines."""
    n = int(SIZE / RAIN_CELL)
    nx = int((TOWN[1] - TOWN[0]) / RAIN_CELL) - n + 1
    ny = int((TOWN[3] - TOWN[2]) / RAIN_CELL) - n + 1
    return [(TOWN[0] + i * RAIN_CELL, TOWN[2] + j * RAIN_CELL)
            for j in range(ny) for i in range(nx)]


def square_mean(grid, x0, y0):
    """The mean rate over the square at (x0, y0), mm/h: the mean of its cells."""
    n = int(SIZE / RAIN_CELL)
    return math.fsum(grid.rate(x0 + (i + 0.5) * RAIN_CELL, y0 + (j + 0.5) * RAIN_CELL)
                     for j in range(n) for i in range(n)) / (n * n)


def window_depth(times, means, k0):
    """The depth over [times[k0], times[k0] + WINDOW] of the step rates, mm."""
    t1 = times[k0] + WINDOW
    parts = []
    for k in range(k0, len(times)):
        if times[k] >= t1:
            break
        end = min(times[k + 1], t1) if k + 1 < len(times) else t1
        parts.append(means[k] * (end - times[k]) / 3600.0)
    return math.fsum(parts)


def choose(index):
    """The square and window of most rain: (depth in mm, x0, y0, k0). The
    candidates are ranked with sums over the town's cells; the depth returned
    is that of window_depth, summed with math.fsum."""
    rows = read_index(index)
    base = Path(index).resolve().parent
    times = [t for t, _ in rows]
    grids = [Grid(base / p) for _, p in rows]
    n = int(SIZE / RAIN_CELL)
    nx = int((TOWN[1] - TOWN[0]) / RAIN_CELL)
    ny = int((TOWN[3] - TOWN[2]) / RAIN_CELL)
    corners = squares()
    # depth[c][k]: the mean depth over square c of entry k, mm
    depth = [[0.0] * (len(times) - 1) for _ in corners]
    for k in range(len(times) - 1):
        g = grids[k]
        cells = [[g.rate(TOWN[0] + (i + 0.5) * RAIN_CELL, TOWN[2] + (j + 0.5) * RAIN_CELL)
                  for i in range(nx)] for j in range(ny)]
        if not any(any(r) for r in cells):
            continue
        dt = (times[k + 1] - times[k]) / 3600.0
        for c, (x0, y0) in enumerate(corners):
            i0 = int(round((x0 - TOWN[0]) / RAIN_CELL))
            j0 = int(round((y0 - TOWN[2]) / RAIN_CELL))
            s = sum(sum(cells[j][i0:i0 + n]) for j in range(j0, j0 + n))
            depth[c][k] = s / (n * n) * dt
    best = None
    for c, (x0, y0) in enumerate(corners):
        d = depth[c]
        for k0 in range(len(times) - 1):
            if times[k0] + WINDOW > times[-1]:
                break
            k1 = bisect.bisect_left(times, times[k0] + WINDOW)
            v = sum(d[k0:k1])
            t1 = times[k0] + WINDOW
            if times[k1] > t1:  # the last entry holds past the window
                v -= d[k1 - 1] * (times[k1] - t1) / (times[k1] - times[k1 - 1])
            if best is None or v > best[0]:
                best = (v, x0, y0, k0)
    if best is None or best[0] <= 0.0:
        fail(f"{index}: no rain falls on the town in any window of {WINDOW:g} s")
    _, x0, y0, k0 = best
    means = [square_mean(g, x0, y0) for g in grids]
    return rows, base, (window_depth(times, means, k0), x0, y0, k0)


def height(i, j):
    """The height of the building of block (i, j), m: 10 to 40 in steps of 5."""
    return 10.0 + 5.0 * ((3 * i + 5 * j) % 7)


def buildings():
    """The footprints, (x0, x1, y0, y1, height) relative to the square's corner."""
    n = int(SIZE / BLOCK)
    return [(i * BLOCK + FOOT[0], i * BLOCK + FOOT[1], j * BLOCK + FOOT[0], j * BLOCK + FOOT[1],
             height(i, j)) for j in range(n) for i in range(n)]


def write_erf(path, x0, y0):
    """ERF's terrain-text table on its nodes: nx, ny, the x nodes, the y nodes,
    then z[ix * ny + iy], one per line (ERF_ProbCommon.H, read_custom_terrain)."""
    n = int(SIZE / ERF_DX) + 1
    feet = buildings()

    def z(px, py):
        for a, b, c, d, h in feet:
            if a <= px <= b and c <= py <= d:
                return h
        return 0.0

    out = [str(n), str(n)]
    out += [f"{x0 + i * ERF_DX:.1f}" for i in range(n)]
    out += [f"{y0 + j * ERF_DX:.1f}" for j in range(n)]
    out += [f"{z(i * ERF_DX, j * ERF_DX):.1f}" for i in range(n) for j in range(n)]
    Path(path).write_text("\n".join(out) + "\n")


def write_flume(path, x0, y0):
    """FLUME's ESRI grid of building heights on the flood cells."""
    n = int(SIZE / FLOOD_DX)
    feet = buildings()

    def z(px, py):
        for a, b, c, d, h in feet:
            if a < px < b and c < py < d:
                return h
        return 0.0

    out = [f"ncols {n}", f"nrows {n}", f"xllcorner {x0:.1f}", f"yllcorner {y0:.1f}",
           f"cellsize {FLOOD_DX}", "NODATA_value -9999"]
    for j in reversed(range(n)):
        py = (j + 0.5) * FLOOD_DX
        out.append(" ".join(f"{z((i + 0.5) * FLOOD_DX, py):g}" for i in range(n)))
    Path(path).write_text("\n".join(out) + "\n")


def main():
    if len(sys.argv) != 3:
        fail("usage: prepare_town.py <index> <directory>")
    index = Path(sys.argv[1])
    out = Path(sys.argv[2])
    if not index.is_file():
        fail(f"no index {index}")
    rows, base, (depth, x0, y0, k0) = choose(index)
    out.mkdir(parents=True, exist_ok=True)
    t0 = rows[k0][0]
    zero = None
    for t, p in rows:
        g = Grid(base / p)
        if all(v == 0.0 for r in g.rows for v in r):
            zero = p
            break
    if zero is None:
        fail(f"{index}: no grid of no rain to close the window with")
    lines = [f"# The rain of {index.resolve()} over x {x0:g} to {x0 + SIZE:g} m, "
             f"y {y0:g} to {y0 + SIZE:g} m, from {t0:g} s, shifted to start at 0"]
    for t, p in rows[k0:]:
        if t - t0 >= WINDOW:
            break
        lines.append(f"{t - t0:.17g} {base / p}")
    lines.append(f"{WINDOW:.17g} {base / zero}")
    (out / "rain.txt").write_text("\n".join(lines) + "\n")
    write_erf(out / "buildings_erf.txt", x0, y0)
    write_flume(out / "buildings.asc", x0, y0)
    (out / "town.inputs").write_text(
        f"# The square of the town case, written by prepare_town.py\n"
        f"geometry.prob_lo = {x0:.1f} {y0:.1f} 0.0\n"
        f"geometry.prob_hi = {x0 + SIZE:.1f} {y0 + SIZE:.1f} {HEIGHT:.1f}\n"
        f"erf.buildings_file_name = {(out / 'buildings_erf.txt').resolve()}\n"
        f"flume.buildings.file = {(out / 'buildings.asc').resolve()}\n"
        f"flume.rain.file = {(out / 'rain.txt').resolve()}\n")
    (out / "town_flume.inputs").write_text(
        f"# The square of the town case for the stand-alone runs, written by prepare_town.py\n"
        f"flume.prob_lo = {x0:.1f} {y0:.1f}\n"
        f"flume.prob_hi = {x0 + SIZE:.1f} {y0 + SIZE:.1f}\n"
        f"flume.buildings.file = {(out / 'buildings.asc').resolve()}\n"
        f"flume.rain.file = {(out / 'rain.txt').resolve()}\n")
    (out / "town_flume_open.inputs").write_text(
        f"# The square of the town case without its buildings, written by prepare_town.py\n"
        f"flume.prob_lo = {x0:.1f} {y0:.1f}\n"
        f"flume.prob_hi = {x0 + SIZE:.1f} {y0 + SIZE:.1f}\n"
        f"flume.rain.file = {(out / 'rain.txt').resolve()}\n")
    (out / "window.txt").write_text(f"{x0:.1f} {y0:.1f} {SIZE:.1f} {t0:.17g} {WINDOW:.17g} "
                                    f"{depth:.17g}\n")
    print(f"prepare_town.py: x {x0:g} to {x0 + SIZE:g} m, y {y0:g} to {y0 + SIZE:g} m, "
          f"from {t0:g} s for {WINDOW:g} s: {depth:.4g} mm of rain "
          f"({len(lines) - 1} entries); {len(buildings())} buildings")


if __name__ == "__main__":
    main()
