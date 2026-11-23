#!/usr/bin/env python3
"""Write the rain grids and index files of the gridded-rain tests.

    make_rain_grids.py [--dir DIR]

The files are committed; this script is how they were made, and running it
again writes them byte for byte. Every grid is synthetic: no data is
downloaded. See README.md for what each file tests.
"""
import argparse
import math
from pathlib import Path


def write(path, nx, ny, origin, cell, values, nodata=None, centre=False):
    """values[j][i]: column i from the west, row j from the south, mm/h. The
    file lists the northern row first. cell is a number (cellsize) or (dx, dy);
    origin is the lower-left corner, or the lower-left centre with centre=True."""
    lines = [f"ncols {nx}", f"nrows {ny}"]
    tag = "center" if centre else "corner"
    lines.append(f"xll{tag} {origin[0]!r}")
    lines.append(f"yll{tag} {origin[1]!r}")
    if isinstance(cell, tuple):
        lines += [f"dx {cell[0]!r}", f"dy {cell[1]!r}"]
    else:
        lines.append(f"cellsize {cell!r}")
    if nodata is not None:
        lines.append(f"NODATA_value {nodata!r}")
    for j in reversed(range(ny)):
        lines.append(" ".join(repr(values[j][i]) for i in range(nx)))
    path.write_text("\n".join(lines) + "\n")


def noise(i, j, seed):
    """A repeatable number in [0, 1) for cell (i, j)."""
    h = (i * 73856093) ^ (j * 19349663) ^ (seed * 83492791)
    return (h % 10007) / 10007.0


def field(nx, ny, x0, y0, dx, dy, f, seed):
    """f at the centres of the cells, plus up to 10 mm/h of noise, to 3 decimals."""
    return [[round(f(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy) + 10.0 * noise(i, j, seed), 3)
             for i in range(nx)] for j in range(ny)]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=str(Path(__file__).resolve().parent))
    out = Path(p.parse_args().dir)
    bad = out / "bad"
    bad.mkdir(parents=True, exist_ok=True)

    # grid_a, grid_b: 18 x 15 cells of 3.7 x 2.9 m from (-5.3, -4.1), so no
    # raster edge falls on a solver edge of the 2 m test grids, and the grid
    # holds the domain (0, 0)-(48, 32) and the shifted one (10, 5)-(58, 37).
    # One NODATA cell in the south-west corner, outside every domain.
    ax, ay, adx, ady, anx, any_ = -5.3, -4.1, 3.7, 2.9, 18, 15
    a = field(anx, any_, ax, ay, adx, ady,
              lambda x, y: 20.0 + 30.0 * (1.0 + math.sin(0.13 * x + 0.4) * math.cos(0.11 * y)), 1)
    a[0][0] = -9999.0
    write(out / "grid_a.asc", anx, any_, (ax, ay), (adx, ady), a, nodata=-9999.0)
    b = field(anx, any_, ax, ay, adx, ady,
              lambda x, y: 80.0 * math.exp(-((x - 30.0) ** 2 + (y - 20.0) ** 2) / 200.0), 2)
    write(out / "grid_b.asc", anx, any_, (ax, ay), (adx, ady), b)
    # grid_c: a finer grid, 1.3 m cells given by their lower-left centre.
    cx, cy, cd, cnx, cny = -1.3, 0.25, 1.3, 50, 32
    c = field(cnx, cny, cx - 0.5 * cd, cy - 0.5 * cd, cd, cd,
              lambda x, y: 5.0 + 0.8 * x + 0.5 * y, 3)
    write(out / "grid_c.asc", cnx, cny, (cx, cy), cd, c, centre=True)

    (out / "index_steps.txt").write_text(
        "# Rain grids for the gridded-rain tests: time (s) and grid (mm/h).\n"
        "# Paths are relative to this file.\n"
        "0     grid_a.asc\n"
        "300   grid_b.asc   # a peak near (30, 20)\n"
        "900   grid_c.asc\n")
    (out / "index_offset.txt").write_text(
        "# Two grids around t = 0, so the first step of a run interpolates.\n"
        "-300  grid_a.asc\n"
        "600   grid_b.asc\n")
    (out / "index_fine.txt").write_text("0 grid_c.asc\n")
    # For the dry, closed domain of Tests/Skeleton/inputs, (-10, 5)-(86, 69):
    # no rain, then a burst of about 1 mm/s from 100 s to 101 s, then none.
    zero = [[0.0] * 11 for _ in range(8)]
    write(out / "grid_zero.asc", 11, 8, (-20.0, 0.0), 10.0, zero)
    burst = field(11, 8, -20.0, 0.0, 10.0, 10.0, lambda x, y: 3590.0, 4)
    write(out / "grid_burst.asc", 11, 8, (-20.0, 0.0), 10.0, burst)
    (out / "index_burst.txt").write_text(
        "0    grid_zero.asc\n"
        "100  grid_burst.asc\n"
        "101  grid_zero.asc\n")

    # grid_aligned: 1 m cells from (-1, -1) to (16, 16), NODATA in the ring
    # outside 0..14 m. On 25 cells of 14/25 m the solver's last edge, formed
    # in raster cells, rounds to 15.000000000000002: a cell that shares no area
    # with the domain must not count as read.
    al = field(17, 17, -1.0, -1.0, 1.0, 1.0, lambda x, y: 30.0 + 2.0 * x - y, 5)
    for j in range(17):
        for i in range(17):
            if i in (0, 15, 16) or j in (0, 15, 16):
                al[j][i] = -9999.0
    write(out / "grid_aligned.asc", 17, 17, (-1.0, -1.0), 1.0, al, nodata=-9999.0)
    (out / "index_aligned.txt").write_text("0 grid_aligned.asc\n")

    # A grid that stops at x = 40, short of the domain's east side (48 m).
    write(bad / "grid_short.asc", 10, 10, (0.0, 0.0), 4.0, [[10.0] * 10 for _ in range(10)])
    # A rate below 0 in a cell the domain reads.
    neg = [[10.0] * 14 for _ in range(10)]
    neg[4][6] = -2.5
    write(bad / "grid_negative.asc", 14, 10, (-4.0, -4.0), 4.0, neg)
    # NODATA in a cell the domain reads.
    nod = [[10.0] * 14 for _ in range(10)]
    nod[5][3] = -1.0
    write(bad / "grid_nodata.asc", 14, 10, (-4.0, -4.0), 4.0, nod, nodata=-1.0)
    for name, text in [
        ("one_item", "0 grid_short.asc\n300\n"),
        ("three_items", "0 grid_short.asc\n300 grid_short.asc extra\n"),
        ("time_not_a_number", "0 grid_short.asc\nlater grid_short.asc\n"),
        ("time_infinite", "0 grid_short.asc\ninf grid_short.asc\n"),
        ("not_increasing", "0 ../grid_a.asc\n600 ../grid_b.asc\n300 ../grid_a.asc\n"),
        ("empty", "# no entry\n\n"),
        ("missing_grid", "0 ../grid_a.asc\n300 no_such_grid.asc\n"),
        ("short", "0 ../grid_a.asc\n300 grid_short.asc\n"),
        ("negative", "0 grid_negative.asc\n"),
        ("nodata", "0 grid_nodata.asc\n"),
    ]:
        (bad / f"index_{name}.txt").write_text(text)


if __name__ == "__main__":
    main()
