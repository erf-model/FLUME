#!/usr/bin/env python3
"""Write the ESRI ASCII grids of the raster tests (Tests/Raster/*.asc).

    make_rasters.py [--dir DIR]

The files are committed; this script is how they were made, and running it
again writes them byte for byte. Every grid is synthetic: no data is
downloaded. See README.md for what each file tests.
"""
import argparse
from pathlib import Path


def plane(a, b, c, d):
    """z = a + b x + c y + d x y, which bilinear interpolation reproduces exactly."""
    return lambda x, y: a + b * x + c * y + d * x * y


def write(path, nx, ny, origin, cell, values, nodata=None, centre=False):
    """values[j][i]: column i from the west, row j from the south. The file
    lists the northern row first. cell is a number (cellsize) or (dx, dy)."""
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=str(Path(__file__).resolve().parent))
    out = Path(p.parse_args().dir)
    bad = out / "bad"
    bad.mkdir(parents=True, exist_ok=True)

    # 1. Bed on the 48 x 32 m domain of inputs_bed: z = 2 + 0.03 x - 0.02 y
    #    + 0.001 x y at the centres of 5 m cells, given by xllcenter. The
    #    centres span -13 to 52 in x and -11 to 39 in y, so the two westmost
    #    columns lie outside what the domain reads (crop_to_domain drops
    #    them) and every cell centre of the domain is inside the hull of the
    #    raster centres. Column 0 holds a NODATA cell, which no cell reads.
    z = plane(2.0, 0.03, -0.02, 0.001)
    nx, ny, x0, y0, h = 14, 11, -13.0, -11.0, 5.0
    vals = [[round(z(x0 + i * h, y0 + j * h), 6) for i in range(nx)] for j in range(ny)]
    vals[4][0] = -9999.0
    write(out / "bed_plane.asc", nx, ny, (x0, y0), h, vals, nodata=-9999.0, centre=True)

    # 2. Bed on the same domain from a coarse raster that holds the domain
    #    exactly: 6 x 8 cells of 8 by 4 m (dx and dy), corners at 0 0. The
    #    cells of the domain within half a raster cell of its sides take the
    #    edge values. Values: a pattern with no plane in it.
    nx, ny = 6, 8
    vals = [[round(1.0 + 0.5 * ((3 * i + 5 * j) % 7) - 0.25 * (i * j % 3), 6) for i in range(nx)]
            for j in range(ny)]
    write(out / "bed_coarse.asc", nx, ny, (0.0, 0.0), (8.0, 4.0), vals)

    # 3. Streets: a bed falling to the east (2 x 2 cells of 20 m holding
    #    [0, 40]^2, 0.075 m at the western centres and 0.025 m at the eastern
    #    ones), and blocks about 4 x 4 m with 4 m streets east of x = 12, on a
    #    2 m raster whose corner is at -0.75 -0.75: a grid centre lies a
    #    quarter or three quarters of a raster cell from a raster edge, so a
    #    half-cell shift in the sampling moves the footprints. Heights 3, 5 or
    #    8 m by block; 0 in the streets. The 0/1 mask has the same footprints.
    write(out / "streets_bed.asc", 2, 2, (0.0, 0.0), 20.0, [[0.075, 0.025], [0.075, 0.025]])
    n, c0 = 21, -0.75
    heights = [[0.0] * n for _ in range(n)]
    mask = [[0] * n for _ in range(n)]
    for j in range(n):
        for i in range(n):
            x, y = c0 + 2 * i + 1, c0 + 2 * j + 1
            if x > 12 and (x - 12) % 8 < 4 and (y - 4) % 8 < 4 and 4 < y < 36:
                bx, by = int((x - 12) // 8), int((y - 4) // 8)
                heights[j][i] = (3.0, 5.0, 8.0)[(bx + 2 * by) % 3]
                mask[j][i] = 1
    write(out / "streets_heights.asc", n, n, (c0, c0), 2.0, heights, nodata=-1.0)
    write(out / "streets_mask.asc", n, n, (c0, c0), 2.0, mask)

    # 4. Manning map of the 8 x 4 m channel of inputs_manning: one value per
    #    row of 1 m (south to north 0.02, 0.03, 0.045, 0.06), columns of 2 m.
    rows = [0.02, 0.03, 0.045, 0.06]
    write(out / "manning_rows.asc", 4, 4, (0.0, 0.0), (2.0, 1.0), [[r] * 4 for r in rows])
    write(out / "manning_uniform.asc", 4, 4, (0.0, 0.0), (2.0, 1.0), [[0.03] * 4 for _ in rows])
    #    Ties: rows of 0.1 m from y = -2, so every grid centre (0.5, 1.5, 2.5,
    #    3.5) lies on the edge between two raster rows. Computed from the
    #    part of the raster kept (from y = -0.1), the position of each centre
    #    in raster rows is a whole number less a few epsilons, which a plain
    #    floor would put in the row below. The row above each centre holds
    #    the coefficient of that grid row, as do the other rows of its band;
    #    the row below holds 0.09, and so do the rows south of the domain.
    vals = []
    for j in range(60):
        k = (j - 20) // 10
        vals.append([0.09 if j < 20 or j == 24 + 10 * k else rows[k]] * 4)
    write(out / "manning_ties.asc", 4, 60, (0.0, -2.0), (2.0, 0.1), vals)

    # 5. Files that must abort, each on the streets domain [0, 40]^2.
    good = "ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 20\n"
    cases = {
        "unknown_key.asc": good.replace("cellsize", "cellsise") + "0 0\n0 0\n",
        "twice.asc": good + "ncols 2\n0 0\n0 0\n",
        "no_ncols.asc": good.replace("ncols 2\n", "") + "0 0\n0 0\n",
        "ncols_fraction.asc": good.replace("ncols 2", "ncols 2.5") + "0 0\n0 0\n",
        "corner_and_center.asc": good + "xllcenter 10\n0 0\n0 0\n",
        "cellsize_and_dx.asc": good + "dx 20\n0 0\n0 0\n",
        "cellsize_zero.asc": good.replace("cellsize 20", "cellsize 0") + "0 0\n0 0\n",
        "not_a_number.asc": good + "0 0\n0 x\n",
        "too_few.asc": good + "0 0\n0\n",
        "too_many.asc": good + "0 0\n0 0 0\n",
        "small.asc": good.replace("cellsize 20", "cellsize 19") + "0 0\n0 0\n",
        "short_y.asc": good.replace("cellsize 20", "dx 20\ndy 19") + "0 0\n0 0\n",
        "header_not_number.asc": good.replace("cellsize 20", "cellsize twenty") + "0 0\n0 0\n",
        "dx_no_dy.asc": good.replace("cellsize 20", "dx 20") + "0 0\n0 0\n",
        "y_corner_and_center.asc": good + "yllcenter 10\n0 0\n0 0\n",
        "no_nrows.asc": good.replace("nrows 2\n", "") + "0 0\n0 0\n",
        "extent_inf.asc": good.replace("cellsize 20", "cellsize 1e308") + "0 0\n0 0\n",
        "huge_header.asc": good.replace("ncols 2", "ncols 1073741824").replace("nrows 2", "nrows 1073741824")
        + "0 0\n0 0\n",
        "beyond_float.asc": good + "0 0\n1e39 0\n",
        "nodata.asc": good + "NODATA_value -9999\n0 0\n-9999 0\n",
        "negative.asc": good + "0 0\n-0.5 0\n",
    }
    for name, text in cases.items():
        (bad / name).write_text(text)


if __name__ == "__main__":
    main()
