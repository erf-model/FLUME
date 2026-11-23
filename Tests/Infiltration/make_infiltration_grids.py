#!/usr/bin/env python3
"""Write the ESRI ASCII grids of the infiltration tests (Tests/Infiltration/*.asc).

    make_infiltration_grids.py [--dir DIR]

The files are committed; this script is how they were made, and running it
again writes them byte for byte. Every grid is synthetic: no data is
downloaded. See README.md for what each file tests.
"""
import argparse
from pathlib import Path

# Columns of 5 m across [0, 40] m, the same in every row.
# Saturated hydraulic conductivity, mm/h: impervious, then the classes of
# Rossman and Huber (2016, Table 4-7) from clay loam to sand (in/h x 25.4),
# and a second impervious column.
CONDUCTIVITY = [0.0, 1.016, 6.604, 10.922, 29.972, 120.396, 0.0, 3.302]
# Curve numbers: impervious, then 98 to 34.74 (65.31 and 34.74 are the
# lambda = 0.05 equivalents of 75 and 50 in Woodward et al., 2003, Table 4).
CURVE_NUMBER = [100.0, 98.0, 90.0, 75.0, 65.31, 50.0, 34.74, 85.0]


def write(path, nx, ny, origin, cell, values, nodata=None):
    """values[j][i]: column i from the west, row j from the south. The file
    lists the northern row first."""
    lines = [f"ncols {nx}", f"nrows {ny}", f"xllcorner {origin[0]!r}",
             f"yllcorner {origin[1]!r}", f"cellsize {cell!r}"]
    if nodata is not None:
        lines.append(f"NODATA_value {nodata!r}")
    for j in reversed(range(ny)):
        lines.append(" ".join(repr(v) for v in values[j]))
    path.write_text("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dir", type=Path, default=Path(__file__).resolve().parent)
    a = p.parse_args()
    d = a.dir
    (d / "bad").mkdir(parents=True, exist_ok=True)
    write(d / "conductivity_columns.asc", 8, 8, (0.0, 0.0), 5.0, [CONDUCTIVITY] * 8)
    write(d / "curve_number_columns.asc", 8, 8, (0.0, 0.0), 5.0, [CURVE_NUMBER] * 8)
    # For the 8 x 8 cells on [0, 40]^2 m of the abort tests: one quadrant of
    # 4 x 4 cells takes a value out of range.
    write(d / "bad" / "deficit_above_one.asc", 2, 2, (0.0, 0.0), 20.0, [[0.3, 0.3], [1.5, 0.3]])
    write(d / "bad" / "curve_number_zero.asc", 2, 2, (0.0, 0.0), 20.0, [[75.0, 75.0], [0.0, 75.0]])
    write(d / "bad" / "curve_number_above_100.asc", 2, 2, (0.0, 0.0), 20.0,
          [[75.0, 75.0], [101.0, 75.0]])


if __name__ == "__main__":
    main()
