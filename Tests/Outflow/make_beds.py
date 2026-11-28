#!/usr/bin/env python3
"""Write the bed rasters of the outflow tests (ESRI ASCII grids).

    make_beds.py [<directory>]

bed_step.asc    50 x 2 cells of 2 m: the bed is at 0 m in the first column
                (the boundary cells of side xlo) and at 0.1 m in the others.
bed_valley.asc  100 x 2 cells of 1 m: the bed is 0.01 |x - 50 m|, so it
                crests at sides xlo and xhi and falls into the domain from
                both.

The files in the repository were written by this script; run it again only to
change them.
"""
import sys
from pathlib import Path


def write(path, nx, ny, dx, bed):
    with open(path, "w") as out:
        out.write(f"ncols {nx}\nnrows {ny}\nxllcorner 0.0\nyllcorner 0.0\ndx {dx}\ndy {dx}\n")
        for _ in range(ny):
            out.write(" ".join(f"{bed((i + 0.5) * dx, i):.3f}" for i in range(nx)) + "\n")


def main():
    where = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    write(where / "bed_step.asc", 50, 2, 2.0, lambda x, i: 0.0 if i == 0 else 0.1)
    write(where / "bed_valley.asc", 100, 2, 1.0, lambda x, i: 0.01 * abs(x - 50.0))


if __name__ == "__main__":
    main()
