#!/usr/bin/env python3
"""The catchment-mean hyetograph of the rain FLUME took from ERF.

Reads the index of the coupled run's rain (flume.rain.replay_file) and its
ESRI ASCII grids, and writes a hyetograph file (flume.rain.file with
flume.rain.type = hyetograph): for each line of the index, its time and the
mean over the catchment's cells of the grid's rate, in mm/h. The catchment
is a rectangle of whole cells (its default is that of prepare_catchment.py).
The cells are of one size, so the mean rate times the catchment's area is
the rate of the grid summed over it: spread evenly, the hyetograph puts the
coupled run's rain on the catchment, step by step.

Usage: make_hyetograph.py <index> <output> [--catchment x0 x1 y0 y1]
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from prepare_catchment import CATCHMENT  # noqa: E402


def read_grid(path):
    """Header {key: value} and rows of values, the northern row first."""
    head = {}
    rows = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            if s[0][0].isalpha():
                head[s[0].lower()] = float(s[1])
            else:
                rows.append([float(x) for x in s])
    for k in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize"):
        if k not in head:
            sys.exit(f"make_hyetograph.py: {path}: no {k} in the header")
    if len(rows) != int(head["nrows"]) or any(len(r) != int(head["ncols"]) for r in rows):
        sys.exit(f"make_hyetograph.py: {path}: not {int(head['nrows'])} rows of {int(head['ncols'])} values")
    return head, rows


def catchment_mean(path, box):
    head, rows = read_grid(path)
    x0, y0, d = head["xllcorner"], head["yllcorner"], head["cellsize"]
    ny = len(rows)
    total, n = 0.0, 0
    for jr, row in enumerate(rows):
        y = y0 + (ny - jr - 0.5) * d
        if not box[2] <= y < box[3]:
            continue
        for i, v in enumerate(row):
            x = x0 + (i + 0.5) * d
            if box[0] <= x < box[1]:
                total += v
                n += 1
    if n == 0:
        sys.exit(f"make_hyetograph.py: {path}: no cell centre in the catchment {box}")
    return total / n, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index")
    ap.add_argument("output")
    ap.add_argument("--catchment", nargs=4, type=float, default=list(CATCHMENT),
                    metavar=("X0", "X1", "Y0", "Y1"))
    a = ap.parse_args()
    index = pathlib.Path(a.index)
    lines = []
    cells = None
    for raw in index.read_text().splitlines():
        s = raw.split("#")[0].split()
        if not s:
            continue
        if len(s) != 2:
            sys.exit(f"make_hyetograph.py: {index}: cannot read the line {raw!r}")
        rate, n = catchment_mean(index.parent / s[1], a.catchment)
        if cells is not None and n != cells:
            sys.exit(f"make_hyetograph.py: {s[1]}: {n} catchment cells, not {cells}")
        cells = n
        lines.append(f"{s[0]} {rate:.17g}")
    if not lines:
        sys.exit(f"make_hyetograph.py: {index}: no lines")
    with open(a.output, "w") as f:
        f.write(f"# The catchment-mean rain of {index.name} over {cells} cells, mm/h\n")
        f.write("# time (s)  rate (mm/h)\n")
        f.write("\n".join(lines) + "\n")
    print(f"make_hyetograph.py: {len(lines)} lines over {cells} cells in {a.output}")


if __name__ == "__main__":
    main()
