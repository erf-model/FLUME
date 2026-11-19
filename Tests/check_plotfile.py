#!/usr/bin/env python3
"""Check a FLUME plotfile against the numbers a test pins.

    check_plotfile.py <plotfile> --n-cell NX NY --prob-lo X Y --prob-hi X Y
        --dx DX DY --max-grid-size N --n-boxes N --n-files N --var NAME --value V
        [--rtol R]

Reads the plotfile's Header and Level_0/Cell_H (text files; no AMReX or numpy
needed) and checks:

  * one variable with the given name, one level, three dimensions;
  * the domain is NX x NY x 1 cells and spans prob_lo to prob_hi in x and y,
    within R times the extent prob_hi - prob_lo (not R times the coordinate, so
    a shifted corner is seen at large coordinates), and exactly 0 to 1 in z;
  * the cell size in x and y, within R times itself, and exactly 1 in z;
  * the boxes: their number, that none is larger than max_grid_size in x or y,
    that each is one cell thick and inside the domain, and that together they
    hold exactly NX * NY cells without overlap;
  * the number of data files the boxes are written to, which is the number of
    MPI ranks that own a box;
  * the minimum and the maximum of the variable in every box equal V.
"""
import argparse
import re
import sys
from pathlib import Path

BOX = re.compile(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)\)")


class Checker:
    def __init__(self):
        self.errors = []

    def eq(self, what, got, want):
        if got != want:
            self.errors.append(f"{what}: got {got}, expected {want}")

    def near(self, what, got, want, tol):
        """Each got[d] within tol[d] of want[d]; no absolute floor."""
        if len(got) != len(want) or not all(abs(g - w) <= t for g, w, t in zip(got, want, tol)):
            self.errors.append(f"{what}: got {got}, expected {want}")


def read_header(path):
    lines = path.read_text().splitlines()
    h = {"version": lines[0].strip()}
    ncomp = int(lines[1])
    h["vars"] = [s.strip() for s in lines[2:2 + ncomp]]
    i = 2 + ncomp
    h["dim"] = int(lines[i])
    h["time"] = float(lines[i + 1])
    h["finest_level"] = int(lines[i + 2])
    h["prob_lo"] = [float(x) for x in lines[i + 3].split()]
    h["prob_hi"] = [float(x) for x in lines[i + 4].split()]
    # lines[i + 5]: refinement ratios (empty for one level)
    m = BOX.search(lines[i + 6])
    if not m:
        sys.exit(f"{path}: cannot read the domain box from {lines[i + 6]!r}")
    v = [int(x) for x in m.groups()]
    h["domain_lo"], h["domain_hi"] = v[0:3], v[3:6]
    h["step"] = int(lines[i + 7].split()[0])
    h["dx"] = [float(x) for x in lines[i + 8].split()]
    return h


def read_cell_h(path):
    """Boxes and per-box minimum and maximum of each component."""
    text = path.read_text()
    lines = text.splitlines()
    ncomp = int(lines[2])
    start = next(k for k, s in enumerate(lines) if re.match(r"^\(\d+ \d+$", s))
    nbox = int(lines[start][1:].split()[0])
    boxes = []
    for s in lines[start + 1:start + 1 + nbox]:
        m = BOX.fullmatch(s.strip())
        if not m:
            sys.exit(f"{path}: cannot read a box from {s!r}")
        v = [int(x) for x in m.groups()]
        boxes.append((v[0:3], v[3:6]))
    # Two blocks follow the FabOnDisk lines, each headed "<nbox>,<ncomp>":
    # the minima, then the maxima, one line per box.
    files = set()
    for s in lines[start + 1 + nbox:]:
        m = re.match(r"^FabOnDisk: (\S+) \d+$", s)
        if m:
            files.add(m.group(1))
    heads = [k for k, s in enumerate(lines) if s.strip() == f"{nbox},{ncomp}"]
    if len(heads) != 2:
        sys.exit(f"{path}: expected a minimum and a maximum block, found {len(heads)}")
    blocks = []
    for k in heads:
        rows = [[float(x) for x in s.strip().rstrip(",").split(",")] for s in lines[k + 1:k + 1 + nbox]]
        blocks.append(rows)
    return boxes, files, blocks[0], blocks[1]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plotfile", type=Path)
    p.add_argument("--n-cell", type=int, nargs=2, required=True)
    p.add_argument("--prob-lo", type=float, nargs=2, required=True)
    p.add_argument("--prob-hi", type=float, nargs=2, required=True)
    p.add_argument("--dx", type=float, nargs=2, required=True)
    p.add_argument("--max-grid-size", type=int, required=True)
    p.add_argument("--n-boxes", type=int, required=True)
    p.add_argument("--n-files", type=int, required=True)
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--var", required=True)
    p.add_argument("--value", type=float, required=True)
    a = p.parse_args()

    for f in (a.plotfile / "Header", a.plotfile / "Level_0" / "Cell_H"):
        if not f.is_file():
            sys.exit(f"FAIL: {f} does not exist")

    c = Checker()
    h = read_header(a.plotfile / "Header")
    nx, ny = a.n_cell
    c.eq("variables", h["vars"], [a.var])
    c.eq("dimensions", h["dim"], 3)
    c.eq("finest level", h["finest_level"], 0)
    c.eq("time", h["time"], 0.0)
    c.eq("step", h["step"], 0)
    c.eq("domain lower corner", h["domain_lo"], [0, 0, 0])
    c.eq("domain upper corner", h["domain_hi"], [nx - 1, ny - 1, 0])
    extent_tol = [a.rtol * abs(hi - lo) for lo, hi in zip(a.prob_lo, a.prob_hi)]
    c.near("prob_lo (x, y)", h["prob_lo"][0:2], a.prob_lo, extent_tol)
    c.near("prob_hi (x, y)", h["prob_hi"][0:2], a.prob_hi, extent_tol)
    c.near("cell size (x, y)", h["dx"][0:2], a.dx, [a.rtol * abs(d) for d in a.dx])
    c.eq("z extent and cell size", [h["prob_lo"][2:], h["prob_hi"][2:], h["dx"][2:]],
         [[0.0], [1.0], [1.0]])

    boxes, files, mins, maxs = read_cell_h(a.plotfile / "Level_0" / "Cell_H")
    c.eq("number of boxes", len(boxes), a.n_boxes)
    c.eq("number of data files (ranks that own a box)", len(files), a.n_files)
    for f in sorted(files):
        if not (a.plotfile / "Level_0" / f).is_file():
            c.errors.append(f"data file Level_0/{f} does not exist")
    covered = set()
    for lo, hi in boxes:
        size = [hi[d] - lo[d] + 1 for d in range(3)]
        if size[0] > a.max_grid_size or size[1] > a.max_grid_size:
            c.errors.append(f"box {lo} {hi} is larger than max_grid_size {a.max_grid_size}")
        if lo[2] != 0 or hi[2] != 0:
            c.errors.append(f"box {lo} {hi} is not one cell thick at k = 0")
        if lo[0] < 0 or lo[1] < 0 or hi[0] > nx - 1 or hi[1] > ny - 1:
            c.errors.append(f"box {lo} {hi} is outside the domain")
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                if (i, j) in covered:
                    c.errors.append(f"cell ({i}, {j}) is in two boxes")
                covered.add((i, j))
    c.eq("cells covered by the boxes", len(covered), nx * ny)

    for name, rows in (("minimum", mins), ("maximum", maxs)):
        for k, row in enumerate(rows):
            if row != [a.value]:
                c.errors.append(f"{name} of {a.var} in box {k}: got {row}, expected {[a.value]}")

    if c.errors:
        print(f"FAIL: {a.plotfile}")
        for e in c.errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {a.plotfile}: {nx} x {ny} cells in {len(boxes)} box(es), {a.var} = {a.value}")


if __name__ == "__main__":
    main()
