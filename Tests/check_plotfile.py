#!/usr/bin/env python3
"""Check the plotfiles of a FLUME run against the numbers a test pins.

    check_plotfile.py <prefix> --steps N... --times T... --n-cell NX NY
        --prob-lo X Y --prob-hi X Y --dx DX DY --max-grid-size N --n-boxes N
        --n-files N [--bed-elevation E] [--bed-slope SX SY] [--surface ETA]
        [--bed-noise A] [--velocity U V] [--tol-eps N] [--rtol R] [--time-rtol R]

The run directory must hold exactly the plotfiles <prefix><step, 5 digits> of
the steps listed, and no other of that prefix. Each is read (text headers and
binary cell data; no AMReX or numpy needed) and checked:

  * the five variables, one level, three dimensions, the step and its time
    (within R times the time; exactly for time 0);
  * the domain is NX x NY x 1 cells and spans prob_lo to prob_hi in x and y,
    within R times the extent prob_hi - prob_lo (not R times the coordinate, so
    a shifted corner is seen at large coordinates), and exactly 0 to 1 in z;
  * the cell size in x and y, within R times itself, and exactly 1 in z;
  * the boxes: their number, that none is larger than max_grid_size in x or y,
    that each is one cell thick and inside the domain, and that together they
    hold exactly NX * NY cells without overlap;
  * the number of data files the boxes are written to, which is the number of
    MPI ranks that own a box;
  * every cell of every variable against the closed form
        bed_elevation     = E + SX (i + 1/2) DX + SY (j + 1/2) DY + A r(i, j)
        depth             = max(ETA - bed_elevation, 0)   (0 without --surface)
        discharge_x, _y   = depth U, depth V
        surface_elevation = bed_elevation + depth
    within N machine epsilons (of the precision the data are written in;
    N = 4 unless --tol-eps is given) times the larger of |bed_elevation| and
    |ETA|, times max(1, |U|, |V|) for the discharges. r(i, j) is the rough-bed
    number of the solver (bed_noise below; Source/FLUME_Init.cpp holds the same
    hash). The check also prints the largest difference found.

The closed form is the initial state. It is also the state at any later time
when the flow is uniform in a periodic direction or at rest under a level
surface, which is what the tests that use this check run.
"""
import argparse
import re
import struct
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


def bed_noise(i, j):
    """The rough-bed number of cell (i, j), in [-1, 1): add_bed_noise() in Source/FLUME_Init.cpp."""
    m = 0xFFFFFFFF
    n = ((i * 73856093) & m) ^ ((j * 19349663) & m)
    n ^= n >> 13
    n = (n * 0x5BD1E995) & m
    n ^= n >> 15
    return ((n & 2047) - 1024) / 1024


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
    """Boxes, and for each box the data file and the offset of its data in it."""
    lines = path.read_text().splitlines()
    start = next(k for k, s in enumerate(lines) if re.match(r"^\(\d+ \d+$", s))
    nbox = int(lines[start][1:].split()[0])
    boxes = []
    for s in lines[start + 1:start + 1 + nbox]:
        m = BOX.fullmatch(s.strip())
        if not m:
            sys.exit(f"{path}: cannot read a box from {s!r}")
        v = [int(x) for x in m.groups()]
        boxes.append((v[0:3], v[3:6]))
    fab_on_disk = []
    for s in lines[start + 1 + nbox:]:
        m = re.match(r"^FabOnDisk: (\S+) (\d+)$", s)
        if m:
            fab_on_disk.append((m.group(1), int(m.group(2))))
    if len(fab_on_disk) != nbox:
        sys.exit(f"{path}: {nbox} boxes but {len(fab_on_disk)} FabOnDisk lines")
    return boxes, fab_on_disk


FAB_HEAD = re.compile(rb"FAB \(\((\d+), \(([\d ]+)\)\),\((\d+), \(([\d ]+)\)\)\)"
                      rb"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)\) (\d+)\n")
VARS = ["depth", "discharge_x", "discharge_y", "bed_elevation", "surface_elevation"]


def read_cell_data(level_dir, fab_on_disk):
    """Cell values {(i, j): [one value per component]} and the size of a real, in bytes."""
    cells = {}
    real_size = None
    for name, offset in fab_on_disk:
        with open(level_dir / name, "rb") as f:
            f.seek(offset)
            head = f.readline()
            m = FAB_HEAD.fullmatch(head)
            if not m:
                sys.exit(f"{level_dir / name}: cannot read the FAB header {head!r}")
            size = int(m.group(3))
            order = [int(x) for x in m.group(4).split()]
            if size not in (4, 8) or order not in (list(range(size, 0, -1)), list(range(1, size + 1))):
                sys.exit(f"{level_dir / name}: unknown real format in {head!r}")
            real_size = size
            lo = [int(m.group(k)) for k in (5, 6, 7)]
            hi = [int(m.group(k)) for k in (8, 9, 10)]
            ncomp = int(m.group(14))
            ni, nj, nk = (hi[d] - lo[d] + 1 for d in range(3))
            n = ni * nj * nk
            raw = f.read(n * ncomp * size)
            if len(raw) != n * ncomp * size:
                sys.exit(f"{level_dir / name}: the data of box {lo} {hi} are cut short")
            endian = "<" if order[0] == size else ">"
            v = struct.unpack(f"{endian}{n * ncomp}{'d' if size == 8 else 'f'}", raw)
            for j in range(nj):
                for i in range(ni):
                    # x runs fastest; one cell in z
                    cells[(lo[0] + i, lo[1] + j)] = [v[c * n + j * ni + i] for c in range(ncomp)]
    return cells, real_size


def check_one(a, c, plotfile, step, time):
    for f in (plotfile / "Header", plotfile / "Level_0" / "Cell_H"):
        if not f.is_file():
            c.errors.append(f"{f} does not exist")
            return 0.0
    h = read_header(plotfile / "Header")
    nx, ny = a.n_cell
    c.eq("variables", h["vars"], VARS)
    c.eq("dimensions", h["dim"], 3)
    c.eq("finest level", h["finest_level"], 0)
    time_rtol = a.rtol if a.time_rtol is None else a.time_rtol
    c.near("time", [h["time"]], [time], [time_rtol * abs(time)])
    c.eq("step", h["step"], step)
    c.eq("domain lower corner", h["domain_lo"], [0, 0, 0])
    c.eq("domain upper corner", h["domain_hi"], [nx - 1, ny - 1, 0])
    extent_tol = [a.rtol * abs(hi - lo) for lo, hi in zip(a.prob_lo, a.prob_hi)]
    c.near("prob_lo (x, y)", h["prob_lo"][0:2], a.prob_lo, extent_tol)
    c.near("prob_hi (x, y)", h["prob_hi"][0:2], a.prob_hi, extent_tol)
    c.near("cell size (x, y)", h["dx"][0:2], a.dx, [a.rtol * abs(d) for d in a.dx])
    c.eq("z extent and cell size", [h["prob_lo"][2:], h["prob_hi"][2:], h["dx"][2:]],
         [[0.0], [1.0], [1.0]])

    boxes, fab_on_disk = read_cell_h(plotfile / "Level_0" / "Cell_H")
    files = {name for name, _ in fab_on_disk}
    c.eq("number of boxes", len(boxes), a.n_boxes)
    c.eq("number of data files (ranks that own a box)", len(files), a.n_files)
    for f in sorted(files):
        if not (plotfile / "Level_0" / f).is_file():
            c.errors.append(f"data file Level_0/{f} does not exist")
            return 0.0
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
    if c.errors:
        return 0.0

    cells, real_size = read_cell_data(plotfile / "Level_0", fab_on_disk)
    c.eq("cells in the data files", len(cells), nx * ny)
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    eta = a.surface
    u, v = a.velocity
    vscale = max(1.0, abs(u), abs(v))
    worst = 0.0
    n_bad = 0
    for (i, j), got in sorted(cells.items()):
        zb = a.bed_elevation + a.bed_slope[0] * (i + 0.5) * a.dx[0] + a.bed_slope[1] * (j + 0.5) * a.dx[1]
        if a.bed_noise != 0.0:
            zb += a.bed_noise * bed_noise(i, j)
        depth = 0.0 if eta is None else max(eta - zb, 0.0)
        want = [depth, depth * u, depth * v, zb, zb + depth]
        scale = max(abs(zb), 0.0 if eta is None else abs(eta))
        tol = [a.tol_eps * eps * scale * s for s in (1.0, vscale, vscale, 1.0, 1.0)]
        for name, g, w, t in zip(VARS, got, want, tol):
            # "not <=" also catches a NaN
            if not abs(g - w) <= t:
                n_bad += 1
                if n_bad <= 5:
                    c.errors.append(f"{name} in cell ({i}, {j}): got {g!r}, expected {w!r} (tolerance {t:.3g})")
            else:
                worst = max(worst, abs(g - w))
    if n_bad > 5:
        c.errors.append(f"... and {n_bad - 5} more cell values")
    return worst


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("prefix")
    p.add_argument("--steps", type=int, nargs="+", required=True)
    p.add_argument("--times", type=float, nargs="+", required=True)
    p.add_argument("--n-cell", type=int, nargs=2, required=True)
    p.add_argument("--prob-lo", type=float, nargs=2, required=True)
    p.add_argument("--prob-hi", type=float, nargs=2, required=True)
    p.add_argument("--dx", type=float, nargs=2, required=True)
    p.add_argument("--max-grid-size", type=int, required=True)
    p.add_argument("--n-boxes", type=int, required=True)
    p.add_argument("--n-files", type=int, required=True)
    p.add_argument("--bed-elevation", type=float, default=0.0)
    p.add_argument("--bed-slope", type=float, nargs=2, default=[0.0, 0.0])
    p.add_argument("--bed-noise", type=float, default=0.0)
    p.add_argument("--surface", type=float, default=None)
    p.add_argument("--velocity", type=float, nargs=2, default=[0.0, 0.0])
    p.add_argument("--tol-eps", type=float, default=4.0,
                   help="tolerance of the cell values, in machine epsilons of the scale (default 4)")
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--time-rtol", type=float, default=None,
                   help="tolerance of the time only (default: --rtol)")
    a = p.parse_args()
    if len(a.steps) != len(a.times):
        sys.exit("check_plotfile.py: --steps and --times need the same number of values")

    c = Checker()
    want = sorted(f"{a.prefix}{n:05d}" for n in a.steps)
    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    c.eq("plotfiles written", found, want)
    worst = 0.0
    for step, time in zip(a.steps, a.times):
        name = f"{a.prefix}{step:05d}"
        if name not in found:
            continue
        one = Checker()
        worst = max(worst, check_one(a, one, Path(name), step, time))
        c.errors += [f"{name}: {e}" for e in one.errors]

    if c.errors:
        print(f"FAIL: {a.prefix}")
        for e in c.errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {', '.join(want)}: {a.n_cell[0]} x {a.n_cell[1]} cells, every cell matches the"
          f" closed form (largest difference {worst:.3g})")


if __name__ == "__main__":
    main()
