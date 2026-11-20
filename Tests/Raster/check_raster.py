#!/usr/bin/env python3
"""Check the bed a FLUME run built from ESRI ASCII grids.

    check_raster.py <prefix> [--plane A B C D] [--bed FILE]
        [--buildings FILE [--building-height H]] [--dry-buildings]
        [--wet-buildings] [--conserved] [--surface ETA] [--wet-beyond X]
        [--tol-eps K]

Every plotfile <prefix>NNNNN in the run directory is read (at least one).

  --plane A B C D   the bed in every cell is A + B x + C y + D x y at the cell
                    centre: bilinear interpolation reproduces this exactly.
  --bed FILE        the bed in every cell is FILE interpolated bilinearly
                    between its cell centres at the cell centre, the edge
                    values held within half a raster cell of its sides.
                    Computed here, in double precision, from the file.
  --buildings FILE  adds to the bed above, in every cell, the value of the
                    raster cell of FILE that holds the cell centre (a centre
                    on an edge, to within 16 epsilons of its position in
                    raster cells, takes the cell above it); with
                    --building-height H, adds H where that value is above 0.
  --dry-buildings   every cell with a building holds no water (depth exactly
                    0) in every plotfile, and the total depth is conserved
                    (as --conserved).
  --wet-buildings   in some plotfile some cell with a building holds more than
                    1 cm of water: the water went over a building.
  --conserved       the total depth is the same in all plotfiles within K
                    epsilons of the largest total (the sides are walls).
  --surface ETA     the lake at rest: in every plotfile every cell holds the
                    depth max(ETA - bed, 0) and discharges of 0, within K
                    epsilons of ETA plus the largest |bed| plus 1 m (depth)
                    and of that scale times sqrt(g ETA) (discharge, g = 9.81).
  --wet-beyond X    in the last plotfile some cell without a building whose
                    centre is beyond x = X holds more than 1 cm of water: the
                    water did reach the streets.

The bed is compared within K epsilons (--tol-eps, default 64) of the
largest |bed| plus 1 m, with the epsilon of the run's precision.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_dam_break as cd  # noqa: E402


def read_asc(path):
    """(header dict, values[j][i] with row 0 the southern one)."""
    tok = Path(path).read_text().split()
    head = {}
    k = 0
    while tok[k][0].isalpha():
        head[tok[k].lower()] = float(tok[k + 1])
        k += 2
    nx, ny = int(head["ncols"]), int(head["nrows"])
    dx = head.get("cellsize", head.get("dx"))
    dy = head.get("cellsize", head.get("dy"))
    x0 = head["xllcorner"] if "xllcorner" in head else head["xllcenter"] - dx / 2
    y0 = head["yllcorner"] if "yllcorner" in head else head["yllcenter"] - dy / 2
    vals = [float(t) for t in tok[k:]]
    assert len(vals) == nx * ny, path
    rows = [vals[r * nx:(r + 1) * nx] for r in range(ny)]
    rows.reverse()
    return dict(nx=nx, ny=ny, x0=x0, y0=y0, dx=dx, dy=dy), rows


def bilinear(r, v, x, y):
    def axis(p, p0, d, n):
        f = min(max((p - p0) / d - 0.5, 0.0), n - 1.0)
        i0 = min(int(math.floor(f)), max(n - 2, 0))
        return i0, min(i0 + 1, n - 1), f - i0
    i0, i1, wx = axis(x, r["x0"], r["dx"], r["nx"])
    j0, j1, wy = axis(y, r["y0"], r["dy"], r["ny"])
    return ((1 - wx) * (1 - wy) * v[j0][i0] + wx * (1 - wy) * v[j0][i1]
            + (1 - wx) * wy * v[j1][i0] + wx * wy * v[j1][i1])


def holding(q, n, eps):
    """The raster cell holding position q (in raster cells): a q within 16
    epsilons of a whole number is that number, as in FLUME_Raster.cpp."""
    whole = round(q)
    if abs(q - whole) <= 16.0 * eps * max(abs(q), 1.0):
        q = whole
    return min(max(int(math.floor(q)), 0), n - 1)


def nearest(r, v, x, y, eps):
    i = holding((x - r["x0"]) / r["dx"], r["nx"], eps)
    j = holding((y - r["y0"]) / r["dy"], r["ny"], eps)
    return v[j][i]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prefix")
    p.add_argument("--plane", type=float, nargs=4)
    p.add_argument("--bed")
    p.add_argument("--buildings")
    p.add_argument("--building-height", type=float)
    p.add_argument("--dry-buildings", action="store_true")
    p.add_argument("--wet-buildings", action="store_true")
    p.add_argument("--conserved", action="store_true")
    p.add_argument("--surface", type=float)
    p.add_argument("--wet-beyond", type=float)
    p.add_argument("--tol-eps", type=float, default=64.0)
    a = p.parse_args()
    if (a.plane is None) == (a.bed is None):
        sys.exit("FAIL: give exactly one of --plane and --bed")
    if (a.dry_buildings or a.wet_buildings or a.building_height is not None) and not a.buildings:
        sys.exit("FAIL: --dry-buildings, --wet-buildings and --building-height need --buildings")
    if a.dry_buildings and a.wet_buildings:
        sys.exit("FAIL: give at most one of --dry-buildings and --wet-buildings")

    found = sorted(d for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if not found:
        sys.exit(f"FAIL: no plotfile {a.prefix}NNNNN")
    runs = [cd.read(f) for f in found]
    head = runs[0][0]
    eps = 2.0 ** -52 if runs[0][2] == 8 else 2.0 ** -23
    lo, dx = head["prob_lo"], head["dx"]

    def centre(i, j):
        return lo[0] + (i + 0.5) * dx[0], lo[1] + (j + 0.5) * dx[1]

    if a.bed:
        br, bv = read_asc(a.bed)
    if a.buildings:
        hr, hv = read_asc(a.buildings)
    want, building = {}, {}
    for (i, j) in runs[0][1]:
        x, y = centre(i, j)
        if a.plane:
            A, B, C, D = a.plane
            z = A + B * x + C * y + D * x * y
        else:
            z = bilinear(br, bv, x, y)
        b = nearest(hr, hv, x, y, eps) if a.buildings else 0.0
        building[(i, j)] = b > 0.0
        if b > 0.0:
            z += a.building_height if a.building_height is not None else b
        want[(i, j)] = z
    scale = max(abs(z) for z in want.values()) + 1.0
    tol = a.tol_eps * eps * scale
    n_build = sum(building.values())
    print(f"{len(found)} plotfile(s), {len(want)} cells, {n_build} with a building,"
          f" bed tolerance {tol:.3e} m")

    errors = []
    worst = 0.0
    for f, (_, cells, _) in zip(found, runs):
        for c, v in sorted(cells.items()):
            err = abs(v[3] - want[c])
            worst = max(worst, err)
            if not err <= tol:
                errors.append(f"{f.name}, cell {c}: bed {v[3]!r}, expected {want[c]!r}")
            if a.dry_buildings and building[c] and v[0] != 0.0:
                errors.append(f"{f.name}, cell {c}: building cell holds depth {v[0]!r}")
    print(f"largest bed error {worst:.3e} m")

    if a.surface is not None:
        d_scale = a.surface + scale
        d_tol = a.tol_eps * eps * d_scale
        q_tol = d_tol * math.sqrt(9.81 * abs(a.surface))
        worst_h = worst_q = 0.0
        n_dry = sum(1 for c in want if want[c] >= a.surface)
        for f, (_, cells, _) in zip(found, runs):
            for c, v in sorted(cells.items()):
                eh = abs(v[0] - max(a.surface - want[c], 0.0))
                eq = max(abs(v[1]), abs(v[2]))
                worst_h, worst_q = max(worst_h, eh), max(worst_q, eq)
                if not eh <= d_tol:
                    errors.append(f"{f.name}, cell {c}: depth {v[0]!r}, at rest {max(a.surface - want[c], 0.0)!r}")
                if not eq <= q_tol:
                    errors.append(f"{f.name}, cell {c}: discharge {v[1]!r} {v[2]!r}")
        print(f"lake at rest: {n_dry} cell(s) above the surface; largest depth error {worst_h:.3e} m"
              f" (tolerance {d_tol:.3e}), largest discharge {worst_q:.3e} m^2/s (tolerance {q_tol:.3e})")
        if n_dry == 0 or n_dry == len(want):
            errors.append("no cell is dry, or none is wet: the test would not test the wet-dry fronts")
    if a.wet_buildings:
        wettest = max((v[0] for _, cells, _ in runs for c, v in cells.items() if building[c]), default=0.0)
        print(f"deepest water on a building cell: {wettest:.4f} m")
        if not wettest > 0.01:
            errors.append("no building cell holds more than 1 cm of water in any plotfile")
    if (a.dry_buildings or a.wet_buildings) and n_build == 0:
        errors.append("no cell has a building: the test would be vacuous")
    if a.dry_buildings or a.conserved:
        mass = [math.fsum(v[0] for v in cells.values()) for _, cells, _ in runs]
        m_tol = a.tol_eps * eps * max(mass) if mass else 0.0
        if len(mass) < 2:
            errors.append("--dry-buildings and --conserved need at least two plotfiles")
        elif not max(mass) - min(mass) <= m_tol:
            errors.append(f"total depth changes from {mass[0]!r} to {mass[-1]!r}, more than {m_tol:.3e}")
        else:
            print(f"total depth {mass[0]!r}, largest change {max(mass) - min(mass):.3e}")
    if a.wet_beyond is not None:
        last = runs[-1][1]
        deepest = max((v[0] for c, v in last.items()
                       if not building.get(c) and centre(*c)[0] > a.wet_beyond), default=0.0)
        print(f"deepest street cell beyond x = {a.wet_beyond}: {deepest:.4f} m")
        if not deepest > 0.01:
            errors.append(f"no street cell beyond x = {a.wet_beyond} holds more than 1 cm")

    if errors:
        for e in errors[:20]:
            print("FAIL: " + e)
        if len(errors) > 20:
            print(f"FAIL: ... and {len(errors) - 20} more")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
