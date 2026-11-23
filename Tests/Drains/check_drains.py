#!/usr/bin/env python3
"""Check FLUME runs with storm-drain inlets against closed-form numbers.

    check_drains.py <prefix> MODE [options]

Every plotfile <prefix><step> in the directory is read (at least N with
--n-plotfiles). The plotfiles must hold the fields of the state, the bed and
the surface, then drained_depth. Every mode checks, in every plotfile:

  * every depth is >= 0;
  * drained_depth is 0 in every cell that holds no inlet (--inlet I J, one
    per inlet cell) and >= 0 in the inlet cells;
  * the domain is closed and has no rain, so the exact sum over the cells
    of depth plus drained_depth is that of the first plotfile, within
    --ledger-eps epsilons of it: each cell keeps its own depth plus
    drained_depth to a few epsilons (the compensated sum of remove_water).

An inlet of perimeter P, open area A, weir coefficient C_w and orifice
coefficient C_o (--grate P A [C_w C_o]) in a cell of area S drains it at

    dh/dt = -min(a h^(3/2), b h^(1/2)),  a = C_w P / S,  b = C_o A sqrt(2 g) / S,

(FHWA HEC-22, 3rd edition, 2009, grates in a sag). From the depth h0 the
time to reach the depth h is computed here from the inverse of the closed
form, independently of the solver's forward form: above d* = b / a,
t = 2 (sqrt(h0) - sqrt(h)) / b; below it, t = t(d*) + 2 (1 / sqrt(h) -
1 / sqrt(h_s)) / a from the depth h_s the weir part starts at (h0 or d*).

Modes:

--uniform H0
    Every cell holds an inlet with the same grate, the water starts still
    at the depth H0 and stays level and still. Each cell must hold the
    depth whose time is that of the plotfile: the time error, times the
    drain rate at that depth, within --tol-eps epsilons of H0 (each step
    rounds to about an epsilon of the depth it starts from, so the error
    grows with the number of steps; it does not shrink with the depth).
    Discharges
    must be 0.

--basin H0 LO HI
    One inlet in a closed, flat basin of n cells, water still at H0 at the
    start. The basin drains as one cell of its area would, but the inlet
    cell is drawn below the mean depth, so it drains less: the drained
    volume (H0 minus the mean depth, per cell) must be within LO and HI
    times that of the lumped closed form (a and b over n), at every
    plotfile after the first.

--first-step H0
    The water starts still and level at the depth H0 (the film), so no water
    moves in the first step. In the first plotfile after the start, which
    must be that of step 1, each inlet cell (--inlet) must have drained the
    depth of the exact solution over that step: H0 minus drained_depth is
    checked through the inverse of the closed form, as in --uniform, but
    within --tol-eps epsilons of H0, since H0 minus drained_depth cancels. An
    explicit step, capped at the water of the cell, takes it all.
"""
import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_plotfile as cp  # noqa: E402

VARS = cp.VARS + ["drained_depth"]
I_H, I_QX, I_QY, I_D = 0, 1, 2, len(VARS) - 1


def read(plotfile):
    for f in (plotfile / "Header", plotfile / "Level_0" / "Cell_H"):
        if not f.is_file():
            sys.exit(f"FAIL: {f} does not exist")
    h = cp.read_header(plotfile / "Header")
    if h["vars"] != VARS:
        sys.exit(f"FAIL: {plotfile}: variables {h['vars']}, expected {VARS}")
    _, fab_on_disk = cp.read_cell_h(plotfile / "Level_0" / "Cell_H")
    cells, real_size = cp.read_cell_data(plotfile / "Level_0", fab_on_disk)
    nx = h["domain_hi"][0] + 1
    ny = h["domain_hi"][1] + 1
    if len(cells) != nx * ny:
        sys.exit(f"FAIL: {plotfile}: {len(cells)} cells, expected {nx * ny}")
    return h, cells, real_size


def time_to(h, h0, a, b):
    """Time for the depth to fall from h0 to h (0 < h <= h0)."""
    ds = b / a
    t = 0.0
    hs = h0
    if h0 > ds:
        if h >= ds:
            return 2.0 * (math.sqrt(h0) - math.sqrt(h)) / b
        t = 2.0 * (math.sqrt(h0) - math.sqrt(ds)) / b
        hs = ds
    return t + 2.0 * (1.0 / math.sqrt(h) - 1.0 / math.sqrt(hs)) / a


def rate(h, a, b):
    return min(a * h ** 1.5, b * math.sqrt(h))


def lumped_depth(t, h0, a, b):
    """Depth at time t from h0, by bisection on time_to (no forward form)."""
    lo, hi = 0.0, h0
    if t <= 0.0:
        return h0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= 0.0 or time_to(mid, h0, a, b) > t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("prefix")
    p.add_argument("--uniform", type=float, metavar="H0")
    p.add_argument("--basin", type=float, nargs=3, metavar=("H0", "LO", "HI"))
    p.add_argument("--grate", type=float, nargs="+", required=True,
                   help="P A [C_w C_o]")
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--first-step", type=float, metavar="H0")
    p.add_argument("--inlet", type=int, nargs=2, action="append", default=[], metavar=("I", "J"))
    p.add_argument("--all-inlets", action="store_true", help="every cell holds an inlet")
    p.add_argument("--ledger-eps", type=float, default=8.0)
    p.add_argument("--tol-eps", type=float, default=16.0)
    p.add_argument("--n-plotfiles", type=int, default=2)
    a = p.parse_args()
    if len(a.grate) not in (2, 4):
        sys.exit("FAIL: --grate takes P A or P A C_w C_o")
    perim, area = a.grate[0], a.grate[1]
    cw, co = (a.grate[2], a.grate[3]) if len(a.grate) == 4 else (1.66, 0.67)

    plotfiles = sorted((d for d in Path(".").iterdir()
                        if d.is_dir() and d.name.startswith(a.prefix)
                        and d.name[len(a.prefix):].isdigit()),
                       key=lambda d: int(d.name[len(a.prefix):]))
    if len(plotfiles) < a.n_plotfiles:
        sys.exit(f"FAIL: {len(plotfiles)} plotfiles {a.prefix}*, expected at least {a.n_plotfiles}")

    errors = []
    first_total = None
    worst = 0.0
    worst_ratio = [math.inf, -math.inf]
    inlets = {tuple(c) for c in a.inlet}
    first_step_done = False
    for pf in plotfiles:
        head, cells, real_size = read(pf)
        eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
        t = head["time"]
        s_cell = head["dx"][0] * head["dx"][1]
        ncell = len(cells)
        for (i, j), v in cells.items():
            if not v[I_H] >= 0.0:
                errors.append(f"{pf.name} cell ({i}, {j}): depth {v[I_H]!r} below 0")
            has_inlet = a.all_inlets or (i, j) in inlets
            if not has_inlet and v[I_D] != 0.0:
                errors.append(f"{pf.name} cell ({i}, {j}): drained_depth {v[I_D]!r} without an inlet")
            if has_inlet and not v[I_D] >= 0.0:
                errors.append(f"{pf.name} cell ({i}, {j}): drained_depth {v[I_D]!r} below 0")
        # Exact sums of the stored numbers.
        total = math.fsum(v[I_H] for v in cells.values()) + math.fsum(v[I_D] for v in cells.values())
        if first_total is None:
            first_total = total
        elif abs(total - first_total) > a.ledger_eps * eps * first_total:
            errors.append(f"{pf.name}: depth + drained_depth sums to {total!r}, "
                          f"at the start {first_total!r}")

        if a.uniform is not None:
            h0 = a.uniform
            ca = cw * perim / s_cell
            cb = co * area * math.sqrt(2.0 * a.gravity) / s_cell
            for (i, j), v in cells.items():
                h = v[I_H]
                if v[I_QX] != 0.0 or v[I_QY] != 0.0:
                    errors.append(f"{pf.name} cell ({i}, {j}): discharge {v[I_QX]!r} {v[I_QY]!r}, expected 0")
                if not h > 0.0:
                    continue
                dt_err = abs(time_to(h, h0, ca, cb) - t) if h < h0 else t
                err = dt_err * rate(h, ca, cb) / h0
                worst = max(worst, err / eps)
                if err > a.tol_eps * eps:
                    errors.append(f"{pf.name} cell ({i}, {j}) at t = {t}: depth {h!r} is "
                                  f"{err / eps:.1f} eps of H0 from the closed form")
        if a.first_step is not None and t > 0.0 and not first_step_done:
            first_step_done = True
            h0 = a.first_step
            ca = cw * perim / s_cell
            cb = co * area * math.sqrt(2.0 * a.gravity) / s_cell
            if head["step"] != 1:
                errors.append(f"{pf.name}: step {head['step']}, expected the plotfile of step 1")
            for (i, j) in inlets:
                d = cells[(i, j)][I_D]
                h1 = h0 - d
                if not (d > 0.0 and h1 > 0.0):
                    errors.append(f"{pf.name} cell ({i}, {j}): drained {d!r} of {h0} in the first step")
                    continue
                # Of H0: h1 is found as H0 minus the drained depth, which
                # cancels to an error of an epsilon of H0.
                err = abs(time_to(h1, h0, ca, cb) - t) * rate(h1, ca, cb) / h0
                worst = max(worst, err / eps)
                if err > a.tol_eps * eps:
                    errors.append(f"{pf.name} cell ({i}, {j}): drained {d!r} in the first step, "
                                  f"{err / eps:.1f} eps of H0 from the closed form")
        if a.basin is not None and t > 0.0:
            h0, lo, hi = a.basin
            ca = cw * perim / (s_cell * ncell)
            cb = co * area * math.sqrt(2.0 * a.gravity) / (s_cell * ncell)
            mean = math.fsum(v[I_H] for v in cells.values()) / ncell
            ratio = (h0 - mean) / (h0 - lumped_depth(t, h0, ca, cb))
            worst_ratio = [min(worst_ratio[0], ratio), max(worst_ratio[1], ratio)]
            if not lo <= ratio <= hi:
                errors.append(f"{pf.name} at t = {t}: drained volume {ratio:.6f} of the lumped "
                              f"closed form, expected {lo} to {hi}")

    if errors:
        for e in errors[:20]:
            print("FAIL:", e)
        sys.exit(1)
    msg = f"PASS: {len(plotfiles)} plotfiles"
    if a.first_step is not None and not first_step_done:
        sys.exit("FAIL: no plotfile after the start for --first-step")
    if a.uniform is not None or a.first_step is not None:
        msg += f"; largest depth error {worst:.2f} eps of H0"
    if a.basin is not None:
        msg += f"; drained {worst_ratio[0]:.6f} to {worst_ratio[1]:.6f} of the lumped closed form"
    print(msg)


if __name__ == "__main__":
    main()
