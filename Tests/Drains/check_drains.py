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
    With --basin the water moves, and every step rounds the depth of every
    wet cell by up to half an epsilon, with no preferred sign: the bound
    grows by sqrt(n) epsilons of the sum, n the steps from the first
    plotfile (the basin drifts by 17 epsilons over 6095 steps in single
    precision, 1 in double).

An inlet of perimeter P, open area A, weir coefficient C_w and orifice
coefficient C_o (--grate P A [C_w C_o]) in a cell of area S drains it at

    dh/dt = -min(a h^(3/2), b h^(1/2)),  a = C_w P / S,  b = C_o A sqrt(2 g) / S,

(FHWA HEC-22, 3rd edition, 2009, grates in a sag). From the depth h0 the
time to reach the depth h is computed here from the inverse of the closed
form, independently of the solver's forward form: above d* = b / a,
t = 2 (sqrt(h0) - sqrt(h)) / b; below it, t = t(d*) + 2 (1 / sqrt(h) -
1 / sqrt(h_s)) / a from the depth h_s the weir part starts at (h0 or d*).

An inlet of the approach-flow law of area A, coefficient k and Froude
exponent alpha (--approach A [k alpha]) in water of speed u (--velocity UX UY)
drains its cell at

    dh/dt = -c h^p,  c = k A g^p u^(1 + alpha) / S,  p = -alpha / 2,

(Xia et al. 2022 as given by Annis, Badas and Mascaro 2026, equation 3:
Q = k u A Fr^alpha, Fr = u / sqrt(g h)). From the depth h0 the time to reach
the depth h is, with q = 1 - p, t = (h0^q - h^q) / (q c): the cell is empty
after h0^q / (q c). With --approach, --uniform also checks that the discharge
of every cell is its depth times the velocity, that no cell is empty before
the closed form reaches the dry depth (--dry-depth, 1e-6 m by default: the
step takes a film at or below it), and that every cell is empty after the
emptying time.

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

--group A K ALPHA ILO IHI  (with --approach-first-step H0; one per group)
    The cells of columns ILO to IHI hold an approach-flow inlet of area A,
    coefficient K and Froude exponent ALPHA; the water starts level at H0
    and moves at --velocity, so no cell gains or loses water to another in
    the first step. In the plotfile of step 1, each cell of a group must
    hold the depth of the closed form of its own inlet, through the
    inverse, within --tol-eps epsilons of H0. Every cell must be in one
    group.

--still
    No inlet takes water: drained_depth is 0 in every cell of every
    plotfile, and every depth and discharge is that of the first plotfile.

--drained-positive
    Each inlet cell (--inlet) has drained more than 0 by the last plotfile.
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


def approach_c(area, k, alpha, speed, gravity, s_cell):
    return k * area * gravity ** (-0.5 * alpha) * speed ** (1.0 + alpha) / s_cell


def approach_time_to(h, h0, c, alpha):
    """Time for the depth to fall from h0 to h (0 <= h <= h0)."""
    q = 1.0 + 0.5 * alpha
    return (h0 ** q - h ** q) / (q * c)


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
    p.add_argument("--grate", type=float, nargs="+", help="P A [C_w C_o]")
    p.add_argument("--approach", type=float, nargs="+", help="A [k alpha]")
    p.add_argument("--velocity", type=float, nargs=2, default=[0.0, 0.0], metavar=("UX", "UY"))
    p.add_argument("--still", action="store_true")
    p.add_argument("--dry-depth", type=float, default=1.0e-6)
    p.add_argument("--approach-first-step", type=float, metavar="H0")
    p.add_argument("--group", type=float, nargs=5, action="append", default=[],
                   metavar=("A", "K", "ALPHA", "ILO", "IHI"))
    p.add_argument("--drained-positive", action="store_true")
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--first-step", type=float, metavar="H0")
    p.add_argument("--inlet", type=int, nargs=2, action="append", default=[], metavar=("I", "J"))
    p.add_argument("--all-inlets", action="store_true", help="every cell holds an inlet")
    p.add_argument("--ledger-eps", type=float, default=8.0)
    p.add_argument("--tol-eps", type=float, default=16.0)
    p.add_argument("--n-plotfiles", type=int, default=2)
    a = p.parse_args()
    if a.approach_first_step is not None:
        if a.grate is not None or a.approach is not None or not a.group:
            sys.exit("FAIL: --approach-first-step takes --group, not --grate or --approach")
        a.approach = [0.0]
    elif (a.grate is None) == (a.approach is None):
        sys.exit("FAIL: give one of --grate and --approach")
    if a.approach is not None:
        if len(a.approach) not in (1, 3):
            sys.exit("FAIL: --approach takes A or A k alpha")
        if a.basin is not None or a.first_step is not None:
            sys.exit("FAIL: --basin and --first-step are for --grate")
        ap_area = a.approach[0]
        ap_k, ap_alpha = (a.approach[1], a.approach[2]) if len(a.approach) == 3 else (0.302, -0.816)
        speed = math.hypot(a.velocity[0], a.velocity[1])
        perim = area = cw = co = None
    else:
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
    first_cells = None
    last_cells = None
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
        step = int(pf.name[len(a.prefix):])
        if first_total is None:
            first_total = total
            first_step = step
        else:
            ledger = a.ledger_eps + (math.sqrt(step - first_step) if a.basin is not None else 0.0)
            if abs(total - first_total) > ledger * eps * first_total:
                errors.append(f"{pf.name}: depth + drained_depth sums to {total!r}, "
                              f"at the start {first_total!r}")

        last_cells = cells
        if first_cells is None:
            first_cells = cells
        if a.still:
            for (i, j), v in cells.items():
                if v[I_D] != 0.0:
                    errors.append(f"{pf.name} cell ({i}, {j}): drained_depth {v[I_D]!r}, expected 0")
                w = first_cells[(i, j)]
                if (v[I_H], v[I_QX], v[I_QY]) != (w[I_H], w[I_QX], w[I_QY]):
                    errors.append(f"{pf.name} cell ({i}, {j}): depth and discharge "
                                  f"{v[I_H]!r} {v[I_QX]!r} {v[I_QY]!r} are not those of the "
                                  f"first plotfile")
        if a.uniform is not None and a.approach is not None:
            h0 = a.uniform
            c = approach_c(ap_area, ap_k, ap_alpha, speed, a.gravity, s_cell)
            p_exp = -0.5 * ap_alpha
            t_empty = approach_time_to(0.0, h0, c, ap_alpha)
            t_film = approach_time_to(min(a.dry_depth, h0), h0, c, ap_alpha)
            for (i, j), v in cells.items():
                h = v[I_H]
                # The discharge is the depth times the velocity: each step
                # scales it with the depth, rounding by about an epsilon.
                for q, u in ((v[I_QX], a.velocity[0]), (v[I_QY], a.velocity[1])):
                    if abs(q - h * u) > a.tol_eps * eps * h0 * max(speed, 1.0):
                        errors.append(f"{pf.name} cell ({i}, {j}): discharge {q!r}, expected "
                                      f"{h * u!r}")
                if not h > 0.0:
                    if t < t_film * (1.0 - 1e-9):
                        errors.append(f"{pf.name} cell ({i}, {j}) at t = {t}: empty before "
                                      f"{t_film!r} s")
                    continue
                if t > t_empty * (1.0 + 1e-9):
                    errors.append(f"{pf.name} cell ({i}, {j}) at t = {t}: depth {h!r} after "
                                  f"{t_empty!r} s, when the cell is empty")
                    continue
                dt_err = abs(approach_time_to(h, h0, c, ap_alpha) - t) if h < h0 else t
                err = dt_err * c * h ** p_exp / h0
                worst = max(worst, err / eps)
                if err > a.tol_eps * eps:
                    errors.append(f"{pf.name} cell ({i}, {j}) at t = {t}: depth {h!r} is "
                                  f"{err / eps:.1f} eps of H0 from the closed form")
        elif a.uniform is not None:
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
        if a.approach_first_step is not None and t > 0.0 and not first_step_done:
            first_step_done = True
            h0 = a.approach_first_step
            speed = math.hypot(a.velocity[0], a.velocity[1])
            if head["step"] != 1:
                errors.append(f"{pf.name}: step {head['step']}, expected the plotfile of step 1")
            for (i, j), v in cells.items():
                mine = [g for g in a.group if g[3] <= i <= g[4]]
                if len(mine) != 1:
                    errors.append(f"{pf.name} cell ({i}, {j}): in {len(mine)} groups")
                    continue
                g_area, g_k, g_alpha = mine[0][:3]
                c = approach_c(g_area, g_k, g_alpha, speed, a.gravity, s_cell)
                h = v[I_H]
                if not 0.0 < h < h0:
                    errors.append(f"{pf.name} cell ({i}, {j}): depth {h!r} after the first step "
                                  f"from {h0}")
                    continue
                err = abs(approach_time_to(h, h0, c, g_alpha) - t) * c * h ** (-0.5 * g_alpha) / h0
                worst = max(worst, err / eps)
                if err > a.tol_eps * eps:
                    errors.append(f"{pf.name} cell ({i}, {j}): depth {h!r} after the first step, "
                                  f"{err / eps:.1f} eps of H0 from the closed form of its inlet")
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

    if a.drained_positive:
        if not inlets:
            sys.exit("FAIL: --drained-positive needs --inlet")
        for (i, j) in inlets:
            if not last_cells[(i, j)][I_D] > 0.0:
                errors.append(f"{plotfiles[-1].name} cell ({i}, {j}): drained_depth "
                              f"{last_cells[(i, j)][I_D]!r}, expected above 0")
    if errors:
        for e in errors[:20]:
            print("FAIL:", e)
        sys.exit(1)
    msg = f"PASS: {len(plotfiles)} plotfiles"
    if (a.first_step is not None or a.approach_first_step is not None) and not first_step_done:
        sys.exit("FAIL: no plotfile after the start for --first-step")
    if a.uniform is not None or a.first_step is not None or a.approach_first_step is not None:
        msg += f"; largest depth error {worst:.2f} eps of H0"
    if a.basin is not None:
        msg += f"; drained {worst_ratio[0]:.6f} to {worst_ratio[1]:.6f} of the lumped closed form"
    print(msg)


if __name__ == "__main__":
    main()
