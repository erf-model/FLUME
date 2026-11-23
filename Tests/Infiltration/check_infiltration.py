#!/usr/bin/env python3
"""Check FLUME runs with infiltration against closed-form numbers.

    check_infiltration.py <prefix> MODE [options]

Every plotfile <prefix><step> in the directory is read (at least N with
--n-plotfiles). The plotfiles must hold the fields of the state, the bed and
the surface, then infiltrated_depth. Rain rates are in mm/h, conductivity in
mm/h and suction in mm, as in the inputs. Tolerances are K epsilons (of the
precision the data are written in, --tol-eps) times the size of the number
compared, times the step number where the error can grow with the steps.

Modes:

--green-ampt K PSI THETA
    A closed, flat domain where every cell is the same: ponded water of depth
    H0 at t = 0 (--initial-depth) and no rain, or a dry start under constant
    rain R (--rain). With M = PSI THETA, the cumulative infiltration F(t) is
    (Rossman and Huber, 2016, section 4.4):
      * ponded from t = 0: K t = F - M ln(1 + F / M);
      * rain R <= K: F = R t;
      * rain R > K: F = R t up to t_p = F_s / R, F_s = K M / (R - K), then
        K (t - t_p) = (F - F_s) - M ln((F + M) / (F_s + M)).
    F is solved here by bisection, independently of the solver's Newton
    iteration. Every cell must hold infiltrated_depth F(t) and depth
    H0 + R t - F(t) within the tolerance (--rtol-f: an extra relative
    tolerance on F), and no discharge; with --velocity U V (a periodic
    domain whose water moves at a uniform velocity) the discharge h (U, V)
    instead: water that infiltrates takes its share of the discharge.

    --ledger-eps L checks each cell's own balance as well: depth plus
    infiltrated_depth equal to H0 + R t within L epsilons of H0 + R t, a
    tolerance that does not grow with the steps. Also in --curve-number
    with --flat (H0 = 0).

    --time-shift D allows F to lag or lead by D seconds of infiltration at
    the rate K (1 + M / F) as well. Ponded water of depth H0 is first
    infiltrated as rain of rate H0 / dt: in the first step the surface
    saturates at F_s = K M / (H0 / dt - K), reached after about
    K M dt^2 / H0^2 s; exact ponding reaches F_s in half that time, so F
    lags by about K M dt^2 / (2 H0^2) s. From there the integrated
    equation is exact, so F carries that shift of time, and only that.

--curve-number CN
    Constant rain R (--rain) from t = 0, P = R t. With S = 25.4 (1000 / CN -
    10) mm and Ia = lambda S (--ia-ratio, default 0.2), the cumulative loss is
    F = P for P <= Ia and P - (P - Ia)^2 / (P - Ia + S) above (Woodward et
    al., 2003, eq. 1). With --cn-map FILE each cell takes the curve number of
    the raster cell that holds its centre (CN is then ignored). Every cell
    must hold infiltrated_depth F within the tolerance. With --flat the
    domain is closed and flat and stays still: the depth must be P - F and
    the discharges 0. Without it water may run from cell to cell; the loss
    still depends on the rain of the cell only. With --impervious FILE the
    cells it marks (building footprints) must lose nothing. With
    --steady-tol the tolerance does not grow with the step number: the rain
    total and the infiltrated depth are compensated sums, so a run of many
    steps stays within a few epsilons.

--balance
    A closed domain: in every plotfile the sum over the cells of depth plus
    infiltrated_depth, minus the sum of depth in the first plotfile, is the
    number of cells times the rain fallen since t = 0 (--rain R, or
    --hyetograph FILE held between entries; none if neither). The sums are
    exact sums of the stored numbers. The infiltrated volume must be above
    0 in the last plotfile.

--cap
    Every plotfile after the first holds no water: depth exactly 0, no
    discharge, and infiltrated_depth exactly the depth of the cell in the
    first plotfile.

Checks for any mode:

--impervious FILE   Every cell whose raster cell in FILE (taking the cell
                    centre) is above 0 holds infiltrated_depth exactly 0.
--impervious-zero FILE
                    The same for the cells whose raster cell in FILE is 0
                    (a conductivity map).
--min-infiltrated V Some cell of the last plotfile holds at least V m.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_plotfile as cp  # noqa: E402

MM = 1.0e-3
MM_PER_H = 1.0 / 3.6e6
VARS = cp.VARS + ["infiltrated_depth"]
I_H, I_QX, I_QY, I_F = 0, 1, 2, len(VARS) - 1


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


def read_raster(path):
    """A raster as a function of a point: the value of the cell holding it."""
    head, values = {}, []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0][0].isalpha():
            head[parts[0].lower()] = float(parts[1])
        else:
            values.append([float(v) for v in parts])
    nx, ny = int(head["ncols"]), int(head["nrows"])
    x0, y0, cs = head["xllcorner"], head["yllcorner"], head["cellsize"]

    def at(x, y):
        i = min(int(math.floor((x - x0) / cs)), nx - 1)
        j = min(int(math.floor((y - y0) / cs)), ny - 1)
        return values[ny - 1 - j][i]
    return at


def centre(head, i, j):
    return (head["prob_lo"][0] + (i + 0.5) * head["dx"][0],
            head["prob_lo"][1] + (j + 0.5) * head["dx"][1])


def bisect(fun, lo, hi):
    """The root of the increasing fun in [lo, hi], to the last bit."""
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if fun(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def ponded(f1, k, m, tau):
    """F_2 - F_1 under ponding for the time tau, from F_1."""
    if tau <= 0.0 or k == 0.0:
        return 0.0
    if m == 0.0:
        return k * tau
    a = f1 + m
    # K tau <= x <= K tau + sqrt(2 M K tau)
    hi = k * tau + math.sqrt(2.0 * m * k * tau) * 1.01 + 1e-300
    return bisect(lambda x: x - m * math.log1p(x / a) - k * tau, 0.0, hi)


def green_ampt(t, k, m, rain, h0):
    if h0 > 0.0:
        return ponded(0.0, k, m, t)
    if rain <= k:
        return rain * t
    fs = k * m / (rain - k)
    tp = fs / rain
    if t <= tp:
        return rain * t
    return fs + ponded(fs, k, m, t - tp)


def cn_loss(p, cn, ratio):
    s = 25.4 * (1000.0 / cn - 10.0) * MM
    ia = ratio * s
    if p <= ia:
        return p
    return p - (p - ia) ** 2 / (p - ia + s)


def read_table(path):
    time, rate = [], []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].split()
        if line:
            time.append(float(line[0]))
            rate.append(float(line[1]) * MM_PER_H)
    return time, rate


def held_integral(time, rate, t):
    """The integral from 0 to t of the rate held from each entry to the next."""
    def from_first(tt):
        if tt <= time[0]:
            return rate[0] * (tt - time[0])
        total = 0.0
        for k in range(len(time) - 1):
            a, b = time[k], time[k + 1]
            if tt <= a:
                break
            total += rate[k] * (min(tt, b) - a)
        if tt > time[-1]:
            total += rate[-1] * (tt - time[-1])
        return total
    return from_first(t) - from_first(0.0)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("prefix")
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--green-ampt", nargs=3, type=float, metavar=("K", "PSI", "THETA"))
    mode.add_argument("--curve-number", type=float, metavar="CN")
    mode.add_argument("--balance", action="store_true")
    mode.add_argument("--cap", action="store_true")
    p.add_argument("--rain", type=float, default=0.0)
    p.add_argument("--hyetograph")
    p.add_argument("--initial-depth", type=float, default=0.0)
    p.add_argument("--ia-ratio", type=float, default=0.2)
    p.add_argument("--cn-map")
    p.add_argument("--flat", action="store_true")
    p.add_argument("--impervious")
    p.add_argument("--impervious-zero")
    p.add_argument("--min-infiltrated", type=float)
    p.add_argument("--tol-eps", type=float, default=16.0)
    p.add_argument("--rtol-f", type=float, default=0.0)
    p.add_argument("--time-shift", type=float, default=0.0)
    p.add_argument("--velocity", nargs=2, type=float, default=(0.0, 0.0), metavar=("U", "V"))
    p.add_argument("--n-plotfiles", type=int, default=2)
    p.add_argument("--ledger-eps", type=float)
    p.add_argument("--steady-tol", action="store_true")
    a = p.parse_args()

    found = sorted(d for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) < a.n_plotfiles:
        sys.exit(f"FAIL: {len(found)} plotfile(s) {a.prefix}NNNNN, expected at least {a.n_plotfiles}")
    rain = a.rain * MM_PER_H
    errors = []

    def fail(msg):
        if len(errors) < 20:
            errors.append(msg)
        else:
            errors.append(None)

    first_head, first, size = read(found[0])
    eps = 2.0 ** -52 if size == 8 else 2.0 ** -23

    def ledger(f, i, j, v, water_in):
        """The depth plus infiltrated_depth of the cell is the water it got,
        within --ledger-eps epsilons, however many steps the run took."""
        tol = a.ledger_eps * eps * water_in
        if abs(v[I_H] + v[I_F] - water_in) > tol:
            fail(f"{f.name} cell ({i}, {j}): depth {v[I_H]!r} + infiltrated {v[I_F]!r}"
                 f" = {v[I_H] + v[I_F]!r}, water in {water_in!r} (tol {tol:.3g})")
    data = [(f, *read(f)) for f in found]

    if a.green_ampt:
        k = a.green_ampt[0] * MM_PER_H
        m = a.green_ampt[1] * MM * a.green_ampt[2]
        worst = 0.0
        for f, head, cells, _ in data:
            t, step = head["time"], head["step"]
            want_f = green_ampt(t, k, m, rain, a.initial_depth)
            want_h = a.initial_depth + rain * t - want_f
            scale = a.initial_depth + rain * t
            growth = 1 if a.steady_tol else max(step, 1)
            tol = a.tol_eps * eps * scale * growth + a.rtol_f * want_f
            if want_f > 0.0:
                tol += a.time_shift * k * (1.0 + m / want_f)
            for (i, j), v in cells.items():
                if a.ledger_eps is not None:
                    ledger(f, i, j, v, scale)
                err = max(abs(v[I_F] - want_f), abs(v[I_H] - want_h))
                worst = max(worst, err)
                if err > tol:
                    fail(f"{f.name} cell ({i}, {j}): infiltrated {v[I_F]!r}, depth {v[I_H]!r};"
                         f" expected {want_f!r}, {want_h!r} (tol {tol:.3g})")
                u, w = a.velocity
                q_tol = 4.0 * eps * abs(v[I_H]) * max(abs(u), abs(w))
                if abs(v[I_QX] - v[I_H] * u) > q_tol or abs(v[I_QY] - v[I_H] * w) > q_tol:
                    fail(f"{f.name} cell ({i}, {j}): discharge ({v[I_QX]!r}, {v[I_QY]!r}),"
                         f" expected depth x ({u!r}, {w!r})")
            print(f"{f.name}: t = {t} s, F = {want_f * 1e3:.6f} mm, depth {want_h * 1e3:.6f} mm")
        print(f"largest difference {worst:.3g} m")

    elif a.curve_number is not None:
        cn_at = read_raster(a.cn_map) if a.cn_map else (lambda x, y: a.curve_number)
        building_at = read_raster(a.impervious) if a.impervious else (lambda x, y: 0.0)
        worst = 0.0
        for f, head, cells, _ in data:
            t, step = head["time"], head["step"]
            p_rain = rain * t
            growth = 1 if a.steady_tol else max(step, 1)
            tol = a.tol_eps * eps * p_rain * growth + a.rtol_f * p_rain
            for (i, j), v in cells.items():
                x, y = centre(head, i, j)
                want_f = 0.0 if building_at(x, y) > 0.0 else cn_loss(p_rain, cn_at(x, y), a.ia_ratio)
                err = abs(v[I_F] - want_f)
                if a.flat:
                    if a.ledger_eps is not None:
                        ledger(f, i, j, v, p_rain)
                    err = max(err, abs(v[I_H] - (p_rain - want_f)))
                    if v[I_QX] != 0.0 or v[I_QY] != 0.0:
                        fail(f"{f.name} cell ({i}, {j}): discharge ({v[I_QX]!r}, {v[I_QY]!r})")
                worst = max(worst, err)
                if err > tol:
                    fail(f"{f.name} cell ({i}, {j}): infiltrated {v[I_F]!r}, depth {v[I_H]!r};"
                         f" expected loss {want_f!r} of rain {p_rain!r} (tol {tol:.3g})")
            print(f"{f.name}: t = {t} s, rain {p_rain * 1e3:.6f} mm")
        print(f"largest difference {worst:.3g} m")

    elif a.balance:
        if a.hyetograph:
            time, rate = read_table(a.hyetograph)
            fallen = lambda t: held_integral(time, rate, t)  # noqa: E731
        else:
            fallen = lambda t: rain * t  # noqa: E731
        v0 = sum(v[I_H] for v in first.values())
        n = len(first)
        for f, head, cells, _ in data:
            t, step = head["time"], head["step"]
            water = sum(v[I_H] for v in cells.values())
            soil = sum(v[I_F] for v in cells.values())
            want = n * fallen(t)
            tol = a.tol_eps * eps * (v0 + want) * max(step, 1)
            print(f"{f.name}: t = {t} s: surface {water!r} + infiltrated {soil!r} - initial"
                  f" {v0!r} = {water + soil - v0!r}, rain {want!r} (tol {tol:.3g})")
            if abs(water + soil - v0 - want) > tol:
                fail(f"{f.name}: water + infiltrated - initial = {water + soil - v0!r}, rain {want!r}")
        if not sum(v[I_F] for v in data[-1][2].values()) > 0.0:
            fail("no water infiltrated")

    else:  # --cap
        for f, head, cells, _ in data[1:]:
            for (i, j), v in cells.items():
                h0 = first[(i, j)][I_H]
                if v[I_H] != 0.0 or v[I_QX] != 0.0 or v[I_QY] != 0.0 or v[I_F] != h0:
                    fail(f"{f.name} cell ({i}, {j}): depth {v[I_H]!r}, discharge"
                         f" ({v[I_QX]!r}, {v[I_QY]!r}), infiltrated {v[I_F]!r}; expected 0 and {h0!r}")
            print(f"{f.name}: t = {head['time']} s, no water left")

    for path, impervious in ((a.impervious, lambda v: v > 0.0),
                             (a.impervious_zero, lambda v: v == 0.0)):
        if not path:
            continue
        at = read_raster(path)
        n_imp = 0
        for f, head, cells, _ in data:
            for (i, j), v in cells.items():
                if impervious(at(*centre(head, i, j))):
                    n_imp += 1
                    if v[I_F] != 0.0:
                        fail(f"{f.name} cell ({i}, {j}): impervious, infiltrated {v[I_F]!r}")
        if n_imp == 0:
            fail(f"no cell is impervious in {path}")
        print(f"{n_imp} impervious cell values checked ({path})")

    if a.min_infiltrated is not None:
        top = max(v[I_F] for v in data[-1][2].values())
        print(f"largest infiltrated depth {top!r} m")
        if top < a.min_infiltrated:
            fail(f"largest infiltrated depth {top!r} below {a.min_infiltrated!r}")

    if errors:
        for e in errors:
            if e is not None:
                print("FAIL:", e)
        if None in errors:
            print("FAIL: ... more")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
