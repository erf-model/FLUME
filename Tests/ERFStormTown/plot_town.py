#!/usr/bin/env python3
"""Draw the figure of the town case from the directory of its CTest run.

Usage:

    plot_town.py <run directory> <town directory> <figure.png>

(a) the flood depth at the end of the coupled run, with the buildings;
(b) ERF's wind speed at 15 m (the second level of cells, 10 to 20 m), with the solid
cells of terrain_IB_mask; (c) the rain rate over the square and the outflow
through the south side, with and without the buildings. Needs matplotlib
and numpy; not run by CTest.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_town import (ROOF, SLOPE_Y, SquareRain, flood_buildings, outflow,  # noqa: E402
                        read_plotfile, read_window)
from prepare_town import FLOOD_DX  # noqa: E402


def as_array(cells, comp, k=0):
    ni = 1 + max(i for i, _, _ in cells)
    nj = 1 + max(j for _, j, _ in cells)
    a = np.zeros((nj, ni))
    for (i, j, kk), v in cells.items():
        if kk == k:
            a[j, i] = v[comp]
    return a


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: plot_town.py <run directory> <town directory> <figure.png>")
    run, town, png = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    x0, y0, size, t0, window, _ = read_window(town)
    plts = sorted(run.glob("plt?????"))
    floods = sorted(run.glob("flood?????"))
    hf, fl = read_plotfile(floods[-1])
    he, erf = read_plotfile(plts[-1])
    depth = as_array(fl, hf["vars"].index("depth"))
    bld = flood_buildings(fl, hf)
    roof = np.zeros_like(depth)
    for (i, j), h in bld.items():
        roof[j, i] = h
    depth = np.ma.masked_where(roof > ROOF, depth)

    k = 1
    iu = [he["vars"].index(v) for v in ("x_velocity", "y_velocity")]
    u = as_array(erf, iu[0], k)
    v = as_array(erf, iu[1], k)
    mask = as_array(erf, he["vars"].index("terrain_IB_mask"), k)
    speed = np.ma.masked_where(mask >= 1.0, np.hypot(u, v))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    ext = [0, size, 0, size]
    im = ax[0].imshow(depth * 100, origin="lower", extent=ext, cmap="Blues", vmin=0)
    ax[0].imshow(np.ma.masked_where(roof <= ROOF, roof), origin="lower", extent=ext,
                 cmap="Greys", vmin=0, vmax=60)
    fig.colorbar(im, ax=ax[0], label="flood depth, cm", shrink=0.85)
    ax[0].set_title(f"(a) Flood at {hf['time']:g} s (buildings grey)")

    im = ax[1].imshow(speed, origin="lower", extent=ext, cmap="viridis")
    ax[1].imshow(np.ma.masked_where(mask < 1.0, mask), origin="lower", extent=ext,
                 cmap="Greys", vmin=0, vmax=2)
    n = u.shape[0]
    s = slice(2, n, 5)
    xc = (np.arange(n) + 0.5) * size / n
    ax[1].quiver(xc[s], xc[s], u[s, s], v[s, s], color="white", scale=50, width=0.004)
    fig.colorbar(im, ax=ax[1], label="wind speed, m/s", shrink=0.85)
    ax[1].set_title(f"(b) ERF's wind at {10 * k + 5} m, {he['time']:g} s")
    for a in ax[:2]:
        a.set_xlabel("x from the square's corner, m")
        a.set_ylabel("y, m")

    rain = SquareRain(town / "rain.txt", x0, y0, size)
    tr = np.array(rain.times)
    rr = np.array(rain.rates) * 3.6e6
    ax2 = ax[2].twinx()
    ax2.step(tr / 60, rr, where="post", color="0.6", lw=1)
    ax2.set_ylabel("rain over the square, mm/h", color="0.4")
    ax2.set_ylim(0, 2.2 * rr.max())
    ax2.invert_yaxis()
    for where, label, style in ((run, "with buildings (coupled)", "-"),
                                (run / "no_buildings", "without buildings", "--")):
        q, _, _ = outflow(where)
        t = np.array([p[0] for p in q])
        f = np.array([p[1] for p in q])
        w = 40  # a centred mean over 40 FLUME steps of 0.25 s, 10 s
        f = np.convolve(f, np.ones(w) / w, mode="valid")
        t = 0.5 * (t[w // 2 - 1:w // 2 - 1 + len(f)] + t[w // 2:w // 2 + len(f)])
        ax[2].plot(t / 60, f, style, color="C0", label=label)
    ax[2].set_xlabel("time from the window's start, min")
    ax[2].set_ylabel("outflow, south side, m3/s")
    ax[2].set_ylim(bottom=0)
    ax[2].legend(loc="upper left")
    ax[2].set_title("(c) Rain and outflow")
    fig.suptitle(f"Town of 100 buildings under the storm's rain: x {x0:g} to {x0 + size:g} m, "
                 f"y {y0:g} to {y0 + size:g} m, from {t0:g} s; bed rising {100 * SLOPE_Y:g} % "
                 f"to the north, flood cells of {FLOOD_DX:g} m")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
