#!/usr/bin/env python3
"""The figure of the storm case, from the directory of its CTest run.

Usage: plot_storm.py <run directory> <output .png>

The run directory is that of erf_flume_storm (Tests/run/erf_flume_storm):
the coupled run's plt2d plotfiles and water_balance.txt, hyetograph.txt, and
the stand-alone runs in replay/ and hyeto/ (check_storm.py). Three panels:

(a) the rain ERF made in the coupled run (precip_total_accum at the end, mm,
    on ERF's columns of 2 km), the catchment, the town, the channel and the
    track of the cell (the centre of each 5 min of rain);
(b) the catchment-mean rain rate of hyetograph.txt;
(c) the discharge at the outlet of the coupled run (to its end, 4800 s), of
    its replay and of the hyetograph run.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from check_balance import read_log  # noqa: E402
from check_plotfile import read_cell_data, read_cell_h, read_header  # noqa: E402
from prepare_catchment import CATCHMENT, TOWN, VALLEY_AXIS, Y_DIVIDE, Y_LO  # noqa: E402

# The first three slots of the reference categorical palette (dataviz skill),
# which stay apart for readers with colour-vision deficiency.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#0b0b0b", "#52514e"


def load2d(path):
    p = Path(path)
    h = read_header(p / "Header")
    _, fab = read_cell_h(p / "Level_0" / "Cell_H")
    cells, _ = read_cell_data(p / "Level_0", fab)
    nx, ny = h["domain_hi"][0] + 1, h["domain_hi"][1] + 1
    a = np.zeros((nx, ny))
    for (i, j), v in cells.items():
        a[i, j] = v[0]
    return h, a


def outlet(where):
    _, names, rows = read_log(Path(where) / "water_balance.txt")
    t = np.array([r[names.index("time")] for r in rows])
    v = np.array([r[names.index("out_outlet")] for r in rows])
    return 0.5 * (t[1:] + t[:-1]) / 3600.0, np.diff(v) / np.diff(t)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    run = Path(sys.argv[1])
    plotfiles = sorted(run.glob("plt2d?????"))
    h, total = load2d(plotfiles[-1])
    lo, hi = h["prob_lo"], h["prob_hi"]
    km = 1e-3
    xe = np.linspace(lo[0], hi[0], total.shape[0] + 1) * km
    ye = np.linspace(lo[1], hi[1], total.shape[1] + 1) * km
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])

    # The track: the rain-weighted centre of each interval between plt2d files.
    track = []
    prev = None
    for p in plotfiles:
        _, a = load2d(p)
        if prev is not None:
            d = a - prev
            if d.sum() > 0.05 * total.max():
                track.append((np.sum(d.sum(axis=1) * xc) / d.sum(),
                              np.sum(d.sum(axis=0) * yc) / d.sum()))
        prev = a
    track = np.array(track)

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(10, 8.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1])
    ax = fig.add_subplot(gs[0, :])
    m = ax.pcolormesh(xe, ye, total.T, cmap="Blues", shading="flat")
    fig.colorbar(m, ax=ax, label="rain from ERF, mm", shrink=0.8, pad=0.01)
    def rect(box, **kw):
        ax.add_patch(plt.Rectangle((box[0] * km, box[2] * km), (box[1] - box[0]) * km,
                                   (box[3] - box[2]) * km, fill=False, **kw))
    rect(CATCHMENT, ec=INK, lw=1.6)
    rect(TOWN, ec=ORANGE, lw=1.6, ls="--")
    ax.plot([VALLEY_AXIS * km] * 2, [Y_LO * km, Y_DIVIDE * km], color=INK, lw=0.8, ls=":")
    ax.plot(track[:, 0], track[:, 1], color=INK, lw=1.0, marker="o", ms=3)
    ax.annotate("", xy=track[-1], xytext=track[-2],
                arrowprops={"arrowstyle": "->", "color": INK, "lw": 1.0})
    ax.text(CATCHMENT[0] * km, Y_DIVIDE * km + 0.5, "catchment (320 km²)", color=INK)
    ax.text(TOWN[0] * km + 0.3, TOWN[2] * km + 0.4, "town", color=ORANGE)
    ax.text(VALLEY_AXIS * km + 0.3, -15.3, "channel → outlet", color=INK, fontsize=8)
    ax.text(track[0, 0] - 13.0, track[0, 1] - 0.3, "storm track\n(every 5 min)", color=INK, fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("x, km")
    ax.set_ylabel("y, km")
    ax.set_title("(a) Rain of the convective cell at 80 min, coupled ERF + FLUME run",
                 loc="left", fontsize=10)

    hy = np.loadtxt(run / "hyetograph.txt")
    ax = fig.add_subplot(gs[1, 0])
    ax.step(hy[:, 0] / 60.0, hy[:, 1], where="post", color=BLUE, lw=1.5)
    ax.set_xlim(15, 80)
    ax.set_xlabel("time, min")
    ax.set_ylabel("catchment-mean rate, mm/h")
    ax.set_title("(b) Hyetograph: the same rain, spread over the catchment", loc="left",
                 fontsize=10)
    ax.grid(alpha=0.25)

    ax = fig.add_subplot(gs[1, 1])
    for where, label, color, kw in (("replay", "replay (gridded rain)", BLUE, {"lw": 1.5}),
                                    ("hyeto", "hyetograph run", ORANGE, {"lw": 1.5}),
                                    (".", "coupled run (to 80 min)", AQUA,
                                     {"lw": 0, "marker": "o", "ms": 3, "markevery": 60})):
        t, q = outlet(run / where)
        ax.plot(t, q, color=color, label=label, **kw)
    t_end = h["time"] / 3600.0
    ax.axvspan(0, t_end, color=AQUA, alpha=0.12, lw=0)
    ax.text(t_end / 2, 0.97, "coupled\nrun", transform=ax.get_xaxis_transform(), ha="center",
            va="top", color=MUTED, fontsize=8)
    ax.set_xlim(0, 8)
    ax.set_xlabel("time, h")
    ax.set_ylabel("outlet discharge, m³/s")
    ax.set_title("(c) Outlet hydrographs", loc="left", fontsize=10)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.savefig(sys.argv[2], dpi=110)
    print(f"plot_storm.py: wrote {sys.argv[2]}")


if __name__ == "__main__":
    main()
