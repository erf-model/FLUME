#!/usr/bin/env python3
"""Check chosen cells of one FLUME plotfile against given numbers.

    check_cells.py <plotfile> --cell I J H HU HV [--cell ...] [--tol-eps K]
        [--others-unchanged <first plotfile>]

Each --cell gives the indices of a cell and the depth and the two discharges
it must hold, within K machine epsilons (of the precision of the data;
default 64) of the larger of 1 and the magnitude of the expected value. With
--others-unchanged every other cell must hold exactly what it holds in the
given plotfile. For runs whose first steps have a closed form.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_dam_break as cd  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("plotfile")
    p.add_argument("--cell", type=float, nargs=5, action="append", required=True)
    p.add_argument("--tol-eps", type=float, default=64.0)
    p.add_argument("--others-unchanged")
    a = p.parse_args()
    _, cells, real_size = cd.read(Path(a.plotfile))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    errors = []
    worst = 0.0
    listed = set()
    for i, j, *expected in a.cell:
        c = (int(i), int(j))
        listed.add(c)
        if c not in cells:
            sys.exit(f"FAIL: {a.plotfile} has no cell {c}")
        for name, got, want in zip(("depth", "discharge_x", "discharge_y"), cells[c], expected):
            off = abs(got - want) / max(1.0, abs(want))
            worst = max(worst, off)
            if not off <= a.tol_eps * eps:
                errors.append(f"cell {c}: {name} is {got!r}, expected {want!r}"
                              f" ({off / eps:.3g} epsilons, at most {a.tol_eps:g})")
    if a.others_unchanged:
        _, first, _ = cd.read(Path(a.others_unchanged))
        changed = [c for c, v in cells.items() if c not in listed and v[:3] != first[c][:3]]
        if changed:
            errors.append(f"{len(changed)} other cell(s) changed, e.g. {changed[0]}:"
                          f" {first[changed[0]][:3]} -> {cells[changed[0]][:3]}")
    if errors:
        print(f"FAIL: {a.plotfile}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {a.plotfile}: {len(a.cell)} cell(s) within {worst / eps:.3g} epsilons")


if __name__ == "__main__":
    main()
