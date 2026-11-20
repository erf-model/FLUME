#!/usr/bin/env python3
"""Check that the plotfiles of a run are byte for byte those of another run.

Usage (in the directory of the run):
    compare_plotfiles.py <prefix> <reference directory> [--cells]

Every plotfile <prefix>NNNNN in this directory must also be in the reference
directory, and the other way round, with at least two of them, and every file
in each plotfile (header, box layout, data) must be identical to its copy.
Used where two decks must give the same answer to the last bit, for example
one boundary segment and the same cells split into two segments of one type.

With --cells the two runs may split the domain into different boxes, so the
files differ; then every value of every cell (read with
Tests/check_plotfile.py) and the time and step must be the same.
"""
import filecmp
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def plotfiles(where, prefix):
    pattern = re.escape(prefix) + r"\d{5,}"
    return sorted(d.name for d in Path(where).iterdir() if d.is_dir() and re.fullmatch(pattern, d.name))


def files_of(root):
    return sorted(str(p.relative_to(root)) for p in Path(root).rglob("*") if p.is_file())


def main():
    args = sys.argv[1:]
    cells = "--cells" in args
    args = [x for x in args if x != "--cells"]
    if len(args) != 2:
        sys.exit("usage: compare_plotfiles.py <prefix> <reference directory> [--cells]")
    prefix, ref = args[0], Path(args[1])
    if not ref.is_dir():
        sys.exit(f"FAIL: no reference directory {ref}")
    if ref.resolve() == Path(".").resolve():
        sys.exit(f"FAIL: the reference directory {ref} is this directory")
    here = plotfiles(".", prefix)
    there = plotfiles(ref, prefix)
    if here != there:
        sys.exit(f"FAIL: plotfiles here {here}, in {ref} {there}")
    if len(here) < 2:
        sys.exit(f"FAIL: {len(here)} plotfile(s); the comparison needs the first and a later one")
    if cells:
        from check_dam_break import read  # noqa: E402
        n = 0
        for name in here:
            h_mine, mine, _ = read(Path(name))
            h_theirs, theirs, _ = read(ref / name)
            if (h_mine["time"], h_mine["step"]) != (h_theirs["time"], h_theirs["step"]):
                sys.exit(f"FAIL: {name} is of step {h_mine['step']} at {h_mine['time']!r} here and of step"
                         f" {h_theirs['step']} at {h_theirs['time']!r} in {ref}")
            if mine.keys() != theirs.keys():
                sys.exit(f"FAIL: {name} holds other cells here than in {ref}")
            for c, v in mine.items():
                if list(v) != list(theirs[c]):
                    sys.exit(f"FAIL: {name} cell {c}: {list(v)} here, {list(theirs[c])} in {ref}")
                n += 1
        print(f"PASS: {len(here)} plotfiles, {n} cells, every value identical to those in {ref}")
        return
    n = 0
    for name in here:
        mine = files_of(Path(name))
        theirs = files_of(ref / name)
        if mine != theirs:
            sys.exit(f"FAIL: {name} holds {mine} here and {theirs} in {ref}")
        for f in mine:
            n += 1
            if not filecmp.cmp(Path(name) / f, ref / name / f, shallow=False):
                sys.exit(f"FAIL: {name}/{f} differs from its copy in {ref}")
    print(f"PASS: {len(here)} plotfiles, {n} files, identical to those in {ref}")


if __name__ == "__main__":
    main()
