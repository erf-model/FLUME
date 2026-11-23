#!/usr/bin/env python3
"""Check that the balance check finds each term of a log left out.

    check_dropped_terms.py [<log>] [--min-columns N]

Reads a water-balance log (default water_balance.txt) that passes
Tests/check_balance.py, then, for each column between 'stored' and
'residual', forms the residual a run would log if that column were missing
from the balance (a source that does not report its water, or an outflow
segment that is not counted) and checks that Tests/check_balance.py fails
on it: of the pairs of lines between which that column changes, at least
SHARE must exceed the bound of the check, so the term is found throughout the
run and not only where the bound is small (the first lines). Each such column
must move water in the run, and the log must have at least N of them
(default 1), so a run that exercises no term cannot pass.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_balance as cb  # noqa: E402

#: Share of the pairs of lines, among those between which the left-out
#: column changes, on which the check must fail.
SHARE = 0.9


def write_log(path, eps, names, rows, added, removed):
    """Write the header lines check() reads, with every column of removed
    under '# removed:' (the check does not tell removed from out)."""
    with open(path, "w") as f:
        f.write(f"# epsilon {eps!r}\n# added: {' '.join(added)}\n"
                f"# removed: {' '.join(removed)}\n# out:\n# {' '.join(names)}\n")
        for r in rows:
            f.write(" ".join(repr(v) for v in r) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("log", nargs="?", default="water_balance.txt")
    p.add_argument("--min-columns", type=int, default=1)
    a = p.parse_args()

    ok, msg, _ = cb.check(a.log)
    if not ok:
        sys.exit(f"FAIL: the log itself fails: {msg}")
    eps, names, rows = cb.read_log(a.log)
    added, removed = cb.read_terms(a.log)
    i_r = names.index("residual")
    terms = names[3:i_r]
    if sorted(terms) != sorted(added + removed):
        sys.exit(f"FAIL: the columns {terms} are not the terms added {added} and removed {removed}")
    if len(terms) < a.min_columns:
        sys.exit(f"FAIL: the log has the terms {terms}, fewer than {a.min_columns}")
    errors = []
    for name in terms:
        i = names.index(name)
        if all(r[i] == 0.0 for r in rows):
            errors.append(f"{name} moves no water in this run")
            continue
        # A run that does not count the term logs neither its column nor its
        # part of the residual, which subtracts the water added and adds the
        # water removed (both logged as positive volumes).
        sign = -1.0 if name in added else 1.0
        mutant = [r[:i] + r[i + 1:i_r] + [r[i_r] - sign * r[i]] + r[i_r + 1:] for r in rows]
        path = Path(f"dropped_{name}.txt")
        write_log(path, eps, names[:i] + names[i + 1:], mutant,
                  [c for c in added if c != name], [c for c in removed if c != name])
        ok, msg, _ = cb.check(path)
        m_eps, m_names, m_rows = cb.read_log(path)
        moved = failed = 0
        for (ra, rb, change, scale), a_row, b_row in zip(
                cb.scaled_changes(m_eps, m_names, m_rows), rows, rows[1:]):
            if b_row[i] != a_row[i]:
                moved += 1
                failed += change > cb.DEFAULT_FACTOR * m_eps * scale
        share = failed / moved
        print(f"without {name}: {'passes' if ok else 'fails'}: {msg}; fails on {failed} of the"
              f" {moved} pairs of lines where {name} changes")
        if ok:
            errors.append(f"the check passes the log without {name}")
        elif share < SHARE:
            errors.append(f"the log without {name} fails on {share:.1%} of the pairs where it"
                          f" changes, below {SHARE:.0%}")
    if errors:
        print("FAIL")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: the check finds each of {terms} left out")


if __name__ == "__main__":
    main()
