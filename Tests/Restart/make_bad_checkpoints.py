#!/usr/bin/env python3
"""Write copies of a checkpoint that a restart must refuse.

Usage: make_bad_checkpoints.py <checkpoint>

In the current directory, each a copy of <checkpoint> with one change:

    bad_version    the version is 2
    bad_precision  real_bytes is that of the other precision
    bad_duplicate  the key "time" is given twice
    bad_step       the step is -1
    bad_nan        the depth of the first cell of the state file is NaN
"""
import re
import shutil
import struct
import sys
from pathlib import Path


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def copy(src, name):
    dst = Path(name)
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


def edit_header(d, pattern, repl):
    h = d / "Header"
    text = h.read_text()
    new, n = re.subn(pattern, repl, text, count=1, flags=re.M)
    if n != 1:
        fail(f"no line matching {pattern!r} in {h}")
    h.write_text(new)


def main():
    if len(sys.argv) != 2:
        fail("usage: make_bad_checkpoints.py <checkpoint>")
    src = Path(sys.argv[1])
    header = (src / "Header").read_text()
    m = re.search(r"^real_bytes (\d+)$", header, re.M)
    if not header.startswith("FLUME_checkpoint\n") or not m:
        fail(f"{src} is not a FLUME checkpoint")
    real_bytes = int(m.group(1))

    edit_header(copy(src, "bad_version"), r"^version 1$", "version 2")
    edit_header(copy(src, "bad_precision"), r"^real_bytes \d+$",
                f"real_bytes {4 if real_bytes == 8 else 8}")
    t = re.search(r"^time .*$", header, re.M).group(0)
    edit_header(copy(src, "bad_duplicate"), r"^time .*$", t + "\n" + t)
    edit_header(copy(src, "bad_step"), r"^step \d+$", "step -1")

    # The first value of the first FAB of the state: the depth of its first
    # cell. A FAB is a text line "FAB ..." and then the values.
    d = copy(src, "bad_nan")
    data = sorted(d.glob("state_D_*"))
    if not data:
        fail(f"no state_D_* file in {src}")
    raw = bytearray(data[0].read_bytes())
    start = raw.index(b"\n") + 1
    if not raw.startswith(b"FAB "):
        fail(f"{data[0]} does not start with a FAB header")
    nan = struct.pack("<d", float("nan")) if real_bytes == 8 else struct.pack("<f", float("nan"))
    raw[start:start + real_bytes] = nan
    data[0].write_bytes(raw)
    print("PASS: wrote bad_version bad_precision bad_duplicate bad_step bad_nan from", src)


if __name__ == "__main__":
    main()
