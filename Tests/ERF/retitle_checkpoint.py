#!/usr/bin/env python3
"""Copy an ERF checkpoint with another title line in its Header.

    retitle_checkpoint.py <checkpoint> <copy> <title>

ERF skips the title line when it reads the Header; FLUME reads the step of a
restart from that Header and stops on a title other than "Checkpoint file
for ERF" (abort_erf_flume_checkpoint_title).
"""
import shutil
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: retitle_checkpoint.py <checkpoint> <copy> <title>")
    src, dst, title = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    lines = (src / "Header").read_text().splitlines(keepends=True)
    if not lines or lines[0] != "Checkpoint file for ERF\n":
        sys.exit(f"FAIL: {src}/Header does not start with the line \"Checkpoint file for ERF\"")
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    (dst / "Header").write_text(title + "\n" + "".join(lines[1:]))
    print(f"{dst}: {src} with the title {title!r}")


if __name__ == "__main__":
    main()
