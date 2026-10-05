#!/usr/bin/env python3
"""Tally chunked exact-polytope results for the {5,2} certification."""

import glob
import os
import re
import sys
from fractions import Fraction as F

TARGETS = {
    1: F(5, 504),
    2: F(1, 360),
    3: F(13, 2520),
}


def parse(paths):
    vals = {1: {}, 2: {}, 3: {}}
    pat = re.compile(
        r"52 d=(\d+) term (\d+) sign (-?\d+) value (\S+) wall"
    )
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                m = pat.search(line)
                if not m:
                    continue
                d = int(m.group(1))
                term = int(m.group(2))
                sign = int(m.group(3))
                val = F(m.group(4))
                signed = sign * val
                if term in vals[d] and vals[d][term] != signed:
                    raise RuntimeError(
                        f"conflicting duplicate d={d} term={term}: "
                        f"{vals[d][term]} vs {signed}"
                    )
                vals[d][term] = signed
    return vals


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    paths = sorted(glob.glob(os.path.join(root, "**", "*.txt"), recursive=True))
    if not paths:
        raise SystemExit(f"no result files found under {root!r}")

    vals = parse(paths)
    totals = {}
    all_complete = True
    all_match = True

    for d in (1, 2, 3):
        got = vals[d]
        missing = sorted(set(range(150)) - set(got))
        extra = sorted(set(got) - set(range(150)))
        total = sum(got.values(), F(0))
        totals[d] = total
        complete = not missing and not extra
        match = complete and total == TARGETS[d]
        all_complete &= complete
        all_match &= match

        print(
            f"d={d}: {len(got)}/150 terms  total={total} "
            f"= {float(total):+.15f}"
        )
        if missing:
            print("  missing:", missing)
        if extra:
            print("  extra:", extra)
        if complete:
            print(
                f"  target={TARGETS[d]} = {float(TARGETS[d]):+.15f} "
                f"{'MATCH' if match else 'MISMATCH'}"
            )

    if all_complete:
        joint = 7 * sum(totals.values(), F(0))
        print(f"{{5,2}} = 7*(U1+U2+U3) = {joint} = {float(joint):+.15f}")
        print(f"registered target = 1/8 = {float(F(1,8)):+.15f}")
        if joint != F(1, 8):
            all_match = False

    if not all_complete:
        raise SystemExit("INCOMPLETE exact term set")
    if not all_match:
        raise SystemExit("EXACT CANDIDATE MISMATCH")

    print("EXACT {5,2} CERTIFICATION PASSED")


if __name__ == "__main__":
    main()
