#!/usr/bin/env python3
"""Midpoint reconnaissance for the seventh-moment {5,2} joint class.

This is a MODEL-SIDE NUMERICAL FALSIFIER, not a proof.

Definition:
    {5,2} = 7*(U1 + U2 + U3),
where Ud integrates the 7-walk overlap against C2(v) and the
partition-cyclic five-point cumulant C5.

The implementation is independently derived from
notes/spectator_52_spec.md.

Usage:
    python scripts/spectator_52_midpoint.py 0.25
    python scripts/spectator_52_midpoint.py 0.2 --sym

Fine grids are expensive (five dimensions, 150 C5 terms).  The point
of this script is to establish coarse rungs and symmetry/support gates
before an exact polytope implementation is attempted.
"""

import argparse
import itertools
import math
import os
import time

import numpy as np


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def compile_cumulant_terms(b=5):
    """Return (sign, prefix coefficient vectors) for C_b.

    Frequencies are (c1,...,c_{b-1}, cb) with
    cb = -sum_{j<b} cj.  Every subset sum is represented by an
    integer coefficient vector on the b-1 free variables.
    """
    out = []
    free = b - 1
    for part in set_partitions(list(range(b))):
        m = len(part)
        sign = -1.0 if (m - 1) % 2 else 1.0
        first = next(i for i, block in enumerate(part) if 0 in block)
        others = [i for i in range(m) if i != first]
        for perm in itertools.permutations(others):
            order = [first] + list(perm)
            prefixes = []
            current = set()
            for bi in order[:-1]:
                current.update(part[bi])
                coeff = [0] * free
                for idx in current:
                    if idx < free:
                        coeff[idx] += 1
                    else:
                        for j in range(free):
                            coeff[j] -= 1
                prefixes.append(tuple(coeff))
            out.append((sign, tuple(prefixes)))
    return out


TERMS5 = compile_cumulant_terms(5)
assert len(TERMS5) == 150


def overlap_from_positions(positions):
    mx = np.maximum.reduce(positions)
    mn = np.minimum.reduce(positions)
    return np.clip(1.0 - (mx - mn), 0.0, None)


def lincomb(coeff, C):
    out = np.zeros_like(C[0])
    for a, x in zip(coeff, C):
        if a:
            out = out + a * x
    return out


def c5_value(C):
    """Partition-cyclic five-point cumulant on four free arrays."""
    out = np.zeros_like(C[0])
    zero = np.zeros_like(C[0])
    for sign, prefixes in TERMS5:
        if not prefixes:
            out = out + sign
            continue
        positions = [zero]
        for coeff in prefixes:
            positions.append(lincomb(coeff, C))
        out = out + sign * overlap_from_positions(positions)
    return out


def walks(v, c1, c2, c3, c4):
    z = np.zeros_like(c1)
    V = np.full_like(c1, v)

    # d=1: increments v,-v,c1,c2,c3,c4,c5
    w1 = [
        z,
        V,
        z,
        c1,
        c1 + c2,
        c1 + c2 + c3,
        c1 + c2 + c3 + c4,
    ]

    # d=2: increments v,c1,-v,c2,c3,c4,c5
    w2 = [
        z,
        V,
        V + c1,
        c1,
        c1 + c2,
        c1 + c2 + c3,
        c1 + c2 + c3 + c4,
    ]

    # d=3: increments v,c1,c2,-v,c3,c4,c5
    w3 = [
        z,
        V,
        V + c1,
        V + c1 + c2,
        c1 + c2,
        c1 + c2 + c3,
        c1 + c2 + c3 + c4,
    ]
    return w1, w2, w3


def run(dv, use_sym=False, checkpoint=True):
    # The overlap support is contained in [-2,2] in every free
    # frequency; midpoint grid matches the lower-order engines.
    g = np.arange(-2.0 + dv / 2.0, 2.0, dv)
    if use_sym:
        # Global negation makes positive/negative v slices equal.
        vgrid = g[g > 0]
        sym_factor = 2.0
    else:
        vgrid = g
        sym_factor = 1.0

    tag = f"/tmp/increasezerozeta_52_dv{dv}_sym{int(use_sym)}.npz"
    acc = np.zeros(3, dtype=float)
    done = 0
    if checkpoint and os.path.exists(tag):
        d = np.load(tag)
        acc = d["acc"].astype(float)
        done = int(d["done"])

    t0 = time.time()
    # Chunk by v and c1; each mesh is only 3D.
    C2, C3, C4 = np.meshgrid(g, g, g, indexing="ij")
    cell_weight = dv ** 5

    total_slices = len(vgrid) * len(g)
    slice_no = 0
    for iv, v in enumerate(vgrid):
        for i1, c1_scalar in enumerate(g):
            if slice_no < done:
                slice_no += 1
                continue
            C1 = np.full_like(C2, c1_scalar)

            # C5 depends on c5=-c1-c2-c3-c4 internally through the
            # coefficient representation used by c5_value.
            c5 = c5_value([C1, C2, C3, C4])
            pair = min(abs(float(v)), 1.0)
            density = c5 * pair

            w1, w2, w3 = walks(v, C1, C2, C3, C4)
            acc[0] += float(np.sum(overlap_from_positions(w1) * density))
            acc[1] += float(np.sum(overlap_from_positions(w2) * density))
            acc[2] += float(np.sum(overlap_from_positions(w3) * density))

            slice_no += 1
            if checkpoint:
                np.savez(tag, acc=acc, done=slice_no)
            if slice_no % max(1, total_slices // 20) == 0:
                print(
                    f"{slice_no}/{total_slices} slices, "
                    f"partial={sym_factor * 7 * acc.sum() * cell_weight:+.8f}, "
                    f"wall={time.time()-t0:.1f}s",
                    flush=True,
                )

    U = sym_factor * acc * cell_weight
    total = 7.0 * float(U.sum())
    print(f"dv={dv}  U1={U[0]:+.10f} U2={U[1]:+.10f} U3={U[2]:+.10f}")
    print(f"{{5,2}} = 7*(U1+U2+U3) = {total:+.10f}")
    print(f"wall={time.time()-t0:.1f}s")
    return U, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dv", type=float)
    ap.add_argument("--sym", action="store_true")
    ap.add_argument("--no-checkpoint", action="store_true")
    args = ap.parse_args()
    if args.dv <= 0 or args.dv > 1:
        raise SystemExit("choose 0 < dv <= 1")
    run(args.dv, use_sym=args.sym, checkpoint=not args.no_checkpoint)


if __name__ == "__main__":
    main()
