#!/usr/bin/env python3
"""Integer lattice evaluation of the corrected model {5,2} integral.

See notes/spectator_52_exact_lattice.md for the degree and period argument.
Requires NumPy; all sums and interpolation use exact integer arithmetic.
This evaluates the stated model integral, not its transport to zeta moments.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import time

import numpy as np

from spectator_52_reduced_midpoint import TERMS
from spectator_reduction import spectator_integral


def cnum(xs, n):
    zero = np.zeros_like(xs[0], dtype=np.int32)
    forms = {}
    for _, prefixes in TERMS:
        for coeff in prefixes:
            if coeff not in forms:
                forms[coeff] = sum((a*x for a, x in zip(coeff, xs) if a), zero)
    result = np.zeros_like(zero)
    for sign, prefixes in TERMS:
        points = [zero] + [forms[c] for c in prefixes]
        span = np.maximum.reduce(points) - np.minimum.reduce(points)
        result += sign*np.maximum(n-span, 0)
    return result


def jnum(prefixes, d, n):
    fixed, shifted = [prefixes[0]] + prefixes[d-1:], prefixes[:d]
    m, M = np.minimum.reduce(fixed), np.maximum.reduce(fixed)
    a, A = np.minimum.reduce(shifted), np.maximum.reduce(shifted)
    value = (np.abs(M-a-n)**3 + np.abs(n+m-A)**3
             - np.abs(m-a)**3 - np.abs(M-A)**3)
    return np.where(np.maximum(M-m, A-a) < n, value, 0)


def key(xs, n):
    result = np.zeros_like(xs[0], dtype=np.int64)
    for x in xs[:4]:
        result = result*(2*n+1) + x+n
    return result


def cache(n):
    """C5 is permutation invariant; store each sorted zero-sum quintuple once."""
    table = np.full((2*n+1)**4, 32767, dtype=np.int16)
    chunks, size, count = [], 0, 0

    def flush():
        nonlocal chunks, size, count
        if not chunks:
            return
        rows = np.concatenate(chunks)
        values = cnum(list(rows.T), n)
        if np.any(np.abs(values) > 150*n):
            raise AssertionError("C5 numerator bound")
        table[key(list(rows.T), n)] = values
        count += len(rows)
        chunks, size = [], 0

    for a in range(-n, 1):
        for b in range(a, n+1):
            for c in range(b, n+1):
                lo, hi = max(c, -a-b-c-n), min(n, (-a-b-c)//2)
                if lo <= hi:
                    ds = np.arange(lo, hi+1, dtype=np.int32)
                    chunks.append(np.column_stack((np.full_like(ds, a),
                                                   np.full_like(ds, b),
                                                   np.full_like(ds, c), ds)))
                    size += len(ds)
                    if size >= 50000:
                        flush()
    flush()
    return table, count


def run(n, use_cache=True, symmetry=True):
    if not 1 <= n <= 60:
        raise ValueError("1 <= N <= 60 is required by the integer/memory bounds")
    # J <= integral_{-1}^1 |v| dv = 1, |Cnum| <= 150N.
    # Even the sum of absolute values fits signed int64 at these scales.
    bound = (2*n+1)**4 * 900*n**4
    assert bound < 2**63
    assert 4*(5*n)**3 < 2**31
    started = time.monotonic()
    table, count = cache(n) if use_cache else (None, None)
    grid = np.arange(-n, n+1, dtype=np.int32)
    c2, c3, c4 = np.meshgrid(grid, grid, grid, indexing="ij")
    zero = np.zeros_like(c2)
    totals, pure = [0, 0, 0], 0
    for a in (range(0, n+1) if symmetry else range(-n, n+1)):
        c1 = np.full_like(c2, a)
        xs = [c1, c2, c3, c4]
        last = -c1-c2-c3-c4
        valid = np.abs(last) <= n
        if table is None:
            value = cnum(xs, n)
        else:
            ordered = np.sort(np.stack(xs + [last]), axis=0)
            # Out-of-support rows have invalid keys; do not index them.
            value = np.zeros_like(c2)
            value[valid] = table[key(list(ordered[:, valid]), n)]
            if np.any(value[valid] == 32767):
                raise AssertionError("missing sorted tuple")
        prefixes = [zero, c1, c1+c2, c1+c2+c3, c1+c2+c3+c4]
        factor = 2 if symmetry and a else 1
        for d in (1, 2, 3):
            weight = jnum(prefixes, d, n)
            if np.any(weight[~valid]):
                raise AssertionError("support outside the five-frequency cube")
            if np.any(weight < 0) or np.any(weight > 6*n**3):
                raise AssertionError("spectator integral bound")
            totals[d-1] += factor*int(np.sum(weight.astype(np.int64)*value))
        overlap = np.maximum(n-np.maximum.reduce(prefixes)
                             + np.minimum.reduce(prefixes), 0)
        pure += factor*int(np.sum(overlap.astype(np.int64)*value))
    return {"N": n, "S": totals, "pure_S": pure,
            "sorted_tuples": count, "absolute_sum_bound": bound,
            "seconds": time.monotonic()-started}


def determinant(rows):
    return sum((-1)**sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4))
               * math.prod(rows[i][p[i]] for i in range(4))
               for p in itertools.permutations(range(4)))


def gates():
    normals = list(itertools.product((0, 1), repeat=4))[1:]
    maxdet = max(abs(determinant(rows))
                 for rows in itertools.combinations(normals, 4))
    assert maxdet == 3
    # Every difference in each nested C5 prefix chain is +/- binary.
    for _, prefixes in TERMS:
        for p, q in itertools.combinations(((0, 0, 0, 0),)+prefixes, 2):
            diff = {a-b for a, b in zip(p, q)}
            assert diff <= {0, 1} or diff <= {0, -1}
    checked = 0
    for n in (1, 2):
        direct, cached = run(n, False, False), run(n, True, True)
        assert direct["S"] == cached["S"] and direct["pure_S"] == cached["pure_S"]
        for row in itertools.product(range(-n, n+1), repeat=4):
            prefixes = [0] + list(itertools.accumulate(row))
            arrays = [np.array(x, dtype=np.int32) for x in prefixes]
            for d in (1, 2, 3):
                assert F(int(jnum(arrays, d, n)), 6*n**3) == spectator_integral(
                    tuple(F(x, n) for x in prefixes), d)
                checked += 1
    return {"maximum_binary_determinant": maxdet,
            "rational_spectator_cases": checked,
            "cache_and_negation_comparison_scales": [1, 2]}


def leading(values, degree, step=6):
    work = [int(v) for v in values]
    assert len(work) == degree+1
    for _ in range(degree):
        work = [b-a for a, b in zip(work, work[1:])]
    return F(work[0], math.factorial(degree)*step**degree)


def predict(values):
    work, edges = list(values), []
    while work:
        edges.append(work[-1])
        work = [b-a for a, b in zip(work, work[1:])]
    return sum(edges)


def certificate(records):
    by_n = {r["N"]: r for r in records}
    if not all(n in by_n for n in range(6, 61, 6)):
        return None
    us = []
    for d in range(3):
        samples = [by_n[n]["S"][d] for n in range(6, 55, 6)]
        assert predict(samples) == by_n[60]["S"][d], "degree-eight held-out failure"
        us.append(leading(samples, 8)/6)
    pure = [by_n[n]["pure_S"] for n in range(6, 43, 6)]
    for n in (48, 54, 60):
        assert predict(pure) == by_n[n]["pure_S"], "C5 anchor held-out failure"
        pure.append(by_n[n]["pure_S"])
    anchor = leading(pure[:7], 6)
    assert anchor == F(1, 36)
    return {"U": [str(u) for u in us], "J52": str(7*sum(us)),
            "pure_C5_anchor": str(anchor), "held_out_N": 60,
            "scope": "exact stated model integral; analytic moment transport unresolved"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scales", type=int, nargs="+", default=list(range(6, 61, 6)))
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    checks = gates()
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (Path(__file__), Path(__file__).with_name("spectator_reduction.py"),
                        Path(__file__).with_name("spectator_52_reduced_midpoint.py"))}
    checkpoint = json.loads(args.out.read_text()) if args.out.exists() else None
    if checkpoint is not None and checkpoint["source_sha256"] != hashes:
        raise ValueError("checkpoint source hashes differ; choose a fresh output path")
    records = checkpoint["records"] if checkpoint is not None else []
    if len({r["N"] for r in records}) != len(records):
        raise ValueError("duplicate checkpoint scales")
    if args.verify_only:
        result = certificate(records)
        if result is None or checkpoint["certificate"] != result or checkpoint["gates"] != checks:
            raise ValueError("missing or inconsistent certificate")
        print(json.dumps(result), flush=True)
        return
    for n in args.scales:
        if any(r["N"] == n for r in records):
            continue
        result = run(n)
        records.append(result)
        payload = {"source_sha256": hashes, "gates": checks, "records": records,
                   "certificate": certificate(records), "numpy_version": np.__version__,
                   "python_version": platform.python_version()}
        args.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.out.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
        tmp.replace(args.out)
        print(json.dumps(result), flush=True)
    print(json.dumps(certificate(records)), flush=True)


if __name__ == "__main__":
    main()
