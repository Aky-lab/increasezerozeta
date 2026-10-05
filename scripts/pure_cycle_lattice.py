#!/usr/bin/env python3
"""Exact pure-cycle integrals by integral flow polytopes and Ehrhart counts.

See notes/pure_cycle_flow_polytopes.md. NumPy is required. The output
certifies finite continuum integrals; arithmetic transport is separate.
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


def partitions(labels):
    if not labels:
        yield []
        return
    first, rest = labels[0], labels[1:]
    for part in partitions(rest):
        for j in range(len(part)):
            yield part[:j]+[[first]+part[j]]+part[j+1:]
        yield [[first]]+part


def terms(b):
    result = []
    for part in partitions(list(range(b))):
        anchor = next(j for j, block in enumerate(part) if 0 in block)
        others = [j for j in range(len(part)) if j != anchor]
        for order in itertools.permutations(others):
            blocks = tuple(tuple(part[j]) for j in (anchor,)+order)
            prefix, forms = set(), []
            for block in blocks[:-1]:
                prefix.update(block)
                forms.append(tuple(int(j in prefix)-int(b-1 in prefix)
                                   for j in range(b-1)))
            result.append(((-1)**(len(part)-1), tuple(forms), blocks))
    return result


def cumulant(xs, n, compiled):
    zero = np.zeros_like(xs[0], dtype=np.int32)
    values = {}
    for _, forms, _ in compiled:
        for form in forms:
            if form not in values:
                value = zero.copy()
                for a, x in zip(form, xs):
                    if a:
                        value += a*x
                values[form] = value
    total = zero.copy()
    for sign, forms, _ in compiled:
        positions = [zero]+[values[f] for f in forms]
        span = np.maximum.reduce(positions)-np.minimum.reduce(positions)
        total += sign*np.maximum(n-span, 0)
    return total


def encode(xs, radius):
    value = np.zeros_like(xs[0], dtype=np.int64)
    for x in xs:
        value = value*(2*radius+1)+x+radius
    return value


def sorted_cache(b, n, compiled):
    radius = n-1
    rows = []
    for prefix in itertools.combinations_with_replacement(range(-radius, radius+1), b-1):
        last = -sum(prefix)
        if prefix[-1] <= last <= radius:
            rows.append(prefix)
    xs = list(np.array(rows, dtype=np.int32).T)
    keys = encode(xs, radius)
    assert np.all(keys[1:] > keys[:-1])
    values = cumulant(xs, n, compiled)
    assert np.all(np.abs(values) <= len(compiled)*n)
    return keys, values


def run(b, n, cached=True):
    if not 2 <= b <= 7 or not 1 <= n <= 12:
        raise ValueError("require 2<=b<=7 and 1<=n<=12")
    started = time.monotonic()
    compiled = terms(b)
    # Each point contributes at most term_count*n in absolute value.
    bound = len(compiled)*n**(b+1)
    assert bound < 2**63 and len(compiled)*n < 2**31
    keys, values = sorted_cache(b, n, compiled) if cached else (None, None)
    grid = np.arange(n, dtype=np.int32)
    positions = np.meshgrid(*([grid]*(b-1)), indexing="ij")
    total = 0
    for first in range(n):
        xs = [positions[0]-first]
        xs += [positions[j]-positions[j-1] for j in range(1, b-1)]
        if cached:
            last = first-positions[-1]
            ordered = np.sort(np.stack(xs+[last]), axis=0)
            query = encode(list(ordered[:-1]), n-1)
            index = np.searchsorted(keys, query)
            assert np.all(index < len(keys))
            assert np.all(keys[index] == query)
            contribution = values[index]
        else:
            contribution = cumulant(xs, n, compiled)
        total += int(np.sum(contribution, dtype=np.int64))
    return {"b": b, "n": n, "dilation_N": n-1, "signed_count": total,
            "terms": len(compiled), "cube_points": n**b,
            "sorted_frequency_tuples": len(keys) if cached else None,
            "absolute_sum_bound": bound, "seconds": time.monotonic()-started}


def leading(values, degree):
    work = list(values)
    assert len(work) == degree+1
    for _ in range(degree):
        work = [y-x for x, y in zip(work, work[1:])]
    return F(work[0], math.factorial(degree))


def predict(values):
    work, result = list(values), 0
    while work:
        result += work[-1]
        work = [y-x for x, y in zip(work, work[1:])]
    return result


def certificates(records):
    result = {}
    for b in sorted({r["b"] for r in records}):
        by_n = {r["n"]: r["signed_count"] for r in records if r["b"] == b}
        if not all(n in by_n for n in range(1, b+4)):
            continue
        samples = [by_n[n] for n in range(1, b+3)]
        assert predict(samples) == by_n[b+3], "held-out count differs"
        result[str(b)] = {"C": str(leading(samples, b+1)), "degree_bound": b+1,
                          "period": 1, "held_out_n": b+3,
                          "scope": "exact continuum model integral"}
    return result


def structural_gates():
    checked = 0
    for b, expected in ((2, 2), (3, 6), (4, 26), (5, 150), (6, 1082), (7, 9366)):
        compiled = terms(b)
        assert len(compiled) == expected
        for _, forms, blocks in compiled:
            assert sorted(itertools.chain.from_iterable(blocks)) == list(range(b))
            m = len(blocks)
            owner = {label:j for j, block in enumerate(blocks) for label in block}
            columns = []
            for i in range(b):
                col = [0]*m
                col[owner[i]] += 1
                col[owner[(i-1)%b]] -= 1
                columns.append(col)
            for j in range(m):
                col = [0]*m
                col[(j-1)%m] += 1
                col[j] -= 1
                columns.append(col)
            for col in columns:
                assert all(v in (-1, 0, 1) for v in col)
                assert sum(col) == 0 and sum(abs(v) for v in col) in (0, 2)
            # Check the flow equations after eliminating the inner cycle's
            # spanning path, using a nonsymmetric integer test vector.
            x = [(i*i+3*i+2)%11 for i in range(b)]
            c = [x[(i+1)%b]-x[i] for i in range(b)]
            y = [3]+[3+sum(a*v for a,v in zip(f,c[:-1])) for f in forms]
            z = x+y
            assert all(sum(col[j]*v for col,v in zip(columns,z)) == 0 for j in range(m))
            checked += 1
    comparisons = []
    for b in (2, 3, 4, 5, 6, 7):
        for n in (1, 2):
            direct, cached = run(b, n, False), run(b, n, True)
            assert direct["signed_count"] == cached["signed_count"]
            comparisons.append([b, n])
    # Independent permutation-and-cut construction of cyclic partitions.
    for b in range(2, 6):
        cycles = set()
        for perm in itertools.permutations(range(b)):
            for mask in range(1 << (b-1)):
                cuts = [0]+[j+1 for j in range(b-1) if mask & (1<<j)]+[b]
                blocks = [tuple(sorted(perm[lo:hi])) for lo,hi in zip(cuts,cuts[1:])]
                anchor = next(j for j,block in enumerate(blocks) if 0 in block)
                cycles.add(tuple(blocks[anchor:]+blocks[:anchor]))
        recursive = {tuple(tuple(sorted(block)) for block in blocks)
                     for _,_,blocks in terms(b)}
        assert cycles == recursive
    return {"incidence_and_elimination_terms_checked": checked,
            "cached_vs_direct": comparisons, "independent_compiler_orders": [2,3,4,5]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b", type=int, nargs="+", default=[4, 5, 6, 7])
    parser.add_argument("--n", type=int, nargs="+")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    gates = structural_gates()
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    saved = json.loads(args.out.read_text()) if args.out.exists() else None
    if saved is not None and saved["source_sha256"] != source_hash:
        raise ValueError("checkpoint source differs; use a fresh output path")
    records = saved["records"] if saved else []
    if len({(r["b"],r["n"]) for r in records}) != len(records):
        raise ValueError("duplicate scales")
    if args.verify_only:
        derived = certificates(records)
        if not derived or derived != saved["certificates"] or gates != saved["gates"]:
            raise ValueError("incomplete or inconsistent certificate")
        print(json.dumps(derived), flush=True)
        return
    for b in args.b:
        for n in args.n or range(1, b+4):
            if any(r["b"] == b and r["n"] == n for r in records):
                continue
            record = run(b, n)
            records.append(record)
            payload = {"source_sha256": source_hash, "gates": gates, "records": records,
                       "certificates": certificates(records), "python_version": platform.python_version(),
                       "numpy_version": np.__version__}
            args.out.parent.mkdir(parents=True, exist_ok=True)
            tmp = args.out.with_suffix(".tmp")
            tmp.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
            tmp.replace(args.out)
            print(json.dumps(record), flush=True)
    print(json.dumps(certificates(records)), flush=True)


if __name__ == "__main__":
    main()
