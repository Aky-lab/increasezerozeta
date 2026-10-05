#!/usr/bin/env python3
"""Independent direct-block and bounded-flow checks of pure-cycle counts."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

import pure_cycle_lattice as cycle


def direct_term_counts(b, n, selected):
    """Positive per-term Ehrhart counts, without frequency sorting or forms."""
    compiled = cycle.terms(b)
    grid = np.arange(n, dtype=np.int32)
    x = np.meshgrid(*([grid]*(b-1)), indexing="ij")
    totals = {i:0 for i in selected}
    zero = np.zeros_like(x[0])
    for first in range(n):
        c = [x[0]-first]+[x[j]-x[j-1] for j in range(1,b-1)]+[first-x[-1]]
        for i in selected:
            _, _, blocks = compiled[i]
            q, lo, hi = zero.copy(), zero.copy(), zero.copy()
            for block in blocks[:-1]:
                for label in block:
                    q = q+c[label]
                lo, hi = np.minimum(lo,q), np.maximum(hi,q)
            count = np.maximum(n-hi+lo,0)
            totals[i] += int(np.sum(count,dtype=np.int64))
    return totals


def brute_flow_count(b, n, blocks):
    """Enumerate every bounded x and y, checking conservation directly."""
    total = 0
    for x in itertools.product(range(n),repeat=b):
        c = [x[(i+1)%b]-x[i] for i in range(b)]
        block_sums = [sum(c[i] for i in block) for block in blocks]
        for y in itertools.product(range(n),repeat=len(blocks)):
            total += all(y[(j+1)%len(blocks)]-y[j] == block_sums[j]
                         for j in range(len(blocks)))
    return total


def compiler_gate(b):
    independent = set()
    for perm in itertools.permutations(range(b)):
        for mask in range(1<<(b-1)):
            cuts = [0]+[j+1 for j in range(b-1) if mask & (1<<j)]+[b]
            blocks = [tuple(sorted(perm[lo:hi])) for lo,hi in zip(cuts,cuts[1:])]
            first = next(j for j,block in enumerate(blocks) if 0 in block)
            independent.add(tuple(blocks[first:]+blocks[:first]))
    recursive = {tuple(tuple(sorted(block)) for block in blocks)
                 for _,_,blocks in cycle.terms(b)}
    assert independent == recursive
    return len(independent)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    fixture_path = root/"results/pure_cycle_reference_anchors.json"
    fixture = json.loads(fixture_path.read_text())
    started = time.monotonic()
    checks = []
    for b in (5,6):
        anchors = [a for a in fixture["anchors"] if a["b"]==b]
        ids = [a["term"] for a in anchors]
        samples = {i:[] for i in ids}
        for n in range(1,b+4):
            counts = direct_term_counts(b,n,ids)
            for i in ids:
                samples[i].append(counts[i])
        compiled = cycle.terms(b)
        for anchor in anchors:
            i = anchor["term"]
            assert compiled[i][0] == anchor["sign"]
            assert cycle.predict(samples[i][:-1]) == samples[i][-1]
            value = cycle.leading(samples[i][:-1],b+1)
            assert value == F(anchor["value"]), (b,i,value,anchor["value"])
            checks.append({"b":b,"term":i,"value":str(value),
                           "lattice_count_samples":samples[i],"held_out_n":b+3})
        print(f"b={b}: {len(anchors)} per-term reference fractions and held-out counts agree",flush=True)
    compiler_counts = {str(b):compiler_gate(b) for b in (6,7)}
    print("Independent permutation-and-cut compilers agree at orders 6 and 7",flush=True)
    b, selected = 7, [0,1,4,26,500,9365]
    direct = direct_term_counts(b,2,selected)
    for i in selected:
        assert brute_flow_count(b,2,cycle.terms(b)[i][2]) == direct[i]
    print("Six bounded-flow enumerations agree with overlap counts",flush=True)
    uncached, cached = cycle.run(7,4,False), cycle.run(7,4,True)
    assert uncached["signed_count"] == cached["signed_count"]
    saved = json.loads((root/"results/pure_cycle_flow_2026-10-05.json").read_text())
    base = sorted((r for r in saved["records"] if r["b"]==7),key=lambda r:r["n"])
    assert [r["n"] for r in base] == list(range(1,11))
    values = [r["signed_count"] for r in base[:9]]
    for n in (10,11):
        prediction = cycle.predict(values)
        observed = base[-1] if n==10 else cycle.run(7,n)
        assert prediction == observed["signed_count"]
        values.append(prediction)
    print("Uncached n=4 sum and additional n=11 polynomial check agree",flush=True)
    sources = [Path(__file__),Path(cycle.__file__),fixture_path,
               root/"results/pure_cycle_flow_2026-10-05.json"]
    report = {"all_checks_passed":True,"source_sha256":{
                  str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sources},"per_term_checks":checks,
              "independent_compiler_counts":compiler_counts,
              "bounded_flow_term_indices":selected,
              "uncached_n4_signed_count":uncached["signed_count"],
              "additional_held_out":observed,"seconds":time.monotonic()-started}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding="utf-8")
    print(f"Verification record: {args.out}",flush=True)


if __name__=="__main__":
    main()
