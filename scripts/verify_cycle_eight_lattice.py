#!/usr/bin/env python3
"""Independent eighth-order compiler and unsorted full-cube lattice sums.

Requires NumPy. No primary cumulant forms, sorted-frequency cache or
spanning-tree coordinate eliminator is used to compute these counts.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def independent_cycles(b):
    cycles = set()
    def visit(labels):
        if len(labels)<b:
            for label in range(max(labels)+2):
                visit(labels+(label,))
            return
        blocks = tuple(tuple(i for i,j in enumerate(labels) if j==label)
                       for label in range(max(labels)+1))
        # Restricted-growth label zero contains position zero.
        for tail in itertools.permutations(blocks[1:]):
            cycles.add((blocks[0],)+tail)
    visit((0,))
    return cycles


def direct_count(b,n,cycles):
    x = np.indices((n,)*b,dtype=np.int32).reshape(b,-1)
    c = [x[(i+1)%b]-x[i] for i in range(b)]
    # Direct subset sums use every frequency, including the last one.
    subsets = [np.zeros_like(c[0])]
    for mask in range(1,1<<b):
        bit = mask&-mask
        subsets.append(subsets[mask^bit]+c[bit.bit_length()-1])
    total = 0
    for blocks in cycles:
        mask = 0
        low,high = subsets[0].copy(),subsets[0].copy()
        for block in blocks[:-1]:
            mask |= sum(1<<i for i in block)
            low = np.minimum(low,subsets[mask])
            high = np.maximum(high,subsets[mask])
        total += (-1)**(len(blocks)-1)*int(np.sum(np.maximum(n-high+low,0),dtype=np.int64))
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    started = time.monotonic()
    cycles = independent_cycles(8)
    assert len(cycles)==94586
    # Compare compiler sets after the independent compiler is constructed.
    from pure_cycle_lattice import terms
    primary = {tuple(tuple(sorted(block)) for block in blocks) for _,_,blocks in terms(8)}
    assert cycles==primary
    saved_path = ROOT/"results/pure_cycle_eight_2026-10-05.json"
    saved = json.loads(saved_path.read_text())
    expected = {r["n"]:r["signed_count"] for r in saved["records"]}
    checks = []
    for n in (2,3):
        observed = direct_count(8,n,cycles)
        assert observed==expected[n],(n,observed,expected[n])
        checks.append({"n":n,"signed_count":observed,"cube_points":n**8})
        print(f"C8: independent unsorted full cube n={n} agrees: {observed}",flush=True)
    paths = (Path(__file__),ROOT/"scripts/pure_cycle_lattice.py",saved_path)
    report = {"all_checks_passed":True,"independent_cyclic_terms":len(cycles),
              "full_cube_checks":checks,"source_sha256":{
                  p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in paths},"seconds":time.monotonic()-started}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("EIGHTH-ORDER COMPILER AND DIRECT LATTICE CHECKS PASSED",flush=True)


if __name__=="__main__":
    main()
