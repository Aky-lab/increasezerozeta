#!/usr/bin/env python3
"""Integer-count certificates for cycles partitioned entirely into pairs.

See notes/paired_cycle_flow_polytopes.md. Requires NumPy; the result is a
continuum-model integral, with arithmetic transport a separate obligation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from bell8_orbits import orbits

ROOT = Path(__file__).resolve().parents[1]


def leading(values, degree):
    if len(values) != degree+1:
        raise ValueError("wrong interpolation sample count")
    work = list(values)
    for _ in range(degree):
        work = [b-a for a,b in zip(work,work[1:])]
    return F(work[0],math.factorial(degree))


def predict(values):
    work,total = list(values),0
    while work:
        total += work[-1]
        work = [b-a for a,b in zip(work,work[1:])]
    return total


def run(k,n):
    if not 1 <= k <= 4 or not 1 <= n <= 12:
        raise ValueError("require 1<=k<=4 and 1<=n<=12")
    started = time.monotonic()
    representatives = orbits(2*k,(2,)*k)
    placements = sum(size for _,size in representatives)
    assert placements == math.factorial(2*k)//(2**k*math.factorial(k))
    bound = placements*(2*n-1)**k*n*(n-1)**k
    assert bound < 2**63
    grid = np.arange(1-n,n,dtype=np.int32)
    frequencies = np.meshgrid(*([grid]*k),indexing="ij")
    weight = np.ones_like(frequencies[0],dtype=np.int64)
    for v in frequencies:
        weight *= np.abs(v)
    rows = []
    for pairs,size in representatives:
        owner = {label:(j,1 if label==pair[0] else -1)
                 for j,pair in enumerate(pairs) for label in pair}
        current = np.zeros_like(frequencies[0])
        low,high = current.copy(),current.copy()
        for label in range(2*k):
            j,sign = owner[label]
            current = current+sign*frequencies[j]
            low,high = np.minimum(low,current),np.maximum(high,current)
        assert np.all(current==0)
        count = int(np.sum(np.maximum(n-high+low,0)*weight,dtype=np.int64))
        rows.append({"pairs":pairs,"multiplicity":size,"count":count})
    total = sum(r["multiplicity"]*r["count"] for r in rows)
    assert 0 <= total <= bound
    return {"k":k,"n":n,"dilation_N":n-1,"signed_count":total,
            "placements":placements,"orbits":rows,"absolute_sum_bound":bound,
            "seconds":time.monotonic()-started}


def certificates(records):
    out = {}
    for k in sorted({r["k"] for r in records}):
        rows = sorted((r for r in records if r["k"]==k),key=lambda r:r["n"])
        assert [r["n"] for r in rows]==list(range(1,2*k+4))
        values = [r["signed_count"] for r in rows]
        assert predict(values[:-1])==values[-1]
        terms = []
        for i,orbit in enumerate(rows[0]["orbits"]):
            samples = [row["orbits"][i]["count"] for row in rows]
            assert predict(samples[:-1])==samples[-1]
            terms.append({"pairs":orbit["pairs"],"multiplicity":orbit["multiplicity"],
                          "value":str(leading(samples[:-1],2*k+1))})
        out[str(k)] = {"value":str(leading(values[:-1],2*k+1)),
                       "degree_bound":2*k+1,"period":1,"held_out_n":2*k+3,
                       "orbit_values":terms,"scope":"exact continuum pair-model integral"}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    records = []
    for k in (1,2,3,4):
        for n in range(1,2*k+4):
            records.append(run(k,n))
    result = {"schema_version":1,"records":records,"certificates":certificates(records),
              "source_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                               for p in ("scripts/paired_cycle_lattice.py","scripts/bell8_orbits.py")}}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for k,cert in result["certificates"].items():
        print(f"{{2^{k}}} = {cert['value']}; per-orbit and aggregate held-out counts agree")


if __name__=="__main__":
    main()
