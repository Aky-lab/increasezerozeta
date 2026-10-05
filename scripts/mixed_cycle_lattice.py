#!/usr/bin/env python3
"""Exact mixed-cycle integrals from bounded integral network polytopes.

Requires NumPy. Checkpoints pin all evaluator sources. The continuum-model
certificates do not establish arithmetic transport to zeta moments.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
from bell8_orbits import orbits
from pure_cycle_lattice import terms, sorted_cache, encode, leading, predict

ROOT = Path(__file__).resolve().parents[1]
SIGNATURES = ((4,2),(6,2),(4,2,2),(4,4))
SOURCES = ("scripts/mixed_cycle_lattice.py", "scripts/bell8_orbits.py",
           "scripts/pure_cycle_lattice.py")


def signature_name(sig):
    return ",".join(map(str,sig))


@lru_cache(None)
def representatives(sig):
    return orbits(sum(sig),sig)


def chart(blocks):
    """Eliminate a spanning tree of outer incidence columns over Z."""
    b,r = sum(map(len,blocks)),len(blocks)
    owner = {i:j for j,block in enumerate(blocks) for i in block}
    columns = []
    for i in range(b):
        col = [0]*r
        col[owner[i]] += 1
        col[owner[(i-1)%b]] -= 1
        columns.append(col)
    parent = list(range(r))
    def find(j):
        while parent[j]!=j:
            j = parent[j]
        return j
    tree = []
    for i in range(b):
        a,c = find(owner[i]),find(owner[(i-1)%b])
        if a!=c:
            parent[a] = c
            tree.append(i)
    assert len(tree)==r-1
    free = [i for i in range(b) if i not in tree]
    # Reduced incidence tree is unimodular. Solve all free-column RHSs.
    matrix = [[F(columns[i][j]) for i in tree]+
              [F(-columns[i][j]) for i in free] for j in range(r-1)]
    for j in range(r-1):
        pivot = next(i for i in range(j,r-1) if matrix[i][j])
        matrix[j],matrix[pivot] = matrix[pivot],matrix[j]
        divisor = matrix[j][j]
        matrix[j] = [v/divisor for v in matrix[j]]
        for i in range(r-1):
            if i!=j:
                divisor = matrix[i][j]
                matrix[i] = [v-divisor*w for v,w in zip(matrix[i],matrix[j])]
    coefficients = [[v for v in row[r-1:]] for row in matrix]
    assert all(v.denominator==1 for row in coefficients for v in row)
    coefficients = [[int(v) for v in row] for row in coefficients]
    for j in range(r):
        for k,i in enumerate(free):
            assert columns[i][j]+sum(columns[t][j]*co[k]
                                    for t,co in zip(tree,coefficients))==0
    return free,tree,coefficients


def orbit_count(blocks,n,caches,chunk_size=131072):
    b = sum(map(len,blocks))
    free,tree,coefficients = chart(blocks)
    total = 0
    points = n**len(free)
    for start in range(0,points,chunk_size):
        indices = np.arange(start,min(start+chunk_size,points),dtype=np.int64)
        x = [None]*b
        for i in free:
            x[i] = (indices % n).astype(np.int32)
            indices //= n
        valid = np.ones_like(x[free[0]],dtype=bool)
        for i,row in zip(tree,coefficients):
            value = np.zeros_like(x[free[0]])
            for a,j in zip(row,free):
                if a:
                    value += a*x[j]
            x[i] = value
            valid &= (value>=0)&(value<n)
        if not np.any(valid):
            continue
        x = [value[valid] for value in x]
        c = [x[(i+1)%b]-x[i] for i in range(b)]
        weight = np.ones_like(c[0],dtype=np.int64)
        for block in blocks:
            frequencies = [c[i] for i in block]
            assert np.all(sum(frequencies)==0)
            if len(block)==2:
                inner = np.abs(frequencies[0])
            else:
                ordered = np.sort(np.stack(frequencies),axis=0)
                keys,values = caches[len(block)]
                query = encode(list(ordered[:-1]),n-1)
                index = np.searchsorted(keys,query)
                assert np.all(index<len(keys)) and np.all(keys[index]==query)
                inner = values[index]
            weight *= inner
        total += int(np.sum(weight,dtype=np.int64))
    return total


def run(sig,n,chunk_size=131072):
    if sig not in SIGNATURES or not 1<=n<=12 or not 1<=chunk_size<=1048576:
        raise ValueError("unsupported signature, scale or chunk size")
    started = time.monotonic()
    b,r = sum(sig),len(sig)
    reps = representatives(sig)
    placements = sum(size for _,size in reps)
    term_counts = {j:len(terms(j)) for j in set(sig)}
    bound = placements*math.prod(term_counts[j] for j in sig)*n**(b+1)
    assert bound<2**63 and max(term_counts.values())*n<2**31
    caches = {j:sorted_cache(j,n,terms(j)) for j in set(sig) if j>2}
    rows = [{"blocks":blocks,"multiplicity":size,
             "count":orbit_count(blocks,n,caches,chunk_size)} for blocks,size in reps]
    total = sum(row["multiplicity"]*row["count"] for row in rows)
    assert abs(total)<=bound
    return {"signature":list(sig),"n":n,"dilation_N":n-1,"signed_count":total,
            "placements":placements,"outer_free_coordinates":b-r+1,
            "absolute_sum_bound":bound,"orbits":rows,
            "seconds":time.monotonic()-started}


def certificates(records):
    output = {}
    for sig in SIGNATURES:
        b = sum(sig)
        rows = sorted((r for r in records if tuple(r["signature"])==sig),key=lambda r:r["n"])
        if [r["n"] for r in rows]!=list(range(1,b+4)):
            continue
        values = [r["signed_count"] for r in rows]
        assert predict(values[:-1])==values[-1],"held-out aggregate differs"
        result = []
        for i,term in enumerate(rows[0]["orbits"]):
            samples = [r["orbits"][i]["count"] for r in rows]
            assert predict(samples[:-1])==samples[-1],"held-out orbit differs"
            result.append({"blocks":term["blocks"],"multiplicity":term["multiplicity"],
                           "value":str(leading(samples[:-1],b+1))})
        value = leading(values[:-1],b+1)
        assert sum(t["multiplicity"]*F(t["value"]) for t in result)==value
        output[signature_name(sig)] = {"value":str(value),"degree_bound":b+1,"period":1,
                                      "held_out_n":b+3,"orbit_values":result,
                                      "scope":"exact continuum mixed-cycle integral"}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signatures",nargs="+",default=[signature_name(s) for s in SIGNATURES])
    parser.add_argument("--n",nargs="+",type=int)
    parser.add_argument("--chunk-size",type=int,default=131072)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    hashes = {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    saved = json.loads(args.out.read_text()) if args.out.exists() else None
    if saved and saved["source_sha256"]!=hashes:
        raise ValueError("checkpoint sources differ; use a fresh output path")
    records = saved["records"] if saved else []
    assert len({(tuple(r["signature"]),r["n"]) for r in records})==len(records)
    for name in args.signatures:
        sig = tuple(map(int,name.split(",")))
        if sig not in SIGNATURES:
            raise ValueError("unsupported signature")
        for n in args.n or range(1,sum(sig)+4):
            if any(tuple(r["signature"])==sig and r["n"]==n for r in records):
                continue
            record = run(sig,n,args.chunk_size)
            records.append(record)
            payload = {"schema_version":1,"source_sha256":hashes,"records":records,
                       "certificates":certificates(records),"python_version":platform.python_version(),
                       "numpy_version":np.__version__}
            args.out.parent.mkdir(parents=True,exist_ok=True)
            tmp = args.out.with_suffix(".tmp")
            tmp.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
            tmp.replace(args.out)
            print(json.dumps({k:v for k,v in record.items() if k!="orbits"}),flush=True)
    print(json.dumps(certificates(records)),flush=True)


if __name__=="__main__":
    main()
