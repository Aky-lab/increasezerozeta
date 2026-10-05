#!/usr/bin/env python3
"""Reduced exact network certificates via centered Ehrhart reciprocity.

Generation needs NumPy; rational interpolation and stored-record checks do not.
See notes/centered_reciprocity.md for the degree and factorization proof.
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

ROOT = Path(__file__).resolve().parents[1]
SIGNATURES = ((2,), (2,2), (2,2,2), (2,2,2,2), (4,), (5,), (6,),
              (7,), (8,), (4,2), (5,2), (6,2), (4,2,2), (4,4))
SOURCES = ("scripts/reciprocity_certificates.py", "scripts/mixed_cycle_lattice.py",
           "scripts/pure_cycle_lattice.py", "scripts/bell8_orbits.py")


def name(sig):
    return ",".join(map(str, sig))


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def evaluate(coefficients, x):
    total = F(0)
    for a in reversed(coefficients):
        total = total*x+a
    return total


def reduced_fit(b, samples):
    """S(n)=n^a(n^2-1)P(n^2), using floor(b/2) positive samples."""
    q, a = b//2, 1 if b%2==0 else 2
    if set(samples)!=set(range(2,q+2)):
        raise ValueError("incorrect interpolation nodes")
    out = [F(0)]*q
    for n, count in sorted(samples.items()):
        term = [F(count, n**a*(n*n-1))]
        for m in sorted(samples):
            if m!=n:
                term = multiply(term, [F(-m*m, n*n-m*m), F(1,n*n-m*m)])
        out = [x+y for x,y in zip(out,term)]
    # Convert to an ordinary polynomial in n, in ascending coefficient order.
    expanded = [F(0)]*(b+2)
    for j, value in enumerate(out):
        expanded[a+2*j] -= value
        expanded[a+2*j+2] += value
    return out, expanded


def certify(sig, rows):
    b, q = sum(sig), sum(sig)//2
    assert [r["n"] for r in rows]==list(range(1,q+3))
    assert rows[0]["signed_count"]==0
    samples = {r["n"]:r["signed_count"] for r in rows[1:-1]}
    reduced, expanded = reduced_fit(b,samples)
    assert evaluate(expanded,q+2)==rows[-1]["signed_count"]
    assert (expanded[-1]*math.factorial(b+1)).denominator==1
    orbit_values = []
    for i, first in enumerate(rows[0]["orbits"]):
        assert first["count"]==0
        for row in rows:
            assert row["orbits"][i]["blocks"]==first["blocks"]
            assert row["orbits"][i]["multiplicity"]==first["multiplicity"]
        rp, poly = reduced_fit(b,{r["n"]:r["orbits"][i]["count"] for r in rows[1:-1]})
        assert evaluate(poly,q+2)==rows[-1]["orbits"][i]["count"]
        assert (poly[-1]*math.factorial(b+1)).denominator==1
        orbit_values.append({"blocks":first["blocks"],"multiplicity":first["multiplicity"],
                             "value":str(poly[-1]), "polynomial":list(map(str,poly))})
    assert sum(F(o["value"])*o["multiplicity"] for o in orbit_values)==expanded[-1]
    return {"value":str(expanded[-1]), "degree_bound":b+1, "period":1,
            "parity": "odd" if b%2==0 else "even",
            "factor_power":1 if b%2==0 else 2,
            "fit_n":list(range(2,q+2)), "held_out_n":q+2,
            "reduced_polynomial":list(map(str,reduced)),
            "polynomial":list(map(str,expanded)), "orbit_values":orbit_values}


@lru_cache(None)
def compiled(b):
    from pure_cycle_lattice import terms
    return terms(b)


def generate(sig,n):
    from bell8_orbits import orbits
    from mixed_cycle_lattice import orbit_count
    from pure_cycle_lattice import sorted_cache
    started = time.monotonic()
    reps = orbits(sum(sig),sig)
    term_counts = {b:len(compiled(b)) for b in set(sig)}
    placements = sum(m for _,m in reps)
    bound = placements*math.prod(term_counts[b] for b in sig)*n**(sum(sig)+1)
    assert bound<2**63 and max(term_counts.values())*n<2**31
    caches = {b:sorted_cache(b,n,compiled(b)) for b in set(sig) if b>2}
    rows = [{"blocks":blocks, "multiplicity":mult,
             "count":orbit_count(blocks,n,caches)} for blocks,mult in reps]
    total = sum(r["multiplicity"]*r["count"] for r in rows)
    assert abs(total)<=bound
    return {"signature":list(sig), "n":n, "signed_count":total,
            "outer_free_coordinates":sum(sig)-len(sig)+1,
            "absolute_sum_bound":bound, "orbits":rows,
            "seconds":time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    hashes = {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    saved = json.loads(args.out.read_text()) if args.out.exists() else None
    if saved and saved["source_sha256"]!=hashes:
        raise ValueError("checkpoint sources differ; use a fresh output path")
    records = saved["records"] if saved else []
    assert len({(tuple(r["signature"]),r["n"]) for r in records})==len(records)
    certificates = {}
    for sig in SIGNATURES:
        for n in range(1,sum(sig)//2+3):
            if not any(tuple(r["signature"])==sig and r["n"]==n for r in records):
                record = generate(sig,n)
                records.append(record)
                print(f"{name(sig)} n={n}: {record['signed_count']}",flush=True)
            rows = sorted((r for r in records if tuple(r["signature"])==sig),key=lambda r:r["n"])
            if len(rows)==sum(sig)//2+2:
                certificates[name(sig)] = certify(sig,rows)
            import numpy as np
            payload = {"schema_version":1,"source_sha256":hashes,"records":records,
                       "certificates":certificates,"python_version":platform.python_version(),
                       "numpy_version":np.__version__,"scope":"finite network model; arithmetic transport open"}
            args.out.parent.mkdir(parents=True,exist_ok=True)
            tmp = args.out.with_suffix(".tmp")
            tmp.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
            tmp.replace(args.out)
    print(json.dumps({k:v["value"] for k,v in certificates.items()}))


if __name__=="__main__":
    main()
