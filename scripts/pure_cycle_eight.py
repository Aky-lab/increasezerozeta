#!/usr/bin/env python3
"""Local, chunked integer enumeration of the pure eighth-order cycle."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
from mixed_cycle_lattice import orbit_count
from pure_cycle_lattice import terms, sorted_cache, leading, predict

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ("scripts/pure_cycle_eight.py","scripts/mixed_cycle_lattice.py",
           "scripts/pure_cycle_lattice.py","scripts/bell8_orbits.py")


def run(n,compiled,chunk_size=131072):
    if not 1<=n<=12:
        raise ValueError("require 1<=n<=12")
    started = time.monotonic()
    assert len(compiled)==94586
    bound = len(compiled)*n**9
    assert bound<2**63 and len(compiled)*n<2**31
    cache = sorted_cache(8,n,compiled)
    count = orbit_count((tuple(range(8)),),n,{8:cache},chunk_size)
    assert abs(count)<=bound
    return {"b":8,"n":n,"dilation_N":n-1,"signed_count":count,
            "terms":len(compiled),"cube_points":n**8,"sorted_frequency_tuples":len(cache[0]),
            "absolute_sum_bound":bound,"seconds":time.monotonic()-started}


def certificate(records):
    rows = sorted(records,key=lambda r:r["n"])
    if [r["n"] for r in rows]!=list(range(1,12)):
        return None
    values = [r["signed_count"] for r in rows]
    assert predict(values[:-1])==values[-1],"held-out count differs"
    return {"C":str(leading(values[:-1],9)),"degree_bound":9,"period":1,
            "held_out_n":11,"scope":"exact continuum eighth-order pure-cycle integral"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n",type=int,nargs="+")
    parser.add_argument("--chunk-size",type=int,default=131072)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    hashes = {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    saved = json.loads(args.out.read_text()) if args.out.exists() else None
    if saved and saved["source_sha256"]!=hashes:
        raise ValueError("checkpoint sources differ; use a fresh output path")
    records = saved["records"] if saved else []
    assert len({r["n"] for r in records})==len(records)
    compiled = terms(8)
    for n in args.n or range(1,12):
        if any(r["n"]==n for r in records):
            continue
        record = run(n,compiled,args.chunk_size)
        records.append(record)
        payload = {"schema_version":1,"source_sha256":hashes,"records":records,
                   "certificate":certificate(records),"python_version":platform.python_version(),
                   "numpy_version":np.__version__}
        args.out.parent.mkdir(parents=True,exist_ok=True)
        tmp = args.out.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
        tmp.replace(args.out)
        print(json.dumps(record),flush=True)
    print(json.dumps(certificate(records)),flush=True)


if __name__=="__main__":
    main()
