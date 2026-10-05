#!/usr/bin/env python3
"""Standard-library verification using independent scalar frequency sums.

The placement and cyclic-partition enumerators, full-frequency cumulants
and prefix walks do not import the primary evaluator's implementations.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {"scripts/mixed_cycle_lattice.py","scripts/bell8_orbits.py",
           "scripts/pure_cycle_lattice.py"}
EXPECTED = {"4,2":(6,15,3),"6,2":(8,28,4),"4,2,2":(8,210,22),"4,4":(8,35,7)}


def require(condition,message):
    if not condition:
        raise ValueError(message)


def leading(values,degree):
    require(len(values)==degree+1,"wrong number of samples")
    work = list(values)
    for _ in range(degree):
        work = [b-a for a,b in zip(work,work[1:])]
    return F(work[0],math.factorial(degree))


def predict(values):
    total,work = 0,list(values)
    while work:
        total += work[-1]
        work = [b-a for a,b in zip(work,work[1:])]
    return total


@lru_cache(None)
def placements(sig):
    """Restricted-growth strings, rather than recursive block insertion."""
    b = sum(sig)
    output = []
    def visit(labels):
        if len(labels)==b:
            blocks = tuple(tuple(i for i,j in enumerate(labels) if j==label)
                           for label in range(max(labels)+1))
            if tuple(sorted(map(len,blocks),reverse=True))==sig:
                output.append(tuple(sorted(blocks)))
            return
        for label in range(max(labels)+2):
            visit(labels+(label,))
    visit((0,))
    return tuple(sorted(output))


def orbit_inventory(sig):
    b = sum(sig)
    unseen = set(placements(sig))
    output = []
    while unseen:
        rep = min(unseen)
        images = {tuple(sorted(tuple(sorted((sign*i+shift)%b for i in block))
                               for block in rep)) for sign in (-1,1) for shift in range(b)}
        require(images<=set(placements(sig)),"invalid dihedral image")
        output.append((rep,len(images)))
        unseen -= images
    return output


@lru_cache(None)
def independent_terms(b):
    """Permutation-and-cut construction with cyclic rotation deduplication."""
    cycles = set()
    for perm in itertools.permutations(range(b)):
        for cuts in itertools.product((False,True),repeat=b-1):
            ends = [0]+[i+1 for i,cut in enumerate(cuts) if cut]+[b]
            blocks = [tuple(sorted(perm[a:z])) for a,z in zip(ends,ends[1:])]
            anchor = next(i for i,block in enumerate(blocks) if 0 in block)
            cycles.add(tuple(blocks[anchor:]+blocks[:anchor]))
    require(len(cycles)=={4:26,6:1082}[b],"cyclic term inventory differs")
    output = []
    for blocks in sorted(cycles):
        mask,prefixes = 0,[]
        for block in blocks[:-1]:
            mask |= sum(1<<i for i in block)
            prefixes.append(mask)
        output.append(((-1)**(len(blocks)-1),tuple(prefixes)))
    return tuple(output)


def scalar_cumulant(frequencies,n):
    b = len(frequencies)
    require(sum(frequencies)==0,"frequency block is not closed")
    if b==2:
        return abs(frequencies[0])
    subsets = [0]*(1<<b)
    for mask in range(1,1<<b):
        bit = mask&-mask
        subsets[mask] = subsets[mask^bit]+frequencies[bit.bit_length()-1]
    total = 0
    for sign,masks in independent_terms(b):
        positions = [0]+[subsets[m] for m in masks]
        total += sign*max(n-max(positions)+min(positions),0)
    return total


@lru_cache(None)
def frequency_options(b,n):
    options = []
    for free in itertools.product(range(1-n,n),repeat=b-1):
        last = -sum(free)
        if abs(last)>=n:
            continue
        frequencies = free+(last,)
        value = scalar_cumulant(frequencies,n)
        if value:
            options.append((frequencies,value))
    return tuple(options)


def direct_count(sig,n):
    total = 0
    for blocks in placements(sig):
        options = [frequency_options(len(block),n) for block in blocks]
        for chosen in itertools.product(*options):
            increments = [0]*sum(sig)
            for block,(frequencies,_) in zip(blocks,chosen):
                for label,value in zip(block,frequencies):
                    increments[label] = value
            walk = [0]+list(itertools.accumulate(increments))
            outer = max(n-max(walk)+min(walk),0)
            total += outer*math.prod(value for _,value in chosen)
    return total


def main():
    data = json.loads((ROOT/"results/mixed_cycle_flow_2026-10-05.json").read_text())
    require(set(data["source_sha256"])==SOURCES,"source inventory differs")
    for name,expected in data["source_sha256"].items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,"source hash differs: "+name)
    require(set(data["certificates"])==set(EXPECTED),"certificate inventory differs")
    for name,(b,num,norb) in EXPECTED.items():
        sig = tuple(map(int,name.split(",")))
        rows = sorted((r for r in data["records"] if tuple(r["signature"])==sig),key=lambda r:r["n"])
        require([r["n"] for r in rows]==list(range(1,b+4)),"scale inventory differs")
        cert = data["certificates"][name]
        require((cert["degree_bound"],cert["period"],cert["held_out_n"])==(b+1,1,b+3),
                "certificate parameters differ")
        inventory = orbit_inventory(sig)
        require(len(placements(sig))==num and len(inventory)==norb,"placement inventory differs")
        observed = [(tuple(map(tuple,t["blocks"])),t["multiplicity"]) for t in cert["orbit_values"]]
        require(inventory==observed,"independent orbit inventory differs")
        values = [r["signed_count"] for r in rows]
        require(predict(values[:-1])==values[-1],"held-out aggregate differs")
        value = leading(values[:-1],b+1)
        require(str(value)==cert["value"],"aggregate value differs")
        require(sum(t["multiplicity"]*F(t["value"]) for t in cert["orbit_values"])==value,
                "weighted orbit values disagree")
        for row in rows:
            require(row["placements"]==num and row["dilation_N"]==row["n"]-1,
                    "record parameters differ")
            identity = [(tuple(map(tuple,t["blocks"])),t["multiplicity"]) for t in row["orbits"]]
            require(identity==inventory,"record orbit inventory differs")
            require(sum(t["multiplicity"]*t["count"] for t in row["orbits"])==row["signed_count"],
                    "weighted orbit counts disagree")
        for i,term in enumerate(cert["orbit_values"]):
            samples = [r["orbits"][i]["count"] for r in rows]
            require(predict(samples[:-1])==samples[-1],"held-out orbit differs")
            require(str(leading(samples[:-1],b+1))==term["value"],"orbit value differs")
        for n in (2,3):
            require(direct_count(sig,n)==values[n-1],"independent scalar frequency sum differs")
        print(f"{{{name}}} = {value}; exact record and independent scalar checks passed",flush=True)
    from model_moments import MODEL_J42, MODEL_J62, MODEL_J422, MODEL_J44
    shared = dict(zip(EXPECTED,(MODEL_J42,MODEL_J62,MODEL_J422,MODEL_J44)))
    require(all(F(data["certificates"][name]["value"])==value for name,value in shared.items()),
            "shared model class inputs differ")
    print("ALL MIXED-CYCLE CERTIFICATE CHECKS PASSED")
    print("Scope: continuum model; analytic transport remains open.")


if __name__=="__main__":
    main()
