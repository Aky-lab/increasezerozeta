#!/usr/bin/env python3
"""Verify pair-cycle records, independent finite counts and exact cell integrals.

Uses only the standard library. The small direct counts recompute the
integrand without dihedral reduction; full lattice recomputation uses NumPy.
"""
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition,message):
    if not condition:
        raise ValueError(message)


def leading(values,degree):
    require(len(values)==degree+1,"wrong sample count")
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


def matchings(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for i,second in enumerate(labels[1:],1):
        for remaining in matchings(labels[1:i]+labels[i+1:]):
            yield ((first,second),)+remaining


def direct_count(k,n):
    total = 0
    for pairs in matchings(tuple(range(2*k))):
        for values in itertools.product(range(1-n,n),repeat=k):
            weight = math.prod(abs(v) for v in values)
            if not weight:
                continue
            increments = [0]*(2*k)
            for (first,second),v in zip(pairs,values):
                increments[first],increments[second] = v,-v
            positions = [0]+list(itertools.accumulate(increments))
            total += max(n-max(positions)+min(positions),0)*weight
    return total


def simplex_integral(nested):
    # Exact integral of ov(0,v,a,b)*|v*a*(b-a)| for nested=True,
    # or ov(0,v,a,b)*|v*a*b| for nested=False, on 24 sorted cells.
    total = F(0)
    for order in itertools.permutations(range(4)):
        rank = {label:i for i,label in enumerate(order)}
        def coordinate(label):
            return tuple(int(j<rank[label])-int(j<rank[0]) for j in range(3))
        v,a,b = (coordinate(i) for i in (1,2,3))
        forms = [v,a,tuple(y-x for x,y in zip(a,b))] if nested else [v,a,b]
        signs = [1 if rank[1]>rank[0] else -1,1 if rank[2]>rank[0] else -1,
                 (1 if rank[3]>rank[2] else -1) if nested else (1 if rank[3]>rank[0] else -1)]
        polynomial = {(0,0,0):1}
        for form,sign in zip(forms,signs):
            updated = {}
            for exp,coefficient in polynomial.items():
                for j,value in enumerate(form):
                    if value:
                        power = list(exp); power[j]+=1; power=tuple(power)
                        updated[power] = updated.get(power,0)+coefficient*sign*value
            polynomial = updated
        total += sum(F(coefficient*math.prod(math.factorial(e) for e in exp),
                       math.factorial(sum(exp)+4)) for exp,coefficient in polynomial.items())
    return total


def main():
    data = json.loads((ROOT/"results/paired_cycle_flow_2026-10-05.json").read_text())
    require(set(data["source_sha256"])=={"scripts/paired_cycle_lattice.py","scripts/bell8_orbits.py"},
            "source inventory differs")
    for name,expected in data["source_sha256"].items():
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT),"source path outside repository")
        require(hashlib.sha256(path.read_bytes()).hexdigest()==expected,"source hash differs: "+name)
    totals = {}
    for k in (1,2,3,4):
        rows = sorted((r for r in data["records"] if r["k"]==k),key=lambda r:r["n"])
        require([r["n"] for r in rows]==list(range(1,2*k+4)),"missing scales")
        values = [r["signed_count"] for r in rows]
        require(predict(values[:-1])==values[-1],"held-out aggregate differs")
        totals[k] = leading(values[:-1],2*k+1)
        cert = data["certificates"][str(k)]
        require(cert["degree_bound"]==2*k+1 and cert["period"]==1 and cert["held_out_n"]==2*k+3,
                "certificate parameters differ")
        require(str(totals[k])==cert["value"],"aggregate value differs")
        expected_placements = math.factorial(2*k)//(2**k*math.factorial(k))
        terms = cert["orbit_values"]
        require(len(terms)=={1:1,2:2,3:5,4:17}[k],"orbit inventory differs")
        require(sum(t["multiplicity"] for t in terms)==expected_placements,"placement count differs")
        require(sum(t["multiplicity"]*F(t["value"]) for t in terms)==totals[k],
                "weighted orbit values disagree")
        for row in rows:
            require(row["placements"]==expected_placements and row["dilation_N"]==row["n"]-1,
                    "record parameters differ")
            require(len(row["orbits"])==len(terms),"record orbit inventory differs")
            require(sum(t["multiplicity"]*t["count"] for t in row["orbits"])==row["signed_count"],
                    "weighted orbit counts disagree")
            for sample,term in zip(row["orbits"],terms):
                require(sample["pairs"]==term["pairs"] and sample["multiplicity"]==term["multiplicity"],
                        "orbit identity differs")
        for i,term in enumerate(cert["orbit_values"]):
            samples = [row["orbits"][i]["count"] for row in rows]
            require(predict(samples[:-1])==samples[-1],"held-out orbit differs")
            require(str(leading(samples[:-1],2*k+1))==term["value"],"orbit value differs")
        for n in (2,3,4):
            require(direct_count(k,n)==values[n-1],"independent full-placement count differs")
    require(totals=={1:F(1,3),2:F(4,15),3:F(32,105),4:F(1661,3780)},"pair-model values differ")
    adjacent,nested = simplex_integral(False),simplex_integral(True)
    require((adjacent,nested)==(F(3,70),F(17,420)),"independent simplex integrals differ")
    orbit_values = {tuple(map(tuple,row["pairs"])):F(row["value"])
                    for row in data["certificates"]["3"]["orbit_values"]}
    require(orbit_values[((0,1),(2,3),(4,5))]==adjacent,"adjacent lattice/cell disagreement")
    require(orbit_values[((0,1),(2,5),(3,4))]==nested,"nested lattice/cell disagreement")
    require(adjacent!=nested,"incorrect noncrossing equivalence")
    from model_moments import T222, PAIR4, MODEL_MOMENTS, M7_BASE
    require(T222==totals[3] and PAIR4==totals[4],"shared model inputs disagree")
    require(MODEL_MOMENTS[6]==F(640,63) and M7_BASE==F(685,36),"dependent model inputs disagree")
    print("PAIR-MODEL CERTIFICATES, DIRECT COUNTS AND EXACT SIMPLEX CHECKS PASSED")
    print("{2,2,2} = 32/105; {2^4} = 1661/3780")
    print("Scope: finite continuum model; analytic transport remains open.")


if __name__=="__main__":
    main()
