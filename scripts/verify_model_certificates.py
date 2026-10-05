#!/usr/bin/env python3
"""Verify recorded integer certificates and provenance without NumPy.

This checks transcript arithmetic and source identity. Recomputing the
lattice sums and independent term checks uses the NumPy research tools.
"""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def leading(values, degree, step=1):
    require(len(values)==degree+1,"wrong interpolation sample count")
    work = list(values)
    for _ in range(degree):
        work = [b-a for a,b in zip(work,work[1:])]
    return F(work[0],math.factorial(degree)*step**degree)


def predict(values):
    work, result = list(values),0
    while work:
        result += work[-1]
        work = [b-a for a,b in zip(work,work[1:])]
    return result


def hash_matches(path, expected):
    resolved = (ROOT/path).resolve()
    require(resolved.is_relative_to(ROOT),"source path outside repository")
    require(hashlib.sha256(resolved.read_bytes()).hexdigest()==expected,
            f"source hash differs: {path}")


def pure_cycle_value(records,b):
    rows = sorted((r for r in records if r["b"]==b),key=lambda r:r["n"])
    require([r["n"] for r in rows]==list(range(1,b+4)),"missing pure-cycle scales")
    values = [r["signed_count"] for r in rows[:-1]]
    require(predict(values)==rows[-1]["signed_count"],"pure-cycle held-out count differs")
    return leading(values,b+1)


def main():
    joint = json.loads((ROOT/"results/spectator_52_exact_lattice_2026-10-05.json").read_text())
    for name, digest in joint["source_sha256"].items():
        hash_matches("scripts/"+name,digest)
    rows = sorted(joint["records"],key=lambda r:r["N"])
    require([r["N"] for r in rows]==list(range(6,61,6)),"missing joint scales")
    components = []
    for d in range(3):
        values = [r["S"][d] for r in rows[:-1]]
        require(predict(values)==rows[-1]["S"][d],"joint held-out count differs")
        components.append(leading(values,8,6)/6)
    require([str(u) for u in components]==joint["certificate"]["U"],"joint components differ")
    require(7*sum(components)==F(joint["certificate"]["J52"]),"joint total differs")
    pure = [r["pure_S"] for r in rows[:7]]
    for r in rows[7:]:
        require(predict(pure)==r["pure_S"],"joint pure-C5 held-out count differs")
        pure.append(r["pure_S"])
    require(leading(pure[:7],6,6)==F(1,36),"joint pure-C5 anchor differs")

    cycles = json.loads((ROOT/"results/pure_cycle_flow_2026-10-05.json").read_text())
    hash_matches("scripts/pure_cycle_lattice.py",cycles["source_sha256"])
    values = {b:pure_cycle_value(cycles["records"],b) for b in (4,5,6,7)}
    require(values=={4:F(-1,60),5:F(1,36),6:F(-1,126),7:F(-17,360)},"cycle values differ")
    for b,v in values.items():
        require(str(v)==cycles["certificates"][str(b)]["C"],"cycle certificate differs")

    independent = json.loads((ROOT/"results/pure_cycle_independent_checks_2026-10-05.json").read_text())
    for path,digest in independent["source_sha256"].items():
        hash_matches(path,digest)
    require(independent["all_checks_passed"],"independent checks failed")
    fixture = json.loads((ROOT/"results/pure_cycle_reference_anchors.json").read_text())
    anchors = {(a["b"],a["term"]):F(a["value"]) for a in fixture["anchors"]}
    checked = set()
    for row in independent["per_term_checks"]:
        key = row["b"],row["term"]
        samples = row["lattice_count_samples"]
        require(predict(samples[:-1])==samples[-1],"per-term held-out count differs")
        require(leading(samples[:-1],row["b"]+1)==anchors[key],"per-term reference differs")
        checked.add(key)
    require(checked==set(anchors),"missing independent anchor terms")
    seventh = sorted((r for r in cycles["records"] if r["b"]==7),key=lambda r:r["n"])
    next_value = predict([r["signed_count"] for r in seventh[:9]])
    next_two = predict([r["signed_count"] for r in seventh[:9]]+[next_value])
    require(independent["additional_held_out"]["n"]==11 and
            independent["additional_held_out"]["signed_count"]==next_two,
            "additional held-out count differs")
    m7 = F(1717,90)+7*sum(components)+values[7]
    require(m7==F(862,45),"seventh-order model ledger differs")
    print("RECORDED INTEGER CERTIFICATES AND SOURCE HASHES VERIFIED")
    print("J52 =",7*sum(components),"; C7 =",values[7],"; model m7 =",m7)
    print("Scope: transcript arithmetic and provenance; analytic transport remains open.")


if __name__=="__main__":
    main()
