#!/usr/bin/env python3
"""Verify reduced certificates against every earlier network count."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from model_moments import (PAIR1, PAIR2, T222, PAIR4, C4, C5, C6,
                          MODEL_C7, MODEL_C8, MODEL_J42, MODEL_J52,
                          MODEL_J62, MODEL_J422, MODEL_J44)
from reciprocity_certificates import ROOT, SOURCES, SIGNATURES, name, certify, evaluate

RECORD = ROOT/"results/reciprocity_2026-10-05.json"


def canonical(blocks):
    return tuple(sorted(tuple(sorted(block)) for block in blocks))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",type=Path,default=RECORD)
    args = parser.parse_args()
    payload = json.loads(args.record.read_text())
    assert payload["source_sha256"]=={
        p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    assert len({(tuple(r["signature"]),r["n"]) for r in payload["records"]})==len(payload["records"])
    expected = dict(zip(map(name,SIGNATURES),
        (PAIR1,PAIR2,T222,PAIR4,C4,C5,C6,MODEL_C7,MODEL_C8,MODEL_J42,
         MODEL_J52,MODEL_J62,MODEL_J422,MODEL_J44)))
    assert set(payload["certificates"])==set(expected)
    for sig in SIGNATURES:
        rows = sorted((r for r in payload["records"] if tuple(r["signature"])==sig),key=lambda r:r["n"])
        certificate = certify(sig,rows)
        assert certificate==payload["certificates"][name(sig)]
        assert F(certificate["value"])==expected[name(sig)]
        coefficients = list(map(F,certificate["polynomial"]))
        assert evaluate(coefficients,0)==evaluate(coefficients,1)==evaluate(coefficients,-1)==0
        assert all(not a for i,a in enumerate(coefficients) if i%2!=(sum(sig)+1)%2)
        assert (coefficients[-1]*math.factorial(sum(sig)+1)).denominator==1
        print(f"{name(sig)}: {certificate['value']}, {sum(sig)//2} samples + held-out n={sum(sig)//2+2}")

    count = orbit_checks = 0
    files = ("paired_cycle_flow_2026-10-05.json", "pure_cycle_flow_2026-10-05.json",
             "mixed_cycle_flow_2026-10-05.json", "pure_cycle_eight_2026-10-05.json")
    for filename in files:
        old = json.loads((ROOT/"results"/filename).read_text())
        for row in old["records"]:
            sig = tuple(row["signature"]) if "signature" in row else ((2,)*row["k"] if "k" in row else (row["b"],))
            if name(sig) not in payload["certificates"]:
                continue
            cert = payload["certificates"][name(sig)]
            assert evaluate(list(map(F,cert["polynomial"])),row["n"])==row["signed_count"]
            count += 1
            lookup = {canonical(o["blocks"]):o for o in cert["orbit_values"]}
            for orbit in row.get("orbits",[]):
                candidate = lookup[canonical(orbit.get("blocks",orbit.get("pairs")))]
                assert candidate["multiplicity"]==orbit["multiplicity"]
                assert evaluate(list(map(F,candidate["polynomial"])),row["n"])==orbit["count"]
                orbit_checks += 1
    # The spectator method used a distinct weighted lattice and period-six bound.
    # The direct network method identifies each of its distance components.
    distances = {}
    for orbit in payload["certificates"]["5,2"]["orbit_values"]:
        pair = next(block for block in orbit["blocks"] if len(block)==2)
        distance = min((pair[1]-pair[0])%7,(pair[0]-pair[1])%7)
        assert orbit["multiplicity"]==7
        distances[distance] = F(orbit["value"])
    assert distances=={1:F(5,504),2:F(1,360),3:F(13,2520)}
    print(f"All {count} earlier aggregate counts and {orbit_checks} earlier orbit counts agree.")
    print("All fourteen reduced certificates and the three spectator components passed.")


if __name__=="__main__":
    main()
