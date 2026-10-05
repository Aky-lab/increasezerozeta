#!/usr/bin/env python3
"""Independent exact Haar-unitary Gram moments from Weyl constant terms.

Uses integer Laurent polynomials, not cumulants, network charts or orbits.
See notes/cue_gram_model.md. Only the Python standard library is required.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
import json
import math
from pathlib import Path

from model_moments import MODEL_MOMENTS
from reciprocity_certificates import ROOT, evaluate


def finite_polynomial(order, certificates):
    """Assemble singleton deletions in numerator n^(order+1) M_order(n)."""
    result = [F(0)]*(order+2)
    result[-1] = 1
    for signature, cert in certificates.items():
        size = sum(map(int,signature.split(",")))
        if size<=order:
            for i,a in enumerate(map(F,cert["polynomial"])):
                result[i+order-size] += math.comb(order,size)*a
    return result


def sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def weyl_moments(n,max_order):
    """Exact normalized E Tr((V V*)^b), V_aj=z_j^a/sqrt(n)."""
    zero = (0,)*(n-1)
    # Joint eigenangle density: |det(z_j^a)|^2 / n!.
    density = defaultdict(int)
    perms = list(permutations(range(n)))
    for p in perms:
        for q in perms:
            density[tuple(p[j]-q[j] for j in range(n-1))] += sign(p)*sign(q)
    density = {e:c for e,c in density.items() if c}
    assert density[zero]==math.factorial(n)
    entries = []
    for a in range(n):
        row = []
        for b in range(n):
            poly = defaultdict(int)
            for j in range(n):
                exponent = tuple(a-b if k==j else 0 for k in range(n-1))
                poly[exponent] += 1
            row.append(dict(poly))
        entries.append(row)
    power = [[{zero:int(a==b)} for b in range(n)] for a in range(n)]
    moments = [F(1)]
    for order in range(1,max_order+1):
        new = [[defaultdict(int) for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for k in range(n):
                    for e,c in power[a][k].items():
                        if c:
                            for f,d in entries[k][b].items():
                                new[a][b][tuple(x+y for x,y in zip(e,f))] += c*d
        power = [[dict(p) for p in row] for row in new]
        numerator = sum(c*density.get(tuple(-x for x in e),0)
                        for a in range(n) for e,c in power[a][a].items())
        moments.append(F(numerator,math.factorial(n)*n**(order+1)))
    return moments


def cue_two_moment(b):
    # Weyl's two-angle density integrates cos^(2j)(delta/2) to C_j/4^j.
    return sum(F(math.comb(b,2*j)*math.comb(2*j,j), (j+1)*4**j)
               for j in range(b//2+1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",type=Path,default=ROOT/"results/reciprocity_2026-10-05.json")
    args = parser.parse_args()
    payload = json.loads(args.record.read_text())
    certs = payload["certificates"]
    # n=2 closed-form and n=3,4 multivariate constant terms provide independent
    # full-moment checks, including every singleton and joint-class placement.
    independently = {n:weyl_moments(n,8) for n in (1,2,3,4)}
    for b in range(9):
        polynomial = finite_polynomial(b,certs)
        assert polynomial[-1]==MODEL_MOMENTS[b]
        assert all(not a for i,a in enumerate(polynomial) if i%2!=(b+1)%2)
        for n in (1,2,3,4):
            value = evaluate(polynomial,n)/n**(b+1)
            assert value==independently[n][b],(n,b,value,independently[n][b])
        assert independently[2][b]==cue_two_moment(b)
        print(f"b={b}: limit {polynomial[-1]}; exact CUE2 {independently[2][b]}; exact CUE3 {independently[3][b]}; exact CUE4 {independently[4][b]}")
    print("All finite CUE moments through eighth order agree with independent Weyl integration.")
    print("All corrections are even powers of 1/n. Arithmetic transport remains open.")


if __name__=="__main__":
    main()
