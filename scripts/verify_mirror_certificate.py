#!/usr/bin/env python3
"""Exact mirror-certificate identities, pole isolation and residue signs."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import permutations
import json
from math import isqrt
from pathlib import Path

Q = [F(2519), F(-8232), F(7368), F(-1932)]
D = list(map(F, (6345361, 104885808, 86095872, 3732624)))
N = list(map(F, (41472816, 131040168, 28469952)))
P = [a*(-1)**j for j, a in enumerate(D)]
BRACKETS = [(F(1,16), F(1,15)), (F(1), F(2)), (F(20), F(24))]


def add(a, b):
    return [sum((p[j] if j < len(p) else F(0)) for p in (a,b))
            for j in range(max(len(a),len(b)))]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j+k] += x*y
    return out


def scale(a, c):
    return [c*x for x in a]


def evaluate(a, x):
    out = F(0)
    for c in reversed(a):
        out = out*x+c
    return out


def iadd(a,b):
    return a[0]+b[0], a[1]+b[1]


def imul(a,b):
    v = [x*y for x in a for y in b]
    return min(v), max(v)


def idiv(a,b):
    assert b[0]*b[1] > 0
    return imul(a,(1/b[1],1/b[0]))


def ievaluate(a,x):
    out = F(0),F(0)
    for c in reversed(a):
        out = iadd(imul(out,x),(c,c))
    return out


def sqrt_interval(v):
    scale_int = 10**18
    lo = isqrt(v[0].numerator*scale_int**2//v[0].denominator)
    hi = isqrt(v[1].numerator*scale_int**2//v[1].denominator)+1
    assert F(lo,scale_int)**2 <= v[0] <= v[1] <= F(hi,scale_int)**2
    return F(lo,scale_int),F(hi,scale_int)


def poles():
    records = []
    intervals = []
    for lo,hi in BRACKETS:
        original = lo,hi
        signs = evaluate(P,lo),evaluate(P,hi)
        assert signs[0]*signs[1] < 0
        for _ in range(80):
            mid = (lo+hi)/2
            if evaluate(P,lo)*evaluate(P,mid) < 0:
                hi = mid
            else:
                lo = mid
        v = lo,hi
        assert evaluate(P,lo)*evaluate(P,hi) < 0
        numerator = ievaluate([N[0],-N[1],N[2]],v)
        derivative = ievaluate([D[1],-2*D[2],3*D[3]],v)
        weight = idiv(numerator,derivative)
        assert weight[0] > 0
        r = sqrt_interval(v)
        intervals.append((r,weight))
        records.append(dict(initial_bracket=list(map(str,original)),
                            endpoint_values=list(map(str,signs)),
                            squared_pole_interval=list(map(str,v)),
                            pole_interval=list(map(str,r)),
                            weight_interval=list(map(str,weight))))
    # Three disjoint sign-changing intervals exhaust the degree-three P.
    contribution = F(0),F(0)
    for r,w in intervals[1:]:
        r2,r3 = imul(r,r),imul(imul(r,r),r)
        lower = iadd(idiv((F(1),F(1)),r2),
                     idiv((F(-2,3),F(-2,3)),r3))
        contribution = iadd(contribution,imul(w,lower))
    target = idiv(iadd((F(9,10),F(9,10)),(-contribution[1],-contribution[0])),
                  intervals[0][1])
    assert F(104388,100000) < target[0] <= target[1] < F(104389,100000)
    return records,contribution,target


def mmul(a,b):
    return [[sum(a[j][k]*b[k][i] for k in range(3)) for i in range(3)] for j in range(3)]


def minverse(a):
    m = [list(row)+[F(i==j) for j in range(3)] for i,row in enumerate(a)]
    for j in range(3):
        k = next(k for k in range(j,3) if m[k][j])
        m[j],m[k] = m[k],m[j]
        pivot = m[j][j]
        m[j] = [x/pivot for x in m[j]]
        for k in range(3):
            if k!=j:
                c = m[k][j]
                m[k] = [x-c*y for x,y in zip(m[k],m[j])]
    return [row[3:] for row in m]


def matrix_polynomial(a,m):
    out = [[F(0)]*3 for _ in range(3)]
    power = [[F(i==j) for j in range(3)] for i in range(3)]
    for c in a:
        out = [[out[i][j]+c*power[i][j] for j in range(3)] for i in range(3)]
        power = mmul(power,m)
    return out


def partial_fractions():
    # Companion-matrix identity proves the complete rational decomposition,
    # using polynomial coefficients rather than sampled evaluation points.
    m = [[F(0),F(0),-D[0]/D[3]], [F(1),F(0),-D[1]/D[3]],
         [F(0),F(1),-D[2]/D[3]]]
    dp = matrix_polynomial([D[1],2*D[2],3*D[3]],m)
    w = mmul(matrix_polynomial(N,m),minverse(dp))
    a = [[[-m[i][j],F(i==j)] for j in range(3)] for i in range(3)]
    determinant = [F(0)]
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        term = [F(sign)]
        for i in range(3):
            term = mul(term,a[i][perm[i]])
        determinant = add(determinant,term)
    assert determinant == scale(D,1/D[3])
    adj = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rows = [k for k in range(3) if k!=j]
            cols = [k for k in range(3) if k!=i]
            minor = add(mul(a[rows[0]][cols[0]],a[rows[1]][cols[1]]),
                        scale(mul(a[rows[0]][cols[1]],a[rows[1]][cols[0]]),-1))
            adj[i][j] = scale(minor,(-1)**(i+j))
    numerator = [F(0)]
    for i in range(3):
        for j in range(3):
            numerator = add(numerator,scale(adj[j][i],w[i][j]))
    assert numerator == scale(N,1/D[3])
    assert sum(w[i][i] for i in range(3)) == N[2]/D[3]
    return dict(companion_determinant=list(map(str,determinant)),
                trace_adjugate_numerator=list(map(str,numerator)),
                sum_weights=str(N[2]/D[3]))


def scalar_identities():
    square = mul(Q,Q)
    assert square == list(map(F,(6345361,-41472816,104885808,-131040168,
                                86095872,-28469952,3732624)))
    assert square[::2] == D and scale(square[1::2],-1) == N
    f_one = evaluate(square,F(1))/evaluate(D,F(1))
    assert f_one == F(76729,201059665)
    assert f_one-F(1,10) == F(-8011695,80423866)
    # Clearing the positive 2*r^3*(x^2+r^2) denominator gives
    # 2*r^3*x - 2*r*x*(x^2+r^2) + x^2*(x^2+r^2)
    # = x^2*(x-r)^2. Compare every x-polynomial coefficient in Q[r].
    # Independently expand (x^2+r^2)*(x^2-2*r*x)+2*r^3*x.
    expanded = [[],[],[],[],[]]
    for i,ai in enumerate(([F(0),F(0),F(1)],[F(0)],[F(1)])):
        for j,bj in enumerate(([F(0)],[F(0),F(-2)],[F(1)])):
            expanded[i+j] = add(expanded[i+j],mul(ai,bj))
    expanded[1] = add(expanded[1],[F(0),F(0),F(0),F(2)])
    # Independently square x-r, then multiply by x^2.
    rhs = [[F(0)],[F(0)],[F(0)],[F(0)],[F(0)]]
    for i,ai in enumerate(([F(0),F(-1)],[F(1)])):
        for j,bj in enumerate(([F(0),F(-1)],[F(1)])):
            rhs[i+j+2] = add(rhs[i+j+2],mul(ai,bj))
    def trim(p):
        p = list(p)
        while len(p)>1 and p[-1]==0:
            p.pop()
        return p
    assert list(map(trim,expanded)) == list(map(trim,rhs))
    block_cap = sum(weight*evaluate(square,x)/evaluate(D,x*x)
                    for x,weight in ((F(0),F(1,6)),(F(1),F(2,3)),(F(2),F(1,6))))
    assert block_cap >= F(1,6) > F(1,10)
    return dict(square_coefficients=list(map(str,square)),
                value_at_one=str(f_one),prime_removal_threshold=str(f_one-F(1,10)),
                low_moment_block_certificate_value=str(block_cap))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    isolated,contribution,target = poles()
    record = dict(schema_version=1,generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scalar=scalar_identities(),partial_fractions=partial_fractions(),
                  isolated_poles=isolated,
                  large_pole_moment_contribution_interval=list(map(str,contribution)),
                  sufficient_small_pole_target_interval=list(map(str,target)),
                  arithmetic_resolvent_cap_proved=False,
                  scope='Exact polynomial and interval certificates; no prime-phase or covariance estimate.')
    if args.output:
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('PASS: mirror square identities; three simple imaginary pole pairs; positive weights;')
    print('complete rational partial fractions; quadratic lower bound; isolated one-point target.')
    print('The actual arithmetic resolvent inequality remains unproved.')


if __name__=='__main__':
    main()
