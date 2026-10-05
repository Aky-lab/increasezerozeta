#!/usr/bin/env python3
"""Exact finite checks for the conditional actual-matrix counting bridge.

These checks do not estimate the unbounded-height prime-side trace.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
import json
from math import comb
from pathlib import Path
import random

from christoffel_exact import christoffel_at_zero, det, hankel
from model_moments import MODEL_MOMENTS
from verify_progression_profile import closed_profile


def square(q):
    return [sum(q[i]*q[k-i] for i in range(len(q))
                if 0 <= k-i < len(q)) for k in range(2*len(q)-1)]


def shift(a, epsilon):
    return [sum(a[j]*comb(j,k)*(-epsilon)**(j-k)
                for j in range(k,len(a))) for k in range(len(a))]


def evaluate(a, x):
    out = F(0)
    for c in reversed(a):
        out = out*x+c
    return out


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a,b):
    return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]


def scale(a,c):
    return [[c*x for x in row] for row in a]


def multiply(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list,zip(*a)))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def outer(u):
    return [[x*y for y in u] for x in u]


def matrix_poly(q,a):
    out = scale(eye(len(a)),F(0))
    for c in reversed(q):
        out = add(multiply(out,a),scale(eye(len(a)),c))
    return out


def inertia(a):
    """Exact symmetric congruence elimination, including 2-by-2 pivots."""
    n = len(a)
    if not n:
        return (0,0,0)
    pivot = next((i for i in range(n) if a[i][i]),None)
    if pivot is not None:
        order = [pivot]+[i for i in range(n) if i != pivot]
        a = [[a[i][j] for j in order] for i in order]
        p = a[0][0]
        rest = [[a[i][j]-a[i][0]*a[0][j]/p
                 for j in range(1,n)] for i in range(1,n)]
        pos,neg,zero = inertia(rest)
        return (pos+(p>0),neg+(p<0),zero)
    pair = next(((i,j) for i in range(n) for j in range(i+1,n)
                 if a[i][j]),None)
    if pair is None:
        return (0,0,n)
    i,j = pair
    order = [i,j]+[k for k in range(n) if k not in pair]
    a = [[a[i][j] for j in order] for i in order]
    p = a[0][1]
    rest = [[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/p
             for j in range(2,n)] for i in range(2,n)]
    pos,neg,zero = inertia(rest)
    return (pos+1,neg+1,zero)


def matrix_checks(q):
    rng = random.Random(814079)
    # Known diagonal signatures transformed by a nonsingular congruence.
    calibration = 0
    for n in range(1,7):
        for index,signature in enumerate(product((-1,0,1),repeat=n)):
            if n > 3 and index % 13:
                continue
            u = [[F(i == j)+F(rng.randrange(-2,3) if i<j else 0)
                  for j in range(n)] for i in range(n)]
            diagonal = [[F(signature[i] if i==j else 0)
                         for j in range(n)] for i in range(n)]
            congruence = multiply(transpose(u),multiply(diagonal,u))
            assert inertia(congruence) == (signature.count(1),
                                           signature.count(-1),
                                           signature.count(0))
            calibration += 1
    assert inertia([[F(0),F(1)],[F(1),F(0)]]) == (1,1,0)

    cases = 0
    for seed in range(36):
        n = 2+seed%6
        s1,s2,pairs = seed%4,(seed//4)%3,(seed//12)%3
        a = scale(eye(n),F(0))
        multiplicity = s1
        def vector():
            return [F(rng.randrange(-2,3),2) for _ in range(n)]
        for _ in range(s1):
            a = add(a,outer(vector()))
        for _ in range(s2):
            m = rng.randrange(2,5)
            a = add(a,scale(outer(vector()),m))
            multiplicity += m
        for _ in range(pairs):
            m = rng.randrange(1,4)
            u,v = vector(),vector()
            # Include exact dependence and arbitrarily small negative parts.
            if seed%3 == 0:
                v = [x/F(1000) for x in u]
            a = add(a,scale(add(outer(u),scale(outer(v),-1)),2*m))
            multiplicity += 2*m
        pos,neg,zero = inertia(a)
        assert pos <= s1+s2+pairs
        assert pos+neg <= s1+s2+2*pairs
        assert s1 >= 2*pos-multiplicity
        assert s1+s2+2*pairs >= pos
        for epsilon in (F(0),F(1,1000),F(1,5)):
            perturbation = [[epsilon*F(rng.randrange(-1,2),4*n)
                             for _ in range(n)] for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    perturbation[i][j] = perturbation[j][i]
            assert max(sum(map(abs,row)) for row in perturbation) <= epsilon
            b = add(a,perturbation)
            shifted = add(b,scale(eye(n),-epsilon))
            above = inertia(shifted)[0]
            assert above <= pos
            qb = matrix_poly(q,shifted)
            qtrace = trace(multiply(qb,qb))
            assert qtrace == sum(x*x for row in qb for x in row)
            assert qtrace >= n-above
            powers = [eye(n)]
            for _ in range(2*(len(q)-1)):
                powers.append(multiply(powers[-1],b))
            coeff = shift(square(q),epsilon)
            assert sum(coeff[k]*trace(powers[k]) for k in range(len(coeff))) == qtrace
            assert s1 >= 2*(n-qtrace)-multiplicity
            cases += 1
    # Strict threshold must include an eigenvalue equal to epsilon.
    epsilon = F(1,5)
    b = [[epsilon,F(0)],[F(0),F(1)]]
    shifted = add(b,scale(eye(2),-epsilon))
    assert inertia(shifted) == (1,0,1)
    assert evaluate(q,F(0))**2 == 1
    return cases,calibration


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    lam,q,_ = christoffel_at_zero(MODEL_MOMENTS,4)
    d = 162540559
    assert q == [F(1),-F(851656632,d),F(1278096456,d),
                 -F(729534540,d),F(140399280,d)]
    assert lam == F(12241115,d)
    for n in range(5):
        assert det(hankel(MODEL_MOMENTS,n)) > 0
    coeff = square(q)
    assert all(c*(-1)**k > 0 for k,c in enumerate(coeff))
    assert sum(coeff[k]*MODEL_MOMENTS[k] for k in range(9)) == lam
    margin = F(1,10)-lam
    assert margin == F(40129409,1625405590)
    sensitivity = sum(map(abs,coeff[1:]))
    assert sensitivity == F(9973263119729203608,26419433320032481)
    tolerance = margin/sensitivity
    assert tolerance == F(6522656571199631,99732631197292036080)

    lam3,q3,_ = christoffel_at_zero(MODEL_MOMENTS,3)
    assert q3 == [F(1),-F(8232,2519),F(7368,2519),-F(1932,2519)]
    assert lam3 == F(247,2519)
    coeff3 = square(q3)
    assert coeff3 == [F(c,6345361) for c in
                      (6345361,-41472816,104885808,-131040168,
                       86095872,-28469952,3732624)]
    assert sum(coeff3[k]*MODEL_MOMENTS[k] for k in range(7)) == lam3
    assert all(c*(-1)**k > 0 for k,c in enumerate(coeff3))
    assert F(1,10)-lam3 == F(49,25190)
    assert 1-2*lam3 == F(2025,2519)

    shifts = 0
    for epsilon in (F(0),F(1,10000),F(1,100),F(1,4),F(1)):
        shifted = shift(coeff,epsilon)
        assert all(c*(-1)**k > 0 for k,c in enumerate(shifted))
        for x in (F(-3),F(-1,2),F(0),F(1,3),F(2),F(7)):
            assert evaluate(shifted,x) == evaluate(q,x-epsilon)**2
            if x <= epsilon:
                assert evaluate(shifted,x) >= 1
            shifts += 1
        model_shift = sum(shifted[k]*MODEL_MOMENTS[k] for k in range(9))
        # The fixed-point Christoffel extremum also bounds shifted polynomials.
        assert model_shift >= lam*evaluate(q,-epsilon)**2
        bad_moments = {k:MODEL_MOMENTS[k]+(-1)**k*tolerance
                       for k in range(1,9)}
        bad_moments[0] = F(1)
        excess = sum(shifted[k]*(bad_moments[k]-MODEL_MOMENTS[k])
                     for k in range(9))
        assert excess == sum(abs(shifted[k])*tolerance for k in range(1,9))

    finite,calibrations = matrix_checks(q)
    finite3,_ = matrix_checks(q3)
    collar_checks = 0
    for s1,s2,p in product(range(5),repeat=3):
        nprime = s1+2*s2+2*p
        for removed_simple in range(s1+1):
            for removed_multiple in range(s2+1):
                for removed_pairs in range(p+1):
                    collar = removed_simple+2*removed_multiple+2*removed_pairs
                    n = nprime-collar
                    if not n:
                        continue
                    possible_positive = s1+s2+p
                    assert s1-removed_simple >= 2*possible_positive-n-2*collar
                    distinct = s1+s2+2*p-removed_simple-removed_multiple-2*removed_pairs
                    assert distinct >= possible_positive-collar
                    collar_checks += 1

    atoms = ((F(1,9),F(-1,2)),(F(5,9),F(1)),(F(1,3),F(3,2)))
    assert sum(w for w,x in atoms) == 1
    assert sum(w*x for w,x in atoms) == MODEL_MOMENTS[1]
    assert sum(w*x*x for w,x in atoms) == MODEL_MOMENTS[2]
    counterexample_q = sum(w*evaluate(q,x)**2 for w,x in atoms)
    assert counterexample_q >= F(1,9) > lam

    exponent_checks = 0
    for nu in (F(1),F(5,4),F(3,2),F(2)):
        for b1,b2 in product((F(k,24) for k in range(25)),repeat=2):
            h = nu-b1-b2
            slot_range = b1 <= h/3 and b2 <= h/3
            physical_range = b1 <= (nu-b2)/4 and b2 <= (nu-b1)/4
            assert slot_range == physical_range
            exponent_checks += 1
    coverage = closed_profile(F(1,4))
    assert coverage == F(1109,144000)
    tail_exponents = [-1+F(k,24)/2 for k in range(1,25)]
    assert all(x <= F(-1,2) for x in tail_exponents)
    support_cases = 0
    for degree in range(2,9):
        maximum = max(sum(abs(x[i]-x[i-1]) for i in range(degree))
                      for x in product((0,1),repeat=degree))
        assert maximum == 2*(degree//2)
        for bandwidth in (F(1,4),F(1,2),F(1)):
            for vertex in product((-bandwidth/2,bandwidth/2),repeat=degree):
                frequency = [vertex[i]-vertex[i-1] for i in range(degree)]
                assert sum(frequency) == 0
                assert sum(map(abs,frequency)) <= maximum*bandwidth
                support_cases += 1
    interior = [F((-1)**i*3,8) for i in range(8)]
    assert sum(abs(interior[i]-interior[i-1]) for i in range(8)) == 6 > 2
    assert coeff[8] == q[4]**2 > 0

    report = {
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'finite certificate and counting checks; arithmetic trace target unproved',
        'q4_coefficients':[str(c) for c in q],
        'q3_coefficients':[str(c) for c in q3],
        'q3_squared_coefficients':[str(c) for c in coeff3],
        'q3_model_trace':str(lam3),
        'q3_conditional_simple_bound':str(1-2*lam3),
        'q3_excess_margin_for_80_percent':str(F(1,10)-lam3),
        'q3_matrix_perturbation_cases':finite3,
        'squared_coefficients':[str(c) for c in coeff],
        'model_trace':str(lam),
        'conditional_simple_bound':str(1-2*lam),
        'trace_cap_for_80_percent':str(F(1,10)),
        'excess_margin_for_80_percent':str(margin),
        'conservative_uniform_positive_order_moment_tolerance':str(tolerance),
        'matrix_perturbation_cases':finite,
        'inertia_signature_grid_cases':calibrations,
        'shifted_polynomial_cases':shifts,
        'collar_cases':collar_checks,
        'progression_exponent_cases':exponent_checks,
        'cyclic_fourier_support_cases':support_cases,
        'eighth_cycle_unit_bandwidth_support_radius':8,
        'quarter_profile_model_coverage':str(coverage),
        'low_moment_counterexample_trace':str(counterexample_q),
        'limitations':[
            'Finite checks do not establish any actual zeta trace asymptotic.',
            'The polynomial-trace cap and its prime-side transport remain unproved.',
            'The published large-progression input has exceptional moduli and a restricted range.',
            'Model coverage is not a proved fraction of an actual arithmetic error.',
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {finite} exact matrix/perturbation cases, {shifts} shifts,')
    print(f'      {collar_checks} collar cases and {exponent_checks} progression scales.')
    print('Conditional trace cap for 80%: 1/10; model trace:',lam)
    print('Unproved arithmetic target: actual polynomial trace <= that cap asymptotically.')


if __name__ == '__main__':
    main()
