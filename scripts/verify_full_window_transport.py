#!/usr/bin/env python3
"""Exact normalization audits for quantitative full-window prime transport."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import gcd, prod


def cube_sum(a, modulus=None):
    """Choose the base and three adjacent vertices, then recover the cube."""
    result = 0
    for x, p, q, s in product(a, repeat=4):
        nodes = (x, p, q, s, p+q-x, p+s-x, q+s-x, p+q+s-2*x)
        if modulus:
            nodes = tuple(n % modulus for n in nodes)
        result += prod(a.get(n, 0) for n in nodes)
    assert result>=0
    return result


def lock_count(f, g, r, q, h):
    result = defaultdict(int)
    for m, n, k in product(f, g, range(-h, h+1)):
        result[k, r*n-q*m] += f[m]*f.get(m-r*k, 0)*g[n]*g.get(n-q*k, 0)
    return result


def main():
    prefixes = progressions = locks = weights = exponents = 0
    for endpoint in range(2, 10):
        full = {x: ((3*x*x+2*x) % 7)-3 for x in range(1, endpoint+1)}
        size = 8*endpoint+1
        assert cube_sum(full)==cube_sum(full, size)
        for start in range(endpoint):
            lower = {x: v for x, v in full.items() if x<=start}
            window = {x: v for x, v in full.items() if x>start}
            assert all(full.get(x, 0)-lower.get(x, 0)==window.get(x, 0)
                       for x in range(size))
            assert cube_sum(window)<=2**7*(cube_sum(full)+cube_sum(lower))
            assert cube_sum(window)==cube_sum(window, size)
            prefixes += 1
    # Insufficient embedding can introduce genuine cyclic cubes.
    assert cube_sum({0: 1, 1: 1, 2: 1}, 3)>cube_sum({0: 1, 1: 1, 2: 1})
    for modulus, size in ((2, 12), (3, 12), (4, 12), (5, 15)):
        a = {x: ((x*x+2*x+1) % 7)-3 for x in range(size)}
        raw = cube_sum(a, size)
        for c in range(modulus):
            gate = {x: v for x, v in a.items() if x % modulus==c}
            fiber = {(x-c)//modulus: v for x, v in gate.items()}
            selected = cube_sum(gate, size)
            assert selected==cube_sum(fiber, size//modulus)<=raw
            assert F(selected, (size//modulus)**4)<=modulus**4*F(raw, size**4)
            progressions += 1
    for r, q in ((1, 1), (2, 3), (4, 3), (5, 2)):
        h, i, j = 2, r, 2*q
        f = {x: x % 3+1 for x in range(i+1, i+r*h+1)}
        g = {x: (2*x) % 3+1 for x in range(j+1, j+q*h+1)}
        u, v = {x: 1 for x in f}, {x: 1 for x in g}
        actual, model = lock_count(f, g, r, q, h), lock_count(u, v, r, q, h)
        keys = actual.keys() | model.keys()
        error = {key: actual[key]-model[key] for key in keys}
        variance = sum(e*e for e in error.values())
        differences, norms = [], []
        for first, second, step in ((f, u, r), (g, v, q)):
            for c in range(step):
                ff = {(x-c)//step: a for x, a in first.items() if x % step==c}
                uu = {(x-c)//step: a for x, a in second.items() if x % step==c}
                differences.append(cube_sum({x: ff[x]-uu[x] for x in ff}))
                norms.extend((cube_sum(ff), cube_sum(uu)))
        assert variance**4 <= (16*r*q)**4*max(differences)*max(norms)**3
        assert all(abs(k)<=h and abs(lock-(r*j-q*i))<=r*q*h for k, lock in keys)
        for mass in (F(1), F(3,7), F(5)):
            b = {key: F((key[0]+key[1]) % 5-2) for key in keys}
            consumed = sum(b[key]*error[key] for key in keys)
            assert consumed**2 <= mass*variance*sum(value**2/mass for value in b.values())
            weights += 1
        locks += 1
    for c0 in (F(1,100), F(1,3), F(1), F(3)):
        a, kappa = min(c0/2, F(1,4)), min(c0/2, F(1,4))/2
        assert 0<kappa<=F(1,8)<F(1,2)
        assert kappa/2-a<0  # progression discrepancy sqrt(R)*D
        assert kappa-2*a < -kappa  # spare exponent absorbs (log L2)^6
        assert -1 < -kappa  # spare exponent absorbs the head factor
        exponents += 1
    for w in range(2, 20):
        primes = [p for p in range(2, 21) if all(p % d for d in range(2, p))]
        assert [p for p in primes if p<F(w)+F(1,2)]==[p for p in primes if p<=w]
    print(f'PASS: {prefixes} prefix/embedding cases, {progressions} subgroup cases,')
    print(f'      {locks} balanced lock variances, {weights} weighted bounds, {exponents} exponent audits.')


if __name__=='__main__':
    main()
