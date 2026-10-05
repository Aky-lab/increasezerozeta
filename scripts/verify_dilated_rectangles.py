#!/usr/bin/env python3
"""Exact independent audits of heterogeneous dilated rectangle locks."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import gcd, prod


def locked(f, g, r, q):
    out = defaultdict(int)
    extent = max(max(f)-min(f), max(g)-min(g))
    for m, n, k in product(f, g, range(-extent, extent+1)):
        out[k, r*n-q*m] += f[m]*f.get(m-r*k, 0)*g[n]*g.get(n-q*k, 0)
    return out


def cube(f, g, r, q):
    extent = max(max(f)-min(f), max(g)-min(g))
    return sum(f[m]*f.get(m-r*k, 0)*f.get(m-r*t, 0)*f.get(m-r*(k+t), 0)
               *g[n]*g.get(n-q*k, 0)*g.get(n-q*t, 0)*g.get(n-q*(k+t), 0)
               for m, n, k, t in product(f, g, range(-extent, extent+1), range(-extent, extent+1)))


def integer_u3(a):
    if not a:
        return 0
    extent = max(a)-min(a)
    vertices = list(product(range(2), repeat=3))
    value = sum(prod(a.get(x+aa*h+bb*k+cc*t, 0) for aa, bb, cc in vertices)
                for x, h, k, t in product(a, range(-extent, extent+1),
                                           range(-extent, extent+1), range(-extent, extent+1)))
    assert value >= 0
    return value


def check_locked_geometry():
    f = dict(enumerate((-1, 2, 0, 3, 1), start=11))
    g = dict(enumerate((2, 1, -2, 1, 3), start=17))
    u = dict(enumerate((1, 0, -1, 2, 1), start=11))
    v = dict(enumerate((1, 1, 0, -1, 2), start=17))
    count = 0
    for r, q in product(range(1, 6), repeat=2):
        if gcd(r, q) != 1:
            continue
        raw, model = locked(f, g, r, q), locked(u, v, r, q)
        assert sum(value**2 for value in raw.values()) == cube(f, g, r, q)
        error = sum((raw[key]-model[key])**2 for key in raw.keys() | model.keys())
        discrepancies, norms = [], []
        for aa, bb, modulus in ((f, u, r), (g, v, q)):
            for residue in range(modulus):
                first = {(x-residue)//modulus: value for x, value in aa.items() if x % modulus == residue}
                second = {(x-residue)//modulus: value for x, value in bb.items() if x % modulus == residue}
                discrepancies.append(integer_u3({x: first[x]-second[x] for x in first}))
                norms.extend((integer_u3(first), integer_u3(second)))
        # Raise the Gowers inequality to the fourth power: cyclic M factors cancel.
        assert error**4 <= (16*r*q)**4*max(discrepancies)*max(norms)**3
        count += 1
    return count


def cyclic_u3(a):
    size = len(a)
    vertices = list(product(range(2), repeat=3))
    return sum(prod(a[(x+aa*h+bb*k+cc*t) % size] for aa, bb, cc in vertices)
               for x, h, k, t in product(range(size), repeat=4))


def check_progression_extraction():
    cases = 0
    for size, modulus in ((8, 2), (8, 4), (12, 3), (12, 4)):
        h = [((5*x*x+3*x+2) % 7)-3 for x in range(size)]
        full = cyclic_u3(h)
        assert full >= 0
        for residue in range(modulus):
            gated = [value if x % modulus == residue else 0 for x, value in enumerate(h)]
            fiber = h[residue::modulus]
            raw_fiber = cyclic_u3(fiber)
            assert cyclic_u3(gated) == raw_fiber <= full
            assert Q(raw_fiber, (size//modulus)**4) <= modulus**4*Q(full, size**4)
            indicator = [int(x % modulus == residue) for x in range(size)]
            assert Q(cyclic_u3(indicator), size**4) == Q(1, modulus**4)
            assert cyclic_u3([1]*(size//modulus)) == (size//modulus)**4
            cases += 1
    return cases


def admissible(modulus, r, q):
    return sum(all(gcd(value, modulus) == 1 for value in
                   (m, m-r*k, m-r*t, m-r*(k+t), n, n-q*k, n-q*t, n-q*(k+t)))
               for m, n, k, t in product(range(modulus), repeat=4))


def prime_count(p, r, q):
    if r*q % p == 0:
        return (p-1)*(p**3-4*p*p+6*p-3)
    return 1 if p == 2 else p**4-8*p**3+28*p*p-44*p+23


def local(p, r, q, k, lock):
    if r*q % p:
        c = -lock*pow(r*q, -1, p) % p
        density = 1-Q(len({0, k % p, c, (c+k) % p}), p)
    elif lock % p == 0:
        density = Q(0)
    else:
        density = 1-Q(1 if k % p == 0 else 2, p)
    return density/Q(p-1, p)**4


def period_counts(modulus, r, q):
    out = defaultdict(int)
    for m, n, k in product(range(modulus), repeat=3):
        if all(gcd(value, modulus) == 1 for value in (m, m-r*k, n, n-q*k)):
            out[k, (r*n-q*m) % modulus] += 1
    return out


def check_prime_and_crt_budgets():
    field_cases, crt_cases = 0, 0
    for p in (2, 3, 5, 7, 11):
        for r, q in ((1, 1), (2, 3), (p, 1), (1, p), (p*p, p+1)):
            direct = admissible(p, r, q)
            assert direct == prime_count(p, r, q)
            counts = period_counts(p, r, q)
            assert sum(value**2 for value in counts.values()) == direct
            beta = Q(direct*p**4, (p-1)**8)
            assert sum(local(p, r, q, k, lock)**2 for k, lock in product(range(p), repeat=2))/p**2 == beta
            for k, lock in product(range(p), repeat=2):
                assert local(p, r, q, k, lock) == Q(counts[k, lock], p)/Q(p-1, p)**4
            if r*q % p == 0:
                assert beta == Q(p, p-1)**3*(1+Q(1, (p-1)**3))
            field_cases += 1
    for ps in ((2, 3), (2, 5), (3, 5)):
        modulus = prod(ps)
        for r, q in ((1, 1), (2, 3), (4, 9), (5, 2), (7, 5)):
            counts = period_counts(modulus, r, q)
            assert sum(value**2 for value in counts.values()) == prod(prime_count(p, r, q) for p in ps)
            for k, lock in product(range(modulus), repeat=2):
                assert Q(counts[k, lock], modulus) == prod(local(p, r, q, k, lock)*Q(p-1, p)**4 for p in ps)
            crt_cases += 1
    assert admissible(6, 4, 9) == prime_count(2, 4, 9)*prime_count(3, 4, 9)
    return field_cases, crt_cases


def solution(r, q, lock):
    m0 = (-lock*pow(q, -1, r)) % r if r > 1 else 0
    n0 = (lock+q*m0)//r
    assert r*n0-q*m0 == lock
    return m0, n0


def overlap(r, q, k, lock, a, ym, b, yn, translate=0):
    m0, n0 = solution(r, q, lock)
    m0, n0 = m0+r*translate, n0+q*translate
    lower = max((a-m0)//r, (b-n0)//q)
    upper = min((a+ym-m0)//r, (b+yn-n0)//q)
    return max(upper-lower-abs(k), 0)


def check_windows():
    cases = 0
    for r, q in ((1, 1), (2, 3), (4, 9), (5, 2), (7, 5)):
        for a, ym, b, yn in ((-4, 7, 2, 9), (11, 9, 17, 6), (0, 3, 0, 4)):
            f = {m: 1 for m in range(a+1, a+ym+1)}
            g = {n: 1 for n in range(b+1, b+yn+1)}
            raw = locked(f, g, r, q)
            for k, lock in product(range(-3, 4), range(-12, 13)):
                z = overlap(r, q, k, lock, a, ym, b, yn)
                assert raw[k, lock] == z
                assert overlap(r, q, k, lock, a, ym, b, yn, 5) == z
                for ps in ((2, 3), (3, 5)):
                    modulus, c = prod(ps), prod(Q(p, p-1) for p in ps)
                    sifted = sum(1 for m, n in product(f, g)
                                 if r*n-q*m == lock and m-r*k in f and n-q*k in g
                                 and all(gcd(x, modulus) == 1 for x in (m, m-r*k, n, n-q*k)))
                    main = z*prod(local(p, r, q, k, lock) for p in ps)
                    assert abs(c**4*sifted-main) <= c**4*modulus
                cases += 1
    return cases


if __name__ == '__main__':
    print('Resolved squared counts, cubes and exact progression-norm criterion:', check_locked_geometry())
    print('Coset cube identity, sharp extraction loss and selected-progression examples:', check_progression_extraction())
    print('Prime-field and CRT local second-moment budgets:', check_prime_and_crt_budgets())
    print('Bezout-independent raw overlaps and presieved main terms:', check_windows())
    print('Dilated rectangle finite audits passed; growing-modulus prime uniformity remains open.')
