#!/usr/bin/env python3
"""Exact finite audits for weighted, slowly growing dilated prime transport.

No numerical samples prove the qualitative analytic theorem or its rate.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product
from math import gcd, isqrt, lcm, prod


def primes(height):
    return [p for p in range(2, height+1) if all(p % d for d in range(2, isqrt(p)+1))]


def check_translated_divisors():
    cases = 0
    for start, end in ((0, 24), (17, 31), (27, 31), (113, 127)):
        ps = primes(end)
        for exponent in (1, 2, 4, 16):
            coefficients = {1: Q(1)}
            for p in ps:
                g = (1+Q(15, p))**exponent-1
                coefficients.update({d*p: value*g for d, value in list(coefficients.items()) if d*p <= end})
            direct = sum(prod(1+Q(15, p) for p in ps if n % p == 0)**exponent
                         for n in range(start+1, end+1))
            expansion = sum(value*(end//d-start//d) for d, value in coefficients.items())
            main = (end-start)*prod(1+((1+Q(15, p))**exponent-1)/p for p in ps)
            endpoint = prod((1+Q(15, p))**exponent for p in ps)
            assert direct == expansion <= main+endpoint
            assert all(end//d-start//d <= Q(end-start, d)+1 for d in coefficients)
            cases += 1
    assert 2//2-1//2 > Q(1, 2)  # The unshifted floor bound fails after translation.
    return cases


def check_translated_prime_moments():
    cases = 0
    for start, end, cutoff in ((0, 24, 5), (17, 31, 5), (27, 31, 7)):
        ps = [p for p in primes(end) if p > cutoff]
        for exponent in (2, 4):
            direct = sum(sum((Q(1, p) for p in ps if n % p == 0), Q(0))**exponent
                         for n in range(start+1, end+1))
            expanded, density = Q(0), Q(0)
            for indices in product(ps, repeat=exponent):
                modulus = lcm(*indices)
                expanded += Q(end//modulus-start//modulus, prod(indices))
                density += Q(1, modulus*prod(indices))
            endpoint = sum((Q(1, p) for p in ps), Q(0))**exponent
            assert direct == expanded <= (end-start)*density+endpoint
            bell = 2 if exponent == 2 else 15
            assert density <= Q(bell, cutoff**exponent)
            cases += 1
    return cases


def bezout_solution(r, q, lock):
    m0 = (-lock*pow(q, -1, r)) % r if r > 1 else 0
    n0 = (lock+q*m0)//r
    return m0, n0


def check_discriminants_and_geometry():
    local_cases, boxes = 0, 0
    for r, q in ((1, 1), (2, 3), (4, 9), (5, 7), (8, 25)):
        for k, lock in product(range(-6, 7), repeat=2):
            values = (k, lock, r*q*k-lock, r*q*k+lock)
            if not all(values):
                continue
            m0, n0 = bezout_solution(r, q, lock)
            for p in (5, 7, 11, 13):
                if r*q % p == 0:
                    continue
                roots = {z for z in range(p) if any(value % p == 0 for value in
                         (m0+r*z, m0+r*(z-k), n0+q*z, n0+q*(z-k)))}
                assert (len(roots) < 4) == (prod(values) % p == 0)
                generic = Q(p-4, p)/Q(p-1, p)**4
                factor = Q(p-len(roots), p-4)
                local = Q(p-len(roots), p)/Q(p-1, p)**4
                assert local == generic*factor and 1 <= factor <= 1+Q(15, p)
                local_cases += 1
        for height, center in ((2, -100), (3, 17), (3, 100)):
            histograms = [Counter() for _ in range(4)]
            for k, lock in product(range(-height, height+1),
                                   range(center-r*q*height, center+r*q*height+1)):
                values = (k, lock, r*q*k-lock, r*q*k+lock)
                if all(values):
                    for histogram, value in zip(histograms, values):
                        histogram[abs(value)] += 1
            assert max(histograms[0].values(), default=0) <= 2*(2*r*q*height+1)
            assert all(max(h.values(), default=0) <= 2*(2*height+1) for h in histograms[1:])
            boxes += 1
    return local_cases, boxes


def check_carry_progression_splitting():
    cases = 0
    for modulus, r, q in ((2, 2, 3), (3, 4, 9), (6, 5, 2), (6, 8, 25)):
        for mres, nres, kres, tres in product(range(modulus), repeat=4):
            for mm, nn, kk, tt in ((-2, 3, -1, 2), (1, -1, 2, -2), (3, 4, 0, 1)):
                m, n = mres+modulus*mm, nres+modulus*nn
                k, t = kres+modulus*kk, tres+modulus*tt
                aa, bb = mm % r, nn % q
                z, y = (mm-aa)//r, (nn-bb)//q
                for first, second in product(range(2), repeat=2):
                    for base, step, residue, quotient_base, ap_residue, ap_base in (
                            (m, r, mres, mm, aa, z), (n, q, nres, nn, bb, y)):
                        raw = base-step*(first*k+second*t)
                        numerator = residue-step*(first*kres+second*tres)
                        carry, vertex_residue = divmod(numerator, modulus)
                        reconstructed = modulus*(quotient_base-step*(first*kk+second*tt)+carry)+vertex_residue
                        assert raw == reconstructed
                        ap_quotient = ap_residue+step*(ap_base-first*kk-second*tt)+carry
                        assert ap_quotient == quotient_base-step*(first*kk+second*tt)+carry
                        # Carry magnitude can grow with the step, but its support translation after division cannot.
                        assert -3 <= (ap_residue+carry)//step <= 1
                cases += 1
    return cases


def check_weighted_exception_and_consumption():
    cases = 0
    for size in (6, 10, 17):
        for mode in ('diagonal', 'all', 'nonuniform'):
            weights = {}
            for i, j in product(range(size), repeat=2):
                value = int(i == j) if mode == 'diagonal' else 1
                if mode == 'nonuniform':
                    value = (3*i+5*j) % 4
                if value:
                    weights[i, j] = Q(value)
            marginal1, marginal2 = defaultdict(Q), defaultdict(Q)
            for (i, j), weight in weights.items():
                marginal1[i] += weight
                marginal2[j] += weight
            bad1, bad2 = set(range(0, size, 4)), set(range(1, size, 5))
            badmass = sum(weight for (i, j), weight in weights.items() if i in bad1 or j in bad2)
            assert badmass <= max(marginal1.values())*len(bad1)+max(marginal2.values())*len(bad2)
            mass = sum(weights.values())
            if mode in ('diagonal', 'all'):
                assert max(marginal1.values()) == max(marginal2.values()) == mass/size
            errors = {key: Q((2*key[0]+key[1]) % 7-3, 5) for key in weights}
            coefficients = {key: Q(key[0]-2*key[1], 3) for key in weights}
            lhs = sum(coefficients[key]*errors[key] for key in weights)**2
            weighted_variance = sum(weights[key]*errors[key]**2 for key in weights)
            consumption = sum(coefficients[key]**2/weights[key] for key in weights)
            assert lhs <= weighted_variance*consumption
            cases += 1
    return cases


if __name__ == '__main__':
    print('Translated positive-divisor expansions and endpoint bounds:', check_translated_divisors())
    print('Translated repeated-prime moments and lcm endpoint counts:', check_translated_prime_moments())
    print('Dilated discriminant support and translated-box multiplicities:', check_discriminants_and_geometry())
    print('Heterogeneous W-residue, carry and progression splitting:', check_carry_progression_splitting())
    print('Weighted exceptional mass and consumption inequality:', check_weighted_exception_and_consumption())
    print('Weighted dilation finite audits passed; the asymptotic claims use the written proofs and published input.')
