#!/usr/bin/env python3
"""Exact finite audits supporting the averaged rectangle Euler-tail proof.

The infinite estimates and prime-uniformity input are proved/cited in the
notes; these tests independently check their arithmetic bookkeeping.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product, permutations
from math import prod, lcm


def primes(height):
    return [p for p in range(2, height+1)
            if all(p % d for d in range(2, int(p**0.5)+1))]


def forms(k, h):
    return k, h, k-h, k+h


def local(p, k, h):
    distinct = len({a % p for a in (0, -k, h, h-k)})
    return Q(p-distinct, p) / Q(p-1, p)**4


def generic(p):
    return Q(p-4, p) / Q(p-1, p)**4


def check_factors():
    count = 0
    for k, h in product(range(-9, 10), repeat=2):
        values = forms(k, h)
        if not all(values):
            continue
        assert local(2, k, h)*local(3, k, h) <= 27
        for p in primes(31)[2:]:
            distinct = len({a % p for a in (0, -k, h, h-k)})
            assert (distinct < 4) == (prod(values) % p == 0)
            factor = Q(p-distinct, p-4)
            assert 0 < generic(p) < 1
            assert local(p, k, h) == generic(p)*factor
            assert 1 <= factor <= 1+Q(15, p)
            count += 1
        for w, z in ((5, 13), (7, 19), (11, 31)):
            head = prod(local(p, k, h) for p in primes(w))
            full = prod(local(p, k, h) for p in primes(z))
            tail = prod(generic(p)*Q(p-len({a % p for a in (0, -k, h, h-k)}), p-4)
                        for p in primes(z) if p > w)
            assert full == head*tail
            envelope = 27*prod(prod(1+Q(15, p) for p in primes(2*9) if n % p == 0)
                               for n in values)
            assert 0 <= head <= envelope and 0 <= full <= envelope
    return count


def check_divisor_moments():
    count = 0
    for height in (1, 8, 24):
        ps = primes(height)
        for exponent in (1, 2, 4, 16):
            coefficients = {1: Q(1)}
            for p in ps:
                coefficient = (1+Q(15, p))**exponent-1
                coefficients.update({d*p: g*coefficient for d, g in list(coefficients.items())})
            direct = sum(prod(1+Q(15, p) for p in ps if n % p == 0)**exponent
                         for n in range(1, height+1))
            expanded = sum(g*(height//d) for d, g in coefficients.items())
            upper = height*prod(1+((1+Q(15, p))**exponent-1)/p for p in ps)
            assert direct == expanded <= upper
            count += 1
    return count


def partitions(size):
    """Generate set partitions by adding the next labeled slot."""
    if size == 0:
        yield ()
        return
    for partition in partitions(size-1):
        for index in range(len(partition)):
            yield partition[:index]+(partition[index]+(size-1,),)+partition[index+1:]
        yield partition+((size-1,),)


def check_prime_moments():
    count = 0
    for height, w in ((24, 5), (40, 7)):
        ps = [p for p in primes(height) if p > w]
        for exponent, bell in ((1, 1), (2, 2), (4, 15)):
            direct = sum(sum((Q(1, p) for p in ps if n % p == 0), Q(0))**exponent
                         for n in range(1, height+1))
            expanded = sum(Q(height//lcm(*indices), prod(indices))
                           for indices in product(ps, repeat=exponent))
            density = sum(Q(1, prod(indices)*lcm(*indices))
                          for indices in product(ps, repeat=exponent))
            partition_density, unrestricted = Q(0), Q(0)
            parts = list(partitions(exponent))
            assert len(parts) == bell
            for partition in parts:
                sizes = [len(block) for block in partition]
                partition_density += sum(Q(1, prod(p**(b+1) for p, b in zip(indices, sizes)))
                                         for indices in permutations(ps, len(sizes)))
                unrestricted += prod(sum(Q(1, p**(b+1)) for p in ps) for b in sizes)
            assert direct == expanded <= height*density
            assert density == partition_density <= unrestricted <= Q(bell, w**exponent)
            # The all-equal tuple contributes and must not be silently discarded.
            assert density >= sum(Q(1, p**(exponent+1)) for p in ps)
            count += 1
    return count


def check_geometry_and_rare_shifts():
    geometry = 0
    for y in (1, 2, 5, 9):
        counters = [Counter() for _ in range(4)]
        for k, h in product(range(-y, y+1), repeat=2):
            values = forms(k, h)
            if all(values):
                for counter, value in zip(counters, values):
                    counter[abs(value)] += 1
        for counter in counters:
            assert all(1 <= n <= 2*y and multiplicity <= 2*(2*y+1)
                       for n, multiplicity in counter.items())
            geometry += 1
    rare = 0
    for w, z, upper in ((5, 13, 23), (7, 19, 31), (11, 23, 37)):
        k = 6*prod(p for p in primes(z) if p > w)
        h = 2*k
        assert all(forms(k, h))
        assert local(2, k, h)*local(3, k, h) == 27
        for p in primes(upper):
            expected = Q(p, p-1)**3 if p in (2, 3) or w < p <= z else generic(p)
            assert local(p, k, h) == expected
        ratio = prod(local(p, k, h) for p in primes(upper) if p > w)
        predicted = prod(Q(p, p-1)**3 for p in primes(z) if p > w)
        predicted *= prod(generic(p) for p in primes(upper) if p > z)
        assert ratio == predicted
        rare += 1
    return geometry, rare


if __name__ == '__main__':
    print('Exact local collision and discriminant checks:', check_factors())
    print('Positive divisor expansions and endpoint bounds:', check_divisor_moments())
    print('Repeated-prime moment and partition identities:', check_prime_moments())
    print('Nondegenerate geometry and rare-shift decompositions:', check_geometry_and_rare_shifts())
    print('Rectangle Euler-tail finite audits passed; infinite estimates use the written proof.')
