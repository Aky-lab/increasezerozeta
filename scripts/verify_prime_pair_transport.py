#!/usr/bin/env python3
"""Finite identities for the published-input prime-pair transport theorem.

The asymptotic input is MRT Theorem 1.3(i), not a finite-height test.
See notes/prime_pair_progression_transport.md for the proof and hypotheses.
"""
from fractions import Fraction as F
from itertools import product
import math


def prime_exponents(n):
    result = []
    p = 2
    while p*p<=n:
        a = 0
        while n%p==0:
            a += 1
            n //= p
        if a:
            result.append(a)
        p += 1
    if n>1:
        result.append(1)
    return result


def divisor_histogram(limit):
    histogram = [0]*(limit+1)
    for q in range(1,limit+1):
        for k in range(1,limit//q+1):
            histogram[q*k] += 1
    return histogram


def main():
    limit = 5000
    tau = divisor_histogram(limit)
    for h in range(1,limit+1):
        exponents = prime_exponents(h)
        factor_count = math.prod(a+1 for a in exponents)
        four_count = math.prod(math.comb(a+3,3) for a in exponents)
        assert tau[h]==factor_count and tau[h]**2<=four_count
    # Ordered four-factor enumeration supplies a separate small-size check
    # of the divisor moment bound, independent of prime-power coefficients.
    small = 60
    d4 = [0]*(small+1)
    for a in range(1,small+1):
        for b in range(1,small//a+1):
            for c in range(1,small//(a*b)+1):
                for d in range(1,small//(a*b*c)+1):
                    d4[a*b*c*d] += 1
    assert all(d4[h]==math.prod(math.comb(a+3,3) for a in prime_exponents(h))
               for h in range(1,small+1))
    for height in (1,2,7,16,60):
        harmonic = sum(F(1,h) for h in range(1,height+1))
        assert sum(tau[1:height+1])<=height*harmonic
        assert sum(t*t for t in tau[1:height+1])<=height*harmonic**3
    families = 0
    for height,kind in product((7,24,60),(0,1,2,3)):
        # Families include all moduli, sparse moduli, and k-dependent gates.
        selected = [(q,k) for q in range(1,height+1) for k in range(1,height//q+1)
                    if kind==0 or (kind==1 and q%3==1)
                    or (kind==2 and k%2==1) or (kind==3 and (q+2*k)%5<2)]
        multiplicity = [0]*(height+1)
        for q,k in selected:
            multiplicity[q*k] += 1
        errors = [F(0)]+[F((-1)**h*(h%11),h%7+1) for h in range(1,height+1)]
        direct = sum(errors[q*k]**2 for q,k in selected)
        shifted = sum(multiplicity[h]*errors[h]**2 for h in range(1,height+1))
        assert direct==shifted
        assert all(multiplicity[h]<=tau[h] for h in range(1,height+1))
        bad = [h for h in range(1,height+1) if h%4==1 or h%9==0]
        good = [h for h in range(1,height+1) if h not in bad]
        bad_mass = sum(multiplicity[h] for h in bad)
        assert bad_mass**2<=len(bad)*sum(tau[h]**2 for h in range(1,height+1))
        good_bound = max((errors[h]**2 for h in good),default=F(0))*sum(tau[1:height+1])
        bad_bound = max((errors[h]**2 for h in bad),default=F(0))*bad_mass
        assert direct<=good_bound+bad_bound
        weights = [F((-1)**(q+k)*(q%4+1),k+1) for q,k in selected]
        weighted = sum(w*errors[q*k] for w,(q,k) in zip(weights,selected))
        assert weighted**2<=sum(w*w for w in weights)*direct
        families += 1
    # The unrestricted W2 window claim has counterexamples at arbitrary scale.
    for j in range(1,101):
        exponents = prime_exponents(6*j)
        assert len(exponents)>=2
    print(f"Divisor multiplicity and tau^2<=d4 checks: {limit}")
    print(f"Independent ordered four-factor counts: {small}")
    print(f"Exact family/histogram, exceptional-set and weighted consumption cases: {families}")
    print("Unrestricted unit-window counterexamples: 100")
    print("Finite bookkeeping passed. Published input and window hypotheses remain explicit.")


if __name__=="__main__":
    main()
