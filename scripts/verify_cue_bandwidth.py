#!/usr/bin/env python3
"""Exact rectangular CUE checks: Weyl integration versus trace cumulants.

No random samples are used. The all-order limit argument is in
notes/cue_gram_bandwidth.md; finite tests check its normalizations and
second-moment variance, rather than purporting to prove the limit.
"""
from fractions import Fraction as F
from itertools import product
from math import comb

from verify_cue_fluctuations import (
    WeylJoint, partitions, component_count, trace_cumulant,
)


def connected(sizes,m,n):
    length = sum(sizes)
    successors = []
    start = 0
    for size in sizes:
        successors.extend(start+(i+1)%size for i in range(size))
        start += size
    pi_connected = [pi for pi in partitions(length)
                    if component_count(pi,sizes)==1]
    total = 0
    for x in product(range(m),repeat=length):
        frequencies = tuple(x[j]-x[i] for i,j in enumerate(successors))
        for pi in pi_connected:
            term = 1
            for block in pi:
                term *= trace_cumulant(tuple(sorted(frequencies[i] for i in block)),n)
                if not term:
                    break
            total += term
    return F(total,m**len(sizes)*n**length)


def exact_variance(m,n):
    d = 2*m-n
    crossing = comb(d+2,5) if d>=3 else 0
    return F(4,m*m*n**4)*(F(m**5-m,30)-crossing)


def pair_variance(m,n):
    diagonal = sum((m-k)**2*k*k for k in range(1,m))
    crossing = sum((m-k)*(m-l)*max(k+l-n,0)
                   for k,l in product(range(1,m),repeat=2))
    return F(4*(diagonal-crossing),m*m*n**4)


def limit_variance(lam):
    return F(2,15)*lam**3-max(2*lam-1,F(0))**5/(30*lam**2)


def main():
    cases = 0
    for n in range(1,6):
        for m in range(1,n+1):
            weyl = WeylJoint(n,max_order=3,rows=m)
            assert weyl.moment((1,))==1
            assert weyl.moment((2,))==1+F(m*m-1,3*n*n)
            assert weyl.cumulant((2,2))==exact_variance(m,n)
            for sizes in ((2,),(3,),(2,2)):
                assert weyl.cumulant(sizes)==connected(sizes,m,n),(m,n,sizes)
                cases += 1
            if n<=4:
                assert weyl.cumulant((2,3))==connected((2,3),m,n)
                cases += 1
    # Independent finite sum identities, including both sides of 2m=n.
    sum_cases = 0
    for n in range(1,41):
        for m in range(1,n+1):
            assert exact_variance(m,n)==pair_variance(m,n)
            assert exact_variance(m,n)>=0
            if m==1:
                assert exact_variance(m,n)==0
            sum_cases += 1
    assert limit_variance(F(1))==F(1,10)
    assert limit_variance(F(1,2))==F(1,60)
    assert limit_variance(F(1,4))==F(1,480)
    for lam in (F(i,100) for i in range(1,101)):
        assert limit_variance(lam)>=lam**3/10>0
    # At rational aspect ratio use n=b*t,m=a*t. This checks the explicit
    # asymptotic formula against the finite one, with O(1/n) correction.
    for a,b in ((1,4),(1,2),(2,3),(3,4),(1,1)):
        lam = F(a,b)
        for t in (100,200,400):
            m,n = a*t,b*t
            assert abs(n*exact_variance(m,n)-limit_variance(lam))<=F(10,n)
    print('Exact rectangular Weyl/connected checks:',cases)
    print('Finite variance double-sum identities:',sum_cases)
    print('Bandwidth variance anchors, positivity and rational limits passed.')


if __name__=='__main__':
    main()
