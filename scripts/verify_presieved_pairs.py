#!/usr/bin/env python3
"""Exact finite comb, Rankin-factor and autocorrelation audits.

The prime asymptotic theorem is imported from short-interval uniformity.
These tests verify the finite algebra used in notes/presieved_pair_transport.md.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import math


def divisors(n):
    return [d for d in range(1,n+1) if n%d==0]


def mobius(n):
    result,p = 1,2
    while p*p<=n:
        if n%p==0:
            n //= p
            result = -result
            if n%p==0:
                return 0
        p += 1
    return -result if n>1 else result


def phi(n):
    return sum(math.gcd(a,n)==1 for a in range(1,n+1))


def ramanujan(q,h):
    return sum(d*mobius(q//d) for d in divisors(math.gcd(q,h)))


def check_comb():
    densities,intervals = 0,0
    for primes in ((2,3),(2,3,5),(2,3,7)):
        period = math.prod(primes)
        ds = divisors(period)
        c = F(period,phi(period))
        full = [c if math.gcd(n,period)==1 else F(0) for n in range(period)]
        assert sum(full)/period==1
        assert sum(v*v for v in full)/period==c
        for n in range(period):
            expanded = sum(F(mobius(q),phi(q))*ramanujan(q,n) for q in ds)
            assert expanded==full[n]
        for dmax in sorted({1,2,3,5,7,11,period}):
            kept = [d for d in ds if d<=dmax]
            alpha = {q:c*sum(F(mobius(d),d) for d in kept if d%q==0)
                     for q in kept}
            partial = [c*sum(mobius(d) for d in kept if n%d==0) for n in range(period)]
            for n in range(period):
                assert sum(alpha[q]*ramanujan(q,n) for q in kept)==partial[n]
            if dmax>=period:
                assert all(alpha[q]==F(mobius(q),phi(q)) for q in ds)
            if period==30 and dmax==5:
                assert alpha[1]==F(-1,8) and min(partial)<0
            for h in range(-period,period+1):
                direct = sum(partial[n]*partial[(n+h)%period] for n in range(period))/period
                spectral = sum(a*a*ramanujan(q,h) for q,a in alpha.items())
                crt = c*c*sum(F(mobius(d)*mobius(e),math.lcm(d,e))
                              for d,e in product(kept,repeat=2) if h%math.gcd(d,e)==0)
                assert direct==spectral==crt
                full_corr = sum(full[n]*full[(n+h)%period] for n in range(period))/period
                local = math.prod(1+F(ramanujan(p,h),(p-1)**2) for p in primes)
                assert full_corr==local
                if h==0:
                    assert full_corr==c
                elif h%2:
                    assert full_corr==0
                for start,length in ((0,1),(3,5),(period-2,9),(2*period+1,17)):
                    interval = sum(partial[n%period]*partial[(n+h)%period]
                                   for n in range(start+1,start+length+1))
                    assert abs(interval-length*crt)<=c*c*dmax*dmax
                    for d,e in product(kept,repeat=2):
                        count = sum(n%d==0 and (n+h)%e==0
                                    for n in range(start+1,start+length+1))
                        main = F(length,math.lcm(d,e)) if h%math.gcd(d,e)==0 else F(0)
                        assert abs(count-main)<=1
                    intervals += 1
                densities += 1
    return densities,intervals


def check_rankin_factors():
    cases = 0
    for primes in ((2,3),(2,3,5),(2,3,7)):
        period = math.prod(primes)
        ds = divisors(period)
        for values in (tuple(F(3,2) for _ in primes),tuple(F(p,1) for p in primes)):
            weight = {d:math.prod(z for p,z in zip(primes,values) if d%p==0) for d in ds}
            expanded = sum(F(weight[d]*weight[e],math.lcm(d,e)) for d,e in product(ds,repeat=2))
            euler = math.prod(1+(2*z+z*z)/p for p,z in zip(primes,values))
            assert expanded==euler
            for limit in (7,31,100):
                counted = sum(sum(weight[d] for d in ds if n%d==0)**2 for n in range(1,limit+1))
                assert counted<=limit*euler
                cases += 1
        # Sigma=1 gives rational tests of the Rankin cutoff inequality.
        for dmax in (1,2,5,11):
            for n in range(1,101):
                tail = sum(1 for d in ds if d>dmax and n%d==0)
                weighted = sum(d for d in ds if n%d==0)
                assert tail<=F(weighted,dmax)
    return cases


def multiply(a,b):
    result = defaultdict(F)
    for i,c in a.items():
        for j,d in b.items():
            result[i+j] += c*d
    return dict(result)


def correlation(sequence,h):
    return sum(sequence[n]*sequence[n+h] for n in range(len(sequence))
               if 0<=n+h<len(sequence))


def check_parseval():
    cases = 0
    for length,kind in product(range(1,10),range(4)):
        g = [F((n+kind)%5-2,n%3+1) for n in range(length)]
        u = [v+F((-1)**n*(n%4+kind),7) for n,v in enumerate(g)]
        poly_g = {n:v for n,v in enumerate(g)}
        poly_u = {n:v for n,v in enumerate(u)}
        dens_g = multiply(poly_g,{-n:v for n,v in poly_g.items()})
        dens_u = multiply(poly_u,{-n:v for n,v in poly_u.items()})
        difference = {n:dens_g.get(n,F(0))-dens_u.get(n,F(0))
                      for n in set(dens_g)|set(dens_u)}
        norm = sum((correlation(g,h)-correlation(u,h))**2 for h in range(1-length,length))
        assert multiply(difference,difference).get(0,F(0))==norm
        for h in range(1-length,length):
            assert difference[h]==correlation(g,h)-correlation(u,h)
        sup_bound = sum(abs(a-b) for a,b in zip(g,u))
        energy = sum(v*v for v in g+u)
        assert norm<=2*sup_bound**2*energy
        cases += 1
    return cases


def main():
    densities,intervals = check_comb()
    print('Exact Ramanujan/direct-period/CRT pair densities:',densities)
    print('Interval main-term and residue-count bounds:',intervals)
    print('Rankin Euler expansion/count cases:',check_rankin_factors())
    print('Independent Laurent/position autocorrelation Parseval cases:',check_parseval())
    print('Finite comb audits passed; asymptotic prime input remains explicitly cited.')


if __name__=='__main__':
    main()
