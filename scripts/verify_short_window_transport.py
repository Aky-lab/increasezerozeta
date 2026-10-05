#!/usr/bin/env python3
"""Exact finite audits for weighted short-window prime-model transport.

These tests use rational abstract errors and phases 1 or (-1)^n. They
check discrete identities and norm bookkeeping, not the cited analytic
prime-uniformity theorem or an effective asymptotic height.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math


@lru_cache(None)
def progressions(length):
    output = {()}
    for start in range(length):
        for step in range(1,length+1):
            path = tuple(range(start,length,step))
            for end in range(1,len(path)+1):
                output.add(path[:end])
    return tuple(sorted(output))


def variation(values):
    return max(map(abs,values),default=F(0))+sum(abs(b-a) for a,b in zip(values,values[1:]))


def check_abel_and_cells():
    cases = 0
    for x,length in product((3,5,9),range(1,8)):
        for j in range(x,2*x):
            integers = tuple(range(j+1,j+length+1))
            for offset in (F(0),F(1,5),F(1,2),F(9,10)):
                left = F(j)+offset
                right = left+length
                selected = tuple(n for n in range(j,j+length+2) if left<n<=right)
                assert selected==integers
            errors = [F((n%7)-3,n%3+1) for n in integers]
            discrepancy = max((abs(sum(errors[i]*(-1)**(mode*integers[i]) for i in q))
                               for q in progressions(length) for mode in (0,1)),default=F(0))
            weights = ([F(1)]*length,
                       [F((-1)**n*(n%5+1),3) for n in integers],
                       [F(n,7) for n in integers])
            for q,mode,f in product(progressions(length),(0,1),weights):
                values = [errors[i]*(-1)**(mode*integers[i]) for i in q]
                direct = sum(f[i]*v for i,v in zip(q,values))
                prefixes = [sum(values[:k+1]) for k in range(len(values))]
                abel = (f[q[-1]]*prefixes[-1]-sum(prefixes[k]*(f[q[k+1]]-f[q[k]])
                         for k in range(len(q)-1))) if q else F(0)
                assert direct==abel
                assert variation([f[i] for i in q])<=variation(f)
                assert abs(direct)<=variation(f)*discrepancy
                cases += 1
    return cases


def check_coupled_products():
    cases = 0
    x = 9
    errors = [F((j%5)-2,7) for j in range(x)]
    d = [abs(v) for v in errors]
    model = [F(j%4-2,5) for j in range(x)]
    actual = [a+b for a,b in zip(model,errors)]
    amplitude = max(abs(v) for v in actual+model)
    for r,kind in product((2,3,4),range(4)):
        if kind==0:
            starts = [tuple([j]*r) for j in range(x)]  # Fully shared starts.
        elif kind==1:
            starts = [tuple((j+2*i)%x for i in range(r)) for j in range(x)]
        elif kind==2:
            starts = [tuple((j+k*i)%x for i in range(r)) for j,k in product(range(x),range(3))]
        else:
            starts = [tuple([0]*r)]*4+[(1,)*(r-1)+(j,) for j in range(x)]
        histograms = [[F(0)]*x for _ in range(r)]
        direct,telescope = F(0),F(0)
        for index,positions in enumerate(starts):
            coefficient = F((-1)**index*(index%4+1),index%3+1)
            v = [F(1+(index+i)%3,2) for i in range(r)]
            s = [v[i]*actual[j] for i,j in enumerate(positions)]
            t = [v[i]*model[j] for i,j in enumerate(positions)]
            difference = math.prod(s)-math.prod(t)
            expanded = sum((s[i]-t[i])*math.prod(s[:i])*math.prod(t[i+1:]) for i in range(r))
            assert difference==expanded
            direct += coefficient*difference
            telescope += coefficient*expanded
            mass = abs(coefficient)*math.prod(v)
            for i,j in enumerate(positions):
                histograms[i][j] += mass
        assert direct==telescope
        marginal_bound = amplitude**(r-1)*sum(w*d[j] for hist in histograms for j,w in enumerate(hist))
        assert abs(direct)<=marginal_bound
        norm_bound_squared = amplitude**(2*r-2)*r*sum(v*v for v in d)*sum(w*w for hist in histograms for w in hist)
        assert direct**2<=norm_bound_squared
        assert len({sum(hist) for hist in histograms})==1
        cases += 1
    return cases


def check_exceptional_and_concentration():
    for size in range(1,41):
        errors = [F(j%7-3,j%4+1) for j in range(size)]
        good = [j for j in range(size) if j%5]
        bad = [j for j in range(size) if not j%5]
        global_max = max(map(abs,errors))
        good_max = max((abs(errors[j]) for j in good),default=F(0))
        assert sum(e*e for e in errors)<=size*good_max**2+len(bad)*global_max**2
        for occupied in range(1,size+1):
            histogram = [F(1)]*occupied+[F(0)]*(size-occupied)
            total = sum(histogram)
            norm_squared = sum(w*w for w in histogram)
            effective = total**2/norm_squared
            assert effective==occupied
            assert F(size)*norm_squared/total**2==F(size,occupied)
        concentrated = [F(size)]+[F(0)]*(size-1)
        assert sum(concentrated)**2/sum(w*w for w in concentrated)==1
        uniform = [F(1)]*size
        assert sum(uniform)**2/sum(w*w for w in uniform)==size


def main():
    print('Exact progression/Abel checks:',check_abel_and_cells())
    print('Coupled product and marginal histogram cases:',check_coupled_products())
    check_exceptional_and_concentration()
    print('Integer-cell, good/bad and 40-size concentration identities passed.')
    print('Analytic inputs are the cited short-interval uniformity theorems.')


if __name__=='__main__':
    main()
