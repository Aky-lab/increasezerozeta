#!/usr/bin/env python3
"""Exact residue and cube audits for qualitative rectangle transport."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
import math


VERTICES=tuple(product(range(2),repeat=3))
ROWS=tuple((1,)+v for v in VERTICES)


def rank(rows,prime=None):
    a=[list(map(Q,row)) if prime is None else [x%prime for x in row] for row in rows]
    r=0
    for c in range(4):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:
            continue
        a[r],a[pivot]=a[pivot],a[r]
        inverse=1/a[r][c] if prime is None else pow(a[r][c],-1,prime)
        a[r]=[v*inverse if prime is None else v*inverse%prime for v in a[r]]
        for i in range(len(a)):
            if i==r:
                continue
            factor=a[i][c]
            a[i]=[v-factor*w if prime is None else (v-factor*w)%prime for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):
            break
    return r


def determinant(a):
    n=len(a)
    total=0
    for permutation in permutations(range(n)):
        sign=(-1)**sum(permutation[i]>permutation[j] for i in range(n) for j in range(i+1,n))
        total+=sign*math.prod(a[i][permutation[i]] for i in range(n))
    return total


def admissible_count(modulus):
    return sum(all(math.gcd(x+aa*a+bb*b+cc*c,modulus)==1 for aa,bb,cc in VERTICES)
               for x,a,b,c in product(range(modulus),repeat=4))


def bprime(p):
    return 1 if p==2 else p**4-8*p**3+28*p*p-44*p+23


def check_rank_certificate():
    odd_coefficients=[0]*5
    even_coefficients=[0]*5
    for mask in range(256):
        rows=[row for i,row in enumerate(ROWS) if mask&(1<<i)]
        rq=rank(rows)
        assert all(rank(rows,p)==rq for p in (3,5,7))
        odd_coefficients[4-rq]+=(-1)**len(rows)
        even_coefficients[4-rank(rows,2)]+=(-1)**len(rows)
    assert odd_coefficients==[23,-44,28,-8,1]
    assert even_coefficients==[21,-42,28,-8,1]
    minor_count=0
    for size in range(1,5):
        for row_indices,column_indices in product(combinations(range(8),size),combinations(range(4),size)):
            minor=[[ROWS[i][j] for j in column_indices] for i in row_indices]
            assert abs(determinant(minor))<=2
            minor_count+=1
    for p in (2,3,5,7,11):
        assert admissible_count(p)==bprime(p)
    assert Q(bprime(2)*2**4,(2-1)**8)==16
    assert Q(bprime(3)*3**4,(3-1)**8)==Q(81,32)
    assert Q(4,5)**8>Q(4,25)
    for p in (5,7,11,13,17,19,23,31,101):
        numerator=bprime(p)*p**4-(p-1)**8
        assert numerator==12*p**5-47*p**4+56*p**3-28*p*p+8*p-1
        assert numerator<=15*p**5
        assert Q(bprime(p)*p**4,(p-1)**8)<=1+Q(100,p**3)
    return 256,minor_count


def prime_divisors(n):
    return [p for p in range(2,n+1) if n%p==0 and all(p%d for d in range(2,math.isqrt(p)+1))]


def phi(n):
    return sum(math.gcd(i,n)==1 for i in range(n))


def check_crt_budget():
    cases=0
    for period in (2,6,10,30):
        c=Q(period,phi(period))
        count=admissible_count(period)
        assert count==math.prod(bprime(p) for p in prime_divisors(period))
        assert c**8*Q(count,period**4)==math.prod(Q(bprime(p)*p**4,(p-1)**8)
                                                for p in prime_divisors(period))
        cases+=1
    return cases


def check_residue_windows():
    cases=0
    for j,y,w in product(range(11),range(3,10),(2,3,6)):
        start=j//w-1
        length=(y+w-1)//w+6
        common=set(range(start+1,start+length+1))
        for b in range(1,w+1):
            actual={n for n in common if j<w*n+b<=j+y}
            direct={n for n in range(start-2,start+length+3) if j<w*n+b<=j+y}
            assert actual==direct
            assert len(common-actual)<=8
            if actual:
                assert len(actual)==max(actual)-min(actual)+1
            cases+=1
    return cases


def cube_sum(functions,window):
    left,right=window
    total=Q(0)
    for x,a,b,c in product(range(left,right+1),range(left-right,right-left+1),
                           range(left-right,right-left+1),range(left-right,right-left+1)):
        positions=[x+aa*a+bb*b+cc*c for aa,bb,cc in VERTICES]
        total+=math.prod(f.get(n,Q(0)) for f,n in zip(functions,positions))
    return total


def check_cube_splitting():
    cases=0
    for j,y,w in ((2,4,2),(5,5,3),(7,5,6)):
        c=Q(w,phi(w))
        prime={n:Q(n%5+1) if math.gcd(n,w)==1 else Q(0) for n in range(j+1,j+y+1)}
        model={n:c if math.gcd(n,w)==1 else Q(0) for n in prime}
        difference={n:prime[n]-model[n] for n in prime}
        functions=[difference if vertex[1:]==(0,0) else (prime if vertex[0]==0 else model)
                   for vertex in VERTICES]
        direct=cube_sum(functions,(j+1,j+y))
        split=defaultdict(Q)
        for x,a,b,cc in product(range(j+1,j+y+1),range(-y,y+1),range(-y,y+1),range(-y,y+1)):
            variables=(x,a,b,cc)
            residues=tuple(v%w for v in variables)
            quotients=tuple((v-r)//w for v,r in zip(variables,residues))
            values=[]
            all_units=True
            for vertex,function in zip(VERTICES,functions):
                weights=(1,)+vertex
                residue_sum=sum(q*r for q,r in zip(weights,residues))
                residue=residue_sum%w
                carry=(residue_sum-residue)//w
                quotient=sum(q*z for q,z in zip(weights,quotients))+carry
                position=sum(q*v for q,v in zip(weights,variables))
                assert position==w*quotient+residue
                assert 0<=carry<=3
                all_units &= math.gcd(residue,w)==1
                values.append(function.get(position,Q(0))/c)
            contribution=math.prod(values)
            if not all_units:
                assert contribution==0
            split[residues]+=contribution
        assert direct==c**8*sum(split.values())
        cases+=1
    return cases


def model_pair(period,k,h):
    c=Q(period,phi(period))
    offsets=(0,-k,h,h-k)
    direct=sum(c**4 for n in range(period) if all(math.gcd(n+t,period)==1 for t in offsets))/period
    local=math.prod((1-Q(len({t%p for t in offsets}),p))/(1-Q(1,p))**4
                    for p in prime_divisors(period))
    assert direct==local
    return direct


def check_four_point_main():
    cases=0
    for period in (6,30):
        c=Q(period,phi(period))
        for k,h in product(range(-4,5),repeat=2):
            density=model_pair(period,k,h)
            for j,y in ((0,5),(3,9),(31,15)):
                offsets=(0,-k,h,h-k)
                positions=[n for n in range(j+1-max(offsets),j+y+1-min(offsets))
                           if all(j<n+t<=j+y for t in offsets)]
                size=max(y-abs(k)-abs(h),0)
                assert len(positions)==size
                actual=sum(c**4 for n in positions if all(math.gcd(n+t,period)==1 for t in offsets))
                assert abs(actual-size*density)<=c**4*period
                cases+=1
    return cases


def main():
    print('Cube subset-rank certificates and square minors:',check_rank_certificate())
    print('Direct CRT and normalization-budget cases:',check_crt_budget())
    print('Residue-window boundary cases:',check_residue_windows())
    print('Exact quotient/carry cube decompositions:',check_cube_splitting())
    print('Four-point period and finite-window main terms:',check_four_point_main())
    print('Finite cube audits passed; qualitative prime input is the cited theorem.')


if __name__=='__main__':
    main()
