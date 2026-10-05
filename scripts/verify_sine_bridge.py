#!/usr/bin/env python3
"""Finite exact checks for the global CUE/sine Gram comparison.

The checker audits Campbell partition normalization in a finite DPP,
exact Fourier gauge cancellation, and a Plancherel second-moment anchor.
It does not replace the ergodic or range-removal proof in the note.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial

from verify_cue_fluctuations import partitions


def determinant(matrix):
    a = [list(row) for row in matrix]
    value = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:
            return F(0)
        if pivot!=j:
            a[pivot],a[j] = a[j],a[pivot]
            value = -value
        d = a[j][j]
        value *= d
        for i in range(j+1,len(a)):
            ratio = a[i][j]/d
            for k in range(j+1,len(a)):
                a[i][k] -= ratio*a[j][k]
    return value


def trace_power(matrix,order):
    n = len(matrix)
    if not n:
        return F(0)
    power = [[F(i==j) for j in range(n)] for i in range(n)]
    for _ in range(order):
        power = [[sum(power[i][k]*matrix[k][j] for k in range(n))
                  for j in range(n)] for i in range(n)]
    return sum(power[i][i] for i in range(n))


def minor(matrix,indices):
    return [[matrix[i][j] for j in indices] for i in indices]


def check_discrete_campbell():
    size = 5
    sites = tuple(range(size))
    subsets = tuple(s for k in range(size+1) for s in combinations(sites,k))
    # A valid contraction by Gershgorin: diagonal 1/2, two neighbours 1/8.
    q = [[F(1,2) if i==j else F(1,8) if (i-j)%size in (1,size-1)
          else F(0) for j in sites] for i in sites]
    inclusion = {s:determinant(minor(q,s)) for s in subsets}
    probability = {}
    for s in subsets:
        missing = tuple(i for i in sites if i not in s)
        probability[s] = sum((-1)**len(t)*inclusion[tuple(sorted(s+t))]
                             for k in range(len(missing)+1)
                             for t in combinations(missing,k))
    assert all(p>=0 for p in probability.values())
    assert sum(probability.values())==1
    for s in subsets:
        assert sum(p for t,p in probability.items() if set(s)<=set(t))==inclusion[s]
    expected_count = sum(len(s)*p for s,p in probability.items())
    assert expected_count==sum(q[i][i] for i in sites)==F(5,2)
    cases = 0
    for cutoff in (1,2):
        kernel = [[F(3,4) if i==j else
                   F(1,2+min((i-j)%size,(j-i)%size))
                   if min((i-j)%size,(j-i)%size)<=cutoff else F(0)
                   for j in sites] for i in sites]
        for order in range(1,6):
            direct = sum(p*trace_power(minor(kernel,s),order)
                         for s,p in probability.items())
            expanded = F(0)
            for pi in partitions(order):
                block_of = {i:b for b,block in enumerate(pi) for i in block}
                for locations in product(sites,repeat=len(pi)):
                    if len(set(locations))!=len(locations):
                        continue
                    weight = inclusion[tuple(sorted(locations))]
                    for i in range(order):
                        weight *= kernel[locations[block_of[i]]][locations[block_of[(i+1)%order]]]
                    expanded += weight
            assert direct==expanded
            assert direct/expected_count==expanded/expected_count
            if order==2:
                ratio_expectation = sum(p*trace_power(minor(kernel,s),order)/len(s)
                                        for s,p in probability.items() if s)
                # Demonstrate why the random denominator cannot be pulled out.
                assert ratio_expectation!=direct/expected_count
            cases += 1
    return cases


ROOTS = ((F(1),F(0)),(F(0),F(1)),(F(-1),F(0)),(F(0),F(-1)))


def multiply(a,b):
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def fourier(count,n,difference):
    return tuple(sum(ROOTS[(a*difference)%4][j] for a in range(count))/n
                 for j in (0,1))


def check_gauge():
    cases = 0
    for n in (1,3,5,7):
        for m in range(1,n+1,2):
            # Locations x_j=n*j/4, so exponentials are exact fourth roots.
            centered = {}
            for d in range(-3,4):
                value = tuple(sum(ROOTS[((a-(m-1)//2)*d)%4][j]
                                  for a in range(m))/n for j in (0,1))
                assert value[1]==0
                centered[d] = value[0]
                assert fourier(m,n,d)==multiply(ROOTS[(((m-1)//2)*d)%4],value)
            for order in range(1,5):
                for walk in product(range(4),repeat=order):
                    complex_product,real_product = ROOTS[0],F(1)
                    for i in range(order):
                        d = walk[i]-walk[(i+1)%order]
                        complex_product = multiply(complex_product,fourier(m,n,d))
                        real_product *= centered[d]
                    assert complex_product==(real_product,F(0))
                    cases += 1
    return cases


def integrate_polynomial(coefficients,upper):
    return sum(c*upper**(j+1)/F(j+1) for j,c in enumerate(coefficients))


def check_second_moment():
    # Fourier transforms of S_lambda^2 and S_1^2 are triangular functions.
    # Their overlap integral is 2*int_0^lambda (lambda-u)(1-u) du.
    for lam in (F(i,20) for i in range(1,21)):
        overlap = 2*integrate_polynomial((lam,-lam-1,F(1)),lam)
        assert overlap==lam**2-lam**3/3
        column_moment = lam**2+lam-overlap
        assert column_moment==lam*(1+lam**2/3)
    # Nonnegative factorial coefficients for (1+N)^p; at Poisson mean one
    # these sum to Bell_(p+1), the growth bound used in the bridge proof.
    for p in range(1,9):
        coefficients = [sum(comb(p,k)*len([pi for pi in partitions(k)
                                           if len(pi)==j]) for k in range(j,p+1))
                        for j in range(p+1)]
        assert sum(coefficients)==len(partitions(p+1))
        for degree in range(10):
            assert sum(c*(factorial(degree)//factorial(degree-j) if j<=degree else 0)
                       for j,c in enumerate(coefficients))==(1+degree)**p


def main():
    print('Discrete DPP Campbell partition checks:',check_discrete_campbell())
    print('Exact closed-walk Fourier gauge checks:',check_gauge())
    check_second_moment()
    print('20 bandwidth second-moment anchors and Palm degree coefficients passed.')


if __name__=='__main__':
    main()
