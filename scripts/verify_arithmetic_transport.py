#!/usr/bin/env python3
"""Exact finite identities supporting the quantitative arithmetic core proof.

This checks local factors, Euler convolution, overlap geometry and rational
integration. The asymptotic estimates are proved in notes/arithmetic_core_limit.md;
the surrounding zeta-moment reduction remains a separate obligation.
"""
from fractions import Fraction as F
from itertools import product
import math


def factors(n):
    result = {}
    p = 2
    while p*p<=n:
        while n%p==0:
            result[p] = result.get(p,0)+1
            n //= p
        p += 1
    if n>1:
        result[n] = 1
    return result


def beta_local(p,v,e1,e2):
    value = F(1)
    for e in (e1,e2):
        a = max(v-e,0)
        if a>=2:
            return F(0)
        if a==1:
            value *= F(-1,p-1)
    return value


def kappa_local(p,v,e1,e2):
    """Independent singular-series density formula, with kappa(-1)=0."""
    if v<0:
        return F(0)
    alpha = F(1,p-1)
    if e1==e2==0:
        return 1-alpha**2 if v==0 else 1+alpha
    if min(e1,e2)==0:
        return 1+alpha if v==0 else F(0)
    e = min(e1,e2)
    if e1==e2:
        if v<e:
            return F(0)
        return p**e*(1-alpha**2 if v==e else 1+alpha)
    return p**e*(1+alpha) if v==e else F(0)


def expected_local(p,v,e1,e2):
    e,alpha = min(e1,e2),F(1,p-1)
    vector = (1-alpha**2,alpha**2) if e1==e2 else (1+alpha,-alpha)
    return vector[v-e] if e<=v<=e+1 else F(0)


def h_coefficient(n):
    value = F(1)
    for p,e in factors(n).items():
        if p==2:
            if e!=1:
                return F(0)
            value *= F(-1,2)
        elif e==1:
            value *= F(2,p*(p-2))
        elif e==2:
            value *= F(-1,p*(p-2))
        else:
            return F(0)
    return value


def a_coefficient(n):
    value = F(1)
    for p,e in factors(n).items():
        if p==2 or e!=1:
            return F(0)
        value *= F(1,p-2)
    return value


def sigma(n):
    return math.prod((p**(e+1)-1)//(p-1) for p,e in factors(n).items())


def cutoff_weight(n):
    return F(mobius(n)**2*sigma(n),phi(n)**2)


def cutoff_h(n):
    value = F(1)
    for p,e in factors(n).items():
        if e==1:
            value *= F(3*p-1,p*(p-1)**2)
        elif e==2:
            value *= F(-p-1,p*(p-1)**2)
        else:
            return F(0)
    return value


def mobius(n):
    fs = factors(n)
    return 0 if any(e>1 for e in fs.values()) else (-1)**len(fs)


def phi(n):
    result = n
    for p in factors(n):
        result = result//p*(p-1)
    return result


def beta_q(q,b1,b2):
    q1,q2 = q//math.gcd(q,b1),q//math.gcd(q,b2)
    return F(mobius(q1)*mobius(q2),phi(q1)*phi(q2))


def cutoff_coefficients(cutoff,b1,b2):
    special = set(factors(2*b1*b2))
    result = {}
    signed_first = F(0)
    absolute_bound = 0
    for q in range(1,cutoff+1):
        if set(factors(q))<=special:
            continue
        beta = beta_q(q,b1,b2)
        assert abs(beta)<=1
        signed_first += beta*phi(q)
        for d in range(1,q+1):
            if q%d==0:
                result[d] = result.get(d,F(0))+mobius(q//d)*beta
                absolute_bound += d
    assert sum(d*c for d,c in result.items())==signed_first
    assert sum(d*abs(c) for d,c in result.items())<=absolute_bound
    return result


def overlap(*positions):
    return max(1-max(positions)+min(positions),F(0))


def original_overlap(x,y,nu):
    total = F(0)
    for a,c in ((x,nu-x),(nu-x,x)):
        for b,d in ((y,nu-y),(nu-y,y)):
            delta = a+c-b-d
            total += overlap(0,-a+delta,c-d,-d)
            total += overlap(0,b-c-d,-c-d,-d)
            total += overlap(0,-b-c+d,-c+d,d)
    return total


def reduced_overlap(x,y,nu):
    pair = 4*(1-nu+min(x,y))
    rest = sum(max(1-max(nu+c-b,b),F(0))
               for c in (x,nu-x) for b in (y,nu-y))
    return pair+2*rest


def multiply(a,b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[i+j] += x*y
    return result


def integrate(polynomial,left,right):
    return sum(a*F(right**(i+1)-left**(i+1),i+1)
               for i,a in enumerate(polynomial))


def main():
    primes = [p for p in range(2,102) if factors(p)=={p:1}]
    local_checks = 0
    for p,e1,e2,v in product(primes,range(5),range(5),range(7)):
        difference = beta_local(p,v,e1,e2)-beta_local(p,v+1,e1,e2)
        density = (kappa_local(p,v,e1,e2)-kappa_local(p,v-1,e1,e2))/p**v
        assert difference==density==expected_local(p,v,e1,e2)
        local_checks += 1
    # In the coprime case the weighted absolute special factors are bounded
    # independently of the exponents of the original prime-power moduli.
    for p,e in product(primes,range(1,6)):
        vector = [beta_local(p,v,e,0)-beta_local(p,v+1,e,0) for v in range(7)]
        assert vector==[1+F(1,p-1),F(-1,p-1)]+[F(0)]*5
        assert sum(p**v*abs(c) for v,c in enumerate(vector))<=4
    for p in primes[1:]:
        z = F(1,p)
        h1,h2 = F(2,p*(p-2)),F(-1,p*(p-2))
        assert 1+h1+h2==F((p-1)**2,p*(p-2))
        assert h1+2*h2==0  # odd-prime contribution to sum h(r) log(r)
        tail = z*(1+z)/(1-z)**2-z
        assert tail<=12*z*z
    # Exact Dirichlet convolution a=(1/id)*h through 2,000 coefficients.
    limit = 2000
    convolution = [F(0)]*(limit+1)
    for r in range(1,limit+1):
        h = h_coefficient(r)
        if h:
            for k in range(1,limit//r+1):
                convolution[r*k] += h/k
    assert all(convolution[n]==a_coefficient(n) for n in range(1,limit+1))
    # A second, independent positive weight controls growing q cutoffs.
    cutoff_convolution = [F(0)]*(limit+1)
    for r in range(1,limit+1):
        h = cutoff_h(r)
        if h:
            for k in range(1,limit//r+1):
                cutoff_convolution[r*k] += h/k
    assert all(cutoff_convolution[n]==cutoff_weight(n) for n in range(1,limit+1))
    for p in primes:
        assert abs(cutoff_h(p))+abs(cutoff_h(p*p))==F(4,(p-1)**2)
    # Compare direct q weights with the separated Euler-factor bound, without
    # using the divisor sawtooth coefficients to generate the left side.
    growing_cases = 0
    cases = ((1,1),(3,5),(81,625),(2**8,7**4),(3**7,5**6),
             (2,2),(2**4,2**7),(3,3),(3**4,3**7),(5**3,5**5))
    for cutoff,(b1,b2) in product((8,32,128),cases):
        ordinary = sum(cutoff_weight(q) for q in range(1,cutoff+1))
        direct = sum(abs(beta_q(q,b1,b2))*sigma(q) for q in range(1,cutoff+1))
        harmonic = sum(F(1,q) for q in range(1,cutoff+1))
        assert ordinary<=37*harmonic
        if math.gcd(b1,b2)==1:
            special_odd = set(factors(b1*b2))-{2}
            multiplier = math.prod(1+F(p+1,p-1) for p in special_odd)
            assert multiplier<=9 and direct<=multiplier*ordinary
        else:
            assert direct<=8*min(b1,b2)*ordinary
        growing_cases += 1
    for p,a,b in product(primes,range(1,8),range(1,8)):
        e = min(a,b)
        weighted = sum(abs(beta_local(p,v,a,b))*sigma(p**v) for v in range(e+2))
        assert weighted<=2*F(p,p-1)**2*p**e<=8*p**e
    geometry_checks = 0
    for i,j,k in product(range(13),range(7),range(7)):
        nu = 1+F(i,12)
        width = 1-nu/2
        x,y = nu-1+width*F(j,6),nu-1+width*F(k,6)
        original,reduced = original_overlap(x,y,nu),reduced_overlap(x,y,nu)
        assert original==reduced>=0
        geometry_checks += 1
    # Q(nu)=(2-nu)^3/6 + (3-2nu)_+^3/6, proved by cube and simplex integration.
    first,second = [F(1)],[F(1)]
    for _ in range(3):
        first = multiply(first,[F(2),F(-1)])
        second = multiply(second,[F(3),F(-2)])
    first,second = [a/6 for a in first],[a/6 for a in second]
    mass = integrate(first,F(1),F(2))+integrate(second,F(1),F(3,2))
    weighted = integrate(multiply(first,[F(-1),F(1)]),F(1),F(2))
    weighted += integrate(multiply(second,[F(-1),F(1)]),F(1),F(3,2))
    assert mass==F(1,16) and weighted==F(1,96) and -2*weighted==F(-1,48)
    ramanujan_checks = 0
    for q,j in product(range(1,41),range(1,41)):
        divisors = sum(d*mobius(q//d) for d in range(1,q+1) if q%d==0 and j%d==0)
        reduced = q//math.gcd(q,j)
        assert divisors==F(mobius(reduced)*phi(q),phi(reduced))
        ramanujan_checks += 1
    cutoff_cases = 0
    for cutoff,(b1,b2) in product((1,2,8,40),((1,1),(3,5),(9,25),(8,27),(5,5),(25,125))):
        coefficients = cutoff_coefficients(cutoff,b1,b2)
        # Check finite Ramanujan decomposition before using the sine transform.
        special = set(factors(2*b1*b2))
        for j in range(1,41):
            original = F(0)
            for q in range(1,cutoff+1):
                if not set(factors(q))<=special:
                    reduced = q//math.gcd(q,j)
                    original += beta_q(q,b1,b2)*F(mobius(reduced)*phi(q),phi(reduced))
            assert original==sum(d*c for d,c in coefficients.items() if j%d==0)
        cutoff_cases += 1
    print(f"Independent beta/density/shifted-vector checks: {local_checks}")
    print(f"Exact Euler-convolution coefficients: {limit}")
    print(f"Growing-cutoff convolution coefficients: {limit}; separated q-weight cases: {growing_cases}")
    print(f"Independent twelve-overlap/reduced-geometry checks: {geometry_checks}")
    print(f"Exact geometric mass {mass}; weighted mass {weighted}; model core {-2*weighted}")
    print(f"Exact Ramanujan identities: {ramanujan_checks}; independent cutoff decompositions: {cutoff_cases}")
    print("Finite identities passed. Asymptotic proofs and zeta-moment reduction are distinct.")


if __name__=="__main__":
    main()
