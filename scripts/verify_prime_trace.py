#!/usr/bin/env python3
"""Finite exact algebra and numerical calibrations of prime trace formulas.

The cyclic-frame checks calibrate algebra and Fourier factors. They do not
prove any asymptotic off-balance prime estimate or certify analytic proofs.
"""
import argparse
import cmath
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
import json
from math import comb, factorial, log, pi, sin, sqrt
from pathlib import Path
import random

from christoffel_exact import christoffel_at_zero
from model_moments import MODEL_MOMENTS
from verify_spectral_bridge import (add, eye, matrix_poly, multiply, scale,
                                    square, trace)


@dataclass(frozen=True, slots=True)
class G:
    """Gaussian rationals; all finite Fourier computations are exact."""
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other):
        if not isinstance(other,G):
            other = G(F(other))
        return G(self.re+other.re,self.im+other.im)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.re,-self.im)

    def __sub__(self, other):
        return self+-other if isinstance(other,G) else self+G(-F(other))

    def __mul__(self, other):
        if not isinstance(other,G):
            other = G(F(other))
        return G(self.re*other.re-self.im*other.im,
                 self.re*other.im+self.im*other.re)

    __rmul__ = __mul__

    def conjugate(self):
        return G(self.re,-self.im)


ZERO,ONE = G(),G(F(1))
ROOTS = (ONE,G(F(0),F(1)),-ONE,G(F(0),-F(1)))
PHI = (F(1),F(1,2),F(0),F(1,2))
N = 4


def phase(k):
    return ROOTS[k%N]


def convolution_frame(indices,s):
    return [[sum((PHI[u]*PHI[(s-u)%N]*phase(-a*u-b*(s-u))
                  for u in range(N)),ZERO) for b in indices] for a in indices]


def frames(indices):
    return {s:convolution_frame(indices,s) for s in range(-3,4)}


def zero_matrix(n):
    return [[ZERO for _ in range(n)] for _ in range(n)]


def identity(n):
    return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]


def finite_fourier_checks():
    indices = (0,1,2)
    frame = frames(indices)
    values = [[sum((PHI[u]*phase((t-a)*u) for u in range(N)),ZERO)
               for t in range(N)] for a in indices]
    assert all(v.im==0 for row in values for v in row)
    cases = 0
    for s in (-3,-2,-1,0,1,2,3):
        for a,b in product(range(len(indices)),repeat=2):
            direct = sum((values[a][t]*values[b][t]
                          *F(phase(t*s).re) for t in range(N)),ZERO)
            assert direct == (frame[s][a][b]+frame[-s][a][b])*F(N,2)
            assert frame[-s][a][b] == frame[s][b][a].conjugate()
            assert frame[s][a][b] == frame[s][b][a]
            cases += 1
    # Rational weights only calibrate the algebra, not actual Lambda values.
    atoms = ((2,1,(1,0),F(1,2)),(4,2,(2,0),F(1,4)),
             (3,3,(0,1),F(3,8)))
    h = zero_matrix(len(indices))
    word_atoms = []
    for n,s,exponents,weight in atoms:
        for sign in (-1,1):
            term = scale(frame[sign*s],weight/N)
            h = add(h,term)
            word_atoms.append((tuple(sign*x for x in exponents),term))
    alternate = [[sum((values[a][t]*values[b][t]
                     *sum(weight*phase(t*s).re for n,s,e,weight in atoms)
                     *F(2,N*N) for t in range(N)),ZERO)
                  for b in range(len(indices))] for a in range(len(indices))]
    assert h == alternate
    assert all(h[a][b].im==0 and h[a][b]==h[b][a]
               for a,b in product(range(len(indices)),repeat=2))

    ledger = []
    states = {(0,0):identity(len(indices))}
    power = identity(len(indices))
    for degree in range(1,7):
        new = {}
        for balance,word in states.items():
            for delta,atom in word_atoms:
                key = tuple(a+b for a,b in zip(balance,delta))
                contribution = multiply(word,atom)
                new[key] = add(new.get(key,zero_matrix(len(indices))),contribution)
        states = new
        power = multiply(power,h)
        total = sum((trace(m) for m in states.values()),ZERO)*F(1,len(indices))
        diagonal = trace(states.get((0,0),zero_matrix(len(indices))))*F(1,len(indices))
        off = sum((trace(m) for key,m in states.items() if key!=(0,0)),ZERO)*F(1,len(indices))
        assert total == trace(power)*F(1,len(indices)) == diagonal+off
        assert total.im==diagonal.im==off.im==0
        ledger.append({'degree':degree,'balance_cells':len(states),
                       'moment':str(total.re),'balanced':str(diagonal.re),
                       'off_balance':str(off.re)})
    return cases,ledger


def closed_loop_checks():
    frame = frames(tuple(range(N)))
    cases = 0
    for degree in range(2,7):
        selected = 0
        for shifts in product((-2,-1,1,2),repeat=degree):
            if sum(shifts):
                continue
            if degree>4 and selected>=64:
                break
            selected += 1
            word = identity(N)
            for s in shifts:
                word = multiply(word,scale(frame[s],F(1,N)))
            actual = trace(word)*F(1,N)
            prefixes,p = [0],0
            for s in shifts[:-1]:
                p += s
                prefixes.append(p)
            overlap = sum(product_weight(PHI[(u+p)%N]**2 for p in prefixes)
                          for u in range(N))/N
            assert actual == G(overlap)
            cases += 1
    return cases


def product_weight(values):
    out = F(1)
    for v in values:
        out *= v
    return out


def prime_word_classification():
    # Independent integer-product enumeration: exact balance, not congruence.
    atoms = ((2,2,1),(4,2,2),(3,3,1))
    signed = tuple((n,p,a,s) for n,p,a in atoms for s in (-1,1))
    counts = []
    for degree in range(1,7):
        balanced,leading = 0,0
        for word in product(signed,repeat=degree):
            plus = product_weight(n for n,p,a,s in word if s==1)
            minus = product_weight(n for n,p,a,s in word if s==-1)
            if plus!=minus:
                continue
            balanced += 1
            groups = {}
            for n,p,a,s in word:
                groups.setdefault(p,[]).append((a,s))
            for group in groups.values():
                assert len(group)>=2
                assert sum(a*s for a,s in group)==0
                assert {s for a,s in group}=={-1,1}
                if len(group)==2:
                    assert group[0][0]==group[1][0]
                else:
                    common = sum(a for a,s in group if s==1)
                    assert common >= (len(group)+1)//2
            if degree%2:
                assert any(len(group)>=3 for group in groups.values())
            if all(len(group)==2 and all(a==1 for a,s in group)
                   for group in groups.values()):
                leading += 1
                assert len(groups)==degree//2
        counts.append({'degree':degree,'balanced_words':balanced,
                       'distinct_ordinary_prime_pair_words':leading})
    assert counts[0]['balanced_words']==0
    return counts


def telescoping_checks(q):
    rng = random.Random(260510)
    for _ in range(24):
        b = [[F(rng.randrange(-4,5),3) for _ in range(3)] for _ in range(3)]
        r = [[F(rng.randrange(-2,3),9) for _ in range(3)] for _ in range(3)]
        for i,j in product(range(3),repeat=2):
            if i<j:
                b[j][i],r[j][i] = b[i][j],r[i][j]
        m = add(b,r)
        bp,mp = [eye(3)],[eye(3)]
        for _ in range(len(q)-1):
            bp.append(multiply(bp[-1],b))
            mp.append(multiply(mp[-1],m))
        telescoped = scale(eye(3),F(0))
        for k in range(1,len(q)):
            for h in range(k):
                telescoped = add(telescoped,scale(
                    multiply(mp[h],multiply(r,bp[k-1-h])),q[k]))
        assert add(matrix_poly(q,m),scale(matrix_poly(q,b),-1))==telescoped
    return 24


def double_leakage_checks():
    """Independent full-frame products versus compressed two-factor products."""
    frame = frames(tuple(range(N)))
    operators = {s:scale(a,F(1,N)) for s,a in frame.items()}
    cases = 0
    for d in range(1,N):
        for s,t in product(range(-3,4),repeat=2):
            a,b = operators[s],operators[t]
            sub_a = [row[:d] for row in a[:d]]
            sub_b = [row[:d] for row in b[:d]]
            compressed = trace(multiply(sub_a,sub_b))
            full = sum(multiply(a,b)[i][i] for i in range(d))
            crossing = sum((a[i][k]*b[k][i]
                            for i in range(d) for k in range(d,N)),ZERO)
            assert full-compressed == crossing
            left = sum(a[i][k].re**2+a[i][k].im**2
                       for i in range(d) for k in range(d,N))
            right = sum(b[k][i].re**2+b[k][i].im**2
                        for i in range(d) for k in range(d,N))
            assert crossing.re**2+crossing.im**2 <= left*right
            amplitude = sum(PHI[u]*PHI[(u+s+t)%N]*PHI[(u+s)%N]**2
                            for u in range(N))/N
            geometric = sum((phase(-i*(s+t)) for i in range(d)),ZERO)*F(1,d)
            assert full*F(1,d) == amplitude*geometric
            cases += 1
    return cases


def logarithmic_mean_value_calibrations():
    """Numerical finite calibrations; not a proof of the cited inequality."""
    prime_powers = ((2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),
                    (11,11),(13,13),(16,2),(17,17),(19,19),(23,23),
                    (25,5),(27,3),(29,29),(31,31),(32,2),(37,37))
    length = log(40)
    frequency = [log(n) for n,p in prime_powers]
    for i,(n,p) in enumerate(prime_powers):
        circular = min(min(abs(frequency[i]-s),length-abs(frequency[i]-s))
                       for j,s in enumerate(frequency) if i!=j)
        assert circular+1e-14 >= 1/(2*n)
    mean_cases,geometric_cases = 0,0
    for start,d in product((0.0,10.0,123.5),(1,7,31,80)):
        kernel = []
        for i,s in enumerate(frequency):
            row = []
            for j,t in enumerate(frequency):
                delta = s-t
                direct = sum(cmath.exp(-1j*(start+2*pi*a/length)*delta)
                             for a in range(d))/d
                if i!=j:
                    endpoint = (cmath.exp(-1j*(start-pi/length)*delta)
                                -cmath.exp(-1j*(start+(2*d-1)*pi/length)*delta))
                    endpoint /= 2j*d*sin(pi*delta/length)
                    assert abs(direct-endpoint)<1e-10
                    geometric_cases += 1
                row.append(direct)
            kernel.append(row)
        for v in (length/8,length/3,length/2):
            coefficients = [log(p)/sqrt(n)*max(0.0,1-4*(v-s)**2/length**2)
                            *cmath.exp(0.31j*n)
                            for (n,p),s in zip(prime_powers,frequency)]
            actual = sum(abs(sum(x*cmath.exp(-1j*(start+2*pi*a/length)*s)
                                 for x,s in zip(coefficients,frequency)))**2
                         for a in range(d))/d
            diagonal = sum(abs(x)**2 for x in coefficients)
            off = sum(coefficients[i]*coefficients[j].conjugate()*kernel[i][j]
                      for i in range(len(coefficients))
                      for j in range(len(coefficients)) if i!=j)
            assert abs(actual-diagonal-off)<1e-9
            energy = sum(n*abs(x)**2 for (n,p),x in zip(prime_powers,coefficients))
            assert abs(off)<=3*length*energy/d+1e-9
            mean_cases += 1
    alias_cases = 0
    for length in (3.0,5.0,10.0,100.0):
        for fraction in (0.0,0.25,0.5,0.9,0.99,0.99999):
            s = log(4)+(length-log(4))*fraction
            ratio = (1-s/length)/sin(pi*s/length)
            assert ratio<=length/(2*log(4))+1e-9
            alias_cases += 1
    exponents = []
    for degree in range(2,9):
        for bandwidth in (F(1),F(2,degree),F(1,4)):
            exponent = degree*bandwidth/2-1
            assert (exponent<=0)==(degree*bandwidth<=2)
            exponents.append({'degree':degree,'bandwidth':str(bandwidth),
                              'power_exponent':str(exponent)})
    return mean_cases,geometric_cases,alias_cases,exponents


def fixed_taper_geometry():
    """Exact taper calculus, smoothing moments and nonperiodic gauge checks."""
    ramp = [F(0),F(0),F(0),F(10),F(-15),F(6)]
    def derivative(p):
        return [j*p[j] for j in range(1,len(p))]
    def value(p,x):
        return sum(c*x**j for j,c in enumerate(p))
    def integral(p,lo,hi):
        return sum(c*(hi**(j+1)-lo**(j+1))/F(j+1) for j,c in enumerate(p))
    density = derivative(ramp)
    second = derivative(density)
    assert density == [F(0),F(0),F(30),F(-60),F(30)]
    assert second == [F(0),F(60),F(-180),F(120)]
    assert value(ramp,F(0))==0 and value(ramp,F(1))==1
    assert all(value(p,x)==0 for p in (density,second) for x in (F(0),F(1)))
    assert integral(ramp,F(0),F(1))==F(1,2)
    assert integral(density,F(0),F(1))==1
    assert value(density,F(1,2))==F(15,8)
    assert integral(second,F(0),F(1,2))-integral(second,F(1,2),F(1))==F(15,4)
    # Two disjoint transition edges of phi for L>2.
    first_l1,first_linf,second_l1 = F(2),F(15,8),F(15,2)
    product_second = 2*second_l1+2*first_linf*first_l1
    assert product_second==F(45,2)
    # Independently integrate centered density versus its closed moment law.
    centered_density = [F(15,8),F(0),F(-15),F(0),F(30)]
    moments = []
    for j in range(13):
        moment = integral([F(0)]*(2*j)+centered_density,F(-1,2),F(1,2))
        expected = F(15,2**(2*j)*(2*j+1)*(2*j+3)*(2*j+5))
        assert moment==expected
        moments.append(str(moment))
        # Series coefficient of 120*((12-w^2)*sin(w/2)-6*w*cos(w/2))/w^5.
        n = 2*j+5
        sine = lambda k: F((-1)**((k-1)//2),2**k*factorial(k))
        cosine = lambda k: F((-1)**(k//2),2**k*factorial(k))
        numerator = 120*(12*sine(n)-sine(n-2)-6*cosine(n-1))
        assert numerator==(-1)**j*moment/F(factorial(2*j))
    for n in (1,3):
        sine = lambda k: F((-1)**((k-1)//2),2**k*factorial(k))
        cosine = lambda k: F((-1)**(k//2),2**k*factorial(k))
        assert 12*sine(n)-(sine(n-2) if n>=3 else 0)-6*cosine(n-1)==0
    # Rational envelope: a/k and b/k^2 bounds imply a finite harmonic
    # contribution plus a uniformly bounded cubic-series tail.
    envelope_cases = 0
    for a,b,d in product((F(1,3),F(2),F(7,3)),(F(1,2),F(3),F(17)),(1,7,31)):
        ratio = b/a
        cutoff = max(1,(ratio.numerator+ratio.denominator-1)//ratio.denominator)
        harmonic = sum(F(1,k) for k in range(1,cutoff+1))
        bound = 2*a*a*harmonic+b*b/F(cutoff**2)
        direct = 2*sum(min(d,k)*min(a*a/F(k*k),b*b/F(k**4))
                       for k in range(1,257))
        assert direct<=bound and b*b/F(cutoff**2)<=a*a
        envelope_cases += 1
    for k in range(2,130):
        assert F(1,2)*(F(1,(k-1)**2)-F(1,k*k))-F(1,k**3)==F(3*k-2,2*k**3*(k-1)**2)
    # A noncommensurate rational unit phase, with zero-extended gates.
    unit = G(F(3,5),F(4,5))
    def unit_power(k):
        base = unit if k>=0 else unit.conjugate()
        out = ONE
        for _ in range(abs(k)):
            out = out*base
        return out
    taper = [F(0),F(1,2),F(1),F(1,2),F(0)]
    gauge_cases = 0
    for s in range(-4,5):
        operator = [[G(taper[u]*taper[v]) if v==u+s else ZERO
                     for v in range(5)] for u in range(5)]
        for u,v in product(range(5),repeat=2):
            transformed = unit_power(u)*operator[u][v]*unit_power(-v)
            assert transformed==unit_power(-s)*operator[u][v]
        gauge_cases += 1
    # Allowing surviving wraps would violate that identity at this phase.
    assert unit_power(3)*G(taper[3]*taper[1])*unit_power(-1)!=unit_power(-3)*G(taper[3]*taper[1])
    # Direct midpoint integration of the two shifted tapers is independent
    # of the centered-density/convolution derivation of the closed kernel.
    kernel_cases,max_kernel_error = 0,0.0
    def ramp_value(x):
        if x<=0:
            return 0.0
        if x>=1:
            return 1.0
        return 10*x**3-15*x**4+6*x**5
    for length in (4.0,6.0,10.0):
        def phi_value(x):
            return ramp_value(length/2+x)*ramp_value(length/2-x)
        grid = [-length/2+length*(j+0.5)/8192 for j in range(8192)]
        for s in (1.0,length-2):
            overlap = [phi_value(v+s/2)*phi_value(v-s/2) for v in grid]
            for k in (0,1,3,9):
                w = 2*pi*k/length
                direct = sum(h*cmath.cos(w*v).real for h,v in zip(overlap,grid))/8192
                if k==0:
                    closed = 1-(s+1)/length
                else:
                    if abs(w)<1:
                        psi = sum((-1)**j*w**(2*j)/factorial(2*j)
                                  *float(F(15,2**(2*j)*(2*j+1)*(2*j+3)*(2*j+5)))
                                  for j in range(11))
                    else:
                        psi = 120*((12-w*w)*sin(w/2)-6*w*cmath.cos(w/2).real)/w**5
                    closed = psi*sin(pi*k*(1-(s+1)/length))/(pi*k)
                error = abs(direct-closed)
                assert error<1e-8
                max_kernel_error = max(max_kernel_error,error)
                kernel_cases += 1
    return dict(ramp_coefficients=list(map(str,ramp)),
                ramp_integral='1/2',phi_first_derivative_l1=str(first_l1),
                phi_first_derivative_linf=str(first_linf),
                phi_second_derivative_l1=str(second_l1),
                overlap_second_derivative_l1_upper_bound=str(product_second),
                centered_smoothing_even_moments=moments,
                exact_smoothing_series_coefficients=13,
                rational_two_decay_envelope_cases=envelope_cases,
                exact_noncommensurate_gauge_cases=gauge_cases,
                wrapped_gate_counterexample=True,
                numerical_direct_overlap_kernel_cases=kernel_cases,
                numerical_direct_overlap_max_error=max_kernel_error,
                continuous_leakage_proof='notes/fixed_taper_projection.md')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    q = [F(1),-F(8232,2519),F(7368,2519),-F(1932,2519)]
    p = [sum(q[k]*comb(k,j)*(-1)**j for k in range(j,4)) for j in range(4)]
    assert p == [F(c,2519) for c in (-277,-708,1572,1932)]
    coeff = square(p)
    numerator = [76729,392232,-369624,-3296280,-264528,6074208,3732624]
    assert coeff == [F(c,6345361) for c in numerator]
    centered = [sum((-1)**k*comb(j,k)*MODEL_MOMENTS[k]
                    for k in range(j+1)) for j in range(7)]
    assert centered == [F(1),F(0),F(1,3),F(0),F(1,4),-F(1,36),F(61,252)]
    diagonal = [F(1),F(0),F(1,3),F(0),F(4,15),F(0),F(32,105)]
    off = [h-d for h,d in zip(centered,diagonal)]
    assert off == [F(0),F(0),F(0),F(0),-F(1,60),-F(1,36),-F(79,1260)]
    model = sum(c*m for c,m in zip(coeff,centered))
    balanced = sum(c*m for c,m in zip(coeff,diagonal))
    off_model = sum(c*m for c,m in zip(coeff,off))
    assert model == F(247,2519)
    assert balanced == F(5102709,31726805)
    assert off_model == -F(1991744,31726805)
    off_cap = F(1,10)-balanced
    assert off_cap == -F(3860057,63453610)
    assert off_cap-off_model == F(49,25190)
    assert balanced > F(1,10)

    lam4,q4,_ = christoffel_at_zero(MODEL_MOMENTS,4)
    p4 = [sum(q4[k]*comb(k,j)*(-1)**j for k in range(j,5)) for j in range(5)]
    coeff4 = square(p4)
    diagonal4 = diagonal+[F(0),F(1661,3780)]
    balanced4 = sum(c*m for c,m in zip(coeff4,diagonal4))
    off_cap4 = F(1,10)-balanced4
    off_model4 = lam4-balanced4
    assert balanced4 == F(38601698857019973,132097166600162405)
    assert off_cap4 == -F(10156792878801493,52838866640064962)
    assert off_model4 == -F(28653310482603548,132097166600162405)
    assert off_cap4-off_model4 == F(40129409,1625405590)

    fourier,ledger = finite_fourier_checks()
    loops = closed_loop_checks()
    classifications = prime_word_classification()
    telescope = telescoping_checks(q)
    double_leakage = double_leakage_checks()
    taper_geometry = fixed_taper_geometry()
    means,geometric,aliases,envelopes = logarithmic_mean_value_calibrations()
    leakage = 0
    for d in range(1,12):
        for shift in range(-17,18):
            assert sum(not 0<=a+shift<d for a in range(d))==min(d,abs(shift))
            leakage += 1
    report = {
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'finite algebra and numerical calibrations checked; analytic proof draft; higher off-balance target unproved',
        'p3_coefficients':[str(x) for x in p],
        'squared_prime_polynomial_coefficients':[str(x) for x in coeff],
        'unit_bandwidth_balanced_limits':[str(x) for x in diagonal],
        'centered_model_moments':[str(x) for x in centered],
        'model_off_balance_targets':[str(x) for x in off],
        'balanced_polynomial_trace_limit':str(balanced),
        'model_off_balance_polynomial_trace':str(off_model),
        'off_balance_cap_for_80_percent':str(off_cap),
        'excess_margin':str(off_cap-off_model),
        'q4_balanced_polynomial_trace_limit':str(balanced4),
        'q4_model_off_balance_polynomial_trace':str(off_model4),
        'q4_off_balance_cap_for_80_percent':str(off_cap4),
        'q4_excess_margin':str(off_cap4-off_model4),
        'exact_cyclic_fourier_entry_cases':fourier,
        'exact_finite_frame_word_ledger':ledger,
        'exact_periodic_closed_loop_cases':loops,
        'integer_product_classifications':classifications,
        'noncommuting_telescoping_cases':telescope,
        'projection_leakage_count_cases':leakage,
        'exact_double_projection_and_open_word_cases':double_leakage,
        'fixed_taper_geometry':taper_geometry,
        'numerical_logarithmic_sampled_mean_value_cases':means,
        'numerical_geometric_endpoint_cases':geometric,
        'numerical_same_sign_alias_cases':aliases,
        'absolute_frame_envelope_exponents':envelopes,
        'analytic_second_off_balance_status':'published-input deduction in proof draft; O2=o(1)',
        'remaining_cubic_off_balance_orders':[3,4,5,6],
        'limitations':[
            'Finite cyclic frames calibrate Fourier algebra, not continuous asymptotics.',
            'The analytic Schatten and balanced-prime proofs require independent review.',
            'No actual off-balance prime estimate is established by these checks.',
            'Higher model off-balance constants are targets, not actual prime moment evaluations.',
            'No new unconditional zeta-zero counting bound follows.',
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {fourier} exact Fourier entries, six word-ledger degrees,')
    print(f'      {loops} periodic closed loops and {telescope} noncommuting telescopes.')
    print(f'      {double_leakage} exact double-leakage cases, {means} sampled mean-value calibrations.')
    print('Balanced trace limit:',balanced)
    print('Unproved off-balance cap for 80%:',off_cap)


if __name__=='__main__':
    main()
