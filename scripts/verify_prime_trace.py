#!/usr/bin/env python3
"""Exact finite checks of the actual-prime trace decomposition.

The cyclic-frame checks calibrate algebra and Fourier factors. They do not
prove any asymptotic off-balance prime estimate or certify analytic proofs.
"""
import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
import json
from math import comb
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
    leakage = 0
    for d in range(1,12):
        for shift in range(-17,18):
            assert sum(not 0<=a+shift<d for a in range(d))==min(d,abs(shift))
            leakage += 1
    report = {
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'finite algebra checked; analytic proof draft; off-balance target unproved',
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
        'limitations':[
            'Finite cyclic frames calibrate Fourier algebra, not continuous asymptotics.',
            'The analytic Schatten and balanced-prime proofs require independent review.',
            'No actual off-balance prime estimate is established by these checks.',
            'Model off-balance constants are targets, not actual prime moment evaluations.',
            'No new unconditional zeta-zero counting bound follows.',
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {fourier} exact Fourier entries, six word-ledger degrees,')
    print(f'      {loops} periodic closed loops and {telescope} noncommuting telescopes.')
    print('Balanced trace limit:',balanced)
    print('Unproved off-balance cap for 80%:',off_cap)


if __name__=='__main__':
    main()
