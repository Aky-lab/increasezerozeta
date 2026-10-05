#!/usr/bin/env python3
"""Exact bounded-certificate algebra; no arithmetic resolvent estimate."""
import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
import json

from verify_prime_trace import G
from verify_spectral_bridge import add, eye, multiply, scale, square, trace

A = F(1932,2519)
B = 1/A  # r^3=B, r is the positive real root.
Q = [F(1),-F(8232,2519),F(7368,2519),-A]
N = square(Q)


@dataclass(frozen=True)
class E:
    """Exact Q[r,i] coordinates, with r^3=B and i^2=-1."""
    v: tuple = (F(0),)*6

    @staticmethod
    def lift(x):
        return x if isinstance(x,E) else E((F(x),)+(F(0),)*5)

    def __add__(self,x):
        x=E.lift(x)
        return E(tuple(a+b for a,b in zip(self.v,x.v)))

    __radd__=__add__

    def __neg__(self):
        return E(tuple(-x for x in self.v))

    def __sub__(self,x):
        return self+-E.lift(x)

    def __rsub__(self,x):
        return E.lift(x)+-self

    def __mul__(self,x):
        x=E.lift(x)
        out=[F(0)]*6
        for j,a in enumerate(self.v):
            for k,b in enumerate(x.v):
                degree=j%3+k%3
                imag=j//3+k//3
                out[degree%3+3*(imag%2)] += a*b*B**(degree//3)*(-1)**(imag//2)
        return E(tuple(out))

    __rmul__=__mul__

    def __pow__(self,n):
        assert n>=0
        out=E.lift(1)
        for _ in range(n):
            out=out*self
        return out

    def conjugate(self):
        return E(self.v[:3]+tuple(-x for x in self.v[3:]))

    def record(self):
        return [str(x) for x in self.v]


R=E((F(0),F(1),F(0),F(0),F(0),F(0)))
I=E((F(0),F(0),F(0),F(1),F(0),F(0)))
Z=I*R
RINV=R**2*F(1,B)
ZINV=-I*RINV


def pmul(p,q):
    out=[E.lift(0)]*(len(p)+len(q)-1)
    for j,x in enumerate(p):
        for k,y in enumerate(q):
            out[j+k]=out[j+k]+x*y
    return out


def ppow(p,n):
    out=[E.lift(1)]
    for _ in range(n):
        out=pmul(out,p)
    return out


def peval(p,x):
    return sum((E.lift(c)*x**j for j,c in enumerate(p)),E.lift(0))


def derivative(p):
    return [j*p[j] for j in range(1,len(p))]


def partial_fractions():
    assert R**3==E.lift(B) and Z*ZINV==E.lift(1)
    u=ZINV*F(1,2)
    nz,n1,n2=(peval(p,Z) for p in (N,derivative(N),derivative(derivative(N))))
    coeff=[R**6*(n2*u**3*F(1,2)-3*n1*u**4+6*nz*u**5),
           R**6*(n1*u**3-3*nz*u**4),R**6*nz*u**3]
    minus,plus=[-Z,E.lift(1)],[Z,E.lift(1)]
    denominator=ppow([R**2,E.lift(0),E.lift(1)],3)
    lhs=[R**6*c-denominator[j] for j,c in enumerate(N)]
    rhs=[E.lift(0)]*7
    for k,c in enumerate(coeff,1):
        p=pmul(ppow(minus,3-k),ppow(plus,3))
        pc=pmul(ppow(plus,3-k),ppow(minus,3))
        for j in range(len(p)):
            rhs[j]=rhs[j]+c*p[j]+c.conjugate()*pc[j]
    assert lhs==rhs
    tinv=(1+Z)*(1-R**2+R**4)*F(1,1+B*B)
    assert (1-Z)*tinv==E.lift(1)
    closing=[coeff[0]*tinv**2+2*coeff[1]*tinv**3+3*coeff[2]*tinv**4,
             coeff[1]*tinv**2+2*coeff[2]*tinv**3,coeff[2]*tinv**2]
    dinv=(1-R**2+R**4)*F(1,1+B*B)
    f1=R**6*sum(Q)**2*dinv**3
    value=sum((c*tinv**k for k,c in enumerate(coeff,1)),E.lift(0))
    assert 1+value+value.conjugate()==f1
    return coeff,closing,f1,len(lhs)


def ginv(x):
    if not isinstance(x,G):
        x=G(F(x))
    norm=x.re*x.re+x.im*x.im
    assert norm
    return G(x.re/norm,-x.im/norm)


def inverse(m):
    n=len(m)
    rows=[[x if isinstance(x,G) else G(F(x)) for x in row]
          +[G(F(j==k)) for k in range(n)] for j,row in enumerate(m)]
    for j in range(n):
        k=next(k for k in range(j,n) if rows[k][j]!=G())
        rows[j],rows[k]=rows[k],rows[j]
        rows[j]=[ginv(rows[j][j])*x for x in rows[j]]
        for k in range(n):
            if k!=j:
                factor=rows[k][j]
                rows[k]=[x-factor*y for x,y in zip(rows[k],rows[j])]
    assert [row[:n] for row in rows]==[[G(F(j==k)) for k in range(n)] for j in range(n)]
    return [row[n:] for row in rows]


def power(m,n):
    out=eye(len(m))
    for _ in range(n):
        out=multiply(out,m)
    return out


def block(m,rows,cols):
    return [[m[j][k] for k in cols] for j in rows]


def transpose(m):
    return [list(row) for row in zip(*m)]


def norm_squared(z):
    return z.re*z.re+z.im*z.im


def resolvent_checks():
    # Large omitted diagonal entries test that no ||J|| factor is inserted.
    cases,insertion=0,0
    for outside,z in product((F(-3),F(5),F(10**6)),(G(F(0),F(1)),G(F(1),F(2)))):
        j=[[F(-1),F(1,3),F(2,5),F(-1,7)],
           [F(1,3),F(2),F(3,7),F(1,4)],
           [F(2,5),F(3,7),outside,F(1,2)],
           [F(-1,7),F(1,4),F(1,2),-outside]]
        p,q=range(2),range(2,4)
        am,dm,v=block(j,p,p),block(j,q,q),block(j,q,p)
        ra,rd,rj=(inverse(add(m,scale(eye(len(m)),-z))) for m in (am,dm,j))
        vh=transpose(v)
        hs2=sum(x*x for row in v for x in row)
        km=add(eye(len(j)),scale(j,-1))
        first=trace(block(km,p,p))*F(1,2)
        ti=ginv(G(F(1))-z)
        xis=[trace(block(multiply(multiply(km,power(rj,k)),km),p,p))*F(1,2)
             for k in (1,2,3)]
        formula=[ti+ti*ti*(first+xis[0]),
                 ti*ti+2*ti*ti*ti*(first+xis[0])+ti*ti*xis[1],
                 ti*ti*ti+3*ti*ti*ti*ti*(first+xis[0])+2*ti*ti*ti*xis[1]+ti*ti*xis[2]]
        for k,value in enumerate(formula,1):
            assert value==trace(block(power(rj,k),p,p))*F(1,2)
            insertion+=1
        for k in (1,2,3):
            lhs=add(block(power(rj,k),p,p),scale(power(ra,k),-1))
            rhs=[[G(),G()],[G(),G()]]
            terms=0
            for b in range(k):
                for c in range(k-b):
                    e=k-1-b-c
                    term=multiply(multiply(multiply(multiply(power(ra,b+1),vh),power(rd,c+1)),v),block(power(rj,e+1),p,p))
                    rhs=add(rhs,term)
                    terms+=1
            assert terms==k*(k+1)//2 and lhs==rhs
            bound=terms*hs2/z.im**(k+2)
            assert norm_squared(trace(lhs))<=bound*bound
            cases+=1
    return cases,insertion


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=str)
    args=parser.parse_args()
    expected=[6345361,-41472816,104885808,-131040168,86095872,-28469952,3732624]
    assert N==[F(x,6345361) for x in expected]
    assert 0<A<1 and N[2]>3 and N[4]>3 and N[6]==A*A
    assert all((-1)**j*N[j]>0 for j in range(7))
    assert sum(Q)**2==F(76729,6345361)<F(1,10)
    coeff,closing,f1,identities=partial_fractions()
    blocks,insertion=resolvent_checks()
    # Known first two moments cannot entail the required cap.
    assert F(1,4)*0+F(3,4)*F(4,3)==1
    assert F(3,4)*F(4,3)**2==F(4,3)
    assert F(1,4)>F(1,10)
    # The stronger block example also keeps unit underlying vector norms.
    block_law=((F(1,6),F(0)),(F(2,3),F(1)),(F(1,6),F(2)))
    assert sum(w*x for w,x in block_law)==1
    assert sum(w*x*x for w,x in block_law)==F(4,3)
    assert F(2,3)+2*F(1,6)==1 and F(1,6)>F(1,10)
    # Endpoint leakage counts underlying the continuous Fourier bound.
    leakage=0
    for d,s in product(range(1,12),range(-17,18)):
        assert sum(not 0<=j+s<d for j in range(d))==min(d,abs(s))
        leakage+=1
    record={
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'exact finite algebra checked; analytic reduction draft; arithmetic resolvent target unproved',
        'q_coefficients':[str(x) for x in Q],
        'r_cubed':str(B),
        'field_basis':['1','r','r^2','i','i*r','i*r^2'],
        'partial_fraction_coefficients_A1_A2_A3':[c.record() for c in coeff],
        'closing_insertion_coefficients_C1_C2_C3':[c.record() for c in closing],
        'f_at_one':f1.record(),
        'exact_partial_fraction_coefficient_identities':identities,
        'exact_block_resolvent_power_and_trace_bound_cases':blocks,
        'exact_two_prime_resolvent_insertion_cases':insertion,
        'projection_leakage_count_cases':leakage,
        'model_trace_upper_bound':'247/2519',
        'actual_trace_cap_required':'1/10',
        'actual_signed_resolvent_combination_required':'-9/20',
        'first_two_moment_block_counterexample':{'law':[[str(w),str(x)] for w,x in block_law],'simple_fraction':'2/3','f_trace_lower_bound':'1/6'},
        'limitations':[
            'Finite matrix checks do not prove the continuous leakage estimate or asymptotic transfer.',
            'Model trace bounds are not estimates for the actual arithmetic matrix.',
            'The three actual resolvent powers and their signed combination remain unestimated.',
            'No unconditional new zeta-zero counting bound is established.',
        ],
    }
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {identities} exact partial-fraction coefficients, {blocks} exact block resolvent-power cases, {insertion} exact insertion identities, {leakage} leakage counts.')
    print('Unproved signed actual-prime resolvent target: -9/20.')


if __name__=='__main__':
    main()
