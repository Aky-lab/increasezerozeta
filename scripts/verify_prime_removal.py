#!/usr/bin/env python3
"""Exact prime-removal algebra and finite covariance geometry.

The prime phases themselves are not estimated by these finite examples.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

from verify_bounded_certificate import block, ginv, inverse, norm_squared, power
from verify_prime_trace import G, ROOTS
from verify_spectral_bridge import add, eye, multiply, scale, trace


def zero(n):
    return [[G() for _ in range(n)] for _ in range(n)]


def adjoint(m):
    return [[m[k][j].conjugate() if isinstance(m[k][j],G) else m[k][j]
             for k in range(len(m))] for j in range(len(m))]


def prod(*matrices):
    out=matrices[0]
    for m in matrices[1:]:
        out=multiply(out,m)
    return out


def real_bound(m):
    assert all(not isinstance(x,G) or x.im==0 for row in m for x in row)
    return max(sum(abs(x.re if isinstance(x,G) else x) for x in row) for row in m)


def tau(m):
    return trace(block(m,range(2),range(2)))*F(1,2)


def grouped_expansion():
    groups=[[[F(0),F(1,3),F(1,7)],
             [F(1,3),F(-1,5),F(0)],
             [F(1,7),F(0),F(1,5)]],
            [[F(1,4),F(0),F(1,9)],
             [F(0),F(-1,4),F(2,7)],
             [F(1,9),F(2,7),F(0)]],
            [[F(1,11),F(1,13),F(0)],
             [F(1,13),F(0),F(-1,17)],
             [F(0),F(-1,17),F(-1,11)]]]
    km=zero(3)
    for m in groups:
        km=add(km,m)
    count=0
    for z in (G(F(0),F(1)),G(F(1),F(2)),G(F(-1),F(3))):
        t=G(F(1))-z
        rm=inverse(add(scale(eye(3),t),scale(km,-1)))
        removed=[inverse(add(add(scale(eye(3),t),scale(km,-1)),ap)) for ap in groups]
        for k in (1,2,3):
            zk,qk,ek=G(),G(),G()
            for ap,rp in zip(groups,removed):
                zk+=tau(prod(ap,power(rp,k)))
                for u in range(1,k+1):
                    v=k+1-u
                    qk+=tau(prod(ap,power(rp,u),ap,power(rp,v)))
                for u in range(1,k+1):
                    for v in range(1,k+2-u):
                        w=k+2-u-v
                        ek+=tau(prod(ap,power(rp,u),ap,power(rp,v),ap,power(rm,w)))
                # Check the differentiated two-term resolvent expansion itself.
                first=zero(3)
                last=zero(3)
                for u in range(1,k+1):
                    first=add(first,prod(power(rp,u),ap,power(rp,k+1-u)))
                    for v in range(1,k+2-u):
                        last=add(last,prod(power(rp,u),ap,power(rp,v),ap,power(rm,k+2-u-v)))
                assert power(rm,k)==add(add(power(rp,k),first),last)
            lhs=t*tau(power(rm,k))
            if k>1:
                lhs-=tau(power(rm,k-1))
            assert lhs==(G(F(1)) if k==1 else G())+zk+qk+ek
            bound=F(k*(k+1),2)*sum(real_bound(ap)**3 for ap in groups)/z.im**(k+2)
            assert norm_squared(ek)<=bound*bound
            count+=1
    # Uniformity in a large outside background for a small increment.
    outside=0
    ap=groups[0]
    for b in (F(0),F(10**6),F(-10**6)):
        background=[[b,F(0),F(0)],[F(0),-b,F(0)],[F(0),F(0),F(1,2)]]
        z=G(F(0),F(1))
        rp=inverse(add(background,scale(eye(3),-z)))
        rm=inverse(add(add(background,scale(ap,-1)),scale(eye(3),-z)))
        assert rm==add(add(rp,prod(rp,ap,rp)),prod(rp,ap,rp,ap,rm))
        outside+=1
    return count,outside


def phase_separation():
    u=[[F(0),F(1,2),F(0)],[F(0),F(0),F(3,4)],[F(0),F(0),F(0)]]
    ua=adjoint(u)
    h=[[F(1,3),F(1,5),F(0)],[F(1,5),F(-1,4),F(2,7)],[F(0),F(2,7),F(1,6)]]
    r=inverse(add(h,scale(eye(3),-G(F(0),F(1)))))
    cases=0
    for phase,k in product(ROOTS,(1,2,3)):
        ap=add(scale(u,phase),scale(ua,phase.conjugate()))
        assert adjoint(ap)==ap
        direct=zero(3)
        alternate=zero(3)
        for b in range(1,k+1):
            c=k+1-b
            rb,rc=power(r,b),power(r,c)
            direct=add(direct,prod(ap,rb,ap,rc))
            balanced=add(prod(u,rb,ua,rc),prod(ua,rb,u,rc))
            unbalanced=add(scale(prod(u,rb,u,rc),phase*phase),
                           scale(prod(ua,rb,ua,rc),phase.conjugate()*phase.conjugate()))
            alternate=add(alternate,add(balanced,unbalanced))
        assert direct==alternate
        cases+=1
    return cases


def covariance_geometry():
    cases=0
    for n in range(2,18):
        diagonal=[]
        for j in range(n):
            direct=sum(F(s,n*n) for s in range(1,n) if 0<=j+s<n)
            direct+=sum(F(s,n*n) for s in range(1,n) if 0<=j-s<n)
            expected=F(j*(j+1)+(n-1-j)*(n-j),2*n*n)
            assert direct==expected
            diagonal.append(direct)
            cases+=1
        assert max(diagonal)==F(n-1,2*n)<F(1,2)
        assert sum(diagonal)/n==F(n*n-1,3*n*n)
    return cases


def small_increment_obstruction():
    t=G(F(1),F(-1))
    limit=(ginv(t-G(F(1)))-ginv(t+G(F(1))))*F(1,6)
    assert limit==G(-F(1,15),F(2,15))
    for n in range(1,25):
        q=F(n-1,n)
        value=(ginv(t-G(q))-ginv(t+G(q)))*F(1,6)
        assert norm_squared(value-limit)<=F(1,n*n)
    assert F(1,6)+F(1,6)==F(1,3)  # second moment of diag(1,-1,0,0,0,0).
    return 24


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    groups,outside=grouped_expansion()
    phase=phase_separation()
    covariance=covariance_geometry()
    obstruction=small_increment_obstruction()
    report={
        'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'finite exact identities and bounds checked; analytic prime-removal draft; signed arithmetic target unproved',
        'exact_grouped_derivative_trace_identities':groups,
        'exact_large_outside_background_identities':outside,
        'exact_balanced_unbalanced_phase_separations':phase,
        'exact_finite_covariance_profile_entries':covariance,
        'exact_small_increment_obstruction_cases':obstruction,
        'small_increment_limit_at_z_i':['-1/15','2/15'],
        'analytic_cubic_remainder_envelope':'O_z((log T)^(-3))',
        'analytic_quadratic_prime_power_replacement_envelope':'O_z((log T)^(-2))',
        'covariance_map_norm_upper_limit':'1/2',
        'covariance_map_normalized_identity_trace_limit':'1/3',
        'remaining_statistics':['linear prime-power conditional phases','ordinary-prime second conditional harmonics','nonlinear balanced covariance trace'],
        'limitations':[
            'Rational matrices calibrate noncommutative identities; they are not the actual prime operator.',
            'Finite covariance profiles do not prove the continuous Mertens passage.',
            'Prime-removal remainder bounds do not imply independence or cancellation of actual prime phases.',
            'The nonlinear covariance trace has not been factorized or evaluated.',
            'No new unconditional zeta-zero bound follows.',
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {groups} grouped derivative traces, {outside} outside-background identities, {phase} phase separations, {covariance} covariance entries, {obstruction} obstruction cases.')


if __name__=='__main__':
    main()
