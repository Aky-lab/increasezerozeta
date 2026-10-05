#!/usr/bin/env python3
"""Exact certificates for prime-position coefficient profiles of the model core."""
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
import argparse
import json
from pathlib import Path

from verify_dilation_coverage import clip, integrate_affine, add, sub, twelve_positions, multiply, integrate, coverage


def solve(rows):
    matrix = [[F(x) for x in row[:3]]+[F(-row[3])] for row in rows]
    for c in range(3):
        pivot = next((r for r in range(c,3) if matrix[r][c]), None)
        if pivot is None:
            return None
        matrix[c],matrix[pivot] = matrix[pivot],matrix[c]
        scale = matrix[c][c]
        matrix[c] = [x/scale for x in matrix[c]]
        for r in range(3):
            if r!=c:
                scale = matrix[r][c]
                matrix[r] = [x-scale*y for x,y in zip(matrix[r],matrix[c])]
    return tuple(row[3] for row in matrix)


def positions3():
    """Lift the original affine overlap positions into (u,x,y,constant)."""
    first,last = list(twelve_positions(F(0))),list(twelve_positions(F(1)))
    return [tuple((b[2]-a[2],a[0],a[1],a[2]) for a,b in zip(p,q))
            for p,q in zip(first,last)]


def difference(a,b,constant=0):
    return tuple(x-y for x,y in zip(a,b))[:3]+(a[3]-b[3]+constant,)


def canonical(row):
    row = tuple(F(x) for x in row)
    scale = next((x for x in row[:3] if x),None)
    return tuple(x/scale for x in row) if scale else None


def base(tau):
    return [(1,0,0,0),(-1,0,0,tau),(0,1,0,0),(0,0,1,0),
            (-1,-2,0,1),(-1,0,-2,1),(-1,-1,-tau,tau),(-1,-tau,-1,tau)]


def breaks(tau):
    """Include all partition-polyhedron vertex levels in the compact base."""
    constraints = base(tau)
    planes = {canonical(row) for row in constraints}
    for positions in positions3():
        for p,q in combinations(positions,2):
            planes.add(canonical(difference(p,q)))
            planes.add(canonical(difference(p,q,1)))
            planes.add(canonical(difference(p,q,-1)))
    planes.discard(None)
    levels = {F(0),tau}
    for rows in combinations(sorted(planes),3):
        point = solve(rows)
        if point and all(sum(a*b for a,b in zip(row[:3],point))+row[3]>=0 for row in constraints):
            levels.add(point[0])
    return sorted(levels),len(planes)


@lru_cache(maxsize=None)
def geometry(tau,u):
    t = (1-u)/2
    polygon = [(F(0),F(0)),(t,F(0)),(t,t),(F(0),t)]
    polygon = clip(polygon,(-1,-tau,tau-u))
    polygon = clip(polygon,(-tau,-1,tau-u))
    total = F(0)
    for positions in twelve_positions(u):
        positions = tuple(dict.fromkeys(positions))
        for lo in positions:
            for hi in positions:
                cell = polygon
                for p in positions:
                    cell = clip(cell,sub(p,lo))
                    cell = clip(cell,sub(hi,p))
                height = add((0,0,1),sub(lo,hi))
                cell = clip(cell,height)
                total += integrate_affine(cell,height)
    return total


def interpolate(xs,ys):
    polynomial = [F(0)]*len(xs)
    for x,y in zip(xs,ys):
        term, denominator = [F(1)],F(1)
        for other in xs:
            if other!=x:
                term = multiply(term,[-other,1])
                denominator *= x-other
        for i,c in enumerate(term):
            polynomial[i] += y*c/denominator
    return polynomial


def evaluate(p,x):
    return sum(c*x**i for i,c in enumerate(p))


def profile(tau):
    """Certify cubic slices after the geometric degree/partition proof."""
    if not 0<=tau<=1:
        raise ValueError('The coefficient profile must lie in [0,1].')
    levels,planes = breaks(tau)
    result, slices = F(0), []
    for a,b in zip(levels,levels[1:]):
        xs = [a+(b-a)*F(j,5) for j in range(1,5)]
        ys = [geometry(tau,u) for u in xs]
        polynomial = interpolate(xs,ys)
        heldout = (a+(b-a)/2,a+(b-a)*F(3,7),a,b)
        for point in heldout:
            assert geometry(tau,point)==evaluate(polynomial,point),(tau,a,b,point)
        result += integrate(multiply([0,1],polynomial),a,b)
        slices.append({'interval':[str(a),str(b)],
                       'fit_points':[str(x) for x in xs],
                       'coefficients_ascending':[str(c) for c in polynomial],
                       'heldout_points':[str(x) for x in heldout]})
    return {'tau':str(tau),'leading_fraction':str(96*result),
            'partition_planes':planes,'breakpoints':list(map(str,levels)),
            'slices':slices}


def closed_profile(tau):
    if not 0<=tau<=F(1,2):
        raise ValueError('The closed expression is proved only for tau<=1/2.')
    return 2*tau**5*(59+40*tau+5*tau**2)/(5*(1+tau)**2*(2+tau))


def closed_geometry(tau,u):
    return ((7+tau)*(tau-u)**3/(3*(1+tau)**2)
            +(3+tau)*max(tau-2*u,F(0))**3/(3*(1+tau)*(2+tau)))


def finite_scale_checks():
    count = 0
    for tau in (F(1,4),F(1,2),F(2,3),F(3,4),F(9,10),F(1)):
        for u in (F(j,8) for j in range(9)):
            nu,t = 1+u,(1-u)/2
            for x in (t*j/8 for j in range(9)):
                for y in (t*j/8 for j in range(9)):
                    b1,b2 = u+x,u+y
                    physical = b1<=tau*(nu-b2) and b2<=tau*(nu-b1)
                    clipped = x+tau*y<=tau-u and tau*x+y<=tau-u
                    assert physical==clipped
                    if clipped:
                        assert max(b1,b2)<=tau
                        assert nu-b1-b2>=nu*(1-tau)/(1+tau)
                    if max(b1,b2)<=tau/(1+tau):
                        assert clipped
                    count += 1
    for b1,b2,g in ((3,5,1),(9,27,9),(4,16,4),(4,9,1)):
        product = F(10000)
        m,n,r,q = product/b2,product/b1,b1//g,b2//g
        assert m/r==n/q==product*g/(b1*b2)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    scales = finite_scale_checks()
    closed_checks = 0
    for tau in (F(1,8),F(1,4),F(3,8),F(1,2)):
        for u in (F(0),tau/4,tau/2,3*tau/4,tau):
            assert geometry(tau,u)==closed_geometry(tau,u),(tau,u)
            closed_checks += 1
        # Exact symbolic-u integration of the two derived cubic terms.
        s=[tau,-1]
        first=multiply(multiply(s,s),s)
        s=[tau,-2]
        second=multiply(multiply(s,s),s)
        exact=96*((7+tau)*integrate(multiply([0,1],first),0,tau)/(3*(1+tau)**2)
                  +(3+tau)*integrate(multiply([0,1],second),0,tau/2)/(3*(1+tau)*(2+tau)))
        assert exact==closed_profile(tau)
    expected = ((F(1,2),F(107,600)),(F(2,3),F(73,150)),
                (F(3,4),F(10243,15750)),(F(9,10),F(256636537,277981875)),(F(1),F(1)))
    records = []
    for tau,fraction in expected:
        record=profile(tau)
        assert F(record['leading_fraction'])==fraction,(tau,record)
        assert coverage(tau/(1+tau))<=fraction<=coverage(tau)
        records.append(record)
    report={'schema_version':1,'scope':'Exact leading-model coefficient-profile coverage; no prime theorem.',
            'degree_certificate':'Affine slice vertices give quadratic areas times affine heights: degree at most three.',
            'profiles':records}
    record_path=Path(__file__).resolve().parents[1]/'results/progression_profile_2026-10-05.json'
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(record_path.read_text(encoding='utf-8'))==report
    slices=sum(len(record['slices']) for record in records)
    print(f'PASS: {len(records)} exact profile fractions, {slices} proved-degree slice cubics,')
    print(f'      {4*slices} independent held-out/endpoint evaluations, {closed_checks} closed-form geometries,')
    print(f'      {scales} coefficient/progression scale cases and exact coverage containment bounds.')
    print('Profile fractions:',[(r['tau'],r['leading_fraction']) for r in records])


if __name__=='__main__':
    main()


