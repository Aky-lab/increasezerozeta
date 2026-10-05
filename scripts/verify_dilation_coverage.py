#!/usr/bin/env python3
"""Independent exact geometry for the full dilation model-core distribution."""
from fractions import Fraction as F
from itertools import product
from math import gcd


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def value(a, p):
    return a[0]*p[0]+a[1]*p[1]+a[2]


def clip(poly, affine):
    """Intersect a counterclockwise polygon with affine>=0."""
    if not poly:
        return []
    result = []
    for p, q in zip(poly, poly[1:]+poly[:1]):
        vp, vq = value(affine, p), value(affine, q)
        if vp>=0:
            result.append(p)
        if (vp<0<vq) or (vq<0<vp):
            f = vp/(vp-vq)
            result.append(tuple(x+f*(y-x) for x, y in zip(p, q)))
    return result


def integrate_affine(poly, affine):
    result = F(0)
    if len(poly)<3:
        return result
    p = poly[0]
    for q, r in zip(poly[1:], poly[2:]):
        twice_area = (q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])
        assert twice_area>=0
        result += twice_area*sum(value(affine, s) for s in (p, q, r))/6
    return result


def overlap_integral(positions, t):
    """Partition by actual minimum and maximum; integrate each affine cell."""
    positions = tuple(dict.fromkeys(positions))  # avoid full-dimensional ties
    square = [(F(0), F(0)), (t, F(0)), (t, t), (F(0), t)]
    total = F(0)
    for lo, hi in product(positions, repeat=2):
        cell = square
        for a in positions:
            cell = clip(cell, sub(a, lo))
            cell = clip(cell, sub(hi, a))
        height = add((F(0), F(0), F(1)), sub(lo, hi))
        cell = clip(cell, height)
        total += integrate_affine(cell, height)
    return total


def twelve_positions(u):
    z = (F(0), F(0), F(0))
    nu = (F(0), F(0), 1+u)
    x, y = (F(1), F(0), u), (F(0), F(1), u)
    for a, c in ((x, sub(nu, x)), (sub(nu, x), x)):
        for b, d in ((y, sub(nu, y)), (sub(nu, y), y)):
            delta = sub(add(a, c), add(b, d))
            yield (z, add(neg(a), delta), sub(c, d), neg(d))
            yield (z, sub(sub(b, c), d), neg(add(c, d)), neg(d))
            yield (z, add(neg(add(b, c)), d), sub(d, c), d)


def raw_overlap(x, y, u):
    return sum(max(1-max(v)+min(v), F(0))
               for v in (tuple(value(a, (x, y)) for a in positions)
                         for positions in twelve_positions(u)))


def restricted_overlap(x, y, u):
    return 4*min(x, y)+2*(max(y-x-u, F(0))+y+max(min(x-y-u, y), F(0)))


def geometry_integral(delta, u):
    t = min(delta-u, (1-u)/2)
    return sum(overlap_integral(p, t) for p in twelve_positions(u))


def full_geometry(delta, u):
    t, s = min(delta-u, (1-u)/2), 1-2*u
    cube = sum(c*max(s-j*t, F(0))**3 for j, c in enumerate((1, -3, 3, -1)))/6
    return F(4,3)*t**3+F(1,2)*max(t-u, F(0))**3+cube


def core_coefficient(delta):
    if not 0<=delta<=1:
        raise ValueError('The cutoff exponent must lie in [0,1].')
    if delta<=F(1,2):
        return -F(59,240)*delta**5+max(3*delta-1, F(0))**5/60
    eta = 1-delta
    return (-F(1,48)+F(2,3)*eta**4-F(14,15)*eta**5
            +(11*delta-4)*max(2-3*delta, F(0))**4/80)


def coverage(delta):
    return -48*core_coefficient(delta)


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def integrate(poly, left, right):
    return sum(c*(right**(i+1)-left**(i+1))/(i+1) for i, c in enumerate(poly))


def combine(*terms):
    result = [F(0)]*max(len(p) for _, p in terms)
    for coefficient, polynomial in terms:
        for i, c in enumerate(polynomial):
            result[i] += coefficient*c
    return result


def cube(p):
    return multiply(multiply(p, p), p)


def power(p, exponent):
    result = [F(1)]
    for _ in range(exponent):
        result = multiply(result, p)
    return result


def derivative(p):
    return [i*c for i, c in enumerate(p) if i]


def evaluate(p, x):
    return sum(c*x**i for i, c in enumerate(p))


def check_branch_densities():
    eta, x, a, b = [1, -1], [0, 1], [-1, 3], [2, -3]
    first = combine((F(-59,240), power(x, 5)))
    second = combine((1, first), (F(1,60), power(a, 5)))
    last = combine((1, [F(-1,48)]), (F(2,3), power(eta, 4)),
                   (F(-14,15), power(eta, 5)))
    third = combine((1, last), (F(1,80), multiply([-4, 11], power(b, 4))))
    density_last = multiply(power(eta, 3), [-96, 224])
    densities = (combine((59, power(x, 4))),
                 combine((59, power(x, 4)), (-12, power(a, 4))),
                 combine((1, density_last), (F(3,5), multiply(power(b, 3), [-70, 165]))),
                 density_last)
    branches = (first, second, third, last)
    for polynomial, density in zip(branches, densities):
        assert combine((-48, derivative(polynomial)))==density
    for left, right, boundary in zip(branches, branches[1:], (F(1,3), F(1,2), F(2,3))):
        assert evaluate(left, boundary)==evaluate(right, boundary)
        assert evaluate(derivative(left), boundary)==evaluate(derivative(right), boundary)
    return len(branches)


def independent_core_integral(delta):
    """Integrate the clipped-cube geometry, splitting at actual affine roots."""
    total = F(0)
    split = max(2*delta-1, F(0))
    for left, right, t in ((F(0), split, [F(1,2), F(-1,2)]),
                           (split, delta, [delta, F(-1)])):
        if left==right:
            continue
        expressions = [(F(1,2), combine((1, t), (-1, [0, 1])))]
        expressions += [(F(c,6), combine((1, [1, -2]), (-j, t)))
                        for j, c in enumerate((1, -3, 3, -1))]
        breaks = {left, right}
        for _, (constant, slope) in expressions:
            if slope:
                root = -constant/slope
                if left<root<right:
                    breaks.add(root)
        ends = sorted(breaks)
        for a, b in zip(ends, ends[1:]):
            midpoint = (a+b)/2
            terms = [(F(4,3), cube(t))]
            terms += [(c, cube(p)) for c, p in expressions if p[0]+p[1]*midpoint>0]
            total += integrate(multiply([0, 1], combine(*terms)), a, b)
    return -2*total


def quantile_bracket(target, denominator=10**6):
    """Exact integer bisection; report neighboring rational grid points."""
    lo, hi = 0, denominator
    while hi-lo>1:
        mid = (lo+hi)//2
        if coverage(F(mid, denominator))<target:
            lo = mid
        else:
            hi = mid
    assert coverage(F(lo, denominator))<target<=coverage(F(hi, denominator))
    return lo, hi


def main():
    polygons = points = integrals = classification = 0
    for delta in (F(1,12), F(1,6), F(1,4), F(1,3)):
        for u in (F(0), delta/4, delta/2, 3*delta/4, delta):
            t = delta-u
            expected = F(7,3)*t**3+F(1,2)*max(delta-2*u, F(0))**3
            assert geometry_integral(delta, u)==expected, (delta, u)
            polygons += 1
            for i, j in product(range(9), repeat=2):
                x, y = t*i/8, t*j/8
                assert raw_overlap(x, y, u)==restricted_overlap(x, y, u)
                points += 1
        p = multiply(multiply([delta, -1], [delta, -1]), [delta, -1])
        q = multiply(multiply([delta, -2], [delta, -2]), [delta, -2])
        coefficient = -2*(F(7,3)*integrate(multiply([0, 1], p), 0, delta)
                          +F(1,2)*integrate(multiply([0, 1], q), 0, delta/2))
        assert coefficient==F(-59,240)*delta**5
        assert coefficient/F(-1,48)==F(59,5)*delta**5
        integrals += 1
    assert geometry_integral(F(1,2), F(0))==F(1,3)
    assert F(7,3)*F(1,2)**3+F(1,2)*F(1,2)**3==F(17,48)!=F(1,3)
    assert F(59,5)*F(1,3)**5==F(59,1215)
    full_polygons = 0
    for delta in (F(1,12), F(1,3), F(2,5), F(1,2), F(3,5), F(2,3),
                  F(3,4), F(4,5), F(9,10), F(1)):
        starts = {F(0), delta/4, delta/2, 3*delta/4, delta}
        starts |= {u for u in (F(1,2), F(1,3), 2*delta-1, delta/2, 1-delta, 3*delta-1)
                   if 0<=u<=delta}
        for u in sorted(starts):
            assert geometry_integral(delta, u)==full_geometry(delta, u), (delta, u)
            full_polygons += 1
    for delta in (F(j,60) for j in range(61)):
        assert independent_core_integral(delta)==core_coefficient(delta), delta
        assert 0<=coverage(delta)<=1
        if delta:
            assert coverage(delta)>coverage(delta-F(1,60))
    assert coverage(F(1,2))==F(11,32)
    assert coverage(F(3,4))==F(147,160)
    assert coverage(F(4,5))==F(15049,15625)
    assert coverage(F(9,10))==F(15582,15625)
    assert coverage(F(1))==1
    print('Coverage values:', [(str(d), str(coverage(d)))
                               for d in (F(1,3), F(1,2), F(2,3), F(3,4), F(4,5), F(9,10))])
    print('Quantile brackets (millionths):', [(str(target), quantile_bracket(target))
                                             for target in (F(1,2), F(9,10), F(99,100))])
    density_checks = check_branch_densities()
    atoms = [(p**a, p) for p in (2, 3, 5, 7) for a in range(1,5)]
    for (b1, p1), (b2, p2), cutoff in product(atoms, atoms, (1, 2, 3, 8, 20, 100)):
        g = gcd(b1, b2)
        full = max(b1, b2)<=cutoff
        reduced = max(b1//g, b2//g)<=cutoff
        if p1!=p2:
            assert g==1 and full==reduced
        if full!=reduced:
            assert p1==p2 and reduced and not full
        classification += 1
    print(f'PASS: {polygons} original-overlap polygon integrals, {points} pointwise overlaps,')
    print(f'      {integrals} exact core integrals, threshold countercheck and {classification} cutoff cases.')
    print(f'      {full_polygons} full-range polygon integrals, 61 piecewise core integrations,')
    print(f'      {density_checks} symbolic density identities with branch joins and exact quantile brackets.')


if __name__=='__main__':
    main()
