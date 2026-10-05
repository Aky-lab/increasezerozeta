#!/usr/bin/env python3
"""Independent exact geometry for the small-dilation model-core cutoff."""
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
    return sum(overlap_integral(p, delta-u) for p in twelve_positions(u))


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def integrate(poly, left, right):
    return sum(c*(right**(i+1)-left**(i+1))/(i+1) for i, c in enumerate(poly))


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


if __name__=='__main__':
    main()
