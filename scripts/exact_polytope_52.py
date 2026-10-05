#!/usr/bin/env python3
"""Exact rational polytope evaluator for the {5,2} joint class.

Adapted independently from the polytope-volume identity used in the
reference lower-order certifications.

For one partition-cyclic C5 term,

    integral ov(A) * ov(B) * min(|v|,1) dx

has outer support |v| <= 1, so min(|v|,1)=|v|.  The two overlap
factors are lifted with variables t,s and |v| with one extra y
variable.  All inequalities have rational coefficients, hence every
term is an exact Fraction.

Usage:
    python scripts/exact_polytope_52.py 1 --terms 0:1
    python scripts/exact_polytope_52.py 2 --terms 0:10 --out part.txt
    python scripts/exact_polytope_52.py 3

A full component contains 150 signed C5 terms.  Expected candidate
components are registered separately; this script does not bake them
into the computation.
"""

import argparse
import itertools
import time
from fractions import Fraction as F

Z = F(0)
O = F(1)


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def compile_terms(b=5):
    out = []
    for part in set_partitions(list(range(b))):
        m = len(part)
        sign = -1 if (m - 1) % 2 else 1
        first = next(i for i, block in enumerate(part) if 0 in block)
        others = [i for i in range(m) if i != first]
        for perm in itertools.permutations(others):
            order = [first] + list(perm)
            masks = []
            cur = set()
            for bi in order[:-1]:
                cur.update(part[bi])
                masks.append(tuple(sorted(cur)))
            out.append((sign, tuple(masks)))
    return out


TERMS = compile_terms(5)
assert len(TERMS) == 150


def solve(M, rhs, n):
    A = [list(M[i]) + [rhs[i]] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fac = A[r][c]
                A[r] = [a - fac * b for a, b in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def vertices(rows, n):
    out = set()
    for comb in itertools.combinations(range(len(rows)), n):
        x = solve([rows[i][0] for i in comb],
                  [rows[i][1] for i in comb], n)
        if x is None:
            continue
        good = True
        for a, b in rows:
            if sum(ai * xi for ai, xi in zip(a, x) if ai) > b:
                good = False
                break
        if good:
            out.add(tuple(x))
    return out


def tight(rows, verts):
    ans = []
    for i, (a, b) in enumerate(rows):
        if all(sum(ai * xi for ai, xi in zip(a, v) if ai) == b
               for v in verts):
            ans.append(i)
    return frozenset(ans)


def volume(rows, n, verts=None, memo=None):
    if memo is None:
        memo = {}
    if verts is None:
        verts = vertices(rows, n)
    if not verts:
        return Z
    if n == 0:
        return O
    if n == 1:
        xs = [v[0] for v in verts]
        return max(xs) - min(xs)

    verts = tuple(sorted(verts))
    key = (n, tight(rows, verts), frozenset(verts))
    if key in memo:
        return memo[key]

    centre = [sum(v[i] for v in verts) / len(verts) for i in range(n)]
    total = Z
    seen = set()

    for a, b in rows:
        face = tuple(v for v in verts
                     if sum(ai * xi for ai, xi in zip(a, v) if ai) == b)
        if len(face) < n or face in seen:
            continue
        seen.add(face)

        k = next((j for j in range(n) if a[j] != 0), None)
        if k is None:
            continue
        ak = a[k]

        projected_rows = []
        for a2, b2 in rows:
            if a2[k] == 0:
                co = tuple(a2[j] for j in range(n) if j != k)
                if any(co):
                    projected_rows.append((co, b2))
            else:
                fac = a2[k] / ak
                co = tuple(a2[j] - fac * a[j]
                           for j in range(n) if j != k)
                rb = b2 - fac * b
                if any(co):
                    projected_rows.append((co, rb))

        projected_verts = set(
            tuple(v[j] for j in range(n) if j != k) for v in face
        )
        sub = volume(projected_rows, n - 1, projected_verts, memo)
        if sub == 0:
            continue

        height = b - sum(ai * ci for ai, ci in zip(a, centre) if ai)
        total += height * sub / abs(ak)

    ans = total / n
    memo[key] = ans
    return ans


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def spec52_positions(distance, masks):
    """Return outer 7-walk A and C5-term walk B in 5 free variables.

    Variables: (v,c1,c2,c3,c4); c5=-(c1+c2+c3+c4).
    """
    d = 5
    V = [O, Z, Z, Z, Z]
    C1 = [Z, O, Z, Z, Z]
    C2 = [Z, Z, O, Z, Z]
    C3 = [Z, Z, Z, O, Z]
    C4 = [Z, Z, Z, Z, O]
    S1 = C1
    S2 = add(C1, C2)
    S3 = add(S2, C3)
    S4 = add(S3, C4)
    zero = [Z] * d

    outer = {
        1: [zero, V, zero, S1, S2, S3, S4],
        2: [zero, V, add(V, S1), S1, S2, S3, S4],
        3: [zero, V, add(V, S1), add(V, S2), S2, S3, S4],
    }[distance]

    W = [
        C1,
        C2,
        C3,
        C4,
        [-x for x in S4],
    ]

    merged = [zero]
    for mask in masks:
        co = zero[:]
        for j in mask:
            co = add(co, W[j])
        merged.append(co)

    return d, outer, merged


def hrep(A, B, d):
    """Lift ov(A)*ov(B) to a polytope in variables (x,t,s)."""
    n = d + 2
    rows = []
    for positions, aux in ((A, d), (B, d + 1)):
        for p in positions:
            c = list(p) + [Z, Z]
            c[aux] = -O
            rows.append((tuple(-v for v in c), Z))
            rows.append((tuple(c), O))
    return n, rows


def weighted_v(rows, n):
    """Exact integral of |v| over a lifted overlap polytope."""
    total = Z
    for sign in (O, -O):
        rs = [(tuple(list(a) + [Z]), b) for a, b in rows]

        # sign*v >= 0
        r = [Z] * (n + 1)
        r[0] = -sign
        rs.append((tuple(r), Z))

        # y >= 0
        r = [Z] * (n + 1)
        r[n] = -O
        rs.append((tuple(r), Z))

        # y <= sign*v
        r = [Z] * (n + 1)
        r[0] = -sign
        r[n] = O
        rs.append((tuple(r), Z))

        total += volume(rs, n + 1)
    return total


def term_value(distance, masks):
    d, A, B = spec52_positions(distance, masks)
    n, rows = hrep(A, B, d)
    return weighted_v(rows, n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("distance", type=int, choices=(1, 2, 3))
    ap.add_argument("--terms", default="0:150")
    ap.add_argument("--out")
    args = ap.parse_args()

    lo, hi = [int(x) for x in args.terms.split(":")]
    hi = min(hi, len(TERMS))
    total = Z
    fh = open(args.out, "a") if args.out else None
    start = time.time()

    for i in range(lo, hi):
        sign, masks = TERMS[i]
        t0 = time.time()
        val = term_value(args.distance, masks)
        total += sign * val
        line = (
            f"52 d={args.distance} term {i} sign {sign} "
            f"value {val} wall {time.time()-t0:.3f}"
        )
        print(line, flush=True)
        if fh:
            fh.write(line + "\n")
            fh.flush()

    if fh:
        fh.close()
    print(
        f"# partial d={args.distance} [{lo},{hi}) = {total} "
        f"= {float(total):+.15f} wall={time.time()-start:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
