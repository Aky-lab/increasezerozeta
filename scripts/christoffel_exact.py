#!/usr/bin/env python3
"""Exact Christoffel consumption for finite moment sequences.

Pure standard-library Fraction arithmetic.  The script reproduces the
fourth-moment 13/18 certificate and gives the canonical exact six-moment
certificate.

The zeta-function interpretation remains conditional on the surrounding
spectral/counting interface.  This file certifies only the moment algebra.
"""

from fractions import Fraction as F


def solve_linear(A, b):
    """Gauss-Jordan solve Ax=b over Fractions."""
    n = len(A)
    aug = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if aug[i][col] != 0)
        aug[col], aug[pivot] = aug[pivot], aug[col]
        piv = aug[col][col]
        aug[col] = [x / piv for x in aug[col]]
        for i in range(n):
            if i == col:
                continue
            fac = aug[i][col]
            if fac:
                aug[i] = [aug[i][j] - fac * aug[col][j]
                          for j in range(n + 1)]
    return [aug[i][-1] for i in range(n)]


def det(A):
    """Determinant by elimination over Fractions."""
    A = [list(r) for r in A]
    n = len(A)
    out = F(1)
    sign = 1
    for col in range(n):
        pivot = next((i for i in range(col, n) if A[i][col] != 0), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            sign *= -1
        piv = A[col][col]
        out *= piv
        for i in range(col + 1, n):
            fac = A[i][col] / piv
            for j in range(col + 1, n):
                A[i][j] -= fac * A[col][j]
    return out * sign


def hankel(moments, n):
    return [[moments[i + j] for j in range(n + 1)]
            for i in range(n + 1)]


def christoffel_at_zero(moments, n):
    """Return lambda_n(0) and q coefficients with q(0)=1."""
    H = hankel(moments, n)
    u = solve_linear(H, [F(1)] + [F(0)] * n)  # H^{-1} e0
    K00 = u[0]
    lam = 1 / K00
    q = [x / K00 for x in u]
    assert q[0] == 1

    # Check E[q^2] exactly.
    eq2 = F(0)
    for i, ci in enumerate(q):
        for j, cj in enumerate(q):
            eq2 += ci * cj * moments[i + j]
    assert eq2 == lam
    return lam, q, det(H)


def poly_string(q):
    bits = []
    for k, c in enumerate(q):
        if not c:
            continue
        if k == 0:
            bits.append(str(c))
        else:
            bits.append(f"({c})*x^{k}")
    return " + ".join(bits)


def eval_poly(q, x):
    return sum(c * x**k for k, c in enumerate(q))


def alternating_negative_halfline(q):
    """Sufficient exact condition q(-t)>=1 for all t>=0."""
    if q[0] != 1:
        return False
    # q(-t) = 1 + sum_{k>=1} q_k (-1)^k t^k.
    return all(q[k] * ((-1) ** k) >= 0 for k in range(1, len(q)))


def main():
    moments = {
        0: F(1),
        1: F(1),
        2: F(4, 3),
        3: F(2),
        4: F(13, 4),
        5: F(101, 18),
        6: F(12809, 1260),
    }

    # Degree 2: recover 13/18.
    lam2, q2, det2 = christoffel_at_zero(moments, 2)
    assert lam2 == F(5, 36)
    assert q2 == [F(1), -F(7, 4), F(2, 3)]
    assert alternating_negative_halfline(q2)
    assert 1 - 2 * lam2 == F(13, 18)

    # Degree 3: exact six-moment certificate.
    lam3, q3, det3 = christoffel_at_zero(moments, 3)
    assert det3 == F(283, 108864)
    assert lam3 == F(1415, 13891)
    assert q3 == [
        F(1),
        -F(43428, 13891),
        F(37704, 13891),
        -F(9660, 13891),
    ]
    assert alternating_negative_halfline(q3)

    simple = 1 - 2 * lam3
    distinct = 1 - lam3
    assert simple == F(11061, 13891)
    assert distinct == F(12476, 13891)

    print("degree 2:")
    print("  det(H2) =", det2)
    print("  lambda_2(0) =", lam2)
    print("  q2(x) =", poly_string(q2))
    print("  simple >=", float(1 - 2 * lam2))

    print("\ndegree 3:")
    print("  det(H3) =", det3)
    print("  lambda_3(0) =", lam3, "=", float(lam3))
    print("  q3(x) =", poly_string(q3))
    print("  q3(-t)>=1 sign test:", alternating_negative_halfline(q3))
    print("  simple >=", simple, "=", float(simple))
    print("  distinct >=", distinct, "=", float(distinct))

    # Compare the exact polynomial with the reference rationalised roots.
    roots_approx = [F(5323, 10000), F(6561, 5000), F(10293, 5000)]
    vals = [float(eval_poly(q3, r)) for r in roots_approx]
    print("\nq3 at reference rationalised roots:", vals)

    print("\nALL EXACT CHRISTOFFEL CHECKS PASSED")


if __name__ == "__main__":
    main()
