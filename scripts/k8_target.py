#!/usr/bin/env python3
"""Exact k=8 target geometry from the pinned m0..m6 moments.

This script uses only Fraction arithmetic.  It derives:
- the exact shifted-Hankel floor for m7;
- the ordinary-Hankel floor m8_min(m7);
- the exact degree-4 Christoffel function lambda_4(0);
- target m8 curves for desired simple-zero proportions.
"""

from fractions import Fraction as F


M7_STAR = F(25866469, 1352400)
LAMBDA3 = F(1415, 13891)
C_GAP = F(18663120, 3931153)


def m8_floor(m7):
    return (
        F(13335840000) * m7 * m7
        - F(504286322400) * m7
        + F(4794396454037)
    ) / F(748818000)


def lambda4(m7, m8):
    num = (
        F(13335840000) * m7 * m7
        - F(504286322400) * m7
        - F(748818000) * m8
        + F(4794396454037)
    )
    den = F(4) * (
        F(24004512000) * m7 * m7
        - F(903891063600) * m7
        - F(1837779300) * m8
        + F(8574904106497)
    )
    return num / den


def lambda4_gap(m7, m8):
    g = m7 - M7_STAR
    delta = m8 - m8_floor(m7)
    return LAMBDA3 * delta / (delta + C_GAP * g * g)


def delta_budget(simple_target):
    L = (F(1) - simple_target) / 2
    if not (F(0) < L < LAMBDA3):
        raise ValueError("target must improve on the degree-3 bound")
    return L * C_GAP / (LAMBDA3 - L)


def m8_target(m7, simple_target):
    g = m7 - M7_STAR
    return m8_floor(m7) + delta_budget(simple_target) * g * g


def main():
    print("m7 Stieltjes floor:")
    print(" ", M7_STAR, "=", float(M7_STAR))

    targets = [
        ("80%", F(4, 5)),
        ("80.2%", F(401, 500)),
        ("81%", F(81, 100)),
        ("90%", F(9, 10)),
    ]
    print("\ndelta <= coefficient * (m7-m7*)^2:")
    for name, s in targets:
        c = delta_budget(s)
        print(f"  {name:>5}: {c} = {float(c):.12f}")

    for m7 in (F(19130, 1000), F(19136, 1000), F(862, 45)):
        floor = m8_floor(m7)
        b80 = m8_target(m7, F(4, 5))
        print(f"\nm7={float(m7):.6f}")
        print(f"  m8 floor = {float(floor):.12f}")
        print(f"  m8 max for 80% = {float(b80):.12f}")
        # Verify target boundary exactly.
        assert lambda4(m7, b80) == F(1, 10)
        assert lambda4_gap(m7, b80) == F(1, 10)

    # Algebraic identity check at a generic rational point.
    m7 = F(19133, 1000)
    floor = m8_floor(m7)
    m8 = floor + F(1, 100)
    assert lambda4(m7, m8) == lambda4_gap(m7, m8)

    print("\nALL EXACT K=8 TARGET CHECKS PASSED")


if __name__ == "__main__":
    main()
