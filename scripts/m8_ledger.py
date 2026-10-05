#!/usr/bin/env python3
"""Exact Bell(8) ledger and aggregate target calculator.

The five genuinely new eighth-order classes are collected into A8:
    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

This script uses the currently pre-registered seventh-order candidates
{5,2}=7/72 and C7=-17/360.  Therefore the resulting inherited baseline
is candidate-level until those inputs are certified.
"""

from fractions import Fraction as F

T_ADJ = F(7, 60)
T_OPP = F(1, 30)
PHI4 = F(-1, 60)
T222 = F(131, 420)
J42 = F(-23, 420)
C5 = F(1, 36)
C6 = F(-1, 126)

J52 = F(7, 72)      # active exact-v numerical candidate
C7 = F(-17, 360)    # identified in source, not yet exact-certified

M7 = F(3443, 180)   # candidate assembled from the same seventh inputs
M7_STAR = F(25866469, 1352400)
LAMBDA3 = F(1415, 13891)
C_GAP = F(18663120, 3931153)


def m8_floor(m7):
    return (
        F(13335840000) * m7 * m7
        - F(504286322400) * m7
        + F(4794396454037)
    ) / F(748818000)


def delta_budget(simple_target):
    L = (F(1) - simple_target) / 2
    return L * C_GAP / (LAMBDA3 - L)


def m8_target(m7, simple_target):
    g = m7 - M7_STAR
    return m8_floor(m7) + delta_budget(simple_target) * g * g


def main():
    base = (
        F(1)
        + F(28, 3)
        + 140 * T_ADJ
        + 70 * T_OPP
        + 70 * PHI4
    )
    assert base == F(167, 6)

    inherited = (
        base
        + 28 * T222
        + 28 * J42
        + 56 * C5
        + 8 * J52
        + 28 * C6
        + 8 * C7
    )
    assert inherited == F(1103, 30)

    print("Bell(8) base pair/four layer =", base, "=", float(base))
    print("candidate inherited m8 baseline =", inherited, "=", float(inherited))
    print("m8 = 1103/30 + A8")
    print("A8 := {2^4}+{4,2,2}+{4,4}+{6,2}+C8")

    floor = m8_floor(M7)
    print("\nAt candidate m7=3443/180:")
    print("  m8 Stieltjes floor =", floor, "=", float(floor))
    print("  A8 floor =", floor - inherited, "=", float(floor - inherited))

    for label, target in [
        ("80%", F(4, 5)),
        ("80.2%", F(401, 500)),
        ("81%", F(81, 100)),
        ("90%", F(9, 10)),
    ]:
        cap = m8_target(M7, target)
        a8cap = cap - inherited
        print(
            f"  {label:>5}: m8 <= {float(cap):.12f}; "
            f"A8 <= {a8cap} = {float(a8cap):.12f}"
        )

    # Exact 80% aggregate cap.
    cap80 = m8_target(M7, F(4, 5)) - inherited
    assert cap80 == F(18073331, 68531400)

    print("\nCANDIDATE BELL(8) LEDGER CHECKS PASSED")


if __name__ == "__main__":
    main()
