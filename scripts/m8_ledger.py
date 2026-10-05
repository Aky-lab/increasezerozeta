#!/usr/bin/env python3
"""Exact Bell(8) ledger and aggregate target calculator.

The five genuinely new eighth-order classes are collected into A8:
    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

The default seventh-order model inputs are {5,2}=1/8 and C7=-17/360.
Explicit overrides calculate alternative model scenarios. Eighth-order
class values and arithmetic transport remain unresolved.
"""

import argparse
from fractions import Fraction as F

T_ADJ = F(7, 60)
T_OPP = F(1, 30)
PHI4 = F(-1, 60)
T222 = F(131, 420)
J42 = F(-23, 420)
C5 = F(1, 36)
C6 = F(-1, 126)
MODEL_J52 = F(1, 8)
MODEL_C7 = F(-17, 360)


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--j52", type=F, default=MODEL_J52, help="joint model input; default 1/8")
    parser.add_argument("--c7", type=F, default=MODEL_C7, help="pure-cycle input; default -17/360")
    args = parser.parse_args()
    J52, C7 = args.j52, args.c7
    M7 = F(1717, 90) + J52 + C7
    if M7 < M7_STAR:
        parser.error("scenario m7 is below the pinned Stieltjes floor")
    print("CONTINUUM MODEL TARGETS: arithmetic transport remains open.")
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
    if J52 == F(7, 72) and C7 == F(-17, 360):
        assert inherited == F(1103, 30)
    if (J52, C7) == (MODEL_J52, MODEL_C7):
        assert M7 == F(862, 45) and inherited == F(3329, 90)

    print("Bell(8) base pair/four layer =", base, "=", float(base))
    print("model inherited m8 baseline =", inherited, "=", float(inherited))
    print(f"m8 = {inherited} + A8")
    print("A8 := {2^4}+{4,2,2}+{4,4}+{6,2}+C8")

    floor = m8_floor(M7)
    print(f"\nAt model m7={M7}:")
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
    if J52 == F(7, 72) and C7 == F(-17, 360):
        assert cap80 == F(18073331, 68531400)

    print("\nMODEL BELL(8) ARITHMETIC CHECKS PASSED")


if __name__ == "__main__":
    main()
