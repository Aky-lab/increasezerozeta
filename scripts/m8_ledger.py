#!/usr/bin/env python3
"""Exact Bell(8) ledger and aggregate target calculator.

The five genuinely new eighth-order classes are collected into A8:
    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

This script requires explicit --j52 and --c7 scenario inputs. The former
7/72 preference is quarantined after the actual-walk audit. No supplied
fraction is implicitly promoted to a certified research input.
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
    parser = argparse.ArgumentParser(description="Unverified-input scenario calculator")
    parser.add_argument("--j52", type=F, help="explicit, unverified joint input")
    parser.add_argument("--c7", type=F, help="explicit, unverified C7 input")
    args = parser.parse_args()
    if args.j52 is None or args.c7 is None:
        if args.j52 is not None or args.c7 is not None:
            parser.error("provide both --j52 and --c7")
        print("UNRESOLVED: no certified J52 or C7 is supplied.")
        print("The historical 7/72 preference is quarantined; see actual-walk audit.")
        print("Pass --j52 and --c7 explicitly for scenario arithmetic only.")
        return
    J52, C7 = args.j52, args.c7
    M7 = F(1717, 90) + J52 + C7
    if M7 < M7_STAR:
        parser.error("scenario m7 is below the pinned Stieltjes floor")
    print("SCENARIO ONLY: supplied fractions are not certified research inputs.")
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

    print("Bell(8) base pair/four layer =", base, "=", float(base))
    print("scenario inherited m8 baseline =", inherited, "=", float(inherited))
    print(f"m8 = {inherited} + A8")
    print("A8 := {2^4}+{4,2,2}+{4,4}+{6,2}+C8")

    floor = m8_floor(M7)
    print(f"\nAt scenario m7={M7}:")
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

    print("\nSCENARIO BELL(8) ARITHMETIC CHECKS PASSED")


if __name__ == "__main__":
    main()
