#!/usr/bin/env python3
"""Exact Bell(7) ledger for the seventh moment.

This file separates:
- identities already certified at lower order;
- the pre-registered {5,2}=1/8 candidate;
- the C7=-17/360 identification candidate.

It therefore computes a *candidate* m7, not yet a theorem-level input.
"""

from fractions import Fraction as F

T_ADJ = F(7, 60)
T_OPP = F(1, 30)
PHI4 = F(-1, 60)
T222 = F(131, 420)
C5 = F(1, 36)
J42 = F(-23, 420)
C6 = F(-1, 126)

# k=7 research inputs not yet symbolically certified in this repo.
J52 = F(1, 8)        # pre-registered candidate
C7 = F(-17, 360)     # source identification candidate

M7_STAR = F(25866469, 1352400)


def main():
    # Pair/four-cycle layer:
    # 1 singleton class;
    # C(7,2)=21 one-pair classes, each 1/3;
    # C(7,4)=35 choices of four pair endpoints, with two
    # noncrossing + one crossing pairing per choice;
    # C(7,4)=35 pure four-block placements.
    pair4 = (
        F(1)
        + F(21, 3)
        + 70 * T_ADJ
        + 35 * T_OPP
        + 35 * PHI4
    )
    assert pair4 == F(67, 4)

    # Frozen singleton multiplicities:
    # (5,1,1): 21 copies of C5
    # (2,2,2,1): 7 copies of the six-cycle {2,2,2} aggregate
    # (4,2,1): 7 copies of {4,2}
    # (6,1): 7 copies of C6
    m7 = (
        pair4
        + 21 * C5
        + 7 * T222
        + 7 * J42
        + J52
        + 7 * C6
        + C7
    )

    assert m7 == F(862, 45)

    print("pair/four layer =", pair4, "=", float(pair4))
    print("candidate m7 =", m7, "=", float(m7))
    print("exact Stieltjes floor =", M7_STAR, "=", float(M7_STAR))
    print("gap =", m7 - M7_STAR, "=", float(m7 - M7_STAR))

    # The source outlook's historical lower endpoint 19.123 is
    # already below the exact floor implied by pinned m0..m6.
    assert F(19123, 1000) < M7_STAR
    print("NOTE: source reconnaissance lower endpoint 19.123 is stale")
    print("CANDIDATE M7 LEDGER CHECKS PASSED")


if __name__ == "__main__":
    main()
