#!/usr/bin/env python3
"""Exact Bell(7) ledger for the seventh moment.

This file separates:
- identities already certified at lower order;
- the historical, quarantined {5,2}=7/72 scenario;
- the C7=-17/360 identification candidate.

It computes only explicitly requested scenarios; no moment input is certified.
"""

import argparse
from fractions import Fraction as F

T_ADJ = F(7, 60)
T_OPP = F(1, 30)
PHI4 = F(-1, 60)
T222 = F(131, 420)
C5 = F(1, 36)
J42 = F(-23, 420)
C6 = F(-1, 126)

# k=7 research inputs not yet symbolically certified in this repo.

M7_STAR = F(25866469, 1352400)


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
    print("SCENARIO ONLY: supplied fractions are not certified research inputs.")
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

    assert m7 == M7

    print("pair/four layer =", pair4, "=", float(pair4))
    print("scenario m7 =", m7, "=", float(m7))
    print("exact Stieltjes floor =", M7_STAR, "=", float(M7_STAR))
    print("gap =", m7 - M7_STAR, "=", float(m7 - M7_STAR))

    # The source outlook's historical lower endpoint 19.123 is
    # already below the exact floor implied by pinned m0..m6.
    assert F(19123, 1000) < M7_STAR
    print("NOTE: source reconnaissance lower endpoint 19.123 is stale")
    print("SCENARIO M7 ARITHMETIC CHECKS PASSED")


if __name__ == "__main__":
    main()
