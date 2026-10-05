#!/usr/bin/env python3
"""Bell(7) continuum ledger using exact finite model inputs.

Defaults: {5,2}=1/8, C7=-17/360 and audited {2,2,2}=32/105.
The frozen-singleton identities and arithmetic transport remain analytic
obligations. Explicit overrides calculate alternative model scenarios.
"""

import argparse
from fractions import Fraction as F
from model_moments import T222, M7_BASE, MODEL_J52, MODEL_C7, MODEL_M7
from k8_target import M7_STAR

T_ADJ = F(7, 60)
T_OPP = F(1, 30)
PHI4 = F(-1, 60)
C5 = F(1, 36)
J42 = F(-23, 420)
C6 = F(-1, 126)

# Exact finite model evaluations; see the lattice and flow-polytope notes.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--j52", type=F, default=MODEL_J52, help="joint model input; default 1/8")
    parser.add_argument("--c7", type=F, default=MODEL_C7, help="pure-cycle input; default -17/360")
    args = parser.parse_args()
    J52, C7 = args.j52, args.c7
    M7 = M7_BASE + J52 + C7
    print("CONTINUUM MODEL LEDGER: arithmetic transport remains open.")
    if (J52, C7) == (MODEL_J52, MODEL_C7):
        assert M7 == MODEL_M7
    else:
        print("Using explicit alternative model inputs.")
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
    print("model m7 =", m7, "=", float(m7))
    print("exact Stieltjes floor =", M7_STAR, "=", float(M7_STAR))
    print("gap =", m7 - M7_STAR, "=", float(m7 - M7_STAR))

    print("MODEL M7 ARITHMETIC CHECKS PASSED")


if __name__ == "__main__":
    main()
