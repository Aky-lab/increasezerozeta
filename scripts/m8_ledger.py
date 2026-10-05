#!/usr/bin/env python3
"""Exact Bell(8) continuum-model ledger and aggregate target calculator."""
import argparse
from fractions import Fraction as F
from model_moments import T222, M7_BASE, MODEL_J52, MODEL_C7, MODEL_M7, MODEL_M8_INHERITED, PAIR4
from k8_target import M7_STAR, LAMBDA3, m8_floor, m8_target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--j52",type=F,default=MODEL_J52)
    parser.add_argument("--c7",type=F,default=MODEL_C7)
    args = parser.parse_args()
    m7 = M7_BASE+args.j52+args.c7
    if m7<M7_STAR:
        parser.error("scenario m7 is below the model Stieltjes floor")
    base = F(1)+F(28,3)+140*F(7,60)+70*F(1,30)+70*F(-1,60)
    assert base==F(167,6)
    inherited = base+28*T222+28*F(-23,420)+56*F(1,36)+8*args.j52+28*F(-1,126)+8*args.c7
    if (args.j52,args.c7)==(MODEL_J52,MODEL_C7):
        assert m7==MODEL_M7 and inherited==MODEL_M8_INHERITED
    print("CONTINUUM MODEL TARGETS: arithmetic transport remains open.")
    print("model m7 =",m7,"=",float(m7))
    print("m8 =",inherited,"+ A8")
    print("A8 := {2^4}+{4,2,2}+{4,4}+{6,2}+C8")
    print("exact {2^4} =",PAIR4)
    floor = m8_floor(m7)
    print("m8 moment-cone floor =",floor,"=",float(floor))
    print("A8 floor =",floor-inherited,"=",float(floor-inherited))
    for label,target in (("80%",F(4,5)),("80.2%",F(401,500)),("81%",F(81,100)),("90%",F(9,10))):
        cap = m8_target(m7,target)
        if cap is None:
            print(label+": already reached by the degree-three model conversion")
        else:
            print(label+": A8 <=",cap-inherited,"=",float(cap-inherited))
    print("MODEL BELL(8) ARITHMETIC CHECKS PASSED")


if __name__=="__main__":
    main()
