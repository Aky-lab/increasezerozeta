#!/usr/bin/env python3
"""Exact Bell(8) continuum-model ledger and aggregate target calculator."""
import argparse
from fractions import Fraction as F
from model_moments import (T222, M7_BASE, MODEL_J52, MODEL_C7, MODEL_M7,
                          MODEL_M8_INHERITED, PAIR4, MODEL_J422, MODEL_J44,
                          MODEL_J62, MODEL_C8, MODEL_A8, MODEL_M8, MODEL_MOMENTS)
from k8_target import M7_STAR, LAMBDA3, m8_floor, m8_target
from christoffel_exact import christoffel_at_zero, alternating_negative_halfline


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
    print("exact {4,2,2} =",MODEL_J422,"; {4,4} =",MODEL_J44)
    print("exact {6,2} =",MODEL_J62,"; C8 =",MODEL_C8)
    print("exact A8 =",MODEL_A8)
    floor = m8_floor(m7)
    print("m8 moment-cone floor =",floor,"=",float(floor))
    print("A8 floor =",floor-inherited,"=",float(floor-inherited))
    m8 = inherited+MODEL_A8
    if m8<=floor:
        parser.error("scenario m8 is not inside the positive-definite moment cone")
    extended = dict(MODEL_MOMENTS)
    extended[7],extended[8] = m7,m8
    lam,q,determinant = christoffel_at_zero(extended,4)
    assert determinant>0 and alternating_negative_halfline(q)
    if (args.j52,args.c7)==(MODEL_J52,MODEL_C7):
        assert m8==MODEL_M8==F(747361,20160)
        assert lam==F(12241115,162540559)
        assert 1-2*lam==F(138058329,162540559)
    print("model m8 =",m8,"=",float(m8))
    print("degree-four lambda(0) =",lam,"=",float(lam))
    print("conditional simple-zero conversion =",1-2*lam,"=",float(1-2*lam))
    for label,target in (("80%",F(4,5)),("80.2%",F(401,500)),("81%",F(81,100)),("90%",F(9,10))):
        cap = m8_target(m7,target)
        if cap is None:
            print(label+": already reached by the degree-three model conversion")
        else:
            print(label+": A8 <=",cap-inherited,"=",float(cap-inherited))
    print("MODEL BELL(8) ARITHMETIC CHECKS PASSED")


if __name__=="__main__":
    main()
