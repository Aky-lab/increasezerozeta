#!/usr/bin/env python3
"""Exact moment targets derived from the shared audited continuum inputs."""
from fractions import Fraction as F
from model_moments import MODEL_MOMENTS, MODEL_M7
from christoffel_exact import solve_linear, hankel, christoffel_at_zero, alternating_negative_halfline

H3 = hankel(MODEL_MOMENTS,3)
INVERSE = [solve_linear(H3,[F(i==j) for i in range(4)]) for j in range(4)]
SHIFTED = [[MODEL_MOMENTS[i+j+1] for j in range(3)] for i in range(3)]
TAIL = [MODEL_MOMENTS[i+4] for i in range(3)]
M7_STAR = sum(x*y for x,y in zip(TAIL,solve_linear(SHIFTED,TAIL)))
LAMBDA3 = 1/INVERSE[0][0]
C_GAP = LAMBDA3*INVERSE[0][3]**2
assert M7_STAR==F(1031677,54096) and LAMBDA3==F(247,2519)
assert C_GAP==F(3732624,622193)


def m8_floor(m7):
    b = TAIL+[m7]
    return sum(x*y for x,y in zip(b,solve_linear(H3,b)))


def lambda4(m7,m8):
    b = TAIL+[m7]
    delta = m8-m8_floor(m7)
    coupling = sum(x*y for x,y in zip(INVERSE[0],b))
    return delta/(INVERSE[0][0]*delta+coupling*coupling)


def lambda4_gap(m7,m8):
    gap,delta = m7-M7_STAR,m8-m8_floor(m7)
    return LAMBDA3*delta/(delta+C_GAP*gap*gap)


def delta_budget(simple_target):
    L = (F(1)-simple_target)/2
    if not F(0)<L<LAMBDA3:
        raise ValueError("target must improve on the degree-three bound")
    return L*C_GAP/(LAMBDA3-L)


def m8_target(m7,simple_target):
    if not F(0)<simple_target<F(1):
        raise ValueError("target must be between zero and one")
    if simple_target<=1-2*LAMBDA3:
        return None  # Every admissible positive-definite extension meets it.
    return m8_floor(m7)+delta_budget(simple_target)*(m7-M7_STAR)**2


def main():
    print("CONTINUUM MODEL TARGETS: analytic transport remains open.")
    print("degree-three simple-zero conversion =",1-2*LAMBDA3,"=",float(1-2*LAMBDA3))
    print("m7 Stieltjes floor =",M7_STAR,"=",float(M7_STAR))
    for label,target in (("80%",F(4,5)),("80.2%",F(401,500)),("81%",F(81,100)),("90%",F(9,10))):
        cap = m8_target(MODEL_M7,target)
        if cap is None:
            print(label+": already reached by the degree-three model conversion")
            continue
        expected = (1-target)/2
        assert lambda4(MODEL_M7,cap)==lambda4_gap(MODEL_M7,cap)==expected
        extended = dict(MODEL_MOMENTS)
        extended[7],extended[8] = MODEL_M7,cap
        direct,q,determinant = christoffel_at_zero(extended,4)
        assert direct==expected and determinant>0 and alternating_negative_halfline(q)
        print(label+": m8 <=",cap,"=",float(cap))
    for m7 in (MODEL_M7,M7_STAR+F(1,100),M7_STAR+F(1,50)):
        m8 = m8_floor(m7)+F(1,100)
        assert lambda4(m7,m8)==lambda4_gap(m7,m8)
        extended = dict(MODEL_MOMENTS)
        extended[7],extended[8] = m7,m8
        direct,_,determinant = christoffel_at_zero(extended,4)
        assert direct==lambda4(m7,m8) and determinant>0
    print("ALL EXACT TARGET AND DIRECT HANKEL CHECKS PASSED")


if __name__=="__main__":
    main()
