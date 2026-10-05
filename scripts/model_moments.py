"""Certified finite continuum-model moments, with analytic transport open."""
from fractions import Fraction as F

# Direct perfect-matching model; see the pairing audit and exact certificates.
T222 = F(32,105)
PAIR4 = F(1661,3780)
MODEL_J52 = F(1,8)
MODEL_C7 = F(-17,360)
MODEL_J42 = F(-23,420)
MODEL_J62 = F(-563,11340)
MODEL_J422 = F(-127,840)
MODEL_J44 = F(23,4536)
MODEL_C8 = F(157,4032)
PAIR1 = F(1,3)
PAIR2 = F(4,15)
C4, C5, C6 = F(-1,60), F(1,36), F(-1,126)

# Singleton deletion and three-block vanishing are proved in the finite
# model; see notes/model_local_identities.md. All class inputs are certified.
MODEL_MOMENTS = {0:F(1),1:F(1),2:1+PAIR1,3:1+3*PAIR1,
                 4:1+6*PAIR1+PAIR2+C4,
                 5:1+10*PAIR1+5*(PAIR2+C4)+C5,
                 6:1+15*PAIR1+15*(PAIR2+C4)+T222+MODEL_J42+6*C5+C6}
assert MODEL_MOMENTS[6]==F(640,63)

M7_BASE = F(67,4)+21*C5+7*T222+7*MODEL_J42+7*C6
M8_BASE = F(167,6)+28*T222+28*MODEL_J42+56*C5+28*C6
MODEL_M7 = M7_BASE+MODEL_J52+MODEL_C7
MODEL_M8_INHERITED = M8_BASE+8*MODEL_J52+8*MODEL_C7
MODEL_A8 = PAIR4+MODEL_J422+MODEL_J44+MODEL_J62+MODEL_C8
MODEL_M8 = MODEL_M8_INHERITED+MODEL_A8
MODEL_MOMENTS[7],MODEL_MOMENTS[8] = MODEL_M7,MODEL_M8
assert M7_BASE==F(685,36) and MODEL_M7==F(3439,180)
assert MODEL_M8_INHERITED==F(3311,90)
