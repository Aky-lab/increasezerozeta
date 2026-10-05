"""Shared continuum-model inputs, with analytic transport still open."""
from fractions import Fraction as F

# Direct perfect-matching model; see the pairing audit and exact certificates.
T222 = F(32,105)
PAIR4 = F(1661,3780)
MODEL_J52 = F(1,8)
MODEL_C7 = F(-17,360)

# The reference sixth-order ledger's other class contributions are retained.
REFERENCE_T222 = F(131,420)
REFERENCE_M6 = F(12809,1260)
M6_OTHER = REFERENCE_M6-REFERENCE_T222
MODEL_MOMENTS = {0:F(1),1:F(1),2:F(4,3),3:F(2),4:F(13,4),5:F(101,18),
                 6:M6_OTHER+T222}
assert MODEL_MOMENTS[6]==F(640,63)

M7_BASE = F(67,4)+21*F(1,36)+7*T222+7*F(-23,420)+7*F(-1,126)
M8_BASE = F(167,6)+28*T222+28*F(-23,420)+56*F(1,36)+28*F(-1,126)
MODEL_M7 = M7_BASE+MODEL_J52+MODEL_C7
MODEL_M8_INHERITED = M8_BASE+8*MODEL_J52+8*MODEL_C7
assert M7_BASE==F(685,36) and MODEL_M7==F(3439,180)
assert MODEL_M8_INHERITED==F(3311,90)
