#!/usr/bin/env python3
"""Exact spectator integration for the actual three seven-cycle walks.

This integrates the spectator variable only. It does not evaluate the
remaining C5 integral, identify {5,2}, or establish arithmetic transport.
"""

from fractions import Fraction as F


def prefix_sets(prefixes, distance):
    if len(prefixes) != 5 or prefixes[0] != 0:
        raise ValueError("expected (0,s1,s2,s3,s4)")
    if distance not in (1, 2, 3):
        raise ValueError("distance must be 1, 2, or 3")
    fixed = (prefixes[0],) + tuple(prefixes[distance - 1:])
    shifted = tuple(prefixes[:distance])
    return fixed, shifted


def spectator_integral(prefixes, distance):
    """Return integral C2(v)*O7(W_distance(v)) dv as a Fraction."""
    fixed, shifted = prefix_sets(tuple(F(x) for x in prefixes), distance)
    m, M = min(fixed), max(fixed)
    a, A = min(shifted), max(shifted)
    if max(M - m, A - a) >= 1:
        return F(0)
    lo, hi = M - a - 1, 1 + m - A
    b, c = m - a, M - A
    # O'' = delta_lo - delta_b - delta_c + delta_hi.
    # Since (|v|^3/6)'' = |v|, integration by parts gives J below.
    # Both prefix sets contain 0, so support is within [-1,1],
    # where C2(v)=min(|v|,1)=|v|.
    return (abs(lo)**3 + abs(hi)**3 - abs(b)**3 - abs(c)**3) / 6
