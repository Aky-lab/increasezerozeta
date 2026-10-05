#!/usr/bin/env python3
"""Independent finite rational checks of the exact-v spectator reduction.

The integral is computed directly from the piecewise-linear overlap,
then compared to the proposed cubic. This checks a reduction used by
the frontier branch, not the remaining C5 integral or a zeta theorem.
"""

from fractions import Fraction as F
from itertools import product


def overlap(prefixes, shifted, v):
    walk = list(prefixes) + [v + x for x in shifted]
    return max(F(0), F(1) - max(walk) + min(walk))


def direct_integral(prefixes, shifted):
    # All changes of the maximum/minimum occur at a prefix difference.
    # Include every possible zero of 1 - (x - y), and split |v| at 0.
    knots = {F(0)}
    for x in prefixes:
        for y in shifted:
            knots.update((x - y, x - y - 1, x - y + 1))
    knots = sorted(knots)
    total = F(0)
    for lo, hi in zip(knots, knots[1:]):
        left = overlap(prefixes, shifted, lo)
        right = overlap(prefixes, shifted, hi)
        slope = (right - left) / (hi - lo)
        intercept = left - slope * lo
        sign = 1 if lo >= 0 else -1
        total += sign * (
            slope * (hi**3 - lo**3) / 3
            + intercept * (hi**2 - lo**2) / 2
        )
    return total


def reduced_integral(prefixes, shifted):
    low, high = min(prefixes), max(prefixes)
    span = high - low
    if span >= 1:
        return F(0)
    alpha = min(shifted) - low
    beta = high - max(shifted)
    return (1 - span) / 6 * (
        2 * span**2 - 3 * span * (alpha + beta) - 4 * span
        + 3 * alpha**2 + 3 * alpha + 3 * beta**2 + 3 * beta + 2
    )


def main():
    grid = [F(i, 4) for i in range(-3, 4)]
    checked = nonzero = 0
    for suffix in product(grid, repeat=4):
        prefixes = (F(0),) + suffix
        for distance in (1, 2, 3):
            shifted = prefixes[:distance]
            direct = direct_integral(prefixes, shifted)
            reduced = reduced_integral(prefixes, shifted)
            if direct != reduced:
                raise AssertionError((prefixes, distance, direct, reduced))
            checked += 1
            nonzero += int(direct != 0)
    print(f"exact rational spectator cases checked: {checked}")
    print(f"cases with nonzero integral: {nonzero}")
    print("ALL FINITE SPECTATOR REDUCTION CHECKS PASSED")
    print("Scope: finite checks of J_d only; no certification of {5,2}=7/72.")


if __name__ == "__main__":
    main()
