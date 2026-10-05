#!/usr/bin/env python3
"""Finite rational gates against actual walks, not a shared prefix assumption."""

from fractions import Fraction as F
from itertools import product

from spectator_reduction import spectator_integral


def walk_sets(prefixes, distance):
    # Build W_d independently from its increments v,c1,...,-v,...,c5.
    frequencies = [prefixes[i + 1] - prefixes[i] for i in range(4)]
    frequencies.append(-prefixes[-1])
    increments = ([(1, F(0))]
                  + [(0, x) for x in frequencies[:distance - 1]]
                  + [(-1, F(0))]
                  + [(0, x) for x in frequencies[distance - 1:]])
    positions = [(0, F(0))]
    for coefficient, constant in increments[:-1]:
        previous = positions[-1]
        positions.append((previous[0] + coefficient, previous[1] + constant))
    if len(positions) != 7 or any(a not in (0, 1) for a, _ in positions):
        raise AssertionError(positions)
    return (tuple(b for a, b in positions if a == 0),
            tuple(b for a, b in positions if a == 1))


def overlap(fixed, shifted, v):
    walk = list(fixed) + [v + x for x in shifted]
    return max(F(0), F(1) - max(walk) + min(walk))


def direct_integral(fixed, shifted):
    # Split at every possible max/min switch and overlap zero.
    knots = {F(-1), F(0), F(1)}
    for x in fixed:
        for y in shifted:
            knots.update((x - y, x - y - 1, x - y + 1))
    knots = sorted(knots)
    total = F(0)
    for lo, hi in zip(knots, knots[1:]):
        left = overlap(fixed, shifted, lo)
        right = overlap(fixed, shifted, hi)
        slope = (right - left) / (hi - lo)
        intercept = left - slope * lo
        sign = 1 if lo >= 0 else -1
        # Integral |v|*(slope*v+intercept), computed independently.
        total += sign * (
            slope * (hi**3 - lo**3) / 3
            + intercept * (hi**2 - lo**2) / 2
        )
    return total


def historical_cubic(prefixes, distance):
    low, high = min(prefixes), max(prefixes)
    span = high - low
    if span >= 1:
        return F(0)
    shifted = prefixes[:distance]
    alpha = min(shifted) - low
    beta = high - max(shifted)
    return (1 - span) / 6 * (
        2 * span**2 - 3 * span * (alpha + beta) - 4 * span
        + 3 * alpha**2 + 3 * alpha + 3 * beta**2 + 3 * beta + 2
    )


def main():
    counterexample = (F(0), F(4, 5), F(1, 5), F(3, 10), F(2, 5))
    fixed, shifted = walk_sets(counterexample, 3)
    if direct_integral(fixed, shifted) != F(2, 75):
        raise AssertionError("counterexample integration")
    if historical_cubic(counterexample, 3) != F(1, 375):
        raise AssertionError("historical regression")
    if spectator_integral(counterexample, 3) != F(2, 75):
        raise AssertionError("corrected regression")
    grid = [F(i, 4) for i in range(-3, 4)]
    checked = nonzero = historical_mismatches = 0
    for suffix in product(grid, repeat=4):
        prefixes = (F(0),) + suffix
        for distance in (1, 2, 3):
            fixed, shifted = walk_sets(prefixes, distance)
            direct = direct_integral(fixed, shifted)
            reduced = spectator_integral(prefixes, distance)
            if direct != reduced:
                raise AssertionError((prefixes, distance, direct, reduced))
            if historical_cubic(prefixes, distance) != direct:
                if distance != 3:
                    raise AssertionError("unexpected historical failure")
                historical_mismatches += 1
            checked += 1
            nonzero += int(direct != 0)
    if not historical_mismatches:
        raise AssertionError("historical bug was not detected")
    print(f"actual-walk rational cases checked: {checked}")
    print(f"cases with nonzero integral: {nonzero}")
    print(f"historical distance-three mismatches: {historical_mismatches}")
    print("counterexample: actual 2/75; historical 1/375")
    print("ALL FINITE ACTUAL-WALK SPECTATOR CHECKS PASSED")
    print("Scope: J_d only; no certification of the remaining C5 integral.")


if __name__ == "__main__":
    main()
