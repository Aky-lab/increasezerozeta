#!/usr/bin/env python3
"""Four-dimensional midpoint audit of corrected and historical spectator weights.

Numerical reconnaissance only: neither values nor Richardson probes are
rigorous enclosures. NumPy is required for this optional research script.
"""

import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import time

import numpy as np

from spectator_reduction import spectator_integral


def partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def cumulant_terms():
    terms = []
    for part in partitions(list(range(5))):
        sign = (-1) ** (len(part) - 1)
        first = next(i for i, block in enumerate(part) if 0 in block)
        others = [i for i in range(len(part)) if i != first]
        for order in itertools.permutations(others):
            current, prefixes = set(), []
            for index in (first,) + order[:-1]:
                # For a one-block term there are no partial sums.
                if len(part) == 1:
                    break
                current.update(part[index])
                prefixes.append(tuple(
                    int(j in current) - int(4 in current) for j in range(4)
                ))
            terms.append((sign, tuple(prefixes)))
    if len(terms) != 150:
        raise AssertionError("C5 term count")
    return terms


TERMS = cumulant_terms()


def cumulant_value(frequencies):
    zero = np.zeros_like(frequencies[0])
    forms = {}
    for _, prefixes in TERMS:
        for coeff in prefixes:
            if coeff not in forms:
                value = np.zeros_like(zero)
                for a, x in zip(coeff, frequencies):
                    if a:
                        value = value + a * x
                forms[coeff] = value
    result = np.zeros_like(zero)
    for sign, prefixes in TERMS:
        positions = [zero] + [forms[c] for c in prefixes]
        span = np.maximum.reduce(positions) - np.minimum.reduce(positions)
        result += sign * np.maximum(1 - span, 0)
    return result


def weight(fixed, shifted):
    m, M = np.minimum.reduce(fixed), np.maximum.reduce(fixed)
    a, A = np.minimum.reduce(shifted), np.maximum.reduce(shifted)
    lo, hi = M - a - 1, 1 + m - A
    b, c = m - a, M - A
    value = (np.abs(lo)**3 + np.abs(hi)**3
             - np.abs(b)**3 - np.abs(c)**3) / 6
    return np.where(np.maximum(M - m, A - a) < 1, value, 0)


def calibration():
    rng = np.random.default_rng(5203)
    samples = rng.integers(-8, 9, size=(100, 4))
    for row in samples:
        prefixes = (F(0),) + tuple(F(int(i), 10) for i in row)
        arrays = [np.array([float(x)]) for x in prefixes]
        for d in (1, 2, 3):
            fixed = [arrays[0]] + arrays[d - 1:]
            numeric = float(weight(fixed, arrays[:d])[0])
            exact = float(spectator_integral(prefixes, d))
            if abs(numeric - exact) > 1e-12:
                raise AssertionError((prefixes, d, numeric, exact))


def run(step):
    ratio = F(2) / step
    if ratio.denominator != 1:
        raise ValueError("step must divide the interval length 2 exactly")
    count = ratio.numerator
    if not 2 <= count <= 80:
        raise ValueError("choose 2 <= 2/step <= 80")
    grid = -1 + (np.arange(count) + 0.5) * float(step)
    c2, c3, c4 = np.meshgrid(grid, grid, grid, indexing="ij")
    zero = np.zeros_like(c2)
    actual, historical, pure = np.zeros(3), np.zeros(3), 0.0
    started = time.monotonic()
    for value in grid:
        c1 = np.full_like(c2, value)
        s2, s3, s4 = c1 + c2, c1 + c2 + c3, c1 + c2 + c3 + c4
        prefixes = [zero, c1, s2, s3, s4]
        c5 = cumulant_value([c1, c2, c3, c4])
        pure_overlap = np.maximum(
            1 - np.maximum.reduce(prefixes) + np.minimum.reduce(prefixes), 0
        )
        pure += float(np.sum(pure_overlap * c5))
        for d in (1, 2, 3):
            shifted = prefixes[:d]
            fixed = [zero] + prefixes[d - 1:]
            actual[d - 1] += float(np.sum(weight(fixed, shifted) * c5))
            historical[d - 1] += float(np.sum(weight(prefixes, shifted) * c5))
    scale = float(step)**4
    actual, historical = actual * scale, historical * scale
    return {
        "step": str(step), "grid_points_per_axis": count,
        "corrected_U": actual.tolist(),
        "historical_U": historical.tolist(),
        "corrected_joint": 7 * float(actual.sum()),
        "historical_joint": 7 * float(historical.sum()),
        "pure_C5_anchor": pure * scale,
        "anchor_target": "1/36",
        "wall_seconds": time.monotonic() - started,
        "status": "numerical reconnaissance; not a rigorous enclosure",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("steps", nargs="+", type=F)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    calibration()
    results = []
    for step in args.steps:
        result = run(step)
        results.append(result)
        print(json.dumps(result), flush=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({
            "scope": "corrected actual walks; outer support in [-1,1]^4",
            "float_weight_vs_rational_cases": 300,
            "numpy_version": np.__version__,
            "results": results,
        }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
