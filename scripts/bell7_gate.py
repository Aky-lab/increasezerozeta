#!/usr/bin/env python3
"""Exact Bell(7) bookkeeping gate for the k=8 research track.

Pure combinatorics:
- Bell(7)=877 and block-size class counts;
- the {5,2} class has 21 placements;
- those placements split 7/7/7 across cyclic distances 1,2,3.

This is the combinatorial specification for the new seventh-moment
joint constant {5,2}.
"""

from collections import Counter


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def signature(part):
    return tuple(sorted((len(block) for block in part), reverse=True))


def cyclic_distance(pair, n=7):
    a, b = sorted(pair)
    d = b - a
    return min(d, n - d)


def main():
    parts = list(set_partitions(list(range(7))))
    counts = Counter(signature(p) for p in parts)

    expected = {
        (1, 1, 1, 1, 1, 1, 1): 1,
        (2, 1, 1, 1, 1, 1): 21,
        (2, 2, 1, 1, 1): 105,
        (2, 2, 2, 1): 105,
        (3, 1, 1, 1, 1): 35,
        (3, 2, 1, 1): 210,
        (3, 2, 2): 105,
        (3, 3, 1): 70,
        (4, 1, 1, 1): 35,
        (4, 2, 1): 105,
        (4, 3): 35,
        (5, 1, 1): 21,
        (5, 2): 21,
        (6, 1): 7,
        (7,): 1,
    }

    assert len(parts) == 877
    assert counts == expected

    dists = Counter()
    for p in parts:
        if signature(p) != (5, 2):
            continue
        pair = next(block for block in p if len(block) == 2)
        dists[cyclic_distance(pair)] += 1

    assert dists == {1: 7, 2: 7, 3: 7}

    print("Bell(7) =", len(parts))
    print("block-size classes:")
    for sig in sorted(counts, key=lambda s: (len(s), s), reverse=True):
        print(" ", sig, counts[sig])
    print("{5,2} cyclic-distance multiplicities:", dict(sorted(dists.items())))
    print("ALL BELL(7) CHECKS PASSED")


if __name__ == "__main__":
    main()
