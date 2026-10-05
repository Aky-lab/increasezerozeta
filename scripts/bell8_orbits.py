#!/usr/bin/env python3
"""Exact dihedral-orbit gates for the new Bell(8) joint classes."""

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
    return tuple(sorted((len(b) for b in part), reverse=True))


def canon(part):
    return tuple(sorted(tuple(sorted(b)) for b in part))


def transform(part, n, shift, reflect):
    out = []
    for block in part:
        b = []
        for x in block:
            y = -x if reflect else x
            b.append((y + shift) % n)
        out.append(tuple(sorted(b)))
    return tuple(sorted(out))


def orbits(n, sig):
    parts = {canon(p) for p in set_partitions(list(range(n)))
             if signature(p) == sig}
    unseen = set(parts)
    out = []
    while unseen:
        p = min(unseen)
        orb = set()
        for reflect in (False, True):
            for shift in range(n):
                orb.add(transform(p, n, shift, reflect))
        orb &= parts
        out.append((p, len(orb)))
        unseen -= orb
    return sorted(out)


EXPECTED = {
    (6, 2): (28, 4, [4, 8, 8, 8]),
    (4, 4): (35, 7, [1, 2, 4, 4, 8, 8, 8]),
    (4, 2, 2): (
        210, 22,
        [2, 4, 4, 4, 4,
         8, 8, 8, 8, 8, 8, 8, 8, 8,
         16, 16, 16, 16, 16, 16, 16, 16],
    ),
    (2, 2, 2, 2): (
        105, 17,
        [1, 2, 2, 4, 4, 4, 4, 4,
         8, 8, 8, 8, 8, 8, 8, 8, 16],
    ),
}


def main():
    for sig, (total, norb, sizes) in EXPECTED.items():
        oo = orbits(8, sig)
        got_sizes = sorted(s for _, s in oo)
        assert sum(got_sizes) == total
        assert len(oo) == norb
        assert got_sizes == sizes

        print(f"{sig}: {total} placements -> {norb} dihedral orbits")
        print("  orbit sizes:", got_sizes)
        for i, (rep, size) in enumerate(oo):
            print(f"  {i:2d}: size={size:2d} rep={rep}")

    print("ALL BELL(8) DIHEDRAL ORBIT GATES PASSED")


if __name__ == "__main__":
    main()
