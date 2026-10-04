#!/usr/bin/env python3
"""Exact-arithmetic checks for notes/class_subtracted_universality.md.

This is not a numerical proof of the surrounding zeta-function argument.
It checks the finite algebra used in the local universality lemma and the
rational arithmetic in the continuum-core consumption.
"""

from fractions import Fraction
from itertools import combinations


def primes_upto(n):
    out = []
    for x in range(2, n + 1):
        if all(x % p for p in range(2, int(x ** 0.5) + 1)):
            out.append(x)
    return out


def l1(v):
    return sum(abs(x) for x in v)


def generic_vector(p):
    a = Fraction(1, p - 1)
    return (1 - a * a, a * a)


def special_vector(p):
    a = Fraction(1, p - 1)
    return (1 + a, -a)


def tensor(v, w):
    return tuple(x * y for x in v for y in w)


def sub(v, w):
    return tuple(x - y for x, y in zip(v, w))


def local_data(p, q):
    alpha = Fraction(1, p - 1)
    beta = Fraction(1, q - 1)

    gp, gq = generic_vector(p), generic_vector(q)
    sp, sq = special_vector(p), special_vector(q)

    # Exact local identities.
    assert sum(gp) == 1 and sum(gq) == 1
    assert sum(sp) == 1 and sum(sq) == 1
    assert l1(sub(sp, gp)) == 2 * alpha * (1 + alpha)
    assert l1(sub(sq, gq)) == 2 * beta * (1 + beta)
    assert l1(sp) == 1 + 2 * alpha
    assert l1(sq) == 1 + 2 * beta

    # Full-tensor comparison, with the common positive Euler tail factored out.
    full_exact = l1(sub(tensor(sp, sq), tensor(gp, gq)))

    # Symmetric telescoping upper bound used in the note.
    full_bound = 2 * (1 + alpha) * (1 + beta) * (alpha + beta)
    assert full_exact <= full_bound

    # The actual class-d replacement subtracts s_p tensor s_q;
    # the (1,1) convention subtracts delta_(0,0).
    delta = (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
    class_exact = l1(sub(tensor(sp, sq), delta))
    class_formula = 2 * alpha + 2 * beta + 4 * alpha * beta
    assert class_exact == class_formula

    total_bound = full_bound + class_bound(p, q)
    coarse_pminus1 = 8 * (alpha + beta)
    coarse_p = 12 * (Fraction(1, p) + Fraction(1, q))

    assert total_bound <= coarse_pminus1
    assert coarse_pminus1 <= coarse_p

    return {
        "full_exact": full_exact,
        "full_bound": full_bound,
        "class_exact": class_exact,
        "total_bound": total_bound,
        "coarse_pminus1": coarse_pminus1,
        "coarse_p": coarse_p,
    }


def class_bound(p, q):
    alpha = Fraction(1, p - 1)
    beta = Fraction(1, q - 1)
    return 2 * alpha + 2 * beta + 4 * alpha * beta


def verify_prime_pairs(limit=101):
    ps = [p for p in primes_upto(limit) if p % 2]
    checked = 0
    worst_ratio = Fraction(0)
    worst_pair = None

    for p, q in combinations(ps, 2):
        data = local_data(p, q)
        ratio = data["total_bound"] / data["coarse_pminus1"]
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_pair = (p, q)
        checked += 1

    return checked, worst_pair, worst_ratio


def verify_core_arithmetic():
    # Integral of (nu-1) * (2-nu)^3 / 6 on [1,2]:
    first = Fraction(1, 120)

    # Integral of (nu-1) * (4/3)*(3/2-nu)^3 on [1,3/2]:
    second = Fraction(1, 480)

    weighted_q = first + second
    assert weighted_q == Fraction(1, 96)

    c_core = -2 * weighted_q
    assert c_core == Fraction(-1, 48)

    # Integral Q(nu) dnu = 1/24 + 1/48 = 1/16.
    q_mass = Fraction(1, 24) + Fraction(1, 48)
    assert q_mass == Fraction(1, 16)

    # Keep the reference candidate's separate 0.0111 sawtooth-tail charge.
    e_tail = Fraction(111, 10_000)
    remainder = c_core + e_tail
    assert remainder == Fraction(-73, 7_500)

    delta = remainder + Fraction(1, 30)
    assert delta == Fraction(59, 2_500)

    lam2 = (
        Fraction(5, 108) + delta / 3
    ) / (
        Fraction(1, 3) + 4 * delta / 3
    )
    assert lam2 == Fraction(457, 3_078)

    simple = 1 - 2 * lam2
    assert simple == Fraction(1082, 1539)

    return q_mass, weighted_q, c_core, remainder, lam2, simple


def main():
    checked, worst_pair, worst_ratio = verify_prime_pairs()
    q_mass, weighted_q, c_core, remainder, lam2, simple = verify_core_arithmetic()

    print(f"prime pairs checked: {checked}")
    print(
        "largest (proved-bound)/(8/(p-1)+8/(q-1)) ratio: "
        f"{float(worst_ratio):.6f} at {worst_pair}"
    )
    print(f"integral Q = {q_mass} = {float(q_mass):.12f}")
    print(
        "integral (nu-1)Q = "
        f"{weighted_q} = {float(weighted_q):.12f}"
    )
    print(f"C_core = {c_core} = {float(c_core):.12f}")
    print(
        "C_core + 0.0111 = "
        f"{remainder} = {float(remainder):.12f}"
    )
    print(f"Lambda_2(0) = {lam2} = {float(lam2):.12f}")
    print(
        "simple-critical proportion (conditional consumption) = "
        f"{simple} = {float(simple):.12%}"
    )
    print("ALL EXACT-ARITHMETIC CHECKS PASSED")


if __name__ == "__main__":
    main()
