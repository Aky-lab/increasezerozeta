# Complete eighth-order continuum-model certificate

The finite model uses the cumulant/prefix-walk conventions in the
[research overview](../RESEARCH.md). Every class required by the Bell(8)
ledger has an exact value or an exact local reduction.

## New eighth-order classes

| Class | Exact value | Certificate |
|---|---:|---|
| {2^4} | 1661/3780 | [Pair-cycle lift](paired_cycle_flow_polytopes.md) |
| {4,2,2} | -127/840 | [Mixed-cycle lift](mixed_cycle_flow_polytopes.md) |
| {4,4} | 23/4536 | [Mixed-cycle lift](mixed_cycle_flow_polytopes.md) |
| {6,2} | -563/11340 | [Mixed-cycle lift](mixed_cycle_flow_polytopes.md) |
| C8 | 157/4032 | Pure-cycle lift below |

Their sum is A8=633/2240. Inherited classes are handled by
[singleton deletion and three-block vanishing](model_local_identities.md),
with the audited lower-order values. The resulting model moments are

    m7 = 3439/180,
    m8 = 3311/90 + A8 = 747361/20160.

## Pure eighth-order cycle

The [pure-cycle network proof](pure_cycle_flow_polytopes.md) applies
at b=8 without a change to its dimension or lattice normalization.
Each of the 94,586 signed inner cyclic-partition terms is an integral
bounded circulation polytope of dimension nine. Thus S8(n) is a
period-one polynomial of degree at most nine, with leading coefficient C8.

The counts n=1,...,11 are

    0, 12, 666, 9656, 73972, 386136, 1555050, 5188584,
    15005682, 38780544, 91520352.

The first ten determine the polynomial. The last is held out and
matches its prediction. The ninth forward difference divided by 9!
is 157/4032. The [integer record](../results/pure_cycle_eight_2026-10-05.json)
includes every count, integer bounds and pinned source hashes.

`scripts/pure_cycle_eight.py` uses bounded chunks of outer positions
and cached sorted frequency tuples. At n<=12 the absolute sum is
bounded by 94,586*n^9<2^63, and every point's cumulant fits int32.
Products and final sums use int64; interpolation uses exact Fractions.
The existing lower-order evaluator and its certificate source bytes
are preserved.

The [independent record](../results/pure_cycle_eight_independent_2026-10-05.json)
checks all 94,586 cyclic partitions using a restricted-growth-string
compiler. It also sums the entire bounded outer cube at n=2 and n=3,
using all frequencies and direct subset sums. It does not use the
primary compiler's coefficient forms, sorted-frequency cache or outer
coordinate eliminator to compute these counts. They recover 12 and 666.

## Exact consumption and positivity

For the complete model sequence, the fourth-degree Hankel matrix has

    det(H4) = 2448223/46088663040 > 0,
    m8 - m8_floor(m7) = 2448223/104569920 > 0.

Its H3 block is positive definite and m7 exceeds the shifted Stieltjes
floor, so the ordinary and shifted matrices satisfy the required
positivity constraints. Direct rational inversion gives

    lambda4(0) = 12241115/162540559,
    1-2*lambda4(0) = 138058329/162540559 = 0.8493777174... .

The minimizing polynomial has coefficients, in increasing powers,

    (1, -851656632/162540559, 1278096456/162540559,
        -729534540/162540559, 140399280/162540559).

The signs alternate, hence q4(-t)>=1 for t>=0. Its square therefore
controls mass on the nonpositive half-line as well as at zero.

The approximately 84.938% number is a **conditional counting conversion**
for this supplied model. It is not an established zeta-zero bound.
Proving the arithmetic moment identities and spectral/counting interface,
and independently reviewing the finite certificates, remain necessary.

## Reproduction

```sh
python scripts/verify_project.py
python scripts/m8_ledger.py
```

These standard-library checks verify the exact algebra, recorded counts,
source hashes and independent scalar mixed-class sums. To repeat the
larger NumPy enumerations:

```sh
python scripts/mixed_cycle_lattice.py --out results/mixed_cycle_reproduction.json
python scripts/pure_cycle_eight.py --out results/cycle_eight_reproduction.json
python scripts/verify_cycle_eight_lattice.py --out results/cycle_eight_checks.json
```

The full evaluators resume their own checkpoints only when source hashes
match. The independent checker compares against the committed C8 record.
