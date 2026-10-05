# Bell(8) ledger and the 80% target

The combinatorial ledger is conditional on the frozen-singleton and 3-block-vanishing identities in the lower-order model. Its analytic transport remains open.

## New eighth-order classes

Define

    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

The inherited contributions give

    m8 = 167/6 + 28{2,2,2} + 28{4,2} + 56C5
         + 8{5,2} + 28C6 + 8C7 + A8.

Using {2,2,2}=131/420, {4,2}=-23/420, C5=1/36 and C6=-1/126 reduces this to

    m8 = 1091/30 + 8{5,2} + 8C7 + A8.

The corrected exact model input {5,2}=1/8 therefore gives

    m8 = 1121/30 + 8C7 + A8.

The [exact flow-polytope value](pure_cycle_flow_polytopes.md) C7=-17/360 gives the model ledger m7=862/45 and m8=3329/90+A8. Arithmetic transport remains open.

## Candidate class values

| Class | Numerical candidate |
|---|---:|
| {2^4} | 1661/3780 |
| {4,2,2} | -127/840 |
| {6,2} | -1/20 |

These values require exact computation or rigorous enclosure. Their sum is 1801/7560. Thus, under those additional assumptions,

    A8 = 1801/7560 + T8,    T8 = {4,4} + C8.

## Target geometry

The degree-four Christoffel bound reaches a conditional simple-zero proportion of 80% when lambda_4(0)<=1/10. The permitted m8 interval depends on m7 and the lower moments.

[The generic target derivation](k8_target_geometry.md) and `scripts/k8_target.py` give the Hankel floor and exact target curves. The calculator `scripts/m8_ledger.py` uses the exact seventh-order model inputs by default and accepts explicit overrides.

At m7=862/45 and m8=3329/90+A8, the moment-cone floor and conditional
80% target give

    204409997/748818000 <= A8 <= 33596531/68531400.

This is approximately 0.272977 <= A8 <= 0.490236. The interval is a
constraint on supplied model moments; it does not certify any of the
five eighth-order classes or their arithmetic transport.

Earlier thresholds based on {5,2}=7/72 use an erroneous distance-three prefix set and do not apply to the corrected model. Eighth-order class evaluations and arithmetic transport are needed before a zeta bound follows.
