# Bell(8) ledger and conditional counting targets

The ledger assumes frozen-singleton and 3-block-vanishing identities.
The other retained reference class inputs and arithmetic transport remain
subject to analytic review.

## Eighth-order classes

Define

    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

The inherited contributions are

    m8 = 167/6 + 28{2,2,2} + 28{4,2} + 56C5
         + 8{5,2} + 28C6 + 8C7 + A8.

Using the audited {2,2,2}=32/105, retained {4,2}=-23/420,
C5=1/36 and C6=-1/126 gives

    m8 = 217/6 + 8{5,2} + 8C7 + A8.

The exact {5,2}=1/8 and C7=-17/360 certificates then give

    m7 = 3439/180,    m8 = 3311/90 + A8.

See [the pairing audit](../results/pairing_model_audit.md) for the
lower-order correction and [the seventh-order ledger](m7_ledger_derivation.md).

## Class evaluation status

| Class | Value | Status |
|---|---:|---|
| {2^4} | 1661/3780 | Exact continuum-model certificate |
| {4,2,2} | -127/840 | Numerical candidate |
| {6,2} | -1/20 | Numerical candidate |
| {4,4}, C8 | Unresolved | Evaluation needed |

The [pair-cycle proof](paired_cycle_flow_polytopes.md) certifies the
first row. If the two mixed-class candidates are also certified, their
combined contribution with {2^4} is 1801/7560, so

    A8 = 1801/7560 + T8,    T8 = {4,4} + C8.

## Target geometry

The audited degree-three certificate already gives the conditional
counting conversion 2025/2519, approximately 80.389%. The calculator
therefore imposes no extra eighth-moment cap for 80% or 80.2%.

At the default seventh-order inputs, the Hankel floor requires

    A8 >= 6775529/26142480.

For an 81% conditional target the permitted interval is

    6775529/26142480 <= A8 <= 8128408/16967475.

The corresponding 90% upper cap is 28456753/106766100.
These bounds constrain supplied model moments; they do not establish
the remaining class values or any unconditional zeta-zero result.

[The Schur-complement derivation](k8_target_geometry.md) gives the generic
curves. `scripts/m8_ledger.py` uses the exact seventh-order model inputs
by default and accepts explicit overrides for alternative scenarios.
