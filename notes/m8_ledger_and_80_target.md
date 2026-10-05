# Bell(8) ledger and conditional counting targets

The [local-identity proof](model_local_identities.md) establishes singleton
deletion and three-block vanishing for the finite continuum model.
Arithmetic transport and the spectral/counting interface remain open.

## Class assembly

Define

    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

Bell enumeration and singleton deletion give

    m8 = 167/6 + 28{2,2,2} + 28{4,2} + 56C5
         + 8{5,2} + 28C6 + 8C7 + A8.

The initial 167/6 includes all-singleton, one-pair, two-pair and pure
four-block contributions. Three-block classes vanish pointwise.
Every other inherited aggregate is multiplied by the number of choices
of its non-singleton endpoint set.

The certified {2,2,2}=32/105, {4,2}=-23/420, C5=1/36 and C6=-1/126
reduce the inherited ledger to 217/6+8{5,2}+8C7. Using {5,2}=1/8
and C7=-17/360 gives

    m7 = 3439/180,    m8 = 3311/90 + A8.

## Exact new classes

| Class | Exact continuum-model value |
|---|---:|
| {2^4} | 1661/3780 |
| {4,2,2} | -127/840 |
| {4,4} | 23/4536 |
| {6,2} | -563/11340 |
| C8 | 157/4032 |

The [pair-cycle](paired_cycle_flow_polytopes.md),
[mixed-cycle](mixed_cycle_flow_polytopes.md) and
[pure eighth-order](eighth_order_certificate.md) certificates give

    A8 = 633/2240,
    m8 = 747361/20160.

The {6,2}=-1/20 numerical candidate is rejected, rather than included
in this assembly. The [pairing audit](../results/pairing_model_audit.md)
documents the earlier lower-order grouping correction.

## Consumption and target geometry

The [exact Schur-complement derivation](k8_target_geometry.md) gives the
Hankel floor and generic target curves. At m7=3439/180,

    A8 >= 6775529/26142480

is required by the moment cone. The exact A8=633/2240 lies above this
floor. Direct fourth-degree rational Hankel inversion gives

    lambda4(0) = 12241115/162540559,
    1-2*lambda4(0) = 138058329/162540559 = 0.8493777174... .

The approximately 84.938% value is conditional on the analytic counting
interface. It does not establish a zeta-zero bound. The polynomial's
alternating coefficients verify its nonpositive-half-line certificate.

For comparison, the 81% conditional target has A8 cap
8128408/16967475 and the 90% cap is 28456753/106766100. The exact
aggregate meets the former and exceeds the latter. The degree-three
model certificate already meets 80% and 80.2% without an extra cap.

`scripts/m8_ledger.py` assembles the exact finite model and checks the
direct Hankel inverse, positivity and coefficient signs. Explicit
seventh-order overrides are available for alternative supplied scenarios.
