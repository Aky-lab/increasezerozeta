# Research overview

## Definitions and analytic scope

The continuum cumulants and outer prefix walks follow
[JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8),
revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`. A class signature
denotes the aggregate over all placements of those block sizes on a cycle.
Each block's frequencies sum to zero, with one frequency eliminated.

The [finite-model identities](notes/model_local_identities.md) prove
singleton deletion and three-block vanishing on outer-overlap support.
These identities assemble the finite model ledger. The class-d replacement,
arithmetic-to-continuum transport and spectral/counting interface remain
separate analytic obligations.

## Exact model evaluations

| Class | Exact value | Method |
|---|---:|---|
| {2} | 1/3 | Pair-cycle lift |
| {2,2} | 4/15 | Pair-cycle lift |
| {2,2,2} | 32/105 | Pair-cycle lift and simplex audit |
| {2^4} | 1661/3780 | Pair-cycle lift |
| C4 | -1/60 | Pure-cycle lift |
| C5 | 1/36 | Pure-cycle lift and spectator anchor |
| C6 | -1/126 | Pure-cycle lift |
| C7 | -17/360 | Pure-cycle lift |
| C8 | 157/4032 | Pure-cycle lift |
| {4,2} | -23/420 | Mixed-cycle lift |
| {5,2} | 1/8 | Weighted spectator lattice calculation |
| {6,2} | -563/11340 | Mixed-cycle lift |
| {4,2,2} | -127/840 | Mixed-cycle lift |
| {4,4} | 23/4536 | Mixed-cycle lift |

The [pair-cycle proof](notes/paired_cycle_flow_polytopes.md),
[pure-cycle proof](notes/pure_cycle_flow_polytopes.md),
[mixed-cycle proof](notes/mixed_cycle_flow_polytopes.md) and
[spectator proof](notes/spectator_52_exact_lattice.md) fix the volume
normalization and polynomial-degree bounds before interpolation.
Their records include source hashes and held-out integer counts.

The {5,2} distance components are U1=5/504, U2=1/360 and U3=13/2520,
with {5,2}=7(U1+U2+U3). The defining distance-three fixed prefix set is
{0,s2,s3,s4}. Adding s1 changes the integral; see
[the spectator audit](results/spectator_actual_walk_audit.md).

The [pairing audit](results/pairing_model_audit.md) separates adjacent
and nested noncrossing patterns, whose integrals are 3/70 and 17/420.
Exact integration on 24 simplex cells independently confirms the corrected
three-pair aggregate 32/105. The {6,2} certificate also rejects the earlier
-1/20 candidate, which differs from the exact value by 1/2835.

## Complete eighth-order ledger and consumption

The exact model sequence through eighth order is

    (m0,...,m8) = (1,1,4/3,2,13/4,101/18,640/63,3439/180,747361/20160).

The new eighth-order aggregate is

    A8 = {2^4}+{4,2,2}+{4,4}+{6,2}+C8 = 633/2240,
    m8 = 3311/90+A8.

[The full certificate](notes/eighth_order_certificate.md) includes the
pure C8 counts, an independent compiler and unsorted cube checks.
The [seventh-order](notes/m7_ledger_derivation.md) and
[eighth-order](notes/m8_ledger_and_80_target.md) notes give Bell-class
multiplicities and exact assembly.

For H_n=(m_(i+j)), the origin-mass bound is

    lambda_n(0) = 1/(e0^T H_n^-1 e0).

The [Christoffel derivation](notes/christoffel_tower.md) gives the
degree-three value 247/2519. The complete eighth-order model gives
lambda4(0)=12241115/162540559. Under the analytic counting interface,
1-2*lambda4(0)=138058329/162540559, approximately 84.938%.
This is a conditional conversion, not an established zeta-zero bound.
The [target geometry](notes/k8_target_geometry.md) separately describes
the positive moment cone and alternative supplied inputs.

## Open research

1. Independently review the finite certificates, normalization and grouping corrections.
2. Prove the arithmetic moment identities and their continuum transport through eighth order.
3. Establish the spectral/counting interface for these supplied moments.
4. Resolve the fixed-P tail obligation in the separate [fourth-moment route](notes/class_subtracted_universality.md) and [arithmetic analysis](notes/arithmetic_to_continuum.md).
5. Assess further model orders after the analytic interface is validated.

A certified finite integral, an assembled model moment and an established
arithmetic zeta moment are distinct stages. The remaining analytic work is
essential to any unconditional simple-zero bound.
