# Research overview

## Model conventions

The higher-moment conventions follow [JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8), revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`. The class-d replacement, cumulant normalization and frozen-singleton identities are part of the analytic framework requiring review.

The notation {5,2} denotes the aggregate over all 21 placements of a five-block and a pair in a seven-cycle. It equals 7(U1+U2+U3), where the three U terms correspond to cyclic distances one, two and three.

## Exact model calculation

The [spectator reduction](notes/spectator_52_eliminate_v.md) integrates the pair frequency exactly. In distance three the fixed prefix set is {0,s2,s3,s4}; adding s1 changes the integral.

[Weighted integer lattice evaluation](notes/spectator_52_exact_lattice.md) then gives:

| Quantity | Exact value |
|---|---:|
| U1 | 5/504 |
| U2 | 1/360 |
| U3 | 13/2520 |
| {5,2} | 1/8 |
| Pure C5 anchor | 1/36 |

The chamber normals have binary entries up to sign. Their determinant bound gives period dividing six, and the weighted sum has degree at most eight on each residue class. Nine samples determine the leading coefficient; a tenth sample checks the prediction. The [raw record](results/spectator_52_exact_lattice_2026-10-05.json) contains exact sums and source hashes.

The [pure-cycle flow-polytope method](notes/pure_cycle_flow_polytopes.md) gives C7=-17/360. Each lifted term is an integral bounded circulation polytope, so its lattice count is polynomial with period one. The calculation recovers C4=-1/60, C5=1/36 and C6=-1/126; independent checks reproduce 150 C5 term certificates and twelve C6 term certificates. See [the cycle counts](results/pure_cycle_flow_2026-10-05.json) and [independent checks](results/pure_cycle_independent_checks_2026-10-05.json).

## Moment consumption

For a moment matrix H_n, the mass-at-origin bound is

    lambda_n(0) = 1 / (e0^T H_n^-1 e0).

The supplied six-moment sequence gives lambda_3(0)=1415/13891. Under the spectral/counting interface, this yields the conditional simple-zero proportion 11061/13891. See [the Christoffel derivation](notes/christoffel_tower.md) and [eighth-order target geometry](notes/k8_target_geometry.md).

The [fourth-moment reduction](notes/class_subtracted_universality.md) and [arithmetic-to-continuum analysis](notes/arithmetic_to_continuum.md) form a separate analytic route, with a fixed-P tail obligation still open.

## Higher-moment frontier

The [seventh-order ledger](notes/m7_ledger_derivation.md), [eighth-order ledger](notes/m8_ledger_and_80_target.md), Bell(7)/Bell(8) enumerators and candidate computations are included in this repository.

The conditional seventh-order ledger is

    m7 = 1717/90 + {5,2} + C7 = 6913/360 + C7.

Both new model inputs are exact: {5,2}=1/8 and C7=-17/360. The ledger therefore gives m7=862/45. Its interpretation as a zeta moment still requires the frozen-singleton, vanishing and arithmetic transport statements.

The new eighth-order aggregate is

    A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8.

The model candidates {2^4}=1661/3780, {4,2,2}=-127/840 and {6,2}=-1/20 still require exact evaluation or rigorous enclosure. Bounds on {4,4}+C8 must use the current seventh-order inputs.

## Open problems

1. Independently check the {5,2} certificate and reference normalization.
2. Independently review the C7 flow-polytope certificate and volume normalization.
3. Prove the arithmetic transport for the seventh-order classes.
4. Certify the eighth-order classes and bound their remaining aggregate.
5. Prove eighth-order transport and apply the exact moment-consumption engine.

A model integral, numerical candidate and analytically established zeta moment have different roles. The outstanding transports are essential to any unconditional simple-zero bound.
