# The {5,2} model integral

The corrected integral is exactly **{5,2}=1/8**, with U1=5/504, U2=1/360 and U3=13/2520. See [the lattice proof](spectator_52_exact_lattice.md) and [integer record](../results/spectator_52_exact_lattice_2026-10-05.json).

## Numerical discrepancy

An earlier spectator reduction used {0,s1,s2,s3,s4} as the distance-three fixed set. The defining walk uses {0,s2,s3,s4}. Including s1 produces a different integral and the numerical ladder near 7/72.

At (s1,s2,s3,s4)=(4/5,1/5,3/10,2/5), the actual spectator integral is 2/75; the extra-prefix formula gives 1/375. The [reduction proof](spectator_52_eliminate_v.md) identifies the error, and the [numerical audit](../results/spectator_actual_walk_audit.md) compares both definitions on the same grids.

The initial numerical convergence toward 1/8 was consistent with the corrected model. It was not a proof; exact lattice evaluation supplies the value.

## Seventh-order ledger

The conditional ledger is

    m7 = 685/36 + {5,2} + C7 = 3439/180.

The separate [flow-polytope certificate](pure_cycle_flow_polytopes.md) gives C7=-17/360 exactly. The resulting model ledger is m7=3439/180; arithmetic transport remains open.

The lower-order base uses the audited three-pair value 32/105; see
[the pairing audit](../results/pairing_model_audit.md) and the
[finite-model local identities](model_local_identities.md) for its assembly.
