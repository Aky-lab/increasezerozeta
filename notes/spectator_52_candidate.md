# The {5,2} model integral

The corrected integral is exactly **{5,2}=1/8**, with U1=5/504, U2=1/360 and U3=13/2520. See [the lattice proof](spectator_52_exact_lattice.md) and [integer record](../results/spectator_52_exact_lattice_2026-10-05.json).

## Numerical discrepancy

An earlier spectator reduction used {0,s1,s2,s3,s4} as the distance-three fixed set. The defining walk uses {0,s2,s3,s4}. Including s1 produces a different integral and the numerical ladder near 7/72.

At (s1,s2,s3,s4)=(4/5,1/5,3/10,2/5), the actual spectator integral is 2/75; the extra-prefix formula gives 1/375. The [reduction proof](spectator_52_eliminate_v.md) identifies the error, and the [numerical audit](../results/spectator_actual_walk_audit.md) compares both definitions on the same grids.

The initial numerical convergence toward 1/8 was consistent with the corrected model. It was not a proof; exact lattice evaluation supplies the value.

## Seventh-order ledger

The conditional ledger is

    m7 = 1717/90 + {5,2} + C7 = 6913/360 + C7.

C7=-17/360 remains a numerical candidate. That assumption would give m7=862/45. Certification of C7 and arithmetic transport are separate requirements.
