# Exact four-pair class {2^4}

The class sums the continuum-model integral over all 105 perfect
matchings of eight cyclic positions. Each matching has four frequencies,
outer prefix-walk overlap and weight |v1*v2*v3*v4|.

The exact value is **1661/3780**. The
[integral-flow proof](paired_cycle_flow_polytopes.md) establishes that
the integer count is a period-one polynomial of degree at most nine.
Ten samples determine each of the 17 dihedral orbit integrals; an
eleventh sample checks each prediction and their aggregate.

The [record](../results/paired_cycle_flow_2026-10-05.json) stores all
integer counts, multiplicities, orbit values and source hashes.
`scripts/verify_pairing_model.py` independently enumerates every matching
on small scalar grids without using dihedral reduction.

Earlier midpoint calculations suggested this fraction. Their claimed
three-pair calibration used an incorrect noncrossing grouping. The
[pairing audit](../results/pairing_model_audit.md) explains why the direct
three-pair value is 32/105 and independently integrates the adjacent
and nested patterns. The exact four-pair value follows from its own
certificate rather than that calibration.

This certifies a finite continuum-model class. Its arithmetic
interpretation and the remaining eighth-order classes are still open.
