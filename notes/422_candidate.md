# Exact {4,2,2} class

The class consists of one connected four-block and two pairs on an
eight-cycle. Its 210 placements split into 22 dihedral orbits.
The [mixed-cycle certificate](mixed_cycle_flow_polytopes.md) gives

    {4,2,2} = -127/840.

Each lifted signed term is an integral network polytope of dimension
nine. Ten integer counts determine each orbit's degree-nine polynomial;
the eleventh checks every orbit and the aggregate. The
[raw record](../results/mixed_cycle_flow_2026-10-05.json) includes all
counts, multiplicities, source hashes and exact orbit values.

The independent standard-library verifier generates all placements,
uses scalar closed frequencies and reconstructs the outer walk. Its
inner cumulant compiler uses permutations and cuts, independent of the
primary coordinate elimination and sorted-frequency cache.

This is an exact finite continuum-model evaluation. Its arithmetic
interpretation is addressed by the remaining transport obligations.
