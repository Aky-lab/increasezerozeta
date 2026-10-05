# Verification records

- `baseline-2026-10-05.json` is a historical record of three algebra checks and an earlier spectator checker. Its lower-moment inputs predate the pairing audit, and its spectator checker used an extra distance-three prefix. It does not describe the current model.
- `actual-walk-2026-10-05.json` records the independent increment-based checker, exact counterexample and source hashes.
- [The lattice record](../../results/spectator_52_exact_lattice_2026-10-05.json) contains the complete corrected model evaluation {5,2}=1/8 and the C5=1/36 anchor.

- [The current verification record](../../results/project_verification_2026-10-05.json) contains twenty-seven local checks, including exact model certificates, independent Weyl integration, prime-transport normalization and rational overlap integration. Each executed script and all script sources are hashed; the record states the limits of finite verification.

The [spectator audit](../../results/spectator_actual_walk_audit.md) explains the prefix correction, and the [pairing audit](../../results/pairing_model_audit.md) explains the nested-pair correction. Analytic transport to zeta moments remains open.
