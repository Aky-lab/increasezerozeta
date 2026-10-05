# increasezerozeta

Exact computations, proof notes and formal verification for research on simple zeros of the Riemann zeta function.

The project studies a continuum moment model and the arithmetic estimates needed to connect it to zeros on the critical line. It contains exact class evaluations through eighth order, random-matrix interpretations, and a conditional spectral counting argument. **The arithmetic estimate needed for a new zeta-zero bound remains open.**

## Start here

| Material | Contents |
|---|---|
| [Research guide](RESEARCH.md) | Definitions, exact values, proof status and a map of the notes |
| [Counting bridge](notes/spectral_counting_bridge.md) | The actual-matrix hypothesis sufficient for an 80% simple-zero bound |
| [Bounded resolvent criterion](notes/bounded_resolvent_bridge.md) | An alternative arithmetic target using three resolvent powers |
| [Lean proofs](formal/README.md) | 18 kernel-checked scalar and finite counting theorems, with their assumptions |
| [Computational records](results/project_verification_2026-10-05.json) | Local verification results, source hashes and reproduction details |

The finite model has moments

```text
(m0,...,m8) = (1, 1, 4/3, 2, 13/4, 101/18, 640/63, 3439/180, 747361/20160).
```

They yield an exact degree-four origin-mass certificate of `12241115/162540559`. Its counting conversion is approximately 84.938%, conditional on arithmetic transport. The [eighth-order certificate](notes/eighth_order_certificate.md) supplies the class values and their derivation.

## Reproduce the checks

From the repository root, using Python 3.12 or later:

```sh
python scripts/verify_project.py
```

This runs 32 checks using the Python standard library. With Lean 4.33.0 installed, include the formal proof audit:

```sh
python scripts/verify_project.py --with-lean --output verification.json
```

The Lean file uses `Std`; no Mathlib download is required. The audit checks every theorem's axiom dependencies and tests rejection of invalid proofs. Finite computations and formal proofs verify the stated identities and implications; analytic assumptions are documented in the corresponding notes.

To regenerate all fourteen network class certificates through eighth order, install NumPy and run:

```sh
python -m pip install -r requirements-lattice.txt
python scripts/reciprocity_certificates.py --out results/reciprocity_reproduction.json
python scripts/verify_reciprocity.py --record results/reciprocity_reproduction.json
python scripts/verify_cue_model.py --record results/reciprocity_reproduction.json
```

Individual evaluators and their independent checks are linked from the [research guide](RESEARCH.md).

## Review and support

Mathematical review is especially useful for the model normalization, the CUE/sine comparison and the connection between the spectral counting argument and the actual prime operator. Review comments should identify the statement, its assumptions and the step requiring justification.

The [support proposal](sponsorship/BRIEF.md) sets out a six-month research plan beginning with a six-week pilot. AI assists with derivation and implementation; accepted calculations require reproducible verification, and analytic proof notes remain subject to mathematical review.

The [literature map](notes/literature_map.md) records source attribution and related work. Reference conventions are pinned to [JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8).
