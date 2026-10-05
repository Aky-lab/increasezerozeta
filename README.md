# increasezerozeta

Exact computation and analytic research on lower bounds for the proportion of simple zeros of the Riemann zeta function on the critical line.

## Results and scope

The corrected continuum model integral **{5,2}=1/8** has an exact integer-lattice evaluation, with components 5/504, 1/360 and 13/2520. The [proof](notes/spectator_52_exact_lattice.md) establishes the polynomial degree and period used in the calculation; the [integer record](results/spectator_52_exact_lattice_2026-10-05.json) includes source hashes and held-out checks.

The repository also contains exact Christoffel/Hankel moment consumption and seventh/eighth-moment target geometry. The six-moment calculation gives lambda_3(0)=1415/13891 and the conditional simple-zero bound 11061/13891, approximately 79.627%.

These are model and algebraic results. Certification of C7 and the arithmetic transport to zeta moments remain open. They do not yet establish a new unconditional zeta-zero bound.

See [RESEARCH.md](RESEARCH.md) for definitions, references and open problems.

## Reproduction

The baseline checks use the Python standard library:

```sh
python scripts/verify_lemmas.py
python scripts/christoffel_exact.py
python scripts/k8_target.py
python scripts/verify_spectator_reduction.py
```

The exact lattice evaluator requires NumPy. The recorded calculation used Python 3.12.12 and NumPy 2.3.5:

```sh
python scripts/spectator_52_exact_lattice.py --out results/lattice_reproduction.json
```

To check the stored certificate, source hashes and rational calibration gates without repeating the large lattice sums:

```sh
python scripts/spectator_52_exact_lattice.py --out results/spectator_52_exact_lattice_2026-10-05.json --verify-only
```

The [spectator audit](results/spectator_actual_walk_audit.md) documents the distance-three prefix correction and the numerical comparison.

## Research support

The [research proposal](sponsorship/BRIEF.md) describes a six-week pilot for certificate review, higher-moment computation and analytic transport.
