# increasezerozeta

Exact computation and analytic research on lower bounds for the proportion of simple zeros of the Riemann zeta function on the critical line.

The project develops reproducible mathematical calculations toward stronger bounds. Its exact model results are available below; the remaining research is to evaluate the eighth-order classes and prove their connection to zeta moments.

## Results and scope

The corrected continuum model integral **{5,2}=1/8** has an exact integer-lattice evaluation, with components 5/504, 1/360 and 13/2520. The [proof](notes/spectator_52_exact_lattice.md) establishes the polynomial degree and period used in the calculation; the [integer record](results/spectator_52_exact_lattice_2026-10-05.json) includes source hashes and held-out checks.

The pure-cycle integral **C7=-17/360** is also exact. Its [flow-polytope proof](notes/pure_cycle_flow_polytopes.md) reduces the calculation to nine integer lattice counts and recovers the known lower-order constants.

The [pair-cycle calculation](notes/paired_cycle_flow_polytopes.md) gives **{2^4}=1661/3780**. Its lower-order checks exposed a grouping error: the direct three-pair aggregate is **32/105**, because adjacent and nested noncrossing patterns have different integrals. The [pairing audit](results/pairing_model_audit.md) includes an independent exact integration of both patterns.

Retaining the framework's other class inputs and local identities, the corrected ledger gives m6=640/63 and m7=3439/180. Exact Christoffel/Hankel algebra then gives lambda_3(0)=247/2519 and the conditional counting conversion 2025/2519, approximately 80.389%. These percentages depend on the unresolved analytic framework; they are not established zeta bounds.

These are model and algebraic results. The arithmetic transport to zeta moments and eighth-order class evaluations remain open. They do not yet establish a new unconditional zeta-zero bound.

See [RESEARCH.md](RESEARCH.md) for definitions, references and open problems.

## Reproduction

Run all ten project checks from the repository root with Python 3.12 or later. They use only the standard library:

```sh
python scripts/verify_project.py
```

This checks finite algebra, exact moment targets, actual-walk examples, certificate records and source hashes, independent pairing counts and cell integrals, Bell(7)/Bell(8) enumeration, and both moment ledgers. Use `--output verification.json` to save the full results.

The seventh- and eighth-order calculators are available directly:

```sh
python scripts/m7_ledger.py
python scripts/m8_ledger.py
```

The exact lattice evaluator requires NumPy. The recorded calculation used Python 3.12.12 and NumPy 2.3.5:

```sh
python -m pip install -r requirements-lattice.txt
python scripts/spectator_52_exact_lattice.py --out results/lattice_reproduction.json
```

To check the stored certificate, source hashes and rational calibration gates without repeating the large lattice sums:

```sh
python scripts/spectator_52_exact_lattice.py --out results/spectator_52_exact_lattice_2026-10-05.json --verify-only
```

The [spectator audit](results/spectator_actual_walk_audit.md) documents the distance-three prefix correction and the numerical comparison.

Recompute the pure-cycle counts and their independent term checks with NumPy:

```sh
python scripts/pure_cycle_lattice.py --out results/pure_cycle_reproduction.json
python scripts/verify_pure_cycle_lattice.py --out results/pure_cycle_checks.json
```

Recompute the pair-cycle counts through four pairs with NumPy:

```sh
python scripts/paired_cycle_lattice.py --out results/paired_cycle_reproduction.json
python scripts/verify_pairing_model.py
```

The standard-library certificate checks verify recorded integer differences, held-out counts and source hashes. The NumPy tools perform the full lattice enumeration. Shared conditional moment inputs are defined in `scripts/model_moments.py`.

## Project files

- [RESEARCH.md](RESEARCH.md): model conventions, references and open mathematical problems.
- [Seventh-order ledger](notes/m7_ledger_derivation.md) and [eighth-order targets](notes/m8_ledger_and_80_target.md): the current calculations and the unresolved classes.
- `scripts/`: executable checks, exact evaluators and exploratory numerical tools.
- `results/`: integer records and computational audits supporting the mathematical notes.

## Research support

The [research proposal](sponsorship/BRIEF.md) describes a six-week pilot for certificate review, higher-moment computation and analytic transport.
