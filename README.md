# increasezerozeta

Exact computation and analytic research on lower bounds for the proportion of simple zeros of the Riemann zeta function on the critical line.

The project develops reproducible mathematical calculations toward stronger bounds. All class evaluations needed for the finite continuum model through eighth order are certified. The model also has a Haar-unitary Gram realization, almost sure convergence to a unique limiting spectral distribution, and a joint Gaussian limit for polynomial spectral statistics. Arithmetic transport to zeta moments remains a central open problem.

## Results and scope

The corrected continuum model integral **{5,2}=1/8** has an exact integer-lattice evaluation, with components 5/504, 1/360 and 13/2520. The [proof](notes/spectator_52_exact_lattice.md) establishes the polynomial degree and period used in the calculation; the [integer record](results/spectator_52_exact_lattice_2026-10-05.json) includes source hashes and held-out checks.

The pure-cycle integral **C7=-17/360** is also exact. Its [flow-polytope proof](notes/pure_cycle_flow_polytopes.md) reduces the calculation to nine integer lattice counts and recovers the known lower-order constants.

The [pair-cycle calculation](notes/paired_cycle_flow_polytopes.md) gives **{2^4}=1661/3780**. Its lower-order checks exposed a grouping error: the direct three-pair aggregate is **32/105**, because adjacent and nested noncrossing patterns have different integrals. The [pairing audit](results/pairing_model_audit.md) includes an independent exact integration of both patterns.

The [complete eighth-order certificate](notes/eighth_order_certificate.md) gives {4,2,2}=-127/840, {4,4}=23/4536, {6,2}=-563/11340 and C8=157/4032. It rejects the earlier {6,2}=-1/20 guess. All mixed-class orbits and the pure cycle pass held-out integer-count predictions and independent small-grid checks.

With [finite-model local identities](notes/model_local_identities.md), the exact ledgers give m6=640/63, m7=3439/180 and m8=747361/20160. Rational Hankel inversion gives lambda_4(0)=12241115/162540559 and a **conditional counting conversion of approximately 84.938%**. This depends on the unresolved analytic interface and is not an established zeta bound.

These are finite model and algebraic results. Arithmetic transport to zeta moments remains open; the calculations do not yet establish a new unconditional zeta-zero bound.

See [RESEARCH.md](RESEARCH.md) for definitions, references and open problems.

## Structural research

[Centered Ehrhart reciprocity](notes/centered_reciprocity.md) reduces an eighth-order class certificate to four fit counts and one held-out count. The new method reproduces all fourteen class integrals, including a separate period-one certificate for {5,2}. Its polynomials match all 119 earlier aggregate counts and 641 orbit counts.

[The Haar-unitary Gram identity](notes/cue_gram_model.md) identifies the finite network model with exact expected moments of a positive random matrix. A separate integer Laurent-polynomial calculation checks every moment through eighth order at matrix sizes 1, 2, 3 and 4.

[An occupancy and moment-growth argument](notes/cue_limit_determinacy.md) proves that the expected spectral measures converge to a unique distribution with the continuum moments at every order. These derivations use established reciprocity, CUE and moment-problem results; priority for their specific application remains to be assessed.

[Connected multi-cycle counts](notes/cue_gram_fluctuations.md) strengthen this to almost sure convergence of the random empirical spectrum and joint Gaussian fluctuations of fixed polynomial statistics. They give exact finite-size joint cumulants. The order-two/order-three limiting covariance matrix is [[1/10,1/3],[1/3,79/70]], with positive determinant 11/6300. Independent Weyl integration checks the joint certificates at sizes 1 through 5.

[Logarithmic control](notes/cue_gram_logdet.md) proves that the limiting model measure has no atom at zero. The exact identity E[log det(G_n)]/n=H_n-1-log n gives a uniform bound on small-eigenvalue mass, which survives the weak limit. This all-order model conclusion is stronger than the finite-moment origin-mass certificate; arithmetic counting transport remains open.

[The bandwidth extension](notes/cue_gram_bandwidth.md) proves strong spectral convergence, polynomial fluctuations and logarithmic control for rectangular CUE Gram matrices with m/n tending to any fixed 0<lambda<=1. The row law has no zero atom; the column law has zero mass exactly 1-lambda. Its second-moment fluctuation variance is 2*lambda^3/15-(2*lambda-1)_+^5/(30*lambda^2), with an overlap threshold at lambda=1/2. Independent exact Weyl checks cover every 1<=m<=n<=5.

[The global sine-process comparison](notes/sine_gram_identification.md) identifies the growing-window sine Gram law with the CUE column law, with almost sure convergence of all moments. At unit bandwidth it has the certified continuum moments; at smaller bandwidth its limiting zero atom is exactly 1-lambda, even though every finite sine Gram matrix is invertible. The proof controls long interactions and specifies the point/frequency normalization.

[The literature comparison](notes/literature_map.md) credits the earlier sine-process Gram motivation, independent-phase Vandermonde work, CUE pair statistics and stationary determinantal-process background. Arithmetic transport and independent review of the model proofs remain open.

[The quantitative arithmetic-core proof](notes/arithmetic_core_limit.md) establishes the explicitly defined class-subtracted prime-power model as -1/48+O(1/log T). It proves the universal coefficient asymptotic, controls exceptional moduli uniformly, and derives the overlap integral from cube and simplex volumes.

[A cutoff obstruction](notes/fixed_cutoff_obstruction.md) proves the uniform bound C_ell,P=O((1+log P)/log T). Every subpower cutoff, including any fixed power of log T, leaves the entire leading model contribution in its tail. An absolute envelope for that same tail has lower limit at least 1/48, so the older 0.0111 finite-height charge cannot be promoted to such an asymptotic bound.

[Prime-pair progression transport](notes/prime_pair_progression_transport.md) uses the published Matomaki–Radziwill–Tao theorem to control arbitrary divisor subfamilies of shifts h=qk. It supplies a weighted error budget and resolves the full-dyadic, one-chain aggregation step. Coupled prime chains and the actual zeta-moment reduction remain open.

[Short-window Fourier transport](notes/short_window_fourier_transport.md) applies published higher-uniformity results to replace single-prime sums by a presieved model. It derives an integer-position mean-square bound and controls coupled products through weighted start histograms. The averaged range begins at X^(1/3+epsilon); a longer pointwise range begins at j^(5/8+epsilon).

[The finite presieved comb](notes/presieved_pair_transport.md) gives explicit Fourier coefficients for a power-sized divisor approximation and identifies its pair main term uniformly. Parseval then yields an actual prime-pair estimate averaged over short-window positions and shifts, with arbitrary logarithmic savings and a divisor-family consumption budget. Prescribed individual windows, higher prime correlations and the zeta normalization still require separate arguments. This is a deduction from published uniformity input; priority for the specific formulation has not been established.

[The resolved lock frame](notes/resolved_lock_frame.md) derives the correct two-dimensional Fourier transform for four-point locks. Exact nonnegative examples show that separate autocorrelations can agree while resolved squared counts differ. A full-lock stability theorem gives an averaged prime-model replacement, but its restriction to rectangle locks loses a window-length factor. The note identifies the cubic uniformity norm required to repair that step and keeps the published qualitative and quantitative inputs distinct.

[Qualitative rectangle transport](notes/averaged_rectangle_transport.md) uses published cubic uniformity to replace actual four-prime rectangle counts by a finite presieved model, averaged over integer starts and both shifts. It proves a vanishing squared error at the natural Y^4 scale and gives an explicit four-point local product. An exact admissible-residue budget avoids a growing normalization loss. The saving is qualitative; weighted, dilated and prescribed-start zeta transport remains open.

[The averaged Euler-tail estimate](notes/rectangle_singular_series_tail.md) replaces that finite main term by the full four-point singular series on nondegenerate shifts. For every fixed integer r>=1, the r-th absolute tail moment over the rectangle box is O_r(Y^2/w^r), uniformly in height and cutoff. Its mean-square case gives the full Hardy--Littlewood main term with the same qualitative prime saving. The proof uses positive divisor moments and repeated-prime partitions; a uniform pointwise tail bound would fail.

## Reproduction

Run all twenty-three project checks from the repository root with Python 3.12 or later. They use only the standard library:

```sh
python scripts/verify_project.py
```

This checks finite algebra, Euler convolution and prime-power coefficients, overlap geometry and Ramanujan cutoff decompositions, divisor-family error bookkeeping, exact moment targets, actual-walk examples, certificate records and source hashes, independent pairing counts and cell integrals, scalar mixed-class counts, the reduced reciprocity certificates, independent Weyl moments and joint cumulants, connected multi-cycle identities, Bell(7)/Bell(8) enumeration, and both moment ledgers. Use `--output verification.json` to save the full results.

The faster route to recompute all fourteen network certificates requires NumPy:

```sh
python -m pip install -r requirements-lattice.txt
python scripts/reciprocity_certificates.py --out results/reciprocity_reproduction.json
python scripts/verify_reciprocity.py --record results/reciprocity_reproduction.json
python scripts/verify_cue_model.py --record results/reciprocity_reproduction.json
```

The generation step took about nine seconds in the recorded environment. These commands check the regenerated certificate; without `--record` the verifiers check the committed one. The independent Gram verifier computes its Weyl integrals afresh. Older evaluators below retain the original, larger sample ranges for comparison.

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

Recompute the mixed classes and the pure eighth-order cycle with NumPy:

```sh
python scripts/mixed_cycle_lattice.py --out results/mixed_cycle_reproduction.json
python scripts/pure_cycle_eight.py --out results/cycle_eight_reproduction.json
python scripts/verify_cycle_eight_lattice.py --out results/cycle_eight_checks.json
```

The standard-library certificate checks verify recorded integer differences, held-out counts and source hashes. The NumPy tools perform the full lattice enumeration. Shared finite model moments are defined in `scripts/model_moments.py`. The [verification record](results/project_verification_2026-10-05.json) identifies all twenty-three executed checks.

## Project files

- [RESEARCH.md](RESEARCH.md): model conventions, references and open mathematical problems.
- [Seventh-order ledger](notes/m7_ledger_derivation.md) and [eighth-order targets](notes/m8_ledger_and_80_target.md): the exact class assembly and conditional counting targets.
- `scripts/`: executable checks, exact evaluators and exploratory numerical tools.
- `results/`: integer records and computational audits supporting the mathematical notes.

## Research support

The [research proposal](sponsorship/BRIEF.md) describes a six-week pilot for certificate review, higher-moment computation and analytic transport.
