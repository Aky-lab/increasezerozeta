# Research guide

## Question and conventions

The objective is to improve lower bounds for the proportion of simple zeros of the Riemann zeta function on the critical line. The principal open step is an arithmetic trace estimate for the actual compressed Weil matrix. Exact continuum-model moments provide candidate certificates for that estimate.

Cumulants and outer prefix walks follow [JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8), revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`. A class signature denotes the sum over placements of its block sizes on a cycle. Each block's frequencies sum to zero. The [local identities](notes/model_local_identities.md) derive singleton deletion and three-block vanishing on outer-overlap support.

## Exact finite model

| Class | Value | Derivation |
|---|---:|---|
| {2} | 1/3 | [Pair-cycle lift](notes/paired_cycle_flow_polytopes.md) |
| {2,2} | 4/15 | Pair-cycle lift |
| {2,2,2} | 32/105 | [Pair-cycle integration](results/pairing_model_audit.md) |
| {2^4} | 1661/3780 | Pair-cycle lift |
| C4 | -1/60 | [Pure-cycle lift](notes/pure_cycle_flow_polytopes.md) |
| C5 | 1/36 | Pure-cycle lift |
| C6 | -1/126 | Pure-cycle lift |
| C7 | -17/360 | Pure-cycle lift |
| C8 | 157/4032 | [Eighth-order certificate](notes/eighth_order_certificate.md) |
| {4,2} | -23/420 | [Mixed-cycle lift](notes/mixed_cycle_flow_polytopes.md) |
| {5,2} | 1/8 | [Weighted lattice evaluation](notes/spectator_52_exact_lattice.md) |
| {4,2,2} | -127/840 | Mixed-cycle lift |
| {4,4} | 23/4536 | Mixed-cycle lift |
| {6,2} | -563/11340 | Mixed-cycle lift |

The [seventh-order ledger](notes/m7_ledger_derivation.md) and [eighth-order ledger](notes/m8_ledger_and_80_target.md) assemble

```text
(m0,...,m8) = (1, 1, 4/3, 2, 13/4, 101/18, 640/63, 3439/180, 747361/20160).
A8 = {2^4}+{4,2,2}+{4,4}+{6,2}+C8 = 633/2240.
```

For the Hankel matrix `H_n=(m_(i+j))`, the origin-mass certificate is `lambda_n(0)=1/(e0^T H_n^-1 e0)`. The [Christoffel calculation](notes/christoffel_tower.md) gives `lambda_3(0)=247/2519`; degree four gives `lambda_4(0)=12241115/162540559`. The associated counting conversions are approximately 80.389% and 84.938%. Applying them to zeta zeros requires the unresolved arithmetic hypotheses.

The certificates use polynomial lattice counts with justified degree and period bounds. Held-out counts check predictions; separate scalar enumeration and exact integration check the construction. [Centered Ehrhart reciprocity](notes/centered_reciprocity.md) reduces the number of fit counts. The [verification record](results/project_verification_2026-10-05.json) identifies executed checks and source hashes.

## Counting argument and arithmetic target

| Note | Content and status |
|---|---|
| [Spectral counting bridge](notes/spectral_counting_bridge.md) | Finite inertia/collar argument and zero-tail reduction. A trace cap of 1/10 is sufficient for the 80% target; the actual cap is open. |
| [Lean formalization](formal/README.md) | 31 kernel-checked theorems for scalar majorization, boundedness, resolvent minorization, the rational-pole sum of squares, exact constants and conditional finite counting. Matrix identification, inertia and arithmetic estimates are supplied hypotheses. |
| [Actual prime trace](notes/actual_prime_trace.md) | Proof draft reducing the Weil matrix to a prime matrix and evaluating its multiplicatively balanced terms. The required signed off-balance estimate is open. |
| [Second prime moment](notes/prime_second_moment.md) | Deduction from sampled Montgomery–Vaughan input, with finite-frame errors. Higher-word boundary costs remain. |
| [Bounded resolvent bridge](notes/bounded_resolvent_bridge.md) | Bounded rational certificate, partial fractions and compression estimate. The signed arithmetic resolvent combination remains open. |
| [Prime-removal expansion](notes/prime_removal_resolvent.md) | Deterministic remainder and covariance-map bounds. Conditional-phase cancellation and the nonlinear covariance trace remain open. |
| [Reflected-square certificate](notes/mirror_resolvent_certificate.md) | A certificate in [0,2] and an exact rational majorant. A lower bound of 0.91 at 11i/40 would imply approximately 80.13%; the model supports the threshold, while the actual arithmetic estimate remains open. |
| [Fixed-taper projection](notes/fixed_taper_projection.md) | A dimension-independent leakage bound, exact Fourier kernel, high-frequency gauge and a whole-space norm obstruction. The projected arithmetic statistic remains unestimated. |

At unit bandwidth, the balanced second, fourth and sixth prime moments are `1/3`, `4/15` and `32/105`. The cubic polynomial certificate's balanced trace is approximately 0.160833; meeting its 0.1 cap requires a signed off-balance contribution at most -0.060833. Small prime increments and operator-norm bounds alone do not establish this cancellation.

The bounded alternative uses `f(x)=q3(x)^2/(1+(1932/2519)^(2/3)*x^2)^3`. It majorizes the nonpositive half-line and reduces the target to three resolvent powers at one nonreal point. This route avoids assuming a sixth-moment limit, while retaining an explicit arithmetic estimate to prove.

The [rational one-point criterion](notes/mirror_resolvent_certificate.md#5-a-rational-one-point-criterion) applies a globally verified majorant directly to the actual symmetric matrix. Its first two moments and resolvent transfer show that `liminf Re m_T(11i/40) >= 113102/124503` suffices for 80%. The slightly stronger bound at 0.91 gives approximately 80.13%. The scalar proof uses an exact degree-ten Sturm certificate; it assumes no higher arithmetic moment limits. Independently, a Bessel disk from the eight certified model moments gives a value above 0.9112, confirming model compatibility. The [moment-disk analysis](notes/mirror_resolvent_certificate.md#6-model-moment-compatibility) also quantifies a positive-model obstruction to requiring a value of one at i/4.

## Random-matrix model

These proof notes concern the model ensemble. Their identification with arithmetic zeta moments is a separate problem; the analytic notes have not undergone independent peer review.

| Note | Result developed |
|---|---|
| [Haar-unitary Gram model](notes/cue_gram_model.md) | Finite expected-moment identity with the network model |
| [Moment determinacy](notes/cue_limit_determinacy.md) | Occupancy bounds and a unique limiting spectral measure |
| [Polynomial fluctuations](notes/cue_gram_fluctuations.md) | Connected multi-cycle counts, almost sure convergence and joint Gaussian limits |
| [Log-determinant control](notes/cue_gram_logdet.md) | Small-eigenvalue control and absence of a zero atom |
| [Variable bandwidth](notes/cue_gram_bandwidth.md) | Rectangular Gram laws and explicit fluctuation variance |
| [Sine-process comparison](notes/sine_gram_identification.md) | Growing-window sine Gram law and its CUE column-law identification |
| [Literature map](notes/literature_map.md) | Primary sources, related results and attribution |

## A separate fourth-moment transport route

The [class-subtracted framework](notes/class_subtracted_universality.md) and [arithmetic-to-continuum specification](notes/arithmetic_to_continuum.md) study prime-power lock families. Its explicitly defined arithmetic model has a [quantitative core limit](notes/arithmetic_core_limit.md) of `-1/48+O(1/log T)`. The [cutoff obstruction](notes/fixed_cutoff_obstruction.md) shows that a subpower cutoff retains none of the leading contribution.

| Stage | Notes |
|---|---|
| Prime pairs and weighted windows | [Progression transport](notes/prime_pair_progression_transport.md), [short-window Fourier transport](notes/short_window_fourier_transport.md), [presieved pairs](notes/presieved_pair_transport.md) |
| Four-point locks | [Resolved Fourier frame](notes/resolved_lock_frame.md), [averaged rectangle transport](notes/averaged_rectangle_transport.md), [singular-series tails](notes/rectangle_singular_series_tail.md) |
| Dilation and prescribed starts | [Dilated geometry](notes/dilated_rectangle_geometry.md), [weighted dilation transport](notes/weighted_dilated_prime_transport.md), [full-window transport](notes/full_window_prime_transport.md) |
| Remaining arithmetic range | [Dilation distribution](notes/dilation_core_coverage.md), [physical coefficient profiles](notes/progression_range_profile.md), [manuscript source](papers/dilation_core_distribution.tex) |

The published-input deductions cover specified averaged windows or limited dilation ranges. Positive-power weighted lock transport and the final zeta-moment reduction remain open. The exact square-root coefficient profile covers `107/600` of the leading model term; this quantifies the range gap rather than closing it.

## Research priorities

1. Review the actual-matrix counting bridge and its normalization, including the correspondence with the formal finite theorem.
2. Estimate the resolvent-inserted two-prime statistic at 11i/40, using the prime-removal phase and covariance identities or another arithmetic approach.
3. Obtain independent mathematical review of the CUE/sine comparison and continuum certificate definitions.
4. Resolve the weighted positive-power transport gap in the fourth-moment route.

The [support proposal](sponsorship/BRIEF.md) scopes an initial six-week pilot around these tasks.
