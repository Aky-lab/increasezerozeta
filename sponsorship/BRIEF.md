# Research support proposal

## Project

Develop auditable exact computations for higher moments in research on simple zeros of the Riemann zeta function. A six-week pilot seeks **$1,000 in AI API credits** for independent derivations, computation and analytic proof development.

The immediate objective is independent review of the exact model and its structural proofs, further mathematical research, and progress on the lemmas connecting it to arithmetic zeta moments.

## Existing results

The corrected four-dimensional model integral {5,2} is exactly 1/8. Its components are 5/504, 1/360 and 13/2520. The [proof](../notes/spectator_52_exact_lattice.md) bounds the degree and period of the weighted lattice sums; the [integer record](../results/spectator_52_exact_lattice_2026-10-05.json) includes held-out checks and source hashes. A separate C5 anchor evaluates to 1/36.

The [pure-cycle calculation](../notes/pure_cycle_flow_polytopes.md) gives C7=-17/360 exactly using integral network-flow polytopes. Independent checks recover all 150 C5 term values from the reference certificates and twelve C6 term values.

The [pair-cycle certificate](../notes/paired_cycle_flow_polytopes.md) gives {2^4}=1661/3780. It also exposes a lower-order grouping error: adjacent and nested noncrossing matchings have integrals 3/70 and 17/420 respectively. Exact integration on 24 simplex cells independently confirms the corrected three-pair aggregate 32/105; see [the audit](../results/pairing_model_audit.md).

The repository provides a standard-library engine for exact Christoffel/Hankel moment consumption and eighth-order target geometry. An independent spectator checker passes 7,203 rational cases and detects 358 failures of an earlier distance-three formula. The [audit](../results/spectator_actual_walk_audit.md) identifies the missing-prefix error.

The [eighth-order certificate](../notes/eighth_order_certificate.md) completes the remaining finite classes: {4,2,2}=-127/840, {4,4}=23/4536, {6,2}=-563/11340 and C8=157/4032. It rejects the earlier {6,2}=-1/20 guess. Independent scalar frequency sums check the mixed classes; a separate compiler and unsorted cube enumeration check the pure eighth-order cycle.

The assembled model moments are m6=640/63, m7=3439/180 and m8=747361/20160. Exact degree-four Hankel algebra gives a conditional counting conversion of approximately 84.938%. The finite-model local identities have written proofs, and the [direct counting bridge](../notes/spectral_counting_bridge.md) supplies an explicit sufficient polynomial-trace condition on the actual Weil matrix. The required arithmetic trace bound remains unproved. These calculations do not establish a new unconditional zeta-zero bound.

See [the research overview](../RESEARCH.md) for conventions and proof obligations.

The [actual prime-matrix proof draft](../notes/actual_prime_trace.md)
now connects the pair-cycle certificates to the multiplicatively balanced
terms of the actual arithmetic matrix. It removes archimedean and pole
terms in a normalized Schatten norm and evaluates the balanced moments
through degree six as 1/3, 4/15 and 32/105 at unit bandwidth. The cubic
80% target becomes one precisely specified signed off-balance estimate.
Its balanced contribution is approximately 0.160833, so at least
0.060833 of negative off-balance contribution is needed. Developing and
independently reviewing that cancellation is a concrete support milestone;
the estimate itself remains open.

The [finite-frame second-moment proof](../notes/prime_second_moment.md)
now recovers the known arithmetic second moment in these conventions,
with explicit projection and sampled mean-value errors. Its second
off-balance contribution vanishes, leaving orders three through six
in the cubic target. The proof tracks the boundary cost that higher
moments still need to overcome and attributes its classical input.

## Structural research and computational value

The [bounded resolvent reduction](../notes/bounded_resolvent_bridge.md)
offers a second analytic milestone: prove one signed combination of
three actual-prime resolvent powers at a single nonreal point. Its bounded
certificate avoids high-moment outlier control, and the draft proves that
finite-mode compression has vanishing cost for this statistic. Exact
algebraic-field and block-matrix checks verify the reduction's identities;
the arithmetic estimate itself remains open.

The [prime-removal proof draft](../notes/prime_removal_resolvent.md)
now controls the local cubic remainder and quadratic higher-power cost,
with a geometric covariance-map bound. This makes the outstanding
conditional-phase and nonlinear resolvent statistics explicit. Exact
matrix checks audit the identities, and a counterexample rejects the
shortcut from small perturbations to phase cancellation. Evaluating the
remaining signed statistic is still an analytic support milestone.

[Centered reciprocity](../notes/centered_reciprocity.md) reduces the number of
fit counts for an order-B class to floor(B/2). At eighth order the largest
outer cube shrinks from 214358881 to 1679616 points. All fourteen class
certificates regenerate in about nine seconds in the recorded environment;
every earlier network count agrees with the reduced polynomials. This also
provides a separate period-one certificate for the spectator integral.

[The Haar-unitary Gram realization](../notes/cue_gram_model.md) gives a
positive random-matrix foundation for the model. Independent exact Weyl
integration checks the moments through order eight at four matrix sizes.
[The occupancy proof](../notes/cue_limit_determinacy.md) bounds moments at
every order, establishes moment determinacy and proves convergence of the
expected spectral measures to a unique limit. These are mathematical model
results beyond the individual class evaluations. Their application to zeta
zeros and priority in the literature remain separate questions.

[Connected multi-cycle counts](../notes/cue_gram_fluctuations.md) further
establish almost sure convergence of the random empirical spectrum and a
joint Gaussian limit for fixed polynomial statistics. Exact joint cumulant
certificates include a nondegenerate order-two/order-three covariance matrix,
independently checked by Weyl integration at five matrix sizes. The argument
extends the network method beyond individual expected moments.

[Logarithmic control](../notes/cue_gram_logdet.md) proves that the limiting
Gram measure has no atom at zero and bounds its small-eigenvalue mass.
The proof uses an exact expected log-determinant identity, independently
checked from the full Weyl density at five sizes. Arithmetic transport
remains essential to any zeta-zero application of this model result.

[The bandwidth theorem](../notes/cue_gram_bandwidth.md) extends these model
results to rectangular CUE Gram matrices. It distinguishes the row law
from the dimension-forced zero atom in the column law and derives an
explicit fluctuation variance across the half-bandwidth overlap threshold.
Independent exact checks compare both matrix representations. The
[literature map](../notes/literature_map.md) credits prior Gram motivation,
Vandermonde moment methods and CUE pair statistics, and specifies what
requires further novelty assessment and mathematical review.

[The global sine-process identification](../notes/sine_gram_identification.md)
proves that the CUE column law is also the almost sure growing-window
sine Gram law, with all moments converging. The proof supplies the missing
global limit comparison using interaction truncation, ergodic averages
and a uniform squared-kernel tail bound. It specifies the zero atom and
the point/frequency normalization required by the model interface.
The arithmetic zeta identities remain a separate research objective.

[The quantitative arithmetic-core proof](../notes/arithmetic_core_limit.md)
establishes the explicit prime-power model as -1/48+O(1/log T), including
coefficient asymptotics, exceptional-class control and exact geometry.
[A cutoff obstruction](../notes/fixed_cutoff_obstruction.md) identifies a
failure of the older absolute-tail strategy: every subpower cutoff,
including any fixed power of log T, leaves the leading model contribution
in its tail, so the finite-height 0.0111 charge
cannot certify that same asymptotic envelope. This directs the analytic
work toward signed estimates and actual-minus-model discrepancies.

[Prime-pair progression transport](../notes/prime_pair_progression_transport.md)
resolves the full-dyadic, one-chain divisor-family aggregation step using a
published prime-correlation theorem. Its weighted bound states the precise
consumption budget. Coupled prime chains retain separate scale and
correlation obligations.

[Short-window Fourier transport](../notes/short_window_fourier_transport.md)
derives a newer position-averaged prime-model replacement and controls
coupled products with explicit marginal concentration budgets. A separate
longer-window pointwise option is available. The work replaces an invalid
unrestricted endpoint deduction with verified ranges and a presieved main
term whose consumption remains an analytic research target.

[The finite presieved comb](../notes/presieved_pair_transport.md) identifies
an explicit divisor model's Fourier coefficients and pair main term. A
Parseval deduction gives actual short-window prime pairs averaged over
positions and shifts, with arbitrary logarithmic savings and a weighted
divisor-family budget. This advances the arithmetic transport program
without assuming independent prime chains. Prescribed individual windows,
weighted and dilated four-prime locks and the final zeta reduction remain open. Independent
review and priority assessment are still needed.

[The resolved lock audit](../notes/resolved_lock_frame.md) identifies an
information loss in a proposed Fourier reduction: separate power spectra
can agree while resolved four-point variances differ. It supplies exact
nonnegative counterexamples, the correct two-dimensional transform, a
valid full-lock stability theorem and a cubic-uniformity replacement
criterion. This narrows the missing arithmetic input and prevents an
aggregated statistic from being consumed as a resolved variance.

[Qualitative rectangle transport](../notes/averaged_rectangle_transport.md)
advances this restricted problem using published cubic uniformity. It
replaces actual four-prime counts averaged over short-window starts and
both shifts by an explicit finite local main term, with vanishing error at
the natural squared scale. An exact cube-residue calculation preserves the
normalization budget. Its qualitative rate and undilated range are stated
explicitly; the weighted zeta interface still needs further estimates.

[The rectangle Euler-tail theorem](../notes/rectangle_singular_series_tail.md)
supplies uniform absolute tail moments O_r(Y^2/w^r) on nondegenerate
shifts. Its mean-square case upgrades the preceding averaged prime
estimate to the full four-point singular-series main term. Positive
divisor moments and repeated-prime partitions avoid a pointwise tail
bound that would fail for rare primorial shifts. The remaining weighted
and dilated zeta interface retains separate proof obligations.

[The dilated lock calculation](../notes/dilated_rectangle_geometry.md)
derives the physical heterogeneous cube, exact raw overlap and finite
presieved main term. It proves a progression-norm replacement criterion
at the balanced raw scale and evaluates the exceptional-prime factors.
Their complete-period second-moment budget stays uniformly bounded for
prime-power moduli of any exponent. The remaining large-progression prime
norms and weighted consumption are explicit research targets.

[Weighted dilated transport](../notes/weighted_dilated_prime_transport.md)
supplies actual prime replacement for fixed coefficients and a sufficiently
slowly growing prime-power family. It permits coupled nonuniform starts
under an explicit marginal-density condition. A translated short-box tail
estimate identifies the full singular-series main term on nondegenerate
locks. The proof retains endpoint errors and the progression extraction
loss; its cutoff and saving are qualitative. Positive-power moduli,
prescribed windows and the final weighted zeta conversion remain open.

[Prescribed full-window transport](../notes/full_window_prime_transport.md)
gives a quantitative double-logarithmic saving at every start for windows
comparable to their positions, with a small power of log log H allowed
for arbitrary coprime dilations. It identifies the full singular-series
main term on nondegenerate locks and covers the reference's fixed-dilation
default windows. Genuine short prescribed windows and the positive-power
dilation range remain research targets.

[The complete dilation distribution](../notes/dilation_core_coverage.md)
measures that gap with an exact piecewise polynomial for every cutoff
R=exp(delta ell), 0<=delta<=1. The cube-root cutoff captures about 4.856%
of the leading model term; 99% coverage requires an exponent in
(0.859563,0.859564]. Subpower cutoffs capture none in the limit.
The [standalone review manuscript source](../papers/dilation_core_distribution.tex)
defines the model, derives the distribution from a clipped-cube volume,
states the uniform arithmetic error and provides independent reproduction
instructions. It gives reviewers a focused result with explicit analytic
limits and directs further research toward a substantial arithmetic range.

[The physical progression profile](../notes/progression_range_profile.md)
maps a coefficient range to the actual prime positions and progression
length. It derives a closed coverage formula through the square-root
profile, whose fraction is exactly 107/600, and gives exact slice
certificates for larger ranges. The manuscript includes this scale map
and its uniform model theorem. This sharpens the missing arithmetic
specification while keeping a candidate coefficient range separate from
a proved four-prime cubic estimate.

## Pilot milestones

| Period | Work | Deliverable |
|---|---|---|
| Week 1 | Review definitions, reference conventions and normalization | Reviewed specification and anchor cases |
| Weeks 2–3 | Review reciprocity, Gram realization, CUE/sine limit comparison, fluctuations and class certificates | Reviewed proofs and certificates, corrections or counterexamples, literature comparison |
| Week 4 | Review the quantitative core and cutoff obstruction; audit moment reduction | Reviewed proofs, precise transport specification and missing lemmas |
| Week 5 | Develop arithmetic-to-model transport and signed tail estimates | Proofs, counterexamples and a precise list of remaining lemmas |
| Week 6 | Evaluate further moments and covariance coefficients; study wider fluctuation tests | Further results or explicit obstacles, reproduction bundle and resource report |

## Resources and evaluation

Proposed credit allocation: $400 for computation and debugging, $300 for independent derivations and falsification, $200 for transport analysis, and $100 for documentation and reproduction.

Evaluation measures include accepted exact certificates, anchor recovery, rejected candidates, remaining proof obligations and model usage per result. AI-generated derivations and code are checked with exact arithmetic and mathematical review.

Cash support can fund a separately scoped implementation or independent-review milestone.
