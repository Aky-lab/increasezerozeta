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

The assembled model moments are m6=640/63, m7=3439/180 and m8=747361/20160. Exact degree-four Hankel algebra gives a conditional counting conversion of approximately 84.938%. The finite-model local identities are proved, while their arithmetic transport and the spectral/counting interface remain open. These calculations do not establish a new unconditional zeta-zero bound.

See [the research overview](../RESEARCH.md) for conventions and proof obligations.

## Structural research and computational value

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

## Pilot milestones

| Period | Work | Deliverable |
|---|---|---|
| Week 1 | Review definitions, reference conventions and normalization | Reviewed specification and anchor cases |
| Weeks 2–3 | Review reciprocity, Gram realization, determinacy and class certificates | Reviewed proofs and certificates, corrections or counterexamples, literature comparison |
| Week 4 | Audit the analytic moment-production framework | Precise transport specification and missing lemmas |
| Week 5 | Develop arithmetic-to-model transport | Proofs and a precise list of remaining lemmas |
| Week 6 | Evaluate further moments and spectral fluctuations; package computations | Further results or explicit obstacles, reproduction bundle and resource report |

## Resources and evaluation

Proposed credit allocation: $400 for computation and debugging, $300 for independent derivations and falsification, $200 for transport analysis, and $100 for documentation and reproduction.

Evaluation measures include accepted exact certificates, anchor recovery, rejected candidates, remaining proof obligations and model usage per result. AI-generated derivations and code are checked with exact arithmetic and mathematical review.

Cash support can fund a separately scoped implementation or independent-review milestone.
