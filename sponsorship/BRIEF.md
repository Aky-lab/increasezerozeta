# Research support proposal

## Project

Develop auditable exact computations for higher moments in research on simple zeros of the Riemann zeta function. A six-week pilot seeks **$1,000 in AI API credits** for independent derivations, computation and analytic proof development.

The immediate objective is a reproducible route from model definitions to exact integrals and clearly specified transport lemmas.

## Existing results

The corrected four-dimensional model integral {5,2} is exactly 1/8. Its components are 5/504, 1/360 and 13/2520. The [proof](../notes/spectator_52_exact_lattice.md) bounds the degree and period of the weighted lattice sums; the [integer record](../results/spectator_52_exact_lattice_2026-10-05.json) includes held-out checks and source hashes. A separate C5 anchor evaluates to 1/36.

The [pure-cycle calculation](../notes/pure_cycle_flow_polytopes.md) gives C7=-17/360 exactly using integral network-flow polytopes. Independent checks recover all 150 C5 term values from the reference certificates and twelve C6 term values. The two seventh-order model inputs together give m7=862/45.

The repository provides a standard-library engine for exact Christoffel/Hankel moment consumption and eighth-order target geometry. An independent spectator checker passes 7,203 rational cases and detects 358 failures of an earlier distance-three formula. The [audit](../results/spectator_actual_walk_audit.md) identifies the missing-prefix error.

The six-moment algebra yields lambda_3(0)=1415/13891. Its zeta interpretation depends on the analytic moment framework. Eighth-order class evaluation and higher-order arithmetic transport remain open; the model calculations alone do not establish a new unconditional bound.

See [the research overview](../RESEARCH.md) for conventions and proof obligations.

## Pilot milestones

| Period | Work | Deliverable |
|---|---|---|
| Week 1 | Review definitions, reference conventions and normalization | Reviewed specification and anchor cases |
| Weeks 2–3 | Independently review the seventh-order certificates and extend the method to eighth-order classes | Certificate review and exact class values or documented obstructions |
| Week 4 | Audit remaining higher-moment candidates | Exact values, enclosures or counterexamples |
| Week 5 | Develop arithmetic-to-model transport | Proofs and a precise list of remaining lemmas |
| Week 6 | Package computations and evaluate the next stage | Reproduction bundle, resource report and mathematical review |

## Resources and evaluation

Proposed credit allocation: $400 for computation and debugging, $300 for independent derivations and falsification, $200 for transport analysis, and $100 for documentation and reproduction.

Evaluation measures include accepted exact certificates, anchor recovery, rejected candidates, remaining proof obligations and model usage per result. AI-generated derivations and code are checked with exact arithmetic and mathematical review.

Cash support can fund a separately scoped implementation or independent-review milestone.
