# Research sponsorship brief

Prepared 5 October 2026. Draft for private review; no application has been submitted and no sponsor has committed support.

## Proposal

Support a six-week pilot to turn numerical higher-moment candidates in zeta-function research into auditable exact computations or explicit counterexamples. Proposed initial request: **$1,000 in AI API credits**. This is a requested pilot budget, not a quotation or a promised award.

The scientific goal is to investigate improved lower bounds on the proportion of simple zeros of the Riemann zeta function on the critical line. The nearer deliverable is a reusable workflow for exact rational integration, candidate falsification, and traceable mathematical verification.

## What exists

The repository has a standard-library Python engine for exact Christoffel/Hankel moment consumption, exact target geometry for the seventh and eighth moments, and a written ledger of unresolved analytic obligations. The active frontier is `research/k8-moments`, draft PR #6. Read [HANDOFF.md](../HANDOFF.md) for the research state.

The existing six-moment algebra gives `lambda_3(0) = 1415/13891`. Its zeta interpretation depends on the surrounding candidate analytic framework. It is **not an independently established unconditional record**. The 80% target remains a research target.

A later actual-walk audit found a structural error in the spectator elimination used to prefer `7/72` over `1/8`. That preference is quarantined. The corrected model integral is now evaluated exactly as **{5,2}=1/8** by integer lattice sums with proved polynomial-degree and period bounds. Its components are 5/504, 1/360 and 13/2520; the independent C5 anchor is 1/36. All held-out checks pass. The [proof and reproduction note](../notes/spectator_52_exact_lattice.md) and [raw integer record](../results/spectator_52_exact_lattice_2026-10-05.json) provide concrete review materials. C7 and the arithmetic transport remain unresolved.

## Evidence produced in this preparation

All three existing main-branch check scripts passed locally. The corrected independent checker constructs the actual walks from increments, integrates their piecewise-linear overlap, and compares with the corrected formula over 7,203 rational cases; 2,539 have nonzero integrals. All corrected comparisons agree exactly; the historical formula fails 358 distance-three cases. An explicit counterexample gives 2/75 instead of 1/375.

The standard-library checks cover finite algebra and the spectator reduction. The separate NumPy lattice certificate evaluates the stated model integral exactly; neither layer proves arithmetic transport or a new zeta theorem.

Reproduce with Python 3.12 or later, without third-party packages:

```sh
python scripts/sponsorship_evidence.py --output sponsorship/evidence/local.json
```

The [current evidence](evidence/actual-walk-2026-10-05.json) includes script hashes, outputs, Python version, return codes, and limitations. It identifies executed files; it does not claim a complete upstream checkout. The historical frontier formula was inspected at `b90e747da2fed0395402396d5b755bdd3a21ad05`; see the actual-walk audit for its correction. The earlier evidence report is retained with a superseded-scope warning.

## Six-week pilot and acceptance gates

| Period | Work | Reviewable output |
|---|---|---|
| Week 1 | Freeze definitions, conventions, and anchor cases | Versioned specification and normalization ledger |
| Weeks 2–3 | Review the completed `{5,2}` certificate and extend exact computation to C7 | Independent certificate review and C7 computation or obstruction |
| Week 4 | Audit normalizations and the remaining candidate history | Exact results/enclosures, or documented candidate rejections |
| Week 5 | Review the arithmetic-to-model transport | Explicit proof obligations with resolved/unresolved status |
| Week 6 | Package evidence and assess the next stage | Reproduction bundle, resource report, and review memorandum |

These are proposed milestones. The pilot can succeed by falsifying a candidate or documenting an obstruction. No target theorem is promised. Independent review must be arranged; no reviewer has yet committed.

Proposed credit allocation: $400 for integration code and debugging, $300 for independent derivations and adversarial checks, $200 for transport analysis, $100 for documentation and reruns. Record actual model, token usage, cost, attempts, and accepted certificates. Only request a larger allocation after demonstrating useful results.

## Sponsor value

The pilot would provide a concrete account of where AI helps research mathematics and where exact verification rejects plausible output. Proposed reports measure cost per accepted certificate, anchor recovery, failed candidates, reproducibility, and remaining proof gaps. Such measurements are a pilot deliverable; they have not yet been collected.

An alternative cash proposal is a separately negotiated $2,000 pilot for implementation time and independent mathematical review. This is an estimate requiring a scoped agreement and named reviewer, not a request already sent or a funded commitment.

## Disclosure and ownership

The repository remains private and currently has no license. External review should use a maintainer-approved summary or controlled access. Public release, licensing, acknowledgment wording, and any case study require a separate decision. Funding would not establish authorship, endorsement, or correctness. The maintainer retains ownership and research judgment.
