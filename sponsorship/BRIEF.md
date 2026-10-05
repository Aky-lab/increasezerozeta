# Research support proposal

## Objective

Support a six-month investigation of the arithmetic estimate needed to turn an exact continuum moment model into a stronger bound on simple critical-line zeros of the Riemann zeta function. The plan begins with a six-week pilot. The requested budget is **$5,000 in AI credits**, with research subscription access where available. Deterministic computations run locally.

The immediate target is a signed trace estimate for the actual prime operator. Exact finite-model certificates and a conditional counting implication are available; the central arithmetic estimate remains open.

## Evidence

| Material | Available evidence |
|---|---|
| [Exact model certificates](../notes/eighth_order_certificate.md) | Fourteen class evaluations through eighth order, integer records, justified interpolation bounds and independent checks |
| [Model and source guide](../RESEARCH.md) | Definitions, exact moment ledger, Gram-model proof notes and primary references |
| [Counting bridge](../notes/spectral_counting_bridge.md) | A sufficient actual-matrix trace hypothesis for the 80% target |
| [Bounded resolvent reduction](../notes/bounded_resolvent_bridge.md) | An alternative target with controlled compression error and three resolvent powers |
| [Prime-removal expansion](../notes/prime_removal_resolvent.md) | Local remainder bounds and explicit conditional-phase/covariance terms to estimate |
| [Positive-weight resolvent criterion](../notes/mirror_resolvent_certificate.md) | A bounded certificate with three first resolvents and a sufficient single-point estimate |
| [Lean formalization](../formal/README.md) | 26 kernel-checked scalar and conditional finite counting theorems, with an axiom audit |
| [Local verification record](../results/project_verification_2026-10-05.json) | 33 Python checks and the separate Lean proof audit, including source hashes |

The eighth-order model certificate has origin-mass value `12241115/162540559`, corresponding to a conditional counting conversion of approximately 84.938%. It does not establish a new unconditional zeta-zero bound. Analytic proof notes require independent mathematical review.

## Six-week pilot

| Period | Work | Deliverable |
|---|---|---|
| Week 1 | Review actual-matrix definitions, the counting bridge and formal assumptions | A precise statement of the arithmetic target and normalization checks |
| Week 2 | Check the bounded certificate, compression estimate and prime-removal expansion through separate derivations | Reviewed identities and reproducible checks |
| Weeks 3–4 | Investigate conditional prime-phase cancellation and the nonlinear covariance trace | Proved lemmas, rigorous counterexamples or a quantified obstruction |
| Week 5 | Extend formal verification where it clarifies the analytic interface; assemble material for expert review | Formal lemmas and a focused review manuscript |
| Week 6 | Reproduce the evidence and assess further research | Reproduction bundle, remaining proof obligations and usage report |

Months 2–4 pursue the most promising estimates and address mathematical review findings. Months 5–6 consolidate proved results, formal statements and the final reproduction report. Scope is reassessed after the pilot.

A documented obstruction is a useful outcome. Success is measured by the mathematical statements resolved and the quality of their verification; a new zeta theorem is a research objective, not a promised deliverable.

## Use of AI and resources

AI will assist with proof development, alternative derivations, counterexample search and Lean implementation. Exact rational checks, separate constructions and kernel verification test accepted finite results. Human mathematical review is sought for analytic arguments and their hypotheses.

Proposed credit allocation: $2,000 for analytic derivations, $1,500 for adversarial mathematical checks, $1,000 for formalization and computational implementation, and $500 for reproduction and reporting. An initial $1,000 tranche can fund the six-week pilot, with further allocation tied to its results. The report will record usage and cost per accepted result where the service exposes those measurements.

Support for independent mathematical review can also be scoped as a separate cash-funded milestone.
