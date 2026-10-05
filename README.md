# increasezerozeta

Private research repository for rigorous improvements to unconditional lower bounds on the proportion of simple zeros of the Riemann zeta function on the critical line.

## Start here

Read **HANDOFF.md** first. It is the authoritative continuation guide and records:

- what is exact versus candidate-level;
- the active research branch;
- retired values that must not be reused;
- the current 80% target;
- the exact next mathematical steps.

## Stable baseline on main

The merged baseline contains:

- the class-subtracted universality/core reduction, including the candidate exact continuum value C_core = -1/48;
- the exact Christoffel/Hankel consumption framework;
- the exact six-moment origin-mass value lambda_3(0)=1415/13891, giving 11061/13891 = 79.62709668...% within the same candidate analytic framework as the reference higher-moment programme;
- exact k=8 target geometry in terms of m7 and m8;
- reproducible exact-arithmetic checks.

## Active frontier

The active higher-moment work is on:

research/k8-moments

with draft PR #6.

Important: the early {5,2}=1/8 midpoint candidate is retired. The active corrected candidate after exact spectator-frequency elimination is {5,2}=7/72. See HANDOFF.md before using any research-branch number.

## Research support

The [sponsorship brief](sponsorship/BRIEF.md) proposes a modest pilot to certify or refute higher-moment candidates. The [funding shortlist](sponsorship/PROGRAMS.md) records program fit and unresolved eligibility; the [application draft](sponsorship/APPLICATION.md) has not been submitted.

Run `python scripts/sponsorship_evidence.py` for the baseline checks and an independent finite check of the spectator-frequency reduction. Passing these checks does not establish the surrounding analytic transport or a new zeta theorem.

## Mathematical status

This is research work in progress. Nothing here should be advertised as an established new unconditional record until the analytic moment-production/transport chain has been independently checked and all candidate constants consumed by a headline result have theorem-level certification.

Keep the repository private until there is something concrete enough to circulate.

