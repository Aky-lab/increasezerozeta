# Research handoff

This file is the entry point for anyone continuing the project without access to the original chat history.

**Repository visibility:** keep this repository private until there is a theorem-level result worth external circulation.

## 1. Stable baseline on main

The material on main contains exact algebraic/combinatorial reductions and reproducible checks. It does not establish a new unconditional zeta-zero record by itself.

### Fourth-moment/core track

- class-subtracted l1 universality lemma;
- exact continuum-core candidate C_core = -1/48;
- explicit remainder ledger separating the fixed-P tail from true o(1) terms.

See notes/class_subtracted_universality.md, notes/arithmetic_to_continuum.md and scripts/verify_lemmas.py.

### Moment-consumption track

The moment-to-origin-mass optimisation is canonical:

lambda_n(0) = 1 / (e0^T H_n^{-1} e0).

For the exact moments through m6,

lambda_3(0) = 1415/13891,

so under the same spectral/counting interface as the reference higher-moment programme,

1 - 2 lambda_3(0) = 11061/13891 = 0.7962709668...

See notes/christoffel_tower.md, notes/k8_target_geometry.md, scripts/christoffel_exact.py and scripts/k8_target.py.

The consumption side is therefore essentially solved for future exact even moment towers. The bottleneck is producing and transporting higher moments.

## 2. External reference implementation inspected

The higher-moment work was derived and audited against:

JoshuaHKU/zeta-0.7947-reproduction

The reference revision used during the original audit was commit:

d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8

Do not silently mix conventions from another revision without checking the class-d replacement, normalisations and Bell/frozen-slot bookkeeping.

## 3. Current frontier branch

The active research branch is:

research/k8-moments

Draft PR: #6.

This branch contains candidate-level seventh/eighth-moment work.

### IMPORTANT correction: retired 1/8 candidate

An early five-dimensional midpoint computation suggested

{5,2} = 1/8.

That candidate is RETIRED.

The exact spectator-frequency elimination exposed the drift as a quadrature artefact. The active corrected numerical candidate is

{5,2} = 7/72.

Authoritative files:

- notes/spectator_52_eliminate_v.md
- results/spectator_52_coarse.md
- notes/spectator_52_candidate.md

If any older branch or script says 1/8, treat it as stale.

### Seventh moment

The exact Bell(7) ledger is

m7 = 1717/90 + {5,2} + C7.

Using the active candidate {5,2}=7/72 and the reference identification candidate C7=-17/360 gives

m7 = 3443/180 = 19.1277777...

(candidate only).

The exact Stieltjes floor implied by m0,...,m6 is

m7* = 25866469/1352400 = 19.126345016...

so the candidate is only about 0.00143 above the floor.

See notes/m7_ledger_derivation.md and scripts/m7_ledger.py.

### Eighth moment

The Bell(8) ledger reduces the new information to

A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8

with

m8 = 1103/30 + A8

under the active seventh-order candidate.

Current pre-registered numerical candidates:

{2^4} = 1661/3780

{4,2,2} = -127/840

{6,2} = -1/20

These are NOT theorem-level inputs.

See notes/pair4_candidate.md, notes/422_candidate.md, notes/62_candidate.md, notes/m8_ledger_and_80_target.md, scripts/bell8_orbits.py and scripts/m8_ledger.py.

Conditional on those three candidates, define

T8 = {4,4} + C8.

The exact target for 80% is

T8 <= 291211/11421900 = 0.02549584570...

while the moment-cone floor corresponds to approximately

T8 >= 0.02497315369.

Thus the decisive model-side k=8 problem is a narrow bound or evaluation of {4,4}+C8, after the earlier candidates are certified.

## 4. Obsolete/experimental branch warning

compute/exact52-batch was created when 1/8 was still the pre-registered {5,2} target.

It is historical/experimental and should NOT be treated as authoritative without updating it to the exact-v reduction and active 7/72 target.

The authoritative seventh-order state is on research/k8-moments.

## 5. What must be proved next

In recommended order:

1. Exact-certify {5,2} using the exact spectator-frequency elimination followed by four-dimensional rational piecewise/polytope integration.
2. Exact-certify C7=-17/360 rather than relying on numerical identification.
3. Prove the arithmetic transport for all seventh-order classes, yielding theorem-grade m7.
4. Exact-certify or rigorously enclose {2^4}, {4,2,2}, {6,2}.
5. Evaluate or rigorously upper-bound {4,4}+C8 strongly enough to meet the target in notes/m8_ledger_and_80_target.md.
6. Prove the eighth-order arithmetic transports.
7. Feed the resulting m7,m8 into the exact Christoffel engine.
8. Only then discuss an 80%+ theorem claim.

## 6. Longer-term track

Issue #5 records the path towards 90%: automate the arithmetic moment-production side so higher Bell/frozen-slot ledgers and transports are generated rather than hand-written.

The Christoffel consumption layer should stay independent and generic.

## 7. Open issues

- #2: theorem-level fixed-P tail bound for the simpler fourth-moment route.
- #4: 80% target: pin m7, control m8.
- #5: automate higher-moment arithmetic transport towards 90%.

For the current ambitious programme, prioritise #4, then #5. Issue #2 remains useful as a simpler independent theorem route.

## 8. Research discipline

Every numerical quantity must be labelled as one of:

- exact/proved algebra;
- exact computer-assisted rational evaluation;
- rigorous enclosure;
- numerical identification candidate;
- retired/falsified candidate.

Do not promote rational reconstruction from numerical data into a theorem input without an independent exact/enclosure gate.

When a falsifier fires, preserve the old value in the audit history and clearly mark it retired. The 1/8 -> 7/72 correction is the model example.

## 9. Minimal resume checklist

A new researcher should:

1. read this file;
2. run the exact checks on main;
3. checkout research/k8-moments;
4. read results/spectator_52_coarse.md before any older {5,2} note;
5. run the Bell(7)/Bell(8) gates and ledgers;
6. work on exact {5,2} certification first.

No mathematical step should require private chat context after following this handoff. If a future decision depends on something not written here or in a cited note/script, add it to the repository before relying on it.
