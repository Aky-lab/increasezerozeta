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

### IMPORTANT: actual-walk audit supersedes the candidate history

The historical spectator elimination added an absent s1 prefix to the distance-three walk. An exact counterexample gives J3=2/75 against the historical 1/375. The corrected reduction and regression gates are on main in notes/spectator_52_eliminate_v.md and scripts/verify_spectator_reduction.py.

The earlier 7/72 preference is quarantined. The stated reason for retiring 1/8 is withdrawn. Corrected numerical refinement approaches 1/8, but neither fraction is certified. Read results/spectator_actual_walk_audit.md before older candidate notes.

### Seventh moment

The exact Bell(7) ledger is

m7 = 1717/90 + {5,2} + C7.

For historical scenario analysis only, substituting {5,2}=7/72 and C7=-17/360 gives

m7 = 3443/180 = 19.1277777...

(quarantined historical scenario only).

The exact Stieltjes floor implied by m0,...,m6 is

m7* = 25866469/1352400 = 19.126345016...

The historical scenario is about 0.00143 above the floor; this does not establish the actual seventh moment.

See notes/m7_ledger_derivation.md and scripts/m7_ledger.py.

### Eighth moment

The Bell(8) ledger reduces the new information to

A8 = {2^4} + {4,2,2} + {4,4} + {6,2} + C8

with

m8 = 1103/30 + A8

under the quarantined 7/72 seventh-order scenario only.

Current pre-registered numerical candidates:

{2^4} = 1661/3780

{4,2,2} = -127/840

{6,2} = -1/20

These are NOT theorem-level inputs.

See notes/pair4_candidate.md, notes/422_candidate.md, notes/62_candidate.md, notes/m8_ledger_and_80_target.md, scripts/bell8_orbits.py and scripts/m8_ledger.py.

Conditional on those three candidates, define

T8 = {4,4} + C8.

For that historical scenario, the exact target for 80% is

T8 <= 291211/11421900 = 0.02549584570...

while the moment-cone floor corresponds to approximately

T8 >= 0.02497315369.

These scenario-specific bounds must be recomputed after the actual seventh-order inputs are certified. The generic moment-space geometry in notes/k8_target_geometry.md remains independent of this correction.

## 4. Obsolete/experimental branch warning

compute/exact52-batch was created when 1/8 was still the pre-registered {5,2} target.

It is historical/experimental. Both this branch and research/k8-moments must be interpreted through the actual-walk correction; no fraction is currently certified.

The latest correction is recorded on main; older research-branch candidate notes are superseded until repaired.

## 5. What must be proved next

In recommended order:

1. Use the corrected actual-walk spectator elimination, then exact-certify {5,2} by rational piecewise/polytope integration, without baking either candidate into the result.
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

When a falsifier fires, preserve the old value in the audit history and clearly mark it retired. The later actual-walk audit overturned the evidence for the 1/8 -> 7/72 change; preserve both stages of that history.

## 9. Minimal resume checklist

A new researcher should:

1. read this file;
2. run the exact checks on main;
3. checkout research/k8-moments;
4. read results/spectator_actual_walk_audit.md and the corrected elimination note before any older {5,2} note;
5. run the Bell(7)/Bell(8) gates and ledgers;
6. work on exact {5,2} certification first.

No mathematical step should require private chat context after following this handoff. If a future decision depends on something not written here or in a cited note/script, add it to the repository before relying on it.

