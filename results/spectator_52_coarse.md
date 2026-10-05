# {5,2} numerical reconnaissance

**Status: model-side falsifier only. Nothing here is consumed as proof.**

The first five-dimensional midpoint ladder discretised the spectator
frequency v and showed a strong but misleading drift towards 1/8.
The exact elimination of v in
`notes/spectator_52_eliminate_v.md` revealed that this was a
quadrature artefact. The 1/8 candidate is therefore **retired**.

## Calibration: pure C5

The same partition-cyclic five-point cumulant evaluator should recover

    C5 = 1/36 = 0.027777777777...

Observed midpoint values:

| dc | pure C5 |
|---:|---:|
| 0.500 | 0.0312500000 |
| 0.250 | 0.0292968750 |
| 0.200 | 0.0288000000 |
| 0.125 | 0.0281982422 |
| 0.100 | 0.0280500000 |

This is consistent with the certified anchor.

## Corrected ladder: spectator frequency integrated exactly

For each four-dimensional midpoint cell, the v-integral is evaluated
analytically using the cubic formula in
`notes/spectator_52_eliminate_v.md`.

| dc | U1 | U2 | U3 | {5,2}=7(U1+U2+U3) |
|---:|---:|---:|---:|---:|
| 0.2500 | 0.0103759766 | 0.0029830933 | 0.0012702942 | 0.1024055481 |
| 0.2000 | 0.0102400000 | 0.0029200000 | 0.0012489600 | 0.1008627200 |
| 0.1250 | 0.0100574493 | 0.0028380156 | 0.0012164712 | 0.0987835526 |
| 0.1000 | 0.0100100000 | 0.0028170312 | 0.0012075812 | 0.0982422875 |
| 0.0800 | 0.0098889730 | 0.0027761171 | 0.0011957809 | 0.0970260964 |
| 0.0625 | 0.0099563102 | 0.0027934096 | 0.0011973564 | 0.0976295331 |
| 0.0500 | 0.0099435807 | 0.0027878263 | 0.0011949089 | 0.0974842118 |

The remaining four-dimensional midpoint rule still crosses rational
kink hyperplanes, so non-monotonicity is expected and these rungs
must not be treated as enclosures.

Two genuinely dyadic Richardson probes are:

- dc 0.125 -> 0.0625: 0.09724485998;
- dc 0.100 -> 0.0500: 0.09723151986.

Both are close to

    7/72 = 0.09722222222...

so **7/72 is now the active candidate to falsify**, not an accepted
identification.

The deviations of the two Richardson probes from 7/72 are about
2.26e-5 and 9.30e-6 respectively. That is encouraging but nowhere
near enough for a proof.

## Next gate

The exact-v reduction turns each partition-cyclic term into a
four-dimensional piecewise polynomial of degree at most four on a
rational hyperplane arrangement. The next decisive step is exact
rational cell/polytope integration. A rational reconstruction is not
to be consumed until that independent exact computation agrees.
