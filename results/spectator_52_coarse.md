# Coarse {5,2} midpoint reconnaissance

**Status: numerical model-side reconnaissance only. Not consumed as proof.**

These values were produced with the independent definition in
`notes/spectator_52_spec.md`.  The implementation uses:

- 150 partition-cyclic C5 terms;
- 15 cached unique prefix linear forms;
- global negation symmetry in the spectator frequency;
- three pair-placement distances with multiplicities 7/7/7.

## Calibration: pure C5

The same C5 evaluator, integrated against the pure five-cycle overlap,
should converge to the already certified value

[
C_5=rac1{36}=0.027777777777ldots.
]

Observed coarse midpoint values:

| dv | pure C5 |
|---:|---:|
| 0.500 | 0.0312500000 |
| 0.250 | 0.0292968750 |
| 0.200 | 0.0288000000 |
| 0.125 | 0.0281982422 |
| 0.100 | 0.0280500000 |

The convergence direction and scale are consistent with the certified
anchor.

## {5,2} ladder

[
{5,2}=7(U_1+U_2+U_3).
]

| dv | U1 | U2 | U3 | {5,2} |
|---:|---:|---:|---:|---:|
| 0.500 | 0.0117187500 | 0.0029296875 | 0.0097656250 | 0.1708984375 |
| 0.250 | 0.0108032227 | 0.0028228760 | 0.0064392090 | 0.1404571533 |
| 0.200 | 0.0105216000 | 0.0028064000 | 0.0059840000 | 0.1351840000 |
| 0.125 | 0.0101709366 | 0.0027887821 | 0.0054831505 | 0.1291000843 |
| 0.100 | 0.0100831500 | 0.0027847875 | 0.0053666250 | 0.1276419375 |
| 0.080 | 0.0098569944 | 0.0028121542 | 0.0051747091 | 0.1249070039 |

The total is not monotone at these coarse grids, so no rational
identification is justified yet.

Two dyadic-style Richardson probes happen to land near (1/8):

- from (dv=0.25,0.125): (0.1253143946);
- from (dv=0.20,0.10): (0.1251279167).

This makes (1/8) a **candidate worth falsifying**, not a result.
The (dv=0.16,0.08) pair behaves much less cleanly, which is a useful
warning against premature reconstruction.

## Next gates

1. produce finer rungs with the cached-prefix implementation;
2. separate the three (U_d) convergence profiles;
3. verify exact negation symmetry by paired slices;
4. implement support pruning and confirm bitwise agreement with the
   unpruned coarse engine;
5. only after stable (h^2)-type convergence, attempt rational
   reconstruction;
6. independently certify any candidate with exact polytope integration.
