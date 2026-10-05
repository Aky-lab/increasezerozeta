# Spectator integral: exact checks and numerical comparison

The model integral is exactly {5,2}=1/8; see [the lattice proof](../notes/spectator_52_exact_lattice.md) and [integer record](spectator_52_exact_lattice_2026-10-05.json). The comparison below distinguishes the defining walk from a walk with an additional prefix.

## Finding

The defining seven-cycle walk has `{0,s2,s3,s4}` as its distance-three
fixed set. Enlarging it to `{0,s1,s2,s3,s4}` changes the overlap.

At `(s1,s2,s3,s4)=(4/5,1/5,3/10,2/5)`, direct exact integration gives
`J3=2/75`; the extra-prefix formula gives `1/375`.
The defining formula and its derivation are in
[the elimination note](../notes/spectator_52_eliminate_v.md).

The defining checker derives the seven positions from increments
independently. It passes 7,203 rational cases and detects 358
distance-three failures of the extra-prefix formula.

## Numerical reproduction

The optional NumPy script evaluates the same signed 150-term C5
definition with defining and extra-prefix spectator weights side by side.
The extra-prefix column measures the effect of enlarging the distance-three support.

| Mesh | Defining actual-walk total | Extra-prefix superset total | Pure C5 anchor |
|---|---:|---:|---:|
| 1/4 | 0.134422302246 | 0.102405548096 | 0.029296875000 |
| 1/8 | 0.127610564232 | 0.098783552647 | 0.028198242188 |
| 1/10 | 0.126690506250 | 0.098242287500 | 0.028050000000 |
| 1/16 | 0.125668724999 | 0.097629533149 | 0.027885437012 |
| 1/20 | 0.125429221777 | 0.097484211768 | 0.027846875000 |
| 1/32 | 0.125168188795 | 0.097325117109 | 0.027804851532 |

The inherited pure-C5 anchor is `1/36=0.0277777777...`.
Three Richardson probes `(4*fine-coarse)/3` for the defining total:

- 1/8 to 1/16: 0.125021445254;
- 1/16 to 1/32: 0.125001343394;
- 1/10 to 1/20: 0.125008793620.

These numerical probes approach the exact value 1/8 but provide no error enclosure. The proof of the full signed integral comes from the separate exact lattice calculation.

## Reproduction and scope

With NumPy available:

```sh
python scripts/spectator_52_reduced_midpoint.py 1/4 1/8 1/10 1/16 1/20 1/32 --out results/local-midpoint.json
```

NumPy 2.3.5 was used for the recorded run. The fine grid can take
minutes; it is excluded from the standard-library CI suite.
The script also checks 300 floating spectator weights against
the rational implementation. That is a finite numerical gate.

Every free ci is bounded in absolute value by one on outer support:
each is a difference between consecutive positions within either
the shifted or the fixed portion of its walk. Thus [-1,1]^4 covers
the outer support for all three distances. Boundary cells have
zero overlap. The midpoint rule remains a numerical quadrature.

The [raw numerical record](spectator_actual_walk_2026-10-05.json)
contains component values, extra-prefix comparisons, anchor values,
mesh sizes, and timing. The [exact-check report](../sponsorship/evidence/actual-walk-2026-10-05.json)
records executed hashes, including the reduction module.
