# Actual-walk spectator audit — 5 October 2026

**Status:** exact counterexample and elementary spectator reduction;
numerical reconnaissance for the full joint constant. No theorem claim.

**Later exact evaluation:** the subsequent integer-lattice calculation gives
the corrected model constant {5,2}=1/8. See
[the degree/period proof](../notes/spectator_52_exact_lattice.md) and
[the raw exact record](spectator_52_exact_lattice_2026-10-05.json).
The numerical table below is preserved as historical reconnaissance;
the later certificate supersedes its unresolved-value disposition.

## Finding

The historical distance-three elimination used the fixed set
`{0,s1,s2,s3,s4}`. The defining seven-cycle walk instead has
`{0,s2,s3,s4}` as its fixed set. Adding s1 changes the overlap.

At `(s1,s2,s3,s4)=(4/5,1/5,3/10,2/5)`, direct exact integration gives
`J3=2/75`; the historical formula gives `1/375`.
The corrected formula and its derivation are in
[the elimination note](../notes/spectator_52_eliminate_v.md).

The corrected checker derives the seven positions from increments
independently. It passes 7,203 rational cases and detects 358
distance-three failures of the historical formula.

## Numerical reproduction

The optional NumPy script evaluates the same signed 150-term C5
definition with corrected and historical spectator weights side by side.
The historical column reproduces the previous ladder at shared mesh
sizes, isolating the changed distance-three support.

| Mesh | Corrected actual-walk total | Historical superset total | Pure C5 anchor |
|---|---:|---:|---:|
| 1/4 | 0.134422302246 | 0.102405548096 | 0.029296875000 |
| 1/8 | 0.127610564232 | 0.098783552647 | 0.028198242188 |
| 1/10 | 0.126690506250 | 0.098242287500 | 0.028050000000 |
| 1/16 | 0.125668724999 | 0.097629533149 | 0.027885437012 |
| 1/20 | 0.125429221777 | 0.097484211768 | 0.027846875000 |
| 1/32 | 0.125168188795 | 0.097325117109 | 0.027804851532 |

The inherited pure-C5 anchor is `1/36=0.0277777777...`.
Three Richardson probes `(4*fine-coarse)/3` for the corrected total:

- 1/8 to 1/16: 0.125021445254;
- 1/16 to 1/32: 0.125001343394;
- 1/10 to 1/20: 0.125008793620.

These support renewed exploration of 1/8. There is no error enclosure,
and rational reconstruction is not proof. The pointwise counterexample
does not alone prove a value for the signed full integral.

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
contains component values, historical comparisons, anchor values,
mesh sizes, and timing. The [exact-check report](../sponsorship/evidence/actual-walk-2026-10-05.json)
records executed hashes, including the reduction module.

## Disposition

- Withdraw the historical explanation that 1/8 was a spectator
  quadrature artifact.
- Quarantine the historical 7/72 ladder and downstream scenario values.
- Keep neither fraction as an exact-certified input.
- Recompute downstream m7/m8 geometry only after the actual joint
  constant and C7 have been certified.
- Preserve prior notes and raw data with explicit superseded warnings.

The earlier sponsorship checker was also affected: it independently
integrated the overlap of the erroneous superset. Passing that test
did not establish agreement with the defining W3 walk. The new
increment-based checker avoids that shared structural assumption.
