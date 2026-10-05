# Numerical comparison for the extra-prefix spectator formula

The table below records an extra-prefix formula with the distance-three fixed set {0,s1,s2,s3,s4}. That set contains an extra s1 and evaluates a different integral from the defining walk.

The model value is {5,2}=1/8 by [exact lattice evaluation](../notes/spectator_52_exact_lattice.md). The [actual-walk audit](spectator_actual_walk_audit.md) compares defining and extra-prefix calculations.

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

## Extra-prefix ladder

The spectator frequency was integrated analytically, while the four remaining variables used midpoint quadrature. The distance-three weight retained the extra prefix.

| dc | U1 | U2 | U3 | {5,2}=7(U1+U2+U3) |
|---:|---:|---:|---:|---:|
| 0.2500 | 0.0103759766 | 0.0029830933 | 0.0012702942 | 0.1024055481 |
| 0.2000 | 0.0102400000 | 0.0029200000 | 0.0012489600 | 0.1008627200 |
| 0.1250 | 0.0100574493 | 0.0028380156 | 0.0012164712 | 0.0987835526 |
| 0.1000 | 0.0100100000 | 0.0028170312 | 0.0012075812 | 0.0982422875 |
| 0.0800 | 0.0098889730 | 0.0027761171 | 0.0011957809 | 0.0970260964 |
| 0.0625 | 0.0099563102 | 0.0027934096 | 0.0011973564 | 0.0976295331 |
| 0.0500 | 0.0099435807 | 0.0027878263 | 0.0011949089 | 0.0974842118 |

These are numerical quadratures without error enclosures. Their apparent convergence toward 7/72 concerns the extra-prefix formula rather than the defining model.
