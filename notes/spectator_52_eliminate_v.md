# Spectator-frequency elimination for {5,2}

**Status:** elementary exact reduction of the spectator integral only.
The full model integral is evaluated in [the lattice certificate](spectator_52_exact_lattice.md). Arithmetic transport remains unresolved.

## 1. Prefix sets

From the defining increments, with `sj=c1+...+cj`, the prefix sets are:

| Distance | Fixed set B_d | Shifted set A_d |
|---|---|---|
| 1 | {0,s1,s2,s3,s4} | {0} |
| 2 | {0,s1,s2,s3,s4} | {0,s1} |
| 3 | {0,s2,s3,s4} | {0,s1,s2} |

Thus `W_d(v)=B_d union (v+A_d)`, ignoring repeated points.
In distance three, A_d need not be a subset of B_d.

## 2. Exact counterexample

Take `(s1,s2,s3,s4)=(4/5,1/5,3/10,2/5)`.
At `v=-2/5`, the actual overlap is `1/5`; the extra-prefix
superset overlap is zero. Integration over v gives:

- actual distance-three integral: `2/75`;
- extra-prefix cubic: `1/375`;
- discrepancy: `3/125`.

This rules out using the same fixed set for all three distances. It does not by itself
evaluate the signed C5-weighted integral: pointwise errors could
cancel after multiplication by C5.

## 3. General reduction without a nesting assumption

Write `m=min B_d`, `M=max B_d`, `a=min A_d`, `A=max A_d`.
If `max(M-m,A-a)>=1`, the overlap is identically zero.

Otherwise define:

```text
l = M-a-1
u = 1+m-A
b = m-a
c = M-A
r = min(b,c)
t = max(b,c)
h = 1-max(M-m,A-a)
```

Then `l<r<=t<u`. The overlap is `v-l` on [l,r],
`h` on [r,t], `u-v` on [t,u], and zero outside.
This follows by checking the switches in
`max(M,v+A)-min(m,v+a)`; it works for both possible orders
of b and c, without assuming either one straddles zero.

Both prefix sets contain 0, so nonzero overlap implies `|v|<1`.
Consequently `C2(v)=min(|v|,1)=|v|` throughout the support.

The distributional second derivative of the overlap is
`delta_l-delta_b-delta_c+delta_u`.
Since `(|v|^3/6)''=|v|`, integration by parts twice, with
vanishing boundary terms for the compactly supported overlap, yields:

```text
J_d = (|l|^3 + |u|^3 - |b|^3 - |c|^3) / 6.
```

This is implemented with rational arithmetic in
`scripts/spectator_reduction.py`. It is piecewise cubic with rational
coefficients. The remaining definition stays
`U_d=integral J_d*C5 dc1...dc4` and `{5,2}=7*(U1+U2+U3)`.

## 4. Independent finite gates

`scripts/verify_spectator_reduction.py` constructs the walks from
increments independently and integrates the actual piecewise overlap.
On a seven-point rational grid for each of s1,...,s4:

- 7,203 distance/grid cases agree with the reduction;
- 2,539 have nonzero integrals;
- the extra-prefix distance-three cubic fails 358 cases;
- the explicit counterexample is checked exactly.

## 5. Full model integral

The reduction gives {5,2}=1/8 by [exact lattice evaluation](spectator_52_exact_lattice.md). The separate [flow-polytope calculation](pure_cycle_flow_polytopes.md) gives C7=-17/360. Arithmetic transport remains open.
