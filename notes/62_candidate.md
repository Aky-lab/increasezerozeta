# Pre-registered candidate for \(\{6,2\}\)

**Status:** coarse numerical identification candidate, frozen before finer-grid or exact evaluation.

The class consists of one connected six-block and one spectator pair on the 8-cycle.

The 28 pair placements split into four cyclic-distance types with multiplicities

\[
8,\;8,\;8,\;4
\]

for distances \(1,2,3,4\).

An independent midpoint implementation uses:

- the 1082-term partition-cyclic \(C_6\) evaluator;
- pair factor \(C_2(v)=|v|\) on the outer-overlap support;
- the four canonical spectator placements.

## Coarse refinement

| \(h\) | \(\{6,2\}\) |
|---:|---:|
| 0.250000 | -0.03594970703125 |
| 0.200000 | -0.04038451200000 |
| 0.166667 | -0.04309040479094 |

A \(1,h^2,h^4\) extrapolation lands at approximately

\[
-0.0500174,
\]

suggesting the simple exact value

\[
\boxed{
\{6,2\}=-\frac1{20}.
}
\]

The errors relative to \(-1/20\), divided by \(h^2\), are approximately

\[
0.2248,\qquad0.2404,\qquad0.2487,
\]

consistent with an \(O(h^2)\) midpoint expansion.

## Pre-registered falsifier

The candidate is frozen as

\[
\boxed{\{6,2\}=-1/20}.
\]

It must be revised if finer refinement, an independent engine, or exact rational integration disagrees.

No theorem-level \(m_8\) input should consume this value before the exact gate passes.
