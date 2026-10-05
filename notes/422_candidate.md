# Pre-registered candidate for \(\{4,2,2\}\)

**Status:** numerical identification candidate, frozen before exact symbolic evaluation.

The class consists of one connected four-block and two pair blocks on the 8-cycle. The \(210\) placements reduce to \(22\) dihedral orbits.

An independent midpoint implementation uses:

- the 26-term partition-cyclic \(C_4\) evaluator;
- pair factors \(C_2(v)=|v|\) and \(C_2(w)=|w|\) on the overlap support;
- the outer 8-walk overlap;
- exact dihedral orbit multiplicities.

## Refinement ladder

| \(h\) | total |
|---:|---:|
| 0.250 | -0.2377166748046875 |
| 0.200 | -0.207316992000000 |
| 0.125 | -0.1734153628349304 |
| 0.100 | -0.1654571820000001 |

These values strongly indicate

\[
\boxed{
\{4,2,2\}=-\frac{127}{840}
=-0.151190476190476\ldots
}
\]

The error divided by \(h^2\) is:

| \(h\) | \((M_h+127/840)/h^2\) |
|---:|---:|
| 0.250 | -1.38442 |
| 0.200 | -1.40316 |
| 0.125 | -1.42239 |
| 0.100 | -1.42667 |

consistent with a stable midpoint \(O(h^2)\) expansion.

## Pre-registered falsifier

The candidate is now frozen as

\[
\boxed{\{4,2,2\}=-127/840}.
\]

It must be rejected or revised if finer grids, an independent implementation, or exact rational integration disagree.

No theorem-level \(m_8\) calculation should consume this value before exact certification.
