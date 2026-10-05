# Pre-registered candidate for the four-pair class \(\{2^4\}\)

**Status:** numerical identification candidate, frozen before exact symbolic evaluation.

The class \(\{2^4\}\) is the sum over all \(105\) perfect matchings of eight cyclic positions. Each matching carries four pair frequencies \(v_1,\dots,v_4\), outer overlap

\[
(1-\operatorname{spread}(\text{prefix walk}))_+,
\]

and weight

\[
|v_1v_2v_3v_4|.
\]

A direct midpoint implementation reproduces the known six-cycle class \(\{2,2,2\}=131/420\) under refinement, providing the calibration.

## Dihedral reduction

The \(105\) matchings split into only \(17\) orbits under the dihedral group of the 8-cycle, with orbit sizes

\[
1,2,2,4,4,4,4,4,8,8,8,8,8,8,8,8,16.
\]

Evaluating one representative per orbit gives exactly the same midpoint totals as summing all \(105\) matchings.

## Refinement ladder

The orbit-reduced midpoint values are:

| \(h\) | total |
|---:|---:|
| 0.1000 | 0.4540814670000006 |
| 0.0625 | 0.4451579764718190 |
| 0.0500 | 0.4430933711601549 |
| 0.0400 | 0.4417709679358756 |

These strongly indicate

\[
\boxed{
\{2^4\}=\frac{1661}{3780}
=0.439417989417989\ldots
}
\]

The error profile against this fraction is:

| \(h\) | \((M_h-1661/3780)/h^2\) |
|---:|---:|
| 0.1000 | 1.46635 |
| 0.0625 | 1.46944 |
| 0.0500 | 1.47015 |
| 0.0400 | 1.47061 |

which is consistent with a stable midpoint \(O(h^2)\) expansion.

## Pre-registered falsifier

The candidate is now frozen as

\[
\boxed{\{2^4\}=\frac{1661}{3780}}.
\]

It is to be rejected or revised if:

1. finer orbit-reduced midpoint rungs stop converging to this fraction with the expected even-power profile;
2. a direct all-105 implementation disagrees with the orbit-reduced computation;
3. exact rational integration produces a different value;
4. the arithmetic \(m_8\) ledger uses a different four-pair object from the model class defined here.

No theorem-level \(m_8\) bound should consume this value before the exact gate passes.
