# Pre-registered candidate for the seventh-moment joint class \(\{5,2\}\)

**Status:** numerical identification candidate only. This is deliberately recorded before any finer-grid run or exact polytope integration.

The independent midpoint engine defined in notes/spectator_52_spec.md gives the following coarse rungs.

| \(dv\) | \(U_1\) | \(U_2\) | \(U_3\) | \(7(U_1+U_2+U_3)\) |
|---:|---:|---:|---:|---:|
| 0.25 | 0.01080322 | 0.00282288 | 0.00643921 | 0.1404571533 |
| 0.20 | 0.01052160 | 0.00280640 | 0.00598400 | 0.1351840000 |
| 0.125 | 0.01017094 | 0.00278878 | 0.00548315 | 0.1291000843 |

A three-rung fit in \(1,dv^2,dv^4\) extrapolates the total to approximately

\[
0.1249970,
\]

suggesting the exact value

\[
\boxed{\{5,2\}=\frac18}.
\]

More strongly, the three placement components independently point to

\[
\boxed{
U_1=\frac5{504},\qquad
U_2=\frac1{360},\qquad
U_3=\frac{13}{2520}.
}
\]

These fractions satisfy

\[
\frac5{504}+\frac1{360}+\frac{13}{2520}
=
\frac1{56},
\]

and therefore

\[
7(U_1+U_2+U_3)=\frac18.
\]

## Error-profile gate

For the candidate component fractions, the observed errors divided by \(dv^2\) are:

| \(dv\) | \(e_1/dv^2\) | \(e_2/dv^2\) | \(e_3/dv^2\) |
|---:|---:|---:|---:|
| 0.25 | 0.01412 | 0.000722 | 0.02049 |
| 0.20 | 0.01502 | 0.000716 | 0.02063 |
| 0.125 | 0.01602 | 0.000704 | 0.02076 |

This is consistent with the expected midpoint \(O(dv^2)\) convergence and is substantially stronger evidence than identifying only the total.

## Pre-registered falsifier

The candidate is now frozen as

\[
(U_1,U_2,U_3)
=
\left(
\frac5{504},
\frac1{360},
\frac{13}{2520}
\right),
\qquad
\{5,2\}=\frac18.
\]

The candidate should be rejected or revised if any of the following occurs:

1. finer midpoint rungs fail to approach the component fractions with an even-power error profile;
2. an independently written implementation disagrees on any coarse rung;
3. exact polytope integration gives a different component fraction;
4. the arithmetic-side joint-class transport reveals that the model class being evaluated is not the class consumed by the \(m_7\) ledger.

No theorem should consume \(1/8\) until the exact gate passes.
