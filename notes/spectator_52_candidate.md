# Seventh-moment joint class \(\{5,2\}\): candidate history and active target

**Status:** unresolved model-side constant. The original \(1/8\) candidate is **retired**. The active numerical candidate is

\[
\boxed{\{5,2\}=\frac7{72}=0.0972222222\ldots}
\]

and must not be consumed as a theorem-level value until exact rational integration agrees.

## 1. Retired first candidate

The first five-dimensional midpoint discretisation treated the spectator frequency \(v\) numerically together with the four \(C_5\) variables. Coarse rungs drifted towards

\[
\{5,2\}=\frac18.
\]

That candidate was correctly pre-registered before finer work, but it did **not** survive the next structural check.

The exact spectator-frequency elimination in
`notes/spectator_52_eliminate_v.md`
showed that the apparent \(1/8\) convergence was a quadrature artefact caused by crossing rational kink hyperplanes in the \(v\)-direction.

Therefore

\[
\boxed{\{5,2\}=1/8\quad\text{is retired.}}
\]

Any script or note still using \(1/8\) is stale and should be corrected rather than interpreted as an alternative candidate.

## 2. Exact-\(v\) reduction

For fixed five-block frequencies, the \(v\)-integral

\[
\int |v|\,O_7(v)\,dv
\]

is evaluated analytically. This reduces the problem from a five-dimensional midpoint quadrature to a four-dimensional piecewise-polynomial integral.

The exact formula and derivation are in
`notes/spectator_52_eliminate_v.md`.

## 3. Corrected numerical ladder

With \(v\) eliminated exactly, the corrected total
\(\{5,2\}=7(U_1+U_2+U_3)\) is:

| mesh \(dc\) | corrected total |
|---:|---:|
| 0.2500 | 0.1024055481 |
| 0.2000 | 0.1008627200 |
| 0.1250 | 0.0987835526 |
| 0.1000 | 0.0982422875 |
| 0.0800 | 0.0970260964 |
| 0.0625 | 0.0976295331 |
| 0.0500 | 0.0974842118 |

Because the remaining four-dimensional midpoint rule still crosses rational kink hyperplanes, these values are **not enclosures** and need not be monotone.

Two genuinely dyadic Richardson probes are:

\[
0.125\to0.0625:\quad0.09724485998,
\]

\[
0.100\to0.0500:\quad0.09723151986.
\]

Both are close to

\[
\frac7{72}=0.09722222222\ldots.
\]

Thus the active candidate to falsify is

\[
\boxed{\{5,2\}=\frac7{72}}.
\]

## 4. Required decisive gate

The exact-\(v\) reduction leaves a four-dimensional rational piecewise-polynomial integral. The decisive next step is exact rational cell/polytope integration.

The active candidate must be rejected or revised if that exact computation disagrees.

## 5. Downstream consequence if certified

Together with the source identification

\[
C_7=-\frac{17}{360},
\]

the exact Bell(7) ledger would give

\[
m_7
=
\frac{1717}{90}
+\frac7{72}
-\frac{17}{360}
=
\boxed{\frac{3443}{180}}
=
19.1277777\ldots.
\]

This lies only slightly above the exact Stieltjes floor from \(m_0,\ldots,m_6\), so the \(k=8\) target geometry is correspondingly tight.
