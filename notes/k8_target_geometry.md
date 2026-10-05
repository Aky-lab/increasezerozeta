# Exact eighth-order target geometry

These are algebraic constraints on supplied moments of a positive measure
on [0,infinity). The model inputs and their counting conversion remain
conditional on the analytic framework.

Use the audited lower moments

    (m0,...,m6) = (1,1,4/3,2,13/4,101/18,640/63).

Their reconstruction is in the [local-identity proof](model_local_identities.md),
with the correction documented in the [pairing audit](../results/pairing_model_audit.md).

## Stieltjes floor

Put G=(m_(i+j+1)) for 0<=i,j<=2 and h=(m4,m5,m6)^T.
The shifted four-by-four moment matrix is positive semidefinite only if

\[
m_7\ge m_7^*:=h^TG^{-1}h=\frac{1031677}{54096}
=19.0712252292\ldots.
\]

The current model ledger m7=3439/180 lies strictly above this floor.

## Eighth-moment floor and Schur complement

Put H=(m_(i+j)) for 0<=i,j<=3, b=(m4,m5,m6,m7)^T and

\[
m_8^{\min}(m_7)=b^TH^{-1}b,\qquad
\delta=m_8-m_8^{\min}(m_7),\qquad g=m_7-m_7^*.
\]

The ordinary five-by-five Hankel matrix is positive semidefinite only
if delta>=0. For delta>0 its block inverse gives

\[
\lambda_4(0)=\frac{\delta}{K\delta+a^2},\quad
K=e_0^TH^{-1}e_0=\frac{2519}{247},\quad
a=e_0^TH^{-1}b.
\]

The linear expression a vanishes at m7=m7*; its slope is
(H^-1)_(0,3). Consequently

\[
\lambda_4(0)=\lambda_3(0)\frac{\delta}{\delta+c g^2},\qquad
\lambda_3(0)=\frac{247}{2519},\qquad c=\frac{3732624}{622193}.
\]

The inverse formula is used in the positive-definite interior. Singular
boundary points require a separate moment analysis, rather than evaluating
an undefined 0/0 expression.

## Conditional counting targets

For a target S put L=(1-S)/2. If L<lambda3, the inequality
1-2*lambda4>=S is equivalent to

\[
\delta\le\frac{Lc}{\lambda_3(0)-L}g^2.
\]

If L>=lambda3, the degree-three certificate already supplies the target
under the same conditional conversion. Thus 80% and 80.2% impose no
additional eighth-moment cap for these inputs.

At m7=3439/180, write m8=3311/90+A8. The moment floor gives

    A8 >= 6775529/26142480 = 0.2591769794... .

The stronger conditional targets require

| Target | Upper bound on A8 |
|---|---:|
| 81% | 8128408/16967475 = 0.4790581981... |
| 90% | 28456753/106766100 = 0.2665336001... |

These are target constraints. The [complete finite certificate](eighth_order_certificate.md)
gives A8=633/2240, inside the 81% interval and above the 90% cap.
Arithmetic transport remains unresolved.

`scripts/k8_target.py` derives all constants from the shared moment
sequence. It checks the gap formula against direct five-by-five rational
matrix inversion, verifies the target equality and checks alternating
coefficient signs for the resulting half-line certificates.
