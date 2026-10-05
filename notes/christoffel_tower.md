# Exact Christoffel consumption for the moment tower

The mass bound below is exact algebra for a supplied positive moment
sequence. Its application to zeta zeros additionally requires the
spectral/counting interface and analytic moment identities.

## 1. General theorem

Let \(\nu\) be a positive probability measure on \([0,\infty)\) with finite moments

\[
m_k=\int x^k\,d\nu(x),\qquad 0\le k\le2n,
\]

and suppose its Hankel moment matrix

\[
H_n=(m_{i+j})_{0\le i,j\le n}
\]

is positive definite.

Let \(e_0=(1,0,\dots,0)^T\), and define

\[
\lambda_n(0)
=
\frac{1}{e_0^T H_n^{-1}e_0}.
\]

Then

\[
\boxed{\nu(\{0\})\le\lambda_n(0).}
\]

Moreover the bound is the optimal one obtainable from the moments \(m_0,\dots,m_{2n}\).

### Proof

Among polynomials \(q\) of degree at most \(n\) satisfying \(q(0)=1\), minimise

\[
\int q(x)^2\,d\nu(x).
\]

Writing \(q(x)=c_0+c_1x+\cdots+c_nx^n\), the objective is \(c^TH_nc\), and the constraint is \(e_0^Tc=1\). Lagrange multipliers give the unique minimiser

\[
c_*=
\frac{H_n^{-1}e_0}{e_0^TH_n^{-1}e_0},
\]

and the minimum is

\[
\int q_*(x)^2\,d\nu(x)
=
\lambda_n(0).
\]

Since \(q_*(0)=1\),

\[
\nu(\{0\})
\le
\int q_*(x)^2\,d\nu(x)
=
\lambda_n(0).
\]

Optimality is the standard truncated moment/Christoffel extremal statement.

There is also a useful structural fact. The first-variation condition against every polynomial \(r\) with \(r(0)=0\) gives

\[
\int q_*(x)\,x\,s(x)\,d\nu(x)=0
\]

for every \(s\) of degree at most \(n-1\). Hence \(q_*\) is, after normalisation at zero, the degree-\(n\) orthogonal polynomial for the positive measure \(x\,d\nu(x)\). Its zeros therefore lie in \((0,\infty)\) when the measure is nondegenerate. Thus

\[
q_*(x)=\prod_{j=1}^n\left(1-\frac{x}{r_j}\right),\qquad r_j>0,
\]

and in particular

\[
q_*(x)\ge1\qquad(x\le0).
\]

So the same square \(q_*^2\) also controls mass on the non-positive half-line if that is the form needed by a particular spectral-counting interface.

## 2. Recovery of the fourth-moment certificate

Using

\[
m_0=1,\quad
m_1=1,\quad
m_2=\frac43,\quad
m_3=2,\quad
m_4=\frac{13}{4},
\]

the degree-two Hankel calculation gives

\[
\lambda_2(0)=\frac5{36}.
\]

The minimiser is

\[
q_2(x)
=
1-\frac{21}{12}x+\frac{8}{12}x^2
=
\frac{8x^2-21x+12}{12}.
\]

Therefore

\[
q_2(x)^2
=
\frac{\left(x^2-\frac{21}{8}x+\frac32\right)^2}{(3/2)^2},
\]

exactly the familiar \(13/18\) certificate.

Thus that certificate is already a Christoffel polynomial in disguise.

## 3. Audited six-moment inputs

The current conditional model sequence is

    (m0,...,m6) = (1,1,4/3,2,13/4,101/18,640/63).

The sixth-order value is assembled from certified class values and the
[finite-model local identities](model_local_identities.md). The
[pairing audit](../results/pairing_model_audit.md) explains its difference
from the separately specified reference moment. Arithmetic transport and the
spectral/counting interface remain open.

For this sequence,

\[
H_3=\begin{pmatrix}
1&1&4/3&2\\
1&4/3&2&13/4\\
4/3&2&13/4&101/18\\
2&13/4&101/18&640/63
\end{pmatrix},\qquad
\det H_3=\frac{247}{108864}>0.
\]

Exact inversion gives

\[
\lambda_3(0)=\frac{247}{2519},\qquad
q_3(x)=1-\frac{8232}{2519}x+\frac{7368}{2519}x^2
-\frac{1932}{2519}x^3.
\]

Every coefficient of q_3(-t) is nonnegative, with constant term one.
Thus q_3(x)^2>=1 for x<=0, and the same certificate bounds mass on
the nonpositive half-line.

If the analytic framework supplies these moments and the counting
conversion, the resulting simple-zero fraction is

    1-2*lambda_3(0) = 2025/2519 = 0.8038904327...,

and the distinct-zero fraction is 2272/2519. These are conditional
conversions, not established bounds for zeta zeros.

## 4. Higher moments

For moments through m_(2n), form H_n, solve H_n c=e0 exactly and
normalize c at zero. The reciprocal e0^T H_n^-1 e0 gives the
Christoffel value. The remaining research concerns producing valid
moments and proving their analytic interface.

[Eighth-order target geometry](k8_target_geometry.md) derives the
degree-four extension with the audited inputs. Run
`python scripts/christoffel_exact.py` to check the degree-two and
degree-three calculations with rational arithmetic.
