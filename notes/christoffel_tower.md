# Exact Christoffel consumption for the moment tower

**Status:** rigorous algebraic reduction, conditional only on the same spectral-measure/counting interface used by the reference higher-moment programme.

The current six-moment candidate uses a rationalised cubic whose roots are approximately
\(0.5323,1.3122,2.0586\). Those numbers are not mysterious: they are approximations to the exact roots of the degree-three Christoffel minimiser at the origin for the exact moment sequence.

The point of this note is broader than the small numerical improvement. Once exact moments through \(m_{2n}\) are known, the optimal origin-mass certificate is obtained automatically from the Hankel matrix. No bespoke LP atom search is needed.

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

## 3. Exact six-moment certificate

Now use the exact moments

\[
m_0=1,\qquad
m_1=1,\qquad
m_2=\frac43,\qquad
m_3=2,\qquad
m_4=\frac{13}{4},
\]

\[
m_5=\frac{101}{18},\qquad
m_6=\frac{12809}{1260}.
\]

The \(4\times4\) Hankel matrix is

\[
H_3=
\begin{pmatrix}
1&1&4/3&2\\
1&4/3&2&13/4\\
4/3&2&13/4&101/18\\
2&13/4&101/18&12809/1260
\end{pmatrix}.
\]

Its determinant is

\[
\det H_3=\frac{283}{108864}>0.
\]

A direct exact inversion gives

\[
e_0^TH_3^{-1}e_0=\frac{13891}{1415},
\]

hence

\[
\boxed{
\lambda_3(0)=\frac{1415}{13891}.
}
\]

The exact minimising cubic is

\[
\boxed{
q_3(x)
=
1-\frac{43428}{13891}x
+\frac{37704}{13891}x^2
-\frac{9660}{13891}x^3.
}
\]

Equivalently,

\[
q_3(x)
=
\frac{13891-43428x+37704x^2-9660x^3}{13891}.
\]

Its three roots are approximately

\[
0.5323430232,\qquad
1.3122059467,\qquad
2.0585566202.
\]

These are precisely the values that the reference candidate rationalised to approximately

\[
0.5323,\qquad1.3122,\qquad2.0586.
\]

Because the coefficients alternate,

\[
q_3(-t)
=
1+\frac{43428}{13891}t
+\frac{37704}{13891}t^2
+\frac{9660}{13891}t^3
\ge1
\qquad(t\ge0),
\]

so \(q_3^2\ge1\) on \((-\infty,0]\) as well.

The exact origin-mass bound is therefore

\[
w_0\le\frac{1415}{13891}.
\]

Under the same simple-zero counting conversion used in the reference moment programme,

\[
\boxed{
\frac{N_0^s}{N}
\ge
1-2\frac{1415}{13891}
=
\frac{11061}{13891}
=
0.7962709668\ldots
}
\]

and

\[
\boxed{
\frac{N_d}{N}
\ge
1-\frac{1415}{13891}
=
\frac{12476}{13891}
=
0.8981354834\ldots.
}
\]

This is slightly stronger than the rounded \(0.7962/0.8981\) headline, but the main gain is conceptual: the certificate is exact and canonical.

## 4. Why this matters for higher moments

Suppose moments through \(m_{2n}\) are pinned exactly.

The consumption step becomes:

1. form \(H_n=(m_{i+j})\);
2. solve \(H_nc=e_0\) exactly;
3. normalise \(c\) so that \(q(0)=1\);
4. compute
   \[
   \lambda_n(0)=1/(e_0^TH_n^{-1}e_0);
   \]
5. feed \(1-2\lambda_n(0)\) into the same spectral counting interface.

No atom fitting and no numerical LP are needed.

This is especially relevant to the proposed \(k=8\) rung. Once exact \(m_7,m_8\) are supplied, the degree-four certificate is automatic. The open difficulty is therefore entirely on the moment-production side, not the consumption side.

## 5. Research direction

This suggests separating the programme into two engines.

### Moment engine

Prove exact or one-sided bounds for

\[
m_7,m_8,m_9,\dots
\]

from the arithmetic/compressed-matrix side.

### Consumption engine

Use exact Hankel/Christoffel algebra to convert any pinned even moment tower into the optimal moment-only mass-at-zero bound.

The second engine is now essentially solved for exact moments.

A longer-term structural question is whether the resulting Christoffel sequence

\[
\lambda_n(0)
\]

for the sine-model moment sequence can be analysed asymptotically without calculating every moment individually. If so, that could expose the rate at which the moment tower approaches \(100\%\).
