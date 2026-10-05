# Exact elimination of the spectator frequency in \(\{5,2\}\)

**Status:** exact elementary reduction.  This removes one integration
dimension from the model-side \(\{5,2\}\) problem.

Let the five-block prefix set be

\[
B=\{0,s_1,s_2,s_3,s_4\},
\qquad
s_j=c_1+\cdots+c_j,
\]

and write

\[
m=\min B,\qquad M=\max B,\qquad R=M-m.
\]

For pair distance \(d=1,2,3\), the pair-shifted prefix subset is

\[
A_1=\{0\},\qquad
A_2=\{0,s_1\},\qquad
A_3=\{0,s_1,s_2\}.
\]

Write

\[
a=\min A_d,\qquad A=\max A_d,
\]

and define the margins

\[
\alpha=a-m\ge0,\qquad
\beta=M-A\ge0.
\]

Because \(A_d\subseteq B\),

\[
A-a\le M-m=R.
\]

If \(R\ge1\), the outer seven-cycle overlap is identically zero for
all spectator frequencies \(v\).  Assume \(R<1\).

The full outer prefix set is

\[
B\cup(v+A_d).
\]

Its overlap is

\[
O_7(v)
=
\left[
1-
\left(
\max(M,v+A)-\min(m,v+a)
\right)
\right]_+.
\]

Since \(A_d\subseteq B\), the four breakpoints are ordered as

\[
\ell=M-a-1
<
b=m-a\le0\le
c=M-A
<
u=1+m-A.
\]

Thus

\[
O_7(v)=
\begin{cases}
v-\ell,&\ell\le v\le b,\\
1-R,&b\le v\le c,\\
u-v,&c\le v\le u,\\
0,&\text{otherwise}.
\end{cases}
\]

The pair cumulant is \(C_2(v)=|v|\) on this support.  Therefore the
spectator frequency can be integrated exactly:

\[
J_d(c_1,\ldots,c_4)
:=
\int_{\mathbb R}|v|O_7(v)\,dv.
\]

Splitting at \(b,0,c\) and integrating elementary quadratics gives

\[
\boxed{
J_d
=
\frac{1-R}{6}
\left[
2R^2
-3R(\alpha+\beta)
-4R
+3\alpha^2+3\alpha
+3\beta^2+3\beta
+2
\right].
}
\]

Hence

\[
\boxed{
U_d
=
\int_{\mathbb R^4}
J_d(c_1,\ldots,c_4)
\,C_5(c_1,\ldots,c_5)
\,dc_1\cdots dc_4,
}
\]

with \(c_5=-c_1-c_2-c_3-c_4\), and

\[
\boxed{\{5,2\}=7(U_1+U_2+U_3).}
\]

## Why this matters computationally

The original definition is a five-dimensional weighted integral.
After this reduction:

- the spectator variable is eliminated exactly;
- the remaining domain is four-dimensional;
- on every order chamber of the five-block prefix walk,
  \(J_d\) is a cubic polynomial;
- each partition-cyclic term of \(C_5\) contributes one further
  linear overlap factor.

Thus every term is a degree-at-most-four polynomial on a rational
polyhedral cell.  The kink hyperplanes are the same small-integer
prefix-difference family already used in the exact lower-moment
integrators.

This makes exact rational integration of \(\{5,2\}\) substantially
closer in complexity to the already-certified four-dimensional
\(C_5\) computation than to a new five-dimensional weighted problem.
