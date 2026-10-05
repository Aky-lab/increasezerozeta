# Seventh-moment joint class \(\{5,2\}\): independent specification

**Status:** model-side specification and combinatorial gate. The arithmetic transport is not yet proved.

The \(k=8\) route first needs a rigorous \(m_7\) input. Among the genuinely new seventh-moment objects, the key mixed class is \(\{5,2\}\).

This note defines that class independently of the reference implementation.

## 1. Placement combinatorics

A seven-cycle with one 5-block and one 2-block has

\[
\binom72=21
\]

placements of the 2-block.

Up to cyclic distance, there are three classes:

\[
d=1,\qquad d=2,\qquad d=3,
\]

each with multiplicity \(7\). Hence

\[
\boxed{
\{5,2\}=7U_1+7U_2+7U_3.
}
\]

The exact enumeration is checked by scripts/bell7_gate.py.

## 2. Frequency variables

Let the 2-block carry frequencies

\[
(v,-v),
\]

and let the 5-block carry

\[
(c_1,c_2,c_3,c_4,c_5),
\qquad
c_5=-c_1-c_2-c_3-c_4.
\]

The connected 5-block factor is the sine-model five-point cumulant

\[
C_5(c_1,\dots,c_5),
\]

defined by the partition-cyclic formula.

The pair factor is

\[
C_2(v)=\min(|v|,1).
\]

For every placement class, the integrand is

\[
O_7(\text{full prefix walk})\,
C_5(c_1,\dots,c_5)\,
C_2(v),
\]

where

\[
O_7(p_0,\dots,p_6)
=
\bigl(1-(\max p_j-\min p_j)\bigr)_+.
\]

## 3. Prefix walks

Use the cyclic order in which the first pair frequency is \(v\).

### Distance \(d=1\)

Increment sequence:

\[
(v,-v,c_1,c_2,c_3,c_4,c_5).
\]

Prefix walk:

\[
0,\;
v,\;
0,\;
c_1,\;
c_1+c_2,\;
c_1+c_2+c_3,\;
c_1+c_2+c_3+c_4.
\]

### Distance \(d=2\)

Increment sequence:

\[
(v,c_1,-v,c_2,c_3,c_4,c_5).
\]

Prefix walk:

\[
0,\;
v,\;
v+c_1,\;
c_1,\;
c_1+c_2,\;
c_1+c_2+c_3,\;
c_1+c_2+c_3+c_4.
\]

### Distance \(d=3\)

Increment sequence:

\[
(v,c_1,c_2,-v,c_3,c_4,c_5).
\]

Prefix walk:

\[
0,\;
v,\;
v+c_1,\;
v+c_1+c_2,\;
c_1+c_2,\;
c_1+c_2+c_3,\;
c_1+c_2+c_3+c_4.
\]

These are the only three cyclic-distance types because \(7\) is odd.

## 4. Integral definition

For \(d=1,2,3\), define

\[
U_d
=
\int_{\mathbb R^5}
O_7(W_d(v,c_1,c_2,c_3,c_4))
\,C_2(v)\,
C_5(c_1,c_2,c_3,c_4,-c_1-c_2-c_3-c_4)
\,dv\,dc_1\cdots dc_4.
\]

The overlap support makes the integral compact.

Then

\[
\boxed{
\{5,2\}=7(U_1+U_2+U_3).
}
\]

## 5. Exact-rationality expectation

Each partition-cyclic term of \(C_5\) is a signed overlap volume. Multiplying by the outer seven-walk overlap and \(C_2(v)\) gives a piecewise polynomial with rational coefficients on a finite arrangement of rational hyperplanes.

Therefore the same polytope-volume argument used for lower connected constants predicts

\[
\boxed{\{5,2\}\in\mathbb Q.}
\]

This should be proved by an exact polytope integrator, not merely inferred from high-precision numerics.

## 6. Required gates

Any numerical or symbolic implementation should pass the following before its output is consumed:

1. **Bell gate:** \(21\) placements split \(7/7/7\).
2. **C5 gate:** the internal five-point cumulant evaluator reproduces the already certified \(C_5=1/36\) pure-cycle constant under the appropriate pure-cycle integration.
3. **Negation symmetry:** the full integrand is invariant under
   \[
   (v,c_1,\dots,c_4)\mapsto-(v,c_1,\dots,c_4).
   \]
4. **Support gate:** direct full-grid evaluation and support-pruned evaluation agree.
5. **Refinement gate:** midpoint errors exhibit the expected even-power convergence before any rational reconstruction is attempted.
6. **Exact gate:** the final rational candidate is independently recovered by exact polytope integration.

## 7. Arithmetic task after model evaluation

A model value for \(\{5,2\}\) is not enough.

The corresponding prime-side joint class must be transported through the same truncated arithmetic/spectral framework as the lower moments. Only after that transport and the \(C_7\) term are certified can the resulting \(m_7\) be fed into the exact \(k=8\) target geometry.

The immediate goal is therefore:

\[
\boxed{
\text{model-side exact }\{5,2\}
\quad+\quad
\text{arithmetic transport}
\quad\Longrightarrow\quad
\text{rigorous }m_7>m_7^*.
}
\]
