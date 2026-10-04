# Arithmetic-to-continuum passage and remainder ledger

**Status:** research proof skeleton for external audit.

This note isolates the limiting argument needed to turn the local class-subtracted universality estimate into the exact continuum constant. It deliberately separates that limit from the independent sawtooth-tail envelope and from the inherited covered-zone/band/glue analysis.

## 1. Normalised modulus measure

Let

\[
\ell=\log\frac{T}{2\pi}
\]

and define the positive measure

\[
\mu_\ell
=
\frac1\ell
\sum_{n\le e^\ell}
\frac{\Lambda(n)}n\,
\delta_{\log n/\ell}.
\]

For every continuous \(f\) on \([0,1]\),

\[
\int f\,d\mu_\ell
\longrightarrow
\int_0^1 f(\beta)\,d\beta.
\tag{1.1}
\]

A standard proof uses partial summation and

\[
\sum_{n\le x}\frac{\Lambda(n)}n
=
\log x+O(1).
\]

No quantitative rate is required for the asymptotic proportion argument.

The higher-prime-power part is negligible:

\[
\frac1\ell
\sum_{\substack{p^k\le e^\ell\\k\ge2}}
\frac{\log p}{p^k}
=
O(\ell^{-1}).
\tag{1.2}
\]

Likewise, the diagonal prime pairs have product-measure mass

\[
\frac1{\ell^2}
\sum_p\frac{(\log p)^2}{p^2}
=
O(\ell^{-2}).
\tag{1.3}
\]

Thus the limiting two-modulus measure is carried by distinct ordinary primes.

## 2. The scaled geometric region

Set

\[
\theta
=
2\pi e^{-(\nu-1)\ell},
\qquad 1\le\nu\le2.
\tag{2.1}
\]

The source-code conditions

\[
\theta\le\frac{T}{\max(b_1,b_2)^2}
\]

and, in the \(N>X\) zone,

\[
\min(b_1,b_2)\ge\frac{2\pi}{\theta},
\]

become, after

\[
\beta_i=\frac{\log b_i}{\ell},
\]

the polygon

\[
D_\nu
=
\left\{
(\beta_1,\beta_2):
\nu-1\le\beta_i\le\frac\nu2
\right\}.
\tag{2.2}
\]

The boundaries are finitely many line segments and hence have Lebesgue measure zero.

The overlap function in the source is homogeneous in the logarithmic variables: if \(Q_\ell\) denotes the finite-height overlap weight, then on the scaled variables

\[
Q_\ell=\ell H(\beta_1,\beta_2,\nu)
\tag{2.3}
\]

with \(H\) bounded, continuous, and piecewise linear on the compact limiting region.

Integrating \(H\) over \(D_\nu\) gives

\[
Q(\nu)
=
\iint_{D_\nu}H(\beta_1,\beta_2,\nu)
\,d\beta_1\,d\beta_2
=
\frac{(2-\nu)^3}{6}
+
\frac43\left(\frac32-\nu\right)_+^3.
\tag{2.4}
\]

## 3. Uniform sawtooth law for ordinary distinct primes

For distinct odd primes \(p,q\), let

\[
\mathcal G_{p,q}(\theta)
=
\sum_{d\ge1}\gamma_d^{p,q}G_0(d\theta),
\]

where the coefficients use the **class-subtracted convention** of the reference implementation.

The proof in notes/class_subtracted_universality.md gives

\[
\sum_d
|\gamma_d^{p,q}-\gamma_d^{1,1}|
\le
12\left(\frac1p+\frac1q\right)
\tag{3.1}
\]

and \(|G_0|\le\pi\). In the contributing region,

\[
\min(p,q)\ge\frac{2\pi}{\theta},
\]

so

\[
|\mathcal G_{p,q}(\theta)-\mathcal G_{1,1}(\theta)|
\le12\theta.
\tag{3.2}
\]

For the universal class,

\[
\mathcal G_{1,1}(\theta)
=
-\theta\log\frac1\theta+O(\theta).
\tag{3.3}
\]

Therefore, uniformly over ordinary distinct prime pairs in the contributing region,

\[
\frac{\mathcal G_{p,q}(\theta)}{\theta\ell}
=
-(\nu-1)+O(\ell^{-1})
\tag{3.4}
\]

whenever \(\nu\ge1+\delta\), with fixed \(\delta>0\).

A global bound of the shape

\[
\left|
\frac{\mathcal G_{p,q}(\theta)}{\theta\ell}
\right|
\ll 1+\nu-1
\tag{3.5}
\]

is enough to remove the \(\nu\in[1,1+\delta]\) strip after first taking \(\ell\to\infty\) and then \(\delta\to0\). One obtains such a bound by combining (3.2) with the partial-summation proof of (3.3); when \(\theta\) stays bounded away from zero, absolute summability of the universal coefficient sequence gives the required \(O(\ell^{-1})\) normalised bound.

## 4. Weak convergence with moving polygonal boundaries

Consider the three-dimensional product measures

\[
d\nu\otimes\mu_\ell\otimes\mu_\ell.
\]

By (1.1),

\[
d\nu\otimes\mu_\ell\otimes\mu_\ell
\Longrightarrow
d\nu\,d\beta_1\,d\beta_2
\tag{4.1}
\]

on the compact box \([1,2]\times[0,1]^2\).

The indicator of

\[
D=
\left\{
(\nu,\beta_1,\beta_2):
1\le\nu\le2,\;
\nu-1\le\beta_i\le\nu/2
\right\}
\tag{4.2}
\]

is discontinuous only on a finite union of planes, a set of limiting measure zero. Since \(H\) is bounded and continuous, the standard continuity-set form of weak convergence gives

\[
\int_D H\,
d\nu\,d\mu_\ell(\beta_1)\,d\mu_\ell(\beta_2)
\longrightarrow
\int_D H\,
d\nu\,d\beta_1\,d\beta_2.
\tag{4.3}
\]

If one wants an elementary epsilon proof rather than invoking the continuity-set theorem, remove a \(\delta\)-neighbourhood of the boundary planes, apply uniform weak convergence on the resulting finite union of compact cells, and then bound the discarded strip by \(O(\delta)\).

## 5. Normalisation audit

The finite model has the normalisation

\[
\frac1{\pi\,d\,\ell_1^4},
\qquad
d=\ell X,
\qquad
\ell_1=\ell+2\log2-1,
\qquad
X\sim\frac{T}{2\pi}.
\tag{5.1}
\]

After the change (2.1),

\[
d\theta=-\ell\theta\,d\nu.
\]

The modulus sums contribute two factors of \(\ell\) from \(\mu_\ell\), the overlap contributes one factor of \(\ell\) by (2.3), and the transformed sawtooth term contributes

\[
\ell\,\frac{\mathcal G_{p,q}(\theta)}{\theta}
=
\ell^2
\frac{\mathcal G_{p,q}(\theta)}{\theta\ell}.
\]

Consequently the normalised core can be written, up to a factor tending to one, as

\[
C_\ell
=
2\left(\frac{\ell}{\ell_1}\right)^4
\int_D
\frac{\mathcal G_{b_1,b_2}(\theta)}{\theta\ell}
H(\beta_1,\beta_2,\nu)
\,d\nu\,d\mu_\ell(\beta_1)\,d\mu_\ell(\beta_2)
+
o(1),
\tag{5.2}
\]

after removing the higher-prime-power and diagonal classes using (1.2)-(1.3).

Using (3.4), (4.3), and then removing the \(\nu=1\) strip,

\[
C_\ell
\longrightarrow
-2
\int_1^2
(\nu-1)Q(\nu)\,d\nu.
\tag{5.3}
\]

Because

\[
\int_1^2(\nu-1)Q(\nu)\,d\nu=\frac1{96},
\]

we obtain

\[
\boxed{
C_{\mathrm{core}}=-\frac1{48}.
}
\tag{5.4}
\]

The factor \(2\) in (5.2) agrees with the reproduction code's 2*pi/(pi*l*ell1**4) normalisation after the scaled modulus sums and is the principal factor-of-two check in this passage.

## 6. Remainder ledger

A theorem-level statement should not hide all unresolved terms under a single symbol. The proposed assembly should be recorded as

\[
R(1)
\le
C_{\mathrm{core}}
+
E_{\mathrm{tail}}
+
E_{\mathrm{fixed}}
+
o(1).
\tag{6.1}
\]

The terms mean:

### \(C_{\mathrm{core}}\)

The arithmetic-to-continuum deterministic core. The proposed exact value is

\[
C_{\mathrm{core}}=-\frac1{48}.
\]

The former numerical continuum-convergence allowance should **not** be retained once (5.4) is proved.

### \(E_{\mathrm{tail}}\)

The fixed-\(P\) sawtooth/frequency-tail envelope from the decomposition

\[
\mathfrak S=\bar D+g_P+g_{\mathrm{tail}}.
\]

This is distinct from the logarithmically growing singular-series truncation used elsewhere in the fourth-moment argument.

The reference candidate quotes a charge \(0.0111\) at its operating point. A new proof must state exactly whether it is using:

- a proven asymptotic bound
  \[
  \limsup_{T\to\infty}E_{\mathrm{tail}}(T)\le0.0111,
  \]
  or
- only finite-height evidence.

Those are not interchangeable.

### \(E_{\mathrm{fixed}}\)

Any other non-vanishing deterministic allowance that survives the limit. At present no such term should be inserted without naming its source.

### \(o(1)\)

Only terms actually proved to vanish may be placed here: covered-zone errors, band/glue deviations, smoothing collars, low-zone terms, or analogous quantities must each have a cited estimate showing decay after the final normalisation.

## 7. Target for a full one-percentage-point improvement

If all non-vanishing remainder terms are combined *after being individually identified* into

\[
E=E_{\mathrm{tail}}+E_{\mathrm{fixed}},
\]

then

\[
R(1)\le-\frac1{48}+E+o(1).
\]

The degree-two consumption gives

\[
P(E)
=
\frac{720E+269}{18(80E+21)}.
\tag{7.1}
\]

To exceed \(0.683008528\), it suffices that

\[
E<0.0410681242\ldots.
\tag{7.2}
\]

This is the correct coarse target for the next analytic stage. It does **not** mean that every unresolved issue may be renamed \(E\); only non-vanishing, explicitly bounded remainder terms belong in (7.2).
