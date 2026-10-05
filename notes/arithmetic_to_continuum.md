# Arithmetic-to-continuum passage and remainder ledger

**Status:** model-core limit proved with a quantitative rate; the surrounding arithmetic moment reduction remains open.

The self-contained [quantitative core proof](arithmetic_core_limit.md) establishes the explicitly defined class-subtracted prime-power functional as \(-1/48+O(\ell^{-1})\), including the universal coefficient asymptotic, uniform exceptional-class bounds and exact overlap integration. This ledger separates that proved model limit from the fixed-P sawtooth-tail envelope and the covered-zone/band/glue analysis. The latter steps are not covered by the model theorem.

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

uniformly on the entire contributing region \(1\le\nu\le2\). The universal estimate extends to the bounded initial \(\theta\)-range by absolute summability, so a fixed \(\delta\)-strip is unnecessary.

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

for the explicitly defined model functional. Removing exceptional classes requires more than (1.2)-(1.3): the [quantitative proof](arithmetic_core_limit.md) supplies a uniform normalized transform bound for coprime prime powers and an integrated \(O(\ell^{-2})\) bound for shared-base pairs. The total exceptional error is \(O(\ell^{-1})\). Any further error introduced by reducing an actual arithmetic moment to this model must be proved separately.

Using the uniform estimate (3.4) and the quantitative product-measure passage proved in [the core theorem](arithmetic_core_limit.md),

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

For the model functional specified exactly in that theorem, the stronger conclusion is \(C_\ell=-1/48+O(\ell^{-1})\). The measure discrepancy is \(O(\ell^{-1})\) by Mertens' estimate; integration by parts in the two modulus variables controls the moving interval endpoints because the overlap is continuous and piecewise affine.

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

The deterministic model core, whose limit is proved for the explicitly defined functional. Its exact limiting value is

\[
C_{\mathrm{core}}=-\frac1{48}.
\]

The model-core convergence error is proved to vanish. It does not require a non-vanishing numerical allowance; errors in an unproved arithmetic-to-model reduction remain separate obligations.

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

For the precise Ramanujan cutoff in the [obstruction theorem](fixed_cutoff_obstruction.md), the finite core is uniformly \(O((1+\log P)/\ell)\). Every subpower cutoff, including any fixed power of \(\log T\), therefore leaves an absolute complementary-tail envelope with lower limit at least \(1/48>0.0111\). Thus the first option cannot hold for that same absolute model tail. A viable arithmetic estimate must instead control a signed tail, use a sufficiently different cutoff strategy, or bound a separately defined actual-minus-model discrepancy. The full core theorem already includes its infinite model tail.

### \(E_{\mathrm{fixed}}\)

Any other non-vanishing deterministic allowance that survives the limit. At present no such term should be inserted without naming its source.

### \(o(1)\)

Only terms actually proved to vanish may be placed here: covered-zone errors, band/glue deviations, smoothing collars, low-zone terms, or analogous quantities must each have a cited estimate showing decay after the final normalisation.

The [prime-pair progression theorem](prime_pair_progression_transport.md)
supplies arbitrary logarithmic mean-square savings for every divisor
subfamily of full-dyadic, one-chain shifts \(h=qk\), in the published MRT
range. Its weighted bound makes the required consumption norm explicit.
It does not supply the short-position-window or multiple-chain estimates,
and an application must still verify the shift range and weight budget.

The [short-window Fourier transport](short_window_fourier_transport.md)
adds a position-averaged replacement of single-prime sums by Lambda_sharp
for Y>=X^(1/3+epsilon), and a longer-window pointwise option at
Y>=j^(5/8+epsilon). Coupled products are controlled by actual weighted
start histograms, rather than an assumed variance factorization. To consume
this input, the dilated lattice scales and concentration budget must be
verified, and the presieved main term must be identified with the proposed
comb. These steps remain part of the analytic reduction.

The [finite presieved pair theorem](presieved_pair_transport.md) evaluates
the divisor model's Fourier coefficients and pair main term. It gives
actual prime-pair variance O_A(X Y^3 log(3X)^(-A)) averaged jointly over
integer short-window starts and full-window shifts, also with divisor
multiplicity. Its consumption norm retains the position and window scales.
It does not replace a common-position four-prime lock by independent
chains, and it does not establish the final zeta normalization.

The [resolved frame audit](resolved_lock_frame.md) separates a convolution
of autocorrelations from the two-coordinate four-prime lock variance.
Its nonnegative counterexamples rule out a general identity based only
on separate power spectra. The valid full-lock theorem loses a factor Y
when restricted to the raw rectangle scale; cubic or equivalent relative
uniformity must be supplied before that gap is charged as a vanishing error.

[Qualitative rectangle transport](averaged_rectangle_transport.md) supplies
a replacement at the raw Y^4 squared scale averaged over starts and two
undilated shifts. It retains the other prime factors' U^3 norms, controls
the admissible residue normalization and identifies a finite four-point
main term. Its unspecified little-o does not absorb arbitrary logarithmic
or divisor consumption weights, and it does not yet cover the reference's
actual dilated lattice family. The subsequent
[rectangle Euler-tail theorem](rectangle_singular_series_tail.md) replaces
the finite local product by the full singular series in the unweighted
nondegenerate joint average. It does not bound the class-subtracted
frequency tail above or justify singular consumption weights.

[The dilated geometry theorem](dilated_rectangle_geometry.md) gives the
exact raw overlap and heterogeneous cube for the physical locks. It proves
that the local cube normalization stays uniformly bounded for prime-power
moduli and is exactly the finite singular-series second moment over a
complete period. Its deterministic progression U^3 criterion has the
correct balanced raw scale. Uniform prime norms on the shorter progression
windows and the consumed weight ranges are still needed; an ordinary
W-tricked bound is not silently promoted to a growing-modulus estimate.

[Weighted dilation transport](weighted_dilated_prime_transport.md)
supplies the prime replacement and full nondegenerate singular-series
main term for fixed coefficients and a sufficiently slowly growing
prime-power family. Its coupled start-weight condition and consumption
norm are explicit. This resolves that restricted progression transfer,
including translated short-box Euler tails. The positive-power modulus
region, prescribed starts, remaining singular weights and the reference's
consumption normalization still need separate estimates.

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
