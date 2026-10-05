# Class-subtracted \(\ell^1\) universality and the exact fourth-moment core

**Research status:** local coefficient proof. This note proves a local arithmetic lemma for the coefficient convention used by the reference reproduction code. The separate [quantitative core proof](arithmetic_core_limit.md) supplies the universal coefficient asymptotic, exceptional-class bounds and exact continuum geometry. The surrounding fourth-moment analytic reduction remains unproved here.

Reference convention, pinned to revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`:

- [tail_bound.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/tail_bound.py)
- [mains_envelope.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/mains_envelope.py)
- [m1_suite.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/m1_suite.py)

After constructing the Euler tensor, `gamma_signed_array` performs the class-\(d\) exact replacement. The proposition below includes that replacement exactly.

## 1. Local vectors

Fix distinct odd primes \(p,q\). For an ordinary odd prime \(r\) in the \((1,1)\) class, write

\[
g_r=(g_r(0),g_r(1))
=
\left(
1-\frac{1}{(r-1)^2},
\frac{1}{(r-1)^2}
\right).
\]

These entries are non-negative and sum to \(1\).

If \(p\) divides exactly one modulus, the code's local \(\beta_p\)-difference gives

\[
s_p=(s_p(0),s_p(1))
=
\left(
\frac{p}{p-1},
-\frac{1}{p-1}
\right).
\]

Put

\[
\alpha=\frac1{p-1},\qquad \beta=\frac1{q-1}.
\]

Then

\[
s_p=(1+\alpha,-\alpha),\qquad
g_p=(1-\alpha^2,\alpha^2),
\]

and similarly at \(q\). Hence

\[
s_p-g_p
=
(\alpha+\alpha^2,-\alpha-\alpha^2),
\]

so

\[
\|s_p-g_p\|_1=2\alpha(1+\alpha),
\qquad
\|s_p\|_1=1+2\alpha,
\qquad
\|g_p\|_1=1.
\]

At \(2\), because \(p,q\) are odd, both modulus pairs have the same local vector and it forces \(v_2(d)=1\). Its coefficient mass is one; the common forced divisor \(2\) does not change the \(\ell^1\) comparison.

## 2. The full Euler tensor

Let \(R_{p,q}\) be the tensor product of all ordinary local vectors \(g_r\) with \(r\ne2,p,q\). Since each \(g_r\) is non-negative with total mass \(1\),

\[
\|R_{p,q}\|_1=1.
\]

Before the class-\(d\) replacement, the two coefficient measures are therefore

\[
\Gamma^{p,q}
=
\delta_2\otimes s_p\otimes s_q\otimes R_{p,q},
\]

and

\[
\Gamma^{1,1}
=
\delta_2\otimes g_p\otimes g_q\otimes R_{p,q}.
\]

Using the two telescoping decompositions

\[
s_p\otimes s_q-g_p\otimes g_q
=
(s_p-g_p)\otimes s_q
+
g_p\otimes(s_q-g_q)
\]

and the same identity with \(p,q\) interchanged, then averaging the resulting upper bounds, gives

\[
\|\Gamma^{p,q}-\Gamma^{1,1}\|_1
\le
2(1+\alpha)(1+\beta)(\alpha+\beta).
\tag{2.1}
\]

This argument can be made first for finite Euler products; the infinite statement follows by absolute convergence.

## 3. The class-\(d\) exact replacement

This is the step that must be included to match the repository's actual convention.

For \(b_1=p,b_2=q\), the special-prime set is

\[
S=\{2,p,q\}.
\]

The function `_gamma_free_exact` replaces the full coefficient whenever all prime divisors of \(d\) lie in \(S\). Because the local \(2\)-vector forces exponent \(1\), and the \(p,q\) local vectors vanish for exponents \(\ge2\), the only non-zero affected coefficients occur at

\[
d\in\{2,2p,2q,2pq\}.
\]

The replacement subtracts exactly the finite signed measure

\[
C^{p,q}
=
\delta_2\otimes s_p\otimes s_q.
\]

For \((1,1)\), the same convention subtracts

\[
C^{1,1}=\delta_2.
\]

This recovers, in particular, the archived \(d=2\) shift from \(C_2\) to \(C_2-1\).

Therefore the *actual* class-subtracted coefficient sequences satisfy

\[
\gamma^{p,q}
=
\Gamma^{p,q}-C^{p,q},
\qquad
\gamma^{1,1}
=
\Gamma^{1,1}-C^{1,1}.
\tag{3.1}
\]

The class-measure difference can be evaluated exactly. The four coefficients of
\(s_p\otimes s_q-\delta_{(0,0)}\) have absolute sum

\[
\|C^{p,q}-C^{1,1}\|_1
=
2\alpha+2\beta+4\alpha\beta.
\tag{3.2}
\]

Combining (2.1), (3.1), and (3.2) gives

\[
\sum_{d\ge1}
\left|
\gamma_d^{p,q}-\gamma_d^{1,1}
\right|
\le
2(1+\alpha)(1+\beta)(\alpha+\beta)
+
2\alpha+2\beta+4\alpha\beta.
\tag{3.3}
\]

Since \(0<\alpha,\beta\le1/2\),

\[
2(1+\alpha)(1+\beta)(\alpha+\beta)
\le \frac92(\alpha+\beta),
\]

while \(4\alpha\beta\le\alpha+\beta\), so

\[
\boxed{
\sum_{d\ge1}
\left|
\gamma_d^{p,q}-\gamma_d^{1,1}
\right|
\le
8\left(\frac1{p-1}+\frac1{q-1}\right).
}
\tag{3.4}
\]

For odd primes \(p,q\ge3\),

\[
\frac1{p-1}\le\frac{3}{2p},
\]

and hence the convenient coarser form

\[
\boxed{
\sum_{d\ge1}
\left|
\gamma_d^{p,q}-\gamma_d^{1,1}
\right|
\le
12\left(\frac1p+\frac1q\right).
}
\tag{3.5}
\]

This is the class-subtracted \(\ell^1\) universality estimate needed in the continuum argument.

## 4. Sawtooth-transform corollary

The reproduction code uses

\[
G_0(x)
=
2\sum_{h\ge1}\frac{\sin(2xh)-\sin(xh)}{h},
\]

with the sawtooth closed form and the bound

\[
|G_0(x)|\le\pi.
\]

Define

\[
\mathcal G_{p,q}(\theta)
=
\sum_{d\ge1}\gamma_d^{p,q}G_0(d\theta).
\]

Then (3.5) gives immediately

\[
|\mathcal G_{p,q}(\theta)-\mathcal G_{1,1}(\theta)|
\le
12\pi\left(\frac1p+\frac1q\right).
\tag{4.1}
\]

In the relevant \(N>X\) zone, the source code imposes

\[
\min(p,q)\ge \frac{2\pi}{\theta}.
\]

Therefore

\[
\frac1p+\frac1q
\le\frac{\theta}{\pi},
\]

and hence

\[
\boxed{
|\mathcal G_{p,q}(\theta)-\mathcal G_{1,1}(\theta)|
\le 12\theta.
}
\tag{4.2}
\]

Thus distinct prime modulus pairs have the same leading sawtooth law, uniformly in the contributing zone.

## 5. Universal \((1,1)\) sawtooth slope

The [Euler-convolution proof](arithmetic_core_limit.md) establishes, in the pinned convention,

\[
A(M)
:=
\sum_{d\le M}d\gamma_d^{1,1}
=
\log M+c_0+o(1),
\qquad c_0=\gamma-2.
\tag{5.1}
\]

For \(0<x<\pi\), the exact sawtooth formula gives

\[
G_0(x)=-x.
\tag{5.2}
\]

Choose a fixed \(c\in(0,\pi)\) and

\[
M=\left\lfloor\frac c\theta\right\rfloor.
\]

Then (5.2) holds for every \(d\le M\), and

\[
\sum_{d\le M}\gamma_d^{1,1}G_0(d\theta)
=
-\theta A(M)
=
-\theta\log\frac1\theta+O(\theta).
\tag{5.3}
\]

For \(d>2\), the \((1,1)\) class-subtracted coefficients are non-negative. Partial summation applied to (5.1) yields

\[
\sum_{d>M}\gamma_d^{1,1}
=
\frac1M+o(M^{-1}),
\tag{5.4}
\]

so the tail contributes \(O(\theta)\) using \(|G_0|\le\pi\). Therefore

\[
\boxed{
\mathcal G_{1,1}(\theta)
=
-\theta\log\frac1\theta+O(\theta).
}
\tag{5.5}
\]

Together with (4.2),

\[
\boxed{
\mathcal G_{p,q}(\theta)
=
-\theta\log\frac1\theta+O(\theta)
}
\tag{5.6}
\]

uniformly for distinct prime moduli in the contributing zone.

## 6. Exceptional modulus classes

The remaining modulus configurations require uniform transform control in addition to small measure mass. The [quantitative core proof](arithmetic_core_limit.md) provides that control and separates coprime pairs from shared-base pairs.

Higher prime powers have total one-dimensional normalised Mertens mass

\[
\frac1\ell
\sum_{k\ge2}\sum_{p^k\le e^\ell}
\frac{\log p}{p^k}
=
O(\ell^{-1}),
\]

For coprime prime-power pairs, the normalized sawtooth transform is uniformly bounded on the contributing region. Thus configurations containing a genuine higher prime power contribute \(O(\ell^{-1})\). The pair \(b=2\) is handled by the same bound and its \(O(\ell^{-1})\) mass.

The diagonal prime class satisfies

\[
\frac1{\ell^2}\sum_p\frac{(\log p)^2}{p^2}
=
O(\ell^{-2}),
\]

This mass estimate alone is insufficient for the transform. For all shared-base pairs \((p^a,p^b)\), the coefficient \(\ell^1\) bound and the region's upper limit on \(\nu\) give a total contribution

\[
O\!\left(\ell^{-4}\sum_{p^a,p^b\le e^\ell}
\frac{(\log p)^2}{p^{\max(a,b)}}\right)=O(\ell^{-2}).
\]

This proves negligibility of ordinary diagonal primes and shared-base higher powers without assuming a uniform normalized transform for those pairs.

Thus ordinary distinct primes carry the limiting modulus mass.

For bounded piecewise-continuous \(f\), standard Mertens/PNT partial summation gives

\[
\frac1\ell
\sum_n\frac{\Lambda(n)}n
f\!\left(\frac{\log n}{\ell}\right)
\longrightarrow
\int f(\beta)\,d\beta.
\tag{6.1}
\]

For the moving polygonal integration region, discard a strip of width \(\delta\) around its finitely many linear boundaries, apply (6.1) uniformly on the remaining compact cells, then let \(\ell\to\infty\) followed by \(\delta\to0\). The boundary strip has \(O(\delta)\) Lebesgue measure and the overlap geometry is bounded and piecewise linear.

This avoids the need for a sharp global \(O(1/\ell)\) "universality-collapse rate"; convergence is enough for an asymptotic proportion theorem.

## 7. Exact continuum geometry

The cube/simplex derivation in the [quantitative core proof](arithmetic_core_limit.md) gives the integrated overlap geometry

\[
\boxed{
Q(\nu)
=
\frac{(2-\nu)^3}{6}
+
\frac43\left(\frac32-\nu\right)_+^3,
\qquad 1\le\nu\le2.
}
\tag{7.1}
\]

The exact integrals are

\[
\int_1^2Q(\nu)\,d\nu=\frac1{16}
\]

and

\[
\int_1^2(\nu-1)Q(\nu)\,d\nu=\frac1{96}.
\tag{7.2}
\]

Because only the logarithmic slope in (5.5) survives the final normalisation, the continuum core is therefore

\[
\boxed{
C_{\mathrm{core}}
=
-2\int_1^2(\nu-1)Q(\nu)\,d\nu
=
-\frac1{48}.
}
\tag{7.3}
\]

Numerically,

\[
-\frac1{48}=-0.020833333333\ldots,
\]

consistent with the independent numerical continuum estimate near \(-0.0209\) in the reference candidate.

## 8. Consequence inside the one-sided fourth-moment framework

This section is a **conditional consumption statement**: it assumes the surrounding one-sided fourth-moment analytic chain from the reference candidate.

For a precisely defined remainder distinct from the obstructed absolute fixed-cutoff model tail, if a separate asymptotic estimate is established as

\[
E_{\mathrm{tail}}\le0.0111,
\]

while the continuum-convergence slack is removed using (7.3), then

\[
R(1)
\le
-\frac1{48}+0.0111+o(1)
=
-\frac{73}{7500}+o(1).
\]

The degree-two Christoffel conversion then gives

\[
\liminf\frac{N_0^s(T)}{N(T)}
\ge
\frac{1082}{1539}
=
0.7030539311\ldots
\]

within that framework.

The [fixed-cutoff theorem](fixed_cutoff_obstruction.md) proves that the same \(0.0111\) cannot be an asymptotic absolute envelope for its specified complementary model tail: that envelope has lower limit at least \(1/48\). Consequently the displayed \(70.3054\%\) remains formal conditional arithmetic, rather than a bound supported by the fixed-P model-tail method.

More generally, if all *non-\(o(1)\)* remainder terms not included in the exact core are collected into a clearly defined quantity \(E\),

\[
R(1)\le-\frac1{48}+E+o(1),
\]

then the same conversion is

\[
P(E)
=
\frac{720E+269}{18(80E+21)}.
\]

To beat \(0.683008528\), it is enough that

\[
E<0.0410681242\ldots.
\]

The notation \(E\) must not silently absorb unrelated analytic gaps: a final theorem should provide an explicit remainder ledger.

## 9. What remains

The local comparison and the explicitly defined model-core limit are proved in the notes. The remaining work concerns the surrounding analytic chain:

1. prove the reduction of the actual arithmetic moment to the specified class-subtracted model functional;
2. audit the normalization of that analytic reduction against (7.3);
3. replace the obstructed fixed-P absolute-tail estimate with a signed or actual-minus-model bound, with a named remainder ledger;
4. prove or independently review the covered-zone/band/glue inputs and their decay after normalization;
5. independently review the new coefficient, geometry and quantitative model-core proofs.

Until those steps are externally checked, the \(70.3054\%\) figure should be described as a consequence **within the candidate framework**, not as an established unconditional record.
