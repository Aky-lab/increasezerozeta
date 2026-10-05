# Weighted short-window Fourier transport to a presieved prime model

Published short-interval uniformity gives an averaged replacement for
single-prime Fourier sums, with a precise position-weight budget. It also
controls products of such sums with arbitrarily coupled parameters. No
variance factorization or independence of window positions is needed.

The replacement is by a presieved function Lambda_sharp, not directly by
the continuum Gram model or a small Ramanujan-denominator truncation.
Identifying and consuming that main term remains a separate problem.

## 1. Published input and its range

Let X,Y be positive integers, X>=3, and fix epsilon>0 with

    X^(1/3+epsilon)<=Y<=X^(1-epsilon).

Put L=log(3X), R=exp((log X)^(1/10)), P(R)=product_(p<R) p, and

    Lambda_sharp,X(n)=[P(R)/phi(P(R))] 1_(gcd(n,P(R))=1).

[Matomaki, Radziwill, Shao, Tao and Teravainen, Higher uniformity of
arithmetic functions in short intervals II, arXiv:2411.05770v2,
Theorem 1.1(ii)](https://arxiv.org/pdf/2411.05770v2)
implies the following specialization. Outside a set of real x in [X,2X]
of measure O_(B,epsilon)(X L^(-B)), for every real alpha and every
arithmetic progression Q contained in (x,x+Y] intersect Z,

    |sum_(n in Q) [Lambda(n)-Lambda_sharp,X(n)] e(alpha n)|
      <= O_(B,epsilon)(Y L^(-B)).                                   (1)

Here e(t)=exp(2 pi i t). Choose the fixed one-dimensional torus, the
fixed Lipschitz test F(t)=e(t), and polynomial sequence g(n)=alpha*n.
Fix a small constant delta accommodating its complexity and norm.
The theorem's supremum over g and its maximal progression norm give
one exceptional set valid for all alpha and Q. No union over frequency
bins or moduli is required. Changing log X to L changes constants only.

The deep short-interval estimate is the published theorem. The discrete,
weighted and coupled-product deductions below are applications developed
here, without a claim of priority.

## 2. From real exceptional measure to integer-position L^2

For integer j in J={X,...,2X-1}, define

    D(j)=sup_(alpha,Q subset (j,j+Y] intersect Z, Q an AP)
          |sum_(n in Q) [Lambda(n)-Lambda_sharp,X(n)] e(alpha n)|.

Define the corresponding D(x) at real x. Because Y is an integer, the
set of integers in (x,x+Y] is {j+1,...,j+Y} for every x in [j,j+1).
The sequence Lambda_sharp,X is fixed across these windows. Consequently
D(x)=D(j) on that entire unit cell, and

    sum_(j in J) D(j)^2=integral_X^(2X) D(x)^2 dx.                    (2)

This cell identity is essential: a small exceptional set of real measure
does not by itself control a prescribed set of integer positions. Here
the discrepancy is constant on full unit cells, so it does.

Mertens' estimate gives P(R)/phi(P(R))=O(log R)=O(L). Also
Lambda(n)<=log(3X)=L for n<=2X+Y<=3X. Therefore D(x)=O(Y L)
pointwise. Split (2) into the good and exceptional sets in (1):

    sum_j D(j)^2
      =O_(B,epsilon)(X Y^2 [L^(-2B)+L^(2-B)]).

Since B is arbitrary, for every fixed A>0,

    sum_(j in J) D(j)^2=O_(A,epsilon)(X Y^2 L^(-A)).                 (3)

The estimates are asymptotic, with inherited theorem constants; no
effective finite-height threshold or numerical constant is asserted.
For nonintegral Y, the cell argument needs adjustment. Rounding a length
changes a single-prime sum by at most O(L) per endpoint, but a consumed
normalization must account for that error. Statements here use integer Y.

## 3. Bounded-variation weights and the one-window budget

For each j choose any progression Q_j in its window, any alpha_j, and
a complex weight f_j on the integer window. Define

    V_j=max_(j<n<=j+Y)|f_j(n)|
         +sum_(n=j+1)^(j+Y-1)|f_j(n+1)-f_j(n)|.

Its variation on an ordered subprogression is at most this full variation.
Discrete partial summation, using initial subprogressions of Q_j, gives

    |sum_(n in Q_j) f_j(n)[Lambda(n)-Lambda_sharp,X(n)]e(alpha_j n)|
       <=V_j D(j).                                                   (4)

All choices may depend on j. For arbitrary complex coefficients a_j,
Cauchy--Schwarz and (3), using twice the desired saving, imply

    |sum_j a_j sum_(n in Q_j) f_j(n)
                         [Lambda(n)-Lambda_sharp,X(n)]e(alpha_j n)|
      =O_(A,epsilon)(sqrt(X) Y L^(-A)
                     [sum_j |a_j|^2 V_j^2]^(1/2)).                  (5)

This is a proved short-position-window replacement when the position
weights meet the displayed norm budget. It retains the actual window
length Y, rather than silently substituting Y for the position scale X
in a full-scale endpoint theorem.

## 4. Coupled products: a marginal concentration criterion

Fix r. Let Omega be any finite parameter family. For each omega and slot
i=1,...,r choose a start j_i(omega) in J, a progression Q_i(omega), a
frequency alpha_i(omega), and a weight f_i(omega,n) with variation size
V_i(omega) defined as above. Use the same common X,Y and model
Lambda_sharp,X. Set

    S_i(omega)=sum_(n in Q_i(omega)) Lambda(n) f_i(omega,n)e(alpha_i n),
    T_i(omega)=sum_(n in Q_i(omega)) Lambda_sharp,X(n) f_i(omega,n)e(alpha_i n).

For arbitrary complex a(omega), define nonnegative weighted start
histograms

    w_i(j)=sum_(omega:j_i(omega)=j) |a(omega)| product_l V_l(omega).

**Theorem.** Uniformly in every such finite family,

    |sum_omega a(omega)[product_i S_i(omega)-product_i T_i(omega)]|
      =O_(A,epsilon,r)(sqrt(X) Y^r L^(-A) sum_i ||w_i||_(ell^2(J))).    (6)

The positions, phases, progressions and weights may be coupled in any
way, including identical starts in every slot. The histograms count all
multiplicities. No restriction on the number of parameter labels is
implicit; it appears through their actual weighted concentration.

### Proof

Both |S_i| and |T_i| are at most C Y L V_i by their pointwise arithmetic
weights. The exact product identity is

    product_i S_i-product_i T_i
      =sum_i (S_i-T_i) product_(l<i) S_l product_(l>i) T_l.

Apply (4) in its difference slot. Taking absolute values and summing
gives a bound

    (C Y L)^(r-1) sum_i sum_j w_i(j)D(j).

Cauchy--Schwarz bounds each last sum by ||w_i||_2 ||D||_2.
Use (3) with enough logarithmic saving to absorb L^(r-1). This proves
(6). It controls deterministic products of single-prime sums. It does
not assert that a multiple-prime correlation on a common position variable
factorizes into those sums.

### When the bound vanishes after normalization

Write W=sum_omega |a(omega)| product_l V_l(omega), so sum_j w_i(j)=W
for each i. If W>0, define M_i=W^2/||w_i||_2^2. This effective number of
occupied starts lies between one and the number of starts supporting w_i.
After division by the scale Y^r W, (6) is

    O_(A,epsilon,r)(L^(-A) sum_i sqrt(X/M_i)).                         (7)

Thus M_i>=X L^(-C) at every slot, for any fixed C, supplies arbitrary
logarithmic saving after choosing A large enough. Uniform weights on M
starts have M_i=M; concentrating all weight at one start has M_i=1.
For a fixed support size X^theta with theta<1, the bound in (7) does not
prove a vanishing error: its power of X cannot be absorbed by a fixed
logarithmic saving. This is a limitation of this deduction, not a theorem
that such particular windows fail.

## 5. A longer-window pointwise alternative

[Matomaki, Shao, Tao and Teravainen, Higher uniformity of arithmetic
functions in short intervals I, arXiv:2204.03754v4,
Theorem 1.1(ii)](https://arxiv.org/pdf/2204.03754v4)
supplies the analogous maximal Fourier discrepancy at every base j when

    j^(5/8+epsilon)<=Y<=j^(1-epsilon).

In this version the model is Lambda_sharp,j, with
R_j=exp((log j)^(1/10)). Apply the theorem separately at each j; do not
silently replace all these models by Lambda_sharp,X. For j in [X,2X],
its uniform bound is O_(B,epsilon)(Y L^(-B)). The same partial-summation
and telescoping argument then gives

    |sum_omega a(omega)[product_i S_i-product_i T_i,local]|
      =O_(A,epsilon,r)(Y^r L^(-A) W),                                (8)

provided every selected window meets the longer-window range. This
pointwise deduction needs no spread over start positions. Replacing the
varying presieved models by a common one requires another estimate.

## 6. Remaining interface

The reference MM leaf seeks a rational major-arc main term, and the
fourth/higher moment transport consumes a particular singular-series
comb and geometry. Equations (3)-(8) only replace Lambda by the stated
presieved function. An application must still establish:

- Actual window lengths are in one of the verified theorem ranges at
  each dilated lattice's own position scale; a common physical length
  does not establish this after a modulus-dependent rescaling.
- The short-window starts meet (7), or the longer pointwise range in (8).
- Weight variation, multiplicities and frequency rescalings have been
  included in the consumption norm.
- The presieved main term has the claimed arithmetic comb and continuum
  normalization, with all relevant frequencies retained.
- Any identity reducing a common-position multiple-prime expression to
  products of single sums holds for the precise weighted family.
- The remaining bands, smoothing collars and actual-minus-model tails
  have separate estimates.

In particular R bounds prime factors of a primorial; it is not a cutoff
q<=R on Ramanujan denominators. This theorem does not overturn the
[subpower denominator-cutoff obstruction](fixed_cutoff_obstruction.md).
It provides a rigorously available route around the invalid unrestricted
short-window endpoint deduction, with the exact restrictions stated.

## Finite verification

`python scripts/verify_short_window_transport.py` checks integer-cell
constancy, discrete partial summation on all small progressions, coupled
product telescoping and weighted start histograms, exceptional-set
bookkeeping and concentration examples. These are finite algebraic audits;
the short-interval analytic input is the cited literature.
