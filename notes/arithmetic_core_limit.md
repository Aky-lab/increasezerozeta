# A quantitative limit for the class-subtracted arithmetic core

This note proves the limit of an explicitly defined prime-power/sawtooth
functional, with an O(1/log T) error. It supplies the coefficient asymptotic,
uniform exceptional-class estimates and exact overlap geometry needed in
[the arithmetic-to-continuum ledger](arithmetic_to_continuum.md).
The reduction of an actual zeta moment to this functional, and the separate
fixed-P tail envelope, remain open obligations.

## 1. Coefficients and the functional

The convention is pinned to
[JoshuaHKU/zeta-0.7947-reproduction, revision d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8).
The definitions are `_beta_loc` in
[mains_envelope.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/mains_envelope.py),
`_loc` and `_gamma_free_exact` in
[tail_bound.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/tail_bound.py),
and the twelve-overlap geometry in
[m1_suite.py](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/m1_suite.py).

For a prime p, let e_i=v_p(b_i). Define beta_p(v) as the product of two
factors, one for each e_i: the factor is 1 if v<=e_i, -1/(p-1) if
v=e_i+1, and zero if v>=e_i+2. The local coefficient is

    f_p(v)=beta_p(v)-beta_p(v+1).

Write e=min(e_1,e_2), alpha=1/(p-1). The only nonzero local coefficients
are at v=e,e+1, with vectors

    e_1=e_2: (1-alpha^2, alpha^2),
    e_1!=e_2: (1+alpha, -alpha).

These statements include p=2 and all exponents. Let S be the prime divisors
of 2b_1b_2, eta the finite tensor over S, and R the tensor of ordinary
positive vectors over odd primes outside S. As multiplicative coefficient
measures the full and class-subtracted sequences are

    Gamma^(b_1,b_2) = eta * R,
    gamma^(b_1,b_2) = eta * R - eta.

The convolution has disjoint prime supports. This subtraction is the exact
class replacement; omitting it changes the convention.

For prime-power moduli, S has at most two odd primes. Each odd local vector
has absolute sum at most two, and the vector at 2 has absolute sum at most
three. The ordinary tensor has total mass one. Therefore, uniformly in all
prime-power pairs and exponents,

    sum_d |gamma_d^(b_1,b_2)| <= 24.                         (1)

Use the sawtooth function G_0 from the reference convention. It satisfies
G_0(x)=-x for 0<x<pi and |G_0(x)|<=pi everywhere. In particular,
|G_0(x)|<=min(x,pi) for x>0. Put

    G_(b_1,b_2)(theta)=sum_d gamma_d^(b_1,b_2) G_0(d theta).

This series converges absolutely by (1).

Let ell=log(T/(2pi)), ell_1=ell+2log 2-1, theta=2pi exp(-(nu-1)ell),
and

    mu_ell = (1/ell) sum_(b<=exp ell) Lambda(b)/b delta_(log b/ell),
    D_nu = {nu-1 <= beta_1,beta_2 <= nu/2},  1<=nu<=2.

For beta_i in D_nu, let H be the sum of the twelve unit-length overlaps
of the source `qvec`, with legs (beta_1,nu-beta_1) and
(beta_2,nu-beta_2). Explicit formulas are proved in section 5. Define

    C_ell = 2 (ell/ell_1)^4 integral_(nu=1)^2
            integral_(D_nu) [G_(b_1,b_2)(theta)/(theta ell)] H
            dmu_ell(beta_1) dmu_ell(beta_2) dnu.             (2)

At the atoms, b_i=exp(ell beta_i). This is a model functional with its
normalization specified exactly, rather than a claim about an unproved
prime-correlation or zeta-moment replacement.

**Theorem.** As ell tends to infinity,

    C_ell = -1/48 + O(1/ell).                              (3)

The constants are independent of T. The proof uses the classical Mertens
estimate sum_(p<=x) log(p)/p=log x+O(1), but no conjecture about prime pairs.

## 2. Prove the universal coefficient asymptotic

Write C_2=product_(p>2)(1-1/(p-1)^2). For b_1=b_2=1, the full sequence
is supported at d=2m, with m odd and squarefree, and

    Gamma_(2m)=C_2 product_(p|m) 1/[p(p-2)].

Set a(m)=product_(p|m)1/(p-2) on odd squarefree m, and zero elsewhere.
There is an exact Dirichlet convolution a=(1/id)*h, where h is multiplicative
and its nonzero local values are

    h(1)=1,
    h(2)=-1/2,
    h(p)=2/[p(p-2)], h(p^2)=-1/[p(p-2)] for p>2.

All other prime-power coefficients vanish. This follows by multiplying the
local series for a by the local factor 1-p^(-s-1). Moreover,

    sum_r |h(r)| <= (3/2) product_(p>2)(1+3/[p(p-2)])
                  <= (3/2) exp(3/2) < 8.

The last inequality uses the telescoping sum over all odd integers
sum_(j>=1) 1/[(2j+1)(2j-1)]=1/2. The same Euler product shows
sum_r |h(r)| r^alpha<infinity for every 0<alpha<1/2, as well as absolute
convergence after weighting by log r.

The sums needed for the constant term are

    sum_r h(r)=1/(2C_2),
    sum_r h(r)log r= -log 2/(2C_2).

For each odd prime, h(p)+2h(p^2)=0, so its logarithmic contribution is zero.
Only h(2) contributes to the second identity. Absolute convergence justifies
factoring and differentiating these products.

With H_k the harmonic numbers, convolution gives

    sum_(m<=X) a(m)=sum_(r<=X) h(r) H_floor(X/r).

Use H_floor(y)=log y+EulerGamma+O(1/y). The weighted absolute convergence
above bounds the harmonic-number error and the omitted tails by
O_alpha(X^(-alpha)(1+log X)). Thus

    sum_(d<=M) d Gamma_d^(1,1)
      = log M + EulerGamma + O_alpha(M^(-alpha)(1+log M)).

The class subtraction is eta=delta_2, hence for M>=2,

    A(M):=sum_(d<=M) d gamma_d^(1,1)
      = log M + EulerGamma - 2
        + O_alpha(M^(-alpha)(1+log M)).                    (4)

This derives the coefficient directly. It is not a fitted
numerical slope. Partial summation also gives

    sum_(d>M) gamma_d^(1,1)
      = 1/M + O_alpha(M^(-1-alpha)(1+log M)).              (5)

For M>2 these tail coefficients are nonnegative. Taking M=floor(1/theta),
using G_0(d theta)=-d theta in the initial range, and bounding the tail by
pi times (5), proves

    G_(1,1)(theta) = -theta log(1/theta)+O(theta).         (6)

The bound extends uniformly to 0<theta<=2pi: on any compact range bounded
away from zero it follows from absolute summability and bounded G_0.

## 3. A uniform estimate for coprime prime powers

When b_1,b_2 are coprime prime powers, every special-prime vector has e=0.
Its coefficients are independent of the positive exponent of the modulus.
The weighted absolute sum at an odd special prime is

    (1+alpha)+p alpha=2p/(p-1)<=3.

At 2 it is at most four. Therefore sum_s s|eta_s|<=36.

The universal full sequence has weighted partial sums bounded by
16(1+log M), from the absolute h-convolution estimate. Removing up to two
odd generic factors increases a coefficient by at most (4/3)^2, because
1-1/(p-1)^2>=3/4. Shifting from d=2m to the ordinary odd tensor gives

    sum_(m<=M) m R_m <= (8/9) sum_(d<=2M) d Gamma_d^(1,1)
                      <= 32(1+log M).

Using gamma=eta*R-eta, for every M>=1,

    sum_(d<=M) d|gamma_d^(b_1,b_2)| <= 1200(1+log M).     (7)

Partial summation applied to the absolute coefficients yields the tail
bound 1200(log M+2)/M. Split the sawtooth series at floor(1/theta), which
is at least 1/(2theta) when theta<=1. With K=1200(1+4pi), this proves

    |G_(b_1,b_2)(theta)| <= K theta(1+log_+(1/theta))      (8)

for all theta>0, uniformly in the coprime prime-power pair. For theta>=1,
the simpler bound (1) proves the same inequality. In the scaled variables,

    |G_(b_1,b_2)(theta)/(theta ell)|
      <= K [nu-1+1/ell],  1<=nu<=2.                     (9)

For distinct ordinary odd primes in D_nu, the
[class-subtracted local comparison](class_subtracted_universality.md)
has error at most 12theta, since min(p,q)>=2pi/theta. Combining this
with (6) gives, uniformly on the whole contributing region,

    G_(p,q)(theta)/(theta ell)=-(nu-1)+O(1/ell).          (10)

No separate fixed-width strip at nu=1 is needed.

## 4. Remove exceptional moduli with the required bounds

Mertens' estimate and summability of higher prime powers imply

    mu_ell([0,1])=1+O(1/ell),
    mu_ell({genuine higher prime powers})=O(1/ell).

The atom b=2 also has mass O(1/ell). For coprime pairs involving either
kind of exceptional modulus, (9) and bounded H therefore give a total
O(1/ell) contribution to (2). Small modulus mass is used only after the
uniform transform estimate has been proved.

Pairs sharing a base prime require a different argument. Their local
vectors are shifted by min(e_1,e_2), so (7) does not hold with a uniform
constant. Use (1) instead. The region imposes

    nu<=1+log(min(b_1,b_2))/ell.

Consequently

    integral_(allowed nu) theta^(-1) dnu
      <= min(b_1,b_2)/(2pi ell).

Multiplying this by the normalized modulus weights, by 24pi/ell from (1),
and by bounded H, gives a shared-base contribution bounded by a constant
times

    ell^(-4) sum_(p^a,p^b<=exp ell)
               (log p)^2 / p^max(a,b).                  (11)

There are 2k-1 ordered exponent pairs with max(a,b)=k. The k=1 sum is
O(ell^2), by sum_(p<=exp ell)(log p)^2/p<=ell sum_(p<=exp ell)log p/p.
For k>=2,

    sum_(k>=2) (2k-1)/p^k <= 12/p^2,

and sum_p(log p)^2/p^2 converges. Thus (11) is O(ell^(-2)). This includes
ordinary diagonal primes and every shared-base higher-power pair.

The limiting core is therefore carried by distinct ordinary odd primes,
with exceptional error O(1/ell). The measure-mass estimates alone would
not justify this conclusion for the unbounded normalized transform.

## 5. Derive the overlap integral

In one source position assignment, write a+c=b+d=nu, with all four legs
in [0,1]. The first of the three overlaps has range max(a,b,c,d).
Summing it over the four swaps of the two leg pairs gives

    4(1-nu+min(beta_1,beta_2)).

The other ranges are max(c+d,b) and max(c+b,d). Swapping b and d exchanges
them. Hence the pointwise ordered geometry is

    H = 4(1-nu+min(beta_1,beta_2))
        + 2 sum_(c in {beta_1,nu-beta_1})
            sum_(b in {beta_2,nu-beta_2}) (1-max(nu+c-b,b))_+.

Put t=(2-nu)/2, beta_i=nu-1+x_i, with x_i in [0,t]. The first term
integrates to 4 integral_[0,t]^2 min(x_1,x_2) dx=4t^3/3.

For the second term, the four swaps cover the square [nu-1,1]^2 in (c,b),
with unit Jacobian. Its integral is twice the volume under
(1-max(nu+c-b,b))_+. Introduce a height u between max(nu+c-b,b) and 1,
then set

    h=c-(nu-1), r=1-b, s=1-u.

The region becomes h>=0, r>=s>=0, h+r+s<=3-2nu. For 1<=nu<=3/2,
the original upper bounds on c and b are redundant. The condition r>=s
bisects a three-dimensional simplex, giving volume (3-2nu)^3/12.
For nu>3/2 the region is empty. Therefore

    Q(nu):=integral_(D_nu) H d beta_1 d beta_2
      = (2-nu)^3/6 + (3-2nu)_+^3/6,
    integral_1^2 Q(nu) dnu=1/16,
    integral_1^2 (nu-1)Q(nu) dnu=1/96.                  (12)

This independently derives the geometry used in the old core ledger.

## 6. Quantitative measure passage and conclusion

Mertens' estimate for von Mangoldt weights gives the uniform discrepancy

    sup_(0<=x<=1) |mu_ell([0,x])-x| = O(1/ell).           (13)

Indeed the numerator is sum_(b<=exp(ell x))Lambda(b)/b=ell x+O(1), uniformly
in x, including the bounded initial range. H is bounded and a finite sum of
piecewise affine overlaps, so its Lipschitz constants in either beta coordinate
are uniformly bounded in nu. For a continuous piecewise affine g on [a,b],
integration by parts against sigma=mu_ell-Lebesgue gives

    |integral_[a,b] g dsigma|
      <= ||F_sigma||_infinity (2||g||_infinity+integral_a^b |g'|).

Using the left limit F_sigma(a-) includes endpoint atoms. Apply this estimate
successively in beta_1 and beta_2 on D_nu, and integrate in nu. Equation (13)
then implies an O(1/ell) error for the product-measure replacement, including
the moving boundaries. Combining (10), the exceptional estimates, and (12),

    C_ell = -2 integral_1^2 (nu-1)Q(nu) dnu + O(1/ell)
           = -1/48 + O(1/ell).

The prefactor (ell/ell_1)^4 differs from one by O(1/ell). This proves (3).

## Verification, attribution and remaining obligations

```sh
python scripts/verify_arithmetic_transport.py
```

The standard-library verifier independently compares the beta-difference
coefficients with the singular-series density formula and shifted local
vectors in 4,550 exact cases. It checks 2,000 Euler-convolution coefficients,
637 rational overlap configurations against the original twelve-term formula,
and exact polynomial integration of (12). These checks support the finite
identities; the asymptotic proof is the argument above.

The classical Mertens estimate is stated and developed in
[MIT 18.785, Fall 2021, Problem Set 9, Problem 1](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_pset9.pdf).
The coefficient convention and geometry are credited to the pinned reference
implementation. The convolution asymptotic, uniform exceptional-class bounds
and quantitative assembly are derived here; priority for this application
requires further literature review.

This theorem removes a model-core convergence allowance. It does not bound
the separate fixed-P frequency tail by 0.0111, prove the covered-zone/band/glue
estimates, or identify the model functional with an actual zeta moment.
Those steps cannot be absorbed into the O(1/ell) term in (3).
