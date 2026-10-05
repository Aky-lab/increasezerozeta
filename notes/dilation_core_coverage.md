# How much of the arithmetic core small dilations cover

The [arithmetic-core theorem](arithmetic_core_limit.md) gives
C_ell=-1/48+O(1/ell) for its explicitly defined prime-power model.
This note computes the contribution from a cutoff on the two moduli,
and then on their reduced dilations. It quantifies the range that a
prime-transport theorem must reach to recover this leading model term.
It does not identify that model with an actual zeta moment.

## 1. A uniform cutoff theorem

Retain all definitions and the class subtraction in equation (2) of the
core note. For B=exp(delta ell), restrict its two modulus atoms to
b_1,b_2<=B, and call the resulting functional C_ell[B]. Uniformly for
0<=delta<=1/3,

    C_ell[B] = -59 delta^5/240 + O(1/ell).                    (1)

Alternatively let g=gcd(b_1,b_2), r=b_1/g and q=b_2/g, and restrict
max(r,q)<=R. Call this functional C_ell^red[R]. For R=exp(delta ell),
uniformly on the same range,

    C_ell^red[R] = -59 delta^5/240 + O(1/ell).                (2)

Thus the ratio of the restricted leading coefficient to -1/48 is

    59 delta^5/5.                                          (3)

At delta=1/3 this ratio is 59/1215, approximately 4.856%. For any
log R=o(ell), the restricted contribution tends to zero. This includes
every slowly growing dilation cutoff in
[weighted short-window transport](weighted_dilated_prime_transport.md).
The complement then retains -1/48+o(1).

This cutoff concerns the physical moduli and reduced dilations. It is
distinct from the Ramanujan denominator cutoff in the
[fixed-cutoff obstruction](fixed_cutoff_obstruction.md). Neither statement
by itself proves a signed actual-minus-model estimate for zeta moments.

## 2. Integrate the original overlap geometry

Write nu=1+u. The contributing modulus region after the cutoff is

    0<=u<=delta,  beta_1=u+x, beta_2=u+y,
    0<=x,y<=t=delta-u.

The original upper bound beta_i<=nu/2 is redundant for delta<=1/3.
Section 5 of the core note derives the ordered twelve-overlap geometry:

    H=4(1-nu+min(beta_1,beta_2))
      +2 sum_(c in {beta_1,nu-beta_1})
         sum_(b in {beta_2,nu-beta_2}) (1-max(nu+c-b,b))_+.

The four terms inside the sum become, respectively,

    (y-x-u)_+,  y,  0,  min(x-y-u,y)_+.

For the second term this uses
2u+x+y<=2delta<=1-delta<=1-y. For the first, nu+c-b
dominates b since x>=0 and u+2y<=2delta<1. The third has
nu+c-b=2-x-y>1. Consequently, throughout the restricted square,

    H=4min(x,y)+2[(y-x-u)_+ + y + min(x-y-u,y)_+].           (4)

Integrating the first two elementary terms gives 4t^3/3 and t^3/2.
The triangular wedge for (y-x-u)_+ has integral (t-u)_+^3/6.
For the last term introduce a height z>=0 with z<=y and
z<=x-y-u. Put y=z+v and x=u+2z+v+h. Then z,v,h>=0 and
2z+v+h<=t-u, a simplex with Jacobian one; its volume is
(t-u)_+^3/12. Thus

    Q_delta(1+u):=integral_[0,t]^2 H dx dy
      = (7/3)(delta-u)^3 + (1/2)(delta-2u)_+^3.             (5)

The limiting transform in the core note is -u, so the restricted
continuum contribution is

    -2 integral_0^delta u Q_delta(1+u) du
      = -2[(7/3)delta^5/20 + (1/2)delta^5/80]
      = -59 delta^5/240.                                  (6)

The formula is asserted only for delta<=1/3. Beyond that threshold
the maximum in the second term of (4) has another branch. For example,
at delta=1/2,u=0 the original full square has integral 1/3, whereas
extending (5) would give 17/48.

## 3. Uniform passage from prime powers to the integral

For distinct ordinary odd primes the core note proves, uniformly on
the whole original region,

    G_(b_1,b_2)(theta)/(theta ell)=-(nu-1)+O(1/ell).

Its coprime higher-power and prime-2 exceptions contribute O(1/ell)
in absolute value. Shared-base pairs contribute O(ell^(-2)) in absolute
value over the full domain. Restricting any of these sums cannot increase
those absolute estimates.

The weighted prime-power measure obeys
sup_x |mu_ell([0,x])-x|=O(1/ell). In each beta coordinate H has uniformly
bounded piecewise affine variation. Truncating the integration interval
at delta adds at most one bounded jump. The integration-by-parts bound
from section 6 of the core note therefore applies with a constant
independent of delta, including endpoint atoms and an empty interval.
Successive replacement of the two measures has error O(1/ell).
The prefactor (ell/ell_1)^4=1+O(1/ell) is also uniform here.
This proves (1) from (6).

For prime-power atoms with different base primes, g=1, so the cutoff
max(r,q)<=R is exactly b_1,b_2<=R. The only difference between the two
restricted domains consists of shared-base pairs. Their full absolute
contribution is O(ell^(-2)), as just noted. This proves (2).

The error term need not be smaller than delta^5 when delta tends to zero.
Equation (1) still proves a vanishing contribution for every subpower
cutoff. It does not assert a relative asymptotic in that regime.

## Verification and interpretation

`python scripts/verify_dilation_coverage.py` integrates the original twelve
overlaps by exact rational polygon clipping, independently of (4)-(5).
It also checks the reduced formula pointwise, integrates (6) exactly,
tests the branch threshold, and checks the reduced prime-power cutoff
classification. These are finite identities; the uniform asymptotic
passage is the proof in section 3.

The reference supplies the coefficient convention and original overlap
geometry; the core note credits that source and the classical Mertens
estimate. The cutoff polynomial and its assembly are derived here.
Independent review and priority assessment remain necessary. The result
provides a concrete range requirement for future arithmetic transport:
subpower reduced dilations cannot recover this model's leading term.
