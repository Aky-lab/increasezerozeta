# The full dilation distribution of the arithmetic core

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
0<=delta<=1,

    C_ell[B] = c(delta) + O(1/ell),                         (1)

where the complete coefficient is

    0<=delta<=1/3:
      c(delta)=-59 delta^5/240;
    1/3<=delta<=1/2:
      c(delta)=-59 delta^5/240+(3delta-1)^5/60;
    1/2<=delta<=2/3:
      c(delta)=-1/48+2(1-delta)^4/3-14(1-delta)^5/15
               +(11delta-4)(2-3delta)^4/80;
    2/3<=delta<=1:
      c(delta)=-1/48+2(1-delta)^4/3-14(1-delta)^5/15.

Alternatively let g=gcd(b_1,b_2), r=b_1/g and q=b_2/g, and restrict
max(r,q)<=R. Call this functional C_ell^red[R]. For R=exp(delta ell),
uniformly on the same range,

    C_ell^red[R] = c(delta) + O(1/ell).                     (2)

Thus the ratio of the restricted leading coefficient to -1/48 is

    F(delta)=-48c(delta),                                  (3)

which equals 59delta^5/5 for delta<=1/3. It is continuous and strictly
increasing on [0,1], from zero to one. Some exact values are:

| Cutoff exponent delta | Leading fraction F(delta) | Approximate percentage |
|---|---|---|
| 1/3 | 59/1215 | 4.856% |
| 1/2 | 11/32 | 34.375% |
| 2/3 | 959/1215 | 78.930% |
| 3/4 | 147/160 | 91.875% |
| 4/5 | 15049/15625 | 96.314% |
| 9/10 | 15582/15625 | 99.725% |
| 1 | 1 | 100% |

Exact rational bisection locates the exponents needed to capture a
specified fraction of the leading model term:

| Fraction | Exponent lies in |
|---|---|
| 50% | (0.551320, 0.551321] |
| 90% | (0.734423, 0.734424] |
| 99% | (0.859563, 0.859564] |

Here the height is exp(ell)=T/(2pi); these exponents concern a cutoff
R=exp(delta ell), not a power of a progression-window length. They are
range requirements for recovering this leading term by covering the
restricted region. They do not rule out a different global analytic
argument that avoids regional prime transport.

At delta=1/3 this ratio is 59/1215, approximately 4.856%. For any
log R=o(ell), the restricted contribution tends to zero. This includes
every slowly growing dilation cutoff in
[weighted short-window transport](weighted_dilated_prime_transport.md).
The complement then retains -1/48+o(1).

This cutoff concerns the physical moduli and reduced dilations. It is
distinct from the Ramanujan denominator cutoff in the
[fixed-cutoff obstruction](fixed_cutoff_obstruction.md). Neither statement
by itself proves a signed actual-minus-model estimate for zeta moments.

## 2. The small-cutoff branch

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

## 3. Integrate all branches without a polynomial fit

For general 0<=delta<=1 and 0<=u<=delta the actual square has side

    t=min(delta-u,(1-u)/2).

The four swapped terms in section 2 have values

    (y-x-u)_+, min(y,1-2u-x-y)_+, 0, min(x-y-u,y)_+.

The first and third maximum choices hold on the entire square because
2t<=1-u. Only the second branch changes from the small-cutoff argument.
Write s=1-2u. Its integral B(s,t) is the volume

    0<=x<=t, 0<=z<=y<=t, x+y+z<=s.

The constraint z<=y selects half the cube [0,t]^3 cut by x+y+z<=s:
the other half is its image under exchange of y,z, and their intersection
has zero volume. Inclusion-exclusion on the three upper bounds gives

    B(s,t)=[s_+^3-3(s-t)_+^3+3(s-2t)_+^3-(s-3t)_+^3]/12.

Thus the restricted overlap integral for every cutoff is

    Q_delta(1+u)=4t^3/3+(t-u)_+^3/2
       +[s_+^3-3(s-t)_+^3+3(s-2t)_+^3-(s-3t)_+^3]/6.      (7)

At the full side t=(1-u)/2 the two terms involving
((1-3u)/2)_+^3 cancel, leaving

    Q_1(1+u)=[(1-u)^3+(1-2u)_+^3]/6,

which agrees with the independently derived full core geometry.
The desired coefficient is always c(delta)=-2 integral_0^delta u Q_delta du.

For delta<=1/2 one has t=delta-u. The cube bracket's affine arguments
are 1-2u, 1-delta-u, 1-2delta and 1-3delta+u. The first three stay
nonnegative throughout 0<=u<=delta. The fourth changes sign only if
delta>1/3. Replacing its full cubic integral by its positive part changes
the weighted integral by (3delta-1)^5/20; its coefficient -1/6 in Q
therefore changes c by +(3delta-1)^5/60. The other term has cutoff
u=delta/2. Elementary integration gives the first two branches of (1).

For delta>=1/2, put eta=1-delta and a=2delta-1. Split at u=a:
the side is full for 0<=u<=a and is delta-u for a<=u<=delta.
The term (1-2u)_+^3/6 has total weighted integral 1/480, independently
of delta. The first terms of the two regions have total weighted integral

    integral_0^a u(1-u)^3/6 du
      + integral_a^delta (4/3)u(delta-u)^3 du
        =1/120-eta^4/3+7eta^5/15.                          (8)

The two remaining nonzero corrections are

    (1/2) integral_a^delta u[(delta-2u)_+^3
                              -(1-delta-u)_+^3] du.

They both vanish if delta>=2/3. Otherwise set b=2-3delta and
u=a+v. Their respective upper limits are v=b/2 and v=b, giving

    -(a b^4)/16-3b^5/160=-(11delta-4)b^4/160.               (9)

Combining (8)-(9) and 1/480, then multiplying by -2, gives the last
two branches of (1). This derivation fixes every branch and coefficient
by integration, independently of any finite sampling.

One can also view F as the distribution function of max(beta_1,beta_2)
under the probability measure 96u H d beta_1 d beta_2 du on the full
core region. The normalization follows from integral uH=1/96.
Positivity explains monotonicity. Strict positivity on every interior
cutoff interval follows either from an open overlap region or from
the derivatives of the explicit branches:

    59delta^4;
    59delta^4-12(3delta-1)^4;
    (1-delta)^3(224delta-96)
       +(3/5)(2-3delta)^3(165delta-70);
    (1-delta)^3(224delta-96),

in their respective ranges. Each is positive in the interior.
The formulas agree at the branch endpoints. For delta>=2/3 the
remaining fraction has the exact fourth-order tail

    1-F(delta)=32(1-delta)^4-(224/5)(1-delta)^5.             (10)

## 4. Uniform passage from prime powers to the integral

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
This proves (1) from (7)-(9), uniformly over the entire cutoff range.

For prime-power atoms with different base primes, g=1, so the cutoff
max(r,q)<=R is exactly b_1,b_2<=R. The only difference between the two
restricted domains consists of shared-base pairs. Their full absolute
contribution is O(ell^(-2)), as just noted. This proves (2).

The error term need not be smaller than delta^5 when delta tends to zero.
Equation (1) still proves a vanishing contribution for every subpower
cutoff. It does not assert a relative asymptotic in that regime.

## Verification and interpretation

`python scripts/verify_dilation_coverage.py` integrates the original twelve
overlaps by exact rational polygon clipping: 20 small-cutoff configurations
and 66 configurations spanning all branches. It checks 1620 pointwise
small-cutoff overlaps, integrates (6) exactly, and checks the old branch
threshold. A separate univariate integrator splits at the actual affine
roots in (7) and compares the result with (1) at 61 rational cutoffs.
It verifies the displayed coverage fractions, rational quantile brackets
and 1536 reduced prime-power cutoff classifications. Four symbolic
polynomial derivative identities verify the displayed density formulas
and exact value/derivative agreement at each branch join. These are finite
identities; the uniform asymptotic passage is the proof in section 4.

The [standalone manuscript source](../papers/dilation_core_distribution.tex)
collects the model definition, arithmetic reduction, complete cutoff proof,
distribution plot and reproduction instructions in one review document.

The reference supplies the coefficient convention and original overlap
geometry; the core note credits that source and the classical Mertens
estimate. The cutoff polynomial and its assembly are derived here.
Independent review and priority assessment remain necessary. The result
provides a concrete range requirement for future arithmetic transport:
subpower reduced dilations cannot recover this model's leading term.
