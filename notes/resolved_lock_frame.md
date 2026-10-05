# The Fourier frame for resolved four-point locks

The one-dimensional convolution of two autocorrelations aggregates
four-point counts. It does not retain their two lock coordinates, and its
Parseval norm is not the resolved lock variance. This note gives the
correct transforms, an exact nonnegative counterexample, a valid
full-lock transport theorem, and the cubic uniformity norm for the
rectangle restriction.

The subsequent [dilated geometry calculation](dilated_rectangle_geometry.md)
evaluates the heterogeneous cube's local normalization, raw overlap and
finite presieved main term, and states a progression Gowers criterion.

The distinction matters for the fourth-moment arithmetic interface. It
does not disprove a prime-specific estimate that might be established by
additional arguments, and it does not affect the exact continuum-model
certificates.

## 1. Physical lattices and the collapsed statistic

Let f,g be real finitely supported functions on the integers. For positive
moduli b1,b2, put l=lcm(b1,b2), r=b1/gcd(b1,b2), q=b2/gcd(b1,b2), and

    F(x)=f(x/b2)1_(b2|x),   G(y)=g(y/b1)1_(b1|y).

The raw locked count is

    N4(k,j)=sum_(b1*n-b2*m=j)
                 f(m)f(m-r*k)g(n)g(n-q*k)
            =sum_x F(x)F(x-l*k)G(x+j)G(x+j-l*k).                    (1)

Write A_F(s)=sum_x F(x)F(x-s), and similarly for G. Summing (1) over
all locks j removes the relation between the two base positions:

    sum_j N4(k,j)=A_F(l*k)A_G(l*k).                                (2)

Every common displacement in the supports of both autocorrelations is a
multiple of l. Thus, for

    rho(t)=sum_s A_F(s)A_G(s-t),

the exact identity with unrestricted signed k is

    sum_k sum_j N4(k,j)=rho(0).                                    (3)

For real sequences, rho has Fourier transform |Fhat|^2|Ghat|^2.
Its one-dimensional Parseval identity is valid for rho. Neither (2) nor
(3) identifies sum_(k,j) N4(k,j)^2 with sum_t rho(t)^2.
Restricting k to a positive range does not restore such an identity.

For example, with b1=b2=1 and f=g=1_{0,1}, direct integer counts give

    sum_(k,j) N4(k,j)^2=8,    sum_t rho(t)^2=70,
    sum_(k,j) N4(k,j)=rho(0)=6.                                   (4)

These are nonperiodic counts with zero extension, so wraparound plays no
role in the discrepancy.

## 2. Equal single-sequence spectra need not give equal lock variances

An example rules out recovering the resolved norm from the two separate
autocorrelations by an unspecified scalar normalization. On positions
0,...,4 take

    f=(1,4,4,1,2),   u=(2,5,2,2,1),   g=1_{0,1}.

Their generating polynomials factor as

    f(z)=(1+2z+z^3)(1+2z),
    u(z)=(1+2z+z^3)(2+z).

The second factors have equal modulus on the unit circle. Hence f and u
have identical autocorrelations and identical Fourier power spectra.
The function g is the same in both constructions, so rho is identical too.
But the positive shift k=1 gives

    sum_j N_(f,g)(1,j)^2=sum_x f(x)^2 f(x-1)^2=292,
    sum_j N_(u,g)(1,j)^2=sum_x u(x)^2 u(x-1)^2=220.                 (5)

All entries are nonnegative. Dividing f,u by 5 gives an example bounded
by one, with the same unequal values divided by 625. Thus the separate
power spectra discard information needed even for a single positive
resolved shift. A prime-specific bridge would have to supply that
information by another argument.

## 3. The correct one- and two-dimensional transforms

First fix a physical displacement s and set

    a_s(x)=F(x)F(x-s),   b_s(y)=G(y)G(y-s),
    N(s,t)=sum_x a_s(x)b_s(x+t).

On a cyclic group of order M, use the unnormalized transform
Hhat(a)=sum_x H(x)e(-a*x/M). Orthogonality gives

    Nhat_t(s,b)=ahat_s(-b)bhat_s(b),
    sum_t |N(s,t)|^2=(1/M)sum_b |ahat_s(-b)bhat_s(b)|^2.             (6)

The factors in (6) are transforms of paired products. They are not
transforms of the original single-position functions.

For the transform in both coordinates, introduce x,u,v,w with
s=x-u, t=v-x and w=v-s. The condition is x+w=u+v. Expanding its
indicator by additive characters gives the exact formula

    Nhat(a,b)=(1/M)sum_c
          Fhat(a-b-c)Fhat(c-a)Ghat(c+b)Ghat(-c),
    sum_(s,t) |N(s,t)|^2=(1/M^2)sum_(a,b) |Nhat(a,b)|^2.            (7)

All frequencies in (6)-(7) are modulo M. The internal convolution in (7)
occurs before the squared absolute value. Replacing it by the product
|Fhat|^2|Ghat|^2 changes the statistic.

Zero-extended finite integer sequences can be embedded without aliasing:
if their supports are contained in [0,H], choose M>4H. Then the equation
x+w=u+v modulo M is the same integer equation, and the signed displacement
and lock have unique representatives in their support ranges. The physical
restriction in (1) selects s=l*k; a restricted k,j range is a coordinate
mask on N, not a one-dimensional convolution identity.

## 4. What single-prime Fourier transport does control

For b>=2 and finitely supported real functions f1,...,fb, retain the entire
(b-1)-dimensional lock vector:

    T_f(v2,...,vb)=sum_x f1(x) product_(i=2..b) fi(x+vi).

On a cyclic group, its transform is

    That_f(a2,...,ab)=fhat1(-sum_i ai) product_(i=2..b) fhati(ai).     (8)

This is a valid factorization with coupled frequencies. Directly in
position space it gives the full-lock product identity

    sum_v T_f(v)T_g(v)=sum_s product_i [sum_x fi(x)gi(x+s)].         (9)

For two families f_i,u_i define

    d_i=max_a |fhati(a)-uhati(a)|,
    E_i=max(sum_x |fi(x)|^2, sum_x |ui(x)|^2).

**Full-lock stability theorem.** Without independence assumptions,

    ||T_f-T_u||_2^2 <= b sum_i d_i^2 product_(l!=i) E_l.             (10)

To prove it, telescope the product in (8). In a single replacement term,
bound the difference factor by d_i. The remaining Fourier energy sums
separate by Parseval. If i is not the base slot, change its frequency
variable to the sum of all frequencies; this is a bijection. That term's
norm is at most d_i^2 times the other energies. Apply
|sum_i z_i|^2<=b sum_i |z_i|^2. No dimension-dependent Fourier grid factor
is missing in (10).

For prime weights and the common presieved model on integer windows of
length Y with starts j_i in {X,...,2X-1}, the published-input
[Fourier transport](short_window_fourier_transport.md) supplies the errors
D(j_i), and E_i=O(Y L^2), L=log(3X). Any positive integer dilation of a
position lattice preserves its energy; its Fourier error is at most the
same maximal D(j_i).

For a finite parameter family Omega with nonnegative weights a(omega),
put H_i(j)=sum_(omega:j_i(omega)=j) a(omega). The pointwise bound
D(j)=O(Y L) and the arbitrary-saving mean-square estimate imply
sum_j D(j)^4=O_A(X Y^4 L^(-A)). Apply Cauchy--Schwarz to (10) and these
histograms to obtain, for every fixed A,

    sum_omega a(omega)||T_Lambda(omega)-T_sharp(omega)||_2^2
      =O_(A,epsilon,b)(sqrt(X) Y^(b+1) L^(-A) sum_i ||H_i||_2).      (11)

The windows must all meet X^(1/3+epsilon)<=Y<=X^(1-epsilon) at their
undilated prime-position scales. A large enough group embeds every finite
support; equivalently (10)-(11) hold for the integer arrays with zero
extension. Dilation does not add a favorable factor to the consumption
norm. If all starts coincide at j and each j is used once, (11) becomes
O_A(X Y^(b+1) L^(-A)).

For b=4 the rectangle locks are the slice

    (v2,v3,v4)=(-s,t,t-s).

Restriction can only decrease the squared norm, so (11) gives a bound on
that slice. But its scale is Y^5, whereas the raw rectangle array has
Y^2 coordinates with counts of size Y and hence natural squared scale
Y^4. A power of Y cannot be absorbed into an arbitrary logarithmic
saving. The restriction must therefore be analyzed, rather than silently
assigned the full-lock bound with one dimension removed.

## 5. The rectangle requires cubic uniformity

For a real function h on a cyclic group of order M,

    sum_(s,t) |sum_x h(x)h(x-s)h(x+t)h(x+t-s)|^2
      =M^4 ||h||_(U^3)^8.                                         (12)

Indeed set the three cube increments to -s,t,x'-x when expanding the
square. The eight factors are precisely the eight cube vertices.
Here the normalized Gowers norm is

    ||h||_(U^3)^8=E_(x,a,b,c) product_(omega in {0,1}^3)
                         h(x+omega1*a+omega2*b+omega3*c).

For complex functions the alternating conjugations belong in both the
rectangle and cube definitions. Repeated Cauchy--Schwarz gives the Gowers
cube inequality; see [Green and Tao, The primes contain arbitrarily long
arithmetic progressions, equation (5.5)](https://people.maths.ox.ac.uk/tillmann/GreenTao.pdf).
Telescoping four factors and applying it yields, for
real h,v bounded in absolute value by B,

    ||N_(h,h)-N_(v,v)||_2^2
      <=16 M^4 B^6 ||h-v||_(U^3)^2.                                (13)

Each telescoping term has h-v at two cube vertices in its squared norm;
the other six Gowers norms are at most B. Thus (13) specifies a sufficient
replacement norm for the raw symmetric rectangle.

The need for more than Fourier uniformity also has a nonnegative example.
For odd primes p tending to infinity let

    h_p(x)=1+cos(2*pi*x^2/p),   v_p(x)=1,   x in Z/pZ.

Then 0<=h_p<=2 and the quadratic Gauss-sum identity gives
max_a |hhat_p(a)-vhat_p(a)|<=sqrt(p)=o(p). Nevertheless,

    p^(-4) ||N_(h_p,h_p)-N_(1,1)||_2^2 -> 1/128.                  (14)

For completeness, away from s=0,t=0,s=t,s=-t, expand each cosine into
two characters. The constant choice contributes 1. The only other
terms with zero quadratic and linear coefficient in x are the two
alternating choices, each with coefficient 1/16. All terms with nonzero
quadratic coefficient are O(p^(-1/2)) after averaging over x, and the
remaining nonconstant linear terms average to zero. Hence uniformly there,

    N_(h_p,h_p)(s,t)/p-1=(1/8)cos(4*pi*s*t/p)+O(p^(-1/2)).

The four excluded lines have O(p) points and bounded normalized counts.
Finally the average of cos(4*pi*s*t/p)^2 over all s,t is
1/2+1/(2p), by character orthogonality. This proves (14). This is an
algebraic obstruction to a general Fourier-only rectangle argument,
not a counterexample about the von Mangoldt function.

## 6. Consequences for the reference reduction

The pinned [reference paper](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/paper.tex),
sections labeled `s:four` and `s:five`, uses the one-dimensional density
rho and a full-lock product identity in its proposed higher-moment
transport. The identities themselves are valid for their respective
statistics. A proof relating their model discrepancies to the restricted
N4(k,j) discrepancy is still needed. Equations (4)-(5) show that this step
cannot be inferred as a general identity from the displayed separate
autocorrelations alone. Any additional archive-specific bridge must state
its maps, weights and normalization explicitly and withstand those examples.

Published cubic uniformity is relevant, but its quantitative form must be
preserved. [MRSTT II, Theorem 1.3(i) and Remark 1.4](https://arxiv.org/pdf/2411.05770v2)
give qualitative small U^s norms for Lambda-Lambda_w as w tends to infinity,
with X large in terms of w and an arbitrarily small logarithmic
exceptional-position measure. They do not state arbitrary logarithmic
U^3 savings for Lambda-Lambda_sharp with R=exp((log X)^(1/10)). The bound
(13), applied with B=O(log X), loses B^6. Its decay cannot be concluded
from that qualitative theorem by dropping the loss or changing models.
Relative generalized von Neumann estimates, pseudorandom majorants and
the actual dilated window/weight ranges are possible additional inputs;
they have not been supplied by the present Fourier-only deduction.

The [averaged pair theorem](presieved_pair_transport.md) remains valid:
it estimates two-point autocorrelations, whose Parseval square has the
correct one-dimensional Fourier form. It does not by itself supply
cubic uniformity, restricted four-prime lock variance, or a new zeta-zero
bound.

## Verification

[Qualitative rectangle transport](averaged_rectangle_transport.md) gives
an additional averaged, undilated prime-model replacement by keeping the
other factors' U^3 norms instead of bounding them by B. Its residue
normalization is uniformly controlled. The saving remains qualitative,
and the actual weighted/dilated zeta consumption is not inferred.

`python scripts/verify_lock_frame.py` uses exact rational and Gaussian
integer arithmetic. It checks independent position counts, the two
resolved transforms, full-lock transforms and stability, cubic identities,
the nonperiodic counterexamples, dilation aggregation, and quadratic
character cancellations. No floating-point FFT identity is used as a
substitute for the target statistic.
