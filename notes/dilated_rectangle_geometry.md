# Dilated rectangle locks: exact geometry and local normalization

The [resolved lock frame](resolved_lock_frame.md) gives the physical
four-point count. This note derives its heterogeneous cube, its exact
presieved main term, and a uniformly bounded local second-moment budget
for prime-power moduli. It also states a deterministic progression-norm
criterion at the correct physical scale. The prime estimate must control
those progression norms and the actual consumption weights. The subsequent
[weighted transport theorem](weighted_dilated_prime_transport.md) supplies
this for fixed coefficients and a sufficiently slowly growing prime-power
family, with explicit start-weight and consumption conditions.

## 1. Remove the lock without collapsing the variance

Let b1,b2 be positive integers, g=gcd(b1,b2), r=b1/g, q=b2/g.
Then gcd(r,q)=1. Write j=gJ and use zero-extended finite real functions
f,g0 on integer intervals I_m=(a,a+Y_m], I_n=(b,b+Y_n]. Define

    N_(f,g0)(k,J)=sum_(r n-q m=J)
                    f(m)f(m-rk)g0(n)g0(n-qk).

The original physical lock is zero unless g divides j. Two solutions
of the same reduced lock differ by (m'-m,n'-n)=(rt,qt). Consequently

    sum_(k,J) N_(f,g0)(k,J)^2
      =sum_(m,n,k,t) f(m)f(m-rk)f(m-rt)f(m-r(k+t))
                      g0(n)g0(n-qk)g0(n-qt)g0(n-q(k+t)).           (1)

All sums are unrestricted integer sums with the functions enforcing the
windows. This is a heterogeneous eight-vertex cube. It retains the
resolved squared count, unlike the collapsed autocorrelation statistic.
The identity also holds on a finite cyclic group when r,q generate the
unit ideal in that group; the continuous parameterization is not needed.

Put

    D=min(floor((Y_m-1)/r),floor((Y_n-1)/q)).

Every nonzero summand has |k|,|t|<=D. Thus the number of possible tuples
in (1) is at most Y_m Y_n (2D+1)^2. In the balanced case
Y_m approximately rH, Y_n approximately qH, its scale is rq H^4.
For the reference's proportional windows Y_n approximately (q/r)Y_m,
H approximately Y_m/r and the scale is q Y_m^4/r^3. These are raw
counts; the reference's consumption-normalized volume is a separate
quantity and cannot be inserted in (1) without its explicit map.

## 2. A deterministic progression Gowers criterion

For residues 0<=c<r, 0<=d<q set

    F_c(z)=f(c+rz),    G_d(z)=g0(d+qz).

Let u,v be replacement sequences supported in the same respective windows
and define U_c,V_d similarly.
Translate each of these progression supports by a common integer for
its prime/model pair into [0,H], where
H=max(ceil(Y_m/r),ceil(Y_n/q))+2. Embed in Z/MZ with M>8H and M=O(H).
Use the normalized cyclic U^3 norm. Define

    delta=max_c,d {||F_c-U_c||_(U^3), ||G_d-V_d||_(U^3)},
    B=max_c,d {||F_c||_(U^3),||U_c||_(U^3),
                 ||G_d||_(U^3),||V_d||_(U^3)}.

Then the exact deterministic bound is

    sum_(k,J) |N_(f,g0)(k,J)-N_(u,v)(k,J)|^2
      <=16 rq M^4 delta^2 B^6.                                  (2)

To prove it, telescope the four factors and square with factor four.
Each of four squared terms has two discrepancy vertices and six ordinary
vertices in (1). Split m=c+rz, n=d+qy. For each of the rq residue pairs,
use cube variables z,-k,-t,y-z; the four F vertices and four G vertices
form the two layers of an ordinary cube. Vertex-dependent conjugations
do not change U^3 norms. The
[Gowers--Cauchy--Schwarz inequality, equation (5.5)](https://people.maths.ox.ac.uk/tillmann/GreenTao.pdf) bounds
each mixed cube by M^4 delta^2 B^6. The choice of M prevents aliasing:
all supported base positions lie in [0,H], each increment is a difference
of such positions, and no sum of the increments wraps into the support.
Summing the four telescoped terms proves (2).

In balanced windows this matches the raw rq H^4 squared scale. If B is
bounded and delta tends to zero, it supplies the required replacement at
that scale. The assertion is conditional on these progression norms;
small ordinary Fourier error does not supply them.

### The cost of extracting an individual progression

There is also a sharp general loss when passing from an ordinary U^3
norm to one progression norm. On Z/MZ with d|M, let
h_c(z)=h(c+dz) on Z/(M/d)Z. The subgroup cube identity is

    ||h 1_(x=c mod d)||_(U^3(Z/MZ))^8
      =d^(-4)||h_c||_(U^3(Z/(M/d)Z))^8.

The coset indicator is an average of d linear characters. Linear
characters preserve U^3 norms, so triangle inequality yields

    max_c ||h_c||_(U^3)<=sqrt(d)||h||_(U^3).                      (2a)

This factor is attained by h=1_(x=c mod d): its global norm is
d^(-1/2) and its selected progression is identically one. Thus global
U^3 smallness alone cannot guarantee uniformly small progression norms
when d grows. The example is a general bounded nonnegative sequence,
not a counterexample to a separately proved prime-specific estimate.
For intervals, the support-density normalizations must additionally be
tracked before applying this extraction.

## 3. Exact cube residue counts

For a prime p, let B_p(r,q) count the (m,n,k,t) modulo p for which the
eight forms in (1) are all nonzero. Coprimality rules out p dividing both
r and q. There are two cases:

    p not dividing rq:
      B_p(r,q)=p^4-8p^3+28p^2-44p+23  (p odd),
      B_2(r,q)=1;

    p dividing rq:
      B_p(r,q)=(p-1)(p^3-4p^2+6p-3).                            (3)

In the first case the invertible changes z=m/r, y=n/q, followed by
x=z, increments -k,-t,y-z, give the ordinary cube count certified in
[averaged rectangle transport](averaged_rectangle_transport.md).
If p divides r, all four f forms equal m. The g forms form a parallelogram
in the three variables n,k,t. Every three coefficient rows are independent
over every prime field. Inclusion-exclusion therefore gives

    C_p=p^3-4p^2+6p-3

for that parallelogram, multiplied by p-1 choices of m. The case p|q is
symmetric. This argument includes p=2 and arbitrary prime-power exponents.

For squarefree W put c_W=W/phi(W) and define

    beta_p(r,q)=B_p(r,q) p^4/(p-1)^8.

The CRT makes the full cube normalization exactly product_(p|W) beta_p.
At an exceptional prime p|rq, (3) simplifies to

    beta_p(r,q)=(1-1/p)^(-3)[1+1/(p-1)^3].                        (4)

The ordinary factors have a convergent uniform upper envelope, and so
do 1+1/(p-1)^3. For example the previous cube bound gives

    product_(p|W) beta_p(r,q)
      <= C [rq/phi(rq)]^3,
    C=(81/2)exp(25/8).                                           (5)

Indeed take the maximum of the ordinary factor and the bracket in (4).
At 2 and 3 that maximum is 16 and 81/32; at p>=5 it is at most
1+100/p^3. Extract the factors (1-1/p)^(-3) at p|rq.

If b1,b2 are prime powers, rq has at most two distinct prime divisors.
Hence rq/phi(rq)<=3, and the budget is at most 27C, uniformly in both
prime-power bases, exponents and W. Large prime-power dilations therefore
do not cause a growing residue normalization loss by themselves.
For unrestricted moduli that conclusion would fail: if r is a primorial
and q=1, (4) makes the budget unbounded as the included primes grow.

## 4. Parameterize the finite presieved main term

Fix any integers m0,n0 with r n0-q m0=J. Every solution is
m=m0+rz, n=n0+qz. Put

    L0=max(floor((a-m0)/r),floor((b-n0)/q)),
    U0=min(floor((a+Y_m-m0)/r),floor((b+Y_n-n0)/q)),
    Z_(r,q)(k,J)=max(U0-L0-|k|,0).                               (6)

The allowed z form the integer interval
(L0+max(k,0), U0+min(k,0)]. Changing the Bezout solution only translates
this interval; its length is invariant. This is the exact raw overlap
cardinality before primality is imposed.

For each p define a_p(k,J) as the proportion of z modulo p for which
m0+rz, m0+r(z-k), n0+qz, n0+q(z-k) are all nonzero. The factor

    s_(p;r,q)(k,J)=a_p(k,J)/(1-1/p)^4

is independent of the Bezout solution. Its explicit formula is:

* If p does not divide rq, let c=-J/(rq) modulo p. Then
  a_p=1-#{0,k,c,c+k}/p.
* If p divides rq and J, then a_p=0.
* If p divides rq but not J, then a_p=1-d/p, where d=1 when p|k
  and d=2 otherwise.

The exceptional formulas follow because one pair is constant, its value
is forced by J, and the other pair has invertible slope. Let
S_(w;r,q)=product_(p<=w) s_(p;r,q). For the finite presieve Lambda_w,
CRT counting on the z interval gives, uniformly in r,q,k,J and windows,

    N_(Lambda_w,Lambda_w)(k,J)
      = Z_(r,q)(k,J) S_(w;r,q)(k,J)+O(c_W^4 W).                   (7)

The error follows by counting each good z residue with error at most one;
there are at most W residues. No infinite-series replacement is asserted
on proportional or repeated affine forms.

## 5. The local main-term second moment is the same budget

For each p, even when p|rq, the equation r n-q m=J is surjective and
each solution set has p elements parameterized by z. Expanding two
copies of its admissible count shows exactly that

    (1/p^2)sum_(k,J mod p) s_(p;r,q)(k,J)^2=beta_p(r,q).           (8)

Thus the CRT second moment of the finite singular series is

    (1/W^2)sum_(k,J mod W) S_(w;r,q)(k,J)^2
      =product_(p|W) beta_p(r,q).                                (9)

For prime-power moduli this is uniformly bounded by (5). These are exact
complete-period identities. Converting them to incomplete weighted
parameter ranges still requires endpoint and weight estimates.

## 6. Remaining analytic obligations

The published [MRSTT II, Theorem 1.3(i)](https://arxiv.org/pdf/2411.05770v2)
controls ordinary and squarefree-W-tricked prime U^3 discrepancies on
almost all short intervals, with qualitative decay. It does not directly
state the uniform estimate in (2) for arbitrary growing r,q and every
progression residue. Such a transfer must track its losses, the progression
lengths and the coupled start weights. Its quantitative Fourier theorem
cannot be substituted for this missing cubic estimate.

As a scale audit, if Y_m=X^theta and r=X^beta with 0<=beta<theta<1,
the progression length and position scales are X^(theta-beta) and
X^(1-beta). Even an analogous one-third-power short-window input would
require theta-beta>(1-beta)/3 with a fixed margin, or
beta<(3theta-1)/2. The numerical scale comparison alone supplies no
prime theorem with a large progression coefficient. A contribution with
only a few surviving k cannot be discarded without its weight bound.

Equations (1)-(9) resolve the exact geometry, raw overlap and local
second-moment normalization. They do not establish the weighted dilated
prime replacement or the reference's zeta consumption normalization.
Independent review and priority assessment remain open.

## Verification

`python scripts/verify_dilated_rectangles.py` independently compares direct
locked squared counts with heterogeneous cubes, direct prime-field counts
with (3), residue-root formulas and their second moments, CRT products,
and raw window overlaps and finite-presieve main terms. It also checks
the deterministic progression criterion on exact rational finite examples,
the subgroup extraction identity and the sharp single-coset examples.
These finite checks audit the identities; the cubic inequality and any
eventual analytic prime input require their mathematical arguments.
