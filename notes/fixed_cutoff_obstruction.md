# Subpower frequency cutoffs leave the leading model core in the tail

For the explicit class-subtracted arithmetic model, every Ramanujan cutoff
with log P=o(log T) captures a contribution tending to zero after normalization.
The complementary tail carries the entire limiting core -1/48. As a result,
an absolute envelope for that same tail cannot have asymptotic upper bound
0.0111: its lower limit is at least 1/48, approximately 0.0208333.

This concerns a precisely specified model tail. A differently defined error,
or an actual-minus-model tail discrepancy, requires its own analysis.
It does not refute a finite-height measurement or establish a zeta-zero bound.

## Definition and finite-cutoff bound

Use the coefficients and model functional of
[the quantitative core theorem](arithmetic_core_limit.md). The pinned
reference [`mains_envelope.py`](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/scripts/mains_envelope.py)
defines

    beta_q(b_1,b_2) = mu(q_1)mu(q_2)/(phi(q_1)phi(q_2)),
    q_i=q/gcd(q,b_i),
    S={prime divisors of 2b_1b_2}.

Here mu is the Möbius function and phi is Euler's totient. Keep q<=P
whose prime divisors are not all in S. From

    c_q(j)=sum_(d|gcd(q,j)) d mu(q/d),

the finite piece has divisor coefficients

    gamma_d^(P) = sum_(q<=P, d|q, q not in the S-class)
                  mu(q/d) beta_q(b_1,b_2).

Its sawtooth transform is

    G_P(theta)=sum_(d<=P) gamma_d^(P) G_0(d theta).

Since |beta_q|<=1 and |G_0(x)|<=x for x>0,

    |G_P(theta)| <= theta sum_(q<=P) sum_(d|q) d
                    = theta B_P,

where B_P is a finite constant depending only on P. This bound is uniform
in every prime-power pair, including shared-base pairs.

Let C_ell,P be the core functional with G replaced by G_P. Its normalized
transform is bounded by B_P/ell, the overlap is bounded, and the product
modulus measure has total mass O(1). Therefore

    C_ell,P = O_P(1/ell).

This proves vanishing for every fixed P. The next section strengthens the
estimate uniformly in P.

## Uniform logarithmic cutoff cost

**Theorem.** With the same finite Ramanujan cutoff and functional, uniformly
for P>=1 as ell tends to infinity,

    |C_ell,P| = O((1+log P)/ell).                           (1)

Consequently P=exp(o(ell)), including every fixed power of log T, leaves
the entire leading model core in the complementary tail. The proof controls
the coprime and shared-base modulus pairs separately; a bound valid only for
generic large primes would not suffice.

Let sigma(n)=sum_(d|n) d and define the nonnegative multiplicative sequence

    w(n)=mu(n)^2 sigma(n)/phi(n)^2.

It is supported on squarefree integers, with
w(p)=(p+1)/(p-1)^2. There is an exact Dirichlet convolution
w=(1/id)*h_sigma. Its nonzero prime-power coefficients are

    h_sigma(p)   = (3p-1)/[p(p-1)^2],
    h_sigma(p^2) = -(p+1)/[p(p-1)^2].

Their absolute sum is 4/(p-1)^2. Thus

    sum_r |h_sigma(r)|
      = product_p (1+4/(p-1)^2)
      <= 5 exp(sum_(j>=1) 1/j^2) < 5e^2 < 37.

The inequality bounds the odd-prime sum by the sum over all odd integers
p=2j+1, and treats p=2 exactly. Convolution with the harmonic sum gives,
for every real P>=1,

    sum_(q<=P) w(q) <= 37(1+log P).                        (2)

The bound is on absolute values, so no cancellation is being assumed.

For coprime prime powers b_1,b_2, beta_q vanishes unless q is squarefree.
At an odd prime dividing one modulus, the local absolute coefficient changes
from 1/(p-1)^2 to 1/(p-1); hence the local weight
sigma(p)|beta_p(1)| is (p+1)/(p-1)<=2. There are at most two such odd
primes. At 2 the absolute weight remains sigma(2)=3, as in w.
Expand over the subsets of those at most two odd primes, and bound each
remaining positive ordinary sum by (2). Each expansion factor is at most
1+2=3. Therefore, uniformly in the exponents and primes,

    sum_(q<=P) |beta_q(b_1,b_2)| sigma(q) <= 333(1+log P).

Restricting to the free q-class can only decrease this positive sum. Since
|G_0(d theta)|<=d theta, the finite divisor expansion proves

    |G_P(theta)| <= 333 theta(1+log P)                     (3)

for every coprime pair and every theta>0.

For a shared-base pair (p^a,p^b), put e=min(a,b)>=1. At the prime p,
beta_p(v)=1 for 0<=v<=e, |beta_p(e+1)|<=1/(p-1), and higher terms vanish.
Using sigma(p^v)<=p^(v+1)/(p-1) shows

    sum_(v>=0) |beta_p(v)| sigma(p^v)
      <= 2 [p/(p-1)]^2 p^e <= 8 p^e.

All other primes have the ordinary weights; remove p and apply (2). This
also covers p=2. Thus

    |G_P(theta)| <= 296 theta p^min(a,b)(1+log P).          (4)

For the coprime contribution, divide (3) by theta ell and use the bounded
overlap H and total modulus measure mu_ell([0,1])=O(1). This is
O((1+log P)/ell).

For the shared-base contribution, divide (4) by theta ell and include both
modulus weights Lambda(b_i)/(ell b_i). Enlarging the allowed nu interval
to length one bounds its absolute value by a fixed multiple of

    (1+log P)/ell^3
      sum_(p, a,b>=1; p^a,p^b<=exp ell) (log p)^2/p^max(a,b).

There are 2k-1 ordered exponent pairs with max(a,b)=k. The k=1 sum is
O(ell^2) by Mertens' estimate and partial summation. The k>=2 sum is O(1):

    sum_(k>=2) (2k-1)/p^k
      = (3p-1)/[p(p-1)^2] <= 12/p^2,

and sum_p (log p)^2/p^2 converges. This proves (1). All constants are
independent of P, including when P grows with T.

## The complementary tail

Define G_tail,P=G-G_P and let C_ell,tail,P be its signed core functional.
The quantitative core theorem gives

    C_ell,tail,P = C_ell-C_ell,P
                  = -1/48 + O((1+log P)/ell).

This tends to -1/48 whenever log P=o(ell). For P=(log T)^B with fixed
B>0, the error is O_B(log log T/log T).

For large ell the normalization prefactor is positive. The modulus weights
and overlap H are nonnegative. Hence the absolute tail functional

    A_ell,tail,P = 2(ell/ell_1)^4 integral
                    |G_tail,P(theta)/(theta ell)| H
                    dnu dmu_ell dmu_ell

satisfies A_ell,tail,P>=|C_ell,tail,P|. In particular,

    liminf_(ell->infinity) A_ell,tail,P >= 1/48 > 0.0111
    whenever log P=o(ell).

Any pointwise or integrated absolute envelope dominating this same model
tail inherits the lower limit. Thus a fixed P=40 finite-height estimate
cannot certify a uniform asymptotic absolute-tail charge of 0.0111 for
this convention. A signed estimate can exploit the negative tail; an
absolute estimate loses that cancellation. Replacing 40 by a fixed power
of log T does not repair this obstruction. If the finite piece is to tend
to the full nonzero core, (1) requires log P to be at least a positive
constant times ell for all sufficiently large ell. This is a necessary
growth condition, not a sufficient approximation theorem.

## Relationship to the full frequency series

The coefficient definition also makes the infinite model transform rigorous.
Let Gamma be the full Euler coefficients, eta the class coefficients and
gamma=Gamma-eta. Their divisor expansions are

    D(j)=sum_(d|j) d Gamma_d,
    Dbar(j)=sum_(d|j) d eta_d.

The local identity sum_(v<=v_p(j)) p^v f_p(v)=kappa_p(v_p(j)) verifies the
singular-series density formula. Thus D-Dbar is the class-subtracted density.
For each finite J, rearranging the finite divisor sum gives

    2 sum_(j<=J) [D(j)-Dbar(j)]
      [sin(2theta j)-sin(theta j)]/j
      = sum_(d<=J) gamma_d 2 sum_(h<=J/d)
          [sin(2d theta h)-sin(d theta h)]/h.

Partial sine-harmonic sums are bounded uniformly in their length and angle.
For completeness, reduce the angle to t in [0,pi]. The initial h<=1/t
terms have absolute sum at most one. Geometric partial sums in the remaining
range are at most pi/t; summation by parts bounds the weighted remainder
by 2pi. At t=0 the sum is zero. Consequently each inner difference above
has a fixed uniform bound. Since sum_d|gamma_d|<=24, dominated convergence
then proves that the ordinary J->infinity limit equals sum_d gamma_d G_0(d theta).
The finite q piece admits the same calculation. The obstruction therefore
concerns the full model frequency tail, rather than a formal rearrangement
of an unproved conditionally convergent sum.

## Implication for the remainder ledger

The full model-core theorem already includes its infinite frequency tail.
An arithmetic proof that subtracts this full model should bound the mismatch
between actual and modeled tails, or retain a rigorously controlled signed
tail. It must identify the error being charged. Adding an absolute envelope
for a fixed model tail to a separate full model-core calculation does not
make a 0.0111 asymptotic assertion valid.

The algebraic 70.3054% conversion from an assumed allowance of 0.0111 remains
a conditional calculation. It is not a consequence of an established tail
bound for the cutoff defined here. The fixed-P route needs a changed estimate
or a precisely different remainder before it can support that application.

```sh
python scripts/verify_arithmetic_transport.py
```

The standard-library checker verifies independent Ramanujan identities,
finite-cutoff decompositions, the new harmonic convolution and the exact
weighted coefficient bounds for coprime and shared-base prime powers.
The asymptotic obstruction follows from the proof above, not from
finite-height numerics.
