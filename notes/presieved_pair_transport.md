# A finite presieved comb and averaged short-window prime pairs

A power-sized divisor sum approximates the presieved prime model in global
mean square. Its Fourier coefficients and pair-density main term are
explicit. Combining these facts with
[short-window Fourier transport](short_window_fourier_transport.md) gives
an averaged prime-pair theorem on integer short windows, with the usual
singular series and arbitrary logarithmic savings.

The average is over both position and shift. The theorem does not give
every prescribed short window, a pointwise twin-prime asymptotic, or the
four-prime lock variance needed in the zeta-moment reduction. Its weight
budget must still be checked in any such application.

## 1. Statements

Let X,Y be integers, X>=3, and fix epsilon>0 and 0<delta<1/6. Assume

    X^(1/3+epsilon)<=Y<=X^(1-epsilon).

Put L=log(3X), R=exp((log X)^(1/10)), P=product_(p<R) p,
c_R=P/phi(P), D=floor(X^delta), and

    f_R(n)=c_R 1_(gcd(n,P)=1),
    F_D(n)=c_R sum_(d|P, d<=D, d|n) mu(d).

F_D is an explicit finite divisor sum. It need not be nonnegative.
For j in J={X,...,2X-1}, define the window sequence
g_j(n)=Lambda(n)1_(j<n<=j+Y) and its autocorrelation

    C_j(h)=sum_n g_j(n)g_j(n+h).

For nonzero h, let S(h) be the twin-prime singular series from
[the progression theorem](prime_pair_progression_transport.md), with
S(-h)=S(h). Define

    Delta_j(h)=C_j(h)-S(h)(Y-|h|),  0<|h|<Y.

**Theorem 1 (finite comb).** For every fixed A>0, uniformly in integer
j and nonzero |h|<=Y,

    sum_(j<n<=j+Y) F_D(n)F_D(n+h)
      =Y S(h)+O_(A,delta)(Y L^(-A))+O(c_R^2 D^2).                    (1)

The expression with both positions restricted to the window has Y-|h|
in place of Y and the same error bound. In particular this is a uniform
statement about the explicit divisor model, not about prime pairs.

**Theorem 2 (actual prime pairs, position/shift average).** For every
fixed A>0,

    sum_(j in J) sum_(0<|h|<Y) |Delta_j(h)|^2
      =O_(A,epsilon,delta)(X Y^3 L^(-A)).                             (2)

The diagonal h=0 is excluded from the singular-series replacement and
requires its own arithmetic main term.

For every choice of families E_j of positive pairs (q,k) with qk<Y,
uniformly in all those choices,

    sum_(j in J) sum_((q,k) in E_j) |Delta_j(qk)|^2
      =O_(A,epsilon,delta)(X Y^3 L^(-A)).                             (3)

Repeated shifts are counted with their divisor multiplicity. Arbitrary
complex consumption weights a_(j,q,k) therefore satisfy

    |sum_j sum_((q,k) in E_j) a_(j,q,k) Delta_j(qk)|
      =O_(A,epsilon,delta)(sqrt(X) Y^(3/2) L^(-A) ||a||_2).           (4)

These are asymptotic deductions from the published uniformity input.
Constants and effective height bounds are not computed. Restricting the
shifts to h<=H<Y preserves (2)-(4) with Y^3 on the right; no improvement
to H Y^2 is asserted without another argument.

## 2. Exact Fourier coefficients

The full presieved function has the finite Ramanujan expansion

    f_R(n)=sum_(q|P) mu(q)/phi(q) c_q(n).                             (5)

For a prime p, its local factor is
p/(p-1)1_(p not dividing n)=1-c_p(n)/(p-1). Multiplication over primes
proves (5). Although the prime-factor cutoff is R, denominators in (5)
can be as large as the whole primorial P.

Use 1_(d|n)=(1/d)sum_(q|d)c_q(n) to expand the finite divisor model:

    F_D(n)=sum_(q|P, q<=D) alpha_q(D,R)c_q(n),
    alpha_q=c_R sum_(d|P, d<=D, q|d) mu(d)/d
            =c_R mu(q)/q
               sum_(s|P/q, s<=D/q) mu(s)/s.                         (6)

These coefficients generally differ from mu(q)/phi(q). Substituting
the latter coefficients at finite D would change the model. When D>=P,
the Euler product in (6) gives alpha_q=mu(q)/phi(q), as a calibration.

Period averages over n modulo P satisfy

    E_P[c_q(n)c_s(n+h)]=1_(q=s)c_q(h).

This follows by orthogonality of the primitive additive characters.
Hence the exact pair density of F_D is

    K_D(h)=E_P[F_D(n)F_D(n+h)]
           =sum_(q|P, q<=D) alpha_q^2 c_q(h)
           =c_R^2 sum_(d,e|P, d,e<=D, gcd(d,e)|h)
                         mu(d)mu(e)/lcm(d,e).                       (7)

The second version counts the common congruence n=0 mod d,
n=-h mod e. It directly yields the interval estimate

    sum_(j<n<=j+Z) F_D(n)F_D(n+h)=Z K_D(h)+O(c_R^2 D^2),              (8)

uniformly for integer j,Z>=0 and h: each compatible residue class has
count Z/lcm(d,e)+O(1), and there are at most D^2 ordered pairs.
For the same-window autocorrelation, the intersection interval has
integer length Z=Y-|h|.

## 3. A Rankin bound for the divisor approximation

Write E_D=f_R-F_D. Exact Mobius inversion gives

    |E_D(n)|<=c_R sum_(d|P, d>D, d|n) 1
             <=c_R D^(-sigma) sum_(d|gcd(n,P)) d^sigma,
    sigma=1/log R.

Expand the square and sum over 1<=n<=M. Since the number of multiples
of lcm(d,e) is at most M/lcm(d,e),

    sum_(n<=M) |E_D(n)|^2
      <=M c_R^2 D^(-2 sigma)
             product_(p<R) [1+(2p^sigma+p^(2 sigma))/p].              (9)

The four local choices are p dividing neither divisor, only one of the
two divisors, or both. Since p^sigma<=e, 2e+e^2<13, and Mertens gives
c_R=O(log R) and sum_(p<R)1/p<=log log R+O(1), define a bound eta by

    eta^2=O((log R)^15 D^(-2/log R)).

Then (9) is at most M eta^2. The same bound holds for the exact period
mean E_P[|E_D|^2], by replacing the multiple count with its exact density.
For D=floor(X^delta), eta decays faster than every negative power of L:
the negative exponential in its logarithm has size
-delta*(log X)^(9/10), while the prefactor contributes O(log log X).
Thus eta=O_(A,delta)(L^(-A)) for every fixed A.

This is an elementary global mean-square approximation. It is not a
uniform L^1 estimate on every power-short interval; making that inference
would lose the position/window scale factor again.

## 4. The finite pair density tends uniformly to the singular series

The full presieved period correlation is exactly

    S_R(h)=E_P[f_R(n)f_R(n+h)]
           =product_(p<R) [1+c_p(h)/(p-1)^2],
    E_P[f_R^2]=c_R.                                                  (10)

For even nonzero h, separate the factor at 2:

    S_R(h)=2 C_2(R) product_(2<p<R, p|h) (p-1)/(p-2),
    C_2(R)=product_(2<p<R) [1-1/(p-1)^2].

For odd h both S_R(h) and S(h) are zero once R>2. For even h,

    |log(S_R(h)/S(h))|
      =O(1/R+log(3|h|)/(R log R)).                                  (11)

Indeed the omitted ordinary Euler factors have total logarithm O(1/R).
There are at most log|h|/log R prime divisors of h at least R, each
contributing O(1/R). With |h|<=Y<=X, (11) decays faster than every
negative power of L. The coarse standard bound S(h)=O(log(3|h|))
therefore gives an absolute O_A(L^(-A)) error uniformly in these shifts.

Period translation is an isometry. Cauchy--Schwarz, (9) and (10) show

    |K_D(h)-S_R(h)|<=2 sqrt(c_R) eta+eta^2,                           (12)

uniformly for all integer h. Therefore K_D(h)=S(h)+O_(A,delta)(L^(-A))
for the nonzero shifts above. Equations (8) and (12) prove Theorem 1.
At h=0, instead K_D(0)=c_R+O(sqrt(c_R)eta+eta^2); replacing this by
the nonzero-shift singular series would be incorrect.

## 5. Averaged Fourier replacement by F_D

Let A_j be the supremum over real alpha of

    |sum_(j<n<=j+Y) [Lambda(n)-F_D(n)] e(alpha n)|.

The [short-window theorem](short_window_fourier_transport.md) bounds the
corresponding maximal Lambda-f_R discrepancy in mean square by
O_A(XY^2 L^(-A)), using
[MRSTT II, Theorem 1.1(ii)](https://arxiv.org/pdf/2411.05770v2).
The E_D contribution has supremum at most sum_(window)|E_D(n)|.
Cauchy--Schwarz in each window, followed by counting the at most Y
windows containing each n, gives

    sum_(j in J) [sum_(j<n<=j+Y)|E_D(n)|]^2
      <=Y^2 sum_(n<=3X)|E_D(n)|^2 <=3X Y^2 eta^2.

Thus, for every fixed A,

    sum_j A_j^2=O_(A,epsilon,delta)(X Y^2 L^(-A)).                   (13)

The same reasoning covers the maximal progression version and bounded-
variation weights. Moreover, for D<=Y, direct divisor counting bounds
any weighted Fourier sum of F_D by

    c_R V sum_(d<=D) (Y/d+1)
      <=c_R V [Y(1+log D)+D]=O(Y L^2 V).

Consequently the coupled-product deduction in the short-window note
also holds with F_D as the main term and the same histogram norm budget:
the extra fixed power of L is absorbed in the arbitrary saving. The finite
main term is explicit, although the lock geometry still needs analysis.

## 6. Parseval proves the actual windowed prime-pair estimate

Let u_j(n)=F_D(n)1_(j<n<=j+Y), and use the Fourier transforms
G_j(alpha)=sum_n g_j(n)e(alpha n), U_j(alpha)=sum_n u_j(n)e(alpha n).
Equation (8) at h=0 and the bound on K_D(0) give

    sum_n u_j(n)^2=O(Y L+c_R^2 D^2)=O(Y L^2),
    sum_n g_j(n)^2<=Y L^2.

Here D^2=o(Y) because 2delta<1/3 and Y>=X^(1/3+epsilon).
Parseval for the two autocorrelations gives

    sum_h |C_j(h)-sum_n u_j(n)u_j(n+h)|^2
      =integral_0^1 (|G_j|^2-|U_j|^2)^2 d alpha
      <=A_j^2 integral_0^1 (|G_j|+|U_j|)^2 d alpha
      =O(Y L^2 A_j^2).                                              (14)

Sum (14) over j and use (13). This is O_A(XY^3 L^(-A)).
For 0<|h|<Y, Theorem 1 approximates the model autocorrelation by
S(h)(Y-|h|), with error O_A(Y L^(-A)+c_R^2D^2). Summing its square
over O(XY) pairs gives

    O_A(XY^3 L^(-A))+O(XY L^4 D^4).

The last error is also O_A(XY^3 L^(-A)), since
D^4/Y^2<=X^(4delta-2/3-2epsilon) has a fixed negative exponent.
The squared triangle inequality proves (2).

## 7. Divisor families and consumption

Any E_j in (3) repeats h with multiplicity at most tau(h). To preserve
the arbitrary saving, split the pairs (j,h) at the threshold
|Delta_j(h)|=Y L^(-B). From (2) with exponent C, the number of bad
pairs is O(XY L^(2B-C)). The divisor bounds from the progression note
are sum_(h<Y)tau(h)<=Y(1+log Y) and
sum_(h<Y)tau(h)^2<=Y(1+log Y)^3. Hence the bad multiplicity mass is
O(XY L^(B-C/2+3/2)) by Cauchy--Schwarz.

The pointwise bound |Delta_j(h)|=O(Y L^2) follows by bounding Lambda
and S(h). The good contribution is O(XY^3 L^(1-2B)); the bad one is
O(XY^3 L^(B-C/2+11/2)). Choose B and then C large enough for any
desired A. This proves (3). Cauchy--Schwarz with twice the desired
saving proves (4).

The theorem averages over X integer starts and shifts spanning the full
window scale. It cannot automatically give the natural bound on a sparse
set of starts or a shorter shift interval. Large powers hidden in a zeta
consumption weight also cannot be absorbed by a fixed logarithmic saving.
This explicit budget is part of the application, not an optional condition.

## Verification and scope

`python scripts/verify_presieved_pairs.py` checks exact full-period
Ramanujan coefficients, partial-divisor coefficients, independent direct
period correlations and CRT densities, interval discrepancy bounds and
the local Euler factors of the Rankin square. Finite autocorrelation and
Fourier identities are independently checked on rational toy sequences.
These checks audit the discrete calculations; the asymptotic prime input
is the cited theorem, and the estimates are the proofs above.

The divisor cutoff D is a positive power of X. Its coefficients are (6),
and they are not those of a naive denominator-truncated twin series.
Neither the approximation nor (2) contradicts the earlier subpower
cutoff obstruction. The actual higher-moment lock variance and final
zeta counting interface remain open.

The [resolved frame audit](resolved_lock_frame.md) gives exact counterexamples
to inferring a four-prime lock variance from the separate pair spectra.
It states the valid full-lock replacement and the extra cubic-uniformity
criterion needed for the rectangle restriction.
