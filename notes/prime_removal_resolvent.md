# Deterministic removal of one prime from the resolvent

**Status:** analytic expansion in proof-draft form, with exact finite
checks. It bounds the expansion remainder but does not establish the
signed arithmetic estimate needed by the
[bounded counting certificate](bounded_resolvent_bridge.md).

Use that note's unit-bandwidth circle operator K_T, finite-rank P,
J_T=I-K_T and z=i*r. More generally the identities below hold for
every nonreal z, with eta=|Im z| and t=1-z.

## 1. Group all powers of the same prime

For an ordinary prime p<=X define the bounded self-adjoint operator

    A_p=(log p/ell1) sum_(j>=1; p^j<=X) p^(-j/2)
          *(exp(-i*T*j*log p)*K_(j*log p)
            +exp(i*T*j*log p)*K_(-j*log p)).                 (1)

Then K_T=sum_p A_p. The contraction of each gated translation gives

    ||A_p|| <=alpha_p=2*log p/(ell1*(sqrt(p)-1)),
    sum_p alpha_p^3 <=C/ell1^3.                              (2)

The constant is finite: for p>=2, sqrt(p)-1 is at least
(1-1/sqrt(2))*sqrt(p), and sum_(n>=2) (log n)^3/n^(3/2)
converges. Thus the sum of cubes tends to zero even though the
sum of norms may grow. Grouping prime powers is essential when
discussing removal of a prime phase.

Set

    R=(tI-K_T)^(-1), R_p=(tI-(K_T-A_p))^(-1).

The removed operator retains every other prime base. Both resolvent
norms are bounded by eta^(-1), independently of ||K_T||.

## 2. A finite expansion with a uniform remainder

The resolvent identity, iterated twice, is the exact formula

    R=R_p+R_p*A_p*R_p+R_p*A_p*R_p*A_p*R.                    (3)

Multiply tR=I+K_T*R by P on both sides and take its normalized
finite-dimensional trace. With m_k=Tr(P*R^k*P)/d, (3) gives

    t*m1=1+Z1+Q1+E1,
    Z1=sum_p Tr(P*A_p*R_p*P)/d,
    Q1=sum_p Tr(P*A_p*R_p*A_p*R_p*P)/d,
    E1=sum_p Tr(P*A_p*R_p*A_p*R_p*A_p*R*P)/d.                (4)

The trace of a compression has absolute value at most d times
the operator norm. Consequently

    |E1|<=eta^(-3)*sum_p alpha_p^3=O(ell^(-3)).               (5)

No infinite series or smallness of ||K_T|| is used.

For k=1,2,3 differentiate (4) k-1 times in z and divide by (k-1)!.
Write

    Z_k=sum_p Tr(P*A_p*R_p^k*P)/d,
    Q_k=sum_p sum_(u,v>=1; u+v=k+1)
                         Tr(P*A_p*R_p^u*A_p*R_p^v*P)/d,
    E_k=sum_p sum_(u,v,w>=1; u+v+w=k+2)
                    Tr(P*A_p*R_p^u*A_p*R_p^v*A_p*R^w*P)/d.

The resulting exact identities and bounds are

    t*m1=1+Z1+Q1+E1,
    t*m2-m1=Z2+Q2+E2,
    t*m3-m2=Z3+Q3+E3,
    |E_k|<=k*(k+1)/2 *eta^(-k-2)*sum_p alpha_p^3.             (6)

The differentiation has k(k+1)/2 ordered triples in the remainder,
each with coefficient one. Operators depend on T, but the derivatives
are in z at fixed T; no differentiation of the prime cutoff is involved.

## 3. The quadratic term reduces to ordinary primes

Put s_p=log p, c_p=log p/(ell1*sqrt(p)), U_p=K_(s_p), and

    a_p=c_p*(exp(-i*T*s_p)*U_p+exp(i*T*s_p)*U_p^*),
    beta_p=2*c_p,
    e_p=2*log p/(ell1*p*(1-p^(-1/2))).

Then ||a_p||<=beta_p and ||A_p-a_p||<=e_p. The convergent
integer series with decay p^(-3/2) and p^(-2) show that

    sum_p (2*beta_p*e_p+e_p^2)=O(ell^(-2)),
    sum_p c_p^2*alpha_p=O(ell^(-3)).                          (7)

In Q_k replacing both occurrences of A_p by a_p costs at most

    k*eta^(-k-1)*sum_p(2*beta_p*e_p+e_p^2).                  (8)

This replacement is only in the quadratic term. The linear term Z_k
still contains all prime powers; its higher-power contributions have
not been shown to vanish.

Expanding the two factors of a_p yields balanced and unbalanced terms.
Define

    B_k=sum_p c_p^2 sum_(u,v>=1; u+v=k+1)
       Tr(P*(U_p*R^u*U_p^**R^v+U_p^**R^u*U_p*R^v)*P)/d,
    G_k=sum_p c_p^2 sum_(u,v>=1; u+v=k+1)
       Tr(P*(exp(-2*i*T*s_p)*U_p*R_p^u*U_p*R_p^v
             +exp(2*i*T*s_p)*U_p^**R_p^u*U_p^**R_p^v)*P)/d.

In B_k the removed resolvents have also been restored to R. The
bound ||R_p^u-R^u||<=u*alpha_p*eta^(-u-1) follows by telescoping
the powers and using the resolvent identity. Thus restoring them
in the balanced terms costs at most

    2*k*(k+1)*eta^(-k-2)*sum_p c_p^2*alpha_p=O(ell^(-3)).     (9)

Combining (5)--(9) proves

    Q_k=B_k+G_k+O_z(ell^(-2)),
    t*m1=1+Y1+O_z(ell^(-2)),
    t*m2-m1=Y2+O_z(ell^(-2)),
    t*m3-m2=Y3+O_z(ell^(-2)), Y_k=Z_k+B_k+G_k.              (10)

All traces still mean finite-rank compressions. The balanced term B1
can be written Tr(P*S(R)*R*P)/d, where

    S(M)=sum_p c_p^2*(U_p*M*U_p^*+U_p^**M*U_p).

This is an operator covariance map. Its trace has not been factorized,
and no closed scalar or matrix Dyson equation has been established.

There is a geometric norm bound stronger than simply summing all
prime variances. The map S is completely positive, since each summand
has the form U*M*U^*. In particular ||S||=||S(I)||: apply S to the
positive two-by-two operator block with identity diagonal and a
contraction M off the diagonal to obtain ||S(M)||<=||S(I)||; equality
holds for M=I. Direct multiplication
of the gated translations gives

    S(I)=multiplication by
      phi(u)^2 sum_p c_p^2*(phi(u+s_p)^2+phi(u-s_p)^2).       (11)

For |u|<=L/2, surviving shifts have s_p<=L/2-u or s_p<=L/2+u.
The classical [first Mertens estimate](https://terrytao.wordpress.com/2013/12/11/mertens-theorems/),
followed by partial summation, gives
uniformly for Y>=2

    sum_(p<=Y) (log p)^2/p=(log Y)^2/2+O(log(2Y)).

Thus (11) is bounded above by

    [(L/2-u)^2+(L/2+u)^2]/(2*ell1^2)+O(ell^(-1)).

It follows that

    ||S||<=1/2+O(ell^(-1)), Tr(P*S(I)*P)/d ->1/3.            (12)

For the trace limit, Tr(P*M_g*P)/d is exactly the interval average
of g. In scaled coordinates v=u/L, the fixed-width ramp changes that
average by O(1/L); away from the ramp (11) tends to 1/4+v^2.
Its average on [-1/2,1/2] is 1/3. Complete positivity therefore gives
the additional size envelope |B_k|<=k*(1/2+o(1))*eta^(-k-1).
This controls the covariance term's size; it does not evaluate its
signed trace or justify a factorization of correlated resolvents.

## 4. The precise remaining phase and resolvent statistics

Let A1,A2,A3 be the fixed partial-fraction coefficients of the bounded
certificate, not the prime operators A_p. Define

    D1=A1/t+A2/t^2+A3/t^3,
    D2=A2/t+A3/t^2, D3=A3/t.

Solving the triangular relations (10) and using the previous bridge
gives the actual-matrix expansion

    Tr f(B_T)/d=f(1)+2*Re(D1*Y1+D2*Y2+D3*Y3)+o(1).          (13)

The o(1) includes the established explicit-formula and compression
errors; the new prime-removal remainder is O_z(ell^(-2)). Hence
another exact sufficient closing estimate is

    limsup Re sum_(k=1..3) D_k*(Z_k+B_k+G_k)
                              <=(1/10-f(1))/2.              (14)

Z_k consists of first phase harmonics for each prime power against a
resolvent with that whole prime base removed. G_k consists of second
harmonics of ordinary primes against removed resolvents. These
coefficients have no explicit occurrence of their removed base prime
inside R_p, but this is an algebraic observation, not a statistical
independence theorem. The remaining phases are deterministic functions
of the same T. Neither Z_k=o(1) nor G_k=o(1) has been proved here.
Moreover B_k still contains correlated resolvents and has not been
evaluated. Showing phase cancellation alone would not evaluate B_k.
Restoring R_p in the linear term by an absolute resolvent estimate
would cost sum_p alpha_p^2=O(1), rather than o(1). Similarly the
absolute sum of the p^2 coefficients in the linear term has order
one. Neither of these replacements has been used.

Small individual increments do not imply that the conditional linear
term vanishes. On a six-dimensional space take
K0=diag(1,-1,0,0,0,0) and n identical increments A_j=K0/n.
Their sum is K0, its normalized first moment is zero, and its second
moment is 1/3. Nevertheless sum_j ||A_j||^3=1/n^2 while

    sum_j Tr(A_j*(tI-K0+A_j)^(-1))/6
      ->Tr(K0*(tI-K0)^(-1))/6=1/[3*(t^2-1)],

which is nonzero at t=1-i. This is an abstract counterexample to a
small-increment cancellation argument. It is not an actual-prime
example; the needed prime-phase arithmetic has not been supplied by
the norm estimates.

This is the distinction from a probabilistic replacement argument.
The general Lindeberg method compares smooth observables under stated
independence or weak-dependence assumptions; see
[Chatterjee](https://arxiv.org/abs/math/0508519v3).
Those assumptions are not supplied for these actual prime phases.
The expansion (3) is elementary and deterministic; no novelty claim
for the resolvent identity or the replacement method is intended.

## Reproduction

Run `python scripts/verify_prime_removal.py`. It checks the exact
noncommuting finite expansion, differentiated grouped identities,
balanced/unbalanced separation and norm envelopes over rational
matrices. The [record](../results/prime_removal_2026-10-05.json) separates
these finite identities from the asymptotic bounds and the unproved
arithmetic target. Passing finite checks does not evaluate (14).
