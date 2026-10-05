# A reflected-square certificate with positive resolvent weights

**Status:** scalar identities and conditional counting reduction. The
actual arithmetic inequalities (5), (7) and (9) below remain unproved.
The scalar majorization and finite counting implication are
[formalized in Lean](../formal/README.md).

Use the actual matrix B_T and unit-bandwidth prime operator J_T=I-K_T,
with Fourier projection P of rank d, as defined in the
[bounded bridge](bounded_resolvent_bridge.md). All normalized traces
below are finite projected traces. No trace over the full infinite
Fourier space is used.

## 1. Scalar majorization and comparison

Set

    Q(x)=2519-8232*x+7368*x^2-1932*x^3,
    D(x)=6345361+104885808*x^2+86095872*x^4+3732624*x^6,
    f_m(x)=Q(x)^2/D(x).

Direct expansion gives

    2*D(x)=Q(x)^2+Q(-x)^2,
    Q(x)^2=D(x)-x*(41472816+131040168*x^2+28469952*x^4).

Since D(x)>=6345361 for real x, these identities prove

    0<=f_m(x)<=2,  f_m(-x)=2-f_m(x),
    f_m(x)>=1 for x<=0.                                    (1)

The function is smooth, globally Lipschitz and tends to one at both
infinities: its denominator has no real root and its derivative is a
continuous rational function tending to zero. These facts allow the
same threshold shift and normalized Schatten-two transfer used in the
bounded bridge. In particular, the condition

    limsup Tr f_m(B_T)/d <=1/10                             (2)

is sufficient for the conditional 80% simple-critical-line conclusion.
The finite implication retains all inertia, multiplicity, collar and
dimension hypotheses; none is replaced by a model moment.

For 0<=alpha<=1 with 6345361*alpha^3<=3732624,

    D(x)-6345361*(1+alpha*x^2)^3
      =(104885808-19036083*alpha)*x^2
       +(86095872-19036083*alpha^2)*x^4
       +(3732624-6345361*alpha^3)*x^6 >=0.

Thus f_m is pointwise at most the bounded certificate f_alpha, and
both are at most (Q/2519)^2. The sufficient cap (2) is consequently
weaker than the f_alpha cap. The supplied continuum model law has

    integral f_m dnu <=247/2519 <1/10.

This is an upper bound for a model integral, not an arithmetic estimate.

## 2. Three simple pole pairs and positive weights

Let v_1<v_2<v_3 be the roots of

    P(v)=6345361-104885808*v+86095872*v^2-3732624*v^3.

Exact endpoint evaluations have signs +,-,-,+,+,- at
1/16,1/15,1,2,20,24 respectively. The intermediate value theorem
therefore puts the roots in the three disjoint intervals

    v_1 in (1/16,1/15), v_2 in (1,2), v_3 in (20,24).

A cubic with three distinct roots has no additional or repeated roots.
Put r_j=sqrt(v_j)>0 and

    w_j=(41472816-131040168*v_j+28469952*v_j^2)
          /(104885808-172191744*v_j+11197872*v_j^2).

The numerator factors as
2*(2519-7368*v_j)*(8232-1932*v_j). Its signs on the three
isolating intervals are +,-,+. The derivative denominator has the
same signs, either from the ordered simple roots or direct interval
bounds. Hence every w_j is positive. In particular no pole cancels.

Polynomial interpolation at the three roots of
d(u)=6345361+104885808*u+86095872*u^2+3732624*u^3 gives

    f_m(x)=1-sum_(j=1..3) w_j*x/(x^2+v_j)
          =1-sum_(j=1..3) w_j*Re (x-i*r_j)^(-1).             (3)

Indeed n(u)=41472816+131040168*u+28469952*u^2 has
n(u)/d(u)=sum_j [n(-v_j)/d'(-v_j)]/(u+v_j); multiplying
through by d(u) yields degree-two polynomials agreeing at three
distinct points. The residues at both poles +/-i*r_j are -w_j/2.
For orientation, rounded values are

    r=(0.252652463,1.105801246,4.666813380),
    w=(0.353674096,0.856816637,6.416838460).

The exact definitions and rational isolating intervals, rather than
these decimals, specify the certificate.

## 3. Actual resolvent and prime-removal targets

Define m_T(z)=Tr[P*(J_T-zI)^(-1)*P]/d. At each fixed pole the
resolvent compression error is O(ell^-2), by the block estimate in
the bounded bridge. The smallest imaginary part exceeds 1/4, so
the constants are finite and independent of T. Formula (3) and
the established explicit-formula transfer give

    Tr f_m(B_T)/d=1-sum_j w_j*Re m_T(i*r_j)+o(1).            (4)

The precise sufficient arithmetic target is

    liminf sum_j w_j*Re m_T(i*r_j) >=9/10.                  (5)

Only first resolvents at three points occur. There are no derivatives
and no complex weights.

In the [prime-removal notation](prime_removal_resolvent.md), put
t_j=1-i*r_j and Y_T(z)=Z_1(z)+B_1(z)+G_1(z). Its deterministic
expansion gives t_j*m_T(i*r_j)=1+Y_T(i*r_j)+O(ell^-2).
Since f_m(1)=76729/201059665, (4) is equivalently

    Tr f_m(B_T)/d=f_m(1)-sum_j w_j*Re[Y_T(i*r_j)/t_j]+o(1). (6)

Thus a sufficient closing estimate is

    liminf sum_j w_j*Re[(Z_1+B_1+G_1)(i*r_j)/t_j]
                              >=-8011695/80423866.         (7)

The exact two-prime insertion identity also gives (6) with
Y_T(i*r_j)/t_j replaced by Xi_1,T(i*r_j)/t_j^2; its linear
first moment is o(1). This is another equivalent closing statistic.

## 4. A sufficient estimate at one point

For real x and r>0,

    x/(x^2+r^2)-x/r^2+x^2/(2*r^3)
       =x^2*(x-r)^2/[2*r^3*(x^2+r^2)] >=0.                 (8)

The established first two projected moments of J_T are 1+o(1)
and 4/3+o(1). The second includes the O(ell^-2) projection-leakage
correction to (I-H_T)^2. Apply spectral functional calculus to (8):

    liminf Re m_T(i*r) >= r^(-2)-2/(3*r^3).

The positive weights let this bound handle r_2 and r_3. Define

    C_large=sum_(j=2,3) w_j*(r_j^(-2)-2/(3*r_j^3)),
    theta=(9/10-C_large)/w_1.

The single actual estimate

    liminf Re m_T(i*r_1) >=theta                            (9)

would therefore imply (5). Exact interval computations give
1.04388<theta<1.04389; C_large is approximately 0.530805482.
Condition (9) is stronger than the three-point condition (5).
It is a sufficient research target, not a proved bound.

The first two moments alone cannot establish (5): the probability
law (1/6)*delta_0+(2/3)*delta_1+(1/6)*delta_2 has those moments,
but its f_m integral is at least 1/6>1/10 by (1). This also rules
out deriving (9) solely from these moment inputs. In (7), the
conditional prime-phase terms and correlated covariance trace still
require arithmetic estimates; removing a prime base does not establish
statistical independence.

## Reproduction and scope

Run `python scripts/verify_mirror_certificate.py`. The
[exact record](../results/mirror_certificate_2026-10-05.json) contains
the square coefficients, endpoint signs, refined rational pole and
weight intervals, and the isolated one-point target. A separate
companion-matrix determinant/adjugate identity verifies the complete
partial fractions by coefficient comparison. No sampled numerical
evaluation is used as a polynomial proof.

The [Lean file](../formal/CountingBridge.lean) proves denominator
positivity, the bounds in (1), denominator comparison, threshold
majorization and the conditional finite counting conversion. Pole
construction, functional calculus, compression and arithmetic cap
are outside that formalization. The symmetrization and partial-fraction
steps are elementary; no general novelty claim for them is made.
