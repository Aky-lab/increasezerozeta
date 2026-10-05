# A reflected-square certificate with positive resolvent weights

**Status:** exact scalar certificates and conditional counting reductions.
The actual arithmetic trace and resolvent estimates remain unproved.
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
resolvent compression error is O(log(ell)/ell^3), by the block estimate
and the [fixed-taper leakage bound](fixed_taper_projection.md).
The smallest imaginary part exceeds 1/4, so
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

For r>0 define S=(r^2+4/3)/(r+1), P=r^2-r*S and
C=1/(2*r*S^2). The following identity gives a global quadratic
minorant on the entire real line:

    x/(x^2+r^2)-C*(-(r-S)^2+2*S*x-x^2)
       =C*(x^2-S*x+P)^2/(x^2+r^2) >=0.                     (8)

The cleared identity and scalar inequality are proved in Lean for
every positive r,S; no assumption x>=0 is used.

The established first two projected moments of J_T are 1+o(1)
and 4/3+o(1). The second includes the O(log(ell)/ell^3) projection-leakage
correction to (I-H_T)^2. Apply spectral functional calculus to (8):

    liminf Re m_T(i*r) >= h(r),
    h(r)=(1-1/(6*r))/(r^2+4/3).

Here r,S are fixed before taking the limit. The minorant expectation
simplifies using S*(r+1)=r^2+4/3. Thus the moment errors contribute
o(1), with no unknown high moments.

This bound is sharp among probability measures with moments 1,4/3.
The polynomial x^2-S*x+P has value -1/3 at x=1, so its two real
roots x_-<1<x_+ are distinct. Give them weights
(x_+-1)/(x_+-x_-) and (1-x_-)/(x_+-x_-). These positive weights
give mean one and second moment S-P=4/3. The remainder in (8)
vanishes on their support.

There is an independent Cauchy--Schwarz derivation. For a probability
measure of mean a and second moment b, put sigma^2=b-a^2 and
m=integral (x-i*r)^(-1) dmu. Centering both factors gives

    |1+(i*r-a)*m|^2 <= sigma^2*(Im(m)/r-|m|^2).

Completing the square puts m in the disk with center
a/(b+r^2)+i*(r+sigma^2/(2*r))/(b+r^2) and radius
sigma^2/[2*r*(b+r^2)]. Its leftmost real coordinate is
(a-sigma^2/(2*r))/(b+r^2), which recovers h(r) when a=1,b=4/3.

The positive weights let this bound handle r_2 and r_3. Define

    C_large=sum_(j=2,3) w_j*h(r_j),
    theta=(9/10-C_large)/w_1.

The single actual estimate

    liminf Re m_T(i*r_1) >=theta                            (9)

would therefore imply (5). Exact interval computations give
0.98282<theta<0.98283; C_large is approximately 0.552399195.
Condition (9) is stronger than the three-point condition (5).
It is a sufficient research target, not a proved bound.
Each individual two-moment lower bound is sharp; no joint optimality
for their weighted sum is asserted.

The first two moments alone cannot establish (5): the probability
law (1/6)*delta_0+(2/3)*delta_1+(1/6)*delta_2 has those moments,
but its f_m integral is at least 1/6>1/10 by (1). This also rules
out deriving (9) solely from these moment inputs. In (7), the
conditional prime-phase terms and correlated covariance trace still
require arithmetic estimates; removing a prime base does not establish
statistical independence.

## 5. A rational one-point criterion

The integer polynomials Q and D from section 1 satisfy the global bound

    f_m(x) <= A(x)-(41501/100000)*x/(x^2+121/1600),        (10)
    A(x)=1.00912-0.83810*x+0.22949*x^2.

The terminating decimals here denote exact rational numbers.

This bound has a direct exact proof, with no algebraic pole construction.
Both D(x) and x^2+121/1600 are positive. Define the integer polynomial

    Z(x)=160000000*[A(x)*D(x)*(x^2+121/1600)
                 -(41501/100000)*x*D(x)-Q(x)^2*(x^2+121/1600)].

Its ascending coefficients are

    (700223277072, 16130581267790, 38453491895885,
     -657884432134880, 686791679385776, 656143236526080,
     4216355565563136, -7275624514820640, 3177114169954896,
     -500529947904000, 137055981081600).

The exact signed Euclidean remainder sequence for Z,Z' has degrees
10,9,...,0. Positive rational rescaling to primitive integer coefficients
does not change its signs. At minus and plus infinity these are

    -infinity: (+,-,-,+,-,-,+,-,-,-,-),
    +infinity: (+,+,-,-,-,+,+,+,-,+,-).

Both have five sign variations. The final nonzero constant remainder
also proves that Z is squarefree. By
[Sturm's theorem](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Algebra/Polynomial/Sturm/Sequence.html),
Z has no real root. Since Z(0)=700223277072>0, Z is positive everywhere,
which proves (10). The complete integer chain is supplied in the
[exact record](../results/mirror_certificate_2026-10-05.json).
The checker verifies every Euclidean division identity and checks the
root counter on 30 explicit factored examples, including repeated roots.

Apply (10) directly to the finite actual symmetric matrix B_T. Write
tau_d=Tr/d and A_T=I-H_T. The established prime moments and normalized
Schatten-two approximation give

    tau_d(B_T)=1+o(1),   tau_d(B_T^2)=4/3+o(1).

For the second identity, use
|tau_d(B_T^2-A_T^2)|<=||B_T-A_T||_(2,d)*
(||B_T||_(2,d)+||A_T||_(2,d)); the norms remain bounded by the
second prime moment and the approximation. Spectral functional calculus
therefore yields

    tau_d f_m(B_T) <=71551/150000
       -(41501/100000)*Re tau_d(B_T-(11i/40)I)^(-1)+o(1).

The resolvent identity transfers this single rational point:

    |tau_d(B_T-(11i/40)I)^(-1)-tau_d(A_T-(11i/40)I)^(-1)|
        <=(1600/121)*||B_T-A_T||_(2,d)=o(1).

The [block compression bound](bounded_resolvent_bridge.md) and fixed-taper
leakage give tau_d(A_T-(11i/40)I)^(-1)=m_T(11i/40)+o(1). No trace of
an infinite-dimensional identity and no compression assertion for f_m
itself is needed for this deduction. Consequently

    limsup tau_d f_m(B_T)
       <=71551/150000-(41501/100000)*liminf Re m_T(11i/40). (11)

The entirely rational sufficient arithmetic target for 80% is

    liminf Re m_T(11i/40) >=113102/124503.                 (12)

The cleaner, slightly stronger hypothesis liminf Re m_T(11i/40)>=91/100
would give a trace cap of 2980427/30000000 and, by the counting bridge,

    liminf N_simple,on-line(T,2T)/N(T,2T)
        >=12019573/15000000, approximately 0.801304867.    (13)

Here d/N tends to one, the collar loss vanishes, and the globally
bounded derivative of f_m handles the vanishing spectral threshold.
These are conditional implications. Neither (12) nor its stronger
version at 91/100 has been established for the actual prime operator.

The arithmetic content can be stated in the prime-removal variables.
Set t=1-11i/40 and Y_T=Z_1+B_1+G_1 at z=11i/40. The deterministic
identity t*m_T=1+Y_T+O(ell^-2) gives

    Re m_T(11i/40)=1600/1721+Re(Y_T/t)+o(1).

Thus the sufficient bound at 91/100 is equivalent to the closing estimate

    liminf Re[(Z_1+B_1+G_1)(11i/40)/(1-11i/40)]
        >=-3389/172100.                                   (14)

Alternatively, write Xi_T=Tr[P*K_T*(J_T-(11i/40)I)^(-1)*K_T*P]/d.
The exact two-prime insertion identity and the vanishing first prime
moment give m_T=1/t+Xi_T/t^2+o(1). Since
1/t^2=(2366400+1408000*i)/2961841, (14) is also equivalent to

    liminf Re[t^(-2)*Xi_T] >=-3389/172100.

For the positive measure kappa_T(E)=Tr[P*K_T*1_E(J_T)*K_T*P]/d,
this statistic is the integral of
(2366400*x-387200)/[2961841*(x^2+121/1600)].
Its total mass tends to 1/3, which alone does not determine this signed
integral. These identities isolate the phase and covariance estimate
still needed; they do not assume independence of the prime-removal terms.

## 6. Model-moment compatibility

The [continuum Gram law](cue_gram_model.md) supports the target in (12)
without substituting its moments for actual prime moments. Resolvent
disks derived from moment constraints are classical; see di Dio,
[Weyl Circles for one-dimensional Moment Problems](https://arxiv.org/abs/1506.06589).
The argument below specializes a direct Bessel projection to the
certified moment ledger. More generally,
let mu be any probability measure on the real line with supplied moments
m_0,...,m_(2n), and let H=(m_(j+k))_(j,k=0..n) be positive definite.
For z=i*r, r>0, write m(z)=integral (x-z)^(-1) dmu and

    q_j=z^j,  s_0=0,
    s_j=sum_(k=0..j-1) z^k*m_(j-1-k),
    v_j=integral x^j/(x-z) dmu=q_j*m(z)+s_j.

Bessel's inequality for projection onto 1,x,...,x^n gives

    v^*H^(-1)v <=integral |x-z|^(-2) dmu=Im m(z)/r.

Set a=q^*H^(-1)q, b=q^*H^(-1)s and c=s^*H^(-1)s. Then
a>0, c is real, and the inequality is

    a*|m|^2+2*Re(conj(m)*b)+c <=Im(m)/r.

Completing the square gives a disk centered at
(-Re(b)/a,(1/(2*r)-Im(b))/a). Its radius squared is
[Re(b)^2+(Im(b)-1/(2*r))^2-a*c]/a^2.
All these quantities are exact rational numbers at a rational imaginary
pole with rational moments.

For the certified model moments through eighth order, H is positive
definite: its successive LDL pivots are

    (1, 1/3, 5/36, 247/5040, 2448223/104569920).

At r=11/40, n=4, the disk therefore proves

    Re m_model(11i/40)
       >=23741200782961777600/26054764615856959631
       >0.9112>91/100.                                    (15)

This bound needs no assumption of nonnegative spectral support. It
applies to every probability measure on the real line with those
eight moments. For the actual matrix only the first two moment limits
are established here; (15) supplies model evidence, not the missing
arithmetic estimate.

Support restrictions can also reject unsuitable candidate bounds.
For a measure supported on [0,infinity), use the weighted inner product
integral x*conj(f)*g dmu and the shifted matrix (m_(j+k+1)). With
v_j=integral x^(j+1)/(x-z) dmu, Bessel's inequality becomes
v^*H_shift^(-1)v<=Re m(z). At z=i/4 and degree three, exact completion
of the square with the model moments through seventh order gives

    Re m_model(i/4)<=742724016/753835003<1.                (16)

Thus a sufficient condition requiring a value of one at i/4 would
be incompatible with convergence to this positive model.
The actual matrix may have negative spectral support; (16) does not
give an upper bound for its resolvent from the first two moments alone.

An exact scalar stability estimate at i/4 is

    24255*(16*x^2+1)^2-295680*x
      =55*(21-128*x-96*x^2)^2
       +32*(55*x-384*x^2)^2+983808*x^4 >=0.

Apply the identity at x and -x and divide by the positive denominator
to obtain |x|/(x^2+1/16)^2<=21. Thus for v>=1/16,

    |x/(x^2+v)-x/(x^2+1/16)|<=21*(v-1/16).

Lean checks the cleared identity, both signs of its inequality and
the exact rational constants and comparisons in (11)--(16). The Sturm root count and
the analytic operator transfers remain outside that formalization.

## Reproduction and scope

Run `python scripts/verify_mirror_certificate.py`. The
[exact record](../results/mirror_certificate_2026-10-05.json) contains
the square coefficients, endpoint signs, refined rational pole and
weight intervals, the isolated one-point target and the rational
majorant's complete Sturm certificate, model moment matrices and
Bessel-disk coefficients. A separate
companion-matrix determinant/adjugate identity verifies the complete
partial fractions by coefficient comparison. No sampled numerical
evaluation is used as a polynomial proof.

The [Lean file](../formal/CountingBridge.lean) proves denominator
positivity, the bounds in (1), denominator comparison, threshold
majorization, the rational-pole sum of squares, exact conditional
constants and the conditional finite counting conversion. Pole
construction, functional calculus, compression and arithmetic cap
are outside that formalization. The symmetrization and partial-fraction
steps are elementary; no general novelty claim for them is made.
