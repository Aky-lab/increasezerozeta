# A direct spectral certificate for actual zeta-zero counts

**Status:** proof draft and finite computational checks. The actual arithmetic
trace estimate remains open.

This note isolates a sufficient arithmetic estimate for the project's
conditional counting conversion. The finite counting implication is proved
below, including negative eigenvalues, a nonzero threshold and boundary
zeros. The needed trace estimate for the actual zeta matrix is **not proved**.
The model trace is not substituted for that arithmetic trace.

[Lean kernel proofs](../formal/README.md) cover scalar majorization
and the finite integer count/collar conversion conditional on the
explicit spectral inputs. The inertia argument, actual matrix
identification and arithmetic trace bound remain outside that formal file.

## 1. A finite theorem with no positive-matrix assumption

Let A be a real symmetric d-by-d matrix formed from s1 simple on-line
blocks, s2 multiple on-line blocks, and p off-line pairs. An on-line block
is m*u*u^T, where u is real and its positive integer multiplicity m is
one for a simple block and at least two for a multiple block. An off-line
pair is 2*m*(a*a^T-b*b^T), with real vectors a,b and positive integer m.
No independence of these vectors is assumed.

Let Nprime be the total multiplicity of these blocks. Then

    n_positive(A) <= s1+s2+p,
    rank(A) <= s1+s2+2p,
    Nprime >= s1+2*s2+2*p.

Indeed A=P-Q, where P,Q are positive semidefinite and P is a sum of
s1+s2+p rank-one matrices. On the kernel of P the quadratic form of A
is nonpositive. Its positive subspace therefore has dimension at most
rank(P). The rank bound follows directly from the block expansion.

Write B=A+E, with E symmetric and operator norm at most epsilon>=0.
For any real polynomial q satisfying q(x)>=1 for x<=0, set

    Q(B,epsilon) = Tr(q(B-epsilon*I)^2)/d.

The spectral theorem gives

    #{eigenvalues of B <= epsilon} <= d*Q(B,epsilon),
    n_positive(A) >= #{eigenvalues of B > epsilon}
                  >= d*(1-Q(B,epsilon)).                         (1)

For the second inequality, on a subspace where B-epsilon*I is positive
definite, A=B-E is positive definite because E<=epsilon*I.
The strict threshold in (1) includes equality on its excluded side.
Combining (1) with the multiplicity bound gives

    s1 >= 2*d*(1-Q(B,epsilon))-Nprime,
    s1+s2+2*p >= d*(1-Q(B,epsilon)).                            (2)

Suppose the blocks cover an enlarged height interval Iprime, the target
interval I contains N zeros with multiplicity, and its boundary collar
contains Bdry zeros with multiplicity. Thus Nprime=N+Bdry. Removing the
collar costs at most Bdry simple or distinct zeros. Consequently

    N_simple_on_line(I)/N >= 2*(d/N)*(1-Q)-1-2*Bdry/N,
    N_distinct(I)/N >= (d/N)*(1-Q)-Bdry/N.                       (3)

Negative right-hand sides may be replaced by zero. Equations (1)-(3) are
deterministic; neither a positive limiting measure nor convergence of all
moments is required for them.

## 2. The actual Weil matrix and an independent tail bound

Use the test family and explicit-formula normalization in
[Alpoge--Furman, sections 2 and 4](https://arxiv.org/abs/2608.13637).
Put ell=log(T/(2*pi)), ell1=ell+2*log(2)-1, L=lambda*ell,
X=exp(L), d=floor(L*T/(2*pi)), and tau_k=T+2*pi*k/L for 0<=k<d.
Here lambda is the bandwidth, not a Christoffel bound.
Fix a smooth step function chi taking values in [0,1], equal to zero
on (-infinity,0] and one on [1,infinity), and use the unit-height taper

    phi_T(u)=chi(L/2+u)*chi(L/2-u).

For L>2 this is real, even, smooth, supported in [-L/2,L/2], and equal
to one on [-L/2+1,L/2-1]. Its second derivative has uniformly bounded
L1 norm. The unit-height flat top is part of the arithmetic target:
arbitrary C2 tapers, including zero or constant rescalings, do not have
the same proposed model normalization.

For rho=beta+i*gamma, its spectral argument is

    z_rho=(rho-1/2)/i=gamma-i*(beta-1/2),
    u_rho=(hat(phi)(z_rho-tau_k))_(0<=k<d).

Define the normalized actual matrix

    B_T = (L*ell1)^(-1) sum_rho m_rho*u_rho*u_rho^T.             (4)

The sum runs over distinct nontrivial zeros; m_rho supplies their
multiplicity exactly once.

The transpose in (4) is not a conjugate transpose. Replacing it by a
conjugate transpose would make every summand positive and lose the
off-line information. Reflection rho->1-conjugate(rho) sends z to its
conjugate and u to its conjugate. Each reflected pair therefore contributes
2*m*(Re(u)*Re(u)^T-Im(u)*Im(u)^T), as required in section 1. On-line
vectors are real. Thus B_T is real symmetric without assuming RH.

There is a direct tail estimate avoiding an infinite sampling identity.
Choose D=T^(2/3), I=(T,2*T], Iprime=(T-D,2*T+D], and split (4) as A_T+E_T
by the real ordinate gamma. Two integrations by parts give

    |hat(phi)(z_rho-tau_k)| <= C*X^(1/4)/|gamma-tau_k|^2.

For a tail zero put dist=distance(gamma,[T,2*T])>=D. The norm of its
rank-one contribution is at most C^2*d*X^(1/2)/dist^4. The classical
unit-interval zero-count bound O(log(|gamma|+3)), with multiplicity,
gives sum_tail m_rho/dist^4=O(log(T)/D^3), including negative ordinates.
Summing this absolutely convergent bound and dividing by L*ell1 yields

    ||E_T||_op <= Ctail*T*X^(1/2)/D^3
               = O(T^(-1+lambda/2)),
    Bdry/N = O(T^(-1/3)),        d/N -> lambda.                 (5)

The last ratio uses the Riemann--von Mangoldt formula. These estimates
hold for fixed 0<lambda<=1. Choose epsilon_T to be a nonnegative upper
bound in (5); epsilon_T=T^(-1/4) is admissible for all sufficiently large
T and tends to zero. This weaker collar is adequate for (3).
The proof does not require a trace-norm bound or transfer of the model's
positive definiteness to the actual matrix.

### Removing the vanishing threshold from the arithmetic target

Let q be any fixed real polynomial of degree r>=1 with nonzero leading
coefficient, and put f=q^2. Polynomial growth gives constants depending
only on q such that, for real x and |e|<=1,

    |x|^(2*r) <= C*(1+f(x)),
    |f(x-e)-f(x)| <= C*|e|*(1+|x|^(2*r)).

The first inequality follows from the positive leading coefficient of f
and a bound on a compact interval. For the second, apply the mean-value
theorem to f and bound f' on [x-1,x+1]. Applying both inequalities to
the real eigenvalues of any symmetric B gives

    |Tr(f(B-e*I))/d-Tr(f(B))/d|
      <= C*|e|*(1+Tr(f(B))/d).

Applying the same argument to B-e*I with shift -e gives the reverse
bound in terms of its shifted trace. Thus, if epsilon_T tends to zero,
either finite limsup cap implies the same cap for the other statistic.
No separately assumed convergence of individual moments is necessary.
Consequently the shifted targets below are equivalent to their unshifted
forms whenever the proposed cap is finite. This simplification does not
prove either cap.

The underlying Weil/inertia construction is prior work. This note supplies
an explicit application to the repository's certificate and
states the outstanding arithmetic estimate separately. It claims no
priority for the general counting mechanism.

## 3. The exact degree-four polynomial and the missing estimate

Let D4=162540559 and take

    q4(x)=1-(851656632/D4)*x+(1278096456/D4)*x^2
             -(729534540/D4)*x^3+(140399280/D4)*x^4.

Every nonconstant coefficient of q4(-t) is positive for t>=0. Hence
q4(x)>=1 on x<=0, so section 1 applies even to indefinite B_T.
For the [certified model moments](../scripts/model_moments.py),

    integral q4(x)^2 dnu_model = L4 = 12241115/162540559.

At bandwidth lambda=1 the single **unproved arithmetic target**

    limsup_(T->infinity) Tr(q4(B_T-epsilon_T*I)^2)/d <= L4       (6)

would imply, by (3)-(5),

    liminf N_simple_on_line(T,2*T)/N(T,2*T)
      >= 1-2*L4 = 138058329/162540559 = 0.8493777174...,
    liminf N_distinct(T,2*T)/N(T,2*T) >= 1-L4.                  (7)

For the weaker target of 80%, it suffices to replace L4 on the right of
(6) by 1/10. The permitted excess over the model trace is exactly

    1/10-L4 = 40129409/1625405590 = 0.0246888587....             (8)

For fixed lambda<1, (3) instead gives 2*lambda*(1-Qlimit)-1.
To reach the endpoint by taking lambda->1 after T->infinity, a trace bound
in that same order of limits is needed. Endpoint model moments alone do
not justify moving the two limits.

### A degree-three alternative for the 80% target

The [degree-three model certificate](christoffel_tower.md) is

    q3(x)=1-(8232/2519)*x+(7368/2519)*x^2-(1932/2519)*x^3.

Again q3(x)>=1 on x<=0. Its model squared norm is L3=247/2519.
At unit bandwidth, the single **unproved** estimate

    limsup_(T->infinity) Tr(q3(B_T)^2)/d <= 1/10

therefore suffices for an 80% simple-on-line bound. It involves only
moments through order six. With gamma_j=Tr(B_T^j)/d, gamma_0=1, its
exact unshifted statistic is

    [6345361-41472816*gamma_1+104885808*gamma_2
       -131040168*gamma_3+86095872*gamma_4
       -28469952*gamma_5+3732624*gamma_6]/6345361.

The allowable signed excess above the model is

    1/10-L3 = 49/25190 = 0.00194521635....

Matching L3 would give 2025/2519, approximately 80.3890%, in (3).
This route avoids seventh and eighth moments but has a smaller error
budget than (8); no comparison with actual prime weights is asserted.

### Exact prime-side expression for the target

Let v(t)=(hat(phi)(t-tau_k))_k and K_T(t,s)=v(t)^T*v(s). The explicit
formula gives B_T=(L*ell1)^(-1) integral v(t)*v(t)^T*nu_X(t) dt, where

    nu_X(t)=mu(t)+Pi_X(t)+P_X(t),
    mu(t)=Re[Gamma'/Gamma(1/4+i*t/2)]/(2*pi)-log(pi)/(2*pi),
    Pi_X(t)=Re[X^(1/2+i*t)/(1/2+i*t)]/pi,
    P_X(t)=-(1/pi) sum_(n<=X) Lambda(n)/sqrt(n)*cos(t*log(n)).

For j>=1, multiplication of these finite matrices gives the exact formula

    Tr(B_T^j)=(L*ell1)^(-j) integral_(R^j)
          product_i nu_X(t_i)*K_T(t_i,t_(i+1)) dt_1...dt_j,      (9)

with cyclic indices. These integrals are absolutely convergent: for each
fixed T the finite prime sum is bounded and each matrix-entry integrand
has integrable taper decay against the logarithmic archimedean density.
If q4(x)^2=sum_(j=0..8) a_j*x^j, define

    b_k(epsilon)=sum_(j=k..8) a_j*binomial(j,k)*(-epsilon)^(j-k).

Then the exact scalar to estimate is

    Q_T=b_0(epsilon_T)+sum_(k=1..8)
                      b_k(epsilon_T)*Tr(B_T^k)/d.              (10)

Equations (9)-(10) retain the finite frame kernel and the actual prime
weights. They make no replacement by a sine kernel, presieved primes,
separate autocorrelations or a continuum volume. Bounding their signed
combination is an alternative to evaluating eight limits separately;
it is not automatically easier than doing so.
For q3 the same identities hold with degree six in place of eight.

The [actual prime-matrix derivation](actual_prime_trace.md) gives the
entries explicitly, removes the archimedean and pole terms in a normalized
Schatten norm, and evaluates the multiplicatively balanced contributions.
It reduces the degree-three target to one specified signed off-balance
estimate, which remains unproved.

## 4. Error directions, and why finite checks cannot establish (6)

Write gamma_k=Tr(B_T^k)/d, gamma_0=1. At epsilon=0 the squared-polynomial
coefficients satisfy (-1)^k*a_k>0. The same sign property holds for b_k at
every epsilon>=0. If even gamma_k have supplied upper bounds m_k+e_k and
odd gamma_k have supplied lower bounds m_k-e_k, then

    Q_T <= sum_k b_k(epsilon)*m_k
             +sum_(k=1..8) |b_k(epsilon)|*e_k.                 (11)

Thus one-sided moment estimates suffice. No positive-measure assumption
on the actual spectrum is needed. At epsilon=0 the coefficient sum for
the unknown positive-order moments is

    sum_(k=1..8) |a_k| = 9973263119729203608/26419433320032481.

With all errors bounded by the same e, negligible collar/threshold losses,
and unit bandwidth, (8) is met if

    e <= 6522656571199631/99732631197292036080
       = approximately 0.00006540.

This is a conservative independent-error budget. Cancellation in (10)
could allow a much larger correlated error; it must be proved for the
actual trace rather than inferred from model agreement.

Even matching the first two model moments does not force (6). The positive
probability measure with masses 1/9 at -1/2, 5/9 at 1 and 1/3 at 3/2 has
first moment 1 and second moment 4/3, yet its negative mass 1/9 exceeds L4.
Its q4-square integral must therefore exceed L4. The counterexample is
about the insufficiency of supplied low moments, not about zeta zeros.

## 5. A published large-progression input and its remaining mismatch

[Shao--Teravainen, The Bombieri--Vinogradov theorem for nilsequences,
Theorem 2.9](https://arxiv.org/pdf/2006.05954v2) supplies qualitative
Gowers uniformity for normalized primes in almost all progressions:
the modulus Q is at most x^(1/3-epsilon), **where x is the number of
progression slots**, and the modulus includes the small-prime product.
The theorem has a logarithmically sparse exceptional set. It is relevant
to the [dilated lock criterion](dilated_rectangle_geometry.md), unlike
a theorem controlling only linear prime sums.

For a progression with physical endpoint M and modulus r, its slot count
is approximately H=M/r. Ignoring only fixed small-prime factors at the
exponent level, r<=H^(1/3) means r<=M^(1/4), not r<=M^(1/3).
In the model's two-coefficient coordinates this gives the limiting
tau=1/4 profile; a fixed positive exponent margin gives a smaller range.
The leading model coverage ceiling is exactly

    P_profile(1/4)=1109/144000 = approximately 0.00770139.         (12)

The logarithmic exceptional count can potentially be charged against
log(p)/p weights in a fixed-scale dyadic prime-modulus block. Nevertheless,
the actual varying progression lengths, coupled weights, non-invertible
residues, local factors and signed consumption in (9)-(10) must be mapped
explicitly before invoking it. No full zeta trace estimate is asserted
from that input here. In particular, (12) is a range calculation, not a
proved percentage of the actual remainder. The rest of the model region
cannot be dropped, and a global argument may avoid this regional route.

[Shao's earlier progression theorem](https://arxiv.org/pdf/1607.01814v2)
has a larger range but concerns Mobius/Liouville, not the prime weights
in (9). Its introduction explicitly explains this limitation.

## 6. Why restricted higher correlations do not supply the trace cap

There is another plausible source of a model-to-zero comparison:
[Rudnick--Sarnak's restricted n-level correlation theorem](https://www.math.tau.ac.il/~rudnick/papers/RudnickSarnakCRAS1994.pdf).
Its Fourier support condition is sum_i |xi_i|<2. This is not the full
correlation conjecture. One must also retain its smoothing and complex-zero
conventions before making an unconditional application.

Even the ideal translation-invariant cycle kernel has a support obstacle.
With S_lambda(t)=sin(pi*lambda*t)/(pi*t), a B-cycle test is

    F_B(t_1,...,t_B)=product_i S_lambda(t_i-t_(i+1)).

Write each factor as an integral of exp(2*pi*i*x_i*(t_i-t_(i+1)))
over x_i in [-lambda/2,lambda/2]. Its Fourier support is the image

    xi_i=x_i-x_(i-1),  sum_i xi_i=0,
    max sum_i |xi_i|=2*floor(B/2)*lambda.                       (13)

The maximum follows because cyclic total variation is a convex function
of each coordinate, so a maximum occurs at box vertices. The number of
changes between the two endpoints is at most 2*floor(B/2), attained by
alternating endpoints. Interior points close to these vertices give a
nonzero Fourier density beyond sum_i |xi_i|=2 whenever the maximum exceeds
2; this is a support issue rather than an endpoint of measure zero.

For B=8, the restricted theorem therefore directly covers this entire
test only when lambda<1/4. At such a bandwidth, even the impossible best
case Q=0 in (3) gives 2*lambda-1<0. It yields no positive simple-zero
bound through this certificate. At unit bandwidth the eight-point support
reaches 8, well beyond the restricted theorem's range.
The degree-three alternative has a six-cycle reaching 6*lambda and
requires lambda<1/3 for direct coverage of its entire ideal cycle.
That bandwidth also gives no positive bound through (3).

The single signed trace target does not automatically cancel this
obstacle: its eight-distinct-point cycle has coefficient a_8=q4_4^2>0,
and traces of orders below eight cannot contain eight distinct vertices.
A lower-degree term therefore cannot cancel that component. This does
not rule out an inequality exploiting additional arithmetic or operator
structure. It rules out directly replacing the whole cycle statistic by
its sine-process value using that restricted-support theorem alone.

## Verification and scope

`python scripts/verify_spectral_bridge.py` checks the rational certificate,
shifted coefficient expansion, one-sided error accounting, exact matrix
inertia under rank collapse and bounded perturbations, collar counts,
the low-moment counterexample, progression exponent conversion and cyclic
Fourier support calculation.
The [record](../results/spectral_bridge_2026-10-05.json) saves the exact
polynomial, scalar target and error margin. These finite checks do not
prove the unbounded-height estimate (6), its weaker 80% form, or novelty.
