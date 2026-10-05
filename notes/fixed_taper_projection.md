# Fixed-taper projection bounds and the high-frequency gauge

**Status:** analytic lemmas with exact finite calibrations. These bounds
control projection errors and describe the prime operator; they do not
establish the arithmetic resolvent cap. Independent mathematical review
of the continuous arguments remains necessary.

Use the unit-height taper, gated translations K_s, prime operator K_T,
and Fourier projection P from the [bounded bridge](bounded_resolvent_bridge.md).
Write Q=I-P. The normalization is ell=log(T/(2*pi)),
L=lambda*ell, X=exp(L), d=floor(L*T/(2*pi)). In particular d/N tends
to lambda, where N counts zeros in (T,2T] with multiplicity; the
logarithm is taken after division by 2*pi.

## 1. A boundary estimate independent of the number of modes

For a fixed C2 step chi, the norms ||phi'||_1, ||phi'||_infinity
and ||phi''||_1 are bounded independently of L>2. Set

    C1=2*||phi'||_1,
    C2=2*||phi''||_1+2*||phi'||_infinity*||phi'||_1.

The multiplier psi_s(u)=phi(u)*phi(u+s), extended by zero, has
||psi_s'||_1<=C1 and ||psi_s''||_1<=C2, uniformly in s. Its
normalized Fourier coefficient satisfies, for k!=0,

    |c_k(s)| <= min(C1/(2*pi*|k|), C2*L/(4*pi^2*k^2)).       (1)

Both integrations by parts have zero boundary terms. The coefficients
of K_s in the Fourier basis are c_(a-b)(s) times unit phases.
For a fixed mode difference k, exactly min(d,|k|) modes cross P's
edge in each indicated direction. Parseval therefore gives

    ||Q*K_s*P||_HS^2
      =sum_(k!=0) min(d,|k|)*|c_k(s)|^2.

Let M=ceil(L) and H_M=sum_(k=1..M)1/k. For k<=M use the first
bound in (1), and for k>M use the second. Since
sum_(k>M)k^(-3)<=1/(2*M^2),

    ||Q*K_s*P||_HS^2
       <= C1^2/(2*pi^2)*H_M+C2^2*L^2/(16*pi^4*M^2)
       <= C1^2/(2*pi^2)*H_M+C2^2/(16*pi^4)
       =O_chi(log(2+L)).                                   (2)

The reverse leakage obeys the same bound. The estimate is uniform
in s and d. Summing the prime operator by Hilbert--Schmidt triangle
inequality and using sum_(n<=X)Lambda(n)/sqrt(n)=O(sqrt(X)) gives

    ||Q*K_T*P||_HS^2/d
       =O_chi(X*log(2+L)/(d*ell1^2)).                       (3)

At unit bandwidth, (3) is O(log(ell)/ell^3). The block-resolvent
identity consequently gives, for fixed nonreal z and integer k>=1,

    |Tr[P*(J_T-z)^(-k)*P-(P*J_T*P-z)^(-k)]|/d
       <=k*(k+1)/2 * ||Q*J_T*P||_HS^2/(d*|Im z|^(k+2)).    (4)

Here the second inverse acts on P's range. Equation (4) improves the
compression error only; the explicit-formula transfer and prime-removal
remainders have their own bounds.

For products of a fixed number j of gated translations, the commutator
identity [P,A*B]=[P,A]*B+A*[P,B] bounds the leakage of a product by
the sum of the individual commutator Hilbert--Schmidt norms. Each
factor is a contraction. The two-boundary trace argument in the
[second-moment note](prime_second_moment.md) then gives an
O_(j,chi)(log(2+L)) unnormalized finite-frame error. Summing all prime
words absolutely yields only

    O_(j,chi)(X^(j/2)*log(2+L)/(d*ell1^j))
      =O_(j,chi)(T^(j/2-1)*log(ell)/ell^(j+1))              (5)

at unit bandwidth. For j>=3 this envelope still fails to tend to zero.
The improved boundary estimate does not resolve higher off-balance
moments or the [positive-weight resolvent target](mirror_resolvent_certificate.md).

## 2. Exact constants and a Fourier kernel

A concrete C2 ramp is

    chi(v)=0 for v<=0,
           10*v^3-15*v^4+6*v^5 for 0<=v<=1,
           1 for v>=1.

Its first two derivatives match at both endpoints. For L>2 the taper
is chi(L/2-|u|), constant near zero. This example has

    ||phi'||_1=2, ||phi'||_infinity=15/8, ||phi''||_1=15/2,
    C1=4, C2=45/2.

Indeed chi'=30*v^2*(1-v)^2 is nonnegative and has maximum 15/8
at v=1/2; chi''=60*v*(1-v)*(1-2*v) changes sign only there.
Thus integral_0^1|chi''|=2*chi'(1/2)=15/4. With two disjoint
transition edges these values give the stated norms. Formula (2) is
then the explicit bound

    ||Q*K_s*P||_HS^2 <=8/pi^2*H_(ceil L)+2025/(64*pi^4).     (6)

The quintic example is C2, not C-infinity. The general estimate (2)
also applies to the fixed smooth step used in the counting bridge.
No particular derivative constants are asserted for an unspecified chi.

Centering the defining integral for the actual F matrix gives exactly

    F_ab(s)/L=c_(a-b)(s)*exp(-i*T*s-i*pi*(a+b)*s/L),
    c_k(s)=1/L integral phi(v+s/2)*phi(v-s/2)
                                  *cos(2*pi*k*v/L) dv.      (7)

The coefficients are real and even in k and s. Hence

    (H_T)_ab=(2/ell1) sum_(n<=X) Lambda(n)/sqrt(n)
       *c_(a-b)(log n)*cos((T+pi*(a+b)/L)*log n).            (8)

This real symmetric matrix formula retains the exact finite kernel.
For the quintic ramp, when 1<=|s|<=L-2, the overlap in (7) is the
convolution of an interval of length L-|s|-1 with the centered density
30*(1/4-v^2)^2 on [-1/2,1/2]. Consequently, for k!=0,

    c_k(s)=Psi(2*pi*k/L)
             *sin(pi*k*(1-(|s|+1)/L))/(pi*k),
    c_0(s)=1-(|s|+1)/L,                                   (9)
    Psi(w)=30 integral_(-1/2)^(1/2)(1/4-v^2)^2*cos(w*v) dv.

Repeated integration by parts evaluates the last integral:

    Psi(w)=120*((12-w^2)*sin(w/2)-6*w*cos(w/2))/w^5
                                                    for w!=0,
    Psi(0)=1.                                             (10)

For small w the entire series avoids cancellation in (10):

    Psi(w)=sum_(j>=0) (-1)^j*w^(2*j)/(2*j)!
                *15/[2^(2*j)*(2*j+1)*(2*j+3)*(2*j+5)].     (11)

The integral form gives the Taylor remainder bound after degree 2*J:

    |remainder| <= |w|^(2*J+2)/(2*J+2)!
        *15/[2^(2*J+2)*(2*J+3)*(2*J+5)*(2*J+7)].

Outside the range in (9), the overlap remains piecewise polynomial
of degree at most ten and its Fourier integral can be computed piece
by piece. These formulas supply finite-matrix building blocks, without
asserting convergence to a model law or a computed arithmetic cap.

## 3. The phases move into the projection

On the fundamental interval define U_Tg(u)=exp(i*T*u)*g(u).
This is a measurable unitary on the circle's L2 space even when
T*L/(2*pi) is not an integer. Every surviving K_s path stays inside
the interval, so

    U_T*K_s*U_T^*=exp(-i*T*s)*K_s,
    K_T=U_T*K_0*U_T^*,
    K_0=ell1^(-1) sum_(n<=X) Lambda(n)/sqrt(n)
                                     *(K_(log n)+K_(-log n)).

The transformed projection P_T=U_T^*P*U_T spans
L^(-1/2)*exp(-i*(T+2*pi*a/L)*u), 0<=a<d. Thus the precise
resolvent statistic becomes

    m_T(z)=Tr[P_T*(I-K_0-zI)^(-1)*P_T]/d.                  (12)

The arithmetic remains in this high-frequency projection. Treating
the whole-space operator alone as the target loses that dependence.
Allowing surviving wrap-around gates would also invalidate this
gauge at noncommensurate T; zero extension is essential.

The phase-free operator preserves nonnegative functions but need not
be positive semidefinite. Its norm in fact diverges as L=lambda*ell
grows with fixed lambda>0. Let g be the constant unit vector L^(-1/2).
Then

    <g,K_0*g>=(2/ell1) sum_(n<=X) Lambda(n)/sqrt(n)*c_0(log n).

For any unit-height taper, if L>5 and s in [L-4,L-3], the common
flat core has length at least one, so c_0(s)>=1/L. Put
alpha=exp(-4), beta=exp(-3) and restrict the sum to primes
alpha*X<=p<=beta*X. The
[prime number theorem](https://dlmf.nist.gov/27.12#E4) gives
pi(beta*X)-pi(alpha*X)~(beta-alpha)*X/log X. Since log p
is asymptotic to log X on this interval and p<=beta*X,

    sum_(alpha*X<=p<=beta*X) log p/sqrt(p) >=c*sqrt(X)

for all sufficiently large X, with some fixed c>0. Therefore

    sup spectrum(K_0) >=<g,K_0*g>
                       >=c'*sqrt(X)/(ell1*L) ->infinity.   (13)

A norm-convergent Neumann expansion at t=1-i*r_1, or a positive
Laplace integral requiring Re t>sup spectrum(K_0), is unavailable.
This is a structural obstruction to that shortcut, not an obstruction
to proving the projected estimate (12) by other means.

Kernel positivity and the first two moments alone are insufficient:
a nonnegative matrix with one 2-by-2 adjacency block and four zero
directions has I-K eigenvalues 0,2,1,1,1,1. Its normalized moments
are 1 and 4/3, but its reflected-square certificate trace is at least
1/6>1/10. This is an abstract example, not an actual prime operator.

## Reproduction and proof scope

`python scripts/verify_prime_trace.py` checks the quintic polynomial,
endpoint derivative matching, exact integrals and derivative norms.
It independently compares thirteen centered-density moments with the
series coefficients of (10), checks rational two-decay envelopes,
and checks the noncommensurate gauge on zero-extended finite grids.
A surviving-wrap counterexample distinguishes the required gates.
Twenty-four direct numerical overlap integrals independently calibrate
the Fourier kernel (9); they are numerical checks, not analytic proofs.
The [record](../results/prime_trace_2026-10-05.json) identifies these
calibrations separately from the continuous estimates (2)-(5) and
the PNT deduction (13). They are not Lean proofs of those estimates.
