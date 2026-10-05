# A bounded certificate and a single-point resolvent target

**Status:** analytic reduction in proof-draft form. The arithmetic
resolvent estimate in section 5 is unproved. No new zeta-zero bound
follows from the finite checks or from this reduction alone.

Retain the fixed flat-top taper, normalization and actual matrix B_T
from the [counting bridge](spectral_counting_bridge.md), and the prime
matrix H_T from the [explicit-formula derivation](actual_prime_trace.md).
The discussion below fixes unit bandwidth L=ell, X=exp(L),
d=floor(L*T/(2*pi)), ell1=ell+2*log(2)-1.

## 1. A bounded majorant

Let

    q(x)=1-(8232/2519)*x+(7368/2519)*x^2-(1932/2519)*x^3,
    a=1932/2519, r=a^(-1/3), alpha=r^(-2),
    f(x)=q(x)^2/(1+alpha*x^2)^3.                              (1)

The positive real cube root is intended, so r>1 and 0<alpha<1.
The function f is nonnegative on the real line. For t>=0, the
coefficients of q(-t)^2, in ascending order, are

    (6345361,41472816,104885808,131040168,
     86095872,28469952,3732624)/6345361.

Subtract (1+alpha*t^2)^3. The constant coefficient vanishes; the
sixth coefficient vanishes because alpha^3=a^2. The second and
fourth coefficients remain positive: 104885808/6345361>3>3*alpha
and 86095872/6345361>3>3*alpha^2. The odd coefficients are positive.
Thus

    f(x)>=1 for x<=0.                                        (2)

The denominator has no real zeros. Both leading coefficients in (1)
are a^2, hence f(x)->1 at both infinities. Therefore f is bounded,
smooth and globally Lipschitz. Moreover

    f(x)<=q(x)^2/(1+a^2*x^6)<=q(x)^2.                         (3)

The first inequality uses (1+alpha*x^2)^3>=1+alpha^3*x^6.
For the continuum model law nu, the already computed polynomial
certificate gives

    integral f dnu <=247/2519 <1/10.                          (4)

This is an upper bound, not an evaluation of the rational model integral.

## 2. Counting and archimedean transfer

For any real symmetric d-by-d matrix B and epsilon>=0, (2) gives

    #{eigenvalues of B <=epsilon} <=Tr f(B-epsilon*I).

If M=sup |f'|, the normalized traces with and without this shift differ
by at most M*epsilon. The [collar and inertia argument](spectral_counting_bridge.md)
therefore applies with Q=Tr f(B_T)/d+o(1). In particular,

    limsup Tr f(B_T)/d <=1/10                                (5)

suffices for the same conditional 80% simple-critical-line counting
conclusion. No sixth-moment bound is needed for the shift.

The explicit-formula derivation proves
||B_T-(I-H_T)||_(2,d)->0. Hoffman--Wielandt and Cauchy--Schwarz imply

    |Tr f(B_T)/d-Tr f(I-H_T)/d|
      <=M*||B_T-(I-H_T)||_(2,d) ->0.                         (6)

The eigenvalues are paired in increasing order. The Hermitian
Frobenius perturbation inequality is standard; see
[Li--Mathias, equation (1.5)](https://cklixx.people.wm.edu/lmw.pdf).
Unlike polynomial transfer, (6) requires no high-moment coercivity.

This is a weaker sufficient condition than the polynomial cap: (3)
shows that the latter implies it. Conversely, a vanishing fraction of
arbitrarily large eigenvalues has vanishing cost for bounded f, but can
make normalized sixth moments diverge. For example take d=n^4, d-1
eigenvalues equal to 1 and one equal to n. Its first two moments tend
to 1, its sixth moment diverges, and its f-trace tends to
f(1)<=q(1)^2=76729/6345361<1/10. This is an abstract matrix example,
not a claim about the zeta matrix or its actual second moment.

## 3. Whole-circle prime operator and leakage

On the circle of length L use the orthonormal Fourier basis
e_k(u)=L^(-1/2)*exp(-2*pi*i*k*u/L), k in Z. Represent a point u by
[-L/2,L/2), extend phi by zero on the real line, and define

    (K_s g)(u)=phi(u)*phi(u+s)*g(u+s).

The function g is periodic, but a surviving term requires both real
points u and u+s in the fundamental interval. There is no surviving
wrap-around path. Thus K_s^*=K_(-s) and ||K_s||_op<=1.
Direct integration gives its matrix entries as exp(i*T*s)*F_ab(s)/L.
Let P project onto Fourier modes 0,...,d-1 and Q=I-P. Set

    K_T=ell1^(-1) sum_(n<=X) Lambda(n)/sqrt(n)
          *(exp(-i*T*log n)*K_(log n)
             +exp(i*T*log n)*K_(-log n)),
    J_T=I-K_T, A_T=P*J_T*P|_(P H)=I-H_T.                     (7)

For each T these are bounded self-adjoint operators; their norms
need not be bounded uniformly in T.

The function phi(u)*phi(u+s) has uniformly bounded variation,
independent of s and L. Its nonzero Fourier coefficients satisfy
|c_m|<=C/|m|. For a shift m the number of modes crossing P's edge is
min(d,|m|). Hence

    ||Q*K_s*P||_HS^2
       <=C sum_(m!=0) min(d,|m|)/m^2 <=C'*log(2d).

The same bound holds for the reverse leakage. Partial summation of
Chebyshev's estimate gives sum_(n<=X) Lambda(n)/sqrt(n)=O(sqrt(X)).
Triangle inequality in Hilbert--Schmidt norm yields, with
V=Q*J_T*P=-Q*K_T*P,

    ||V||_HS^2/d = O(X*log(2d)/(d*ell1^2))=O(ell^(-2)).       (8)

This estimate uses no arithmetic cancellation.

## 4. Resolvent powers remove compression

Suppress T and put D=Q*J*Q|_(Q H). For nonreal z set
R_J=(J-zI)^(-1), R_A=(A-zI)^(-1), R_D=(D-zI)^(-1).
All three operator norms are at most |Im z|^(-1), independently of
the norms of J, A and D. Block multiplication gives the identity

    P*R_J*P-R_A = R_A*V^**R_D*V*P*R_J*P.                    (9)

Differentiate k-1 times with respect to z and divide by (k-1)!.
Since the jth derivative of a resolvent is j! times its (j+1)st
power, the result is

    P*R_J^k*P-R_A^k
      =sum_(b+c+e=k-1) R_A^(b+1)*V^**R_D^(c+1)
                                 *V*P*R_J^(e+1)*P.          (10)

All indices are nonnegative; there are k*(k+1)/2 summands, each
with coefficient one. The trace ideal inequality then gives

    |Tr(P*R_J^k*P-R_A^k)|/d
       <=k*(k+1)/2 * ||V||_HS^2/(d*|Im z|^(k+2)).            (11)

The traces are on the finite-dimensional space P H. Boundedness
for each T avoids unbounded-operator domain issues even though Q H
is infinite-dimensional. Equations (9)--(11) hold for every fixed
positive integer k. For the certificate only k=1,2,3 are needed.

## 5. One complex point and the unresolved arithmetic estimate

Put z=i*r, N(x)=q(x)^2. Define

    A3=r^6*N(z)/(2*z)^3,
    A2=r^6*(N'(z)/(2*z)^3-3*N(z)/(2*z)^4),
    A1=r^6*(N''(z)/(2*(2*z)^3)
            -3*N'(z)/(2*z)^4+6*N(z)/(2*z)^5).                (12)

These are the principal-part coefficients at the triple pole z.
Expanding r^6*N(x)/(x+z)^3 at x=z proves, for real x,

    f(x)=1+2*Re sum_(k=1..3) A_k/(x-z)^k.                   (13)

The pole at -z has conjugate coefficients, and the polynomial part
is 1. Functional calculus, (8), (11) and (6) yield

    Tr f(B_T)/d
      =Tr(P*f(J_T)*P)/d+o(1)
      =1+2*Re sum_(k=1..3) A_k*m_(k,T)(i*r)+o(1),            (14)
    m_(k,T)(z)=Tr(P*(J_T-zI)^(-k)*P)/d.

The compression error is O(ell^(-2)); the other o(1) terms come
from the explicit-formula reduction and counting collar. Write the
compressed trace Tr(P*f(J_T)*P), never a full trace Tr f(J_T): the
full space is infinite-dimensional and f tends to 1.

Thus a sufficient arithmetic target is exactly

    limsup_(T->infinity)
      Re sum_(k=1..3) A_k*m_(k,T)(i*r) <=-9/20.              (15)

Only this real combination is required. Separate convergence of
the three resolvent powers is a stronger possible approach.

Neither the known first two moments nor the balanced prime moments
prove (15). The measure (1/4)*delta_0+(3/4)*delta_(4/3) has first
moment 1 and second moment 4/3, yet its f-integral is at least 1/4.
This demonstrates insufficient information, not a possible actual
zeta law. Nor may a Neumann series in K_T be used at the fixed pole
without convergence: the available bound ||K_T||=O(sqrt(X)/ell)
grows. A formal series would conceal uncontrolled prime words of
all lengths. Replacing the actual resolvent by a CUE or sine-model
resolvent is precisely an additional arithmetic assertion.

The reduction avoids high-moment outlier control and the absolute
higher-word boundary cost described in the
[second-moment note](prime_second_moment.md). It leaves the signed
actual-prime resolvent bound (15) as a separate, unresolved problem.

## 6. An exact two-prime sum with a resolvent insertion

There is a finite resolvent identity that uses no Neumann series. Put
t=1-z and K=K_T. Since J-zI=tI-K,

    (J-zI)^(-1)=t^(-1)*I+t^(-2)*K
                          +t^(-2)*K*(J-zI)^(-1)*K.          (16)

Define

    Xi_(k,T)(z)=Tr(P*K*(J-zI)^(-k)*K*P)/d.

Let M1=Tr(P*K*P)/d=Tr H_T/d. Equation (16) and its first two
derivatives give exact identities

    m1=t^(-1)+t^(-2)*M1+t^(-2)*Xi1,
    m2=t^(-2)+2*t^(-3)*M1+2*t^(-3)*Xi1+t^(-2)*Xi2,
    m3=t^(-3)+3*t^(-4)*M1+3*t^(-4)*Xi1
                              +2*t^(-3)*Xi2+t^(-2)*Xi3.    (17)

The [second-moment deduction](prime_second_moment.md) supplies
M1=o(1), and

    ||K*P||_HS^2/d=Tr H_T^2/d+||Q*K*P||_HS^2/d ->1/3.

At the same fixed z=i*r, set

    C1=A1/t^2+2*A2/t^3+3*A3/t^4,
    C2=A2/t^2+2*A3/t^3, C3=A3/t^2.

Equations (13), (14) and (17) imply

    Tr f(B_T)/d
      =f(1)+2*Re(C1*Xi1+C2*Xi2+C3*Xi3)+o(1).               (18)

Consequently the exact alternative closing lemma is

    limsup Re sum_(k=1..3) C_k*Xi_(k,T)(i*r)
                         <=(1/10-f(1))/2.                   (19)

Its arithmetic content is explicit. Write w_n=Lambda(n)/sqrt(n),
s_(n,sigma)=sigma*log(n), and U_(n,sigma)=K_(s_(n,sigma)). Then

    Xi_(k,T)(z)=1/(d*ell1^2)
      *sum_(n,m<=X; sigma,rho in {-1,+1}) w_n*w_m
         *exp(-i*T*(sigma*log(n)+rho*log(m)))
         *Tr(P*U_(n,sigma)*(J_T-zI)^(-k)*U_(m,rho)*P).        (20)

For each T this is a finite two-prime-power sum, not an infinite
word expansion. The inserted resolvent still depends on the entire
prime sum K_T. Standard separable sampled mean-value estimates do
not by themselves control these correlated operator coefficients.
The elementary bound |Xi_k|<=(1/3+o(1))*r^(-k) gives no required
signed cancellation. Indeed define the positive measure

    kappa_T(E)=Tr(P*K_T*1_E(J_T)*K_T*P)/d.

Then Xi_k=integral (x-z)^(-k) d kappa_T(x), and kappa_T has total
mass tending to 1/3. Positivity and this total mass constrain the size
of the remainder, not the required spectral distribution or its sign.
Equations (19)--(20) isolate the needed additional arithmetic estimate;
they do not prove it.

Even the zero-block trace constraints do not repair the information
gap. For d divisible by six, use 2d/3 orthogonal simple blocks
e_j*e_j^T, d/6 orthogonal double blocks 2*e_j*e_j^T, and d/6 unused
directions. All underlying vectors have squared norm one; the total
multiplicity is N=d and the simple fraction is 2/3. The spectral law
is (1/6)*delta_0+(2/3)*delta_1+(1/6)*delta_2. Its first two moments
are 1 and 4/3, but its f-integral is at least 1/6>1/10. This is an
abstract information obstruction, not a zeta example. Additional
arithmetic input is essential even when the block constraints are kept.

## Reproduction and limits

Run `python scripts/verify_bounded_certificate.py`. Its checks use
exact arithmetic for the negative-side coefficients, partial fractions
over Q(r,i), block resolvent-power identities, the finite two-prime
insertion identities, and finite trace bounds.
Use `--output results/bounded_certificate_2026-10-05.json` to save the
record. Finite examples verify algebra; they do not prove (8), the
asymptotic transfer, (15), or a new counting theorem. The analytic
arguments and their normalization still require independent review.
