# Actual prime matrix: balanced moments and the remaining signed target

**Status:** analytic proof draft, with finite algebraic checks. The
off-balance estimate stated below is unproved. No new zeta-zero bound
is claimed; independent mathematical review remains necessary.

This note derives the prime matrix underlying the degree-three target in
the [counting bridge](spectral_counting_bridge.md). It removes the
archimedean and pole terms in a normalized Schatten norm, evaluates the
exact multiplicatively balanced contributions through degree six, and
states the remaining estimate without assuming a four-prime reduction.

## 1. Conventions and exact entries

Retain the bridge's fixed unit-height smooth taper, ell, ell1, L=lambda*ell,
X=exp(L), d=floor(L*T/(2*pi)), tau_a=T+2*pi*a/L and C=L*ell1.
Fix 0<lambda<=1 and let T tend to infinity. The Fourier convention is

    hat(phi)(t)=integral phi(u)*exp(-i*t*u) du.

Since phi is even, its transform is unchanged by reversing that sign.
For real s define the complex symmetric matrix

    F_ab(s)=integral phi(u)*phi(s-u)
                    *exp(-i*tau_a*u-i*tau_b*(s-u)) du.

Then F(-s)=conjugate(F(s))=F(s)^*, and F(s)=0 for |s|>=L.
Fourier inversion gives, with v_a(t)=hat(phi)(t-tau_a),

    integral v_a(t)*v_b(t)*cos(s*t) dt
      = pi*(F_ab(s)+F_ab(-s)).                                (1)

The choice of sign in F corresponds to writing v_a(t) as
integral phi(u)*exp(i*t*u-i*tau_a*u) du, which is permissible by evenness.
In particular the minus sign and the factor two in the prime density are
both retained. Define the real symmetric matrix

    H_T = C^(-1) sum_(2<=n<=X) Lambda(n)/sqrt(n)
                                  *(F(log n)+F(-log n)).       (2)

Only prime powers contribute. Weil's formula is the exact identity

    B_T=D_mu+Z-H_T,
    D_mu=C^(-1) integral v(t)*v(t)^T*mu(t) dt,
    Z=C^(-1)*(v(i/2)*v(i/2)^T+v(-i/2)*v(-i/2)^T).             (3)

The explicit formula, pole absorption and taper estimates are established
background; see [Alpoge--Furman, section 2](https://arxiv.org/html/2608.13637v1#S2).
Here B_T is normalized by C, rather than by that paper's a*L^2.

There is a useful contraction bound. Let
f_a(u)=L^(-1/2)*phi(u)*exp(-i*tau_a*u). Orthogonality of the unweighted
exponentials on [-L/2,L/2] and 0<=phi<=1 show that the synthesis map
A:x->sum_a x_a*f_a has operator norm at most one. Changing u to -u in F
shows F(s)/L=A^* T_s A, where (T_s f)(u)=f(u+s) on the real line.
Consequently

    ||F(s)||_op <= L.                                         (4)

The individual summands of H_T need not be positive.

## 2. Removing the archimedean term without assuming high moments

Write ||M||_(p,d)=(Tr|M|^p/d)^(1/p). For every fixed p>=1,

    ||B_T-(I-H_T)||_(p,d) -> 0.                               (5)

Here is a proof with the dependence on the taper explicit. Freeze mu(t)
at ell1/(2*pi), and let

    A0_ab=L^(-1) integral phi(u)^2
                         *exp(i*(tau_a-tau_b)*u) du.

For t in [T/2,3T], Stirling gives mu(t)-ell1/(2*pi)=O(1).
Parseval and exponential orthogonality give

    integral |sum_a x_a*v_a(t)|^2 dt <= 2*pi*L*||x||_2^2.

Outside that interval, C2 decay and Cauchy--Schwarz give an integrable
majorant d*||x||_2^2/dist(t,[T,2T])^4. Against
|mu(t)-ell1/(2*pi)| this integrates to O(d*log(T)/T^3).
After division by C these two estimates imply

    ||D_mu-A0||_op=O(ell^(-1)+T^(-2)).

Also 0<=A0<=I and Tr(I-A0)/d=1-L^(-1)*integral phi^2=O(L^(-1)).
Thus ||A0-I||_(p,d)=O(L^(-1/p)). The pole evaluation and C2 decay give

    ||Z||_op=O(X^(1/2)/(ell1*T^3)).

This proves (5), including p=6. The taper correction need not have
vanishing operator norm; its normalized Schatten norm does vanish.

For completeness, (5) preserves any finite polynomial-trace cap of the
appropriate degree. Let q have degree r>=1 and nonzero leading
coefficient, let B,M be real symmetric, and suppose
||B-M||_(2*r,d)->0. Polynomial coercivity implies

    ||B||_(2*r,d)^(2*r) <= Cq*(1+Tr(q(B)^2)/d).

If either polynomial trace is bounded, both matrix norms are bounded
by this inequality and the Schatten triangle inequality. Telescoping
each noncommuting power and using normalized Schatten Holder gives

    ||q(B)-q(M)||_(2,d)
      <= ||B-M||_(2*r,d)
          *sum_(k=1..r) |q_k| sum_(h=0..k-1)
                 ||B||_(2*r,d)^h*||M||_(2*r,d)^(k-1-h)
      =o(1).

The product estimate first has exponent 2*r/k>=2; normalized Schatten
monotonicity then bounds its exponent-two norm. Squared norms therefore
differ by o(1). Interchanging B,M proves the reverse implication.
Thus a finite limsup cap for Tr(q3(B_T)^2)/d is equivalent to the same
cap for Tr(q3(I-H_T)^2)/d. This does not assume separate moment limits,
and it does not assert a trace difference of o(1) for unbounded statistics.

## 3. Exact cyclic formula and multiplicative balance

For j>=1 put M_j=Tr(H_T^j)/d. Expanding (2) gives

    M_j=1/(d*C^j) sum_(n_1,...,n_j<=X) product_i Lambda(n_i)/sqrt(n_i)
         *sum_(sigma_i in {-1,1})
                Tr(product_i F(sigma_i*log(n_i))).             (6)

No limiting kernel is substituted in (6). With s_i=sigma_i*log(n_i),
u_i the integration variable in edge i and cyclic indices, the trace is

    integral product_i phi(u_i)*phi(s_i-u_i)
                      *D_T(u_i+s_(i-1)-u_(i-1)) du_1...du_j,
    D_T(w)=exp(-i*T*w)*sum_(a=0..d-1) exp(-2*pi*i*a*w/L).       (7)

In particular (7) contains the global factor

    exp(-i*T*sum_i s_i)
      =exp(-i*T*log(product_(sigma_i=1) n_i/
                       product_(sigma_i=-1) n_i)).

Separate (6) exactly into D_j+O_j, where D_j contains precisely the
tuples with equal integer products on the two sign sides, and O_j
contains the remaining tuples. Reversing every sign conjugates the
trace, so both totals are real. The symbol O_j denotes the off-balance
sum, not big-O notation. Neither sector is asserted to be nonnegative.
There are configurations with up to six prime-power factors; a reduction
to four-prime locks would require a further proof.

## 4. Uniform balanced-word trace and its finite-frame error

For real s_i with sum_i s_i=0 and fixed j,

    Tr(product_i F(s_i))/(d*L^j)
      =1/L integral product_(i=0..j-1) phi(u+p_i)^2 du
            +O_j(sqrt(log(2*d)/d)),                           (8)
    p_0=0, p_i=s_1+...+s_i.

The error is uniform in all |s_i|<=L, which is essential when summing
over prime powers. One way to prove (8) is to work on the circle of
length L with orthonormal basis exp(-2*pi*i*a*u/L)/sqrt(L), a in Z.
Let K_s be the gated translation

    (K_s f)(u)=phi(u)*phi(u+s)*f(u+s),

where phi is extended by zero on the real line, and f is periodic.
Its matrix is F(s)/L with the scalar exp(-i*T*s) removed. Each K_s is
a contraction. Its Fourier coefficients are those of
psi_s(u)=phi(u)*phi(u+s), followed by a diagonal unitary. Uniformly in s,
||psi_s'||_1<=2*||phi'||_1=O(1), so |hat(psi_s)_m|<=C/|m| for m!=0.
Here the Fourier coefficient includes the factor 1/L.

If P_d projects onto a=0,...,d-1, this implies

    ||(I-P_d)*K_s*P_d||_HS^2
      <= C*sum_(m!=0) min(d,|m|)/m^2 = O(log(2*d)).

A telescoping insertion of P_d between the j factors, the contraction
bounds, and Hilbert--Schmidt Cauchy--Schwarz bound the normalized trace
error by O_j(sqrt(log(2*d)/d)). Products of several K_s have leakage at
most the sum of the individual leakage norms.

Because a surviving path never exits the zero-extended taper interval,
the closed gated translation is multiplication by
product_i phi(u+p_i)^2; no surviving path wraps around the interval.
Every Fourier basis vector has squared modulus 1/L, so the trace of its
P_d compression divided by d is the integral in (8). The scalar phases
cancel when sum_i s_i=0. This proves (8) without replacing the finite
matrix kernel inside an unbalanced prime sum.

Finally the fixed-width unit ramps imply, uniformly in the closed walk,

    1/L integral product_i phi(u+p_i)^2 du
      =(1-(max_i p_i-min_i p_i)/L)_+ +O_j(L^(-1)).             (9)

The intersection of the translated length-L intervals has that length;
the union of their transition collars has total length O_j(1).

## 5. Balanced prime moments through degree six

Unique factorization, (4), (8)-(9) and Mertens estimates prove

    D_1=0,
    D_(2r+1)=o(1) for r>=1,
    D_(2r)->lambda^(2r)*{2^r} for fixed r>=1.                 (10)

The pair integral {2^r} is exactly the integral defined in
[the pair-cycle proof](paired_cycle_flow_polytopes.md). Here is the
arithmetic argument, including the exceptional configurations.

Partition a balanced word by its underlying prime bases. Every block
has at least two slots and both signs occur. For a two-slot block the
exponents agree, giving weight (log p)^2/p^a. The a=1 terms have total
O(L^2); the a>=2 terms have uniformly bounded total.

For a block of b>=3 slots, the common sum of positive and negative
exponents is at least ceil(b/2)>=2. Counting compositions of that sum
shows that its total positive weight, summed over p and all exponents,
is bounded by a constant depending only on b: it is dominated by
C_b*sum_p (log p)^b/p^2. Thus any such block loses at least three powers
of L relative to the normalization ell1^j. A higher-power two-slot block
loses two powers. These estimates, the finitely many slot partitions
and |Tr(product F)|/d<=L^j show that all configurations except pairs of
distinct ordinary primes contribute o(1). This also proves the odd claim.

For each remaining pairing, orient the first endpoint by sigma=+/-1.
The normalized signed prime measure satisfies

    L^(-2) sum_(p<=exp(L)) (log p)^2/p
                *[delta_(log(p)/L)+delta_(-log(p)/L)]
       -> |x| dx on [-1,1].

Partial summation of sum_(n<=y) Lambda(n)/n=log y+O(1) proves this:
the difference between Lambda(n)^2 and Lambda(n)*log n has bounded
weighted sum on higher prime powers. Coincidences of two pair bases
have weight O(L^(-4))*sum_p (log p)^4/p^2 and vanish. Apply the product
measure limit to the bounded continuous overlap (9). The uniform error
in (8) is harmless because the total balanced weight is O_j(L^j).
The factor (L/ell1)^j tends to lambda^j. This proves (10).

In particular, at unit bandwidth the balanced limits through degree six
are

    (D_0,...,D_6)->(1,0,1/3,0,4/15,0,32/105).                (11)

This is an evaluation of a specified part of the actual arithmetic
expression. It is not an evaluation of its full moments.

Also M_1=o(1), by a direct geometric-sum bound. The diagonal of F(s) is
exp(-i*tau_a*s)*g_phi(s), with 0<=g_phi(s)<=L-s for 0<s<L.
For s>=log 2,

    g_phi(s)*|sum_(a=0..d-1) exp(-2*pi*i*a*s/L)| <= C*L^2.

Indeed |sin(pi*s/L)|>=2*min(s,L-s)/L. Chebyshev partial summation
gives sum_(n<=X) Lambda(n)/sqrt(n)=O(sqrt(X)). Division by C*d yields
M_1=O(sqrt(X)/(ell1*T))=o(1). Since D_1=0, O_1=o(1) as well.

The [finite-frame second-moment derivation](prime_second_moment.md)
also proves O_2=O_(lambda,chi)(X/(T*ell))=o(1), using the published
weighted cosecant inequality. It supplies a stronger O_j(log(2*d)/d)
open-word frame estimate and retains all errors when summing the
two-prime expression. Thus the full second moment is lambda^2/3.
This recovers a known arithmetic input; higher off-balance moments
remain separate obligations.

## 6. The precise unproved off-balance target

With q3 from the bridge, put p3(h)=q3(1-h). Then

    p3(h)=(-277-708*h+1572*h^2+1932*h^3)/2519,
    p3(h)^2=sum_(j=0..6) c_j*h^j/6345361,
    (c_0,...,c_6)=(76729,392232,-369624,-3296280,
                                  -264528,6074208,3732624).

For the model, binomial centering of its certified moments gives

    (h_0,...,h_6)=(1,0,1/3,0,1/4,-1/36,61/252).

This is only a model calculation. Subtracting (11) gives its suggested
off-balance values

    (o_2,...,o_6)=(0,0,-1/60,-1/36,-79/1260),

The second off-balance value is established in the linked proof draft;
the higher values have not been established for actual primes. The balanced part
of the squared certificate tends exactly to

    A_diag=(c_0+c_2/3+4*c_4/15+32*c_6/105)/6345361
          =5102709/31726805 = 0.16083274064....

It exceeds the desired cap 0.1. The required cancellation is therefore
substantial: ignoring all off-balance terms cannot prove the cap.
By (5), (10), O_1=o(1) and O_2=o(1), the actual degree-three 80% target is equivalent
to the following **unproved signed arithmetic estimate**:

    limsup_(T->infinity) [sum_(j=3..6) c_j*O_j(T)]/6345361
      <= 1/10-A_diag
       =-3860057/63453610 = -0.06083274064....                 (12)

The full model predicts -1991744/31726805, approximately -0.062777957,
leaving the same exact excess allowance 49/25190 as the bridge.
Equivalently one may bound sum c_j*(O_j-o_j)/6345361 by that allowance.
The sum can begin at j=3 because its j=2 term tends to zero.
No absolute-value bound or separate autocorrelation theorem supplies
the negative upper bound in (12) automatically.

The next analytic obligation is now (12), with O_j defined by (6)-(7).
Any regional progression replacement must first derive its coefficients
from these sums and account for all configurations. The current
four-prime transport results are not silently used for six-prime words.

### Degree-four comparison

The same argument applies to q4 using p=8 in (5) and r=4 in (10).
With p4(h)=q4(1-h), its balanced squared trace tends to

    A_diag4=38601698857019973/132097166600162405
           =0.29222200484....

The 80% cap is therefore equivalent to a signed off-balance cap

    1/10-A_diag4=-10156792878801493/52838866640064962
                =-0.19222200484....

Its model predicts -28653310482603548/132097166600162405,
approximately -0.216910864, leaving the larger excess allowance
40129409/1625405590. Matching that model value would give the bridge's
84.938% conditional conversion. These figures compare the arithmetic
requirements of the two certificates; they prove neither off-balance
cap. The quartic route needs words through length eight and permits
more error, while the cubic route needs words only through length six.
The second off-balance term vanishes in both routes. The finite-frame
note quantifies why its absolute boundary-error argument cannot be
promoted to higher unbalanced moments at unit bandwidth.

## Verification and scope

`python scripts/verify_prime_trace.py` independently checks Fourier and
convolution factors on an exact finite cyclic frame, the balanced/off-
balance word decomposition, closed-loop overlaps, noncommuting polynomial
telescoping, and the rational centered-model and signed budgets.
Its [record](../results/prime_trace_2026-10-05.json) states this finite scope.
The analytic proofs above require mathematical review; finite checks do
not establish the unbounded-height off-balance estimate (12).
