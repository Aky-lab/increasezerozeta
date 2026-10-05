# The actual prime second moment with finite-frame errors

**Status:** analytic proof draft and finite algebraic calibrations.
This reconstructs a known second-moment input in the normalization of
the [actual prime-matrix decomposition](actual_prime_trace.md).
It does not prove a new zero-counting result or the higher-order cap.
Independent mathematical review remains necessary.

Retain phi, ell, ell1, L=lambda*ell, X=exp(L), d, F, H_T, D_j and O_j
from that note. Fix 0<lambda<=1. The result proved below is

    O_2(T)=O_(lambda,chi)(X/(T*ell))=o(1),
    Tr(H_T^2)/d -> lambda^2/3.                                (1)

The published prime-side second-moment result is prior work; see
[Alpoge--Furman, section 5](https://arxiv.org/html/2608.13637v1#S5).
The purpose here is to connect that input explicitly to this project's
balanced/off-balance sums, rather than to assume a model moment.

## 1. A sharper finite-frame word estimate

Let K_s be the gated circle translation from the actual-prime note,
P=P_d its d-mode projection, and eta_d=log(2*d)/d. Both leakage norms

    ||(I-P)*K_s*P||_HS, ||P*K_s*(I-P)||_HS

are O(sqrt(log(2*d))), uniformly in s, because the Fourier coefficients
of phi(u)*phi(u+s) are O(1/|m|). Products of a fixed number of K_s have
leakage bounded by the sum of the individual leakage norms.

Inserting P between j factors and telescoping gives trace errors of the
form Tr(P*K_1*P*...*K_i*(I-P)*K_(i+1)*...*K_j*P). The left factor
crosses the projection boundary once; the right factor also crosses it.
Their Hilbert--Schmidt norms are both O_j(sqrt(log(2*d))). Therefore

    |Tr(P*K_1*...*K_j*P)
       -Tr(P*K_1*P*...*P*K_j*P)| <= C_j*log(2*d).             (2)

This improves the earlier one-leakage bound and holds for open words as
well as balanced words. No prime sum has yet been estimated.

For s_1,...,s_j put p_0=0, p_i=sum_(k=1..i) s_k, S=p_j, and

    A(s_1,...,s_j)=1/L integral phi(u)*phi(u+S)
                               *product_(i=1..j-1) phi(u+p_i)^2 du,
    E_d(S)=1/d sum_(a=0..d-1) exp(-i*tau_a*S).

Zero extension of phi ensures that a surviving path has every u+p_i
inside [-L/2,L/2]. Thus the full operator product is the indicated
multiplier followed by translation by S; no surviving path wraps around
the taper interval. Taking its finite-mode diagonal gives the uniform
formula

    Tr(product_i F(s_i))/(d*L^j)
      =A(s_1,...,s_j)*E_d(S)+O_j(eta_d).                     (3)

In particular A=0 if |S|>=L, and

    0<=A<=max(1-spread(p_0,...,p_j)/L,0).

For balanced words E_d(0)=1 and (3) recovers the actual-prime note's
overlap formula with the stronger O_j(eta_d) error.

## 2. The sampled weighted mean-value inequality

For 2<=n<=X let s_n=log n. These frequencies are distinct modulo L.
Their circular separation is at least 1/(2*n): ordinary neighboring
integer logarithms have that spacing, and a wrap-around gap is at least
log 2 because the frequency zero from n=1 is absent.

[Montgomery--Vaughan, Theorem 1, equation (1.3)](https://personal.science.psu.edu/rcv4/personal/Publications/s2-8-1-73.pdf)
is a weighted Hilbert inequality for the cosecant kernel at distinct
frequencies modulo one. Apply it to s_n/L, choosing the permitted
separation weights 1/(2*n*L). For any complex coefficients x_n it yields

    |sum_(n!=m) x_n*conjugate(x_m)
                   /sin(pi*(s_n-s_m)/L)| <= C_MV*L*sum_n n*|x_n|^2.

Multiplying x_n by any unit phase leaves this bound unchanged.
The exact geometric-sum identity is

    sum_(a=0..d-1) exp(-i*tau_a*Delta)
      =[exp(-i*(T-pi/L)*Delta)
          -exp(-i*(T+(2*d-1)*pi/L)*Delta)]
                         /[2*i*sin(pi*Delta/L)].

Using the cosecant bound for both endpoint phases proves

    |sum_(n!=m) x_n*conjugate(x_m)*E_d(s_n-s_m)|
      <= C_MV*(L/d)*sum_n n*|x_n|^2.                         (4)

This is a signed estimate, not an absolute sum of reciprocal gaps. The
primary theorem is used directly on the sampling lattice; no integral
mean-value estimate is silently substituted for a discrete one.

## 3. Opposite signs in the second trace

Put w_n=Lambda(n)/sqrt(n). Applying (3) to s=log n and -t=-log m,
then changing v=u+s, gives the separable open-word amplitude

    A(s,-t)=1/L integral phi(v)^2*phi(v-s)*phi(v-t) dv.         (5)

For each fixed v take x_n=w_n*phi(v-log n) in (4). After multiplying
by phi(v)^2/L and integrating, the n!=m contribution is at most

    C_MV*(L/d)*sum_(n<=X) Lambda(n)^2.

The elementary Chebyshev bound gives sum_(n<=X) Lambda(n)^2=O(X*L).
Division by ell1^2, including both sign orientations, charges

    O(L*X*L/(d*ell1^2)).                                    (6)

The diagonal n=m is precisely the balanced sector, apart from the
finite-frame error already tracked in (3).

## 4. Same signs and the endpoint alias

For two positive logarithms S=log(n*m)>=log 4. If S>=L the open-word
amplitude is zero. Otherwise the gate and geometric sum give

    |A(log n,log m)*E_d(S)|
      <= (1-S/L)/[d*|sin(pi*S/L)|] <= C*L/d.                 (7)

Use |sin(pi*S/L)|>=2*min(S,L-S)/L. Near the alias S=L, the overlap
factor 1-S/L cancels the small denominator. Negative signs give the
conjugate estimate. Since sum_(n<=X) w_n=O(sqrt(X)), their total charge
after normalization is

    O(L*X/(d*ell1^2)).                                      (8)

No nonzero same-sign word is multiplicatively balanced.

## 5. Summing the finite-frame error

There are four sign choices per two-factor word. From (3), their full
absolute error is bounded by

    O(eta_d*ell1^(-2)*(sum_n w_n)^2)
      =O(X*log(2*d)/(d*ell1^2)).                             (9)

The balanced proxy differs from D_2 by a smaller error of the same
form with sum_n w_n^2 in place of (sum_n w_n)^2. Thus (6), (8) and (9)
bound the exact O_2=M_2-D_2. As d is comparable to L*T and L=lambda*ell,
these estimates give (1). The balanced limit D_2->lambda^2/3 was derived
in the actual-prime note using Mertens estimates and the pair overlap.

In particular the normalized Schatten-two reduction there now gives

    Tr(B_T)/d -> 1,
    Tr(B_T^2)/d -> 1+lambda^2/3.                             (10)

Indeed H_T has bounded normalized second norm by (1), and
B_T-(I-H_T) tends to zero in that norm. The Hilbert--Schmidt triangle
inequality and Cauchy--Schwarz justify the second-trace transfer.
At unit bandwidth this recovers the supplied moments 1 and 4/3.
Their counting consequences are established background, not a new record.

## 6. What this removes from the remaining cap

The degree-three target in the actual-prime note now needs only

    limsup [c_3*O_3+c_4*O_4+c_5*O_5+c_6*O_6]/6345361
      <= -3860057/63453610,                                 (11)
    (c_3,c_4,c_5,c_6)=(-3296280,-264528,6074208,3732624).

The second off-balance term tends to zero, so it can be removed without
spending any fixed asymptotic budget. Equation (11) remains unproved.
The model's predicted total is -1991744/31726805; its excess allowance
remains 49/25190. Likewise q4 retains only orders three through eight.

There is a precise limit to the boundary argument. Summing (3) absolutely
over all j-factor prime words gives only the envelope

    O_j(X^(j/2)*log(2*d)/(d*ell1^j))
      =O_(j,lambda)(T^(j*lambda/2-1)/ell^j).                 (12)

It tends to zero when j*lambda/2<=1, including the equality case. At
unit bandwidth it justifies j=2. For higher degrees at that bandwidth,
this envelope fails to tend to zero; it does not prove that the actual
boundary error is large. The balanced-word proof has a much smaller
total arithmetic weight and remains valid. Higher unbalanced traces
need a sharper summed boundary estimate or retention of the exact finite
kernel. The second-moment proof cannot simply be repeated with j=6.

## Verification and attribution

`python scripts/verify_prime_trace.py` checks exact two-factor projection
identities on a finite cyclic frame, sampled geometric sums, small
logarithmic circular-separation and mean-value calibrations, and exact
bandwidth-envelope exponents. These checks do not certify (1) or (11)
at unbounded heights. The analytic deduction above uses the cited
published weighted cosecant inequality and classical Chebyshev--Mertens
estimates; its application and normalization await independent review.
