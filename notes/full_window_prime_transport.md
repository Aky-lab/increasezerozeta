# Quantitative prime transport at prescribed full-length windows

The [weighted short-window theorem](weighted_dilated_prime_transport.md)
averages starts and has a qualitative saving. For windows whose lengths
are comparable to their positions, a global published uniformity estimate
instead gives a quantitative result at every prescribed start. This note
derives that result for a double-logarithmic dilation range, including
arbitrary coprime dilations and the full singular-series main term.

## 1. Statement and normalization

Fix K>=1. Let H tend to infinity through integers, put L2=log log H,
and let r,q be coprime positive integers. Use the balanced windows

    I_m=(i,i+rH], I_n=(j,j+qH],
    0<=i<=KrH, 0<=j<=KqH,                                  (1)

with integer starts. Let N_(i,j)(Lambda;k,J) be the reduced physical
lock count in section 1 of the
[dilated geometry note](dilated_rectangle_geometry.md):

    sum_(rn-qm=J) Lambda(m)Lambda(m-rk)Lambda(n)Lambda(n-qk),

with both m slots in I_m and both n slots in I_n. The exact raw overlap
Z_(i,j;r,q)(k,J), finite product S_(w;r,q) and infinite product S_(r,q)
are as in sections 4-5 of that note and section 5 of the weighted note.
Set w=floor(sqrt(L2)), W=product_(p<=w)p and c_W=W/phi(W).

There is a constant kappa>0 such that, uniformly for
1<=r,q<=R=floor(L2^kappa), all starts in (1) satisfy

    sum_(k,J) |N_(i,j)(Lambda;k,J)-Z S_(w;r,q)(k,J)|^2
      <<_K rq H^4 L2^(-kappa).                             (2)

For Delta_(r,q)=kJ(rqk-J)(rqk+J)!=0, one also has

    sum_(k,J:Delta!=0) |N_(i,j)(Lambda;k,J)-Z S_(r,q)(k,J)|^2
      <<_K rq H^4 L2^(-kappa).                             (3)

The finite statement retains degenerate locks. The infinite product is
asserted only on the nondegenerate domain. The scale is the raw lock
scale rqH^4, without a consumption-normalized zeta volume factor.

The exponent depends on a positive exponent in the published cubic
uniformity theorem. No numerical value for that exponent is claimed.
The result gives a double-logarithmic saving, not arbitrary powers of
log H. The dilation range is subpower and much smaller than the range
needed to recover the [leading arithmetic core](dilation_core_coverage.md).

## 2. Published input and restriction by prefix subtraction

[Tao and Teravainen, Quantitative bounds for Gowers uniformity of the
Mobius and von Mangoldt functions, JEMS 27 (2025), Theorem 1.4](https://ems.press/journals/jems/articles/13625437)
proves, for k=3, some c0>0 and 2<=z<=exp((log N)^(1/10)),

    ||Lambda-Lambda_Cramer,z||_(U^3[N])
      << (log log N)^(-c0)+z^(-c0).                        (4)

The [author version](https://arxiv.org/pdf/2107.02158v4) defines the
interval norm by zero extension on the integers divided by the indicator
norm. Its model uses p<z. Taking z=w+1/2 gives exactly our p<=w model
Lambda_w=c_W 1_(gcd(n,W)=1). This statement concerns the ordinary model;
there is no omitted exceptional-character term in (4).

Let a=min(c0/2,1/4), kappa=a/2. Choose a common integer
M>8(K+2)H with M=O_K(H). For s in {r,q}, view zero-extended sequences
on Z/(sM)Z. If sqrt(H)<=N<=(K+1)sH, the source norm conversion gives

    ||(Lambda-Lambda_w)1_[1,N]||_(U^3(Z/(sM)Z))
      <<_K L2^(-a).                                       (5)

Indeed log log N is comparable to L2, the cutoff condition holds, and
the integer indicator has O(N^4) cubes. Dividing by (sM)^4 and taking
an eighth root costs only O_K(1). The group size exceeds eight times
the support diameter, so its supported cyclic cubes are exactly the
integer cubes. This conversion uses a common physical group for both
prefixes, not a Fourier projection onto the window.

For 0<=N<sqrt(H), the elementary pointwise bound O(log H) and the same
cube count instead give O_K(log H H^(-1/4))=o(L2^(-a)). Subtract the
prefixes at i+sH and i (or j+sH and j). Their exact difference is the
window-restricted discrepancy, and the norm triangle inequality proves

    ||(Lambda-Lambda_w)1_(v,v+sH]||_(U^3(Z/(sM)Z))
      <<_K D:=L2^(-a)                                    (6)

uniformly for 0<=v<=KsH. This is where the full-length assumption matters:
for a short interval far from zero, the two large prefixes would lose
their support-density normalization relative to H.

## 3. Extract progressions and apply the physical lock criterion

The coset extraction inequality (2a) of the geometry note gives

    max_(c mod s) ||(Lambda-Lambda_w)(c+s*.)1_window||_(U^3(Z/MZ))
      <<_K sqrt(s)D.

Each progression has exactly H slots. Translation into a common support
interval preserves its norm. The model progression norm is at most c_W
by its pointwise bound, and the prime norm is at most c_W+O_K(sqrt(R)D).
The deterministic criterion (2) of the geometry note therefore yields

    sum_(k,J) |N(Lambda)-N(Lambda_w)|^2
      <<_K rq H^4 R D^2 (c_W+sqrt(R)D)^6.                  (7)

Mertens' product bound gives c_W=O(1+log w)=O(1+log L2).
With R<=L2^(a/2), the normalized right side is
O_K(L2^(-3a/2)(1+log L2)^6)=O_K(L2^(-kappa)).
The extraction loss has been paid explicitly. No uniform boundedness
of the exceptional-prime residue budget for unrestricted r,q is assumed.

For the model, exact CRT counting on the Bezout overlap gives

    N_(i,j)(Lambda_w;k,J)=Z S_(w;r,q)(k,J)+O(c_W^4 W).

The relevant lock box is |k|<=H, |J-J0|<=rqH,
J0=rj-qi, with |J0|<=2KrqH. There are O(rqH^2) locks.
Since W=exp(O(sqrt(L2))), their total squared CRT remainder is
O(rqH^2 c_W^8 W^2)=o(rqH^4 L2^(-kappa)). This proves (2).

## 4. Replace the finite product by the infinite product

The translated positive-moment proof in section 5 of the weighted note
applies to the same lock box with height O_K(R^2H)<=H^2 for large H.
The four nonzero discriminant forms are k,J,rqk-J,rqk+J. Their image
intervals have lengths O(H) and O(rqH), and their preimage multiplicities
are O(rqH) and O(H), respectively. Thus their normalized divisor
moments are bounded uniformly; the endpoint error
O((log H)^(15s)/H) tends to zero for every fixed moment order s.
Their prime-reciprocal moments are O_m(w^(-m)+(1+L2)^m/H).

Here all primes dividing rq are at most R<=w. Write t0=rq/phi(rq).
The small-prime head is bounded by 27t0^4. The original proof's
fourth-moment bound therefore acquires t0^16, and Cauchy--Schwarz
acquires t0^8. Explicitly,

    sum_(nondegenerate box) |S_(r,q)-S_(w;r,q)|^2
      <<_K rqH^2 t0^8 [w^(-2)+(1+L2)^2/sqrt(H)].            (8)

The logarithmic tail ratio and repeated-prime expansion are unchanged:
at ordinary primes they are supported only on the four discriminant
forms. Small-prime-zero cases have both products zero. This accounts for
unrestricted coprime moduli instead of silently using the prime-power
head bound.

Mertens' bound gives t0<<1+log log(3rq)<<1+log L2. Since
kappa<=1/8 and w is comparable to sqrt(L2), the normalized right side
of (8) is O_K(L2^(-kappa)). Multiply by Z^2<=H^2 and combine with
(2) to prove (3).

## 5. Weights, the reference windows and limits

For any nonnegative start weights omega_(i,j) supported in (1), with
total mass A, summing (2) or (3) gives the same upper bound multiplied
by A. No marginal density or exceptional-set condition is needed.
For a fixed start, Cauchy--Schwarz gives the consumption estimate

    |sum_(k,J) b_(k,J) E_(k,J)|
      <<_K sqrt(rq) H^2 L2^(-kappa/2) ||b||_2.             (9)

Weighted consumption likewise uses sqrt(A rq)H^2 L2^(-kappa/2)
times [sum |b|^2/omega]^(1/2), with b zero where omega is zero.

The reference's default I_m=[X,4X] and proportional I_n have full size
at their positions. For fixed r,q, take H near 3X/r; integer endpoint
rounding changes O_(r,q)(1) slots per progression. Its U^3 norm error
is O_(r,q)(log H/sqrt(H)), negligible in (7). Thus prescribed starts
alone are not an obstruction for those fixed-dilation full windows.
For growing r,q the theorem is stated only for the balanced integer
windows (1); additional endpoint conventions must be tracked separately.

Positive-power dilations, genuine short prescribed windows, a signed
actual-minus-model arithmetic-core estimate, and the reference's
consumption normalization remain open. The theorem cannot charge those
obligations to its double-logarithmic saving.

## Verification and attribution

`python scripts/verify_full_window_transport.py` checks sparse integer
and cyclic cube sums, prefix subtraction, subgroup extraction, exact
physical lock variance and consumption weights. It also checks the
cutoff convention and the exponent inequalities. Finite checks do not
prove (4) or the asymptotic Euler bounds. The theorem is a deduction
from published global uniformity and the repository's deterministic
geometry and translated Euler-tail proofs. Priority and independent
review of this specific deduction remain open.
