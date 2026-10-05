# Weighted prime transport for slowly growing dilated locks

The [dilated geometry](dilated_rectangle_geometry.md) identifies the
progression norms and local budget needed for physical four-point locks.
This note supplies a published-input prime replacement for fixed dilations,
and uniformly for a sufficiently slowly growing prime-power family.
It permits coupled, nonuniformly weighted starts and replaces the finite
main term by the full singular series on nondegenerate locks.

The decay and dilation cutoff are qualitative. This is not the positive-
power modulus range or the prescribed-window assertion needed by the
reference's zeta reduction. Priority and independent review remain open.

## 1. The statements and weight condition

Fix epsilon>0 and C>=0. Let X,H be integers with

    X^(1/3+epsilon)<=H<=X^(1-epsilon),    L=log(3X).

For starts i,j in {X,...,2X-1}, use the balanced integer windows

    I_m=(i,i+rH],    I_n=(j,j+qH],    gcd(r,q)=1.

The reduced lock count is

    N_(i,j)(f;k,J)=sum_(r n-q m=J)
          f(m)f(m-rk)f(n)f(n-qk),

with both m slots in I_m and both n slots in I_n. For nonnegative
start weights omega_(i,j), put A=sum omega>0 and assume

    max_i sum_j omega_(i,j), max_j sum_i omega_(i,j)
      <=(A/X)L^C.                                                (1)

Weights, A and the coupled relation between starts may depend on X.
The condition allows all start pairs and a diagonal family. It excludes
a single prescribed start. Finite fixed dilations of the dyadic start
ranges are allowed after a finite dyadic covering, with adjusted constants.

**Fixed-dilation theorem.** For any fixed positive coprime r,q, there is
a sufficiently slow w(X)->infinity such that

    sum_(i,j) omega_(i,j) sum_(k,J)
      |N_(i,j)(Lambda;k,J)-Z_(i,j;r,q)(k,J)S_(w;r,q)(k,J)|^2
        =o_(r,q,epsilon,C)(A rq H^4).                             (2)

Here Z is the exact raw overlap in equation (6) of the geometry note,
and S_w is its explicit four-affine-form finite local product. Repeated
and proportional forms are retained in (2).

On locks with

    Delta_(r,q)(k,J)=k J(rq k-J)(rq k+J)!=0,

the infinite product S_(r,q)=product_p s_(p;r,q) converges. The same
bound holds with the sum restricted to these locks and S_w replaced by
S_(r,q). The singular series is not assigned on the excluded degeneracies.

**Slowly growing prime-power theorem.** There exist w(X),R(X)->infinity
with R<=w<=log log(3X), for which both assertions hold uniformly over
coprime r,q<=R, each a prime power or 1, and all weights satisfying (1).
In particular, arbitrary nonnegative combinations of the resulting
normalized estimates are valid. No explicit power or logarithmic growth
rate for R is asserted.

## 2. Published input and extraction losses

[MRSTT II, Theorem 1.3(i)](https://arxiv.org/pdf/2411.05770v2) supplies
qualitative local U^3 smallness for normalized squarefree-W-tricked
primes on almost all short intervals, with arbitrarily logarithmically
small exceptional measure. Its integer norm is divided by the interval
indicator's integer norm. Use

    W=product_(p<=w)p, c_W=W/phi(W),
    A_b(z)=c_W^(-1)Lambda(Wz+b),  gcd(b,W)=1.

Write d(w)>0 for an upper envelope of the qualitative discrepancy tending
to zero. Increase it to at least 1/w. There is no asserted numerical rate.
For a fixed dilation s in {r,q}, use quotient length

    H_s=ceil(sH/W)+6,  t_i=floor(i/W)-1,

and the analogous t_j. At position scale X/W these lengths are in the
published range with a smaller fixed epsilon, for X sufficiently large
and W=X^o(1). This remains true uniformly for s<=w<=log log(3X).

The integer local discrepancy is constant on real unit cells. The
exceptional real measure therefore controls integer quotient starts.
There are at most W original starts per quotient start. Union over the
at most R different lengths gives, for any fixed B,

    # bad original starts=O_(B,epsilon)(R X L^(-B))+O(RW).         (3)

The constant number of residue-window endpoint trims changes a cyclic
U^3 norm by O(L/sqrt(sM)), where M is the common progression group size
and sM is the pre-extraction group size. Choose M=O(H/W+1) sufficiently
large to prevent aliasing. A single supported value v has norm
|v|/sqrt(sM). Translations preserve the norm.

The exact coset extraction in the geometry note loses sqrt(s). Thus each
extracted discrepancy norm is at most

    eta=O(sqrt(R)d(w)+L sqrt(W/H)).                               (4)

The normalized presieved factor is an interval indicator after extraction,
so its cyclic norm is at most one. Triangle inequality bounds every other
prime factor by 1+eta. This retains the prime norms instead of replacing
them by their pointwise logarithmic bounds.

For fixed arbitrary r,q, the corresponding bound has sqrt(max(r,q))
in place of sqrt(R), with constants allowed to depend on r,q.

## 3. Cube splitting and the uniform local budget

First restrict primes to unit residues modulo W, giving Lambda^*.
Telescope the four factors and square with factor four. Each squared
term is the heterogeneous cube of equation (1) of the geometry note,
with two discrepancy vertices and six ordinary vertices.

Split m,n,k,t into residues modulo W and integer quotients. A typical
f vertex has quotient form

    M0-r omega1 K-r omega2 T+carry,

where its carry has magnitude O(r); the g carries have magnitude O(q).
The actual unit residue fixes its W-tricked function. Absorb the carry
as a translation of that vertex function. Next split M0 modulo r and
N0 modulo q, giving rq classes. In each, division by the step leaves
an ordinary cube with common increments -K,-T and independent bases.
After division the carry translates the support by only O(1), uniformly
in r,q, so all supports fit in the common progression group of size M.
Translate the f and g base coordinates separately to their respective
window origins; this preserves K,T and centers all eight supports.

For each admissible W-residue class the mixed cube is bounded by

    O(c_W^8 rq M^4 eta^2(1+eta)^6).

The number of admissible classes is B_W(r,q). The exact budget is

    c_W^8 B_W(r,q)/W^4=product_(p|W) beta_p(r,q).

It is uniformly bounded when r,q are coprime prime powers or 1, since
rq has at most two prime divisors. For fixed arbitrary r,q it is bounded
by a constant depending on them. Since M=O(H/W), every good start pair
satisfies

    sum_(k,J)|N_(i,j)(Lambda^*;k,J)-N_(i,j)(Lambda_w;k,J)|^2
      <=O(rq H^4 eta^2(1+eta)^6).                                (5)

All integer sums are zero extended, and the common group is chosen large
enough that supported cubes have no wraparound. No independent-start
assumption or variance factorization is used.

The deleted prime values occur only at powers of primes p<=w, at most
O(wL) integers up to 3X. Fixing one bad m vertex leaves O(qH^3) choices
for n,k,t; fixing one bad n vertex leaves O(rH^3). Their total squared
mass is therefore

    O(w(r+q)H^3 L^9)=o(rqH^4).                                  (6)

At bad start pairs the trivial bound is O(rqH^4 L^8). By (1) and (3),
their total weight is at most

    O(A[R L^(C-B)+RW L^C/X]).

Choose B>C+12. As R<=w<=log log(3X), the bad-start contribution is
o(A rqH^4). On good pairs (5) is multiplied by total weight at most A.

Choose R(w)=floor(min(w,1/d(w))), after any needed fixed-size adjustment.
Then R(w)->infinity and sqrt(R)d(w)<=sqrt(d(w))->0. Choose w(X) so
slowly that all published-input thresholds hold, W=X^o(1), and
w<=log log(3X). Equations (4)-(6) prove the replacement by Lambda_w
uniformly in the slowly growing prime-power family. For fixed arbitrary
coefficients omit the union over growing lengths and let w tend slowly
to infinity after it contains their prime divisors.

## 4. The finite main term at the raw scale

The CRT/Bezout formula in the geometry note gives error O(c_W^4 W)
per lock. Only |k|<H and

    |J-(rj-qi)|<=rqH

can contribute, so there are O(rqH^2) relevant shifts and locks. The
total squared main-term error is

    O(A rqH^2 c_W^8 W^2)=o(A rqH^4).

The squared triangle inequality proves (2).

## 5. Euler-tail moments in translated short boxes

The centered tail theorem alone is insufficient when J is centered near
rj-qi. Its divisor method extends with an explicit endpoint remainder.
For n!=0 define

    A0(n)=product_(p divides n)(1+15/p),
    Q_w(n)=sum_(p>w, p divides n)1/p.

Each prime is used once in A0. Put L2=1+log log(3X).
For an integer interval I of length K
within [1,X^2] and fixed s>0, the positive squarefree-divisor expansion
g_s(p)=(1+15/p)^s-1 gives

    sum_(n in I) A0(n)^s <= C_s K+O_s(L^(15s)).                   (7)

Indeed a divisor d occurs at most K/d+1 times. The main sum
sum_d g_s(d)/d has a convergent Euler product; the endpoint sum is at
most product_(p<=X^2)(1+g_s(p))=O_s(L^(15s)) by Mertens' estimate.

For each fixed positive integer m, the repeated-prime expansion similarly
gives

    sum_(n in I) Q_w(n)^m <= Bell(m)K/w^m+O_m(L2^m).              (8)

The main terms are grouped by equality partitions of the prime indices,
as in the centered proof. The +1 endpoint terms sum to at most
(sum_(p<=X^2)1/p)^m=O_m(L2^m). Reflection handles negative intervals,
and zero is excluded before either expansion.

Now restrict |k|<=H and |J-J0|<=rqH, where |J0|<=4RX and R<=log log(3X).
Use the four nonzero forms k,J,rqk-J,rqk+J and exclude Delta_(r,q)=0.
The first form has O(rqH) preimages and image length O(H); the others
have O(H) preimages and image length O(rqH). After normalization by
rqH^2, (7)-(8) therefore give uniformly bounded A0 moments and Q_w
moments O_m(w^(-m)+L2^m/H). Their heights are at most X^2 for large X.

When w contains all prime divisors of rq, the exceptional local head is
at most [rq/phi(rq)]^4, times at most 27 from the ordinary primes 2,3.
For prime-power r,q this is uniformly bounded. All remaining factors
have the same generic v_p and coincidence support as in the centered
proof. Thus S,S_w are bounded by a constant times product_i A0(|L_i|),
and their logarithmic tail ratio is at most

    40/w+15 sum_i Q_w(|L_i|).

Small-prime-zero cases have both products zero. Holder with A0 exponent
16, (7), and H>=X^(1/3+epsilon) bound the fourth moment of S+S_w.
Using (8) with m=4 and Cauchy--Schwarz yields

    sum_(nondegenerate box) |S_(r,q)-S_(w;r,q)|^2
      <=O(rqH^2 [w^(-2)+L2^2/sqrt(H)]).                          (9)

The implied constant is uniform for the slowly growing prime-power
family. For fixed arbitrary r,q allow dependence on them. This is a
translated short-box mean-square estimate, not a uniform pointwise tail.

Since Z<=H, (9), multiplied by A H^2, is o(A rqH^4). This proves the
full singular-series assertions by the squared triangle inequality.

## 6. Consumption and what remains open

Let E denote the actual-minus-main discrepancy on the stated domain.
For arbitrary complex consumption coefficients b_(i,j,k,J) vanishing
when omega_(i,j)=0, weighted Cauchy--Schwarz gives

    |sum b E| <=o(sqrt(A rq)H^2)
                    [sum |b|^2/omega_(i,j)]^(1/2).               (10)

The little-o is qualitative. Further growing consumption norms cannot
be absorbed without a rate or their separate normalized bound. The
start-weight condition permits logarithmic density concentration because
it is charged against the quantitative exceptional measure, not against
the qualitative discrepancy on good windows.

No prescribed start, positive-power dilation cutoff, arbitrary-log cubic
saving or zeta counting conversion is established here. The reference's
weighted modulus region and consumption-normalized volume still need
their explicit transport map. The present slowly growing family does not
cover that region. The prior arithmetic denominator-tail obstruction
remains a different estimate and is unchanged.

## Verification

`python scripts/verify_weighted_dilations.py` checks the heterogeneous
residue/carry/progression splitting, translated divisor and repeated-prime
endpoint expansions, discriminant support, translated-box multiplicities,
weighted exceptional-set counting and consumption Cauchy--Schwarz in exact
rational arithmetic. The qualitative prime input and infinite Euler bounds
are the cited theorem and the proofs above, not consequences of finite
tests. These derivations use established uniformity and singular-series
moment methods; no priority claim is made.
