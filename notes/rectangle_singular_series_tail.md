# Mean-square Euler tails for rectangle singular series

The finite four-point main term in
[averaged rectangle transport](averaged_rectangle_transport.md) can be
replaced by the full singular series on nondegenerate shifts. An elementary
divisor-moment argument gives a uniform mean-square Euler-tail bound. It
does not require a uniform pointwise tail estimate; such an estimate would
fail for a fixed prime-factor cutoff and unrestricted heights.

## 1. Definitions and conclusions

For integer k,h put

    Delta(k,h)=k*h*(k-h)*(k+h),
    D_p(k,h)=#{0,-k,h,h-k} mod p,
    s_p(k,h)=[1-D_p(k,h)/p]/(1-1/p)^4.

For integer w>=5 define K_w=product_(p<=w)s_p. If Delta!=0, all
four offsets are distinct as integers and

    S_4(k,h)=product_p s_p(k,h)

converges. It is zero when one of the factors at 2 or 3 is zero.
Let Omega_Y be the integer pairs with |k|,|h|<=Y and Delta!=0.

**Euler-tail theorem.** For every fixed positive integer r,

    sum_((k,h) in Omega_Y) |S_4(k,h)-K_w(k,h)|^r
      <=C_r Y^2 w^(-r),                                          (1)

uniformly for integer Y>=1 and w>=5. In particular the squared tail is
O(Y^2/w^2). The constants depend on the fixed moment r, not on the height
or cutoff. This is an arithmetic Euler-product estimate; no prime-tuple
conjecture is used.

Combining r=2 with the published-input rectangle replacement gives the
following actual prime statement. For X^(1/3+epsilon)<=Y<=X^(1-epsilon),
let N_j(Lambda;k,h) be the four-prime window count and
Z_Y(k,h)=(Y-|k|-|h|)_+. Then

    sum_(X<=j<2X) sum_((k,h) in Omega_Y)
      |N_j(Lambda;k,h)-Z_Y(k,h)S_4(k,h)|^2=o_(epsilon)(X Y^4).       (2)

The sum may equivalently restrict to |k|+|h|<Y, as both terms vanish
outside it. The saving remains qualitative because the prime-uniformity
input is qualitative. Individual starts, concentrated or logarithmically
growing consumption weights, variable dilations and degenerate offsets
are not supplied by (2).

## 2. Separate the ordinary and coincidence factors

For p>=5 put

    v_p=(1-4/p)/(1-1/p)^4.

Then 0<v_p<1, and the local factor is

    s_p=v_p A_p,
    A_p=(p-D_p)/(p-4).

The offsets collide modulo p exactly when p divides Delta. Thus
A_p=1 otherwise, and

    1<=A_p<=1+15/p  when p|Delta.                                  (3)

The factor 15 follows from (4-D_p)/(p-4)<=3/(p-4)<=15/p for p>=5.
The factors at 2 and 3 are nonnegative and their product is at most 27.

For any nonzero integer n define

    A(n)=product_(p|n)(1+15/p),
    Q_w(n)=sum_(p>w, p|n) 1/p.

Equations (3) and v_p<1 imply, with the four nonzero linear forms
L1=k, L2=h, L3=k-h, L4=k+h,

    S_4(k,h), K_w(k,h) <=27 product_(i=1..4) A(|L_i|).              (4)

If the small-prime factor vanishes, the tail is zero. Otherwise all
factors are positive and the logarithm of the tail ratio is defined.
Expanding -log(1-4/p)+4log(1-1/p), or bounding the remainder after its
linear term, gives |log v_p|<=40/p^2. Therefore

    |log(S_4/K_w)|
      <=40/w+15 sum_i Q_w(|L_i|).                                 (5)

Only finitely many coincidence factors occur, since each L_i is nonzero.
The ordinary factors converge absolutely, justifying (5) for the infinite
product as well. For positive u,v, |u-v|<=max(u,v)|log(u/v)|. Hence

    |S_4-K_w|
      <=C(S_4+K_w)[1/w+sum_i Q_w(|L_i|)].                          (6)

## 3. One-dimensional divisor moments without endpoint loss

For any fixed s>0, expand

    A(n)^s=sum_(d|n, d squarefree) g_s(d),
    g_s(p)=(1+15/p)^s-1,

where g_s is multiplicative on squarefree integers and zero otherwise.
Every coefficient is nonnegative. Since g_s(p)=O_s(1/p),

    (1/H)sum_(1<=n<=H) A(n)^s
      <=sum_d g_s(d)/d
       =product_p [1+g_s(p)/p]=C_s<infinity.                        (7)

This uses floor(H/d)<=H/d, so there is no primorial-period remainder.

For a positive integer m, expand Q_w(n)^m and count multiples of the
least common multiple of its prime indices. Group each prime tuple by
the partition of slots having equal primes. If a block has size b, its
prime sum is at most

    sum_(p>w) p^(-b-1)<=1/(b w^b).

Ignore distinctness between different blocks to get an upper bound.
There are Bell(m) set partitions, and the block sizes sum to m. Thus

    (1/H)sum_(1<=n<=H) Q_w(n)^m<=Bell(m) w^(-m),                   (8)

uniformly in integer H>=1 and w>=5. For m=4 the right side is 15/w^4.
Repeated prime indices are included; pretending that all indices are
distinct would give an incomplete moment bound.

## 4. Transfer to the rectangle box and prove the tail theorem

Each of k,h,k-h,k+h takes values of magnitude at most 2Y on the box
|k|,|h|<=Y. Every nonzero value has O(Y) preimages. Consequently any
nonnegative function F satisfies

    (1/Y^2)sum_((k,h) in Omega_Y) F(|L_i(k,h)|)
      <=(C/Y)sum_(1<=n<=2Y) F(n).                                 (9)

The exclusion Delta=0 can only decrease this nonnegative sum. Holder's
inequality across the four forms, (4), (7) and (9) give

    (1/Y^2)sum_(Omega_Y) (S_4+K_w)^(2r)=O_r(1).                    (10)

Here (7) is used with s=8r, since each of four factors is raised to
2r and Holder uses exponent four. No independence of the forms is assumed.

Similarly (8)-(9), with m=2r, and the elementary inequality for a sum
of five nonnegative terms give

    (1/Y^2)sum_(Omega_Y) [1/w+sum_i Q_w(|L_i|)]^(2r)
      =O_r(w^(-2r)).                                              (11)

Apply Cauchy--Schwarz to the r-th power of (6), using (10)-(11).
This proves (1) for every fixed integer r. The diagonal lines must be
excluded before using nonzero-value multiplicities or Q_w.

## 5. Full singular-series prime transport

Use the sufficiently slow w(X)->infinity from the
[qualitative rectangle theorem](averaged_rectangle_transport.md).
That theorem, restricted to Omega_Y, bounds the actual-minus-finite-main
squared error by o(XY^4). By Z_Y<=Y and (1) with r=2,

    sum_j sum_(Omega_Y) |Z_Y(S_4-K_w)|^2
      <=C X Y^4/w(X)^2=o(XY^4).

The squared triangle inequality proves (2). This supplies the full
Hardy--Littlewood four-point main term for the stated joint average.
It does not imply the four-prime conjecture at prescribed offsets.

## 6. Why a pointwise tail cannot replace this argument

Fix w>=5. For z>w let

    k=6 product_(w<p<=z)p,  h=2k.

These shifts are nondegenerate, and the factors at 2 and 3 have product
27. At every prime in (w,z] all offsets coincide, so D_p=1 and
s_p=(1-1/p)^(-3). At p>z there are no coincidence primes. Therefore

    S_4(k,h)/K_w(k,h)
      =product_(w<p<=z)(1-1/p)^(-3) product_(p>z)v_p.

The second product tends to one; the first grows without bound by the
standard Mertens product law. Also K_w is bounded below here by the
positive constant 27 product_(p>=5)v_p. Thus the absolute pointwise tail
is unbounded at fixed w as these rare shifts grow. Their rarity is exactly
what the box moment argument accommodates.

The tail in this note is a prime-factor Euler-product tail averaged over
rectangle shifts. It does not estimate the class-subtracted Ramanujan
denominator tail in the separate
[arithmetic-core cutoff obstruction](fixed_cutoff_obstruction.md), nor
the singular weights and normalization in that moment reduction.

## Literature and verification

Singular-series averaging and Euler-product moment methods are established
background. [Ford's proof of Gallagher's estimate](https://arxiv.org/pdf/1108.3861)
treats unrestricted distinct tuples, and
[Kowalski's Euler-product averaging framework](https://people.math.ethz.ch/~kowalski/singular-series-distribution.pdf)
studies singular-series moments and distributions. The present note gives
the explicit bound and consumption for the restricted two-parameter
rectangle family; no priority claim is made.

`python scripts/verify_rectangle_tail.py` independently checks local
collision factors, the discriminant support, divisor-moment expansions,
repeated-prime partition bookkeeping and finite tail decompositions in
exact rational arithmetic. The infinite estimates are the proofs above;
the actual-prime input is the cited theorem in the preceding note.
