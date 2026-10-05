# Qualitative prime transport for averaged rectangle locks

The restricted rectangle in the [resolved frame](resolved_lock_frame.md)
can be treated at its natural Y^4 squared scale by using the Gowers norms
of the other factors, rather than their pointwise logarithmic bounds.
Residue-class admissibility compensates for the W-trick normalization.
This gives an actual four-prime replacement averaged over integer starts
and both rectangle shifts. Its saving is qualitative.

This is a deduction from published short-interval Gowers uniformity. It
does not assert a new priority claim, a prescribed-window theorem,
arbitrary logarithmic savings, or the weighted and dilated zeta interface.

## 1. Statement

Fix epsilon>0 and integers X,Y satisfying

    X^(1/3+epsilon)<=Y<=X^(1-epsilon).

There is a function w=w(X) tending to infinity sufficiently slowly such
that, with W=product_(p<=w)p and c_W=W/phi(W), the model

    Lambda_w(n)=c_W 1_(gcd(n,W)=1)

has the following property. For j in J={X,...,2X-1}, define

    N_j(f;k,h)=sum_n f(n)f(n-k)f(n+h)f(n+h-k)
                   product_(t in {0,-k,h,h-k})1_(j<n+t<=j+Y).

Then, for signed integer k,h,

    sum_(j in J) sum_(k,h) |N_j(Lambda;k,h)-N_j(Lambda_w;k,h)|^2
      =o_(epsilon)(X Y^4).                                        (1)

Only |k|+|h|<Y can contribute. The little-o is uniform in the stated
Y range after a sufficiently slow choice of w(X). It has no asserted
logarithmic decay rate or effective height threshold.

Put Z_Y(k,h)=(Y-|k|-|h|)_+ and

    K_w(k,h)=product_(p<=w)
       [1-D_p(k,h)/p]/(1-1/p)^4,
    D_p(k,h)=#{0,-k,h,h-k} mod p.

The explicit finite main term also satisfies

    sum_j sum_(k,h)
       |N_j(Lambda;k,h)-Z_Y(k,h)K_w(k,h)|^2=o_(epsilon)(X Y^4).       (2)

This finite local product includes repeated offsets on k=0, h=0 and
k=+-h. No infinite four-distinct-prime singular series is substituted
on those degeneracies. Uniform pointwise replacement by the infinite series
is not claimed. The subsequent [Euler-tail proof](rectangle_singular_series_tail.md)
does establish its replacement in the joint mean square on nondegenerate
shifts.

## 2. Published input and the norm used

[MRSTT II, Theorem 1.3(i)](https://arxiv.org/pdf/2411.05770v2) supplies,
for the W-tricked functions

    A_b(n)=c_W^(-1) Lambda(Wn+b),  gcd(b,W)=1,

the bound max_b ||A_b-1||_(U^3(I))=o_(w->infinity)(1) on almost all
intervals I=(t,t+H], in the range T^(1/3+epsilon')<=H<=T^(1-epsilon').
The exceptional real starts in [T,2T] have measure
O_(A,epsilon')(T log(3T)^(-A)) for any fixed A. X is large in terms of
w; a slow diagonal choice of w is part of the statement, not a numerical
rate. The prime majorants and inverse argument proving this input are
the published work. No new majorant construction is assumed here.

For a finite interval I, the published local norm is the integer Gowers
norm of the zero-extended function divided by that of 1_I. Equivalently,
use a cyclic group large enough to embed the cubes without aliasing;
the common normalization cancels from this ratio. Consequently

    ||A_b 1_I||_(U^3) <=(1+delta)||1_I||_(U^3),
    ||(A_b-1)1_I||_(U^3) <=delta||1_I||_(U^3),                     (3)

when the local discrepancy is at most delta. Gowers triangle inequality
proves the first bound from the second. Thus the other prime factors have
bounded cubic norms on a good window despite their unbounded pointwise
values. The [Gowers--Cauchy--Schwarz inequality](https://people.maths.ox.ac.uk/tillmann/GreenTao.pdf),
equation (5.5), bounds a mixed cube by the product of its eight U^3 norms.

## 3. A uniformly bounded residue budget

Let B_W be the number of (x,a,b,c) modulo W for which all eight
x+omega1*a+omega2*b+omega3*c, omega in {0,1}^3, are coprime to W.
The normalization budget is

    c_W^8 B_W/W^4=product_(p|W) beta_p.                            (4)

Each cube form has a zero set of probability 1/p. Distinct forms have
linearly independent coefficient vectors, so two zero sets intersect
with probability 1/p^2. Bonferroni already shows
beta_p=1+O(1/p^2) from above for large p, which suffices to bound (4)
uniformly in W.

An exact finite rank calculation sharpens this. The cube coefficient
rows are (1,omega1,omega2,omega3). Inclusion-exclusion gives

    B_p=sum_(S subset {0,1}^3) (-1)^|S| p^(4-rank_p(S)).

All nonzero square minors of these rows have magnitude at most 2.
Hence their ranks are the same over every odd prime field as over Q.
The signed rank coefficients are (1,-8,28,-44,23), giving

    B_p=p^4-8p^3+28p^2-44p+23,  p odd,
    B_2=1,
    beta_p=B_p p^4/(p-1)^8.                                      (5)

The standard-library checker enumerates all 256 subsets with rational
rank and independently enumerates small prime fields. Over F_2 the signed
rank coefficients are (1,-8,28,-42,21), whose evaluation at 2 is 1.
The small factors are beta_2=16 and beta_3=81/32. For p>=5,

    beta_p<=1+100/p^3.

Indeed the numerator minus (p-1)^8 is
12p^5-47p^4+56p^3-28p^2+8p-1, bounded above by 15p^5;
(1-1/p)^8>0.16 for p>=5. Thus, for example, (4) is at most
(81/2)exp(25/8), using sum_(n>=5)n^(-3)<=1/32.
No factor growing with W is lost when summing the admissible classes.

## 4. Integer starts, residue windows and boundary terms

Choose w so slowly that W=X^o(1), and set

    H=ceil(Y/W)+6,  t=floor(j/W)-1.

Every W-tricked coordinate needed by a prime in (j,j+Y] lies in
I_t=(t,t+H]. At most a fixed number of endpoints distinguish its exact
residue-class interval from I_t. At scale T approximately X/W, the length
H lies in the published range with epsilon'=epsilon/2 for large X.
Small dyadic endpoint adjustments affect only O(W) original starts.

Since H is an integer, the local discrepancy in (3) is constant on
real unit cells of t. The published exceptional measure therefore
controls the integer t as well. At most W starts j use any one t, so
there are O_A(X log(3X)^(-A))+O(W) bad original starts.

Embed the residue-coordinate functions on a cyclic group of size
M=O(H), chosen sufficiently large to avoid aliasing after the bounded
carry shifts. A single supported value v has U^3 norm |v|/sqrt(M):
there is only one cube entirely supported at that point. Removing the
bounded number of endpoints changes the norm by O(log(3X)/sqrt(H)).
Translations preserve it. Thus (3) applies to every translated, trimmed
residue function with error

    delta'=delta+O(log(3X)/sqrt(H)).                               (6)

This tends to zero. It does not introduce a pointwise logarithmic factor
on the six other cube vertices.

## 5. Four-factor telescoping and the cube decomposition

First delete prime weights at nonunit residues modulo W. Write

    Lambda^*(n)=Lambda(n)1_(gcd(n,W)=1).

Expand N_j(Lambda^*)-N_j(Lambda_w) by telescoping its four factors,
and apply |sum_(i=1)^4 z_i|^2<=4 sum_i |z_i|^2. Each squared term is
a mixed eight-vertex cube, using increments -k,h,n'-n. The difference
appears at two vertices; the other six are either prime or model factors.

Split the four cube variables into residues modulo W and integer
quotients. Nonadmissible residue classes vanish. In an admissible class,
each vertex equals W times its quotient cube form plus its unit residue.
The carry is a bounded vertex-dependent translation, absorbed into that
vertex's function. After extracting c_W^8, (3), (6) and Gowers--Cauchy--
Schwarz bound each class by

    O(c_W^8 H^4 (delta')^2(1+delta')^6).

Sum over the B_W admissible classes. Equations (4)-(5) and H=O(Y/W)
give, at each good j,

    sum_(k,h)|N_j(Lambda^*;k,h)-N_j(Lambda_w;k,h)|^2
      =O(Y^4 (delta')^2(1+delta')^6).                              (7)

The class sums are finite signed mixed-cube expressions; taking their
absolute values before using the cube inequality is legitimate. No
independence of positions or factorization of resolved variances is used.

The deleted weights can only occur at powers of primes p<=w. There are
O(w log(3X)) such integers in the whole window, even by the coarse global
count up to 3X. A cube containing one specified vertex has O(Y^3)
remaining choices. Bounding its eight weights by O(log(3X)^8) therefore
gives uniformly

    sum_(k,h)|N_j(Lambda;k,h)-N_j(Lambda^*;k,h)|^2
      =O(w Y^3 log(3X)^9).                                        (8)

Here the weights are nonnegative, so expansion of the difference square
is bounded by the cube mass containing at least one deleted vertex.
For sufficiently slow w, (8) is o(Y^4).

At bad starts the trivial bound is O(Y^4 log(3X)^8). Their number from
section 4 makes their total negligible by choosing the exceptional-set
exponent A large enough. Summing (7)-(8) and taking w(X) slowly to
infinity proves (1). The qualitative norm bound and all finitely many
range conditions can be met simultaneously by a diagonal choice of w.

## 6. The finite four-point main term

The four offsets have span |k|+|h|. Their common window is an integer
interval of length Z_Y(k,h). The CRT gives exactly W times the product
of admissible local proportions as its number of good residues per
period. Each residue class has count Z/W+O(1), so, uniformly in j,k,h,

    N_j(Lambda_w;k,h)
      =Z_Y(k,h)K_w(k,h)+O(c_W^4 W).                                (9)

There are O(Y^2) relevant shifts. The squared error in (9), summed over
starts and shifts, is O(X Y^2 c_W^8 W^2)=o(XY^4), since W=X^o(1)
and Y is a fixed positive power of X. The squared triangle inequality
proves (2).

## 7. Consumption and remaining scope

For arbitrary weights a_(j,k,h), Cauchy--Schwarz yields from (2)

    |sum_(j,k,h) a_(j,k,h)
          [N_j(Lambda;k,h)-Z_Y(k,h)K_w(k,h)]|
      <=o(sqrt(X)Y^2)||a||_2.                                    (10)

Bounded normalized weight norms can consume this qualitative saving.
An unspecified little-o cannot absorb arbitrary powers of log X, divisor
multiplicities, or powers of X. Restricting to a shorter shift family
retains the full Y^4 bound; no proportional saving in its size is inferred.
The proof averages over X integer starts and treats undilated rectangles.
The actual two-modulus lattices, concentrated starts, class subtraction
and zeta normalization require further estimates. The subsequent
[averaged Euler-tail theorem](rectangle_singular_series_tail.md) replaces
K_w by the full singular series on nondegenerate shifts, with the same
o(XY^4) error and consumption scale. Degenerate shifts retain the finite
model; the four-distinct-offset infinite product does not apply there.
In particular (1)-(2) do not overturn the separate subpower-cutoff
obstruction for the arithmetic-core tail.

The [resolved frame audit](resolved_lock_frame.md) remains valid: its
pointwise-bound criterion loses B^6. The present deduction avoids that
loss for this cube geometry by retaining the other six U^3 norms and
the admissible residue budget. It does not change the published
qualitative theorem into a logarithmically quantitative one.

## Verification

`python scripts/verify_rectangle_transport.py` checks the cube rank
certificate and minors, direct prime-field and CRT admissibility counts,
exact residue/quotient cube splitting, residue-window endpoints and the
finite four-point CRT main term. The analytic uniformity input is the
cited theorem; finite audits do not prove it or its asymptotic rate.
