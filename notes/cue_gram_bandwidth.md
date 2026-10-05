# Bandwidth-dependent CUE Gram laws and fluctuations

The square Gram model extends to every fixed bandwidth ratio 0<lambda<=1.
The row spectral measure has a determinate almost sure limit, a polynomial
statistic CLT and no atom at zero. The column measure has an atom of mass
exactly 1-lambda. Keeping these two normalizations distinct is essential
when comparing frequency-space and point-space Gram matrices.

The argument extends the [connected network proof](cue_gram_fluctuations.md).
It is a model theorem, without an asserted arithmetic zeta-moment identity
or a claim of priority. See [the literature comparison](literature_map.md).

## 1. Definition and statements

Let z_1,...,z_n be the eigenvalues of a Haar unitary matrix. For 1<=m<=n set

    V_(a,j)=z_j^a/sqrt(n), a=0,...,m-1,
    H_(m,n)=V V*, K_(m,n)=V* V, lambda_n=m/n,
    Y_B(m,n)=Tr(H_(m,n)^B)/m.

H is an m by m principal compression of the square G_n, has diagonal one
and is positive definite almost surely. K has the same nonzero eigenvalues
and n-m exact zero eigenvalues. If mu_row and mu_col denote the empirical
probability measures normalized by m and n respectively, then exactly

    mu_col=(1-lambda_n) delta_0+lambda_n mu_row.                       (1)

**Theorem.** If m/n->lambda in (0,1], there is a unique probability measure
nu_lambda on [0,infinity) such that mu_row converges weakly almost surely
on every coupling of the Haar marginals. All moments converge almost surely
and in L^2. Write m_B(lambda)=integral x^B dnu_lambda. These moments are
locally Lipschitz in lambda on (0,1], and

    E Y_B(m,n)=m_B(lambda_n)+O_B,lambda0(1/n),                         (2)
    m_1(lambda)=1, m_2(lambda)=1+lambda^2/3,
    m_B(lambda)<=42^B Bell_(B+1)/lambda.

For any fixed positive orders B_1,...,B_r with L=sum B_j,

    cum(Y_(B_1),...,Y_(B_r))
      = n^(1-r) [c_(B_1,...,B_r)(lambda_n)+O_(B,lambda0)(1/n)],        (3)

uniformly for lambda_n in [lambda0,1]. The coefficients are continuous
and locally Lipschitz. Finite collections of

    sqrt(n) [Y_B(m,n)-E Y_B(m,n)]

converge jointly to a centered Gaussian vector with covariance
c_(B,C)(lambda). Centering by m_B(lambda_n) gives the same limit. Centering
by m_B(lambda) also works for m=floor(lambda*n); it need not work for
an arbitrary slowly converging sequence m/n->lambda.

For the second moment the nondegenerate limiting variance is explicit:

    c_(2,2)(lambda)=2 lambda^3/15-(2 lambda-1)_+^5/(30 lambda^2)>0.      (4)

At lambda=1 this is 1/10; at lambda=1/2 it is 1/60. The overlap correction
begins at lambda=1/2. This coefficient is four times continuously
differentiable there; its fifth derivative has a jump.

The row limit satisfies

    nu_lambda({0})=0,
    (EulerGamma-1)/lambda <= integral log(x) dnu_lambda <0,
    nu_lambda([0,epsilon])
      <= [(1-EulerGamma)/lambda+lambda/(2 sqrt(3))]/|log epsilon|      (5)

for 0<epsilon<1. Consequently the limiting column measure is
(1-lambda)delta_0+lambda nu_lambda and its zero atom is exactly 1-lambda.
The row laws converge to nu_1 in moments and weakly as lambda->1.
As lambda->0 they converge to delta_1, with

    integral (x-1)^2 dnu_lambda=lambda^2/3.                            (6)

## 2. Rectangular connected network count

Expand the traces using H_(a,b)=Tr(U^(a-b))/n. Outer coordinates now
lie in {0,...,m-1}; the inner Toeplitz shifts still lie in {0,...,n-1}.
The outer cumulant cancellation is unchanged: only slot partitions
connecting all r outer trace cycles survive.

For a fixed surviving partition and fixed anchored cyclic partition inside
each block, use the same lifted incidence equations as in the square proof:

    y_(A,j+1)-y_(A,j)=sum_(i in group(A,j)) (x_next(i)-x_i).

If there are M inner groups, the connected incidence matrix has rank M-1.
There are L+M variables, so the solution space has dimension D=L+1.
The square proof's integral coordinate chart remains valid: contract the
inner cycles, choose a spanning tree of the resulting connected block graph,
solve its tree x variables by incidence elimination, and retain one free
inner offset per slot block. The chart has integer coefficients in both
directions and has D free coordinates. It does not depend on m or n.

Scale the chart by n. Each retained x coordinate lies in [0,lambda_n),
each retained offset lies in [0,1); the derived variables obey the same
bounds. Let Q(lambda) be the closed polytope obtained by imposing
0<=x_i<=lambda and 0<=y_j<=1 on these linear expressions. Volume is
Lebesgue volume in this integer chart, not ambient Euclidean volume.
All x_i=lambda/2 and y_j=1/2 solve the equations strictly, so Q(lambda)
has full dimension for lambda>0.

For this fixed term its lattice count is

    N_(m,n)=n^D vol(Q(lambda_n))+O(n^(D-1)).                           (7)

Here is a direct justification that avoids a two-parameter Ehrhart claim.
The free-coordinate domain is contained in [0,1]^D. Its boundary lies in
finitely many hyperplanes with fixed integer normals and varying intercepts.
Associate each point of the mesh n^(-1)Z^D to its half-open mesh cube.
The discrepancy between their union and Q(lambda_n) lies in slabs of
width C/n around those hyperplanes. Each slab has volume O(1/n), uniformly
in its intercept. The original integer upper bounds m-1 and n-1 differ
from the closed bounds m and n by another such collection of slabs.
This proves (7), with a constant depending only on the term. Zero normal
constraints are vacuous and can be discarded.

Changing lambda by delta changes membership only in analogous slabs of
width C|delta|, so vol(Q(lambda)) is Lipschitz. There are finitely many
terms for fixed orders. Their signed sum, divided by lambda^r, gives
c_(B_1,...,B_r)(lambda). Dividing (7) by the trace normalization
m^r n^L proves (3), since D=L+1. For r=1 it proves (2).
No even-power finite-size correction or polynomial parity is asserted for
the rectangular case; those stronger square identities use centered
reciprocity and do not follow from this mesh argument.

## 3. Determinacy, concentration and the Gaussian limit

Principal-compression interlacing and positivity give

    Tr(H_(m,n)^B)<=Tr(G_n^B),
    E Y_B(m,n)<=lambda_n^(-1) 42^B Bell_(B+1).

The [square occupancy proof](cue_limit_determinacy.md) supplies the last
bound. It gives tightness, uniform integrability of every fixed moment,
and Carleman determinacy, uniformly when lambda_n>=lambda0>0.
Equation (2) gives all limiting moments. Every subsequential expected
measure limit has those moments, hence the expected measures converge to
a unique nu_lambda. The volume continuity and the same moment bounds also
prove weak and moment continuity at lambda=1.

Equation (3) implies Var(Y_B)=O(1/n) and its fourth cumulant is O(1/n^3).
Thus its centered fourth moment is O(1/n^2). Markov's inequality and
Borel--Cantelli prove simultaneous almost sure convergence at all integer
orders, without independence between sizes. Exact first moment one gives
tightness on each such path; the next moment bounds supply uniform
integrability along the path. Determinacy identifies the almost sure
empirical limit. The mean and variance estimates also prove L^2 convergence.

After scaling centered traces by sqrt(n), every joint cumulant of order
r>=3 is O(n^(1-r/2)) and tends to zero. The second cumulants tend to
c_(B,C)(lambda). Moment-cumulant expansion and the determinate Gaussian
moment problem prove the joint CLT, including possibly singular covariance
matrices. The O(1/n) error in (2) justifies the stated alternative centerings.

## 4. Exact variance and the bandwidth threshold

Writing X_k=Tr(U^k), one has for m<=n

    Y_2=1+2/(m n^2) sum_(k=1)^(m-1) (m-k)|X_k|^2,
    E Y_2=1+(m^2-1)/(3n^2).

The [CUE pair covariance formula](cue_gram_fluctuations.md), valid for
1<=k,l<n, is

    Cov(|X_k|^2,|X_l|^2)=1_(k=l) k^2-(k+l-n)_+.

Consequently

    Var(Y_2)=4/(m^2 n^4)
      [sum_(k=1)^(m-1)(m-k)^2 k^2
       -sum_(k,l=1)^(m-1)(m-k)(m-l)(k+l-n)_+]
      =4/(m^2 n^4)[(m^5-m)/30-binomial(2m-n+2,5)_*].                 (8)

The starred binomial is zero when 2m-n<3, and is the ordinary binomial
otherwise. For the second equality put u=m-k,v=m-l,d=2m-n. Since d<=m,
the crossing sum becomes sum_(u,v,w>=1,u+v+w=d) u v w. Its generating
function is t^3/(1-t)^6, giving the stated binomial coefficient. The
diagonal sum is the elementary polynomial (m^5-m)/30.

Taking n times (8) gives (4), because the binomial leading coefficient is
d^5/120. Alternatively its continuum crossing integral is

    integral_(u,v>=0,u+v<d) u v(d-u-v) du dv=d^5/120,
    d=(2 lambda-1)_+.

Since d<=lambda for lambda<=1, c_(2,2)(lambda)>=lambda^3/10>0.
This is a bandwidth-dependent specialization of established CUE pair
covariance theory, rather than a claim of a new pair-covariance technique.

## 5. Logarithmic control survives compression

Almost surely G_n is positive definite. Split it into H_(m,n) and the
remaining principal block. Its Schur complement S is positive definite,
with diagonal entries at most one. The block determinant identity and
Hadamard's inequality imply

    det(G_n)=det(H_(m,n)) det(S)<=det(H_(m,n)).

For m=n the empty Schur complement has determinant one. The square
logarithm is integrable, and log det(G_n)<=log det(H_(m,n))<=0; the
last inequality follows from the unit diagonal and Hadamard's inequality.
Thus the compressed logarithm is integrable as well. The square
[log-determinant identity](cue_gram_logdet.md) therefore gives

    E log det(H_(m,n))/m >= (n/m)(H_n-1-log n).                        (9)

The positive logarithm is bounded by (x-1)_+. The expected row measure
has mean one and second centered moment (m^2-1)/(3n^2). Hence

    E integral log_+(x) dmu_row <= sqrt(m^2-1)/(2 sqrt(3) n).

Combine this with (9) to bound its negative logarithmic integral. Uniform
second moments imply convergence of the positive logarithmic integrals.
For the negative part use bounded continuous truncations min(log_-(x),R),
with value R at x=0, and then monotone convergence. This yields

    integral log_- dnu_lambda
      <= integral log_+ dnu_lambda+(1-EulerGamma)/lambda
      <= lambda/(2 sqrt(3))+(1-EulerGamma)/lambda.

It proves no zero atom and the small-mass bound (5). The lower bound for
the full logarithmic integral follows from the same inequality. Strict
Jensen gives its upper bound <0, since m_2(lambda)>1. Equality with the
finite-size expected logarithm limit is not asserted; negative logarithmic
uniform integrability has not been established.

## 6. Relation to a sine-kernel bandwidth

After a diagonal unitary phase change, the column kernel evaluated at
an angular separation 2 pi t/n is

    sin(pi m t/n)/(n sin(pi t/n)) -> sin(pi lambda t)/(pi t).

Its diagonal tends to lambda. This exact local entry limit alone does not
justify a global spectral identification or an interchange of limits.
[The global comparison proof](sine_gram_identification.md) supplies that
separate argument, using finite-range ergodic averages and uniform
Hilbert--Schmidt tail control. The growing sine-window law equals the
CUE column law in (1), with almost sure convergence of all moments.

Equation (1) specifies the normalization for any such comparison: trace
per sampled point gives column moment lambda*m_B(lambda), whereas trace
per frequency dimension gives row moment m_B(lambda). Both agree at
lambda=1; their zero masses differ when lambda<1. The earlier zeta paper
already proposes a sine-process Gram moment interface. Establishing its
arithmetic hypotheses remains an open task, not a consequence of (1)-(9).

## Reproduction

Run `python scripts/verify_cue_bandwidth.py`. Exact Weyl integration is
performed on the column Gram and compared with the connected trace
expansion on the row Gram at every 1<=m<=n<=5. Separate integer double
sums check (8) for every 1<=m<=n<=40, including the threshold boundary.
The checks cover normalizations and finite identities; the all-order
asymptotic arguments are the proofs above.
