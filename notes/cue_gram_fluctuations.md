# Connected network counts and fluctuations of the Gram spectrum

The Haar-unitary Gram model has almost surely convergent empirical spectral
measures and a joint Gaussian limit for every finite collection of polynomial
spectral statistics. The proof extends the integral network lift from one
trace cycle to several cycles. It does not require independence of Haar
traces or of matrices at different sizes.

These are results about the explicit Gram model. Arithmetic transport to
zeta-zero counting remains a separate problem. The argument is derived here;
priority for this application requires a broader literature comparison and
independent mathematical review.

## Statements

Use [the Gram definition](cue_gram_model.md): U_n is Haar in U(n),
V_(a,j)=z_j^a/sqrt(n), and G_n=V V*. Set

    Y_B(n)=Tr(G_n^B)/n,  B>=1,
    mu_n=(1/n) sum_j delta_(lambda_j(G_n)).

Let m_B be the continuum moments and nu their unique limiting measure from
[the moment-growth proof](cue_limit_determinacy.md). In particular m_1=1.

**Theorem 1 (joint cumulants).** Fix r>=1 and positive orders B_1,...,B_r,
and write L=sum_j B_j. There is a rational polynomial P_(B_1,...,B_r) of
degree at most L+1 such that, for every positive integer n,

    cum(Y_(B_1)(n),...,Y_(B_r)(n)) = P_(B_1,...,B_r)(n)/n^(L+r),          (1)
    P(-n)=(-1)^(L+1) P(n),  P(0)=0.

Consequently

    cum(Y_(B_1),...,Y_(B_r))
      = c_(B_1,...,B_r) n^(1-r) + O(n^(-1-r)),                         (2)

with rational c and (L+1)! c an integer. A coefficient may be zero.
At r=1 this recovers the expected-moment theorem. At r>=2 the polynomial
also vanishes at n=1, since G_1=[1] is deterministic.

**Theorem 2 (strong spectral limit).** At every fixed B,

    Y_B(n) -> m_B in L^2.

On any common probability space carrying the U_n with their Haar marginals,
all these moments converge almost surely, and mu_n converges weakly almost
surely to nu. Without a chosen coupling the convergence statement is weak
convergence in probability. No independence between sizes is required.

**Theorem 3 (polynomial-statistic CLT).** For any fixed orders B_1,...,B_s,

    sqrt(n) (Y_(B_j)(n)-m_(B_j))_(j=1)^s
      -> a centered Gaussian vector with covariance
         Sigma_(i,j)=c_(B_i,B_j).                                      (3)

The matrix Sigma is positive semidefinite and rational; it may be singular.
For example Y_1=1 exactly. The limit extends to any fixed real polynomial
test functions by taking linear combinations of these statistics.

For the second moment there is an explicit nondegenerate example:

    E Y_2(n)=4/3-1/(3n^2),
    Var Y_2(n)=1/(10n)+1/(6n^3)-4/(15n^5),                            (4)
    sqrt(n)(Y_2(n)-4/3) -> N(0,1/10).

For orders two and three the limiting covariance matrix is

    Sigma = [[1/10, 1/3], [1/3, 79/70]],
    det Sigma = 11/6300 > 0.

Thus this two-dimensional example is nondegenerate.
The [separate logarithmic estimate](cue_gram_logdet.md) proves that the
limiting spectral measure has no atom at zero.

## 1. Why only connected slot partitions survive

Put X_k=Tr(U_n^k). Expand each trace as a cycle of B_j outer indices:

    Y_(B_j)=n^(-B_j-1) sum_(x in {0,...,n-1}^(B_j))
              product_(i in cycle j) X_(x_(i+1)-x_i).

Keep slots distinct, including slots with equal frequencies. For a partition
pi of all L slots, say that pi connects the r outer cycles if the following
relation has one equivalence class: two cycles are joined when a block of
pi meets both, and then take the transitive closure.

The cumulant of the r products of slot variables is

    sum_(pi connecting all outer cycles)
       product_(A in pi) cum(X_(k_i):i in A).                           (5)

Here is a direct proof of this product-cumulant identity. The outer cumulant
expands as a sum over partitions rho of the r cycle indices, with coefficient
(-1)^(|rho|-1)(|rho|-1)!, followed by the product of joint moments on its
blocks. Expand those moments into slot cumulants. For a fixed slot partition
pi, let sigma be its induced cycle-component partition. This term occurs
precisely for rho coarser than sigma. Its coefficient is therefore

    sum_(rho coarser than sigma) (-1)^(|rho|-1)(|rho|-1)!.

If sigma has t blocks, this is the Stirling sum
sum_(a=1)^t {t brace a}(-1)^(a-1)(a-1)!, equal to 1 at t=1 and zero
at t>1. Extracting coefficients from log(exp(z))=z proves the identity.
Thus disconnected contributions cancel exactly, rather than being estimated.

The [Toeplitz log-determinant derivation](cue_gram_model.md) evaluates each
slot cumulant as the signed sum over cyclically ordered partitions of A.
A term is zero unless sum_(i in A) k_i=0. When this sum is zero its value is
(-1)^(number of inner groups-1) max(n-range(prefixes),0). This formula includes
zero frequencies and singleton blocks, so (5) needs no genericity condition.

## 2. The connected multi-cycle lift

Fix a connecting slot partition pi and one cyclic-partition term inside each
slot cumulant. Let M be the total number of inner groups. Introduce L outer
coordinates x_i and M inner offset coordinates y_(A,j), all in [0,n-1].
The overlap counts are exactly the integer solutions of

    y_(A,j+1)-y_(A,j)
      = sum_(i in inner group (A,j)) (x_(i+1)-x_i).                     (6)

Successors of x stay within their own outer cycle; successors of y stay
within their own inner cycle. Summing (6) around an inner cycle imposes
the requisite closure sum_(i in A) k_i=0. Conversely, once it holds, the
allowed starting y values are precisely the overlap count. Thus the lift
introduces neither a missing constraint nor an extra factor.

The matrix of (6) is a node-arc incidence matrix on M inner-group nodes.
Each outer x coordinate supplies an arc between the inner owners of adjacent
slots; each inner y coordinate supplies an arc within its inner cycle.
Loops give zero columns. Within one outer cycle its x arcs connect all nodes
that it visits. Within a slot block its y arcs connect all its inner groups.
Since pi connects all outer cycles, the resulting graph is connected.
Its incidence rank is M-1. Consequently the solution subspace has dimension

    (L+M)-(M-1)=L+1.                                                   (7)

This is the dimension gain: r disconnected trace cycles would have dimension
L+r, but a joint cumulant retains only the dimension L+1 terms.

The incidence matrix is totally unimodular, including loops. Intersecting its
kernel with the unit cube yields an integral polytope Q in the lattice
ker(A) intersect Z^(L+M). The half-ones vector satisfies (6) and lies strictly
inside all cube bounds, so Q has the full relative dimension (7). For positive
integer n its count is the Ehrhart polynomial E_Q(n-1), of degree L+1.

This is the same integer lattice as the original outer/offset sum. To make
the normalization explicit, contract each inner y cycle to its slot block.
The remaining coarse outer incidence graph is connected. Delete one redundant
closure equation and choose a coarse spanning tree. Its reduced incidence
minor has determinant +/-1; solve for those outer x coordinates integrally.
The L-|pi|+1 other outer coordinates, together with one starting inner offset
per slot block, form an integer bijective chart of dimension L+1. No
Euclidean covolume enters the volume coefficient.

Summing the finite set of lifted terms with their integer signs gives the
polynomial numerator P in (1). The coefficient of n^(L+1) is a signed sum of
normalized lattice volumes. For integral polytopes (L+1)! times each volume
is an integer, proving the rationality and denominator assertion.

Every lifted equation also satisfies A*1=0. The interior-shift argument of
[centered reciprocity](centered_reciprocity.md) applies verbatim: subtracting
the all-ones vector identifies interior integer points of (n+1)Q with
points of (n-1)Q. Ehrhart reciprocity gives

    E_Q(-n-1)=(-1)^(L+1) E_Q(n-1).

There are no interior integer points in Q, so E_Q(-1)=0. Therefore P has
the parity and zero in Theorem 1. Its possible powers are L+1,L-1,...,
which gives (2) after dividing by n^(L+r).

## 3. Concentration and the strong spectral limit

At fixed B, (2) with r=2 gives Var(Y_B)=O(1/n). At r=1, the expected-moment
identity gives E Y_B=m_B+O(1/n^2). This proves L^2 convergence.

Let Z_B=Y_B-E Y_B. The fourth central moment is exactly

    E Z_B^4 = cum(Y_B,Y_B,Y_B,Y_B)+3 Var(Y_B)^2 = O(1/n^2).

The r=4 term is O(1/n^3) by (2). Markov's inequality bounds the probability
of |Z_B|>eta by O_B(eta^(-4)n^(-2)). This is summable over n. The first
Borel-Cantelli lemma yields Z_B->0 almost surely for any coupling. Taking
the intersection over the countably many orders B proves simultaneous
almost sure moment convergence.

For every realization, mu_n is a probability measure on [0,infinity) with
first moment exactly one, hence mu_n([R,infinity))<=1/R. It is tight. On
the event of simultaneous moment convergence, each next moment Y_(B+1)
is eventually bounded. This makes x^B uniformly integrable along that
realization: its integral over x>R is at most Y_(B+1)/R. Every weak
subsequential limit therefore has moments m_B at all orders. The previously
proved moment determinacy identifies it as nu. All subsequences have the
same limit, proving weak almost sure convergence.

## 4. Gaussian fluctuations

For any fixed real linear combination of centered statistics, define

    W_n=sqrt(n) sum_j a_j (Y_(B_j)-E Y_(B_j)).

Its first cumulant is zero, its second converges to a^T Sigma a, and its
r-th cumulant for fixed r>=3 is O(n^(1-r/2)), by multilinearity and (2).
It tends to zero. Expressing every fixed moment as a finite sum of products
of cumulants leaves exactly the pair partitions in the limit. These are
the moments of a centered Gaussian with variance a^T Sigma a. Bounded
second moments give tightness; bounded higher even moments allow passage
of moments along weak subsequences. The Gaussian moment sequence is
determinate, so the method of moments proves convergence. If the limiting
variance is zero, convergence to zero follows already in L^2.

The covariance matrix is positive semidefinite as a limit of n times the
finite covariance matrices. Cramer-Wold gives (3). Finally,
sqrt(n)(E Y_B-m_B)=O(n^(-3/2)), so centering by m_B gives the same limit.

## 5. An exact variance certificate

Let P_(2,2)(n)=n^6 Var(Y_2(n)). Theorem 1 gives an odd polynomial of degree
at most five, with zeros at n=0,1,-1. Hence it has the form

    P_(2,2)(n)=n(n^2-1)(a n^2+b).

Independent Weyl integration at n=2,3 yields P(2)=4 and P(3)=28. These two
values determine a=1/10 and b=4/15, so

    P_(2,2)(n)=n(n^2-1)(3n^2+8)/30.

Dividing by n^6 proves (4) for every n>=1. Weyl integration at n=4 gives
the held-out value P(4)=112.

There is a separate trace-cumulant check. For 1<=k,l<n, direct collection
of the 26 cyclic-partition terms gives

    cum(X_k,X_-k,X_l,X_-l)=-(k+l-n)_+.

On the cone 0<=k<=l the net signed range multiplicities are

| Prefix range | Net coefficient |
|---|---:|
| l-k | -1 |
| l | 2 |
| k+l | -1 |

Equal ranges coalesce on the cone boundary; symmetry covers l<=k. Since
k,l<n, their signed overlap sum is
2(n-l)-(n-l+k)-(n-k-l)_+=-(k+l-n)_+. The covariance of the squared traces is

    Cov(|X_k|^2,|X_l|^2)=1_(k=l) k^2-(k+l-n)_+.

The Toeplitz expression

    Y_2=1+(2/n^3) sum_(k=1)^(n-1) (n-k)|X_k|^2

therefore gives the finite double-sum check

    n^6 Var(Y_2)
      =4 sum_(k,l=1)^(n-1) (n-k)(n-l)
           [1_(k=l) k^2-(k+l-n)_+].

### Further joint certificates

The same parity and n=1 zero imply P=n^a(n^2-1)R(n^2), with a=1 for
even L and a=2 for odd L. The number of coefficients is floor(L/2).
Exact Weyl values at n=2,...,floor(L/2)+1 therefore determine P before
any held-out check. The additional numerator polynomials are

| Orders | P(n) in (1) | Leading c |
|---|---|---:|
| (2,3) | n^2(n^2-1)(n^2+2)/3 | 1/3 |
| (3,3) | n(n^2-1)(79n^4+107n^2-12)/70 | 79/70 |
| (2,2,2) | n(n^2-1)(2n^4-5n^2+48)/45 | 2/45 |

Their fit-size numerator values at n=2,3,4 are respectively (24,264,1440),
(144,2520,18792), and (8,88,640). The (2,3) row needs only the first two;
the other rows need all three. Independent Weyl integration at n=5 is
held out for every row. This proves the listed exact identities for all n,
using the already established polynomial space rather than an empirical fit.

## Prior results and attribution

The trace-cumulant method and CUE pair-statistic variance formula are
established background. [Soshnikov and Wu, *A Note on Cumulant Technique in
Random Matrix Theory* (2023), Proposition 2 and equations (50)–(51)](https://escholarship.org/content/qt1888k21h/qt1888k21h_noSplash_8e205173bf8abb441168126ac4d9cd1a.pdf)
give the pair variance for general square-integrable even test functions.
Y_2 is the pair statistic with test function F_n/n^2, where F_n is the
Fejer kernel. Its Fourier coefficients are (n-|k|)/n^3 for |k|<n and zero
otherwise. Substitution into that published variance formula also gives
the double sum above. We do not claim priority for the pair covariance
identity or the cumulant method.

The derivation here adds the integral connected multi-cycle lift, its
finite rational correction polynomials, and the resulting all-order
statements for polynomial statistics of this Gram matrix. A complete
comparison with earlier Vandermonde, point-process and higher-body
statistic results is still needed before a novelty claim.

## Independent verification

```sh
python scripts/verify_cue_fluctuations.py
```

The standard-library verifier computes joint moments directly by integer
Laurent constant terms in the Weyl density. It uses the column Gram V*V,
whose trace powers equal those of V V*, and does not use connected network
counts to obtain those values. A separate connected-partition evaluator
checks joint cumulants at small sizes. It also verifies the cycle-component
cancellation and incidence connectivity, the four-cumulant range table on
the entire rational cone, and the exact variance formula at additional sizes.

The checks support the finite identities and variance certificate. The
all-order concentration and CLT follow from the dimension and cumulant
arguments above, rather than a fit to simulated fluctuations. A more
general nonpolynomial-statistic CLT and the arithmetic counting interface
remain further research questions.
