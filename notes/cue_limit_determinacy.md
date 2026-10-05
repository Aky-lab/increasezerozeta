# Moment determinacy and convergence of the expected Gram spectrum

For the Haar-unitary Vandermonde Gram matrices G_n defined in
[the model realization](cue_gram_model.md), the expected empirical spectral
measures converge weakly to a unique probability measure on [0,infinity).
Its moments are the continuum Bell-class moments at every order.

The key estimate is, for every integer B>=1 and every n,

    E[Tr(G_n^B)/n] <= 42^B Bell_(B+1).

Here Bell_d is the number of set partitions of d labels. The constant 42
is deliberately coarse. The proof combines an elementary sampling inequality,
a rank perturbation, and determinantal occupancy bounds. It establishes
convergence of expected measures; fluctuations of the random empirical
measures require an additional argument.

## 1. Sampling a separated set

Use circle coordinates t in R/Z. For a trigonometric polynomial

    p(t)=sum_(a=0)^(n-1) v_a exp(2 pi i a t),

Parseval gives ||p'||_2 <= 2 pi (n-1) ||p||_2. Suppose a set of sample
points has circular separation at least 1/n. Intervals of length 1/n
centered at these points have disjoint interiors. For g=|p|^2 and one
such interval I centered at t_j, the fundamental theorem of calculus,
followed by averaging over I, gives

    g(t_j) <= n integral_I g + integral_I |g'|.

Sum over the intervals and use |g'|<=2|p'||p| and Cauchy--Schwarz:

    sum_j |p(t_j)|^2
      <= [n+4 pi (n-1)] ||p||_2^2.

Thus the sum of the normalized outer products of the Vandermonde columns
at any such separated set has operator norm at most 1+4 pi.

Partition the circle into n half-open intervals of length 1/n. If every
interval contains at most K points, divide the intervals into at most
three color classes so that adjacent intervals on the circle have different
colors. This is a proper coloring of a cycle; for n=2 two colors suffice.
Number the points within each interval from 1 to K. A fixed color and
point number give a separated set. There are at most 3K such sets, hence

    ||G_n|| <= 3(1+4 pi) K < 42 K.

The case n=1 is immediate. The same inequality holds for the sum formed
from any subset of the columns; the normalization remains 1/n.

## 2. Remove crowded intervals

Let X_j be the number of eigenangles in interval j. For an integer K>=0,
remove all points in intervals with X_j>K. The matrix of the remaining
columns has norm at most 42K (and is zero for K=0). The removed columns
give a positive matrix of rank at most

    H_K = sum_j X_j 1_(X_j>K).

By the min-max principle, G_n has at most H_K eigenvalues above 42K.
Writing mu_n for its expected empirical spectral measure, stationarity
of the CUE eigenangles gives

    mu_n((42K,infinity)) <= E[X 1_(X>K)],

where X is the occupancy of any one interval of length 1/n. The factor n
from summing the n identically distributed occupancies cancels the spectral
normalization 1/n. This cancellation avoids using the maximum occupancy,
which would give a bound deteriorating with n.

## 3. Determinantal occupancy moments

With respect to uniform measure dt on R/Z, the CUE kernel is

    K_n(t,s)=sum_(a=0)^(n-1) exp(2 pi i a(t-s)),
    K_n(t,t)=n.

Its k-point correlation function is det(K_n(t_i,t_j)). The matrix is
positive semidefinite. Hadamard's inequality bounds the determinant by n^k.
Consequently the falling-factorial moments of X satisfy

    E[(X)_k] = integral_(I^k) det(K_n(t_i,t_j)) dt_1...dt_k <= 1.

For k>n the factorial moment is zero. Converting factorial moments to
ordinary moments with nonnegative Stirling numbers gives

    E[X^d] = sum_(k=0)^d {d brace k} E[(X)_k] <= Bell_d,  d>=1.

This follows directly from the correlation function; it does not assume
that the occupancies in different intervals are independent.

## 4. Bound the spectral moments

Integrate the spectral tail in successive intervals [42K,42(K+1)].
For B>=1, positivity and the previous rank estimate give

    M_B(n) = B integral_0^infinity x^(B-1) mu_n((x,infinity)) dx
      <= 42^B sum_(K>=0) [(K+1)^B-K^B] E[X 1_(X>K)]
      = 42^B E[X * sum_(K=0)^(X-1) ((K+1)^B-K^B)]
      = 42^B E[X^(B+1)]
      <= 42^B Bell_(B+1).

All summands are nonnegative, so exchanging expectation and summation is
justified. The finite-size network formula proves M_B(n)->m_B at each
fixed B. Taking the limit preserves this bound for m_B.

## 5. Carleman's criterion and the limit

There is an injection from partitions of d labels into functions from
d labels to themselves: send each label to the smallest label in its block.
Thus Bell_d<=d^d. In particular,

    m_(2k)^(1/(2k)) <= 42 (2k+1)^(1+1/(2k)) <= 126 (2k+1).

The final inequality uses 2k+1<=3^(2k), for k>=1. Hence

    sum_(k>=1) m_(2k)^(-1/(2k)) >=
      (1/126) sum_(k>=1) 1/(2k+1) = infinity.

If an even moment is zero the criterion holds directly. Carleman's criterion
therefore makes the moment problem determinate even among measures on R.

The expected measures are tight because their first moment is exactly one.
Along any weakly convergent subsequence, the uniform moment bound at order
B+1 gives uniform integrability at order B. Every subsequential limit has
all moments m_B and is supported on [0,infinity). Existence follows by
tightness; determinacy makes every limit the same measure. Thus the entire
sequence mu_n converges weakly to this unique measure, with convergence
of moments at every order.

## Attribution and scope

The CUE kernel and determinantal correlation functions are standard; see
Elizabeth Meckes, [The Random Matrix Theory of the Classical Compact Groups,
section 3.2, Proposition 3.7](https://case.edu/artsci/math/mwmeckes/elizabeth/Haar_book.pdf).
The moment criterion is classical; see Gwo Dong Lin,
[Recent Developments on the Moment Problem (2017), Theorem 1](https://arxiv.org/pdf/1703.01027).
The occupancy-to-Gram bound and its application to this continuum model
are developed here. This is a proof for the stated model; priority for
this application has not been established.

The theorem resolves moment determinacy for the full model sequence,
without needing higher-order class enumeration. It leaves open almost-sure
spectral convergence, sharper spectral tails, evaluation of later moments,
and the arithmetic-to-model transport needed for a zeta-zero bound.
