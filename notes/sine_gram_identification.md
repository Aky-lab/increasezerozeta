# The global sine-process Gram law equals the CUE column law

The finite CUE bandwidth law has a global sine-process realization.
The proof uses finite interaction range, ergodic averages and a uniform
Hilbert--Schmidt tail estimate. Local kernel convergence alone would not
justify this identification.

For the stationary sine process of intensity one, the Gram matrix on a
growing interval has almost surely convergent empirical measures and all
moments. At bandwidth 0<lambda<=1 its limiting law is

    bar_nu_lambda=(1-lambda) delta_0+lambda nu_lambda,                 (1)

where nu_lambda is the CUE row law from
[the bandwidth theorem](cue_gram_bandwidth.md). In particular the unit-
bandwidth sine Gram has the certified continuum moments at every order.
This is a random-matrix model result; arithmetic transport to zeta zeros
and the priority of this application remain separate questions.

## 1. Definition and normalization

Let Xi be the simple stationary determinantal process on R with kernel

    S_1(t)=sin(pi t)/(pi t), S_1(0)=1,
    rho_q(x_1,...,x_q)=det[S_1(x_i-x_j)].

At bandwidth lambda define S_lambda(t)=sin(pi lambda t)/(pi t), with
S_lambda(0)=lambda. Let I_R=[-R/2,R/2], N_R=#(Xi intersect I_R), and

    A_R=[S_lambda(x-y)]_(x,y in Xi intersect I_R),
    mu_R=(1/N_R) sum_j delta_(eigenvalue_j(A_R)).

The empty-window measure may be defined arbitrarily; N_R/R->1 almost
surely, so windows are eventually nonempty. This is a probability measure
per sampled point, not per frequency dimension.

**Theorem.** As R->infinity, mu_R converges weakly almost surely to (1).
For every integer B>=1,

    Tr(A_R^B)/N_R -> lambda m_B(lambda) almost surely,                 (2)
    limsup_(R->infinity) Tr(A_R^B)/N_R <=42^B Bell_(B+1).

The expected measure per unit length

    eta_R=(1/R) E sum_j delta_(eigenvalue_j(A_R))

is a probability measure and converges weakly and in all moments to (1).
For R>=1 its B-th moment is at most 3*42^B Bell_(B+1).

Every nonempty finite A_R is positive definite, yet its limiting zero atom
is exactly 1-lambda. At lambda=1 there is no zero atom. At smaller bandwidth
finite invertibility does not prevent mass from accumulating at zero.

The positive-order moments per frequency dimension are the row moments:
Tr(A_R^B)/(lambda R)->m_B(lambda). A measure formed by dividing every
spectral mass by lambda R instead of N_R has asymptotic total mass
1/lambda; it is not nu_lambda. Equation (1), including its zero atom,
specifies the relation correctly.

An independent second-moment anchor follows directly from the two-point
correlation rho_2(0,t)=1-S_1(t)^2:

    lambda^2+integral S_lambda(t)^2[1-S_1(t)^2] dt
      =lambda+lambda^3/3=lambda*m_2(lambda).

Plancherel gives integral S_lambda^2=lambda. The remaining product integral
is the overlap of two Fourier triangles, namely
2*integral_0^lambda (lambda-u)(1-u)du=lambda^2-lambda^3/3. This check
uses the sine-process correlations directly, rather than the CUE network.

## 2. Background and the ergodic averaging tool

The sine process exists as a translation-invariant determinantal process.
It is mixing, hence its unit-translation action is ergodic; see
[Soshnikov, Determinantal Random Point Fields (2000), Theorem 7](https://arxiv.org/pdf/math/0002099).
The correlation matrices are positive semidefinite with diagonal one, so
Hadamard's inequality gives 0<=rho_q<=1. Campbell's formula then implies
that interval counts X have factorial moments E[(X)_q]<=|I|^q and ordinary
moments of every order. For a unit interval, E X^d<=Bell_d.

We use ordinary Birkhoff averaging under unit translations in this form.
For a translation-covariant root score F(x,Xi) satisfying

    E sum_(x in Xi intersect [0,1)) |F(x,Xi)| <infinity,

bin the sum into unit intervals. The bin scores are an integrable stationary
ergodic sequence. Their two-sided averages converge almost surely to the
expectation of one bin. For nonnegative scores, enclosing an arbitrary
I_R in complete bins gives the corresponding limsup bound. Integrable
end-bin scores divided by R tend to zero: stationarity and
sum_(j>=1) P(|bin_score|>epsilon*j)<infinity prove this by Borel--Cantelli.
The same applies to the negative direction. In particular N_R/R->1.

For local scores used below, discrepancies involving the two window edges
are controlled by counts in intervals of fixed length. Those counts have
moments of every order. Markov with a second moment and Borel--Cantelli
make the edge contribution divided by integer R vanish almost surely.
Enlarging each edge interval by one covers R between adjacent integers.

## 3. Truncation and a spectral distance bound

For K>=1 retain only entries at separation at most K:

    A_(R,K)(x,y)=S_lambda(x-y) 1_(|x-y|<=K).

This matrix is Hermitian; positivity is not required. For Hermitian matrices
of size d, Hoffman--Wielandt and Cauchy--Schwarz give

    d_BL(mu_A,mu_B)<=W_1(mu_A,mu_B)
      <=[Tr((A-B)^2)/d]^(1/2),                                      (3)

where d_BL is the bounded-Lipschitz distance. The matrix inequality is
standard; see [Meckes, The Random Matrix Theory of the Classical Compact
Groups, Lemma 5.22](https://case.edu/artsci/math/mwmeckes/elizabeth/Haar_book.pdf).

Define the infinite root score

    g_K(x,Xi)=sum_(y in Xi, y!=x)
                  |S_lambda(x-y)|^2 1_(|x-y|>K).

Campbell's formula and rho_2<=1 imply

    e_K=E sum_(x in Xi intersect [0,1)) g_K(x,Xi)
       <=integral_(|t|>K) |S_lambda(t)|^2 dt <=2/(pi^2 K).             (4)

This proves the required root-score integrability, without independence.
The finite-window squared Hilbert--Schmidt error is bounded by the sum
of g_K over its roots. The ergodic tool and N_R/R->1 therefore yield

    limsup_R d_BL(mu_R,mu_(R,K))<=sqrt(e_K)<=sqrt(2)/(pi sqrt(K))       (5)

almost surely, simultaneously for positive integer K.

For comparison, scale the CUE angles to a circle of circumference n.
Its intensity is one and its correlation kernel is

    D_(n,n)(t)=(1/n) sum_(a=0)^(n-1) exp(2 pi i a t/n).

The CUE column Gram kernel at m rows is D_(m,n), the same sum ending at
m-1. With t the signed circular distance in [-n/2,n/2],

    |D_(m,n)(t)|<=min(1,1/(2|t|)).

Indeed its absolute value is |sin(pi m t/n)|/[n|sin(pi t/n)|] and
sin(pi |t|/n)>=2|t|/n. CUE correlation determinants also satisfy rho_q<=1.
For the circular-distance truncation K<n/2, (3) gives

    E d_BL(mu_(n,col),mu_(n,col,K))<=1/sqrt(2K).                      (6)

The squared bound before taking the square root is
integral_(K<|t|<=n/2) 1/(4t^2) dt<=1/(2K). It is uniform in m<=n.

## 4. A common limit for each fixed interaction range

Fix K and an integer order B. A trace of the truncated matrix is a sum
over closed walks of B steps, allowing repeated vertices. Group the B
slots by equal vertices; this gives a set partition pi of those slots.
Order its q blocks by first occurrence and let a(i) be the block index
of slot i. The block containing slot zero is the root, whose coordinate
is t_0=0. The limiting anchored moment is

    M_(B,K)=sum_pi integral_(R^(q-1))
         rho_q(0,t_1,...,t_(q-1))
         product_(i=0)^(B-1)
           [S_lambda(t_(a(i))-t_(a(i+1)))
            1_(|t_(a(i))-t_(a(i+1))|<=K)] dt,                        (7)

with cyclic slot indices. For q=1 the integral is just lambda^B.
Distinctness diagonals have measure zero. Each slot walk is connected,
so a nonzero integrand puts every t_j in [-BK,BK]. Also |S_lambda|<=1
and rho_q<=1. Thus these are finite compact integrals.

For the sine process, the infinite rooted truncated trace score has absolute
value at most the (B-1)-st power of the point count within distance BK
of its root. Campbell and count moments imply bin-score integrability.
The ergodic tool, followed by the local edge estimate, proves

    Tr(A_(R,K)^B)/N_R -> M_(B,K) almost surely.                        (8)

The expected trace divided by R has the same limit. More explicitly the
finite-window version of (7) has an additional factor
(1-range(0,t_1,...,t_(q-1))/R)_+, which tends to one on the compact domain.

For CUE, when n>2BK a nonzero rooted circular walk has a unique coherent
lift to the line: total path length is at most BK<n/2, so it cannot wind
around the circle. All its rooted lifted coordinates are in [-BK,BK].
Uniformly on this compact set,

    D_(m,n)(t)->exp(pi i lambda t) S_lambda(t),
    D_(n,n)(t)->exp(pi i t) S_1(t),

as m/n->lambda. The first convergence is a Riemann sum. The phases cancel
in the closed slot product and in the correlation determinant respectively.
The same partition expansion and dominated convergence show

    E[Tr(K_(m,n,K)^B)/n] -> M_(B,K).                                 (9)

This is a finite-range comparison; no unbounded-domain dominated
convergence or interchange with K->infinity has been assumed.

To pass from moments to measures, a uniform growth bound is needed for
the truncated, possibly indefinite matrices. Let d_x count all points
within distance K of x, including x, and let D be the diagonal matrix
of these degrees. Since each retained entry has absolute value at most
one, the inequality 2|u_x u_y|<=|u_x|^2+|u_y|^2 gives

    -D<=A_(R,K)<=D,
    Tr(|A_(R,K)|^(2b))<=2 sum_x d_x^(2b).                            (10)

The second inequality follows by applying eigenvalue monotonicity to
both A and -A: bound their positive eigenvalues separately by those of D.
It also holds for the truncated CUE matrices. The reduced Palm factorial
moments of the other-point count within distance K are at most (2K)^j,
because rho_(j+1)(0,t_1,...,t_j)<=1 and root intensity is one.
Expanding (1+count)^p in nonnegative falling-factorial coefficients gives

    E^0[d_0^p]<=max(1,2K)^p Bell_(p+1).                              (11)

One way to check the constant is to compare these factorial moments
with Poisson(2K); E(1+Poisson(1))^p=Bell_(p+1). Campbell's formula turns
(10)-(11) into uniform even-moment bounds for both expected measures
normalized per length and per n. Carleman's criterion applies: the
(2b)-th moment root grows at most C_K b. Tightness and uniform integrability
then identify a unique probability measure eta_K with moments M_(B,K).
Equation (8), its next-moment bound and determinacy also prove

    mu_(R,K)->eta_K weakly almost surely,
    E mu_(n,col,K)->eta_K weakly.                                    (12)

For the expected sine per-length measures, total mass is exactly one
because E N_R=R. Their convergence to eta_K follows from the same bound.
For empirical sine probability measures, N_R/R->1 converts the pathwise
trace averages in (8) to probability moments. This distinction avoids
replacing an expectation of a random ratio by a ratio of expectations.

## 5. Removing the range cutoff identifies the global law

The [CUE bandwidth theorem](cue_gram_bandwidth.md) already gives
E mu_(n,col)->bar_nu_lambda. For any bounded one-Lipschitz test, take
n->infinity in (6) and (12). Taking the supremum afterwards gives

    d_BL(eta_K,bar_nu_lambda)<=1/sqrt(2K).                            (13)

On the almost sure event for all positive integer K, use (5), (12)
and (13) to obtain

    limsup_R d_BL(mu_R,bar_nu_lambda)
      <=sqrt(2)/(pi sqrt(K))+1/sqrt(2K).

Let K->infinity. This proves the asserted global weak convergence.
For the expected per-length sine measures, (3), (4) and Cauchy--Schwarz
give the analogous bound on their test integrals, with E N_R/R=1.
Hence those expected measures have the same weak limit.

This proves the comparison on the level of weak spectral laws. A central
limit theorem for the growing sine windows does not follow: the cutoff
error above is not controlled at the fluctuation scale.

## 6. Uniform moments for the untruncated sine Gram

Here the positivity of S_lambda is useful. Write

    S_lambda(x-y)=integral_(-lambda/2)^(lambda/2)
                         exp(2 pi i s(x-y)) ds.

The nonzero spectrum of A_R equals that of the finite-rank positive
operator sum_(x in Xi intersect I_R) v_x v_x* on L^2([-lambda/2,lambda/2]),
where v_x(s)=exp(2 pi i s x). The same representation proves finite
positive definiteness: distinct finite exponentials are linearly
independent on every interval of positive length.

For a band-limited function p(t)=integral f(s)exp(2 pi i st) ds,
Parseval gives ||p'||_2<=pi lambda ||p||_2. A set separated by at least
one has disjoint unit intervals centered at its points. Applying the
fundamental theorem of calculus to |p|^2 and then averaging gives

    sum_x |p(x)|^2 <=||p||_2^2+2||p||||p'||
                   <=(1+2 pi lambda)||f||_2^2.

Partition R into unit cells and let X_j be their full-process occupancies.
If each retained cell has at most J points, two colors of adjacent cells
and J point labels split the set into 2J separated sets. Its sampling
operator therefore has norm at most 2(1+2 pi lambda)J<42J.

Delete all columns in cells with X_j>J. The deleted positive operator
has rank at most sum_(cells intersecting I_R) X_j 1_(X_j>J). The min-max
principle bounds the number of eigenvalues above 42J by this rank. Integrating
the spectral tail in intervals [42J,42(J+1)] and telescoping yields

    Tr(A_R^B)<=42^B sum_(cells intersecting I_R) X_j^(B+1).             (14)

This is the line version of the
[square occupancy argument](cue_limit_determinacy.md). Unit-translation
ergodicity, N_R/R->1, and E X_j^(B+1)<=Bell_(B+1) give the pathwise
limsup in the theorem. A next-order bound supplies pathwise uniform
integrability for each fixed power. Weak convergence therefore strengthens
to (2), simultaneously at every integer order.

Taking expectations in (14) and using at most R+2 cells gives the stated
3*42^B Bell_(B+1) bound for R>=1. These bounds similarly upgrade weak
convergence of eta_R to convergence of all its moments. Positivity puts
the limiting measure on [0,infinity), as also follows from (1).

## Scope and computational checks

This proof identifies a precisely normalized global sine-process Gram
law with the established CUE column law. The known zero mass of the CUE
row law gives exactly the mass in (1), even though finite sine matrices
are invertible. Arithmetic transport and a sine-window fluctuation theorem
are not supplied by this argument.

`python scripts/verify_sine_bridge.py` checks finite kernel and phase
identities, the partition/Campbell normalization in a discrete determinantal
model, and exact two-point moment anchors. These tests audit finite
ingredients; the ergodic and limiting assertions are the arguments above.
