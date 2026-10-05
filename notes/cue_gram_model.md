# A Haar-unitary Gram realization of the continuum model

The certified continuum moments are limits of expected spectral moments
of explicit positive matrices. This supplies a probabilistic foundation
for the finite model and an independent way to check its assembly.
It does not prove arithmetic transport to the Riemann zeta function.

## Matrix definition and exact finite-size identity

Let U be Haar distributed in U(n), with eigenvalues z_j=exp(i theta_j).
Define the square Vandermonde matrix and its Gram matrix by

    V_(a,j) = z_j^a / sqrt(n),  a=0,...,n-1, j=1,...,n,
    G_n = V V*.

G_n is positive semidefinite, its diagonal entries are one, and its
normalized trace is exactly one. Its entries are

    (G_n)_(a,b) = Tr(U^(a-b))/n.

Writing X_k=Tr(U^k), cyclic expansion gives, for every fixed B>=1,

    M_B(n) = E[Tr(G_n^B)/n]
           = n^(-B-1) sum_(x in {0,...,n-1}^B)
             E product_(i=0)^(B-1) X_(x_(i+1)-x_i).

Reversing the cycle accommodates the sign in the matrix-entry formula.
The expectation is invariant under that reversal. All indices are cyclic.

Let P range over set partitions of the B labels. Moment-cumulant expansion
and the trace-cumulant formula below yield the exact identity

    n^(B+1) M_B(n) = sum_P S_P(n),

where S_P is precisely the signed network count, permitting singleton
blocks as well. Each block must have total frequency zero. The inner count
for its cyclic partition pi is max(n-range(prefixes(pi)),0), with coefficient
(-1)^(number_of_blocks(pi)-1).

## Derivation of the cumulants

For formal variables t_i and f(z)=sum_i t_i z^(k_i), Weyl integration
and determinant integration give the finite-dimensional identity

    E exp(sum_i t_i X_(k_i)) = det T_n(exp(f)),

where T_n(h)_(a,b) is the Fourier coefficient of h at a-b and
a,b=0,...,n-1. Since T_n(1)=I, expand

    log det T_n(exp(f))
      = sum_(m>=1) (-1)^(m-1)/m Tr[T_n(exp(f)-1)^m].

To extract the coefficient of a product of distinct t_i, assign every
label to one of the m nonempty factors. Within each factor the coefficient
is z raised to the sum of its frequencies. Ordered nonempty blocks result.
Cyclic invariance of the trace cancels the factor 1/m: equivalently, use
set partitions and cyclic orders anchored at the block containing label zero.
This is exactly the inner compiler's convention.

Each T_n(z^k) is a truncated integer shift. The trace of the product of m
such shifts vanishes unless their sum is zero. If the sum is zero, its
trace counts starting indices such that every cumulative shifted index
lies in {0,...,n-1}. The count is

    max(n - (max(prefixes)-min(prefixes)), 0).

Thus the joint trace cumulant equals the repository's K_b(n,k), including
its sign and lattice normalization. This proves the finite-size identity
without assuming Gaussian independence of the traces.

A singleton has cumulant n when its frequency is zero and zero otherwise.
Deleting its zero increment leaves the outer overlap unchanged and contributes
one factor n. The [local three-block identity](model_local_identities.md)
vanishes on the bounded outer support. Therefore all moments through eighth
order use exactly the fourteen non-singleton, non-three-block signatures
in the reduced certificate, with singleton multiplicities binom(B,k).

## Consequences at every fixed order

For every partition and every term, the network lift is an integral
polytope of dimension B+1, including partitions with singleton blocks.
Consequently T_B(n)=n^(B+1) M_B(n) is a rational polynomial of degree at most
B+1. [Centered reciprocity](centered_reciprocity.md) gives

    T_B(-n)=(-1)^(B+1) T_B(n),
    T_B(0)=0.

Thus the expected Gram moments have only even inverse-size corrections:

    M_B(n) = m_B + a_(B,1)/n^2 + a_(B,2)/n^4 + ... .

The leading coefficients m_B are exactly the continuum Bell-class moments.
They are rational, and (B+1)! m_B is an integer. For example the complete
eighth-order finite-size identity is

    M_8(n) = 747361/20160
             - 916103/(10080 n^2)
             + 46891/(576 n^4)
             - 147341/(5040 n^6)
             + 277/(105 n^8).

There exists a probability measure on [0,infinity) with all these m_B as
moments. To see this, let mu_n be the expected empirical spectral measure
of G_n. It is a probability measure with first moment one, hence is tight.
At every fixed order its moments are bounded uniformly in n by the finite
inverse-size formula. Along a weakly convergent subsequence, the bound on
the moment of order B+1 makes x^B uniformly integrable; truncation then
gives convergence of the B-th moment to m_B. The limit has the claimed
moments at every order. This argument establishes existence. The separate
[occupancy and moment-growth proof](cue_limit_determinacy.md) establishes
moment determinacy and weak convergence of the entire sequence of expected
measures. It bounds M_B(n) by 42^B Bell_(B+1), uniformly in n.

In particular, the finite model's moment and shifted-moment matrices are
positive semidefinite on mathematical grounds, independently of the finite
Hankel inversions. For its unique limiting measure nu, the computed degree-four
Christoffel polynomial gives the valid model bound

    nu({0}) <= 12241115/162540559.

The separate [log-determinant argument](cue_gram_logdet.md) strengthens this
all-order model conclusion to nu({0})=0, with uniform control of mass near
zero. The displayed Christoffel bound remains an exact finite-moment
certificate.

Transferring this measure or bound to a zeta-zero counting problem still
requires the separate analytic interface.

## Independent exact checks

[The verifier](../scripts/verify_cue_model.py) computes Gram moments directly
from Weyl's eigenangle density

    |det(z_j^a)|^2 / n!,

relative to product uniform circle measure. It expands the integer Laurent
polynomial Tr((n G_n)^B), multiplies by the squared Vandermonde, and extracts
the constant coefficient. It does not use cumulant compilation, network
coordinates, frequency sorting or dihedral representatives to compute the
independent values. Rotational invariance allows z_n=1, leaving n-1 variables;
every term is homogeneous of total degree zero, so this preserves the constant
coefficient. Normalization is n! n^(B+1).

At n=2 there is also a closed check. If delta is the angle difference,
its density is (1-cos(delta))/(2 pi), and the Gram eigenvalues are
1 +/- |cos(delta/2)|. Hence

    M_B(2) = sum_(j=0)^(floor(B/2)) binom(B,2j) Catalan_j / 4^j.

All moments B=0,...,8 agree exactly at n=1,2,3,4. At eighth order the independent
values are M_8(2)=2431/128, M_8(3)=549917/19683 and M_8(4)=1038817/32768. The verifier also checks
the parity of every assembled polynomial and its leading model moment.

```sh
python scripts/verify_cue_model.py
```

This requires only the Python standard library and runs in about a second
in the recorded environment. Outputs are included in the
[project verification record](../results/project_verification_2026-10-05.json).

## Attribution and next research

The CUE eigenangle density, determinant functional and trace-cumulant
technique are established background. See Soshnikov and Wu,
[A Note on Cumulant Technique in Random Matrix Theory (2023)](https://escholarship.org/content/qt1888k21h/qt1888k21h_noSplash_8e205173bf8abb441168126ac4d9cd1a.pdf),
especially the discussion of CUE cumulants, and Bourgade's
[thesis, section 4.1](https://math.nyu.edu/~bourgade/papers/PhDThesis.pdf)
for Weyl integration. The finite Gram identity, network interpretation and
explicit correction polynomials are derived here from this background;
priority for this specific combination remains to be assessed.

The [connected multi-cycle proof](cue_gram_fluctuations.md) establishes
almost sure convergence of the empirical spectrum and Gaussian fluctuations
of fixed polynomial statistics. Further questions concern wider fluctuation
test classes, efficient evaluation of later moments and covariance coefficients,
and the uniform estimates connecting the Gram formulation to arithmetic
zeta moments. The counting interface remains a separate obligation.
