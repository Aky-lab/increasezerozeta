# Pure-cycle integrals from integral flow polytopes

The continuum model integrals through order seven are

| Order b | Exact C_b |
|---|---:|
| 4 | -1/60 |
| 5 | 1/36 |
| 6 | -1/126 |
| 7 | -17/360 |

The order-seven value is obtained by integer counts and a proved
polynomial-degree bound. The lower-order results recover the pinned
reference values. These are continuum model evaluations; their arithmetic
transport to zeta moments is a separate analytic problem.

## 1. Definition

Let c0,...,c_(b-1) sum to zero, with b-1 free frequencies, and put
s0=0 and si=c0+...+c_(i-1). For an ordered cyclic partition into
m blocks B0,...,B_(m-1), put q0=0 and

    qj = sum_(i in B0 union ... union B_(j-1)) ci.

Its contribution is

    V_B = integral ov(s0,...,s_(b-1)) ov(q0,...,q_(m-1)) dc,
    ov(p) = (1-max(p)+min(p))_+.

The pure-cycle integral is the sum of V_B with sign (-1)^(m-1),
over cyclic block orders. The block containing label zero anchors
each order. There are 26, 150, 1,082 and 9,366 terms at orders
four, five, six and seven respectively.

These definitions match the [reference implementation](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/gpu/exact_polytope.py),
revision d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8. The cumulant
normalization is unchanged.

## 2. The lifted polytope is a bounded circulation polytope

Lift both overlaps using t<=si<=t+1 and u<=qj<=u+1.
Define xi=si-t and yj=qj-u. Then

    0<=xi<=1,  0<=yj<=1,
    ci=x_(i+1)-xi,
    y_(j+1)-yj = sum_(i in Bj) (x_(i+1)-xi),

with both sets of indices interpreted cyclically. The original
variables (c0,...,c_(b-2),t,u) and free variables
(x0,...,x_(b-1),y0) are related by an integer linear map with an
integer inverse. Its determinant has absolute value one, so it
preserves both lattice points and the integration measure.

Write the block equations as A(x,y)=0. If owner(i) is the block
containing label i, column xi has +1 at owner(i) and -1 at
owner(i-1). Column yj has +1 at block j-1 and -1 at block j.
An edge with the same source and destination gives a zero column.
Thus A is a directed node-arc incidence matrix. The y columns
contain a spanning cycle, so its rank is m-1, including rank zero
when m=1. The affine dimension is b+m-(m-1)=b+1.

For completeness, incidence matrices are totally unimodular:
in a square submatrix, a column with at most one nonzero permits
determinant expansion and induction; if every column has two,
the opposite signs give zero column sums and dependent rows.
Every square determinant is consequently 0, +1 or -1. Appending
unit bound rows and replacing equalities by their two inequalities
preserves this property. Integer right hand sides therefore give
integer vertices. This is the standard network-polyhedron
integrality argument; see [Goemans, Linear Programming and
Polyhedral Combinatorics](https://math.mit.edu/~goemans/18433S11/polyhedral.pdf).

The free-coordinate projection is bijective, with integer inverse
yj=y0+qj(x). Hence each lifted term is an integral polytope of
dimension b+1 in the explicit free-coordinate lattice. This also
fixes the volume normalization: no Euclidean covolume factor is
introduced by the redundant conservation equations.

## 3. Exact lattice counts

Dilate by an integer N>=0. For a fixed x in {0,...,N}^b,
the admissible integer y0 lie in

    -min(q) <= y0 <= N-max(q).

Their number is max(N+1-range(q),0). Put n=N+1. Summing over
terms with signs gives

    S_b(n) = sum_(x in {0,...,n-1}^b)
             sum_B (-1)^(m-1) (n-range(q_B(c(x))))_+.

This is the signed sum of the Ehrhart counts E_B(n-1).
The Ehrhart theorem makes each E_B a polynomial of degree b+1,
with leading coefficient V_B, because the polytope is integral.
Consequently S_b(n) is a polynomial of degree at most b+1 for
all positive integers n, and its leading coefficient is C_b.
The leading-coefficient statement is also the unweighted case of
[Baldoni et al., weighted Ehrhart quasi-polynomials](https://arxiv.org/abs/1011.1602).

No period estimate or rational reconstruction is needed. At order
seven, nine exact counts n=1,...,9 give

    C7 = eighth_forward_difference(S7) / 8! = -17/360.

The additional n=10 count is held out and equals the polynomial
prediction. The resulting polynomial is

    S7(n) = -n^2 (n^2-1) (n^2-4) (17n^2-27) / 360.

This expression is recovered from the counts; it is not supplied
as an input to the lattice evaluator.

## 4. Implementation and checks

`scripts/pure_cycle_lattice.py` uses int32 for cumulant values,
int64 for summation, and Python integers/Fractions for differences.
At its supported b<=7 and n<=12,

    absolute sum <= term_count * n^(b+1) < 2^63,
    absolute point value <= term_count * n < 2^31.

Every frequency satisfies |ci|<=n-1. The cumulant is invariant
under frequency permutations: relabel its set partitions and
cyclic orders. Choosing a different cyclic anchor only translates
the prefix positions, preserving their range. The implementation
therefore evaluates it once per sorted zero-sum frequency tuple,
then uses exact integer keys to retrieve those values.

The structural gates check all term incidence columns and the
eliminated flow equations. A separate permutation-and-cut
enumerator checks the cyclic-partition compiler through order seven.
Cached sums agree with direct full-cube sums, including the
nontrivial order-seven n=4 count. All 150 reference C5 term values
and twelve C6 term values are reproduced using direct block sums,
without frequency sorting or compiled coefficient forms. Six small
bounded-flow enumerations check the overlap-count interpretation.
The additional order-seven n=11 count agrees with the polynomial
prediction. The [independent verification record](../results/pure_cycle_independent_checks_2026-10-05.json)
and [reference term fixtures](../results/pure_cycle_reference_anchors.json)
identify the exact comparisons and pinned sources.
The recorded calculation uses Python 3.12.12 and NumPy 2.3.5.

Reproduce all four orders with:

    python scripts/pure_cycle_lattice.py --out results/pure_cycle_reproduction.json

Check the stored integer differences, held-out counts, structural
gates and source hash without repeating the larger cubes:

    python scripts/pure_cycle_lattice.py --out results/pure_cycle_flow_2026-10-05.json --verify-only

The [raw record](../results/pure_cycle_flow_2026-10-05.json) contains
the signed integer counts, runtime bounds and certificates.

Run the independent checks with NumPy:

    python scripts/verify_pure_cycle_lattice.py --out results/pure_cycle_checks.json

The standard-library `scripts/verify_model_certificates.py` checks
the stored integer arithmetic, held-out predictions and source
hashes for both the pure-cycle and {5,2} certificates. It is part
of the baseline verification workflow and does not replace full
lattice enumeration.

## 5. Seventh-order consequence

Combining C7=-17/360 with the exact corrected {5,2}=1/8 gives

    1717/90 + {5,2} + C7 = 862/45.

This closes those two finite model evaluations. The identification
of that ledger with an arithmetic seventh moment still requires
the frozen-singleton, vanishing and transport statements.
