# Exact mixed-cycle integrals

For a partition P of B cyclic positions into r blocks, assign frequencies
c_i with sum_(i in A) c_i=0 for each block A. One frequency per block is
eliminated. Let s be the full outer prefix walk. The class integral is

    I_P = integral ov(s) product_(A in P) C_|A|(c_A) dc,

where C_b is the partition-cyclic cumulant in the
[pure-cycle convention](pure_cycle_flow_polytopes.md). The free-frequency
dimension is B-r. A signature denotes the sum over all its placements.

## Integral network lift

Expand each inner cumulant into its signed cyclic-partition terms.
Let m_A be the number of inner blocks chosen in outer block A and
M=sum_A m_A. Lift the outer overlap with B variables x_i=s_i-t in
[0,1], and each inner overlap with m_A variables y_(A,j) in [0,1].
The block equations are

    y_(A,j+1)-y_(A,j)
       = sum_(i in inner block (A,j)) (x_(i+1)-x_i).

All indices in each inner cycle and in the outer cycle are cyclic.
Every x column has one +1 and one -1 at the inner blocks owning its
two adjacent labels. Every y column is an edge of an inner cycle.
Zero columns correspond to loops. The outer cycle visits every inner
block, so the full incidence graph is connected and has rank M-1.

Thus the bounded polytope has dimension

    B+M-(M-1) = B+1.

The incidence matrix is totally unimodular, so integer bounds give
integral vertices. This is the same network-integrality argument used
for the pure-cycle and pair-cycle certificates.

To fix the lattice and volume normalization, first sum the inner
equations around each outer block. This yields the r outer closure
equations sum_(i in A)(x_(i+1)-x_i)=0, of rank r-1. Eliminate the
x columns of a spanning tree in this reduced incidence graph. The
remaining B-r+1 x coordinates freely determine all eliminated x
coordinates by integer linear combinations.

For each outer block, one chosen inner y coordinate determines all
its other y coordinates by integer block sums. These r inner offsets
and the B-r+1 free outer coordinates form a bijection with the full
affine integer lattice. They also have an integer bijection with the
original free frequencies, t and the r inner offsets: integrate the
outer increments and use x_0=-t. The coordinate transformation is
unimodular. No Euclidean surface-measure factor is introduced.

## Integer counts

Dilate by N and set n=N+1. For a given bounded outer x satisfying
the closure equations, the signed inner count for block A is

    K_A(n,c_A) = sum_(inner cyclic partitions)
                (-1)^(m_A-1) max(n-range(inner prefixes),0).

For a pair, K_2=|v| because |v|<=n-1. Therefore the direct count is

    S_P(n) = sum_(x in {0,...,n-1}^B, outer closure)
             product_(A in P) K_A(n,c_A).

It is the signed sum of the lifted Ehrhart polynomials, of period one
and degree at most B+1. Its leading coefficient is I_P. The integer
coordinate charts above establish the integration measure as well as
the polynomial degree; interpolation alone is not used as a proof.

At B=8, counts n=1,...,10 determine degree-nine polynomials and n=11
checks every orbit and the aggregate. At B=6, eight samples determine
degree-seven polynomials and n=9 is held out.

## Results

| Signature | Placements | Dihedral orbits | Exact model integral |
|---|---:|---:|---:|
| {4,2} | 15 | 3 | -23/420 |
| {6,2} | 28 | 4 | -563/11340 |
| {4,2,2} | 210 | 22 | -127/840 |
| {4,4} | 35 | 7 | 23/4536 |

The lower-order mixed class recovers the reference -23/420. The
{4,2,2} numerical candidate is confirmed. The earlier {6,2}=-1/20
candidate is rejected: the exact value is larger by 1/2835.

[The record](../results/mixed_cycle_flow_2026-10-05.json) contains all
integer counts, orbit values, bounds, source hashes and held-out checks.

## Independent checks and reproduction

The primary evaluator uses the spanning-tree chart, bounded outer x,
dihedral representatives and cached sorted inner frequency tuples. Its
integer arrays are processed in bounded chunks rather than storing a
full high-dimensional cube. The conservative absolute aggregate bound
is placement_count * product_A(term_count_A) * n^(B+1), below 2^63
for the supported signatures and n<=12. Inner cumulants fit int32;
products and sums use int64, and interpolation uses Python Fractions.

The standard-library verifier instead generates every placement with
restricted-growth strings, assigns independent closed frequencies,
reconstructs the outer walk and evaluates inner cumulants from all
frequencies. Its inner cyclic compiler uses permutations and cuts, with
cyclic rotation deduplication. It uses neither the primary coordinate
eliminator nor its frequency sorting or dihedral summation to compute
the scalar counts. All-placement counts at n=2,3 agree.

```sh
python scripts/mixed_cycle_lattice.py --out results/mixed_cycle_reproduction.json
python scripts/verify_mixed_cycles.py
```

The first command requires NumPy and supports resumable checkpoints
with pinned source hashes. The second uses only the standard library.
These certificates establish finite continuum-model integrals. Their
identification with arithmetic zeta moments is a separate obligation.
