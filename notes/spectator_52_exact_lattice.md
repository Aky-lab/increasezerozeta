# Exact lattice evaluation of the corrected model {5,2}

This note evaluates the finite continuum model defined by the 150-term C5
compiler and the **actual** distance-one, two and three walks. It does not
prove that this model supplies the seventh zeta moment. The arithmetic
transport remains unresolved. The separate C7 model integral is evaluated
in [the flow-polytope certificate](pure_cycle_flow_polytopes.md).

The completed exact evaluation is

    U1 = 5/504,  U2 = 1/360,  U3 = 13/2520,
    J52 = 7(U1+U2+U3) = 1/8.

The held-out N=60 check passes for each component. The separate pure C5
anchor is exactly 1/36 and passes its three held-out scales. The raw record
is `results/spectator_52_exact_lattice_2026-10-05.json`.

## Integral and support

Put c5 = -c1-c2-c3-c4 and s0=0, sj=c1+...+cj. The walks and
spectator elimination are proved in `spectator_52_eliminate_v.md`.
For distance d=1,2,3 use

    B_d = {0,s_(d-1),...,s4},  A_d = {s0,...,s_(d-1)}.

Duplicates have no effect. In particular B3 does not contain s1.
If m=min B, M=max B, a=min A, A=max A, the eliminated spectator weight is

    J_d = (|M-a-1|^3 + |1+m-A|^3 - |m-a|^3 - |M-A|^3)/6

when both spans are less than one, and zero otherwise. The target is

    U_d = integral_[−1,1]^4 J_d(c) C5(c) dc,
    J52 = 7 (U_1+U_2+U_3).

The cube contains the entire support. Each c_i is a difference of two
positions within B_d or A_d: for d=1 all adjacent pairs are in B1;
for d=2 the first is in A2 and the others in B2; for d=3 the first
two are in A3 and the others in B3. The closing difference c5=-s4
is in B_d for all three distances. Therefore a nonzero weight requires
|c_i|<1 for all five frequencies.

C5 is the signed sum over set partitions of five labels and cyclic
orders of their blocks, with sign (-1)^(number of blocks-1). For each
order its contribution is (1-range(block-prefix sums))_+. There are
150 terms. This is the compiler already calibrated against the pinned
reference implementation, not a replacement normalization.

## Why nine integer samples suffice

All break hyperplanes have the form

    sum_(i in T) c_i = 0, 1, or -1.

For C5, two block prefixes differ by the sum over a subset of labels:
the prefix subsets are nested. Comparisons choose the minimum and maximum,
and a difference equal to one controls the positive part. For J_d,
differences of ordinary prefixes have the same property; the four absolute
values and the two span cutoffs introduce only offsets 0 or +/-1.
The cube boundaries are of the same form.

After eliminating c5, each nonzero normal is a binary four-vector or its
negative. Indeed a subset containing label five becomes minus the sum over
its complement. Any vertex solves four independent equations with these
normals and integer right hand sides. Every nonsingular binary 4x4 matrix
has determinant of absolute value 1, 2 or 3. One proof augments a matrix B
to the 5x5 matrix with top row and first column all ones and bottom-right
block 1-2B. Its determinant has absolute value 16|det B|; Hadamard gives
16|det B| <= 5^(5/2) < 56. Thus |det B| <= 3. The script also checks all
1,365 choices of four distinct nonzero binary rows with a separate integer
determinant expansion. Every chamber vertex consequently has denominator
dividing 6.

On each chamber the weight J_d*C5 is a polynomial of degree at most four.
With c=k/N its numerator 6N^4 J_d(k/N)C5(k/N) is a homogeneous polynomial
of total degree four in (k,N). On chamber faces the formulas agree: the
original weight is continuous. A disjoint face decomposition, or finite
inclusion-exclusion of closed chambers and their faces, gives the exact
lattice sum. No double counting of shared walls is permitted in that
argument; the implementation simply visits each lattice point once.

The weighted Ehrhart theorem says that a polynomial-weighted lattice sum
over an N-dilated rational polytope is a quasi-polynomial, with period
dividing a common denominator of its vertices. For a homogeneous weight
of degree j in dimension r its degree is at most r+j, and its top
coefficient is the integral of the weight. See the general rational
polytope statement in [Baldoni, Berline, De Loera, Koeppe and Vergne,
*Computation of the highest coefficients of weighted Ehrhart
quasi-polynomials*](https://arxiv.org/abs/1011.1602).

Expand each chamber numerator as sum_j N^(4-j) h_j(k), with h_j homogeneous
of degree j. Each term has degree at most r+j+4-j <= 8 in N. The periods
all divide six. Lower-dimensional shared faces have degree at most seven
and cannot change the leading coefficient. Hence

    S_d(N) = sum_(k in Z^4) Jnum_d(k,N) Cnum(k,N)

is polynomial of degree at most eight for positive N divisible by six,
and its leading coefficient is 6U_d. Here Jnum=6N^3 J_d(k/N)
and Cnum=N C5(k/N) are integers. No value at N=0 is assumed.

For the nine samples N=6,12,...,54, therefore,

    U_d = eighth_forward_difference(S_d) / (6 * 8! * 6^8).

The additional N=60 is held out and must equal the polynomial prediction.
This is finite exact interpolation with a proved degree and period bound;
it is not rational reconstruction from floating point estimates.

The independent pure C5 anchor uses (1-range(s0,...,s4))_+ C5.
Its numerator has degree two, so the lattice sum has degree at most six
on the same residue class. Seven samples N=6,...,42 recover its integral;
N=48,54,60 are held out. The required reference value is 1/36.

## Implementation bounds and symmetry

`scripts/spectator_52_exact_lattice.py` uses NumPy integer arrays and
Python integers/Fractions for the final differences. It permits N<=60.
The C5 numerator has absolute value at most 150N. Since each overlap
is at most one and its support in the spectator variable lies in [-1,1],
0<=J_d<=integral_[-1,1] |v| dv=1. Thus Jnum<=6N^3 and

    sum of absolute products <= (2N+1)^4 * 900 N^4 < 2^63.

Each cubed intermediate and their signed sum fit int32 under the script's
additional bound 4(5N)^3<2^31. Products and sums use int64. The final
interpolation is arbitrary precision. No approximate arithmetic enters
the integral values; elapsed timings alone are floating point.

C5 is invariant under permutations of the five frequencies: relabel its
partitions and cyclic block orders. Changing which block anchors an order
only translates its prefix positions, since the total frequency is zero;
the range is unchanged. Consequently C5 is evaluated once per sorted
zero-sum quintuple. Global negation preserves all ranges and J_d, allowing
paired summation of the positive/negative c1 slices. Both optimizations
are compared with unoptimized full-grid summation at N=1 and N=2.
At those scales all 2,118 spectator cases are also compared with the
independent Fraction-based reduction checker.

The JSON record contains the exact integer sums, input source hashes,
checks, rational component results, and held-out test. Reproduce with:

    python scripts/spectator_52_exact_lattice.py --out results/lattice_reproduction.json

This calculation evaluates only {5,2}; the C7 evaluation has a separate
certificate. Arithmetic transport and a new unconditional simple-zero
bound require further analytic work.
