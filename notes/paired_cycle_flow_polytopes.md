# Exact pair-cycle integrals

For a perfect matching P of the 2k cyclic positions, put +v_j at the
first endpoint of pair j and -v_j at its second endpoint. Let s be the
resulting prefix walk and ov(s)=(1-max(s)+min(s))_+. Define

    I_P = integral ov(s) product_j |v_j| dv_1 ... dv_k.
    {2^k} = sum_P I_P.

The overlap support forces every |v_j|<=1, so the pair cumulant
min(|v_j|,1) equals |v_j| on this support.

## Integral-flow lift and polynomial degree

Each pair cumulant is 1-ov(0,v_j). Expand the product into 2^k signed
terms. Lift the outer overlap with x_i=s_i-t in [0,1]. For each pair,
lift its chosen inner term with one bounded y variable for the constant
term (both endpoints in one block), or two bounded y variables whose
difference is v_j for the overlap term (one block per endpoint). The two
endpoint frequencies sum to zero in either case.

The conservation equations have one row per inner block. Each outer
x_i column has +1 at the block owning i and -1 at the block owning i-1;
each inner y column is an incidence column in its pair's inner cycle.
The outer cycle visits every block, so the incidence graph is connected.
For M inner blocks the equation rank is M-1 and the dimension is
2k+M-(M-1)=2k+1. Incidence matrices are totally unimodular, and all
bounds are integral. Every lifted polytope therefore has integral
vertices; see [the network-matrix argument](pure_cycle_flow_polytopes.md).

The variables v_1,...,v_k,t and one inner offset per pair provide an
integer coordinate chart with an integer inverse. In particular,
v_j=x_(a_j+1)-x_(a_j), t=-x_0, and an inner offset is a selected y
coordinate. Thus the chart identifies the full affine lattice, and its
coordinate volume is the original integral. No Euclidean surface-measure
factor is inserted.

By Ehrhart's theorem each lifted term has a period-one lattice-count
polynomial of degree at most 2k+1. For dilation N and n=N+1, summing
the signed inner counts gives product_j |v_j|: a constant inner term
has n choices and an overlap inner term has n-|v_j| choices. Consequently

    S_P(n) = sum_(v_j=-(n-1)..n-1)
             max(n-spread(s),0) product_j |v_j|

is a polynomial of degree at most 2k+1, with leading coefficient I_P.
Dihedral rotations and reflections act by integer coordinate changes,
so representative counts may be multiplied by orbit sizes.

## Exact results

| k | Matchings | Dihedral orbits | Aggregate |
|---|---:|---:|---:|
| 1 | 1 | 1 | 1/3 |
| 2 | 3 | 2 | 4/15 |
| 3 | 15 | 5 | 32/105 |
| 4 | 105 | 17 | 1661/3780 |

The two k=2 orbit values recover 7/60 and 1/30. At k=4, samples
n=1,...,10 determine degree-nine polynomials; n=11 checks every orbit
and their aggregate. The [integer record](../results/paired_cycle_flow_2026-10-05.json)
includes all counts, orbit values, integer bounds and source hashes.

The [independent standard-library verifier](../scripts/verify_pairing_model.py)
enumerates all matchings and scalar frequency tuples at n=2,3,4 without
dihedral reduction. It also evaluates adjacent and nested three-pair
integrals by exact integration on 24 ordered simplex cells. Those checks
identify the [three-pair grouping discrepancy](../results/pairing_model_audit.md).

## Reproduction and scope

```sh
python scripts/paired_cycle_lattice.py --out results/paired_cycle_reproduction.json
python scripts/verify_pairing_model.py
```

The first command requires NumPy; the second uses only the standard
library. These are finite continuum-model certificates. The
[complete eighth-order assembly](eighth_order_certificate.md) uses them with
the finite-model local identities. Arithmetic transport and the counting
interface remain necessary before any zeta-zero bound follows.
