# Three-pair grouping audit

The direct perfect-matching model has five dihedral orbits on six cyclic
positions. Grouping all five noncrossing matchings as copies of the same
integral loses the distinction between an adjacent and a nested pattern.

| Representative matching | Multiplicity | Exact integral |
|---|---:|---:|
| (01)(23)(45) | 2 | 3/70 |
| (01)(24)(35) | 6 | 1/90 |
| (01)(25)(34) | 3 | 17/420 |
| (02)(14)(35) | 3 | 1/180 |
| (03)(14)(25) | 1 | 1/70 |

Their sum is 32/105. Replacing the nested value 17/420 by 3/70 gives
131/420, an excess of 3*(1/420)=1/140.

The pinned [reference certification script](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/certification/exact_t222.py)
defines four representative integrals and reports the aggregation
5*T0+6*T1+3*T2+T3=131/420. Its listed T0 walk is the adjacent pattern;
the nested pattern is not a rotation or reflection of it.

## Independent exact calculation of the discrepancy

For adjacent pairs, the integral is

    integral ov(0,v,a,b) |v*a*b| dv da db = 3/70.

For (01)(25)(34), the prefix positions are 0,v,0,w,w+u,w.
Set a=w and b=w+u, an integer change of variables with determinant one:

    integral ov(0,v,a,b) |v*a*(b-a)| dv da db = 17/420.

Sort the four positions {0,v,a,b} in each of their 24 possible orders.
Their three nonnegative consecutive gaps g_1,g_2,g_3 satisfy sum(g)<=1.
On each cell the absolute values have fixed signs and the weight is a
homogeneous degree-three polynomial. Every monomial integrates exactly:

    integral_(g>=0,sum(g)<=1) g^alpha (1-sum(g)) dg
       = product_i alpha_i! / (|alpha|+4)!.

The [verifier](../scripts/verify_pairing_model.py) expands those polynomials
with integer coefficients and sums the rational integrals. It shares no
lattice-count or orbit implementation with the primary evaluator.

## Consequences for the supplied model ledger

Retaining the reference framework's other class inputs and conditional
frozen-singleton identities, the pairing correction gives

    m6 = 12809/1260 - 1/140 = 640/63,
    m7 = 862/45 - 7/140 = 3439/180,
    m8 = 3329/90 - 28/140 + A8 = 3311/90 + A8.

The audited six-moment algebra gives lambda_3(0)=247/2519 and the
conditional counting conversion 1-2*lambda_3(0)=2025/2519, approximately
80.389%. The earlier 1415/13891 Christoffel value is valid algebra for
the separately supplied reference moment 12809/1260; it is not the
value for the direct perfect-matching model.

The {5,2}=1/8 and C7=-17/360 certificates are unchanged. This audit
does not establish arithmetic transport, validate the entire reference
analytic framework, or claim an unconditional zeta-zero bound.
