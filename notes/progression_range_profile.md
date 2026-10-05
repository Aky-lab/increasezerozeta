# Prime-position coefficient ranges and their model coverage

The [complete cutoff law](dilation_core_coverage.md) measures a bound on
the reduced dilations relative to the global height. A prime theorem is
often formulated instead relative to the physical prime positions.
This note maps that different range into the arithmetic model, proves
its exact coverage through the square-root profile, and supplies exact
certificates for several larger profiles. These are geometric and model
results, not a cubic prime theorem with growing coefficients.

## 1. The product cell and progression length

Retain the core definitions and weights in
[the quantitative core proof](arithmetic_core_limit.md). Write

    A=exp(ell)=T/(2pi), u=nu-1, P=T/theta=A^nu,
    beta_i=log(b_i)/ell, g=gcd(b_1,b_2), r=b_1/g, q=b_2/g.

At the central product cell, b_2 m=b_1 n=P. Thus the physical prime
position scales and their common progression scale are

    M=P/b_2, N=P/b_1, H_*=M/r=N/q=Pg/(b_1b_2).             (1)

This orientation agrees with the reduced physical lock rn-qm=J and
with proportional windows I_n=(q/r)I_m. The product P and the physical
prime position N are different quantities. For distinct base primes,
which carry the leading model term, g=1 and

    log M/ell=nu-beta_2, log N/ell=nu-beta_1,
    log H_*/ell=nu-beta_1-beta_2.                           (2)

The actual raw overlap also depends on k,J and integer endpoints.
Equation (1) is a position/length scale audit; it does not replace that
overlap by a consumption-normalized volume.

For 0<=tau<=1, define the coefficient profile by

    r<=M^tau, q<=N^tau.                                    (3)

Let C_ell^prof(tau) restrict the core's modulus/nu integrand to (3).
Off shared-base pairs, its limiting region is exactly

    beta_1+tau beta_2<=tau nu,
    tau beta_1+beta_2<=tau nu.                              (4)

In particular beta_i<=tau, because beta_other>=nu-1. Summing (4) also
gives the progression-length lower bound

    nu-beta_1-beta_2>=nu(1-tau)/(1+tau).                    (5)

For fixed tau<1 this leaves a positive-power progression length.
At tau=1 the profile covers the whole original model region. A theorem
that supplies (3) must still specify its norm, residue uniformity,
modulus averaging, saving and consumption weights. Ordinary distribution
of single-prime counts in progressions does not supply the resolved
four-prime cubic estimate.

## 2. Exact model theorem through the square-root profile

Uniformly for 0<=tau<=1/2,

    C_ell^prof(tau)=-P_prof(tau)/48+O(1/ell),                (6)
    P_prof(tau)=2 tau^5(59+40tau+5tau^2)
                    /[5(1+tau)^2(2+tau)].                (7)

Consequently the profile r<=sqrt(M),q<=sqrt(N) covers precisely

    P_prof(1/2)=107/600 =17.8333... percent                  (8)

of the leading core. Its complement retains 493/600. This is a statement
about the normalized leading model measure; it is not a fraction of
prime tuples or a new zeta-zero proportion.

For all 0<=tau<=1 the profile fraction is continuous and monotone and
obeys the useful bounds

    F(tau/(1+tau))<=P_prof(tau)<=F(tau),                    (9)

where F is the complete dilation CDF. The left containment follows from
beta_i<=tau/(1+tau) and nu>=1; the right follows from beta_i<=tau.
The closed expression (7) is asserted only through tau=1/2.

## 3. Derive the clipped profile geometry

Set beta_1=u+x,beta_2=u+y. Equation (4) becomes

    x+tau y<=s=tau-u, tau x+y<=s,
    x,y>=0, x,y<=(1-u)/2.                                  (10)

Thus 0<=u<=tau. For tau<=1/2 the last two upper bounds are redundant:
x,y<=s<= (1-u)/2. The remaining quadrilateral has vertices

    (0,0), (s,0), (s/(1+tau),s/(1+tau)), (0,s).

On this region the ordered overlap takes the small-branch form

    H=4min(x,y)+2[y+(y-x-u)_+ + min(x-y-u,y)_+].             (11)

To check the branch in (11), the maximum of x+2y on the quadrilateral
is 3s/(1+tau) for tau<=1/2. This is at most 1-2u, since the difference
after multiplying by 1+tau is (2tau-1)(1+u)<=0. Therefore the changed
overlap min(y,1-2u-x-y)_+ equals y throughout. The other three
choices are already valid on the original core region.

By symmetry, integrating min(x,y) on y>=x and its reflected triangle gives

    integral min(x,y)=s^3/[3(1+tau)^2],
    integral y=(3+tau)s^3/[6(1+tau)^2].

For (y-x-u)_+, the relevant triangle has
0<=x<=(s-u)/(1+tau), x+u<=y<=s-tau x. Its integral is
(s-u)_+^3/[6(1+tau)]. For min(x-y-u,y)_+, introduce a height z<=y
and z<=x-y-u. Put y=z+v, x=u+2z+v+h. The restrictive face is

    (2+tau)z+(1+tau)v+h<=s-u, z,v,h>=0.

Its simplex volume is (s-u)_+^3/[6(1+tau)(2+tau)]. Combining the
four terms in (11) yields the exact slice integral

    Q_tau(1+u)=(7+tau)(tau-u)^3/[3(1+tau)^2]
       +(3+tau)(tau-2u)_+^3/[3(1+tau)(2+tau)].               (12)

The leading measure has density 96uH and total mass one. Hence

    P_prof(tau)=96 integral_0^tau u Q_tau(1+u) du.

Use the exact integrals tau^5/20 and tau^5/80 for the two cubics in
(12); simplification gives (7). No polynomial degree or coefficient
is inferred from numerical fitting.

## 4. Uniform arithmetic passage

The ordinary-prime transform is -u+O(1/ell) uniformly on the entire
core domain. Coprime higher powers and prime 2 contribute O(1/ell)
absolutely; shared-base pairs contribute O(ell^(-2)) absolutely.
Restriction to (3) preserves those bounds. For distinct bases g=1,
so (4) is exactly the coefficient restriction after exception removal.

For fixed nu,beta_2 and tau>0 the beta_1 integration interval has
lower endpoint nu-1 and upper endpoint, clamped to the original interval,

    U(beta_2)=min(nu/2, tau(nu-beta_2), nu-beta_2/tau).

This endpoint is monotone in beta_2 and has total variation at most one
after clamping. Replacing the first prime-power measure by Lebesgue
measure costs O(1/ell), uniformly, by the core's CDF discrepancy and
the bounded variation of H. The resulting function
g(beta_2)=integral_(nu-1)^U H d beta_1 has variation bounded by a
constant times 1+Var(U); the bound uses the boundedness and uniform
coordinate Lipschitz constants of H. Its endpoint jumps are bounded.
Replacing the second measure therefore also costs O(1/ell) uniformly
in tau. At tau=0 the coprime domain is empty; the shared-base absolute
bound applies. The prefactor is 1+O(1/ell).

This proves (6), and the same argument gives the uniform model limit
for the full range 0<=tau<=1 with P_prof defined by (10). It does not
assert the small-profile expression (7) for tau>1/2.

## 5. Exact certificates for larger profiles

The independent rational calculation gives:

| tau | P_prof(tau) | Approximate percentage |
|---|---|---|
| 1/2 | 107/600 | 17.833% |
| 2/3 | 73/150 | 48.667% |
| 3/4 | 10243/15750 | 65.035% |
| 9/10 | 256636537/277981875 | 92.321% |
| 1 | 1 | 100% |

The [certificate record](../results/progression_profile_2026-10-05.json)
stores the exact slice breakpoints, cubic coefficients and held-out
checks. The degree bound follows before interpolation: partition each
original overlap by its affine minimum, maximum and zero-height planes.
Intersect with (10). These are bounded polyhedra in (u,x,y). Between
successive u-coordinates of their vertices, each slice polygon has a
fixed edge incidence and vertices affine in u. Triangulation gives
quadratic areas; its affine height gives a cubic slice integral.

The verifier enumerates all independent triples of partition planes,
retains their intersections in the base domain, and splits at their
u-coordinates. This may insert redundant breakpoints but includes every
relevant polyhedron vertex. On each interval four exact rational samples
determine the proved-degree cubic. Two further interior points and the
endpoints are checked independently against the original twelve
overlaps. Integration of u times that cubic is exact. This supplies
finite rational certificates for the displayed larger-profile values,
without asserting a new closed formula for them.

```sh
python scripts/verify_progression_profile.py
```

The verifier also compares (12) with original-overlap polygon integration,
checks (7), profile containment (9), and exact progression-scale identities.
This is an audit of a candidate coefficient range. It does not substitute
for the actual cubic prime input or the signed zeta consumption estimate.
The reference supplies the product-cell convention and overlaps; the
profile mapping, formula and certificates are derived here. Independent
review and priority assessment remain necessary.
