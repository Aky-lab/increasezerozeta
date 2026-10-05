# Research overview

## Definitions and analytic scope

The continuum cumulants and outer prefix walks follow
[JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8),
revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`. A class signature
denotes the aggregate over all placements of those block sizes on a cycle.
Each block's frequencies sum to zero, with one frequency eliminated.

The [finite-model identities](notes/model_local_identities.md) prove
singleton deletion and three-block vanishing on outer-overlap support.
These identities assemble the finite model ledger. The class-d replacement,
arithmetic-to-continuum transport and spectral/counting interface remain
separate analytic obligations.

## Exact model evaluations

| Class | Exact value | Method |
|---|---:|---|
| {2} | 1/3 | Pair-cycle lift |
| {2,2} | 4/15 | Pair-cycle lift |
| {2,2,2} | 32/105 | Pair-cycle lift and simplex audit |
| {2^4} | 1661/3780 | Pair-cycle lift |
| C4 | -1/60 | Pure-cycle lift |
| C5 | 1/36 | Pure-cycle lift and spectator anchor |
| C6 | -1/126 | Pure-cycle lift |
| C7 | -17/360 | Pure-cycle lift |
| C8 | 157/4032 | Pure-cycle lift |
| {4,2} | -23/420 | Mixed-cycle lift |
| {5,2} | 1/8 | Weighted spectator lattice and reduced network certificate |
| {6,2} | -563/11340 | Mixed-cycle lift |
| {4,2,2} | -127/840 | Mixed-cycle lift |
| {4,4} | 23/4536 | Mixed-cycle lift |

The [pair-cycle proof](notes/paired_cycle_flow_polytopes.md),
[pure-cycle proof](notes/pure_cycle_flow_polytopes.md),
[mixed-cycle proof](notes/mixed_cycle_flow_polytopes.md) and
[spectator proof](notes/spectator_52_exact_lattice.md) fix the volume
normalization and polynomial-degree bounds before interpolation.
Their records include source hashes and held-out integer counts.

The {5,2} distance components are U1=5/504, U2=1/360 and U3=13/2520,
with {5,2}=7(U1+U2+U3). The defining distance-three fixed prefix set is
{0,s2,s3,s4}. Adding s1 changes the integral; see
[the spectator audit](results/spectator_actual_walk_audit.md).

The [pairing audit](results/pairing_model_audit.md) separates adjacent
and nested noncrossing patterns, whose integrals are 3/70 and 17/420.
Exact integration on 24 simplex cells independently confirms the corrected
three-pair aggregate 32/105. The {6,2} certificate also rejects the earlier
-1/20 candidate, which differs from the exact value by 1/2835.

## Complete eighth-order ledger and consumption

The exact model sequence through eighth order is

    (m0,...,m8) = (1,1,4/3,2,13/4,101/18,640/63,3439/180,747361/20160).

The new eighth-order aggregate is

    A8 = {2^4}+{4,2,2}+{4,4}+{6,2}+C8 = 633/2240,
    m8 = 3311/90+A8.

[The full certificate](notes/eighth_order_certificate.md) includes the
pure C8 counts, an independent compiler and unsorted cube checks.
The [seventh-order](notes/m7_ledger_derivation.md) and
[eighth-order](notes/m8_ledger_and_80_target.md) notes give Bell-class
multiplicities and exact assembly.

For H_n=(m_(i+j)), the origin-mass bound is

    lambda_n(0) = 1/(e0^T H_n^-1 e0).

The [Christoffel derivation](notes/christoffel_tower.md) gives the
degree-three value 247/2519. The complete eighth-order model gives
lambda4(0)=12241115/162540559. Under the analytic counting interface,
1-2*lambda4(0)=138058329/162540559, approximately 84.938%.
This is a conditional conversion, not an established zeta-zero bound.
The [target geometry](notes/k8_target_geometry.md) separately describes
the positive moment cone and alternative supplied inputs.

## Open research

1. Independently review the structural proofs, finite certificates and normalization.
2. Prove the arithmetic moment identities and their continuum transport through eighth order.
3. Establish the spectral/counting interface for these supplied moments.
4. Replace the obstructed fixed-P absolute-tail estimate with a signed estimate or an actual-minus-model remainder bound in the separate [fourth-moment route](notes/class_subtracted_universality.md).
5. Develop efficient evaluation of further model moments and study fluctuations of the Gram spectrum.
6. Assess prior literature for the network/Gram connection and reduced certificate method.

A certified finite integral, an assembled model moment and an established
arithmetic zeta moment are distinct stages. The remaining analytic work is
essential to any unconditional simple-zero bound.

## Structural results beyond enumeration

For the lifted network polytopes, the all-ones direction and Ehrhart
reciprocity imply S_P(-n)=(-1)^(B+1) S_P(n). Non-singleton cumulants also
force zeros at n=0,1,-1. The resulting factorization needs only floor(B/2)
fit counts; [the proof and reduced certificates](notes/centered_reciprocity.md)
reproduce every class through order eight and every earlier network count.

For Haar U in U(n), form V_(a,j)=z_j^a/sqrt(n) from its eigenvalues, and
G_n=V V*. [The exact Gram identity](notes/cue_gram_model.md) identifies
n^(B+1) E[Tr(G_n^B)/n] with the same Bell-class network sum. It proves that
the model moments are limits of positive spectral moments, and that every
finite-size correction is an even power of 1/n. Independent Weyl integration
checks the entire moment assembly through order eight at n=1,2,3,4.

[The occupancy argument](notes/cue_limit_determinacy.md) gives the all-order
uniform bound E[Tr(G_n^B)/n]<=42^B Bell_(B+1). Carleman's criterion then
proves moment determinacy and weak convergence of the expected empirical
spectral measures to a unique probability measure on [0,infinity).
This resolves existence and uniqueness for the model sequence without
enumerating later class integrals. It does not prove convergence of random
empirical measures or transport to zeta zeros.

The underlying Ehrhart, determinantal-process and moment-problem theorems
are established background cited in the proof notes. The applications are
derived in this repository; a claim of priority requires further literature
review and independent mathematical assessment.

## Quantitative arithmetic model and a cutoff obstruction

[The arithmetic-core theorem](notes/arithmetic_core_limit.md) proves
C_ell=-1/48+O(ell^(-1)) for the explicitly defined class-subtracted
prime-power/sawtooth functional. An absolutely convergent Euler convolution
gives A(M)=log M+EulerGamma-2+O_alpha(M^(-alpha)(1+log M)), for every
0<alpha<1/2. Uniform coprime transform bounds and integrated shared-base
estimates justify exceptional-modulus removal. The twelve-overlap geometry
reduces to cube and simplex volumes, yielding weighted integral 1/96.

[The cutoff theorem](notes/fixed_cutoff_obstruction.md) proves the uniform
bound C_ell,P=O((1+log P)/ell). Every P=exp(o(ell)), including a fixed
power of log T, has a vanishing finite core, while its complementary model
tail tends to -1/48. An absolute envelope
for that same tail therefore has lower limit at least 1/48. This rules out
0.0111 as an asymptotic absolute-tail bound for the specified subpower-cutoff
convention; it does not rule out a differently defined remainder or signed
cancellation. The older 70.3054% figure remains only formal conditional
consumption, with that tail assumption unestablished.

The new standard-library checker independently verifies local coefficients,
Euler convolution, overlap geometry and finite Ramanujan decompositions.
The asymptotic proofs are mathematical arguments in the notes. Neither the
model theorem nor its checks prove the reduction of actual zeta moments to
the model functional.

## A discharged prime-pair aggregation input

[Prime-pair progression transport](notes/prime_pair_progression_transport.md)
derives, from Matomaki–Radziwill–Tao Theorem 1.3(i), the uniform estimate

    sum_((q,k) in E) |Delta_X(qk)|^2 = O_(A,epsilon)(H X^2 log(3X)^(-A))

for every subset E of positive pairs qk<=H, in the published shift range
X^(8/33+epsilon)<=H<=X^(1-epsilon). The proof treats exceptional shifts
with a divisor second-moment bound; it does not assume variance
factorization. A further Cauchy-Schwarz estimate makes the consumption
weight budget explicit.

This resolves the full-dyadic, one-chain divisor-family aggregation step.
The note also states the endpoint error retained when passing to a shorter
position window. Replacing the full-scale error X E(X) by Y E(X) requires
a separate theorem; unrestricted window-mass replacement has elementary
counterexamples. Short-window prime sums and multiple-chain variance
remain distinct obligations. The new checker tests the divisor identities
and aggregation bookkeeping with exact arithmetic.
