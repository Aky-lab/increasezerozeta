# Research overview

## Definitions and analytic scope

The continuum cumulants and outer prefix walks follow
[JoshuaHKU/zeta-0.7947-reproduction](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/tree/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8),
revision `d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8`. A class signature
denotes the aggregate over all placements of those block sizes on a cycle.
Each block's frequencies sum to zero, with one frequency eliminated.

The [finite-model identities](notes/model_local_identities.md) prove
singleton deletion and three-block vanishing on outer-overlap support.
These identities assemble the finite model ledger. The class-d replacement
and arithmetic-to-continuum transport remain separate analytic obligations.
The [direct counting bridge](notes/spectral_counting_bridge.md) derives the
finite spectral implication and an adequate zero-side tail bound; its
actual arithmetic polynomial-trace hypothesis remains unproved.

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
3. Independently review the direct spectral/counting implication, normalized Schatten reduction and balanced-prime evaluation; prove the actual signed off-balance bound, suitable one-sided moment estimates, or the bounded resolvent criterion below.
4. Replace the obstructed fixed-P absolute-tail estimate with a signed estimate or an actual-minus-model remainder bound in the separate [fourth-moment route](notes/class_subtracted_universality.md).
5. Evaluate further moments and connected covariance coefficients efficiently, and extend the polynomial-statistic CLT to wider test-function classes.
6. Extend the [verified literature comparison](notes/literature_map.md) and review the network/Gram application and reduced certificate method.
7. Independently review the global CUE/sine-process comparison, its normalization, boundary terms and removal of the interaction cutoff.
8. Establish a relative cubic-uniformity or equivalent restricted-lock estimate with the actual dilation and consumption budgets in the [resolved frame](notes/resolved_lock_frame.md).

A certified finite integral, an assembled model moment and an established
arithmetic zeta moment are distinct stages. The remaining analytic work is
essential to any unconditional simple-zero bound.

The [bounded resolvent criterion](notes/bounded_resolvent_bridge.md)
uses f(x)=q3(x)^2/(1+(1932/2519)^(2/3)*x^2)^3. It majorizes the
nonpositive half-line and is globally bounded and Lipschitz. Normalized
Schatten-two perturbation and a block resolvent estimate transfer its
trace to the whole-circle prime operator with vanishing compression
error, without assuming a sixth-moment bound. At z=i*(2519/1932)^(1/3),
the sufficient arithmetic estimate is a real linear combination of
resolvent powers one through three bounded above by -9/20.
This analytically weaker alternative avoids the absolute higher-word
boundary cost; it does not evaluate the required arithmetic combination.

The direct bridge reduces the endpoint target to Tr(q3(B_T)^2)/d or
Tr(q4(B_T)^2)/d, where B_T is the actual normalized Weil matrix.
A limsup at most 1/10 suffices for an 80% simple-on-line bound. The
vanishing threshold can be removed by polynomial coercivity, without
assuming separate moment limits. The cubic model value is 247/2519,
with excess allowance 49/25190; the quartic value is 12241115/162540559.
The six- and eight-cycle Fourier supports reach radii 6 and 8 at unit
bandwidth, beyond the restricted correlation theorem's radius below 2.

The [prime-matrix proof draft](notes/actual_prime_trace.md) derives its
exact entries and removes the archimedean and pole terms in normalized
Schatten-six norm. A noncommuting polynomial-stability argument preserves
any finite trace cap. Uniform finite-frame leakage bounds, unique
factorization and Mertens estimates then evaluate the multiplicatively
balanced moments: D_(2r)->lambda^(2r)*{2^r}, while odd balanced moments
vanish. Higher prime powers and repeated prime bases are negligible.
This connects the exact pair-cycle integrals to specified actual prime
terms; it does not evaluate the full moments.

For the cubic certificate, the balanced limit is 5102709/31726805.
The precise remaining requirement is a signed off-balance upper bound
-3860057/63453610, approximately -0.060833. The model predicts
-1991744/31726805, leaving the same 49/25190 budget. The note defines
every off-balance sum with its finite kernel and multiplicative phase;
it does not assume that six-prime configurations reduce to four-prime
locks. These analytic deductions await independent review, and the
off-balance estimate remains unproved.

The [finite-frame second-moment proof](notes/prime_second_moment.md)
discharges O_2=o(1) with an explicit O_(lambda,chi)(X/(T*ell)) bound,
using Montgomery--Vaughan's periodic weighted cosecant inequality.
Opposite signs have a separable taper amplitude, while the same-sign
alias denominator is canceled by its overlap factor. A two-sided
projection-leakage argument improves the word error to O_j(log(d)/d).
This recovers the actual normalized moments 1 and 1+lambda^2/3 and
removes order two from the remaining cubic signed target. The result
is a reconstruction of established second-moment information; its
normalization and proof await independent review. Absolute summation
of the general word error only gives a vanishing envelope when
j*lambda/2<=1. At unit bandwidth this does not cover higher degrees.

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
enumerating later class integrals. The separate connected-cycle argument
below proves convergence of random empirical measures. Transport to zeta
zeros remains open.

The underlying Ehrhart, determinantal-process and moment-problem theorems
are established background cited in the proof notes. The applications are
derived in this repository; a claim of priority requires further literature
review and independent mathematical assessment.

[The connected multi-cycle theorem](notes/cue_gram_fluctuations.md) proves
that a joint cumulant of r normalized Gram traces, of total order L, has
a rational polynomial numerator of degree at most L+1 and normalization
n^(L+r). Only partitions connecting all outer cycles survive. The resulting
O(n^(1-r)) cumulant bound gives almost sure convergence of every spectral
moment and the empirical measure, on any coupling of the Haar marginals,
and a joint Gaussian limit at scale sqrt(n) for fixed polynomial statistics.

Exact Weyl certificates give the covariance coefficients c22=1/10,
c23=1/3 and c33=79/70, and third cumulant coefficient c222=2/45.
The two-dimensional covariance determinant is 11/6300>0. Independent
constant-term integration at n=5 supplies held-out checks. Published CUE
pair-statistic results are credited in the proof; priority for the broader
Gram application remains to be assessed.

[The logarithmic argument](notes/cue_gram_logdet.md) gives the exact identity
E[log det(G_n)]/n=H_n-1-log n and proves nu({0})=0. It also establishes
EulerGamma-1<=integral log(x)dnu<0 and the explicit bound
nu([0,epsilon])<=[1-EulerGamma+1/(2sqrt(3))]/|log epsilon|.
This uses uniform logarithmic control, not only finite-matrix full rank.
The model result does not supply the arithmetic counting interface.

[The rectangular bandwidth theorem](notes/cue_gram_bandwidth.md) extends
the strong limit and polynomial CLT to every fixed m/n->lambda in (0,1].
A fixed integral coordinate chart and boundary-slab mesh estimate replace
the square polynomial-reciprocity argument. The row measure has moments
m1=1 and m2=1+lambda^2/3 and no zero atom; the column measure is
(1-lambda)delta0+lambda*nu_lambda. The second-moment limiting variance is
2*lambda^3/15-(2*lambda-1)_+^5/(30*lambda^2). Exact independent integration
checks rectangular cases, and integer sums verify the finite variance
formula across its threshold. Local kernel entry convergence alone does
not establish a global sine-process spectral identification. A separate
[global comparison](notes/sine_gram_identification.md) establishes that
identification using compact rooted-walk integrals, ergodic averages and
a uniform square-integrable interaction tail. The sine-window empirical
measure converges almost surely in all moments to the CUE column law.
Finite sine matrices are positive definite, while their limiting zero atom
is exactly 1-lambda. At unit bandwidth the global sine law is nu_1.
This resolves the model comparison; it does not prove arithmetic transport
or transfer the finite-CUE CLT to growing sine windows.

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

## A short-window prime-model replacement

[Short-window Fourier transport](notes/short_window_fourier_transport.md)
uses the published short-interval uniformity theorems of Matomaki,
Radziwill, Shao, Tao and Teravainen. The deduction gives
sum_j D(j)^2=O_(A,epsilon)(X Y^2 log(3X)^(-A)) for integer window lengths
X^(1/3+epsilon)<=Y<=X^(1-epsilon), where D is a maximal Fourier/AP
discrepancy between Lambda and the specified presieved Lambda_sharp.
The integer-position passage uses unit-cell constancy, not an assumption
that a real exceptional measure automatically controls integer samples.

Bounded-variation weights and coupled products of single sums inherit
an error budget determined by weighted start histograms. Effective start
counts at least X times an inverse power of log X suffice for normalized
logarithmic savings. Highly concentrated fixed starts need further input.
A longer-window theorem supplies a pointwise alternative with a model
depending on each start. No independence or variance factorization is used.
The actual dilated lattice ranges, weight budgets and arithmetic-to-continuum
identities must still be established in an application.

## An explicit comb and averaged short-window prime pairs

[Presieved pair transport](notes/presieved_pair_transport.md) writes the
finite divisor model with its exact Ramanujan coefficients. An elementary
Rankin mean-square bound controls its difference from the full presieved
model. CRT residue counting and period orthogonality identify its pair
density with the usual singular series, up to an arbitrarily small
logarithmic error and the explicit interval error O(c_R^2 D^2).

Combining this calculation with the published short-window Fourier input
and Parseval gives

    sum_(X<=j<2X) sum_(0<|h|<Y)
      |sum_n Lambda(n)Lambda(n+h)1_(j<n<=j+Y)1_(j<n+h<=j+Y)
         -S(h)(Y-|h|)|^2 = O_(A,epsilon)(X Y^3 log(3X)^(-A))

for X^(1/3+epsilon)<=Y<=X^(1-epsilon). The same bound holds with
divisor multiplicities h=qk and gives an explicit weighted consumption
estimate. The diagonal is excluded. The formulation is an application of
published uniformity, not a claim of a new twin-prime theorem or established
priority. The finite checker independently compares Ramanujan expansions,
direct period sums, CRT counts, Rankin factors and autocorrelation Parseval.

The estimate averages over both positions and shifts. Restricting to a
shorter shift interval retains the Y^3 bound; no H Y^2 saving is inferred.
Prescribed starts, weighted and dilated four-prime locks and the actual zeta
counting interface remain distinct research targets.

## Resolved lock transport and an information-loss obstruction

[The resolved lock frame](notes/resolved_lock_frame.md) proves that the
one-dimensional convolution of autocorrelations gives an aggregated
four-point statistic, while the restricted squared lock count has a
two-dimensional transform. Two nonnegative integer sequences with the
same autocorrelation give positive-shift squared counts 292 and 220
against the same second sequence. Thus a general bridge to the resolved
variance cannot use those separate power spectra alone.

A full-lock stability theorem is valid for every fixed number of slots,
with maximal Fourier errors and the other slots' energies. Published
short-window input therefore gives an averaged full-lock prime-model
replacement, including dilated lattices and coupled start histograms.
Restricting a four-slot array to rectangle locks retains a Y^5 error
scale, whereas the raw rectangle's natural squared scale is Y^4.

The symmetric rectangle norm is exactly a cubic Gowers norm. A bounded
nonnegative quadratic-phase example has vanishing normalized maximal
Fourier error but nonvanishing normalized rectangle discrepancy. This
specifies why a generic Fourier-only repair is insufficient. The cited
prime U^3 input is qualitative and uses Lambda_w; it must not be changed
to an arbitrary logarithmic saving or to Lambda_sharp without proof.
These are algebraic and conditional transport results, not a falsification
of a separately proved prime-specific theorem or a new zeta-zero bound.

## Qualitative four-prime rectangle replacement

[Averaged rectangle transport](notes/averaged_rectangle_transport.md) now
uses the published W-tricked U^3 estimate and the other factors' Gowers
norms to control the restricted rectangle at its natural squared scale.
For a sufficiently slow w(X)->infinity, it proves

    sum_(X<=j<2X) sum_(k,h)
      |N_j(Lambda;k,h)-Z_Y(k,h)K_w(k,h)|^2=o(X Y^4)

in the short-window range X^(1/3+epsilon)<=Y<=X^(1-epsilon).
Here Z_Y=(Y-|k|-|h|)_+ and K_w is the explicit finite four-point local
product. It includes repeated-offset degeneracies instead of assigning
them an infinite distinct-prime singular series.

The proof keeps the admissible residue count: its normalization factor
is uniformly bounded. The exact odd-prime cube count is
p^4-8p^3+28p^2-44p+23, with B_2=1, certified by all 256 subset ranks and
494 square minors, and checked by direct finite-field enumeration.
Residue carries, integer exceptional positions, boundary trims and deleted
small-prime powers are accounted for. This avoids the earlier crude
logarithmic pointwise loss without asserting a quantitative U^3 theorem.

The result is a published-input deduction with qualitative decay. It
cannot automatically absorb logarithmic consumption weights, divisor
multiplicity, sparse starts or variable dilations. Those remain needed for
the actual zeta reduction. Independent review and priority assessment remain
open.

## Averaged rectangle singular-series tails

[The Euler-tail theorem](notes/rectangle_singular_series_tail.md) proves,
for every fixed positive integer r and integer w>=5,

    sum_(|k|,|h|<=Y, kh(k-h)(k+h)!=0) |S_4(k,h)-K_w(k,h)|^r
      <=C_r Y^2/w^r.

Positive squarefree-divisor expansions give bounded moments for the
coincidence factors. Expanding the reciprocal-prime tail and grouping
equal prime indices by set partitions gives uniform moments of order
w^(-r). Linear-form multiplicities transfer these estimates to the
two-parameter rectangle; no independence or period approximation is used.
Rare primorial shifts show why a height-uniform pointwise tail would fail.

The mean-square case replaces K_w by S_4 in the preceding actual-prime
estimate on nondegenerate shifts. Thus its joint position/shift average
has the full four-point Hardy--Littlewood main term with o(XY^4) error.
The prime saving is still qualitative. This Euler-product tail differs
from the class-subtracted Ramanujan denominator tail in the arithmetic
core obstruction; singular consumption weights and the zeta normalization
remain unresolved. Singular-series moment methods are established
background; the note states its restricted-family proof without a priority
claim.

## Physical dilations and exceptional-prime normalization

[The dilated rectangle calculation](notes/dilated_rectangle_geometry.md)
retains the raw resolved variance as a heterogeneous eight-vertex cube.
Splitting the two sequences into progressions of steps r and q gives the
deterministic criterion

    sum_(k,J)|N_(f,g)-N_(u,v)|^2<=16 rq M^4 delta^2 B^6,

where delta is the maximal progression U^3 discrepancy and B bounds the
ordinary progression U^3 norms. In balanced windows of lengths rH,qH,
M=O(H), so this matches the raw rq H^4 scale. Supplying those norms
uniformly for growing coefficients remains an analytic obligation.

The general passage from a global cyclic U^3 norm to a selected step-d
progression loses sqrt(d). A single-coset indicator attains that factor:
its global norm tends to zero as d grows while its selected progression
is identically one. This specifies why qualitative global cubic smallness
cannot alone supply the required growing-progression prime estimate.

At a prime dividing rq, the exact admissible cube count is
(p-1)(p^3-4p^2+6p-3). Its normalized factor is
(1-1/p)^(-3)[1+1/(p-1)^3]. The product is bounded by a constant times
[rq/phi(rq)]^3. For the reference's prime-power moduli rq has at most two
distinct prime divisors, so this bound is uniform in their bases and
exponents. It is also exactly the complete-period second moment of the
finite four-point singular series.

The note parameterizes each lock by a Bezout solution and gives the
exact raw window overlap and all local factors of the finite presieve.
It distinguishes that raw cardinality from the consumption-normalized
volume. The large-progression prime estimate, incomplete weighted ranges
and the final zeta normalization remain open. The scale audit records
the progression length before any short-window theorem is invoked.

## Weighted transport for a slowly growing dilation family

[The weighted dilation theorem](notes/weighted_dilated_prime_transport.md)
uses the published W-tricked prime input and the explicit progression
extraction to supply the replacement for fixed coprime r,q, and uniformly
for a sufficiently slowly growing family of prime powers or 1. For
balanced windows of lengths rH,qH and nonnegative start weights of total
mass A whose two marginal maxima are at most (A/X)log(3X)^C, it proves

    sum_(i,j) omega_(i,j) sum_(k,J)
      |N_(i,j)(Lambda;k,J)-Z_(i,j;r,q)(k,J)S_(w;r,q)(k,J)|^2
        =o(A rq H^4).

The range is X^(1/3+epsilon)<=H<=X^(1-epsilon). The two starts may be
coupled, and their weights need not be uniform; the diagonal family
satisfies the marginal condition. The local residue budget stays bounded,
while the norm loss is explicitly sqrt(max(r,q)) times the qualitative
prime discrepancy. The slowly growing cutoff has no asserted rate.

The full singular series replaces the finite main term on
k J(rqk-J)(rqk+J)!=0. A translated short-box Euler-tail estimate proves
this uniformly in the same growing family. Positive divisor expansions
retain the endpoint term O_s(log(3X)^(15s)), and repeated-prime moments
retain O_m((1+log log(3X))^m); the polynomial window length absorbs them.

Weighted Cauchy--Schwarz states the actual consumption norm. Logarithmic
concentration of start weights is charged against the quantitative
exceptional measure. Additional growing consumption norms, prescribed
starts and positive-power dilations are not inferred from the qualitative
saving. Those regions and the reference's zeta normalization remain open.
The result is a deduction from published uniformity input with new explicit
geometry and tail bookkeeping; priority and independent review remain open.

## Prescribed full windows and quantitative dilation coverage

[The full-window theorem](notes/full_window_prime_transport.md) instead
uses Tao--Teravainen's published global quantitative cubic uniformity.
Prefix subtraction gives the needed interval norm at every prescribed
start, provided the window length is comparable to its position. Paying
the sqrt(d) progression loss explicitly yields

    sum_(k,J:Delta!=0) |N_(i,j)(Lambda;k,J)-Z S_(r,q)(k,J)|^2
      <<_K rqH^4 (log log H)^(-kappa)

for balanced windows (i,i+rH],(j,j+qH], starts bounded by KrH,KqH,
and arbitrary coprime r,q<=(log log H)^kappa. The positive exponent
kappa depends on the exponent in the published input. The same finite
product bound retains degenerate locks. Nonnegative weights on these
starts need no marginal concentration condition. The translated Euler
tail proof tracks the growing factor [rq/phi(rq)]^8 for unrestricted
coprime moduli. Fixed-dilation default full windows of the reference
are covered; genuine short prescribed windows are outside this result.

[The model coverage theorem](notes/dilation_core_coverage.md) computes
what a reduced-dilation cutoff actually captures. For R=exp(delta ell),
uniformly for 0<=delta<=1/3,

    C_ell^red[R]=-59 delta^5/240+O(1/ell).

The relative leading fraction is 59 delta^5/5, or 59/1215 at delta=1/3.
The complete extension covers every 0<=delta<=1 with branch transitions
at 1/3, 1/2 and 2/3. Its new geometry identity expresses one overlap
as half a clipped cube, so every coefficient comes from exact integration.
At delta=1/2,2/3,3/4 the leading fractions are 11/32,959/1215,147/160.
For delta>=2/3 the remaining fraction is exactly
32(1-delta)^4-(224/5)(1-delta)^5. Rational bisection bounds the exponents
for 50%,90%,99% coverage by (0.551320,0.551321],
(0.734423,0.734424],(0.859563,0.859564], respectively.
The [standalone manuscript](papers/dilation_core_distribution.tex)
defines the arithmetic model and gives the full proof and reproduction route.
Every subpower reduced-dilation cutoff has a vanishing contribution;
its complement retains the entire -1/48 leading term. Shared-base
prime-power pairs do not change the formula because their full absolute
contribution is O(ell^(-2)). Exact clipping integrates the original
twelve overlaps independently of the reduced geometry.

This is a cutoff on physical dilations, separate from the Ramanujan
denominator obstruction. The quantitative prime theorem and the model
coverage theorem together specify a concrete remaining gap: positive-
power dilation transport and the signed, consumption-normalized zeta
interface. Neither theorem establishes a new zeta-zero counting bound.

## Coefficients relative to physical prime positions

[The progression profile calculation](notes/progression_range_profile.md)
records the product-cell scales P=T/theta, M=P/b_2,N=P/b_1 and
H_*=M/r=N/q=Pg/(b_1b_2). Off shared-base pairs, the coefficient range
r<=M^tau,q<=N^tau is exactly
beta_1+tau beta_2<=tau nu, tau beta_1+beta_2<=tau nu. It leaves
a progression exponent at least nu(1-tau)/(1+tau) for fixed tau<1.
This maps the model region to a candidate physical arithmetic range;
it does not replace the exact raw overlap or the zeta consumption volume.

For 0<=tau<=1/2 the exact leading coverage fraction is

    P_prof(tau)=2tau^5(59+40tau+5tau^2)/[5(1+tau)^2(2+tau)].

The uniform core asymptotic is -P_prof(tau)/48+O(1/ell). The square-
root profile captures only 107/600 of the leading model term. A proof
integrates the original overlap branches on a symmetric quadrilateral;
the Mertens measure passage keeps a uniform moving-endpoint variation bound.
For larger profiles the exact certificates give 73/150 at tau=2/3,
10243/15750 at 3/4 and 256636537/277981875 at 9/10. Their cubic slice
degree follows from affine vertices before interpolation; original-
overlap integration checks additional interior and endpoint samples.

The standalone manuscript now includes the physical scale map and the
closed profile theorem. Ordinary single-prime progression distribution
does not supply the cubic lock norms. Even a square-root cubic profile
would leave most of this leading core beyond its range. Averaging the
actual modulus weights, or a different global signed argument, remains
an analytic research target.
