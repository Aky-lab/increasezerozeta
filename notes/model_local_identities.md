# Local identities in the finite continuum model

These statements concern the explicitly defined cumulant/prefix-walk
model. They do not supply its arithmetic transport to zeta moments.

## Singleton deletion

A singleton block forces its frequency to zero and has C1(0)=1. Removing
the corresponding zero increment only removes a repeated outer prefix
position, so it leaves the overlap unchanged. No free integration
variable is removed: the singleton already had no free frequency.
Deleting all singleton blocks therefore gives the same class integral
on the shorter cycle. The empty reduced partition has value one.

For a class whose non-singleton blocks occupy k positions, there are
binomial(B,k) endpoint sets on a B-cycle. The induced cyclic order
preserves each shorter-cycle placement, so its aggregate is multiplied
by binomial(B,k). This explains the frozen-singleton multiplicities in
the moment ledgers.

## Three-block vanishing

On the support of the outer overlap every |c_i|<=1. For a three-block,
write a+b+c=0. Its six cyclic-partition terms give

    C3(a,b,c) = 1 - ov(0,a) - ov(0,b) - ov(0,c)
               + 2*ov(0,a,-c).

The two three-singleton cyclic orders have the same range. For three
real frequencies summing to zero,

    range(0,a,-c) = (|a|+|b|+|c|)/2 = max(|a|,|b|,|c|).

Under the support bound this range is at most one. Substituting
ov(0,a)=1-|a| and ov(0,a,-c)=1-range gives C3=0 pointwise.
Every outer partition containing a three-block consequently contributes
zero. The integer-count version replaces one by n and has the same
cancellation when |c_i|<=n-1.

## Lower moment reconstruction

Let T2=1/3, T22=4/15, T222=32/105 and J42=-23/420 be the certified
pair and mixed aggregates. The pure-cycle certificates give
C4=-1/60, C5=1/36 and C6=-1/126. Bell enumeration and the identities
above give

    m2 = 1 + T2 = 4/3,
    m3 = 1 + 3*T2 = 2,
    m4 = 1 + 6*T2 + T22 + C4 = 13/4,
    m5 = 1 + 10*T2 + 5*(T22+C4) + C5 = 101/18,
    m6 = 1 + 15*T2 + 15*(T22+C4) + T222 + J42 + 6*C5 + C6
       = 640/63.

Thus the sixth model moment can be assembled directly from certified
class values rather than retaining an unevaluated reference remainder.
The [seventh-order](m7_ledger_derivation.md) and
[eighth-order](m8_ledger_and_80_target.md) notes extend this assembly.
All of these model moments still require the analytic spectral/counting
interface before any conditional conversion becomes a zeta-zero bound.
