# Prime-pair errors on arbitrary divisor subfamilies

This note discharges one arithmetic input using a published prime-correlation
theorem: arbitrary subfamilies of nonzero shifts h=qk inherit mean-square
logarithmic savings, even when every eligible modulus q is retained. It
also gives a weighted consumption bound. The result concerns a single
prime-pair correlation on a full dyadic position interval; a short-window
or multiple-chain replacement requires additional estimates.

## Published input and notation

Let Lambda be the von Mangoldt function. For X>=2 and a positive integer h,
write

    C_X(h)=sum_(X<n<=2X) Lambda(n)Lambda(n+h),
    Delta_X(h)=C_X(h)-S(h)X,

where S(h)=0 for odd h and

    S(h)=2 C_2 product_(p|h, p>2) (p-1)/(p-2)

for even h. Here C_2 is the twin-prime constant, and L=log(3X).

The published input is [Matomaki, Radziwill and Tao, *Correlations of the
von Mangoldt and higher divisor functions I*, arXiv:1707.01315v3,
Theorem 1.3(i)](https://arxiv.org/pdf/1707.01315v3). For fixed epsilon>0,

    X^(8/33+epsilon) <= H <= X^(1-epsilon),

and every B>0, all but O_(B,epsilon)(H L^(-B)) positive integer shifts
h<=H satisfy |Delta_X(h)|=O_(B,epsilon)(X L^(-B)). This is the h_0=0
case of the published statement. We use its exceptional-set formulation
directly rather than assume a variance estimate as an extra theorem.
The zero shift is excluded: the pair singular series is defined for
nonzero shifts, and the diagonal sum requires separate accounting.

## Uniform divisor-family theorem

Let E be any set of pairs (q,k) of positive integers with qk<=H. For every
fixed A>0, uniformly in the choice of E,

    sum_((q,k) in E) |Delta_X(qk)|^2
      = O_(A,epsilon)(H X^2 L^(-A)).                      (1)

No bound on the number of selected moduli is needed. In particular, this
applies to every subset Q of the positive integers and arbitrary restrictions
1<=k<=K_q<=H/q. Repeated shifts must be counted with their multiplicity.

### Proof

Let tau(h) be the divisor count. The multiplicity

    m_E(h)=#{(q,k) in E:qk=h}

is at most tau(h). The elementary divisor bounds are

    sum_(h<=H) tau(h) <= H(1+log H),
    sum_(h<=H) tau(h)^2 <= H(1+log H)^3.                 (2)

For the second inequality, tau(h)^2<=d_4(h) multiplicatively, because
(a+1)^2<=binomial(a+3,3) at every prime-power exponent a>=0. Summing the
four-factor divisor count and bounding the last factor count by H/(abc)
gives H times three harmonic sums. The first inequality is the two-factor
version. These arguments apply to real H>=1, with sums through floor H.

The uniform pointwise bound needed on exceptional shifts is

    |Delta_X(h)| = O(X L^2),  1<=h<=H.                  (3)

Indeed, C_X(h)<=O(X log(3X)^2) follows by bounding each Lambda by log(3X).
For S(h), separate p=3 and use

    1/(p-2)=1/p+2/[p(p-2)]

at p>=5. The second terms have convergent sum; Mertens' reciprocal-prime
estimate bounds the first terms by log log(3H)+O(1). Taking logarithms of
the positive Euler factors gives S(h)=O(log(3H))=O(L). Thus (3) requires
no prime-pair conjecture or upper-bound sieve.

Choose B large and let Z be the exceptional set from the published input.
On the good shifts, (2) gives

    sum_(h not in Z) m_E(h)|Delta_X(h)|^2
      = O_(B,epsilon)(H X^2 L^(1-2B)).

On Z, Cauchy-Schwarz and (2) give

    sum_(h in Z) m_E(h)
      <= |Z|^(1/2) [sum_(h<=H) tau(h)^2]^(1/2)
      = O_(B,epsilon)(H L^(3/2-B/2)).

Multiplying by the squared pointwise bound (3) bounds the bad contribution by

    O_(B,epsilon)(H X^2 L^(11/2-B/2)).

Taking, for example, B=2A+12 makes both contributions O(H X^2 L^(-A)),
which proves (1). The constants inherit the published theorem's dependence
on A and epsilon; this is an asymptotic result, not an effective height bound.

## Weighted consumption

For arbitrary complex coefficients a_(q,k) on E, Cauchy-Schwarz and (1),
applied with twice the desired logarithmic saving, yield

    |sum_((q,k) in E) a_(q,k) Delta_X(qk)|
      = O_(A,epsilon)(X sqrt(H) L^(-A)
                     [sum_((q,k) in E) |a_(q,k)|^2]^(1/2)).              (4)

This makes the required weight budget explicit. A target error o(V) follows
whenever the right side is o(V). A claimed consumption normalization must
supply that budget and must stay within the published H range. An exact
Fourier identity by itself supplies neither the prime-correlation estimate
nor the required size of the weights.

## Position windows require a separate scale check

The [pinned reference paper](https://github.com/JoshuaHKU/zeta-0.7947-reproduction/blob/d85bddfe9d8f12856fba735fc9cb3ca23b48b3a8/paper.tex)
uses a short-window prime-sum approximation in its MM leaf and a window-mass
replacement in W2. Neither follows solely by changing X to the window length
Y in a full-scale prime number theorem or Siegel-Walfisz estimate.

For example, suppose an endpoint estimate has the form

    psi(t)=t+O(X E(X)),  X<=t<=2X.

For an interval W=(u,v] within that range, of length Y=v-u, subtracting
endpoints yields only

    sum_(u<n<=v) Lambda(n)=Y+O(X E(X)),
    relative error = O((X/Y)E(X)).                       (5)

For a continuously differentiable weight f, integration by parts gives
the corrected weighted version

    |sum_(u<n<=v) Lambda(n)f(n)-integral_u^v f(t)dt|
      <= C X E(X) (|f(u)|+|f(v)|+integral_u^v |f'(t)|dt).                (6)

The same deterministic argument applies to an arithmetic-progression
endpoint estimate after subtracting its appropriate main density. For
f(t)=exp(2pi i eta t), the last factor is 2+2pi|eta|Y. A genuine short-window
endpoint discrepancy would replace X E(X) by its own proved bound.

The major-arc result in [MRT, Proposition 4.1](https://arxiv.org/pdf/1707.01315v3)
has position scale X and error X times an arbitrary negative power of log X.
Formula (5) retains that X. For fixed Y=X^theta with theta<1, even the stronger
classical endpoint error E(X)=exp(-c sqrt(log X)) does not give a vanishing
relative error by endpoint subtraction: its ratio is
exp((1-theta)log X-c sqrt(log X)), which grows. This shows an insufficiency
of that deduction, not failure of the true prime number theorem in a
particular power-short interval. No assertion about the optimal known
short-interval exponent is needed here.

An unrestricted window-mass assertion is actually false: for W=(6j-1,6j],
j>=1, the sole integer is 6j, which has at least two distinct prime factors,
so its von Mangoldt weight is zero while the interval length is one. Such
windows occur at arbitrarily large scales. The necessary window hypotheses
must therefore accompany an arithmetic singleton deletion. The finite-model
singleton identity remains valid independently of this arithmetic issue.

## Verification and remaining transport

```sh
python scripts/verify_prime_pair_transport.py
```

The checker independently compares divisor-family aggregation with a
shift histogram, verifies the divisor-moment domination, and tests the
good/bad-set Cauchy-Schwarz bookkeeping with exact rational errors. These
are finite checks of the proof's discrete identities; the asymptotic input
is the cited theorem.

This result closes the full-dyadic, one-chain h=qk aggregation step. It
does not close the short-position-window major/minor arc estimates, the
two-modulus four-prime variance, or the assertion that a multiple-chain
variance factorizes. The exact model moments and their zeta counting
conversion retain those analytic obligations.
