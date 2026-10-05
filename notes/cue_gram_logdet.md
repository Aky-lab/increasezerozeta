# Logarithmic control and no atom at zero in the Gram limit

The unique limiting spectral measure nu of the Haar-unitary Vandermonde
Gram model has no atom at zero. This follows from an exact expected
log-determinant identity, rather than the finite matrices' full rank alone.
The result concerns the Gram model; it does not establish a zeta-zero bound.

## Statements

Let G_n and its random empirical spectral measure mu_n be as in
[the Gram realization](cue_gram_model.md). Write H_n=sum_(k=1)^n 1/k
and let gamma be Euler's constant. Then

    (1/n) E log det G_n = H_n-1-log n -> gamma-1.                       (1)

Let log_+(x)=max(log x,0) and log_-(x)=max(-log x,0). Since G_n is
positive definite almost surely, both finite-matrix integrals are defined.
For every n>=1,

    E integral log_-(x) dmu_n(x)
      <= 1+log n-H_n + (1/2)sqrt((1-n^(-2))/3)
      <= 1+1/(2sqrt(3)).                                              (2)

For the unique limiting measure nu,

    integral log_-(x) dnu(x) < infinity,
    gamma-1 <= integral log(x) dnu(x) < 0,                             (3)
    nu({0})=0,
    nu([0,epsilon]) <= [1-gamma+1/(2sqrt(3))]/|log epsilon|,
    0<epsilon<1.                                                     (4)

The logarithm at zero is interpreted as minus infinity. The strict upper
bound in (3) uses m_1=1 and m_2=4/3, so nu is not a point mass.

## 1. Exact log determinant from the CUE two-point kernel

The eigenvalues z_j=exp(i theta_j) of Haar U_n are distinct almost surely,
by their absolutely continuous Weyl density. The Vandermonde determinant
therefore gives

    det G_n=n^(-n) product_(i<j) |z_i-z_j|^2.

Use normalized circle measure dt/(2pi), and put
D_n(t)=sum_(a=0)^(n-1) exp(i a t). The CUE two-point function relative to
normalized circle measure is

    rho_2(t,s)=n^2-|D_n(t-s)|^2.

This is the standard determinantal kernel, cited in
[the occupancy proof](cue_limit_determinacy.md). Let
g(t)=log|1-exp(i t)|. It is integrable and has Fourier coefficients

    g_hat(0)=0,  g_hat(k)=-1/(2|k|),  k!=0.

For completeness, expand log|1-r exp(i t)| for 0<r<1. Its series is
-sum_(k>=1) r^k cos(k t)/k. As r increases to one, dominated convergence
holds: for r>=1/2,
|1-r exp(i t)|^2=(1-r)^2+4r sin(t/2)^2>=2 sin(t/2)^2,
so the negative logarithm is bounded by an integrable constant plus
|log|sin(t/2)||. The positive logarithm is bounded by log 2. This proves
the claimed coefficients without an exchange of a nonabsolute Fourier
series with the expectation.

Writing L_n=sum_(i<j) 2log|z_i-z_j|, the pair-correlation formula yields

    E L_n = integral g(t) [n^2-|D_n(t)|^2] dt/(2pi).

The logarithmic singularity is integrable against this bounded kernel.
There is no diagonal term: rho_2 counts ordered distinct pairs, and the
factor of two in L_n cancels the unordered-pair factor. The finite expansion

    |D_n(t)|^2=n+2 sum_(k=1)^(n-1) (n-k)cos(k t)

gives

    E L_n=sum_(k=1)^(n-1) (n-k)/k = n(H_n-1).

Subtract n log n and divide by n to prove (1). Its finite value also
justifies the expected log determinant: the positive part of L_n is
bounded and its negative part is integrable.

## 2. Uniform control at small eigenvalues

The trace of G_n is exactly n, so integral x dmu_n=1. Since log_+(x)<=x,

    integral log_+(x) dmu_n <= 1.

Subtracting the normalized log determinant and taking expectation gives

    E integral log_-(x) dmu_n
      = E integral log_+(x) dmu_n - (H_n-1-log n)
      <= 2+log n-H_n.

The harmonic sum is at least log n by integral comparison, giving a coarse
bound of two. To get (2), apply the first two moments of the expected
probability measure bar_mu_n=E mu_n. Its mean is one and its second moment
is 4/3-1/(3n^2). Therefore

    integral log_+(x) dbar_mu_n
      <= integral (x-1)_+ dbar_mu_n
      = (1/2) integral |x-1| dbar_mu_n
      <= (1/2)sqrt((1-n^(-2))/3).

Subtracting (1) gives (2). In particular the expected empirical mass in
[0,epsilon] is bounded by the right side of (2) divided by |log epsilon|.
This is the uniform control that full rank
of each finite matrix, by itself, does not provide.

## 3. Passage to the limiting measure

Let bar_mu_n=E mu_n. These expected measures converge weakly to nu, and
their second moments are uniformly bounded by the exact moment identity.
The functions log_+(x) are consequently uniformly integrable: for R>1,

    integral_(x>R) log_+(x) dbar_mu_n
      <= integral_(x>R) x dbar_mu_n <= (1/R) integral x^2 dbar_mu_n.

Truncation and weak convergence give

    integral log_+(x) dbar_mu_n -> integral log_+(x) dnu.

For each K>0, min(log_-(x),K), extended with value K at zero, is a bounded
continuous function on [0,infinity). Weak convergence, followed by monotone
convergence in K, gives

    integral log_-(x) dnu
      <= liminf_n integral log_-(x) dbar_mu_n
      = integral log_+(x) dnu - (gamma-1).

The right side is finite. Thus nu({0})=0, and subtraction proves the lower
bound in (3). Jensen's inequality and integral x dnu=1 give the upper
bound zero. Strict concavity makes it strict because m_2=4/3 rules out
nu=delta_1. The first two moments of nu similarly give
integral log_+ dnu<=1/(2sqrt(3)), hence
integral log_- dnu<=1-gamma+1/(2sqrt(3)). Markov's inequality for log_-
proves (4).

This argument proves a lower bound for the limiting logarithmic integral,
not equality to gamma-1. Weak convergence could lose some negative logarithmic
mass at unusually small finite-size eigenvalues; (2) alone does not exclude
that loss.

## Verification and scope

```sh
python scripts/verify_cue_fluctuations.py
```

The Weyl constant-term verifier also checks E L_n=n(H_n-1) at n=1,...,5,
using the full eigenangle density and the exact logarithmic Fourier
coefficients. This is separate from the two-point-kernel calculation.
The proof of absence of a zero atom uses the uniform logarithmic estimate
and the established weak limit, rather than finite-size numerics.

The computed Christoffel bound nu({0})<=12241115/162540559 remains a valid
finite-moment certificate; the all-order Gram realization and logarithmic
argument strengthen its model conclusion to zero. An arithmetic argument
must prove the relevant transport before using either statement for zeta
zeros. The CUE kernel and logarithmic Fourier identities are established
background; priority for this application is not asserted.
