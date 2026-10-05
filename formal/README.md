# Formal scope of the counting bridge

[CountingBridge.lean](CountingBridge.lean) contains 28 kernel-checked
theorems using Lean 4.33.0 and its standard library. It requires neither
Mathlib nor a package download. The file checks the scalar certificate
and the finite count implication conditional on the spectral inputs.

## What is proved

| Layer | Formal result |
|---|---|
| Integer accounting | Inertia/multiplicity inputs imply simple and distinct count bounds, including the collar. |
| Scalar certificates | Both denominators are positive and both certificates majorize the nonpositive half-line. The reflected-square certificate lies in [0,2] and its denominator dominates the parameterized denominator. |
| Finite spectral list | Pointwise majorization bounds the number of bad entries by the certificate sum; a one-tenth sum cap implies a one-tenth count cap. |
| Conditional counting conversion | The finite trace-cap input gives `4*N <= 5*S + 10*C + 9*delta`. |
| Resolvent minorant | A cleared polynomial-square identity gives a global quadratic lower bound for `x/(x^2+r^2)`, for positive r and S. |
| Exact rational arithmetic | The supplied cubic model moments give 247/2519, margin 49/25190 and conversion 2025/2519. A rational parameter witness is checked. |

Here N is the target multiplicity, S its simple-on-line count, C the
collar multiplicity, and delta a nonnegative dimension deficit satisfying
N<=d+delta. Dividing the formal conclusion by 5*N, when N>0, gives

    S/N >=4/5-2*C/N-(9/5)*delta/N.

The limit argument that the last two terms vanish is outside this file.

The parameterized scalar definitions clear the common denominator:

    qN(x)=2519-8232*x+7368*x^2-1932*x^3,
    f_alpha(x)=qN(x)^2/[2519^2*(1+alpha*x^2)^3].

The proof works in any linearly ordered field satisfying Lean's stated
field and order interfaces, for 0<=alpha<=1 and
2519^2*alpha^3<=1932^2. The parameter in the
[bounded bridge](../notes/bounded_resolvent_bridge.md) satisfies equality.
The file also checks alpha=4/5 in the rational field, so the admissible
parameter conditions have a concrete witness. The definition at alpha=0
is permitted for the pointwise inequalities; boundedness requires a
positive parameter and is not among the formal theorems.

The [reflected-square certificate](../notes/mirror_resolvent_certificate.md)
has the unconditional definitions

    D(x)=6345361+104885808*x^2+86095872*x^4+3732624*x^6,
    f_m(x)=qN(x)^2/D(x).

Lean checks `2*D(x)=qN(x)^2+qN(-x)^2`, denominator positivity,
`0<=f_m(x)<=2`, negative-side majorization, and the denominator
comparison with the admissible parameter family. The threshold and
finite counting theorems apply directly to f_m, without a parameter
hypothesis. The actual f_m trace cap remains a supplied hypothesis.

The resolvent lemmas check

    2*r*S^2*x-(-(r-S)^2+2*S*x-x^2)*(x^2+r^2)
      =(x^2-S*x+r^2-r*S)^2,
    (-(r-S)^2+2*S*x-x^2)/(2*r*S^2)<=x/(x^2+r^2).

The identity is algebraic; the inequality assumes r>0 and S>0.
Choosing S=(r^2+4/3)/(r+1) and taking the supplied first two
moments gives the [sharp two-moment lower bound](../notes/mirror_resolvent_certificate.md)
in the proof note. The moment limits, sharpness construction and
operator functional calculus are outside this Lean file.

## What is still an input

The list in the formal proof must be identified with the relevant real
matrix eigenvalues. The inertia estimate, zero-block multiplicities,
threshold perturbation and collar removal premises remain explicit
hypotheses. This file does not formalize the matrix spectral theorem,
the Weil explicit formula, zero tails, the archimedean/resolvent transfer,
or real-root construction and its connection to a concrete real-number
implementation. The bound f_m<=2 is formalized; boundedness of f_alpha
for positive alpha and Lipschitz regularity remain in the analytic
proof notes. Pole isolation and partial fractions have exact rational
checks, but are outside this Lean file.

The actual arithmetic certificate cap is also an explicit hypothesis.
The supplied model moments only enter a rational identity. They are
never substituted for actual zeta moments in the formal counting proof.
These proofs therefore do not establish the 80% zeta-zero theorem.

## Reproduction and trust

From the repository root, with Lean 4.33.0 installed, run

```sh
python scripts/verify_lean_counting.py
python scripts/verify_project.py --with-lean
```

The checker can also take `--lean-path /path/to/lean`, or the
`ZEROZETA_LEAN` environment variable. It locates an already installed
elan compiler when available. It never installs packages. To save a
fresh record, pass `--output results/lean_counting_reproduction.json`.

The [recorded audit](../results/lean_counting_2026-10-05.json) identifies
the source, checker and toolchain hashes and every theorem's transitive
axiom dependencies. The checker fails on a compiler error, warning,
missing theorem report, unapproved axiom or native-evaluation dependency.
Only the standard foundational axioms propext, Classical.choice and
Quot.sound are permitted; the printed record lists the subset used by
each theorem. The numerical decisions use kernel reduction.
The gate also runs seven rejection cases, including an actual Lean
compilation of a deliberately false stronger counting bound. The
temporary mutation is removed after the check.

Lean validates the elaborated statement relative to its definitions and
axioms. A reviewer must still check that those statements represent the
intended mathematical claims. See the
[official proof-validation guidance](https://lean-lang.org/doc/reference/latest/ValidatingProofs/).
