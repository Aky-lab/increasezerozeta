# increasezerozeta

Research repository for rigorous improvements to unconditional lower bounds on the proportion of simple zeros of the Riemann zeta function on the critical line.

## Current focus

The immediate goal is to isolate and prove the arithmetic-to-continuum step behind the one-sided fourth-moment route, with particular emphasis on:

- the class-subtracted \(\ell^1\) universality estimate for the actual \(\gamma_d^{\mathrm{free}}\) convention;
- an exact evaluation of the continuum core candidate \(C_{\mathrm{core}}=-1/48\);
- a transparent remainder ledger separating the fixed sawtooth-tail charge from genuine \(o(1)\) terms.

## Current research branch

`research/class-subtracted-universality`

- [Proof note](notes/class_subtracted_universality.md): includes the actual class-(d) replacement used by the reference code and proves an explicit
  [
  \ell^1\text{ bound }\le 8\left(\frac1{p-1}+\frac1{q-1}\right)
  \]
  for distinct odd prime moduli.
- [Exact checker](scripts/verify_lemmas.py): verifies the local tensor algebra, class correction, rational continuum integrals, and the conditional (1082/1539\approx70.3053931\%\) consumption arithmetic.
- CI runs the exact checker on pushes and pull requests.

## Status

This is research work in progress. The local algebra and rational calculations above are suitable for review, but the surrounding fourth-moment analytic chain has not been independently certified here. Nothing in this repository should be described as an established new unconditional theorem until that chain has been checked.
