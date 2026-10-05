# Literature comparison and scope of the Gram results

This comparison records verified background and the specific applications
developed in the repository. It is not an exhaustive novelty assessment.
Independent mathematical review remains necessary for the new proof notes.

## The zeta-moment motivation

The paper [More than two thirds of the zeros of the Riemann zeta function
lie on the critical line (2026), Section 7.2(f)](https://www-cdn.anthropic.com/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf)
already proposes a sine-process Gram moment interface, lists the moments
1, 4/3, 2, 13/4 at unit bandwidth, and discusses conditional all-order
counting consequences. Proposition 7.4 states the frequency-dimension cap.
The random-matrix connection and the aspiration to use higher moments
for counting are therefore existing background. Density one in its
conditional discussion does not establish the Riemann hypothesis.

The repository derives exact finite CUE identities and network certificates,
and supplies proofs of moment determinacy, strong spectral convergence,
polynomial fluctuations and logarithmic control for its stated model.
These are candidate contributions to be assessed against prior work.
An identification with a global infinite-window sine-process law and an
arithmetic zeta-moment transport theorem require separate proofs.

## Random Vandermonde matrices

[Ryan and Debbah, Asymptotic Behaviour of Random Vandermonde Matrices
with Entries on the Unit Circle (2009)](https://arxiv.org/pdf/0802.3570)
study Vandermonde Gram moments and polytope-integral coefficients. Their
definition assumes independent identically distributed phases. Here the
phases are correlated Haar-unitary eigenangles. The distinction matters:
an independent-phase moment formula cannot be imported as a CUE formula.
Neither the general use of Vandermonde matrices nor polytope volumes for
their moments is claimed as new.

## CUE pair statistics and fluctuation methods

[Soshnikov and Wu, A Note on Pair Dependent Linear Statistics with
Trigonometric Test Functions (2023), Proposition 2](https://escholarship.org/content/qt1888k21h/qt1888k21h_noSplash_8e205173bf8abb441168126ac4d9cd1a.pdf)
give finite CUE pair-statistic variance identities. The second Gram moment
is such a statistic with a size-dependent Fejer test function. Its covariance
formula, and cumulant-based Gaussian-limit methods, are established
background. Their fixed-test-function theorems should not automatically
be substituted for the varying tests used here.

The repository's [multi-cycle argument](cue_gram_fluctuations.md) handles
every fixed collection of polynomial Gram statistics by connected incidence
networks. The [bandwidth extension](cue_gram_bandwidth.md) changes the
outer coordinate box, distinguishes row and column laws, and derives an
explicit second-moment variance across the overlap threshold. The exact
finite variance is a specialization of the published pair-statistic formula.

## Review targets

- Compare the finite CUE network representation with existing dependent-phase
  Vandermonde and sine-process spectral literature.
- Review the connected-partition cancellation, integral lattice chart,
  reciprocity and the rectangular mesh-count error.
- Check normalization before relating the row law, column law and the
  frequency-dimensional trace convention of a zeta compression.
- Establish the missing global sine-process identification or delimit the
  comparison to the finite CUE model.
- Keep arithmetic moment identities separate from random-matrix theorems.

No search result or absence of a search result establishes priority.
