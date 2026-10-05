# Centered reciprocity for the network certificates

The network lift has more structure than a generic Ehrhart polynomial.
For a partition P of B cyclic labels into non-singleton blocks, its signed
count satisfies

    S_P(-n) = (-1)^(B+1) S_P(n),
    S_P(0) = S_P(1) = S_P(-1) = 0.

Consequently,

    B even: S_P(n) = n (n^2-1) R_P(n^2),
    B odd:  S_P(n) = n^2 (n^2-1) R_P(n^2),

where R_P has at most floor(B/2) coefficients. Its leading coefficient is
the class integral. This applies to individual placements, dihedral orbits
and signature aggregates. Three-block classes vanish on the outer support;
the evaluator covers the fourteen remaining signatures through order eight.

## Proof

Fix one cyclic-partition term in every inner cumulant. The
[network lift](mixed_cycle_flow_polytopes.md) has the form

    Q = ker(A) intersect [0,1]^L,  L=B+M,
    dim Q = d=B+1.

The lattice is ker(A) intersect Z^L, with the integer chart already proved
in the lift. Total unimodularity gives integral vertices. The equations
only involve coordinate differences, so A*1=0. In particular, the vector
whose every coordinate is 1/2 is strictly inside every cube inequality.
Q is full dimensional in ker(A). Its relative interior consists exactly
of the points satisfying 0<z_i<1 for every coordinate: strict inequalities
give a relative neighborhood, while an equality is on a supporting
hyperplane, since the all-ones direction varies that coordinate.

Write E(t) for the Ehrhart polynomial of Q, and E_int(t) for its relative
interior count at positive integer dilations. The standard
Ehrhart--Macdonald reciprocity theorem gives

    E(-t) = (-1)^d E_int(t).

An interior integer point in (n+1)Q has every coordinate in {1,...,n}.
Subtracting the all-ones vector preserves A*z=0 and gives a bijection with
the integer points in (n-1)Q. Therefore

    E_int(n+1) = E(n-1).

For the shifted polynomial S(n)=E(n-1), reciprocity implies

    S(-n) = E(-n-1) = (-1)^d E(n-1) = (-1)^d S(n).

Equality at every positive integer n proves a polynomial identity.
Also E(-1)=(-1)^d E_int(1)=0, since no integer coordinates lie strictly
between zero and one. Hence S(0)=0. These conclusions survive the signed
sum over inner terms.

At n=1 the only outer assignment is x_i=0. Every frequency is zero. For
a non-singleton block of size b, its signed cumulant is

    K_b(1,0) = sum_(m=1)^b (-1)^(m-1) (m-1)! {b brace m} = 0.

The identity follows by extracting the b-th coefficient of
log(1+(exp(t)-1))=t; {b brace m} is a Stirling number of the second kind.
Thus S_P(1)=0. Parity also gives S_P(-1)=0. For even B an odd polynomial
has a factor n. For odd B an even polynomial vanishing at zero has a
factor n^2. Dividing by these factors and n^2-1 gives the stated degree
bound on R_P. This proves the interpolation space before any computation.

## Exact interpolation

Let q=floor(B/2), and a=1 for even B or a=2 for odd B. Use n=2,...,q+1.
The q distinct nodes n^2 determine R_P uniquely. In particular,

    I_P = sum_(n=2)^(q+1)
          [ S_P(n) / (n^a (n^2-1)) ]
          / product_(m=2,...,q+1; m!=n) (n^2-m^2).

The independent count at n=q+2 is held out. Rational interpolation is
performed for every orbit as well as the aggregate. Since each lift is
an integral polytope in its specified lattice, (B+1)! I_P is an integer;
this is an additional necessary check, not a sufficient certificate.

At B=8 the generic method used n=1,...,10 plus held-out n=11. The reduced
method uses four fit counts n=2,...,5 and held-out n=6, with n=1 separately
checking the zero identity. For C8 the largest outer cube shrinks from
11^8=214358881 to 6^8=1679616 points. Both methods use exact integer counts;
the reduced method reuses the same bounded evaluator rather than claiming
an independent implementation.

## Results and reproducibility

[The reduced record](../results/reciprocity_2026-10-05.json) contains all
fourteen signatures, orbit polynomials, fit nodes, held-out counts,
overflow bounds and source hashes. Generation of this record took about
nine seconds on the recorded local environment; this is a single run,
not a hardware-independent benchmark.

Every one of the 119 earlier aggregate counts and 641 earlier orbit counts
matches the reduced polynomials, including counts beyond the new fit range.
The {5,2} integral now also has a period-one network certificate with just
three fit counts. Its three distance integrals are 5/504, 1/360 and 13/2520,
and its aggregate is 1/8. The older spectator certificate used a different
weighted lattice and a period-dividing-six bound; both give the same integral.

```sh
python scripts/reciprocity_certificates.py --out results/reciprocity_reproduction.json
python scripts/verify_reciprocity.py --record results/reciprocity_reproduction.json
```

Generation requires NumPy. The verifier uses only the standard library.
The method does not identify these moments with arithmetic zeta moments.

## Attribution and research scope

Ehrhart reciprocity and lattice-volume integrality are established results;
see Matthias Beck's [2018 reciprocity lectures, slides 11--12](https://matthbeck.github.io/papers/jcca.pdf).
The contribution developed here is their application to the shifted network
counts, the resulting factorization, and the reduced certificates for these
moment classes. Priority for this particular application has not been
established by a comprehensive literature review.
