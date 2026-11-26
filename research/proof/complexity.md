# Part 2: complexity of finding the smallest elimination template (working notes)

## Problem

**MIN-ELIM.**

- **Input:** polynomials f_1, …, f_m ∈ Q[x_1, …, x_n] in sparse encoding, a hidden
  index h, and an integer k in unary.
- **Question:** is there a valid template (columns C, rows x^t f_i with support in C,
  E ⊆ C ∩ (C − e_h), G = C \ E, with M[:, G] of full column rank) such that |G| ≤ k?

**Variant MIN-ELIM-N:** additionally require |E| = N_min, the minimal eigen size.

## Membership in NP

After squarifying, a valid template with |G| ≤ k has exactly |G| rows. A certificate
lists, for each row, the equation index and the multiplier t (binary exponents), plus
the touched columns, of which there are at most k · max_i |supp f_i|. Checking it
needs supports, set membership, and one rank computation of a |G| × |G| rational
matrix, all polynomial time. So **MIN-ELIM ∈ NP** when k is given in unary.

Columns of E that no row touches can be dropped together with their shift chains,
which keeps the certificate polynomial. With k in binary, templates may be
exponentially large and membership is unclear.

## Linear systems and "flat" templates

For affine-linear equations f_i = a_i · x + β_i, the natural templates are flat:
C = {1} ∪ {x_j : j ∈ V}, E = {1}, G = V ∋ x_h. The rows are the equations supported
inside V. Validity says these equations have rank |V| on V. So flat MIN-ELIM becomes:

> **MIN-CLOSED-SPAN.** Given column vectors a_c ∈ Q^p with supports S_c, and a target
> set B, find the smallest U ⊇ B such that {a_c : S_c ⊆ U} spans Q^U.

(This is the transposed form. U ranges over sets of variables, and B = {h} together
with the variables of the equation that contains h.)

The family of good sets U is closed under union, but there is no unique minimal one,
which hints at hardness.

**Easy special cases.** These are informal observations, not yet written proofs:

- **Supports of size ≤ 2, generic coefficients.** Spanning means every connected
  component of (U, available columns) contains a cycle. The minimum is a shortest
  cycle (or a "theta" structure) through B: polynomial time.
- **Difference-type columns e_u − e_v (graphic matroid).** The rank on U is at most
  |U| − (#components), so spanning needs "anchor" columns. The minimum is a shortest
  path to an anchor: polynomial time.
- **Biased graphs** (e_u − e_v + α e_z): the minimum is a shortest unbalanced cycle,
  again polynomial time.

**Why the obvious reductions fail.** Candidate source problems are:

- MINIMUM WEIGHT SOLUTION TO LINEAR EQUATIONS (Garey–Johnson MP5, NP-complete);
- SPARK, the smallest linearly dependent set of columns (NP-hard: Tillmann–Pfetsch,
  IEEE Trans. Inf. Theory 2014; to be verified);
- MINIMUM k-UNION.

A template's cost is the **union of the supports of the rows used**, not the support
of the resulting combination. So encodings of sparse solutions or circuits pay for
every column a row touches, and the cost no longer tracks the sparsity being encoded.
Graph-like (≤ 2-variable) equations are easy in every coefficient regime above.

**Working conjecture.** MIN-CLOSED-SPAN, and hence MIN-ELIM, is NP-hard. Hardness
needs (i) equations involving at least 3 variables and (ii) non-generic rational
coefficients. In the generic, structural case the problem becomes a density/matching
problem on a hypergraph. Whether that case is polynomial is itself open; it is
related to finding the smallest subhypergraph with at least as many edges as vertices
that satisfies Hall's condition.

If both hold, the computer-vision message is: **non-genericity is the source of
hardness**.

## Plan

1. Prove NP-hardness of MIN-CLOSED-SPAN by reduction from SPARK or MP5, with gadgets
   whose supports have size ≥ 3. The main obstacle is making the cost equal the
   encoded sparsity.
2. Prove that flat templates are optimal for affine-linear systems, i.e. that
   multiplying linear equations by monomials never reduces |G|. This lifts hardness
   from MIN-CLOSED-SPAN to MIN-ELIM.
3. Settle the generic (structural) case. If it is polynomial, that confirms the
   message above.
4. Approximation: is MIN-ELIM hard to approximate (cf. Amaldi–Kann 1998 for minimum
   unsatisfied linear relations)?
