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

## Theorem C1 (NP-hardness, flat templates) — proved

**Theorem.** MIN-ELIM restricted to flat templates (columns of degree ≤ 1) is
**strongly NP-complete**. It stays NP-complete for affine-linear systems with a
unique solution, so that the minimal eigen size is 1 and MIN-ELIM-N is NP-complete
too. It holds for coefficients in a small integer range and also for generic
coefficients.

**Reduction from SET COVER.** Take a universe X, sets 𝒮 = {S_1, …, S_m} covering X,
and a budget k. Variables: h, one x per element of X, one y_S per set. The linear
parts of the equations:

    f0       = h − Σ_{x∈X} x
    f_{S,x}  = y_S + r_{S,x} x            for every S ∈ 𝒮 and x ∈ S
    f'_S     = y_S + r'_S x_1(S)          with r'_S ≠ r_{S,x_1(S)}

Here x_1(S) is any fixed element of S. The constants are chosen so that a fixed point
v* is a solution: each f(v) = a·(v − v*). The whole system then has the unique solution
v*, because every x is covered.

**Claim.** The minimal flat template has |G| = 1 + |X| + OPT_SC.

*Proof.* A flat template has E = {1} and G = U ∋ h, and is valid iff the equations
supported in U have rank |U| on U.

- **h forces X.** f0 is the only equation containing h. So the column h is non-zero only
  if f0 is a row, which forces X ⊆ U. Write Y′ = U ∩ {y_S}, so |U| = 1 + |X| + |Y′|.
- **Rank splits.** f0 is the only row with an h-entry, so rank = 1 + rank(gadget rows on
  X ∪ Y′).
- **Gadget rows are a frame matroid.** Each gadget row is 2-sparse, an edge y_S – x with a
  "gain". The rank of 2-sparse vectors on vertex set V is |V| − (#balanced components).
  A component is balanced if every cycle has gain product 1; trees and isolated vertices
  are balanced (Zaslavsky's frame matroids).
- **Every component with a chosen y_S is unbalanced.** The double edge y_S – x_1(S) is a
  2-cycle with gain product r/r′ ≠ 1.
- **Valid ⇔ cover.** If Y′ covers X, every component contains a chosen y_S, so the rank is
  |X| + |Y′| and the template is valid. If some x is uncovered, it is an isolated vertex,
  so the rank is at most |X| + |Y′| − 1 and the template is invalid.

So the valid flat templates are exactly {h} ∪ X ∪ Y′ with Y′ a set cover. ∎

Membership in NP is shown above, so the problem is NP-complete. The coefficients can be
taken in {−1, 1, 2, 3}, so it is strongly NP-hard. With generic coefficients every cycle
is unbalanced, and the same proof applies.

**Brute-force check** (`research/convexity/setcover_reduction.py`): on 40/40 random
instances, the minimal flat |G| equals 1 + |X| + OPT_SC, with small-integer and random
rational coefficients.

**Correction of the earlier working conjecture.** The notes previously guessed that
hardness needs ≥ 3-variable equations and non-generic coefficients. **Both guesses are
wrong.** The reduction uses 2-variable gadget equations, and generic coefficients do not
help. The earlier "easy" cases (shortest cycle or path) only apply when the equations
containing h are sparse. Here f0 involves every element, and that turns the problem into
a covering problem. **The hardness is combinatorial (set cover), not caused by
non-genericity.**

## Plan

1. ~~NP-hardness of the flat linear case~~: done (Theorem C1).
2. **Lifting lemma (open; evidence and proof sketch).** For affine-linear systems, do
   non-flat templates (rows x^t f, columns of degree ≥ 2) ever give smaller |G|? If not,
   Theorem C1 holds for unrestricted MIN-ELIM.

   *Greedy evidence* (`lifting_test.py`). symgen's exact greedy search over degree-≤2
   column sets on 6 gadget instances never beat the flat optimum. It is inconclusive,
   though: on 4 of the 6 the greedy search did not even reach the flat optimum (e.g. 8
   against 6). That shortfall is itself a data point: **symgen's deterministic greedy
   reduction is measurably suboptimal** on instances where the problem is NP-hard.

   *Exhaustive evidence* (`lifting_exhaustive.py`), instance X = {0,1}, sets {0},{1}, flat
   optimum |G| = 5.

   - For a fixed column set C the best eigen basis is the maximal one,
     E = {e ∈ C : e + e_h ∈ C}, since enlarging G only makes validity harder.
   - So a column set is good iff its boundary columns are independent in the matrix of
     all rows fitting in C.
   - With degree ≤ 2, |G| ≤ 4 forces |E| ≤ 6 (only the 6 monomials of degree ≤ 1 have an
     h-shift of degree ≤ 2), so |C| ≤ 10.
   - All such column sets were checked (sizes up to 14). **None has |G| < 5.** So for this
     instance, no degree-≤2 template beats the flat optimum. The harness does find the
     flat optimum when allowed (|G| = 5).

   *Proof sketch.* Take the f0-row x^t f0 of maximal h-degree H − 1.

   - Each shifted column t + e_x (x ∈ X) is either eliminated at level H − 1, or it lies
     in E. In the second case its h-shift t + e_x + e_h lies at the top level H. There it
     can be spanned only by gadget rows based at t + e_h, since an f0-row based at
     t + e_x would contradict maximality.
   - Either way each element x costs one eliminated column. The anchored column t + e_h
     is one more.
   - The gadget rows around t + e_h form a shifted copy of the flat template, with the
     anchored column playing the role of the constant 1. The frame-matroid argument then
     forces the sets used, across the two levels, to cover X.
   - Each set used contributes one eliminated y-column at some level. Total:
     |G| ≥ 1 + |X| + OPT.

   The bookkeeping (an injective map into G across levels, and the rank conditions on the
   restricted blocks) still has to be written out and checked.

3. **Approximation.** The cost has an additive offset of 1 + |X|, so set cover's ln n
   inapproximability does not transfer directly. A padded or gap-preserving variant is
   needed.
4. **Genuine minimal problems.** Hardness for systems with N > 1 and the eigen size fixed
   at N (e.g. by tensoring the gadget with a univariate factor in a new variable).
5. **Lean.** The proof is finite linear algebra plus the frame-matroid rank formula; it is
   formalisable with the existing set-up.
