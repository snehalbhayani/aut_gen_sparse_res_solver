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

## Lemma L (used equations pin down the hidden variable) — proved, any system

Let a template be valid and squarified, i.e. its rows are independent and number |G|.
Let U_eq be the set of equations that occur in its rows. Then x_h takes **at most |E|
distinct values** on V(U_eq) ∩ (K*)^n.

*Proof.*

- **Points give kernel vectors.** Every ξ ∈ V(U_eq) ∩ (K*)^n is a zero of every row, so
  v(ξ) = (ξ^b)_{b∈C} ∈ ker M.
- **The solver matrix exists.** dim ker M = |C| − |G| = |E|. Restriction to E is
  injective (validity), hence bijective. So the solver matrix X, the shift
  v|_E ↦ (v_{e+e_h})_e on kernel vectors, is well defined.
- **Each point is an eigenvector.** X v(ξ)|_E = ξ_h v(ξ)|_E, with v(ξ)|_E ≠ 0.
- **Finitely many eigenvalues.** X has at most |E| eigenvalues. ∎

For affine-linear systems V(U_eq) is an affine subspace through the solution, so
x_h must be constant on it. Equivalently, **e_h lies in the span of the linear parts of
U_eq**: the used equations determine x_h. Also |G| = #rows ≥ |U_eq|.

## Theorem C2 (NP-completeness of unrestricted MIN-ELIM) — proved

**Theorem.** MIN-ELIM is strongly NP-complete with k in unary. This holds for arbitrary
templates (any monomial multiples and column degrees), with the eigen size free or fixed
at its minimum, already for affine-linear systems with a unique solution and coefficients
in {−1, 1, …, 10}.

*Proof.* Use the instance of Theorem C1 with **positive** gadget coefficients
r_{S,x}, r'_S > 0 and r'_S ≠ r_{S,x_1(S)}, and a solution v* with non-zero coordinates.
Take any valid template. Squarifying keeps it valid, so assume it is squarified. By
Lemma L:

1. **h is determined.** e_h ∈ span(lin(U_eq)). Only f0 contains h, so f0 ∈ U_eq. Comparing
   h-coefficients gives Σ_{x∈X} e_x ∈ span(lin(F)), where F is the set of used gadget
   equations.
2. **Span of a gadget component.** View F as a gain graph on X ∪ Y. Each component K
   contributes span = all of K^{V(K)} if K is unbalanced. If K is balanced it contributes
   the hyperplane π^⊥, where π is the potential with π_y + r π_x = 0 on every edge.
3. **Every element is in an unbalanced component.** With r > 0, along each edge the signs
   of π_y and π_x are opposite. In a connected component all element variables therefore
   get potentials of the **same** sign, so Σ_{x∈K∩X} π_x ≠ 0. Hence Σ e_x restricted to K
   is not in π^⊥ whenever K contains an element. An isolated element is balanced with
   span {0}. So every x lies in an unbalanced component.
4. **Counting.** An unbalanced connected component has ≥ |V(K)| edges. Its y-vertices are
   adjacent to its x-vertices, since the graph is bipartite, so together the components'
   y's form a set cover. Hence |F| ≥ |X| + OPT_SC.
5. **Conclusion.** |G| ≥ |U_eq| = 1 + |F| ≥ 1 + |X| + OPT_SC.

The flat template of Theorem C1 achieves equality, with eigen size 1 (the minimum, since
the solution is unique). So for both variants the optimum equals 1 + |X| + OPT_SC, and
SET COVER reduces to MIN-ELIM.

For membership in NP: a certificate gives ≤ k rows as (equation, multiplier) pairs. The
columns are the touched monomials, and E = {c ∈ C : c + e_h ∈ C}. Validity is one rank
test. ∎

**Consistency checks.** The degree-≤2 exhaustive search (`lifting_exhaustive.py`) and the
flat brute force (`setcover_reduction.py`) both agree with C2.

**Take-aways for solver generators.**

- **Finding the smallest template is NP-hard** in general, so every practical generator
  (syzygy reduction, basis sampling, greedy parameter search, the Laurent generator,
  symgen) is necessarily a heuristic unless P = NP.
- **The hardness is combinatorial,** a set cover over which equations to use. Generic
  coefficients do not remove it.
- **Lemma L is a cheap, general pruning rule:** the equations used must have finitely many
  x_h-values on their common zero set.

## Theorem C3 (hardness for genuine minimal problems, any N) — proved

**Theorem.** For every fixed N ≥ 1, MIN-ELIM-N is NP-complete. MIN-ELIM-N is MIN-ELIM
restricted to zero-dimensional systems with exactly N solutions, all with distinct x_h,
with the eigen size fixed at its minimum N.

*Construction.* Keep the SET COVER gadgets of C2 (positive gains, double edges, solution
(x*, y*)). Replace f0 by

    f0′ = Π_{i=1..N} (h − a_i) − ( Σ_{x∈X} x − Σ_{x∈X} x* ),

with the a_i distinct and non-zero. The gadgets force (x, y) = (x*, y*), and then
h ∈ {a_1, …, a_N}. So there are exactly N solutions with distinct h, and the minimal
eigen size is N.

*Upper bound.*

- Take E = {1, h, …, h^(N−1)}, so E + e_h ⊆ C.
- Take G = {h^N} ∪ X ∪ Y′ for a cover Y′. The rows are f0′ and the cover's gadgets.
- f0′ is the only row containing h^N, and `cover_spans` (Lean) shows that the gadget rows
  span X ∪ Y′. So the template is valid with |G| = 1 + |X| + OPT.

*Lower bound.* Take any valid template, squarified.

- **The used gadgets must pin down Σ x.** If they don't, there is a direction d with
  d_h = 0, lin(g)·d = 0 for every used gadget g, and Σ_x d_x = 1.
- **That gives infinitely many solutions of the used rows.** For t ∈ Q, put
  s(t) = Π(t − a_i) and ξ(t) = (x*, y*) + s(t)·d with h = t. Then ξ(t) satisfies f0′
  and every used gadget row.
- **They are torus points with distinct h.** Each coordinate is a non-zero polynomial in t,
  so it vanishes for only finitely many t. The h-values are pairwise distinct.
- **f0′ must be used.** If it is not, take s ≡ 0 and let h = t vary freely.
- **Contradiction with Lemma L** (`card_le`), which bounds the number of distinct
  h-values by |E|.
- **Conclusion.** f0′ is used and Σ e_x ∈ span(used gadgets). The sign-trick lemma
  (`gadget_count`) then gives |F| ≥ |X| + |Y_used| with Y_used a cover, so
  |G| = #rows ≥ 1 + |F| ≥ 1 + |X| + OPT. ∎

*Check* (`research/convexity/c3_check.py`): on 5 random instances with N = 2…4, the
system has exactly N solutions (h = the chosen roots), every cover gives a valid template,
and every non-cover gives an invalid one.

**Lean:** the lower bound is machine-checked as `Hardness.min_elim_lower_bound_poly`
(with `determines_sum_poly`). Achievability uses `cover_spans`.

**Meaning.** Hardness is not an artefact of N = 1. For problems with any number of
solutions, finding the smallest elimination template at the minimal eigen size is
NP-hard.

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
2. ~~**Lifting lemma**~~: **resolved by Theorem C2** (Lemma L gives the lifting directly).
   Earlier evidence, kept for the record: for affine-linear systems, do non-flat templates
   (rows x^t f, columns of degree ≥ 2) ever give smaller |G|?

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

3. **Approximation (open; analysis).** The optimum is 1 + |X| + OPT_SC. The additive
   |X| swamps any multiplicative gap in OPT_SC, so set cover's (1 − ε) ln n
   inapproximability (Dinur–Steurer 2014) does not transfer.

   The natural amplification gives each set a private chain y_S^(1) – … – y_S^(w) of
   positive-gain rows, with the only unbalanced cycle (the double edge) at the end of the
   chain. Element edges get gain 1, so every cycle through shared elements is balanced.
   The sign-trick lemma still applies, since the graph is bipartite with positive gains
   and X lies in one colour class.

   **Loophole:** Lemma L only needs *one* unbalanced cycle per connected component. Sets
   that share elements can share one paid chain, so the cost becomes
   #sets + w · #components, not w · #sets.

   Candidate fixes:
   - private element copies per set, linked to the shared element in a way that keeps
     components separate;
   - reducing from a problem whose cost is naturally per component, such as a group
     Steiner variant.

   Until one works, MIN-ELIM is known to be NP-hard but **not** known to be hard to
   approximate.

   **Structural obstruction (observation).** On the C1–C3 gadget instances, every cover has
   at most |X| sets. So the flat template from *any* cover costs at most 1 + 2|X| ≤ 2·OPT,
   and these instances admit a trivial 2-approximation. Any inapproximability result
   beyond a constant factor therefore needs gadgets in which the "forced" part (here the
   |X| element columns dragged in by f0) is small compared with the optimisation part.

   That points to two directions:
   - **Reduce the forced part:** make the hidden variable's equation sparse and force the
     covering requirement through a chain of determinations instead of one dense
     equation.
   - **Accept the result:** MIN-ELIM may well be constant-factor approximable on natural
     classes, which would be good news for generators.

   Next step: look for an APX-hardness reduction from bounded-occurrence SET COVER or
   VERTEX COVER, with a sparse f0 and a determination chain that provably cannot be
   short-cut.
4. ~~**Genuine minimal problems (N > 1)**~~: **done**, see Theorem C3 below. (The earlier
   z-level construction needed a per-level lemma; C3 avoids it by putting the N solutions
   into the hidden variable itself.)
5. ~~**Lean**~~: **done** for the algebraic core. See `research/lean/Hardness.lean`:
   - `min_elim_lower_bound` (every template ≥ 1 + |X| + |sets used|, with a cover);
   - `cover_spans` (achievability);
   - Lemma L (`card_le`, `determines_h`);
   - `squarify`.

   The proof is also simpler than on paper: the sign-trick lemma replaces the
   frame-matroid / connected-components argument. Not formalised: the template encoding
   of the flat optimum, and the NP-membership / Karp wrapper.
