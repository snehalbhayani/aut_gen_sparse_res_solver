# Part 1: lower bounds on the elimination (matrix-inverse) size

Notation is as in `face_obstruction.md`. A template has columns C, eigen basis
E ⊆ C ∩ (C − e_h), elimination columns G = C \ E, and rows x^t f_i with
t + S_i ⊆ C. It is valid iff M[:, G] has full column rank.

The eliminated block A2 is |G| × |G|, and |G| = |C| − |E|. Two questions:

- **(1a)** min |G| over all valid templates, eigen size free. This is the
  stable, matrix-inverse end of the thesis's spectrum.
- **(1b)** min |G| subject to |E| = N_min, the minimal eigen size (= real
  solutions + unavoidable ghosts, Theorem 3). This is the same as minimising |C|.

## Theorem P1a (face-certificate bound) — proved

For a weight w and the initial system in_w(F), call a finite set L ⊂ Z^n
**complete** for in_w(F) if the multiples x^t in_w(f_i) with support in L span
every monomial of L. Let

    λ_w = min{ |L| : L complete for in_w(F) }   (= ∞ if in_w(F) has a torus solution).

**Theorem.** Every valid template with hidden x_h satisfies

    |G| ≥ max{ λ_w : w_h > 0 }.

*Proof.* Let D = max_{b∈C} ⟨w, b⟩ and let L = {b ∈ C : ⟨w,b⟩ = D} be the top
w-layer. Because w_h > 0, L ⊆ G (as in the proof of Theorem 1). Since M[:, G] is
injective, so is M[:, L].

A row x^t f_i has entries on L only if its top w-degree ⟨w,t⟩ + top_i equals D.
In that case its entries on L are exactly the coefficients of x^t in_w(f_i), and
that polynomial is supported in L. So M[:, L] is the template matrix of in_w(F) on
the column set L. Injectivity means these rows span all of V_L, i.e. L is complete.
Hence |G| ≥ |L| ≥ λ_w. ∎

If in_w(F) has a torus solution, no L is complete. That is Theorem 1 again.

**Lean.** Machine-checked in `research/lean/TopLayer.lean`:
- `layer_disjoint`: L ⊆ G;
- `layer_entry_initial`: the entries on L are the coefficients of in_w;
- `top_layer_card_le`: M[:, L] is injective, so at least |L| rows meet L.

**Two variables (exact).** The top layer lies on a line, and in_w(f_i) becomes a
univariate polynomial g_i of degree d_i in the edge coordinate. Fitting rows satisfy
|L ⊖ S_i| ≤ |L| − d_i (1-D Brunn–Minkowski). A complete L therefore needs
Σ_i (|L| − d_i)_+ ≥ |L|. For two equations this gives |L| ≥ d1 + d2, and the
Sylvester matrix achieves it when the g_i have no common root. So λ_w = d1 + d2.

**Tight on the thesis example** (Ex. 4.2.1 instance, x1 hidden). The right edge
has leading forms of degree 3 and 2, so |G| ≥ 5. The template on the box
[0,3]×[0,4] has **|G| = 5** (eigen size 15), with the right edge as its eliminated
block (`research/convexity/min_elim_2d.py`). So for (1a) the bound is exact on this
example. The optimal template is **not** a Minkowski sum: it is smaller than
NP(f1)+NP(f2) = [0,4]×[0,5].

## P1b (deficiency bound): RETRACTED as stated

An earlier version claimed |G| ≥ s*(N_min) − N_min from the inequality
δ(C) := |C| − Σ_i |C ⊖ S_i| ≤ |E|. The inequality is true, since it is just
rank ≥ |G| with rank ≤ #rows ≤ Σ_i |C ⊖ S_i|. But **it gives no lower bound**:
every small column set (|C| ≤ N) satisfies it trivially. The "first box with
δ ≤ 9 has 20 points" came from scanning only boxes that already contained the Newton
polytopes. Exhaustive enumeration over all Young diagrams gives s* = 1
(`research/convexity/young_bound.py`).

What survives:

- **Compression lemma (correct, kept as a tool).** If every support S_i is a down-set
  (e.g. a full box or a full simplex), pushing the points of C down in each column
  (or left in each row) never decreases |C ⊖ S_i|.
  - Before compression, a column with n points fits at most n − s + 1 vertical segments
    of length s.
  - After compression, the fitting positions in consecutive columns are nested initial
    segments, so their intersection is as large as possible.

  Hence extremal problems about fitting translates reduce to Young diagrams (down-sets
  in N^n).
- **(1b) on the thesis example is now solved: min |G| = 11** at eigen size 9 (computer-assisted
  proof, see below).

**What a correct (1b) bound must use** (constraints the counting argument ignores):

1. The equations that appear in the rows must have finitely many torus solutions, at most
   |E| distinct x_h-values among them. Otherwise the kernel contains a continuum of point
   vectors, which contradicts the eigenvector argument of Theorem 2. In the thesis
   example both f1 and f2 must be used, so C contains a translate of the 4×4 box NP(f1).
2. The face certificates of P1a hold simultaneously for every w with w_h > 0.
3. E + e_h ⊆ C, and the point vectors stay independent on E.

Constraint 1 alone already repairs the counting bound. See "The repaired counting bound" below,
which proves symgen optimal on 7pt (Theorem P1d).

## Theorem P1b′ (exact optimum for the thesis example) — proved, computer-assisted

**Theorem.** For the Ex. 4.2.1 instance with x1 hidden, every valid template with eigen
size 9 (the minimum: 8 solutions + 1 ghost) has **|G| ≥ 11**, and |G| = 11 is achieved on
the box [0,3]×[0,4]. So 11 is the optimal matrix-inverse size.

*Proof.*

1. **Both equations are used.** By Lemma L (`complexity.md`, machine-checked as
   `Hardness.LemmaL.card_le`), x1 takes ≤ 9 values on V(U_eq). V(f1) and V(f2) alone are
   curves, so both f1 and f2 occur in the rows. Hence C contains a translate of
   supp f1 = the full 4×4 box. Templates are translation invariant, so take it to be
   [0,3]².
2. **Locate the right edge.** By P1a, the rightmost column of C (x1 = D) is eliminated and
   contains ≥ 5 points (Sylvester count 3 + 2).
   - If D > 3, that column is disjoint from the box, so |C| ≥ 21 and |G| = |C| − 9 ≥ 12.
   - So D = 3, and the column holds ≥ 1 point outside the box. Already |C| ≥ 17 and |G| ≥ 8.
3. **Reduce to a finite search.** |G| ≤ 10 means |C| ≤ 19, i.e. at most 3 points outside the
   box, all with x1 ≤ 3.
   - Any fitting row that touches an extra point has most of its support in the box (an
     f1 row covers ≥ 13 box points, an f2 row ≥ 3). So extra points that some row reaches
     lie in [−1,3]×[−2,5].
   - An extra point that no row reaches is a zero column. A chain of such points must end
     at a non-shiftable, hence eliminated, column, which is infeasible.
   - So every candidate lies in the region [−4,3]×[−4,7].
4. **Exhaust.** `research/convexity/min_elim_1b_exact.py` checks all 23,132 such column
   sets with **exact rational ranks**. A set is accepted iff its corank is ≤ 9 and its
   boundary columns are eliminable. This is necessary for any valid eigen-size-9 template
   on C, using all fitting rows and the matroid step. None passes. ∎

So on this example the thesis-style box template is **optimal** for the inverse size at
minimal eigen size. It is not the Minkowski sum: [0,3]×[0,4] is smaller than
NP(f1)+NP(f2) = [0,4]×[0,5].

## Theorem P1c (6pt focal: symgen's template is optimal) — proved

**Theorem.** For the 6pt focal relative-pose formulation in `problems/` (a2 hidden), every
valid template with the minimal eigen size 18 has **|G| ≥ 12**. symgen's (and the Maple
generator's) 12 + 18 template is therefore **optimal** in matrix-inverse size at minimal
eigen size.

*Proof.*

1. **A trace equation is used.** By Lemma L (machine-checked as `Hardness.LemmaL.card_le`),
   the equations used must leave only finitely many a2-values on their common torus zeros.
   det F = 0 alone defines a surface, on which a2 takes infinitely many values. So at least
   one of the 9 trace equations occurs in the rows.
2. **So C is large.** Each trace equation has the same 30-monomial support, so C contains 30
   distinct monomials and |C| ≥ 30.
3. **Every valid template has |E| ≥ 18.** By Theorem 2 (`FaceObstruction.spurious_eigenpair`),
   the 15 torus solutions and the 3 face solutions at a3 = ∞ (w = (0,0,1), w_2 = 0) give
   eigenvectors of the solver matrix with eigenvalue a2. For generic data these 18
   a2-values are distinct. Distinctness is a Zariski-open condition, so one instance where
   they are distinct suffices. They are distinct in all 50 tested instances: the smallest
   gap is 2.6e-5, far above the solver's ~1e-12 error. So the
   eigenvectors are independent and |E| ≥ 18.
4. **Conclusion.** At eigen size 18, |G| = |C| − 18 ≥ 12. ∎

Only this formulation is covered. The published 11 × 26 template uses a different input
system.

## The repaired counting bound (P1b fixed by Lemma L)

The retracted P1b used the right inequality with the wrong search space. Every valid template
with eigen size N satisfies

    #rows(C) ≥ rank M_C ≥ |C| − N,                                   (★)

where #rows(C) counts *all* translates t + S_i inside C. Adding every fitting row keeps a
template valid, so all fitting rows may be assumed. P1b failed because tiny column sets satisfy
(★). **Lemma L removes them:** the equations used must have finitely many x_h-values on their
common torus zeros, so C must contain a translate of supp f for some f in every such "finite"
set of equations. Call that translate the **anchor**.

Two facts make the repaired bound computable for a target value c:

- **Unreached points can be deleted.** If no fitting row meets b ∈ C, then C \ {b} has the same
  rows and one point fewer, so (★) only gets easier there. A minimal violating C is therefore
  *anchor ∪ X*, where X is a union of exteriors (row minus anchor) of fitting rows.
- **Finitely many candidates.** |C| ≤ c bounds |X| ≤ c − |anchor|. Enumerating unions of row
  exteriors of that size is a finite search.

`research/convexity/support_bound.py` finds the anchors. It certifies a set of equations as
*infinite* when, at a true solution, its Jacobian has full row rank and its tangent space moves
x_h. That is an open condition, so it holds for generic data. Sets it cannot certify are treated
as finite, which can only weaken the bound.

## Theorem P1d (7pt focal+radial: symgen's template is optimal) — proved, computer-assisted

**Theorem.** For `problem_relpose_7p_fr_1s_partial_elim` (a3 hidden, generic data), every valid
template with the minimal eigen size 19 has **|G| ≥ 22**. symgen's 22 + 19 template is
therefore **optimal**.

*Proof.*

1. **Anchor.** The supports have sizes 6, 26, 35, 38, 38. Every set of equations not
   certified infinite contains f3, f4 or f5. So C contains a translate of S3 (35 points), S4 or
   S5 (38 points). Translate it to the origin.
2. **Minimal eigen size.** The 19 solutions have pairwise distinct a3, with smallest gap 0.07
   against a residual of 1e-11. By Theorem 2 every valid template has |E| ≥ 19.
3. **Counting.** Suppose |C| ≤ 40. By the deletion argument, C = anchor ∪ X with X a union of
   row exteriors and |X| ≤ 5 (anchor S3) or ≤ 2 (S4, S5). There are 71, 5 and 5 such exteriors,
   giving 436 + 9 + 9 candidate column sets. `research/convexity/min_elim_7pt_exact.py` counts
   their rows. Every candidate has #rows(C) ≤ |C| − 21, so (★) fails by at least 2.
   Hence |C| ≥ 41 and |G| ≥ 22. ∎

Steps 1 and 3 are pure counting. The data enter only through "these equation sets are infinite"
and "the 19 a3-values are distinct", both generic conditions.

## Theorem P1e (face-layer bounds for 9pt and opt PnP) — proved, computer-assisted

P1a says that for every w with w_h > 0, the top w-layer L of C lies in G and is *complete*
for in_w(F): the rows x^t in_w(f_i) inside L have rank |L|. The following necessary conditions
hold for every choice of coefficients, provided the coefficients of in_w(F) are non-zero.
`research/convexity/layer_search.py` tests all of them.

- **(a) No empty columns.** Every point of L lies in a row inside L.
- **(b) Face anchor.** Lemma L for the face system: the in_w-equations used in L have no
  common torus zero, so L contains a translate of one of their supports.
- **(c) Rank bound from syzygies.** rank ≤ #rows − #syz. Here #syz counts Koszul syzygies
  in f_j · (x^t in f_i) − in f_i · (x^t in f_j), for i < j and t + S_i′ + S_j′ ⊆ L, that have
  distinct largest rows (i, t + max S_j′). Such relations are triangular, hence independent.
- **(d) P1a inside the layer.** For every direction v, the v-top part L′ of L is met by at
  least |L′| rows of L.

A minimal complete L is *row-connected*: its row-disjoint pieces are complete on their own. So
it suffices to grow L from the face anchor by adding rows that meet it. The search processes the
sets in increasing size and stops at the first size where some set passes.

**Combination with the global anchor.** Lemma L puts a translate A of supp f_a inside C. The top
layer L either misses A or meets it exactly in A's top face. Hence

    |C| ≥ |A| + λ_w − |face_w(supp f_a)|,    where λ_w is the least |L| passing (a)–(d).

**9pt two radial (a3 hidden).**

- Take w = e3. In coordinates (a1, a2, a4) the four initial supports are down-sets:
  - S1′ = S2′ = {a1 + a2 ≤ 1, a1 + a2 + a4 ≤ 2}, with 7 points each;
  - S3′ = the unit simplex, with 4 points;
  - S4′ = {a1 + a2 + a4 ≤ 3, a4 ≤ 2}, with 19 points.
- **(b).** No proper subset of in_{e3}(F) is free of torus zeros. This was checked with a
  Gröbner count at a random instance; the three-equation subsets have positive mixed volume.
  So L contains a translate of S4′.
- **Search.** No row-connected L with |L| ≤ 28 passes (a)–(c), so λ_{e3} ≥ 29. Adding (d)
  does not remove the 29-point set (a full search up to 32 points took 1 hour).
- **Global anchor.** Every proper subset of F is certified infinite (`support_bound.py`), so f4
  is used and A is a translate of S4: 57 points, with a top face of 19.
- **Conclusion.** |C| ≥ 57 + 29 − 19 = 67, so **|G| ≥ 40** at eigen size 27, up from 30.
  27 is the minimal eigen size: the 27 solver eigenvalues are pairwise distinct on 20 random
  instances (smallest relative gap 1.6e-3), and Theorem 2 then gives |E| ≥ 27.
- **Sanity check.** symgen's own top layer has 39 points and passes (a)–(c) with slack 3.
- **The bound is not sharp.** The 29-point set that passes the count is *not* complete: its
  exact GF(p) rank is too low. So λ_{e3} ∈ [29, 39], and an exact-rank search could raise the
  bound further.

**opt PnP (a3 hidden).**

- Take w = (1, 1, 1). The initial forms are the cubic parts g_1, g_2, g_3 of the three
  equations, all supported on the full degree-3 triangle T3 (10 points). The face anchor is
  therefore automatic.
- **Search.** No row-connected L with |L| ≤ 26 passes (a)–(d). So **|G| ≥ 27** for every valid
  template, whatever its eigen size. The anchor bound gave nothing here.
- **A sharper, generic condition (e).** Shift a complete L into N³. Its rows are then
  polynomial multiples of the g_i, so every monomial of L lies in the ideal I = (g_1, g_2, g_3).
  - At three random instances, the degree-6 Macaulay matrix has rank 27 of 28, and no degree-6
    monomial lies in its row space. This is an open condition, so it holds for generic data, and
    then I contains no monomial of degree ≤ 6.
  - So the triangle bounding L must have side ≥ 7. The full side-7 triangle (36 points) is
    complete: its exact rank is 36. The Hilbert function of the complete intersection is
    1, 3, 6, 7, 6, 3, 1.
- **Exact value λ_{(1,1,1)} = 32.** With (e) added, the size-ordered search finds no set of
  ≤ 31 points. The first set passing is a 32-point hexagon, smaller than the 36-point triangle:
  rows a1 = 0…6 holding a2 ∈ [0,5], [−2,4], [−2,3], [−2,2], [−2,1], [−2,0], {−1}. Its rows
  have exact GF(p) rank 32 at two instances, and full rank at one instance implies full rank
  for generic data. So λ_{(1,1,1)} = 32 exactly.
- **Conclusion.** **|G| ≥ 32** for every valid opt PnP template, at any eigen size, for
  generic data. Conditions (a)–(d) need only non-zero coefficients; (e) is the generic
  condition certified above. symgen has |G| = 87.

## Status across the benchmarks (minimal eigen size, matrix-inverse size |G|)

| problem (hidden) | N_min | anchor bound | proved optimum | symgen |G| |
|---|---|---|---|---|
| thesis Ex. 4.2.1 (x1) | 9 | 16 − 9 = 7 | **11** (P1b′) | 11 |
| 6pt focal, `problems/` (a2) | 18 | 30 − 18 = **12** | **12** (P1c) | 12 |
| 7pt focal+radial (a3) | 19 | 35 − 19 = 16 | **22** (P1d) | 22 |
| 9pt two radial (a3) | 27 | 57 − 27 = 30 | open, in **[40, 90]** (P1e) | 90 |
| opt PnP (a3) | 27 | 20 − 27 < 0, vacuous | open, in **[32, 87]** (P1e) | 87 |

So the generator's templates are proven optimal on three of the five systems. For the other two,
P1e raised the lower bounds to 40 and 32. The anchor alone gave only 30 for 9pt and nothing
for opt PnP.

**The randomized search did not improve symgen on 9pt.** `random_reduce.py` ran 150 randomized
greedy reductions from every valid Minkowski start, with and without the recovery constraint.
None got below |G| = 90.

**What would close the gaps.** In both cases one face layer alone already costs 29–32 columns,
yet symgen's templates are about three times larger. A bound that captures that difference has
to use more than one layer, for example:
- the second-highest layer, which must be complete modulo the shifts of E;
- or the Koszul-corrected count (★) over all of C, not just the top layer.

## What this says about the thesis

The thesis's "convex sets fit more translates" (condition C1) is a statement about
the counts |C ⊖ S_i|. Section 2 of `research/convexity/README.md` showed that it is
not what decides *feasibility*; ghosts decide that. The (retracted) P1b shows that
row counts alone do not bound the elimination size either. Whether C1 plays a role in
a correct (1b) bound is open.

## Next steps

1. A rank-aware version of (★) for 9pt and opt PnP. It would subtract the Koszul syzygies
   x^t(f_j e_i − f_i e_j) that fit in C, and require the P1a layers inside G.
2. Compare with the best published templates (73×97 for 9pt, 87×114 for opt PnP). Those use
   other formulations, so the bounds above do not apply to them directly.
3. Lean: P1a is a short extension of `face_kernel`'s restriction argument plus
   `EigenSize.finrank_ker_le`. The counting step (★) is `EigenSize.finrank_ker_le` plus
   rank ≤ #rows.
