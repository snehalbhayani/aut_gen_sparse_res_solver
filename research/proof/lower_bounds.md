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
- **Status of (1b) on the thesis example:** upper bound |G| = 11 (achieved, eigen size 9,
  box [0,3]×[0,4]); lower bound |G| ≥ 5 (P1a). The gap is open.

**What a correct (1b) bound must use** (constraints the counting argument ignores):

1. The equations that appear in the rows must have finitely many torus solutions, at most
   |E| distinct x_h-values among them. Otherwise the kernel contains a continuum of point
   vectors, which contradicts the eigenvector argument of Theorem 2. In the thesis
   example both f1 and f2 must be used, so C contains a translate of the 4×4 box NP(f1).
2. The face certificates of P1a hold simultaneously for every w with w_h > 0.
3. E + e_h ⊆ C, and the point vectors stay independent on E.

Combining 1–3 with compression into a computable bound (an integer-programming
relaxation over down-sets) is the next step for (1b).

## What this says about the thesis

The thesis's "convex sets fit more translates" (condition C1) is a statement about
the counts |C ⊖ S_i|. Section 2 of `research/convexity/README.md` showed that it is
not what decides *feasibility*; ghosts decide that. The (retracted) P1b shows that
row counts alone do not bound the elimination size either. Whether C1 plays a role in
a correct (1b) bound is open.

## Next steps

1. Derive a correct (1b) bound from constraints 1–3 above, using compression to reduce
   to down-sets, and test it against exhaustive search on small instances.
2. Compute the P1a and P1b bounds for the benchmark problems (6pt focal, 9pt radial,
   opt PnP, 7pt) and compare with the best published templates (11×26, 73×97, 87×114,
   34×56). For faces of dimension ≥ 2, λ_w needs a higher-dimensional
   certificate-size bound. Counting is weak when there are many equations (6pt focal
   has 10).
3. Lean: P1a is a short extension of `face_kernel`'s restriction argument plus
   `EigenSize.finrank_ker_le`.
