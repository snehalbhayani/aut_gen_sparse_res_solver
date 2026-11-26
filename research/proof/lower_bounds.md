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

## Theorem P1b (row-counting / deficiency bound) — proved

For any valid template, rank M ≥ |G|, and the rank is at most the number of rows,
which is at most Σ_i |C ⊖ S_i| with C ⊖ S = {t : t + S ⊆ C}. Hence

    δ(C) := |C| − Σ_i |C ⊖ S_i| ≤ |E|.

With |E| = N_min this forces |C| ≥ s*(N_min), where
s*(N) = min{ |C| : δ(C) ≤ N }, and so

    |G| ≥ s*(N_min) − N_min.

**Thesis example.** N_min = 9 (8 real solutions + 1 ghost).

- Among boxes, δ(C) ≤ 9 first holds for [0,3]×[0,4]: 20 points, with 2 + 9 = 11
  fitting rows. So s*_box = 20 and the bound is |G| ≥ 11.
- The best template found has |E| = 9 and **|G| = 11** on exactly that box
  (`min_elim_2d_mineig.py`).

So the bound is tight, **provided boxes minimise the deficiency δ among all sets of
the same size**.

**Conjecture D (a quantitative form of the thesis's §4.2.1 hypothesis).** For
supports whose Newton polytopes are compatible with Q, the minimal deficiency among
column sets of a given size is attained, up to a lower-order term, by lattice points
of shapes compatible with the supports (Minkowski-type shapes).

Proving Conjecture D turns P1b into an explicit, shape-independent lower bound. The
natural tools are the discrete Brunn–Minkowski inequalities:

- (C ⊖ S) + S ⊆ C gives |C| ≥ |(C ⊖ S) + S|;
- this can then be bounded below via Gardner–Gronchi (2001) or Ruzsa-type
  |A+B| ≥ |A| + n|B| − n(n+1)/2.

The 2-D case should be provable this way.

## What this says about the thesis

The thesis's "convex sets fit more translates" (condition C1) is exactly the
statement that Minkowski-type shapes minimise the deficiency δ. Section 2 of
`research/convexity/README.md` showed this is not what decides *feasibility*;
ghosts decide that. It **is** what decides the *minimal elimination size at minimal
eigen size* (P1b). So the thesis intuition has a precise and provable home.

## Next steps

1. Prove Conjecture D in 2-D (lattice polygons; discrete Brunn–Minkowski) and test it
   against exhaustive search on small instances.
2. Compute the P1a and P1b bounds for the benchmark problems (6pt focal, 9pt radial,
   opt PnP, 7pt) and compare with the best published templates (11×26, 73×97, 87×114,
   34×56). For faces of dimension ≥ 2, λ_w needs a higher-dimensional
   certificate-size bound. Counting is weak when there are many equations (6pt focal
   has 10).
3. Lean: P1a is a short extension of `face_kernel`'s restriction argument plus
   `EigenSize.finrank_ker_le`. P1b is a rank ≤ #rows inequality. Both are within reach.
