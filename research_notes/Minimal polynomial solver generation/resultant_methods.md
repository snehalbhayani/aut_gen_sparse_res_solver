# Resultant-based methods for generating minimal polynomial solvers (computer vision)

Abbreviations used below:
- "Syzygy" = Larsson, Åström, Oskarsson, "Efficient solvers for minimal problems by syzygy-based reduction", CVPR 2017. It is the automatic Gröbner basis (GB) generator.
- "GFan" / "Heuristic" = Larsson et al., "Beyond Gröbner bases: Basis selection for minimal solvers", CVPR 2018.
- "Greedy" = Martyushev et al. 2022, the greedy parameter search for elimination templates.
- Template size "r × c" (Bhayani convention): invert or eliminate an r×r block, then solve an eigenvalue problem of size c−r. For example, 11×20 means an 11×11 inverse and a 9×9 eigenproblem.

Primary texts read in full: arXiv 1912.10268 (CVPR 2020 paper plus its supplement), arXiv 2007.10100 (ICPR 2020), arXiv 2301.06443v2 (the JMIV 2024 journal version) and the BMVC 2008 PEP paper.

## Q1. Bhayani, Kukelova, Heikkilä, "A sparse resultant based method for efficient minimal solvers" (CVPR 2020; arXiv 1912.10268)

### Takeaway
The method adds one extra polynomial x_i − λ to the system and hides λ. It builds a sparse resultant matrix from the integer points of a δ-shifted Minkowski sum of a subset of Newton polytopes, then shrinks it by randomized row and column removal. It orders rows and columns so that M = M0 + λM1 with M1 = [[0,0],[−I,0]]. The Schur complement X = A21 − A22·Â12⁻¹·A11 then gives a standard eigenproblem. The resulting templates match or beat the best Gröbner basis generators (Syzygy, GFan, Heuristic) on most of the 15 + 6 problems tested.

### Cited Findings
**Background the paper uses (Sec. 2)**
- The resultant needs n+1 polynomials in n variables. It vanishes iff the system has a common root, and det M of a square full-rank Macaulay-style matrix is a nontrivial multiple of it (after Cox, Little, O'Shea, *Using Algebraic Geometry*) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Hidden variable: treat x_n as a coefficient, so M(x_n)b = 0. Write this as a polynomial eigenvalue problem (M0 + M1·x_n + … + M_l·x_n^l)b = 0, then convert it to a generalized eigenvalue problem (GEP) Ay = x_n·By. The paper says this "relaxes the original problem", because eigenvectors need not respect the monomial structure of b, and it "introduces extra parasitic (zero) eigenvalues leading to slower polynomial solvers" — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- u-resultant: add f_{n+1} = u0 + u1x1 + … + unxn, with u1..un usually random, and hide the u's — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Emiris et al. (ISSAC'93; Canny–Emiris, J. ACM 47(3) 2000) compute sparse resultants by mixed subdivision. Emiris's offline/online solver (arXiv 1201.5810) is described as "usually quite large and not very efficient". [7, 11, 12, 17] all require as many polynomials as unknowns — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)

**The method**
- Augmented system: x′ = [x1..xn, λ], f_{m+1} = x_i − λ. The system can be overdetermined (m ≥ n). Multiplier sets T = {T1..T_{m+1}} and monomial basis B give a tall M·b = 0 with p ≥ ε = |B| rows — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Proposition 3.1 (block partition): split B = Bλ ⊔ Bc with Bλ = B ∩ T_{m+1}.
  - Top rows are the multiples of f1..fm. Bottom rows are x^α(x_i − λ).
  - This gives M = [[A11, A12],[A21, A22]] + λ[[0,0],[−I,0]].
  - Square up Â12 by deleting rows. Then (A21 − A22·Â12⁻¹·A11)·b1 = λ·b1, which is eq. 14.
  - If A12 is not of full column rank, the alternative partition Bλ = {x_m ∈ T_{m+1} | x_i·x_m ∈ B} gives (B21 − B22·Â12⁻¹·A11)·b1 = −(1/λ)·b1, which is eq. 16.
  - Speed depends on |Bλ| and the size of Â12. Accuracy depends mainly on the conditioning of Â12.
  - Source: [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Basis search (Sec. 3.2 and suppl. Alg. 1):
  - For each variable x_i (i = 1..n), for each subset F_sub of the augmented system, and for each displacement δ with entries in {−ε, 0, ε}: compute Q = NP0 (unit simplex) + Σ NP(f), B = Zⁿ ∩ (Q + δ), and T_j = {t : t + supp(f_j) ⊂ B}.
  - Accept B if Σ|T_j| ≥ |B|, min|T_j| > 0, rank(M) = |B|, and the Bc-column / top-row block has full column rank.
  - Keep the smallest such B.
  - The method is derived from Heikkilä's polytope-subset algorithm (ICCV 2017). It does not compute a mixed subdivision.
  - Source: [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Column removal (Sec. 3.3):
  - Try columns "one by one in a random order". For a column, take all rows with nonzeros in it (r1..rk) and all columns touched by those rows (c1..cl), and remove them all.
  - Accept the removal only if every T_i stays non-empty, p−k ≥ ε−l, rank(M_red) = ε−l, and M_red still has the Prop. 3.1 partition.
  - Then delete the extra rows. Rows are deleted first from the λ-block to shrink |Bλ| and hence X; any remaining rows are deleted "randomly chosen from the upper block".
  - Source: [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)

**Table 1 (CVPR 2020, main paper): "Our" vs original vs Syzygy vs GFan (#GB) vs Heuristic** — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)

| Problem (#sols) | Our | Original | Syzygy | GFan (#GB) | Heuristic |
|---|---|---|---|---|---|
| Rel. pose F+λ 8pt‡ (8) | 7×16 | 12×24 | 11×19 | 11×19 (10) | 7×15 |
| Rel. pose E+f 6pt (9) | 11×20 | 21×30 | 21×30 | 11×20 (66) | 11×20 |
| Rel. pose f+E+f 6pt (15) | 12×30 | 31×46 | 31×46 | 31×46 (218) | 21×36 |
| Rel. pose E+λ 6pt (26) | 14×40 | 48×70 | 34×60 | 34×60 (846) | 14×40 |
| Stitching fλ+R+fλ 3pt (18) | 18×36 | 54×77 | 48×66 | 48×66 (26) | 18×36 |
| Abs. pose P4Pfr (16) | 52×68 | 136×152 | 140×156 | 54×70 (1745) | 54×70 |
| P4Pfr (elim. f) (12) | 28×40 | 28×40 | 48×60 | 28×40 (699) | 28×40 |
| Rel. pose λ+E+λ 6pt‡ (52) | 39×95 | 238×290 | 149×201 | – | 53×105 |
| Rel. pose λ1+F+λ2 9pt (24) | 90×117 | 179×203 | 189×213 | 87×111 (6896) | 87×111 |
| Rel. pose E+fλ 7pt (19) | 61×80 | 200×231 | 181×200 | 69×88 (3190) | 69×88 |
| E+fλ 7pt (elim. λ) (19) | 22×41 | – | 52×71 | 37×56 (332) | 24×43 |
| E+fλ 7pt (elim. fλ) (19) | 51×70 | 51×70 | 51×70 | 51×70 (3416) | 51×70 |
| Abs. pose quivers† (20) | 68×92 | 372×386 | 203×223 | – | 68×88 |
| Rel. pose E angle+4pt (20) | – (failed) | 270×290 | 266×286 | – | 183×203 |
| Abs. pose refractive P5P† (16) | 68×93 | 280×399 | 199×215 | 112×128 (8659) | 199×215 |

- Markers: † = inputs pre-reduced by G-J elimination. ‡ = alternative formulation (eq. 16) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- For the 6pt λ+E+λ problem (52 solutions), the paper's text gives a 56×56 eigenproblem plus a 39×39 inverse, against the heuristic solver's 52×52 eigenproblem plus a 53×53 inverse, and says the resultant solver is "slightly faster". This means X can be larger than the number of solutions (spurious eigenvalues) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Supplement Table 2 lists more problems as Our / Original / Syzygy / GFan / Heuristic — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268):
  - Rolling shutter pose (8): 47×55 / 48×56 / 47×55 / 47×55 / 47×55
  - Triangulation from satellite images (27): 87×114 / 93×120 / 88×115 / 88×115 / 88×115
  - Optimal pose 2pt v2 (24): 176×200 / 192×216 / 192×216 / – / 192×216
  - Optimal PnP Cayley (40): 118×158 everywhere
  - Optimal PnP Hesch (27): 87×114 / 93×120 / 88×115 / 88×115 / 88×115
  - Unsynchronized relative pose (16): 150×168 / 633×649 / 467×483 / – / 299×315
- Supplement Table 1 compares against extended Heikkilä hidden-variable GEPs and a simulated Emiris u-resultant, given as (GEP) vs ours (Inv / Eig) vs u-resultant (Inv / Eig) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268):
  - F+λ 8pt: 12² vs 11/9 vs 15/9
  - Stitching: 24² vs 18/18 vs 31/18
  - E+λ 6pt: 30² vs 14/26 vs 44/26
  - quivers: 43² vs 68/24
  - f+E+f: 18² vs 12/18
  - 9pt λ1+F+λ2: 68² vs 90/27
  - E+fλ 7pt: 36² vs 61/19 vs 105/19
  - λ+E+λ: 110² vs 39/56
  - Satellite triangulation: 52² vs 88/27 vs 93/27
  - Unsynch.: 128² vs 150/18
  - Rolling shutter: 18² vs 47/8 vs 48/8
  - The u-resultant generator failed on 5 problems.
- Emiris's computational-biology example (3 equations in 3 unknowns, 15 coefficients, mixed volume 16): the u-resultant solver of arXiv 1201.5810 uses a 56×56 inverse and a 30×30 eigenproblem. The new method uses a 48×48 inverse and a 16×16 eigenproblem — [arXiv 1912.10268 suppl.](https://arxiv.org/abs/1912.10268)
- Stability protocol: 5K random instances; reports mean and median of log10 normalized residual, and fail% = share of instances with any residual > 1e−3. Main-paper Table 2, values as Ours / Syzygy / Heuristic (mean, fail%) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268):
  - f+E+f 6pt: −12.55, 0.52% / −12.09, 2.36% / −12.05, 1.44%
  - P4Pfr (elim. f): −12.86, 0% / −12.59, 0% / −12.73, 0.02%
  - λ+E+λ 6pt: −8.99, 14.66% / −6.92, 25.9% / −8.13, 26.46%
  - E+fλ 7pt: −11.29, 0.36% / −10.69, 7.58% / heuristic failed to extract all variables
  - E+fλ 7pt (elim. λ): −12.53, 2.34% / −11.99, 0.44% / −11.05, 5.70%
  - Refractive P5P: −13.03, 0% / −12.45, 0.10% / −12.23, 0.08%
- Real data: the P4Pfr 28×40 solver in RANSAC on the Rotunda dataset (62 GoPro images) gave focal error avg 0.080% / median 0.063% / max 0.266%, distortion 0.522 / 0.453 / 1.651%, and rotation 0.031 / 0.029 / 0.062°. This is comparable to the 40×50 solver of Larsson et al. 2017 — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268); [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Authors' conclusion: "there is no one general method (GFan/heuristic/resultant), which will provably return the smallest solver for every problem". For complex problems all methods should be tried — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)

### Inferences
- The size gains are largest where the original or Syzygy solvers were poor: f+E+f, quivers, refractive P5P, unsynchronized relative pose. Where GFan or Heuristic already found small templates, the resultant method ties them.
- The method has three stochastic or heuristic parts: the δ choice (3ⁿ shifts), the random column-removal order and the random row removal from the upper block. So template sizes can vary from run to run. The numbers in the JMIV table differ from the CVPR table (see Q3), which fits this.

### Gaps
- The CVPR paper gives no per-solver wall-clock timings; comparisons are by template size only.
- The paper gives no explicit "6pt focal relative pose" row with the x_i − λ variant beyond E+f 6pt (9 solutions) and f+E+f 6pt (15 solutions).
- No P4Pf (without radial distortion) row appears in these tables.

## Q2. Bhayani, Kukelova, Heikkilä, "Computing stable resultant-based minimal solvers by hiding a variable" (ICPR 2020; arXiv 2007.10100)

### Takeaway
This paper goes back to the hidden-variable PEP/GEP, building on Heikkilä's ICCV 2017 algorithm, with four fixes:
- allowing m ≥ n equations;
- an exhaustive search over hidden variables;
- an explicit rank test on M0(x_n);
- removal of the parasitic zero or infinite eigenvalues by row-column elimination.

The point is to avoid any matrix inverse, because the A and B matrices hold only raw input coefficients. The result is a larger and slower GEP that is much more stable in near-degenerate configurations. It was shown on three refractive-pose problems: P5Pr, P6Pfr and P2Pr.

### Cited Findings
- Abstract: the hidden-variable method "results in a larger eigenvalue problem than the action matrix and extra variable resultant-based methods; however, it does not need to compute an inverse or elimination of large matrices that may be numerically unstable". It is "more stable … especially in close to critical configurations" — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
- Prior hidden-variable work in CV: Kukelova et al. and Hartley & Li solved 5-pt and 6-pt relative pose; Hartley & Li used quotient-free Gaussian elimination, Levinson–Durbin iteration and root polishing. Those problems are "simple", because after hiding a variable the matrix is already square. A Dixon-resultant solver (Kasten et al., WACV 2019) was slow. Kukelova's modified Macaulay method "is not general" and can produce dependent or extra rows on sparse systems — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
- Steps of Heikkilä's ICCV 2017 algorithm, as summarized here — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  1. Hide x_n and compute the Newton polytopes in n−1 dimensions.
  2. Compute Minkowski sums of all subsets of the polytopes.
  3. For all δ_k ∈ {−ε, 0, ε}^{n−1}, form B = Z^{n−1} ∩ (Q_j + δ_k).
  4. Find multipliers T_{j,k}; accept the basis if every polynomial contributes and |T| ≥ |B|.
  5. Take the smallest such basis.
- Canny–Emiris instead sums all the Newton polytopes. Drawbacks of the earlier methods listed in the paper — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  - Canny–Emiris and Heikkilä assume square systems.
  - Heikkilä's method can return a rank-deficient M0(x_n).
  - The PEP-to-GEP conversion creates large sparse matrices with parasitic 0 or ∞ eigenvalues.
  - Bhayani CVPR 2020 needs an inverse that can be ill-conditioned.
- Parasitic eigenvalue removal: if column j of A is zero and the matching column of B has one nonzero (row i), delete row i and column j; the nonzero eigenvalues are preserved. This repeats until no such column remains. It is attributed to Kukelova et al. If column j of B has two nonzeros, first apply an elementary row operation G (Gauss step) and then remove — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
- P5Pr (16 solutions) — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  - Haner & Åström's original GB solver: 280×399 template with a 44×44 action matrix (QR basis selection).
  - Syzygy: 199×215 template with a 16×16 eigenproblem. Heuristic: same size.
  - Bhayani CVPR'20: 78×78 inverse with a 25×25 eigenproblem.
  - Canny–Emiris: basis of size 400. Heikkilä: basis of 35 but rank 32, so invalid.
  - Ours: basis of 40, giving an 80×80 GEP that reduces to a 36×36 GEP.
- P6Pfr (18 solutions after symmetry removal) — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  - Original: 648×917 template with a 41×41 action matrix. Syzygy: 636×654. Heuristic: 398×416.
  - Canny–Emiris did not terminate. Heikkilä: basis of 82 but rank 76.
  - CVPR'20: 248×248 inverse with a 52×52 eigenproblem.
  - Ours: basis of 145, giving a 290×290 GEP that reduces to 110×110. The authors note they did not remove the symmetries.
- P2Pr (6 equations in 4 unknowns) — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  - Syzygy: 913×937 template with a 24×24 eigenproblem. Heuristic: 597×621.
  - Haner's hidden-variable solution: 160×160 GEP.
  - CVPR'20: 142×142 inverse with a 32×32 eigenproblem.
  - Ours: a 160×160 resultant that reduces to a 124×124 GEP.
- Timings, Table I (MATLAB, 3.9 GHz i7, 1K runs), as Syzygy / Heuristic / CVPR'20 resultant / ours — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100):
  - P5Pr: 0.494 / 0.668 / 0.334 / 0.474 ms
  - P6Pfr: 4.62 / 2.08 / 1.53 / 5.19 ms
  - P2Pr: 10.57 / 4.96 / 0.729 / 6.56 ms
- Authors' verdict: their solvers "are slower than the fastest ones based on the extra variable resultant-based method", but "significantly numerically more stable" in close-to-degenerate scenes (points near the plane through C and the refractive normal) — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)

### Inferences
- There is a clear speed–stability trade-off between the hidden-variable GEP (no inverse, raw coefficients, larger QZ) and the u−x_i Schur-complement solver (smaller eigenproblem, but conditioning set by Â12).

### Gaps
- Stability is shown only in histograms (Fig. 1–2); no numeric fail% table was extracted.
- Only three refractive problems were tested; there is no evaluation on the standard Larsson benchmark set.

## Q3. Bhayani, Heikkilä, Kukelova, "Sparse resultant based minimal solvers in computer vision and their connection with the action matrix" (arXiv 2301.06443, v2 Sep 2023; J. Math. Imaging Vision 2024, doi 10.1007/s10851-024-01182-1)

### Takeaway
This is the journal extension of CVPR 2020. It adds:
- an iterative scheme over Newton-polytope subsets, motivated by algebraically constrained (non-generic) coefficients;
- two matrix-reduction procedures;
- a refreshed 21-problem comparison that includes Martyushev's greedy search;
- a formal equivalence between action-matrix and extra-polynomial sparse-resultant solvers.

Many templates are smaller than in the CVPR version. Examples: stitching 8×31, unsynchronized relative pose 59×80, quivers 37×71, E+fλ 7pt 34×56, rolling shutter 40×52.

### Cited Findings
- The abstract says sparse resultants were "less successful … primarily because the polytopes do not respect the constraints on the coefficients". The fix is a "simple iterative scheme to test various subsets of the Newton polytopes", plus the extra polynomial and a Schur complement. It also gives conditions under which AM and SRes solvers "are equivalent" — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Implementation detail: the extra polynomial is written f_{m+1} = x_k − u0, with δ entries in {−10⁻³, 0, 10⁻³}. Rank conditions are checked "for a random value of u0 ∈ C". Steps: [offline] favourable-set search → block partition → row-column removal → row removal; [online] fill M(u0), take the Schur complement, solve the eigenproblem of X (eq. 49 or 51) — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Definition 2 (equivalence): an AM solver and an SRes solver are equivalent iff (1) the action matrix M_f equals the Schur complement X, and (2) the G-J elimination of the template Ĉ can be swapped for Â12⁻¹Â11, with Ĉ the same size as the upper block [Â11 Â12]. Propositions C.1–C.3 prove the constructions in both directions. Starting from an AM solver fixes the extra polynomial and the favourable monomial set; starting from an SRes solver fixes the monomial order, the extended set F′ and the quotient basis — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- The authors note that both approaches are instances of normal-form methods (Telen & Van Barel 2018; Telen et al. 2018; Mourrain & Trébuchet 2008), and that "normal forms have yet not been incorporated in the automatic generators for minimal solvers in computer vision" — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Table 1 (JMIV), as Our / Original / Syzygy / GFan (#GB) / Heuristic / Greedy — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2):

| # | Problem (#sols) | Our | Original | Syzygy | GFan (#GB) | Heuristic | Greedy |
|---|---|---|---|---|---|---|---|
|1|F+λ 8pt‡ (8)|7×16|12×24|11×19|11×19 (10)|7×15|7×15|
|2|E+f 6pt (9)|11×20|21×30|19×28|11×20 (66)|11×20|11×20|
|3|f+E+f 6pt* (15)|11×26|31×46|25×40|31×46 (218)|25×40|11×26|
|4|E+λ 6pt (26)|14×40|48×70|34×60|34×60 (846)|14×40|14×40|
|5|Stitching fλ+R+fλ 3pt (18)|8×31|54×77|48×66|48×66 (26)|18×36|18×36|
|6|P4Pfr (16)|52×68|136×152|140×156|54×70 (1745)|54×70|52×68|
|7|P4Pfr elim f (12)|28×40|28×40|58×70|28×40 (699)|28×40|28×40|
|8|λ+E+λ 6pt‡ (52)|39×95|238×290|154×206|–|53×105|39×95|
|9|λ1+F+λ2 9pt (24)|77×104|179×203|165×189|105×129 (6896)|73×97|76×100|
|10|E+fλ 7pt (19)|34×56|200×231|181×200|69×88 (3190)|65×84|55×74|
|11|E+fλ 7pt elim λ (19)|22×41|–|52×71|37×56 (332)|24×43|22×41|
|12|E+fλ 7pt elim fλ (19)|51×70|51×70|51×70|51×70 (3416)|51×70|51×70|
|13|Rolling shutter pose (8)|40×52|48×56|47×55|47×55 (520)|47×55|47×55|
|14|Satellite triangulation (27)|87×114|93×120|88×115|88×115 (837)|88×115|87×114|
|15|Refractive P5P* (16)|57×73|280×399|57×73|57×73 (8659)|57×73|57×73|
|16|Abs. pose quivers (20)|37×71|372×386|156×185|–|68×88|65×85|
|17|Unsynch. rel. pose‡ (16)|59×80|633×649|159×175|–|139×155|139×155|
|18|Optimal PnP Hesch (27)|87×114|93×120|88×115|88×115 (837)|88×115|87×114|
|19|Optimal PnP Cayley (40)|118×158|118×158|118×158|118×158 (2244)|118×158|118×158|
|20|Optimal pose 2pt v2 (24)|112×146|192×216|192×216|–|192×216|139×163†|
|21|E angle+4pt (20)|–|270×290|266×286|–|183×203|99×119†|

- Notes on Table 1 — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2):
  - Markers: * = simplified input system as in Martyushev 2022. ‡ = alternative eigen formulation. † = Greedy template reduced with a Schur complement.
  - The Syzygy column differs from the CVPR table, apparently a re-run; for example E+f 6pt is 19×28 here vs 21×30 in CVPR.
  - Internal inconsistency: quivers is 37×71 in Table 1, but the Fig. 4 histogram legend says "Our (61x86)".
- Table 2 (5K random instances; mean log10 residual and fail%) — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2):
  - The resultant solvers have 0% failures on most problems.
  - λ+E+λ 6pt: −11.86, 0.06%, vs Syzygy −7.50, 5.16%, Heuristic 0.40% and Greedy 0%.
  - Unsynch. relative pose: −13.59, 0% vs Syzygy 3.43%, Heuristic 0.36% and Greedy 0.44%.
  - Optimal pose 2pt v2: the resultant solver is worse at 0.52% vs 0 to 0.05% for the others.
  - Optimal PnP Cayley: mean −12.36 vs Greedy −14.30.
- The conclusion admits that "for some minimal problems, it leads to solvers with larger eigenvalue problems but performing a smaller matrix inverse" — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Code: the generator is "publicly available at (Bhayani et al, 2020a)" = github.com/snehalbhayani/aut_gen_sparse_res_solver. The repo README states it needs MATLAB R2018a+ and Maple 2018+ (through the MATLAB symbolic Maple toolbox) for the offline stage — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2); local repo README (/home/user/aut_gen_sparse_res_solver/README.md)

### Inferences
- The equivalence result suggests that the smallest-template question is the same question for both families: "find the smallest Ĉ whose elimination yields a multiplication/Schur matrix". The difference lies in the search heuristics (polytope subsets vs. GB, syzygy or greedy), not in the online numerics.

### Gaps
- No wall-clock timings were reported in the JMIV tables, and the full Springer version was not compared line by line with arXiv v2.

## Q4. Other resultant / eigenvalue-based solver work (CV and computational algebra)

### Takeaway
Before 2017, resultant solvers in CV were hand-built hidden-variable PEPs for simple problems (5pt, 6pt; Kukelova et al., Hartley & Li). Heikkilä (ICCV 2017) brought in polytope-subset basis selection, and the Bhayani line (2020–2024) automated it. In algebra, Canny–Emiris and Emiris's general solver, Telen–Mourrain–Van Barel truncated normal forms and Bender–Telen toric and "yet another" eigenvalue methods give the theory. Newest (2026): a hidden-variable solver with no inversion, using FFT interpolation of the determinant.

### Cited Findings
- Kukelova, Bujnak, Pajdla, BMVC 2008 — [BMVC 2008](https://www.bmva-archive.org.uk/bmvc/2008/papers/164.pdf):
  - 5-pt relative pose is a cubic PEP (z³C3 + z²C2 + zC1 + C0)v = 0 with 10×10 matrices C_i, built from a 10×20 matrix M. It gives 30 eigenvalues for the 10 real solutions.
  - 6-pt focal is a quadratic eigenvalue problem (QEP) in w = f⁻² with 10×10 matrices, built from a 10×30 matrix. It gives 20 eigenvalues for the 15 solutions.
  - Both are solved with MATLAB polyeig and are described as "fast and somewhat more stable than existing solutions".
- The PAMI 2012 version (vol. 34, no. 7, pp. 1381–1393) adds a "resultant-based method for transforming a system of polynomial equations to a polynomial eigenvalue problem" — [dblp](https://dblp.org/rec/journals/pami/KukelovaBP12)
- Hartley & Li, "An efficient hidden variable approach to minimal-case camera motion estimation", IEEE TPAMI 34(12):2303–2314, 2012 — cited in [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Heikkilä, "Using sparse elimination for solving minimal problems in computer vision", ICCV 2017, pp. 76–84. It offers a new basis-selection algorithm for PEP solvers, more compact than the state of the art, and solves two homography-based planar self-calibration problems (focal length from two homographies of an unknown plane) — [CVF](https://openaccess.thecvf.com/content_iccv_2017/html/Heikkila_Using_Sparse_Elimination_ICCV_2017_paper.html); [Oulu repo](https://urn.fi/URN:NBN:fi-fe202002256379)
- Emiris, "A general solver based on sparse resultants" (arXiv 1201.5810, 2012): reviews the Canny–Emiris constructions, reduces root-finding to an eigenproblem with "a new reduction method … that does not increase the problem's dimension", and gives applications in vision, robotics and structural biology — [arXiv 1201.5810](https://arxiv.org/abs/1201.5810)
  - Emiris applied the method to the 5-point problem, "however the final solver isn't particularly efficient" — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
  - Bhayani et al. state there is "no publicly available code for the method [Emiris]" — [arXiv 1912.10268 suppl.](https://arxiv.org/abs/1912.10268)
- Canny & Emiris, "A subdivision-based algorithm for the sparse resultant", J. ACM 47(3):417–451, 2000; Emiris & Canny, "A practical method for the sparse resultant", ISSAC'93, pp. 183–192 — cited in [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- Telen, Mourrain, Van Barel, "Solving polynomial systems via truncated normal forms", SIAM J. Matrix Anal. Appl. 39(3):1421–1447, 2018. Follow-up: Mourrain, Telen, Van Barel, "Truncated normal forms for solving polynomial systems: generalized and efficient algorithms", J. Symbolic Comput. 102:63–85, 2021 (arXiv 1803.07974) — [arXiv 1803.07974](https://www.arxiv.org/abs/1803.07974)
- Bender & Telen, "Toric eigenvalue methods for solving sparse polynomial systems", Math. Comp. 91 (2022) 2397–2429 (arXiv 2006.10654). It gives an eigenvalue theorem for zero-dimensional subschemes of compact toric varieties (homogeneous ideals in the Cox ring), handles arbitrary multiplicities, and gives regularity and complexity bounds for weighted, multihomogeneous and unmixed systems — [AMS](https://www.ams.org/journals/mcom/2022-91-337/S0025-5718-2022-03744-3/home.html); [arXiv 2006.10654](https://arxiv.org/abs/2006.10654)
- Bender & Telen, "Yet another eigenvalue algorithm for solving polynomial systems" (arXiv 2105.08472, v2 Feb 2022; EMS Press chapter). For overdetermined systems it reduces to an eigenproblem "simpler and considerably faster" than earlier methods and can beat homotopy continuation. Julia code: EigenvalueSolver.jl — [arXiv 2105.08472](https://arxiv.org/abs/2105.08472)
- Wu, Bhayani, Heikkilä, "Solving minimal problems without matrix inversion using FFT-based interpolation" (arXiv 2605.06572, May 2026, listed as CVPR 2026) — [arXiv 2605.06572](https://arxiv.org/abs/2605.06572):
  - It builds on sparse hidden-variable resultants.
  - It reconstructs det M(x_h) from sampled evaluations with inverse FFT instead of symbolic expansion, then solves for the hidden variable.
  - The other unknowns come from rank-1-deficient submatrices via Cramer's rule, with a GCD-based criterion for robustness to noise.
  - Claims: "strong numerical stability and competitive runtime, particularly for small-scale problems".

### Inferences
- The toric and TNF methods are the closest mathematical framework for x_i − λ Schur-complement solvers: both choose a basis from a Macaulay-type matrix and compute a multiplication map. Neither has been packaged as an offline/online CV template generator, which matches the JMIV paper's remark.

### Gaps
- I did not verify the Heikkilä ICCV 2017 template sizes or Telen's thesis contents; the CVF page returned 403.
- I found no published CV benchmark that uses EigenvalueSolver.jl or TNF on minimal problems.
- The FFT paper's quantitative results were not retrieved.

## Q5. Known weaknesses, and reported or available fixes

### Takeaway
The main weaknesses:
- the search is heuristic and randomized: δ shifts, random column-removal order, random upper-block row removal, rank tests at random coefficient values;
- polytopes ignore algebraic constraints on the coefficients;
- the eigenproblem can be larger than the solution count, giving spurious eigenvalues;
- the inverse of Â12 is unstable near degeneracies;
- the hidden-variable GEP has parasitic 0/∞ eigenvalues;
- the reference generator depends on Maple.

Reported fixes so far are the ICPR hidden-variable variant (stability), the JMIV iterative subset scheme, the 2026 FFT interpolation, and a local deterministic reimplementation.

### Cited Findings
- Randomness in the CVPR method: "select columns of M one by one in a random order"; extra upper-block rows are "randomly chosen" — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
- In JMIV, rank conditions are checked "for a random value of u0 ∈ C" — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Sparse resultants were "less successful … primarily because the polytopes do not respect the constraints on the coefficients" — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Failures and sizes:
  - The generator failed on E angle+4pt in both the CVPR and JMIV tables — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
  - λ+E+λ 6pt had 14.66% failures in the CVPR version (the others were worse) — [arXiv 1912.10268](https://arxiv.org/abs/1912.10268)
  - Eigenproblems larger than the solution count are flagged in blue in the JMIV table — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Instability of the Schur complement: the CVPR 2020 solvers "can become unstable if the matrix that is to be inverted becomes close to singular … especially in close-to-degenerate configurations" — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
- Hidden-variable GEP issues: eigenvectors that ignore monomial dependencies, parasitic 0/∞ eigenvalues, and rank-deficient Heikkilä bases (35→32 for P5Pr, 82→76 for P6Pfr). The fixes are explicit rank tests and row-column removal — [arXiv 2007.10100](https://arxiv.org/abs/2007.10100)
- Spurious solutions at infinity or on coordinate hyperplanes: the papers discuss only parasitic zero or ∞ eigenvalues. I found no explicit saturation or torus-specific treatment in the Bhayani papers (see Gaps).
- Tooling: the reference generator requires MATLAB plus Maple — local [README](/home/user/aut_gen_sparse_res_solver/README.md)
- Local fix (this repository, not peer-reviewed): `symgen/` is a Python/SymPy reimplementation. It replaces the random tests with exact rank tests over GF(p), p = 2097143, and uses the exact validity condition rank Top[:, G] = |G|. It finds the minimum eigen size per row set in one matroid-greedy pass, uses deterministic row and column removal orders, and needs no Maple. For relpose_6p_focal it reports a minimum eigen size of 18 across bases of 30–517 monomials (15 real solutions plus 3 spurious) — local [symgen/README.md](/home/user/aut_gen_sparse_res_solver/symgen/README.md)
- FFT interpolation avoids the matrix inversion entirely — [arXiv 2605.06572](https://arxiv.org/abs/2605.06572)

### Inferences
- Because rank tests are run at one random complex or real point, a bad sample could accept or reject a template wrongly. Finite-field exact tests (as in symgen) are the standard remedy for this.
- Spurious eigenvalues from the extra-polynomial construction correspond to points of the relaxed system (for example at infinity or with x_i = 0 on the toric boundary). The toric eigenvalue theory of Bender & Telen is the natural tool for characterizing them, but no CV paper applying it was found.

### Gaps
- No paper reports run-to-run variance of template sizes for the random heuristics.
- I found no reliable source on applying saturation (Larsson ICCV 2017, "Polynomial solvers for saturated ideals") inside the sparse resultant generator; only the local repo commits mention a "saturated 6p focal" experiment.

## Q6. Code availability

### Takeaway
- The Bhayani extra-polynomial generator is open source (MATLAB + Maple).
- The Julia toric and overdetermined eigenvalue solver is EigenvalueSolver.jl.
- Emiris's CV solver has no public code according to Bhayani et al.
- The ICPR hidden-variable and 2026 FFT code was not located.

### Cited Findings
- github.com/snehalbhayani/aut_gen_sparse_res_solver: the generator cited in the JMIV paper. It takes `problems/*.m` configs with fields numOfCoeff, numOfVars, hiddenVarNum, polyComb, noOfRowsToReduce and heurisiticTemplatesize, and outputs `solvers/<name>/solver.m` — [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2); local [README](/home/user/aut_gen_sparse_res_solver/README.md)
- EigenvalueSolver.jl (Bender & Telen) — [arXiv 2105.08472](https://arxiv.org/abs/2105.08472)
- Emiris u-resultant solver: "no publicly available code" — [arXiv 1912.10268 suppl.](https://arxiv.org/abs/1912.10268)

### Gaps
- I did not confirm a public repository for the TNF code (Telen's Python/Julia TNF implementation), the ICPR 2020 hidden-variable generator, or the 2026 FFT solver.
