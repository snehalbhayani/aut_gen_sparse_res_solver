# Recent (2022-2026) automatic minimal solver generation: SOTA, benchmarks, tooling, open problems

Notes compiled 2026-11-26. Template sizes are rows x columns of the elimination template / resultant matrix. For a Gröbner/action-matrix template, columns minus rows is usually the size of the eigenvalue problem (this holds for every row of Martyushev-Bhayani-Pajdla 2023 Table, see Q2). "d" / "sols" = number of solutions (the quotient-ring dimension).

## Q1. Newest papers (2022-2026) on automatic solver generation and new solving paradigms

### Takeaway
After 2022 the elimination-template line is led by Martyushev (greedy parameter search, CVPR 2022, then a Laurent-polynomial generator in 2023/24 that gives the smallest templates on most benchmarks) and by Bhayani-Heikkilä-Kukelova (sparse-resultant generator and its formal equivalence with action-matrix solvers, JMIV 2024). The alternatives are learned or simulated single-path homotopy continuation (Hruby et al. CVPR 2022, best paper; Simulator HC 2024/25), learned online template selection, and, newest, inversion-free FFT-interpolation hidden-variable resultant solvers (CVPR 2026).

### Cited Findings
**Elimination templates / action matrix (Gröbner-free)**
- Martyushev, Vrablikova, Pajdla, "Optimizing Elimination Templates by Greedy Parameter Search", CVPR 2022 (arXiv 2203.14901). The method builds a partial affine parameterization of elimination templates for systems with finitely many solutions, then runs a greedy heuristic over the parameters. Tested on 34 minimal problems, the templates were the same size or smaller than the SOTA, and 2.1, 2.5, 3.8 and 6.6 times smaller on some hard cases — [arXiv 2203.14901](https://arxiv.org/abs/2203.14901); [ar5iv](https://ar5iv.labs.arxiv.org/html/2203.14901)
  - Problem 33 (rel. pose fλ+E+fλ 7pt, 68 sols): 209x277 template, 8.5 ms in Matlab, compared with 886x1011 at ~400 ms for Jiang et al. and 581x659 at ~16 ms for syzygy (which was unstable and needed column pivoting). Oskarsson's C++-core solver was ~3x faster (2.9 ms). About 2% of trials failed. Computing the full Gröbner fan did not finish in ~24 h, so 1000 random weight orders plus 500 heuristic bases were sampled instead — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)
  - The parameter θ_ij is restricted to scalars, not polynomials, so the parameterization is only partial. There is no optimality guarantee. Schur-complement reduction applies only to polynomials with fixed coefficients. Symmetries in problems 15 (2-fold) and 32 (4-fold) were not exploited — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)
- Martyushev, Bhayani, Pajdla, "Automatic Solver Generator for Systems of Laurent Polynomial Equations", arXiv 2307.00320 (July 2023). A later Martyushev paper cites it as a 2024 publication, but I could not confirm the venue. It gives a practical algorithm to test whether a set of Laurent polynomials is enough to build a template, and the generator built on it handles ideals with positive-dimensional components and finds partial p-fold symmetries automatically. It reports new solvers for optimal 3-view triangulation, semi-generalized hybrid pose and time-of-arrival (ToA) self-calibration. The generator is in Maple and the solvers in Matlab, at github.com/martyushev/EliminationTemplates — [arXiv abs](https://arxiv.org/abs/2307.00320); [arXiv html](https://arxiv.org/html/2307.00320)
  - Operating notes: the iteration cap N=10 is enough for all tractable systems. For rel. pose from 9 lines (116 roots), the generator found a 2163x2616 template in ~20 min, and reducing it to 1610x1726 took ~13 h. The paper argues that template size alone is not an adequate efficiency measure, because a smaller template can carry a larger eigenproblem — [arXiv html 2307.00320](https://arxiv.org/html/2307.00320)
- Bhayani, Heikkilä, Kukelova, "Sparse resultant based minimal solvers in computer vision and their connection with the action matrix", J. Math. Imaging and Vision 2024 (arXiv 2301.06443). The method adds an extra polynomial of a special form, uses Newton polytopes / Minkowski sums to build the sparse resultant matrix, reduces it by Schur complement, and converts it to a regular eigenproblem. Section 4 gives step-by-step conversions in both directions between action-matrix and sparse-resultant solvers and the conditions under which they are equivalent. The code is at github.com/snehalbhayani/aut_gen_sparse_res_solver — [Springer JMIV](https://link.springer.com/article/10.1007/s10851-024-01182-1); [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)
- Martyushev, "Forward kinematics of a general Stewart-Gough platform by elimination templates", Mechanism and Machine Theory 215:106170 (2025), arXiv 2505.00634. This is a robotics application of the generator: a 40-solution problem, a sparse 293x362 template, then PLU to a 69x69 generalized eigenproblem pair. The paper says MATLAB/Julia/Python code "will be made publicly available" — [arXiv 2505.00634](https://arxiv.org/abs/2505.00634)
- Macaulay2 package **EliminationTemplates** (v1, dated 2026-04-30; authors Batavia, Chen, Shen, Chlopecki, Duff, Huang, Li). It provides `eliminationTemplate`, `getTemplate`, `getTemplateMatrix`, `getActionMatrix` and `templateSolve`. It references Martyushev et al. CVPR 2022 and Larsson et al. CVPR 2017 (syzygy) — [Macaulay2 docs](https://www.macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EliminationTemplates/html/index.html)

**New solving paradigms**
- Wu, Bhayani, Heikkilä, "Solving Minimal Problems Without Matrix Inversion Using FFT-Based Interpolation", CVPR 2026 (arXiv 2605.06572, May 2026). It uses sparse hidden-variable resultants with no online matrix inversion. The determinant polynomial in the hidden variable is reconstructed from sampled evaluations by inverse FFT interpolation, the remaining unknowns come from rank-1-deficient submatrices via Cramer's rule, and a GCD-based criterion makes this robust to noise. The authors report strong numerical stability and competitive runtime, especially on small problems — [arXiv 2605.06572](https://arxiv.org/abs/2605.06572)
- Hruby, Duff, Leykin, Pajdla, "Learning to Solve Hard Minimal Problems", CVPR 2022 (best paper). It is a pick-and-solve method: a learned picker chooses one start problem-solution pair and a single real homotopy-continuation (HC) path is tracked. Average solve time is under 70 µs per problem — [CVF](https://openaccess.thecvf.com/content/CVPR2022/html/Hruby_Learning_To_Solve_Hard_Minimal_Problems_CVPR_2022_paper.html); [arXiv 2112.03424](https://www.arxiv.org/pdf/2112.03424)
- Zhang, Dai, Xu, Kneip, "Simulator HC: Regression-based Online Simulation of Starting Problem-Solution Pairs for HC in Geometric Vision", arXiv 2411.03745 (v2, March 2025). A regressor predicts a rough solution, an online simulator builds a matching start problem, and one path is tracked back. Applied to generalized camera resectioning and a new generalized relative pose + scale solver. The paper states explicitly that some problems stay intractable for generated templates because of constraint degree or number of unknowns — [arXiv 2411.03745](https://arxiv.org/abs/2411.03745)
- Xu, Hu, Tsakiris, Kneip, "Online Stability Improvement of Groebner Basis Solvers using Deep Learning" (3DV 2019, arXiv posted Jan 2024). Variable reordering amounts to permuting the columns of the coefficient matrix, so one template can be reused in several ways, and a classifier picks the most stable one online — [arXiv 2401.09328](https://arxiv.org/abs/2401.09328)
- GPU-based HC (arXiv 2112.03444) targets problems whose templates are impractical. For one problem the full template is 1866x1975. It is applied to 4-view triangulation and trifocal pose with unknown focal length — [arXiv 2112.03444](https://arxiv.org/html/2112.03444v2)
- Theory and survey context: "Snapshot of Algebraic Vision" (Kileel & Kohn, 2022) says symbolic methods are the most popular way to build minimal solvers, while HC gave the fastest solver for a minimal SfM problem — [arXiv 2210.11443](https://arxiv.org/html/2210.11443). Also "Numerically Computing Galois Groups of Minimal Problems" (arXiv 2507.10407, 2025) — [arXiv](https://arxiv.org/pdf/2507.10407), and PLMP point-line minimal problems for projective SfM (arXiv 2503.04351) — [arXiv](https://arxiv.org/pdf/2503.04351). I saw only listings and reference lists for these two, not their content.
- Generator-using application paper: "Combining Absolute and Semi-Generalized Relative Poses for Visual Localization" (arXiv 2409.14269) — [arXiv](https://arxiv.org/pdf/2409.14269). It came up in searches for generator-built solvers; I did not verify its details.

**Pre-2022 baselines referenced by the above**
- Kukelova, Bujnak, Pajdla, automatic generator, ECCV 2008 — [PDF](https://cmp.felk.cvut.cz/~kukelova/publications/Kukelova-etal-ECCV-2008.pdf). Larsson, Åström, Oskarsson, syzygy-based reduction, CVPR 2017 (template generation in Macaulay2). That paper reports that smaller templates are usually faster but that accuracy does not always correlate with size — [CVF](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf). Larsson et al., "Beyond Gröbner Basis: Basis Selection for Minimal Solvers", CVPR 2018 (Gröbner-fan enumeration plus heuristic non-standard bases) — [arXiv 1803.04360](https://arxiv.org/pdf/1803.04360). Kukelova, Kileel, Sturmfels, Pajdla, "A clever elimination strategy", CVPR 2017: E+fλ 7pt went from 200x231 to 51x70 — [arXiv 1703.05289](https://arxiv.org/pdf/1703.05289). Bhayani et al. sparse-resultant conference version — [arXiv 1912.10268](https://arxiv.org/pdf/1912.10268). GAPS generator (2020) was mentioned in search results — [arXiv 2307.00320 html](https://arxiv.org/html/2307.00320).

### Inferences
- From 2022 to 2024 the main gains came from enlarging the search space: affine parameterization (2022), then Laurent monomials, saturation-free handling of positive-dimensional components and partial symmetries (2023). These cut template sizes 1.1x to 3x relative to Greedy on hard problems (e.g., unsynchronized relative pose 139x155 to 59x79; refractive P6Pf 126x162 to 40x61).
- Between 2024 and 2026 activity moved away from template size towards the online solve: no-inversion FFT resultants, learned or simulated single-path HC, and learned template selection. This suggests template-size minimization is reaching diminishing returns on the standard benchmark set.
- I found no LLM- or AI-agent-based symbolic solver generation in the literature. "Learning" is used only for start-point picking, regression and stability classification.

### Gaps
- I could not confirm the publication venue of arXiv 2307.00320. A citing paper dates it 2024, which suggests a journal version.
- I did not find a separate "Gröbner-free elimination templates" or "template rearrangement" paper under those names. The closest matches are the Greedy (2022) and Laurent (2023) papers.
- Runtimes and success rates for Simulator HC were not in its abstract, and I did not read the full text.
- I did not verify a Ding/Kukelova "fundamental matrix from relative depths" paper or the 2024-2026 multi-camera affine-correspondence papers.
- I could not reach the CVPR/ICCV/ECCV 2022-2025 tutorial pages, so I found no evidence of a 2024/25 minimal-problems tutorial.

## Q2. Consolidated benchmark table (template sizes, solution counts, runtimes)

### Takeaway
On the standard 34-to-36-problem benchmark, Martyushev-Bhayani-Pajdla's Laurent generator (2023) gives the smallest templates on almost all problems. It sometimes does so by accepting extra (spurious or symmetric) roots, i.e. columns minus rows greater than d. Greedy (2022) and Bhayani's sparse-resultant generator (JMIV 2024) are next, and syzygy (2017) and BeyondGB (2018) are clearly larger on hard problems.

### Cited Findings
Sources: [Greedy 2022 Tables 1-2](https://ar5iv.labs.arxiv.org/html/2203.14901), [Laurent 2023 table](https://arxiv.org/html/2307.00320), [Bhayani et al. JMIV 2024 Table 1](https://arxiv.org/pdf/2301.06443v2).
- Greedy = Martyushev et al. CVPR 2022, smallest of its std/nstd basis.
- Laurent = Martyushev-Bhayani-Pajdla 2023; the number in parentheses is the eigenproblem size.
- SparseR (Greedy table) = the 2020 Bhayani resultant solver; SpRes-JMIV = Bhayani-Heikkilä-Kukelova JMIV 2024.
- \* = Schur-complement-reduced; † = positive-dimensional formulation; ‡ = Laurent formulation.
- Timings are Matlab, in ms.

| Problem | d | Syzygy 2017 | BeyondGB 2018 (std / nstd) | SparseR 2020 | SpRes-JMIV 2024 | Greedy 2022 | Laurent 2023 (roots) | Laurent ms / Greedy ms |
|---|---|---|---|---|---|---|---|---|
| Rel. pose F+λ 8pt | 8 | 11x19 | 11x19 / 7x15 | 7x16 | 7x16 | **7x15** | — | — |
| Rel. pose E+f 6pt | 9 | 21x30 (JMIV: 19x28) | 11x20 / 11x20 | 11x20 | 11x20 | **11x20** | — | — |
| Rel. pose f+E+f 6pt | 15 | 31x46 (JMIV: 25x40) | 31x46 / 21x40 | 12x30 | **11x26** | 11x26 | — | — |
| Rel. pose E+λ 6pt | 26 | 34x60 | 34x60 / 14x40 | 14x40 | 14x40 | **14x40** | — | — |
| Stitching fλ+R+fλ 3pt | 18 | 48x66 | 48x66 / 18x36 | 18x36 | **8x31** | 18x36 | 6x30 (24) | 0.2 / 0.2 |
| Abs. pose P4P+fr | 16 | 140x156 | 54x70 | 52x68 | 52x68 | 52x68 | **42x60 (18)** | 0.2 / 0.3 |
| P4P+fr (elim. f) | 12 | 28x40 | 28x40 | 28x40 | 28x40 | 28x40 | — | — |
| Rel. pose λ+E+λ 6pt | 52 | 149x201 | — / 53x105 | 39x95 | 39x95 | **39x95** | — | — |
| Rel. pose λ1+F+λ2 9pt (two-sided radial) | 24 | 165x189 | 87x111 / 87x111 (JMIV lists GFan 105x129, heuristic 73x97) | 90x117 | 77x104 | 76x100 | **73x97 (24)** | 0.3 / 0.4 |
| Rel. pose E+fλ 7pt | 19 | 185x204 (JMIV: 181x200) | 69x88 | 61x80 | **34x56** | 55x74 | 34x56 (22) | 0.2 / 0.6 |
| E+fλ 7pt (elim. λ) | 19 | 52x71 | 37x56 / 24x43 | 22x41 | 22x41 | **22x41** | — | — |
| E+fλ 7pt (elim. fλ) | 19 | 51x70 | 51x70 | 51x70 | 51x70 | 51x70 | — | — |
| Rolling-shutter pose | 8 | 47x55 | 47x55 | 47x55 | 40x52 | 47x55 | **40x52 (12)** | 0.2 / 0.2 |
| Triangulation, satellite | 27 | 88x115 | 88x115 | 87x114 | 87x114 | 87x114 | **74x104 (30)** | 0.4 / 0.8 |
| Abs. pose refractive P5P | 16 | 240x256 | 112x128 / 199x215 | 68x93 | 57x73 | 57x73 | **38x58 (20)** | 0.2 / 0.3 |
| Abs. pose quivers | 20 | 169x189 | — / 68x88 | 68x92 | 37x71 | 65x85 | 56x80 (24) | 0.2 / 0.3 |
| Unsynchronized rel. pose | 16 | 159x175 | — / 299x315 | 150x168 | 59x80 | 139x155 | **59x79 (20)** | 0.3 / 0.6 |
| Optimal PnP (Hesch) | 27 | 88x115 | 88x115 | 87x114 | 87x114 | **87x114** | — | — |
| Optimal PnP (Cayley) | 40 | 118x158 | 118x158 | 118x158 | 118x158 | 118x158 | **78x131 (53)** | 1.1 / 0.9 |
| Optimal pose w/ direction 2pt v2 | 24 | 192x216 | — / 192x216 | 176x200 | 112x146 | 139x163* | **87x120 (33)** | 0.5 / 1.2 |
| Rel. pose E+angle 4pt | 20 | 246x276 | — / 183x249 | — | — | 99x119* | — | — |
| Rel. pose λ+F+λ 8pt | 16 | 32x48 | — | — | — | 31x47 | 19x39 (20) | 0.2 / 0.2 |
| P3.5P + focal | 10 | 20x30 (orig. 20x43) | — | — | — | 18x28 | **12x26 (14)** | 0.1 / 0.1 |
| Gen. P4P + scale | 8 | 47x55 | — | — | — | 47x55 | — | — |
| Rel. pose E+angle 4pt v2 | 20 | 36x56 | — | — | — | 16x36* | **12x36 (24)** | 0.2 / 0.2 |
| Gen. rel. pose E+angle 5pt | 44 | 317x361 | — | — | — | 37x81* | **34x79 (45)** | 0.8 / 0.9 |
| Rel. pose E+fuv+angle 7pt | 6 | 66x72 | — | — | — | 40x46 (orig. cascade 0.4 ms vs 0.6 ms) | 26x33† (7) | 0.1 / 0.2 |
| Rolling shutter R6P | 20 | 204x224 (orig. 196x216) | — | — | — | 120x140 | **66x92 (26)** | 0.3 / 0.9 |
| Opt. pose w/ dir 4pt | 28 | 203x231 | — | — | — | 134x162* | **105x138 (33)** | 0.6 / 1.1 |
| Opt. pose w/ dir 3pt | 48 | 544x592 (orig. 1260x1278) | — | — | — | 385x433* | **297x356 (59)** | 3.0 / 6.8 |
| L2 3-view triangulation (relaxed) | 31 | 231x262 | — | — | — | 217x248 | **190x227 (37)** | 0.6 / 1.4 |
| Refractive P6P + focal | 36 | 636x654 (orig. 648x917) | — | — | — | 126x162 | **40x61 (21, symmetry)** | 0.2 / 0.8 |
| Rel. pose fλ+E+fλ 7pt | 68 | 581x659 (orig. 886x1011) | — | — | — | 209x277 | **138x210 (72)** | 2.3 / 3.0 |
| Rel. pose λ1+E+λ2 7pt | 76 | — | — | — | — | 436x512 | **204x289 (85)** | 4.2 / 27.9 |
| Gen. rel. pose 6pt (Stewenius) | 64 | — | — | — | — | 99x163 | **78x155 (77)** | 2.5 / 2.5 |
| Gen. rel. pose + scale 7pt | 140 | 144x284 | — | — | — | 144x284 | 130x275 (145) | 9.7 / 9.2 |
| Optimal PnP (Zheng-type, [57] in Laurent paper) | 40 | — | — | — | — | 272x312 | **120x175 (55)** | 1.5 / 5.9 |
| Weak PnP | 16 | — | — | — | — | 108x124 | 85x107 (22); with ×2 symmetry 10x20 (10) | 0.3 / 0.6 |
| Rel. pose 9 lines | 116 | — | — | — | — | failed | 1610x1726 | 205 ms |
| Vanishing point est. | 40 | — | — | — | — | 343x383 | 136x194† (58) | 2.1 / 13.8 |
| ToA (4,6) | 38 | — | — | — | — | 863x901 | 427x475† (48) | 5.2 / 17.1 |
| ToA (5,5) | — | — | — | — | — | failed | 772x832† (60) | 17.8 |
| Semi-gen. H51f | — | — | — | — | — | failed | 740x811† (71) | 43.4 |
| Semi-gen. H32f | 26 | — | — | — | — | 258x284 | 174x222† (48) | 2.8 / 5.7 |
| Optimal 3-view triangulation | 47 | — | — | — | — | 401x448 | 69x127‡ (58) | 1.1 / 5.7 |

- Stability (JMIV 2024 Table 2, log10 median residual, noise-free): sparse-resultant solvers mostly reach -12 to -15 with 0% failures. Examples: λ+E+λ 6pt failure rate 0.06% vs 5.16% for syzygy; unsynchronized rel. pose 0% vs 3.43% for syzygy — [arXiv 2301.06443](https://arxiv.org/pdf/2301.06443v2)
- Hand-made solvers listed in the JMIV table's "Author" column include Kuang 2014, Bujnak 2009/2010, Kukelova 2008, Naroditsky 2011, Saurer 2015 and Albl 2017. The "Original" column there gives, for example, P4Pfr 136x152 and E+fλ 7pt 200x231 — [arXiv 2301.06443](https://arxiv.org/pdf/2301.06443v2)
- An HC baseline: Macaulay2 complex HC takes about 105 µs on a benchmark problem versus 5 µs for a dedicated implementation. The pick-and-solve HC solver averages under 70 µs — [CVPR 2022 Hruby et al.](https://openaccess.thecvf.com/content/CVPR2022/html/Hruby_Learning_To_Solve_Hard_Minimal_Problems_CVPR_2022_paper.html)

### Inferences
- Bold marks the smallest entry per row across the cited tables. When comparing, count the eigenproblem size as well: the Laurent generator often trades rows for extra roots (e.g., Optimal PnP Cayley 78x131 has 53 roots against 40 true solutions, and is slightly slower at 1.1 ms vs 0.9 ms).
- The 5pt relative pose (10 sols) and P3P/P4Pf rows are not in the 2022-2024 benchmark tables, because they are saturated and easy. Commonly quoted values (10x20 Gröbner template for 5pt; P4Pf with 10 sols) are not covered by sources fetched here.

### Gaps
- The JMIV-table text extraction for the 9pt λ1+F+λ2 row (heuristic 73x97, GFan 105x129) conflicts with Greedy 2022's BeyondGB column (87x111). This may be a column misalignment in pdftotext output or re-run numbers; check against the published PDF.
- The Syzygy column in JMIV for E+f 6pt (19x28) and f+E+f (25x40) differs from Greedy 2022's syzygy column (21x30, 31x46). The tables likely use different formulations or versions.
- I found no sourced template sizes for 5pt E, P4Pf, hand-eye, refractive relative pose or homography-with-radial-distortion. Larsson's CVPR 2017 PDF returned HTTP 403.
- The Greedy supplementary runtime table (Sec. 13) was not retrieved.

## Q3. Open-source tooling ecosystem

### Takeaway
The working generators are Martyushev's EliminationTemplates (Maple generator, Matlab solvers), Bhayani's sparse-resultant generator (MIT; Matlab, with Maple/Macaulay2 parts in the repo), Larsson's autogen (Macaulay2/Matlab), the new Macaulay2 EliminationTemplates package (2026), and the classic Kukelova generator. On the HC side there are HomotopyContinuation.jl, MINUS and the Hruby learned-HC code. Most symbolic generators depend on proprietary Maple or MATLAB.

### Cited Findings
- **martyushev/EliminationTemplates** (GitHub): Maple implementation of the Laurent generator plus Matlab solvers for all benchmark problems in the 2022 and 2023 papers. The 2022 paper also mentions Gfan (Gröbner fans) and Macaulay2 or Maple for computing H0 — [arXiv html 2307.00320](https://arxiv.org/html/2307.00320); [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)
- **snehalbhayani/aut_gen_sparse_res_solver**: "An automatic generator for sparse resultant based solvers for minimal problems in computer vision", MIT license, last pushed 2026-11-26 (GitHub API). The repo includes `generate_solver.m`, `polysolver.map` (Maple), `convex.mla` and `multires03.txt` (Maple multires resultant code), and `symgen` — local checkout and [JMIV paper](https://arxiv.org/pdf/2301.06443v2)
- **Macaulay2 EliminationTemplates** v1 (2026), which implements the Greedy-2022 and syzygy-2017 approaches in open-source Macaulay2. The docs page states no license — [M2 docs](https://www.macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EliminationTemplates/html/index.html)
- **Kukelova automatic generator (gbsolver, ECCV 2008)** and solver collections (peig5pt/peig6pt, P4Pf, P4Pfr, R6P, refractive) are linked from the CTU AAG "Minimal Problems in Computer Vision" index. The index lists papers up to 2019 and has no tutorials or open-problems section — [aag.ciirc.cvut.cz/minimal](https://aag.ciirc.cvut.cz/minimal/)
- **Syzygy-based reduction (Larsson 2017)**: template generation implemented in Macaulay2 — [arXiv search summary of CVPR 2017](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- **Learned HC (Hruby et al. 2022)**: the supplementary material says the code will be released on GitHub — [CVPR 2022 supp](https://openaccess.thecvf.com/content/CVPR2022/supplemental/Hruby_Learning_To_Solve_CVPR_2022_supplemental.pdf)
- **GPU-HC**: a generic GPU homotopy solver for vision problems — [arXiv 2112.03444](https://arxiv.org/html/2112.03444v2)
- **Stewart-Gough (2025)**: MATLAB/Julia/Python implementations were promised — [arXiv 2505.00634](https://arxiv.org/abs/2505.00634)

### Inferences
- The Maple dependency of the two strongest generators (Martyushev; parts of Bhayani's) is a practical barrier. The 2026 Macaulay2 package and Julia/Python code from Martyushev's 2025 paper point towards open-source replacements.
- No generator emits SIMD or optimized C++ in a way documented in the sources fetched here. Greedy-2022 solvers are pure Matlab and ~3x slower than a C++-core solver on the same problem.

### Gaps
- I could not check licenses or contents of PoseLib, autogen/autogen_v2, HomotopyContinuation.jl, MINUS, learning_minimal or MonodromySolver, because GitHub API access to those repos was blocked in this session. From memory, which is unverified, PoseLib is BSD-3 C++ and HomotopyContinuation.jl is MIT; confirm before citing.
- I did not confirm whether a Python/Julia port of autogen exists.

## Q4. Open research problems stated by authors

### Takeaway
The stated open problems are: no optimality guarantee for template size (generators are heuristic), size being a poor proxy for speed and stability, incomplete parameterizations (scalar vs polynomial multipliers), exploiting symmetries and positive-dimensional components automatically, unifying action-matrix, resultant and normal-form methods, offline generation cost for large problems, and problems still intractable for templates, which motivates HC.

### Cited Findings
- **Heuristic, no optimality**: the greedy row/column reductions are heuristics without an optimality guarantee. Restricting θ to scalars instead of polynomials leaves the parameterization partial — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)
- **Size is not efficiency**: "the size of a template alone is not an appropriate measure of the efficiency", because a smaller template can be slower when its eigendecomposition is larger — [arXiv 2307.00320](https://arxiv.org/html/2307.00320). Larsson 2017 reports no clear correlation between accuracy and template size — [CVPR 2017](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf). The JMIV paper notes that its resultant solvers sometimes have larger eigenproblems but smaller matrix inverses — [arXiv 2301.06443](https://arxiv.org/pdf/2301.06443v2)
- **Symmetry automation**: Greedy 2022 left 2-fold and 4-fold symmetries unexploited. The Laurent generator finds partial p-fold symmetries automatically, e.g. Weak PnP 108x124 to 10x20 and refractive P6Pf 126x162 to 40x61 — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901); [arXiv 2307.00320](https://arxiv.org/html/2307.00320)
- **Saturation / positive-dimensional components**: the Greedy generator applies only to zero-dimensional ideals, so positive-dimensional components had to be saturated by hand. The Laurent generator removes this requirement — [arXiv 2307.00320](https://arxiv.org/html/2307.00320)
- **Unification and normal forms**: the JMIV conclusion asks under what conditions action-matrix and sparse-resultant solvers coincide, and says that "normal forms have yet not been incorporated in the automatic generators for minimal solvers in computer vision. This may be an interesting future topic to investigate" — [arXiv 2301.06443](https://arxiv.org/pdf/2301.06443v2)
- **Scaling / offline cost**: the full Gröbner fan did not finish in 24 h for the 68-solution problem, so the search fell back to random sampling — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901). For 9-line relative pose, template reduction took ~13 h and the resulting 1610x1726 solver runs in 205 ms. Greedy failed outright on 9 lines, ToA (5,5) and H51f — [arXiv 2307.00320](https://arxiv.org/html/2307.00320)
- **Intractable problems → HC**: some problems remain intractable for templates because of degree or number of unknowns — [arXiv 2411.03745](https://arxiv.org/abs/2411.03745). One full template is 1866x1975, which is impractical — [arXiv 2112.03444](https://arxiv.org/html/2112.03444v2)
- **Online numerical cost and stability**: matrix inversion in the online stage is a computational bottleneck, which motivates inversion-free FFT interpolation — [arXiv 2605.06572](https://arxiv.org/abs/2605.06572). A fixed template's accuracy varies with variable ordering per instance, which motivates learned online selection — [arXiv 2401.09328](https://arxiv.org/abs/2401.09328). Greedy's 68-solution solver failed in ~2% of trials, and syzygy solvers can need column pivoting — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)
- **Spurious solutions**: some formulations give extra roots (e.g., problem 8 has 4 spurious solutions, and many Laurent templates have more roots than d), which must be filtered by residuals — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901); [arXiv 2307.00320](https://arxiv.org/html/2307.00320)
- **Implementation speed**: Matlab-only solvers leave a speed-up on the table compared with C++ — [ar5iv 2203.14901](https://ar5iv.labs.arxiv.org/html/2203.14901)

### Inferences
- Open problems the sources imply but do not state formally: (a) a lower bound or certificate of minimal template size, which no paper attempts; (b) a cost model combining template LU, eigenproblem size and stability, to replace size as the objective; (c) open-source (non-Maple) generators with optimized C++ codegen; (d) systematic handling of degenerate and non-generic configurations inside RANSAC. None of the fetched papers addresses (d).

### Gaps
- I found no workshop or tutorial talk with an explicit open-problems list for 2022-2026; the tutorial pages were not reachable.
- The Kileel-Kohn survey's open-problem list was not extracted in detail.
