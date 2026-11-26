# Laurent templates lead; symgen must close the gap

Through 2026, the smallest elimination templates on the standard computer-vision benchmark come from Martyushev, Bhayani and Pajdla's Laurent-polynomial generator (arXiv 2023). It beats the Greedy search of 2022 and the sparse-resultant generator of the JMIV 2024 paper on most hard problems, for example 9pt two-sided radial distortion at 73×97 and refractive P6P+focal at 40×61. It often gets there by accepting a few spurious eigenvalues. The user's sparse-resultant line is competitive, ties or wins on about half the benchmark (6pt shared focal 11×26, E+fλ 7pt 34×56, unsynchronized pose 59×80), and has the best published stability numbers. Its own new step is the inversion-free FFT hidden-variable solver at CVPR 2026. symgen, the deterministic Python reimplementation, is algorithmically cleaner than the Maple generator: every decision is an exact GF(p) rank test, and it finds the minimal eigen size in one matroid pass. But it only reproduces the CVPR 2020 template sizes (6pt focal 12×30, 9pt 90×117), so it lags the best published templates by 10–25% in rows. Outside templates, homotopy continuation (HC) has taken over problems with more than about 100 solutions: GPU all-path HC runs in milliseconds, and learned single-path HC runs in microseconds at a 25–30% success rate per sample. Galois/monodromy analysis now tells you when a formulation cannot be improved. Of the six open problems in the user's thesis, two have seen real progress since 2022 (unified AM/resultant search, and the GEP-to-inverse spectrum). Two remain open. One, normal forms, has not been touched by any computer-vision generator. The most valuable next steps for symgen are mostly engineering: a Greedy/Laurent-style multiplier search on top of its exact tests, input-formulation reuse, C++ code generation, symmetry detection, automated saturation, and a measured cost model that covers the hidden-variable end of the spectrum.

## Template generators shrank 2–20× in a decade, and the race is now about eigen size and constants

The field started from the Kukelova–Bujnak–Pajdla generator (ECCV 2008). It fixed a GRevLex Gröbner basis over Z_p, expanded monomial multiples, and pruned. Representative templates were 31×46 for 6pt shared focal, 179×203 for 9pt two radial and 238×290 for λ+E+λ 6pt ([Kukelova thesis](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)). Every later advance enlarged the space of templates searched:

- **Syzygy reduction (Larsson et al. 2017)** reduced multiplier vectors modulo the syzygy module. It cut λ+E+λ to 154×210 and relative pose with known direction from 411×489 to 40×57 with symmetry, about a 900× estimated speedup. It was publicly released with a 5.7 s median generator time ([Syzygy17](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)).
- **Saturation (Larsson et al. 2017)** extended the action matrix to I : f^∞ without computing generators online, by lifting the basis with f^N. Vanishing-point estimation dropped from 1.8 s to 5 ms ([Sat17](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf)).
- **Beyond Gröbner bases (2018)** sampled non-standard quotient bases. Examples: E+λ 6pt 34×60 → 14×40, P4Pfr 140×156 → 54×70 ([BeyondGB18](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)).
- **Greedy parameter search (Martyushev, Vrablikova, Pajdla, CVPR 2022)** built an affine parameterization of templates from the syzygy view and searched it greedily. On all 34 problems it matched or beat prior work, and refractive P6P+focal shrank 20× ([Greedy22](https://arxiv.org/abs/2203.14901)).
- **Laurent generator (Martyushev, Bhayani, Pajdla, 2023)** tests directly whether a set of Laurent polynomials suffices for a template. It handles positive-dimensional components and finds partial p-fold symmetries automatically ([Laurent23](https://arxiv.org/abs/2307.00320)). A citing paper dates it 2024, but the notes could not confirm a venue.

The sparse-resultant line ran in parallel. The CVPR 2020 generator added x_i − λ to the system, built a basis from shifted Minkowski sums of Newton-polytope subsets, and reduced the matrix to an eigenproblem with a Schur complement. On most problems it tied or beat the Gröbner generators of the time ([SparseR20](https://arxiv.org/abs/1912.10268)). The JMIV 2024 version added an iterative search over polytope subsets, refreshed the table against Greedy, and proved step-by-step equivalence conditions with action-matrix solvers ([JMIV 2024](https://link.springer.com/article/10.1007/s10851-024-01182-1); [arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)).

### Consolidated benchmark (template rows × columns; columns − rows = eigen size unless roots are given)

Sizes use the Bhayani convention: r×c means an r×r elimination or inverse followed by a (c−r)×(c−r) eigenproblem. The sources are [Syzygy17 Table 1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf), [BeyondGB18 Table 1](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf), [SparseR20 Table 1 and supplement](https://arxiv.org/abs/1912.10268), [JMIV Table 1](https://arxiv.org/abs/2301.06443v2), [Greedy22 Tables 1–2](https://arxiv.org/abs/2203.14901) and [Laurent23 Table 1](https://arxiv.org/abs/2307.00320). The symgen column comes from the symgen experiments in this repo. Bold marks the smallest per row.

| Problem | d | Original | Syzygy 2017 | BeyondGB 2018 best | SparseR 2020 | SpRes JMIV 2024 | Greedy 2022 best | Laurent 2023 (roots) | symgen |
|---|---|---|---|---|---|---|---|---|---|
| F+λ 8pt | 8 | 12×24 | 11×20 | 7×16 | 7×16 | 7×16 | **7×15** | — | — |
| E+f 6pt | 9 | 21×30 | 21×30 ^a | **11×20** | **11×20** | **11×20** | **11×20** | — | — |
| f+E+f 6pt | 15 | 31×46 | 31×50 ^a | 21×40 | 12×30 | **11×26** ^b | **11×26** ^b | — | 12×30 (eig 18 = 15+3); 25×40 saturated (eig 15) |
| E+λ 6pt | 26 | 48×70 | 34×60 | **14×40** | **14×40** | **14×40** | **14×40** | — | — |
| Stitching fλ+R+fλ 3pt | 18 | 54×77 | 48×66 | 18×36 | 18×36 | 8×31 | 18×36 | **6×30 (24)** | — |
| P4P+fr | 16 | 136×152 | 140×156 | 54×70 | 52×68 | 52×68 | 52×68 | **42×60 (18)** | — |
| λ+E+λ 6pt | 52 | 238×290 | 154×210 ^a | 53×115 | 39×95 (56 eig) | **39×95** | **39×95** | — | — |
| λ1+F+λ2 9pt | 24 | 179×203 | 165×200 ^a | 84×117 / 87×111 ^c | 90×117 | 77×104 | 76×100 | **73×97 (24)** | 90×117 |
| E+fλ 7pt | 19 | 200×231 | 185×204 ^a | 69×90 | 61×80 | **34×56** | 55×74 | **34×56 (22)** | — |
| E+fλ 7pt (elim. λ) | 19 | — | 52×71 | 24×43 | **22×41** | **22×41** | **22×41** | — | **22×41** |
| Rolling-shutter pose | 8 | 48×56 | 47×55 | 47×55 | 47×55 | **40×52** | 47×55 | **40×52 (12)** | — |
| Refractive P5P | 16 | 280×399 | 240×324 | 157×246 | 68×93 | 57×73 ^b | 57×73 | **38×58 (20)** | — |
| Abs. pose quivers | 20 | 372×386 | 217×253 | 81×119 | 68×92 | **37×71** ^d | 65×85 | 56×80 (24) | — |
| Unsynch. rel. pose | 16 | 633×649 | 159×175 | 299×315 | 150×168 | 59×80 | 139×155 | **59×79 (20)** | — |
| Optimal PnP (Hesch) | 27 | 120×120 | 88×115 | 88×115 | **87×114** | **87×114** | **87×114** | — | **87×114** |
| Optimal PnP (Cayley) | 40 | 118×158 | 118×158 | 118×158 | 118×158 | 118×158 | 118×158 | **78×131 (53)** | — |
| Opt. pose 2pt v2 | 24 | 192×216 | 192×216 | 192×216 | 176×200 | 112×146 | 139×163* | **87×120 (33)** | — |
| Rolling-shutter R6P | 20 | 196×216 | 204×224 | — | — | — | 120×140 | **66×92 (26)** | — |
| Refractive P6P+focal | 36 | 648×917 | 636×851 ^e | — | — | — | 126×162 | **40×61 (21, symmetry)** | — |
| fλ+E+fλ 7pt | 68 | 886×1011 | 581×862 ^a | — | — | — | 209×277 | **138×210 (72)** | — |
| Optimal 3-view triangulation | 47 | — | — | — | — | — | 401×448 | **69×127 (58)** | — |

Footnotes:

- ^a The Syzygy numbers differ between papers. f+E+f is 31×46 in Greedy22 and CVPR 2020, 31×50 in Syzygy17 itself and 25×40 in JMIV. E+f is 19×28 in JMIV. The 9pt problem is 165×200, 165×189 or 189×213. λ+E+λ is 154×210, 149×205, 149×201 or 154×206. fλ+E+fλ is 581×862 or 581×659. These are re-runs or version changes of autogen, so cite the specific table.
- ^b The input system is simplified as in Martyushev 2022 (JMIV marker "*"). The 11×26 for f+E+f is therefore not directly comparable with symgen's 12×30 on the original 10-equation formulation.
- ^c For the 9pt problem, BeyondGB18 gives 84×117 heuristic and 87×111 GFan, while JMIV lists GFan 105×129 and heuristic 73×97. The notes flag this as a possible pdftotext column misalignment.
- ^d The JMIV Table 1 says 37×71 but its Fig. 4 legend says 61×86.
- ^e In Syzygy17 this was 2,196×1,913 → 636×851 with symmetry removal; Greedy22 lists 636×654.
- \* Reduced with a Schur complement.

Runtimes are sparse in these tables. The Laurent paper reports MATLAB times that exclude coefficient construction: 9pt 0.3 ms against 0.4 ms for Greedy, R6P 0.3 against 0.9 ms, λ1+E+λ2 7pt 4.2 against 27.9 ms, and fλ+E+fλ 2.3 against 3.0 ms. It also warns that "template size alone is not an appropriate measure": refractive P5P at 38×58 carries a 20×20 eigenproblem, where Greedy's 57×73 carries 16×16 ([Laurent23](https://arxiv.org/abs/2307.00320)).

The symgen runs used Octave on 200 identical instances (symgen experiments in this repo). Median log10 residual and time per solve went from the Maple generator to symgen as follows:

| Problem | Template | Median log10 residual | Time per solve |
|---|---|---|---|
| 6pt focal | 12+18 | −12.02 → −12.18 | 7.9 → 3.2 ms |
| 7pt E+fλ (elim. λ) | 22+19 | −11.21 → −12.44 | 5.9 → 1.0 ms |
| 9pt two radial | 90+27 | −7.90 → −8.09 | 15.7 → 1.4 ms |

The template sizes are identical. The gains come from leaner code and from better variable recovery, which reads each unknown from the component of the eigenvector with the largest magnitude.

Stability comparisons use incompatible protocols. JMIV reports the mean log10 residual and the failure rate over 5K instances; symgen reports the median over 200. Within JMIV, the resultant solver for λ+E+λ reaches −11.86 with 0.06% failures, against −7.50 and 5.16% for Syzygy, and unsynchronized pose reaches −13.59 with 0% failures, against 3.43% for Syzygy ([arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)).

The main message of the table is that **Laurent23 is now the reference to beat**. It holds the smallest template on most rows where it reports a result, often by trading extra roots for fewer rows. The sparse-resultant generator still holds unique bests on stitching (8×31 eigen 23, against Laurent's 6×30 eigen 24, which is effectively a tie) and on quivers.

## The resultant family spans a speed-stability spectrum that nobody has optimized

Two papers by the user bracket the resultant design space. The extra-polynomial method of CVPR 2020 has a small eigenproblem, but its accuracy is governed by the conditioning of the inverse Â12. The hidden-variable variant of ICPR 2020 removes the inverse entirely: A and B hold raw coefficients, and parasitic zero and infinite eigenvalues are stripped by row and column elimination. The cost is a larger generalized eigenproblem (GEP) ([Hidden20](https://arxiv.org/abs/2007.10100)):

| Problem | GEP size (raw → reduced) | Hidden-variable time | Extra-polynomial time | Syzygy time |
|---|---|---|---|---|
| P5Pr | 80 → 36 | 0.474 ms | 0.334 ms | 0.494 ms |
| P6Pfr | 290 → 110 | 5.19 ms | 1.53 ms | 4.62 ms |
| P2Pr | 160 → 124 | 6.56 ms | 0.729 ms | 10.57 ms |

The hidden-variable solvers were "significantly numerically more stable" near degenerate refractive configurations.

The CVPR 2020 supplement shows the same trade-off on the standard set. Each pair below gives the hidden-variable GEP size first, then the extra-polynomial inverse and eigen sizes ([SparseR20](https://arxiv.org/abs/1912.10268)):

| Problem | Hidden-variable GEP | Extra-polynomial (inverse / eigen) |
|---|---|---|
| f+E+f 6pt | 18² | 12 / 18 |
| E+fλ 7pt | 36² | 61 / 19 |
| 9pt two radial | 68² | 90 / 27 |
| λ+E+λ | 110² | 39 / 56 |

For 6pt focal the GEP is the same size as symgen's eigen basis. The BMVC 2008 quadratic eigenproblem solved 6pt focal from a 10×30 matrix with 20 eigenvalues for 15 solutions ([BMVC 2008](https://www.bmva-archive.org.uk/bmvc/2008/papers/164.pdf)).

The newest point on the spectrum is Wu, Bhayani and Heikkilä's CVPR 2026 solver ([arXiv 2605.06572](https://arxiv.org/abs/2605.06572)). It interpolates det M(x_h) by inverse FFT from sampled evaluations, recovers the other unknowns from rank-1-deficient submatrices by Cramer's rule, and applies a GCD criterion for noise. It claims strong stability and competitive runtime on small problems, though the notes did not extract its numbers. An intermediate point exists in robotics: Martyushev's 2025 Stewart-Gough solver applies a sparse 293×362 template, then PLU, then a 69×69 generalized eigenpair for 40 solutions ([arXiv 2505.00634](https://arxiv.org/abs/2505.00634)). That is a partial elimination ending in a GEP, which is exactly the "in-between" regime the thesis hypothesizes.

On the numerical side the toolbox is unchanged since Byröd et al. Their methods are QR with column pivoting or SVD basis selection over a redundant basis, plus extraction from eigenvalues. They gave up to 10⁷ better accuracy on 3-view triangulation, at 1.4× the cost (QR) or 5× (SVD) ([Byröd 2009](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)). The Laurent paper's triangulation comparison shows how much this still matters:

| Solver | Time | Failure rate |
|---|---|---|
| Laurent | 1.34 ms | 1.11% |
| Byröd standard | 1.54 ms | 3.09% |
| Byröd QR | 2.07 ms | 0.29% |
| Saturated (Sat17) | 1.56 ms | 8.91% |

([Laurent23](https://arxiv.org/abs/2307.00320)). Xu, Hu, Tsakiris and Kneip showed that column permutations of a single template give differently conditioned solvers, and trained a classifier to pick one online ([arXiv 2401.09328](https://arxiv.org/abs/2401.09328)). Truncated normal forms (TNF), from Telen, Mourrain and Van Barel, generalize Byröd's basis choice to cokernels of resultant maps with non-monomial bases ([arXiv 1803.07974](https://www.arxiv.org/abs/1803.07974)). Bender and Telen's toric eigenvalue methods handle zero-dimensional schemes in toric varieties and overdetermined systems, with Julia code in EigenvalueSolver.jl ([arXiv 2006.10654](https://arxiv.org/abs/2006.10654); [arXiv 2105.08472](https://arxiv.org/abs/2105.08472)). No computer-vision generator uses either.

## HC owns problems above about 100 roots, and Galois theory says when to stop optimizing

Homotopy continuation has matured into three regimes:

- **All-path CPU HC.** MINUS solves the 312-solution Chicago problem in 660 ms on average and fails about 1% of the time ([MINUS](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)).
- **All-path GPU HC.** GPU-HC (CVPR 2022) handles trifocal pose with unknown focal length (1784 solutions) in 338 ms. The current README reports 16.7 ms, and 0.88 ms for 3-view triangulation. Systems are capped at 32×32 and must be square ([GPU-HC paper](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf); [repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)).
- **Learned single-path HC.** Hruby et al. pick a start pair with an MLP and track one real path. They solve 26.3% of Scranton (272 solutions) instances in 16.3 µs, or 61.6 µs per correct solution, and 5pt in 26.1 µs per correct solution against about 5 µs for the algebraic solver ([Hruby 2022](https://arxiv.org/abs/2112.03424)). Simulator HC regresses a solution and simulates a matching start problem. It reports 100% success on generalized P4P in 0.316 ms against 1.254 ms for a Gröbner solver, and about 17× speedup on the 140-solution generalized relative pose with scale ([Simulator HC](https://arxiv.org/html/2411.03745)).

GPU-HC's own table puts the crossover clearly. Templates are 3–30× faster for P3P and 5pt. GPU-HC wins 38× on 3-view triangulation (94 solutions) and 5× on optimal PnP (128). Template construction runs out of memory above about 160 solutions. The Laurent generator's 205 ms, 1610×1726 template for 9-line relative pose (116 roots), which took about 13 h to reduce, marks the practical ceiling for templates ([Laurent23](https://arxiv.org/abs/2307.00320)).

Algebraic structure analysis decides how far reformulation can go. Duff, Korotynskiy, Pajdla and Regan compute Galois/monodromy groups numerically, then use group theory to detect decompositions ([arXiv 2105.04460](https://arxiv.org/abs/2105.04460)):

- 5pt is S2≀S10∩A20, which explains Nistér's degree-10 formulation.
- A pair of calibrated homographies drops from degree 64 to 16.
- Problems with the full symmetric group, such as S216 and S312, are already optimal.

Hruby et al.'s four-view radial geometry tracks 28 paths to recover 3584 solutions ([CVPR 2023](https://openaccess.thecvf.com/content/CVPR2023/html/Hruby_Four-View_Geometry_With_Unknown_Radial_Distortion_CVPR_2023_paper.html)). Duff's 2025 "Galois width" gives 10 for 5pt and 8 for P3L ([arXiv 2507.10407](https://arxiv.org/html/2507.10407v1)). Duff and Lee (2026) certify monodromy actions with certified path tracking ([arXiv 2603.17288](https://arxiv.org/abs/2603.17288)). Monodromy plus capture-recapture estimates give true solution counts when BKK bounds overshoot ([Hauenstein–Sherman](https://arxiv.org/abs/2005.00327)).

In template generators, symmetry is the lever that pays off most often. The Laurent generator's automatic partial p-fold detection cuts Weak PnP to 10×20 and refractive P6P+focal to 40×61 ([Laurent23](https://arxiv.org/abs/2307.00320)). autogen detects variable-aligned symmetries by exponent residues ([Syzygy17](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)). The ICPR 2020 P6Pfr solver explicitly did not remove its symmetry ([Hidden20](https://arxiv.org/abs/2007.10100)).

Code availability is good but tied to proprietary tools:

| Tool | Generator stack | Output |
|---|---|---|
| Martyushev EliminationTemplates | Maple | MATLAB solvers |
| Bhayani aut_gen_sparse_res_solver | MATLAB + Maple, MIT licence | MATLAB solvers |
| autogen / GAPS | Macaulay2 + MATLAB | MATLAB solvers |
| Macaulay2 EliminationTemplates package (v1, 2026) | Open-source Macaulay2, Greedy and syzygy approaches | — |
| MINUS, GPU-HC, learning_minimal | — | C++ / CUDA solvers |

Sources: [GAPS](https://arxiv.org/pdf/2004.11765); [M2 package](https://www.macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EliminationTemplates/html/index.html). None of the symbolic generators documents an optimized C++ backend. Greedy22 notes that its MATLAB fλ+E+fλ solver (8.5 ms) is about 3× slower than a C++-core solver (2.9 ms) ([Greedy22](https://arxiv.org/abs/2203.14901)).

## The thesis's six open problems: two advanced, two untouched, two sidestepped

| # | Open problem (user's thesis / JMIV 2024) | What the literature since 2022 has done | Status |
|---|---|---|---|
| 1 | When AM = SRes; normal forms not used in CV generators | JMIV Definition 2 and Props. C.1–C.3 give constructive two-way equivalence ([arXiv 2301.06443](https://arxiv.org/abs/2301.06443v2)). TNF and toric theory matured ([1803.07974](https://www.arxiv.org/abs/1803.07974); [2006.10654](https://arxiv.org/abs/2006.10654)). No CV generator uses them, and the 2026 M2 package implements only Greedy and syzygy. | Equivalence partly answered; normal forms still untouched |
| 2 | Macaulay-style determinant/minor formula for mixed sparse resultants | No CV paper addresses it, and the notes did not survey post-2022 pure-mathematics work, so its current mathematical status is unverified here. | Open (and not on the critical path for solvers) |
| 3 | Optimal resultant matrix for non-generic coefficients | JMIV iterates over polytope subsets. Greedy and Laurent bypass polytopes and build templates from the ideal itself, which takes non-genericity into account by construction. symgen's exact rank tests at the specialized parametric point also respect it. | Sidestepped in practice, no optimality result |
| 4 | Best point on the hidden-GEP ↔ extra-polynomial spectrum | ICPR 2020 and CVPR 2020 data quantify the two ends. FFT interpolation (2026) adds an inversion-free end. Stewart-Gough (2025) is an intermediate PLU→GEP example. Laurent23 states that size ≠ efficiency. No cost model or automated selection exists. | Advanced empirically, not optimized |
| 5 | Unified AM/SRes monomial-set search; border bases; convexity | The Laurent generator, co-authored by the user, searches Laurent multiplier sets independent of Gröbner bases and polytopes. BeyondGB and Greedy already left standard Gröbner bases behind. Border bases are not used in any CV generator. symgen's matroid step gives the minimal eigen size without iterative tests. | Substantially advanced; border bases and convexity untested |
| 6 | Data-driven hybrid RANSAC choosing a solver per iteration | Learned online template-permutation choice ([2401.09328](https://arxiv.org/abs/2401.09328)). Hruby's TRASH class skips hopeless samples. An unverified SSRN GNN solver selector. No peer-reviewed per-iteration solver-switching RANSAC. | Open |

Problem 1 is the most consequential. JMIV proves that once the template and the extra polynomial are fixed, the two families differ only in search heuristics, not in online numerics. TNF adds the missing ingredient: a numerically chosen basis, as in Byröd's QR pivoting, rather than one fixed offline. A generator that emits a TNF-style solver, with a fixed template offline and a pivoted basis choice online, would answer the JMIV remark directly.

Problem 4 now has enough data to be posed as an optimization, but nobody has done so. The symgen sweep shows that within the extra-polynomial family the minimal eigen size is always fastest, because each extra eigen column removes only one column of A2. Yet on 6pt focal the saturated 25+15 solver and the 12+18 solver take nearly the same pure linear-algebra time, 77 µs and 70 µs (symgen experiments in this repo). A naive flop model with eig ≈ 25× LU per n³ would predict about a 1.5× win for 25+15. So at n ≤ 30, constant factors and coefficient evaluation dominate, and any cost model must be calibrated on compiled code, not counted in flops.

Problem 5 has largely been absorbed into the Laurent generator, which is exactly the "unified" search, since its output templates serve both interpretations. What remains distinct for the user's line is the guarantee side. Neither the Greedy nor the Laurent paper certifies minimality. symgen's matroid argument does certify the minimum eigen size for a given row set, which is a small but real theoretical foothold.

## symgen's fastest wins are search breadth, C++ output and symmetry, not new theory

symgen's strengths are real:

- **Determinism.**
- **Exact certificates.** Full rank mod p = 2097143 proves that a minor is a nonzero polynomial in the data, so every emitted template is generically valid.
- **A provable per-row-set eigen minimum.**
- **No Maple.**
- **Better residuals and 2–11× faster Octave solvers at equal size.**

Its weaknesses are just as clear:

- On the problems tested so far it reproduces CVPR 2020 sizes, not JMIV, Greedy or Laurent sizes: 12×30 against 11×26 for 6pt focal, and 90×117 against 73×97 for 9pt.
- Its only reductions are greedy column and row moves in fixed orders.
- Its saturation step is hand-derived.
- It emits MATLAB/Octave only.
- A rank deficiency found mod p is only probabilistic evidence of deficiency over ℚ (Schwartz–Zippel: probability at most about deg/p of a false rejection). Its "cannot go lower" claims are therefore high-confidence, not proofs, unless repeated at several seeds or primes.

The directions below are ordered by expected payoff per unit of risk.

| Priority | Direction | Expected payoff | Risk | How to evaluate |
|---|---|---|---|---|
| 1 | **Multiplier/greedy parameter search on exact GF(p) tests**: Greedy22-style scalar-θ parameterization and multiplier-set growth and shrink, plus randomized-restart greedy over removal orders, all decided by the existing rank Top[:,G]=\|G\| test; also import the "*" simplified input systems and G–J pre-reduction | Close the 10–25% row gap (9pt 90→≤77 rows; 6pt 12×30→11×26); possibly match Laurent on several rows | Low–medium: the test oracle exists; search cost grows with template size | Reproduce JMIV Table 1 and Laurent Table 1 rows on identical formulations; report rows, eigen size, offline time, and run-to-run variance (which no paper reports) |
| 2 | **C++ code generation with structured coefficient evaluation** (intermediate symbols, CSE, fixed-size Eigen/LAPACK, optional sparse LU) | 3×-class speedups suggested by Greedy22's MATLAB vs C++ gap; removes the Octave bottleneck, which is coefficient evaluation | Low | µs/solve and per-correct-solution time against Martyushev MATLAB and the C++ in PoseLib-style libraries on the same data |
| 3 | **Automatic symmetry detection**: exponent-residue (c-weighted p-fold) detection as in autogen, then a quotient eigen basis | 2–10× when symmetry is present (Weak PnP 108×124→10×20; P6Pf 126×162→40×61 in Laurent23) | Low for variable-aligned, medium for partial or unaligned | Run on Weak PnP, refractive P5P/P6Pf, generalized pose; compare with Laurent23 |
| 4 | **Automated saturation via syzygies of leading forms**: generalize the 6pt-focal derivation by computing, over GF(p), syzygies among the leading coefficients in the hidden variable, emitting lower-degree combinations with data-polynomial coefficients, and keeping them only if the eigen count drops to d | Removes spurious roots (18→15) and enables every hidden variable; opposite trade to Laurent23, which accepts extra roots | Medium: on 6pt focal it cost a 25×25 inverse and about 1 digit of accuracy | Compare d against the eigen size across the benchmark; measure residual and time with structured coefficients; verify d independently via monodromy counts or the quotient-ring dimension |
| 5 | **Speed–stability cost model including the hidden-variable GEP end and intermediate PLU→GEP points**: implement GEP templates (Heikkilä/ICPR) and partial Schur in symgen, time them in C++, fit a cost model in (inverse size, eigen size, GEP size, nnz), and select automatically | Directly answers thesis problem 4; principled choice per problem; stability near degeneracies | Medium: needs compiled timings; small-n constants are noisy | ICPR protocol (random plus near-degenerate scenes), fail% at 1e−3, µs; include the FFT solver of CVPR 2026 as a baseline |
| 6 | **Laurent/torus multipliers**: allow negative exponent shifts in T_j and B, i.e. a δ-shift into the torus | Path toward Laurent23 sizes inside the resultant framework | Medium–high: overlaps heavily with Laurent23 itself | Sizes on stitching, P4Pfr, R6P, and opt. pose 2pt against Laurent23 |
| 7 | **TNF / QR-pivoted online basis**: keep a slightly redundant eigen basis and choose it online by column pivoting (Byröd/TNF) | Stability on hard rows (λ+E+λ, 9pt at −8 median) at about 1.4× cost | Medium: no longer pure GF(p) offline; needs float validation | Fail% and median residual on 9pt and λ+E+λ; time overhead |
| 8 | **Galois/monodromy diagnostics** (HomotopyContinuation.jl or M2 MonodromySolver offline) | Tells when reformulation, not template search, is the bottleneck | Low cost, low direct payoff | Flag imprimitive or deck-symmetric benchmark problems |
| 9 | **Hybrid RANSAC solver selection** (thesis problem 6) | Robustness in near-degenerate data | High: data-dependent, no baseline | RANSAC success and time on real datasets with the 12+18 vs 25+15 vs GEP variants |

Directions 1 and 2 are the obvious first quarter of work. Direction 1 attacks the only axis on which symgen is measurably behind, and its exact oracle makes greedy search cheap and reproducible, where the published tables suffer from run-to-run drift. Direction 2 attacks the measured bottleneck. Direction 4 is the most novel: no generator yet automates saturation that avoids an online elimination, and Larsson's saturation lifts by f^N, which grows templates. But its 6pt evidence is mixed, with a larger inverse, a lost digit and equal speed, so it should be framed as a tool for exact eigen size and full hidden-variable choice, not for speed. Direction 5 is where the user's thesis has a unique claim, and the CVPR 2026 FFT work already stakes out the inversion-free end.

## Conclusion

Template-size minimization on the classic benchmark is near saturation: on the small classic problems the best methods now differ by only a few rows, and the frontier has moved to three places: problems above about 100 roots (HC), honest efficiency metrics (eigen size, constants, stability), and open tooling. That reframes symgen's value. It will not win by out-shrinking Laurent23 on 6pt focal. It can become the first open-source, Maple-free, deterministic generator whose decisions are certified mod p, and the first to choose automatically between the inverse, Schur and GEP ends using a calibrated cost model. Neither Greedy, Laurent nor the 2026 Macaulay2 package offers either.

The equivalence theorem in JMIV 2024 has a practical corollary that the literature has not used. Since AM and SRes solvers coincide for a fixed template, the remaining scientific questions concern the online map from template to roots: basis choice (TNF), spectrum position (GEP versus inverse), and spurious-root handling (saturation versus extra roots). These are exactly thesis problems 1, 4 and 5, and they are better posed now than in 2024.
