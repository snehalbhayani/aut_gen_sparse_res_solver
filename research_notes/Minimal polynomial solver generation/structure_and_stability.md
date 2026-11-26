# Algebraic structure analysis of minimal problems and numerical stability/efficiency of minimal solvers

## Q1. What do Galois/monodromy group methods say about decomposable problems and symmetries, and how does decomposition give faster solvers?

### Takeaway
The Galois/monodromy group of the branched cover "problem-solution incidence variety -> data space" shows whether a minimal problem decomposes. An imprimitive group means the problem splits into a chain of smaller subproblems, and a nontrivial deck group means the problem has symmetries. A full symmetric group means the formulation is already algebraically optimal. The classic decompositions are 5pt relative pose (20 = 2 x 10, Nistér's degree-10 formulation), P3P (8 = 2 x 4, sign symmetry) and calibrated homography (12). Since 2021, Galois analysis has produced concrete degree reductions: 64 -> 16 for a pair of calibrated homographies, and 3584 -> 28 + 2 + 4 subproblems for four-view radial-camera geometry. It has also made the mixed point-line absolute pose solvers faster.

### Cited Findings
**Duff, Korotynskiy, Pajdla, Regan, "Galois/monodromy groups for decomposing minimal problems in 3D reconstruction"**
- Published as arXiv 2105.04460 (May 2021) and in SIAM J. Appl. Algebra Geometry (2022). — [arXiv](https://arxiv.org/abs/2105.04460); [dblp](https://dblp.uni-trier.de/rec/journals/siaga/DuffKPR22.html)
- **Method.** They compute the group numerically by monodromy (homotopy continuation), then reduce the questions of decomposability and deck transformations to standard computational group theory. — [arXiv](https://arxiv.org/abs/2105.04460)
- **Software.** Bertini, HomotopyContinuation.jl and Macaulay2 MonodromySolver for the monodromy; GAP for group computations; the Macaulay2 FGLM package for symbolic certification. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **5pt relative pose (degree 20).**
  - The group is S2 ≀ S10 ∩ A20, which is imprimitive.
  - The intermediate essential-matrix cover has monodromy S10 and degree 10, which matches the formulation in Nistér's solver.
  - The "twisted pair" is a nonlinear deck transformation.
  - With the normalization ||t|| = 1 the problem becomes a degree-40 cover with a correspondingly larger group. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **P3P (degree 8).**
  - The group is S2 ≀ S4 ∩ A8 and the sign symmetry is a deck transformation.
  - The restriction to A8 was certified symbolically: the discriminant of the degree-8 eliminant is a square. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **Mixed point/line absolute pose.**
  - (2 points, 1 line): group C2×C2 (degree 4), deck group C2×C2. The FGLM constraints suggest that rotation and translation can be solved separately.
  - (1 point, 2 lines): group S2 ≀ S4 ∩ A8 (degree 8), with a new deck involution (R, t) -> (−R, −t − 2e3).
  - 3 lines: group S8, which is primitive. Any solver must therefore handle a degree-8 polynomial. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **Calibrated two-view homography (degree 12).** The group is C2×C2 ⋊ (S2 ≀ S3 ∩ A6) and the deck group is C2×C2 (the twisted pair plus a reflection-rotation symmetry). — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **Pair of calibrated homographies across three views.** The degree drops from 64 to 16, and the analysis yields new constraints on compatible homographies. — [arXiv](https://arxiv.org/abs/2105.04460)
- **Survey of a table of minimal problems with degree < 1000.**
  - Every problem is either imprimitive or has the full symmetric group.
  - The imprimitive groups act on 12, 16, 20, 32 and 64 points, and one of them (on 32 points) is solvable.
  - Two problems have full symmetric groups S216 and S312. This shows that the earlier homotopy solvers for these problems already track the minimum number of paths. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **No timings.** The paper reports structure and degrees only. — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)

**Duff, "Numerically Computing Galois Groups of Minimal Problems" (ISSAC 2025 tutorial, arXiv 2507.10407)**
- It introduces the "Galois width" as a measure of intrinsic difficulty. Values:

  | Problem | Degree | Galois width |
  |---|---|---|
  | P3P | 8 | 3 |
  | P3P1L | 4 | 2 |
  | P1P2L | 8 | 3 |
  | P3L | 8 | 8 |
  | 5pt relative pose | 20 | 10 |

  — [arXiv html](https://arxiv.org/html/2507.10407v1)
- The width of 10 for 5pt is said to show that the Nistér–Stewénius method is algebraically optimal. Nistér et al. (2007) pioneered computing Galois groups of minimal problems, using arithmetic Galois groups. — [arXiv html](https://arxiv.org/html/2507.10407v1)
- For 5pt, tracking 10 paths suffices to recover all 20 solutions. — [arXiv html](https://arxiv.org/html/2507.10407v1)
- Brysiewicz et al. (2021) are cited for the result that decomposability is equivalent to imprimitivity. Améndola et al. (2016) give the "decomposable monodromy" technique, implemented in MonodromySolver as the Equivalencer option. — [arXiv html](https://arxiv.org/html/2507.10407v1)
- Hruby et al. (2024) "significantly improved" the mixed (2,1) and (1,2) point-line solvers by exploiting the decomposition that the Galois computation suggests. Galois width computations have led "to more efficient minimal solvers than previously known" (Hruby et al. 2023a, 2024). — [arXiv html](https://arxiv.org/html/2507.10407v1)
- **Caveat.** Numerical monodromy only shows that the computed group is a subgroup of the true Galois group, because path tracking is uncertified. — [arXiv html](https://arxiv.org/html/2507.10407v1); [Hauenstein–Rodriguez–Sottile, FoCM 2018](https://arxiv.org/abs/1605.07806)

**Hruby, Korotynskiy, Duff, Oeding, Pollefeys, Pajdla, Larsson, "Four-View Geometry with Unknown Radial Distortion" (CVPR 2023)**
- The calibrated system has 3584 solutions.
- Each Euclidean symmetry class has 16 members, so 224 class representatives are found by solving subproblems with 28, 2 and 4 solutions.
- Computational Galois theory showed the decomposition.
- The minimal case is 13 points in 4 views, or 7 points with 4 upright cameras. — [MSR page](https://www.microsoft.com/en-us/research/?p=947241); [CVF](https://openaccess.thecvf.com/content/CVPR2023/html/Hruby_Four-View_Geometry_With_Unknown_Radial_Distortion_CVPR_2023_paper.html)
- The Galois groups were computed heuristically and numerically, then checked with Macaulay2 and HomotopyContinuation.jl. — [CVPR supplement](https://openaccess.thecvf.com/content/CVPR2023/supplemental/Hruby_Four-View_Geometry_With_CVPR_2023_supplemental.pdf)
- The solver tracks 28 paths to recover all 3584 solutions. — [Duff 2025](https://arxiv.org/html/2507.10407v1)

**Certified Galois groups (2026).** Duff and Lee, "Certifying Galois/monodromy Actions via Homotopy Graphs" (arXiv 2603.17288, March 2026), use certified homotopy path tracking so that the computed monodromy action is guaranteed correct. They used it to certify several notable Galois groups. — [arXiv](https://arxiv.org/abs/2603.17288)

**PLMP and Kileel (classifications of minimal problems)**
- **PLMP (Duff, Kohn, Leykin, Pajdla, ICCV 2019; extended in TPAMI 2024).** Under complete visibility with calibrated cameras there are exactly 30 minimal point-line problems, counted up to adding lines in the two-view case.
  - By number of cameras: none for 7 or more views; 1, 3, 6 and 17 for 6, 5, 4 and 3 views.
  - Each problem has at most 5 points and 6 lines, except for the extra lines allowed in two views. — [arXiv 1903.10008](https://arxiv.org/pdf/1903.10008); [NSF PAR](https://par.nsf.gov/servlets/purl/10530443)
  - A follow-up (IJCV) covers partial visibility with three calibrated cameras. A 2025 work (Kiehn, Ahlbäck, arXiv 2503.04351) extends the setting to projective/uncalibrated SfM and finds 291 minimal problems. — [arXiv 2503.04351](https://arxiv.org/pdf/2503.04351)
- **Kileel, "Minimal problems for the calibrated trifocal variety" (SIAGA 2017).** Gives algebraic degrees for 66 minimal problems over three calibrated views. For example, the 3 PPP + 1 PPL problem has degree 272. — [arXiv 1611.05947](https://arxiv.org/abs/1611.05947); [solved instances page](https://web.ma.utexas.edu/users/jkileel/CalibratedMinimalProblems.html)

**Symmetry removal in template solvers.** Larsson and Åström, "Uncovering symmetries in polynomial systems" (ECCV 2016), detect and remove symmetries of the form p(−x) and reduce the solver size accordingly. Example: P6Pfr has only 18 solutions once symmetries are removed. — cited in [Bhayani et al. arXiv 2007.10100](https://arxiv.org/pdf/2007.10100)

### Inferences
- For a solver generator, the practical consequences are:
  - If the problem has a deck transformation (for example x -> −x), quotienting by it halves the action-matrix size.
  - If the group is imprimitive, it can pay to solve the base cover first and then the fiber subproblems. This is the 5pt 10 x 2 pattern.
  - If the group is the full symmetric group Sn, no reformulation can lower the degree, so engineering effort should go to template size and stability instead.
- The certified methods from 2026 make it possible in principle to prove a decomposition before investing effort in hand-crafting a solver for it.

### Gaps
- I did not verify the "PL1P" paper (point-line problems with one line incidence, about 2020) and its counts in this pass.
- I did not fetch the exact speedups from Hruby et al. 2024 for the mixed point-line absolute pose solvers.
- I did not confirm a "Brazelton & Raman" Galois-groups paper. It may not exist under that title.

## Q2. How are solution counts and degrees determined (BKK/mixed volume, monodromy), and how are spurious solutions handled?

### Takeaway
There are three routes to the solution count:
- **Bounds.** Bézout and BKK (mixed volume) bounds are often far too high for vision problems.
- **Monodromy.** It finds the actual generic count by populating a fiber. Because the expected number of paths grows linearly with the number of solutions, it scales well. Capture-recapture statistics (Lincoln-Petersen/Chapman/Schnabel) combined with a trace test estimate the count and decide when to stop.
- **Template generators.** They count solutions as the dimension of the quotient ring, computed exactly over finite fields. Spurious components are removed by saturation.

### Cited Findings
- **Duff, Hill, Jensen, Lee, Leykin, Sommars, "Solving polynomial systems via homotopy continuation and monodromy"** (arXiv 1609.08722; IMA J. Numer. Anal. 2019).
  - It describes monodromy solvers through "decorated graphs" of homotopies.
  - Assuming the monodromy actions are generated uniformly, the expected number of paths tracked grows linearly with the number of solutions.
  - The software is competitive with the state of the art. — [arXiv](https://arxiv.org/abs/1609.08722)
- **Hauenstein and Sherman, "Using monodromy to statistically estimate the number of solutions"** (arXiv 2005.00327; published 2022 in a Springer volume, DOI 10.1007/978-3-030-91352-6_4). Note the authors: this is not a Duff paper.
  - Each monodromy loop is treated as a capture-recapture experiment.
  - The Lincoln-Petersen estimate is #S·#E/#(S∩E), where S and E are the start and end points of a loop. Chapman's bias-corrected version is (#S+1)(#E+1)/(#(S∩E)+1) − 1, and the Schnabel estimator pools several loops.
  - A trace test serves as the stopping criterion.
  - The method requires a transitive monodromy group.
  - The motivation is that classical upper bounds are often orders of magnitude too large. — [arXiv](https://arxiv.org/abs/2005.00327); [preprint](https://www3.nd.edu/~jhauenst/preprints/hsMonodromyStats.pdf)
- **Degree of the Galois group vs. paths tracked.** The degree of the cover equals the number of complex solutions for generic data. When the group is the full symmetric group, homotopy solvers that track exactly that many paths are optimal (S216 and S312 examples). — [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.04460)
- **Mixed-volume and resultant constructions.**
  - Sparse-resultant generators build the template from Newton polytopes, using the Canny-Emiris subdivision or Heikkilä's algorithm. These constructions can fail: Canny-Emiris did not terminate on P6Pfr and P2Pr, and Heikkilä's matrices were rank deficient (rank 32 of 35 on P5Pr; rank 76 of 82 on P6Pfr).
  - The resultant methods also produce parasitic zero or infinite eigenvalues, which must be removed. Bhayani et al. reduce the generalized eigenvalue problem (GEP) from 80×80 to 36×36 on P5Pr, from 290 to 110 on P6Pfr, and from 160 to 124 on P2Pr. — [Bhayani et al. 2007.10100](https://arxiv.org/pdf/2007.10100)
- **Saturation.** Larsson et al., "Polynomial solvers for saturated ideals" (ICCV 2017), is the standard reference for removing spurious components in template generators. — cited in [Bhayani et al.](https://arxiv.org/pdf/2007.10100)

### Inferences
- A robust generator pipeline would:
  1. Compute the mixed volume as an upper bound.
  2. Use monodromy or an exact finite-field quotient-ring computation for the true count.
  3. Compare the two. A gap signals solutions at infinity or degenerate components that call for saturation or elimination of zero/infinite eigenvalues.

### Gaps
- I did not collect quantitative examples comparing BKK, Bézout and the true degree for specific vision problems (for example, 6pt focal length relative pose: true degree 15).
- I did not cover solutions at infinity beyond the parasitic eigenvalue removal.

## Q3. Which techniques improve the numerical stability of elimination-template solvers?

### Takeaway
The main stability tools are:
- Adaptive basis selection via QR with column pivoting or SVD, over an enlarged "permissible" monomial set (Byröd, Josephson, Åström 2007–2009).
- Extracting solutions from eigenvalues rather than eigenvectors.
- Truncated normal forms with QR/SVD-chosen bases (Telen, Mourrain, Van Barel).
- Hidden-variable resultant GEPs that avoid large matrix inversions (Bhayani et al.).

The reported gains range from about 10x to 10^6x in accuracy, at a cost of up to about 10x in speed.

### Cited Findings
**Byröd, Josephson, Åström**
- **ICCV 2007 ("Improving numerical accuracy of Gröbner basis polynomial equation solvers").**
  - Basis selection: any linearly independent combination of monomials in the permissible set M′ = {x^α ∈ M : x_k·x^α ∈ M} can serve as the basis of R[x]/I.
  - A change of basis C̃nb = Cnb·Q (Q orthogonal, chosen by SVD, not mixing permissible and non-permissible columns) is used to minimize the condition number of the non-basis block.
  - Solutions are extracted from eigenvalues of all the action matrices m_xk, using the shared eigenvectors, rather than from eigenvectors.
  - Reported test problems:

    | Problem | Unknowns | Solutions | Template | Accuracy gain |
    |---|---|---|---|---|
    | 6pt relative pose with unknown focal length | 3 | 15 | 34 eqs × 50 monomials | about 10x |
    | Generalized-camera relative pose | 3 | 64 | 101 eqs × 165 monomials | not given |
    | Optimal 3-view L2 triangulation | 3 | 50 | 225 × 209 | about 10^6x, making the problem solvable in double precision |

    — [Byröd et al. poster, Lund](https://lup.lub.lu.se/search/files/8178280/1244587.pdf)
- **ECCV 2008 ("A column-pivoting based strategy for monomial ordering...").**
  - QR with column pivoting gives a fast basis selection, plus an adaptive truncation of the Gröbner basis.
  - It is about 4x faster than the SVD method with "nearly as good" precision. — [mlanthology](https://mlanthology.org/eccv/2008/byrod2008eccv-column); [LU portal](https://portal.research.lu.se/en/publications/a-column-pivoting-based-strategy-for-monomial-ordering-in-numeric/)
- **IJCV 84(3):237–256, 2009 ("Fast and stable polynomial equation solving and its application to computer vision").** Combines QR column pivoting and SVD. Reports "dramatically improved numerical stability," allowing a larger class of problems to be solved. — [LU record](https://lup.lub.lu.se/record/1389303)
- **Adoption.** Haner and Åström (CVPR 2015) enlarged the action matrix and used QR basis selection for refractive absolute pose:
  - P5Pr: 16 solutions, but a 44×44 action matrix and a 280×399 template.
  - P6Pfr: 41×41 action matrix and a 648×917 template.

  — [Bhayani et al.](https://arxiv.org/pdf/2007.10100)

**Telen, Mourrain, Van Barel, truncated normal forms (TNF)**
- TNFs build the quotient-ring structure from the cokernel of a resultant map. The basis of C[x]/I is chosen for numerical stability rather than fixed by a monomial order.
- The work started with Telen and Van Barel, "A stabilized normal form algorithm for generic systems" (arXiv 1708.07670), and continued with "Solving polynomial systems via a stabilized representation of quotient algebras" (arXiv 1711.04543), later published as "Truncated normal forms for solving polynomial systems" (2018).
- The generalization (arXiv 1803.07974, Mourrain–Telen–Van Barel; later titled "Solving Polynomial Systems Efficiently and Accurately") allows non-monomial bases adapted to the expected root locations and non-generic zero-dimensional systems, with accurate results in double precision. — [arXiv 1708.07670](https://arxiv.org/abs/1708.07670v1); [arXiv 1711.04543](https://arxiv.org/abs/1711.04543v1); [arXiv 1803.07974](https://www.arxiv.org/abs/1803.07974)

**Bhayani, Kukelova, Heikkilä, "Computing stable resultant-based minimal solvers by hiding a variable" (arXiv 2007.10100, ICPR 2020)**
- The hidden-variable sparse resultant gives a GEP directly. It needs no matrix inverse or Gauss-Jordan step, at the cost of a larger eigenvalue problem.
- It removes zero eigenvalues through structured row operations.
- Its stability is comparable on random scenes and "significantly" better near degenerate configurations, where scene points are nearly coplanar with the camera center and the interface normal.
- It evaluates 1K noise-free scenes and 1K near-degenerate scenes, using the solution closest to ground truth. — [arXiv](https://arxiv.org/pdf/2007.10100)
- Sizes and mean times (MATLAB, 3.9 GHz i7) of the main steps:

  | Solver | P5Pr | P6Pfr | P2Pr |
  |---|---|---|---|
  | GB [Larsson 2017] | 199×215 template, 16×16 eig, 0.494 ms | 636×654, 18×18, 4.62 ms | 913×937, 24×24, 10.57 ms |
  | Heuristic GB [Larsson 2018] | 199×215, 16×16, 0.668 ms | 398×416, 18×18, 2.08 ms | 597×621, 24×24, 4.96 ms |
  | Extra-variable resultant [Bhayani CVPR 2020] | 78×93, 25×25, 0.334 ms | 248×300, 52×52, 1.53 ms | 142×174, 32×32, 0.729 ms |
  | Hidden-variable resultant (this paper) | GEP 36×36, 0.474 ms | GEP 110×110, 5.19 ms | GEP 124×124, 6.56 ms |

  The hidden-variable solvers are slower than the extra-variable ones but "significantly numerically more stable." — [arXiv](https://arxiv.org/pdf/2007.10100)
- **Journal version** (J. Math. Imaging Vision 2024, arXiv 2301.06443). Connects the extra-polynomial sparse resultant to the action-matrix method through a Schur complement. For some problems it gives smaller and more stable solvers than Gröbner-basis ones. — [Springer](https://link.springer.com/content/pdf/10.1007/s10851-024-01182-1.pdf); [arXiv](https://arxiv.org/abs/2301.06443v2)

### Inferences
- Every stability method trades template size or eigen-problem size for conditioning. Examples: the 16 -> 44 action matrix in P5Pr, and GEP 36 vs. eigenproblem 16.
- A generator should expose a "stable mode" option that turns on:
  - an enlarged basis with QR pivoting, or
  - a hidden-variable GEP,
  - plus Newton polishing.

### Gaps
- I did not extract median log10 error values from these papers. The Byröd results are histograms, and only the improvement factors (10x, 10^6x) are stated in text.
- I did not cover these techniques with sources in this pass:
  - choice of action variable or random linear form,
  - coefficient scaling and normalization,
  - eigenvalue condition numbers,
  - Newton refinement statistics.

## Q4. What speed engineering has been reported: template size, Sturm sequences, characteristic polynomials, eigen-solver alternatives?

### Takeaway
In many Gröbner-basis solvers the eigendecomposition dominates the runtime. Bujnak, Kukelova and Pajdla (CVPR 2012) replace it with one of two routes, each followed by Sturm-sequence root isolation restricted to a feasible real interval:
- a modified FGLM step that gives a univariate polynomial, or
- fast computation of the characteristic polynomial.

They report significant speedups but give no concrete number in the sources I checked. Template size remains the main lever. Heuristic basis selection (Larsson 2018) and sparse resultants cut template sizes and runtimes by up to about 14x (P2Pr: 10.57 ms -> 0.73 ms).

### Cited Findings
- **Bujnak, Kukelova, Pajdla, "Making minimal solvers fast" (CVPR 2012, DOI 10.1109/CVPR.2012.6247840).**
  - It notes that almost all existing Gröbner-basis solvers spend most of their time in the eigenvalue computation.
  - **Method 1.** A modified matrix FGLM in which polynomial division is replaced by matrix-vector multiplication. This yields a univariate polynomial whose roots are found by Sturm sequences, only on a feasible interval and only to the desired precision.
  - **Method 2.** For larger problems, the characteristic polynomial of the action matrix is computed quickly and then solved with Sturm sequences.
  - It claims to "significantly speed up many existing solvers." — [CMP colloquium abstract](https://cmp.felk.cvut.cz/cmp/events/colloquium-2012.04.12/a_kukelova.html); [mlanthology](https://mlanthology.org/cvpr/2012/bujnak2012cvpr-making)
- **Later use of Sturm sequences.** Kukelova et al. (ICCV 2015) find the roots of a degree-10 univariate polynomial with Sturm sequences in a radial-distortion solver. — [CVF](https://openaccess.thecvf.com/content_iccv_2015/papers/Kukelova_Efficient_Solution_to_ICCV_2015_paper.pdf)
- **Template size vs. runtime.**
  - Larsson et al. (CVPR 2018) heuristic basis selection: P6Pfr goes from a 636×654 to a 398×416 template (4.62 -> 2.08 ms), and P2Pr from 913×937 to 597×621 (10.57 -> 4.96 ms).
  - The extra-variable sparse resultant goes further: P2Pr 0.729 ms, P6Pfr 1.53 ms, P5Pr 0.334 ms. — [Bhayani et al. Table I](https://arxiv.org/pdf/2007.10100)
- **Speed of QR column pivoting.** About 4x faster than SVD-based basis selection. — [ECCV 2008](https://mlanthology.org/eccv/2008/byrod2008eccv-column)
- **Eigenvalue-only extraction.** Computing one eigendecomposition and reusing its eigenvectors to get the eigenvalues of all action matrices (by element-wise division) is both fast and more accurate than reading solutions from eigenvectors. — [Byröd et al. poster](https://lup.lub.lu.se/search/files/8178280/1244587.pdf)
- **Schur complement.** The extra-polynomial resultant solver works on the Schur complement of a submatrix of the resultant matrix. This reduces the problem to a small eigenproblem, but it requires a matrix inverse, which can be ill-conditioned. — [Bhayani et al.](https://arxiv.org/pdf/2007.10100); [JMIV 2024](https://link.springer.com/content/pdf/10.1007/s10851-024-01182-1.pdf)

### Inferences
- In a generated C++/MATLAB solver the costs split into three parts: template elimination, which scales roughly cubically in the number of rows, the eigen step, which is O(n^3) in the number of solutions, and the extraction step.
- For small n (at most about 20), a characteristic polynomial with Sturm sequences or a companion matrix restricted to real roots beats a full eigendecomposition.
- For large templates, sparse LU or structured elimination dominates the gains.

### Gaps
- I found no exact speedup figures for Bujnak 2012 in the accessible sources. The full paper is needed.
- I did not collect sourced numbers for sparse LU and real-eigenvalue-only solvers in this pass.
