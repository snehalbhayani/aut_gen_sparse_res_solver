# Gröbner-basis / action-matrix automatic generators of minimal solvers

Source abbreviations and URLs (all primary; PDFs downloaded and read in full text unless noted):
- [KBP08-thesis] Z. Kukelova, PhD thesis "Algebraic Methods in Computer Vision", CTU Prague 2013 (contains the ECCV 2008 generator as Ch. 6 and its solvers in Ch. 7) — https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf
- [KBP08] Kukelova, Bujnak, Pajdla, "Automatic Generator of Minimal Problem Solvers", ECCV 2008, LNCS 5304, pp. 302–315, DOI 10.1007/978-3-540-88690-7_23 — https://www.doi.org/10.1007/978-3-540-88690-7_23 (I could not get the Springer PDF; the thesis was used instead)
- [Syzygy17] Larsson, Åström, Oskarsson, "Efficient Solvers for Minimal Problems by Syzygy-based Reduction", CVPR 2017 — https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf ; poster: https://openaccess.thecvf.com/content_cvpr_2017/poster/270-poster.pdf
- [Sat17] Larsson, Åström, Oskarsson, "Polynomial Solvers for Saturated Ideals", ICCV 2017 — https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf
- [BeyondGB18] Larsson, Oskarsson, Åström, Wallis, Kukelova, Pajdla, "Beyond Gröbner Bases: Basis Selection for Minimal Solvers", CVPR 2018 — https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf (arXiv 1803.04360)
- [LarssonThesis] V. Larsson, PhD thesis "Computational Methods for Computer Vision: Minimal Solvers and Convex Relaxations", Lund 2018 — https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf
- [Byrod09] Byröd, Josephson, Åström, "Fast and Stable Polynomial Equation Solving and its Application to Computer Vision", IJCV 84(3), 2009, DOI 10.1007/s11263-009-0235-z — author version https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf
- [Greedy22] Martyushev, Vrablikova, Pajdla, "Optimizing Elimination Templates by Greedy Parameter Search", CVPR 2022 — https://arxiv.org/abs/2203.14901 ; CVF: https://openaccess.thecvf.com/content/CVPR2022/papers/Martyushev_Optimizing_Elimination_Templates_by_Greedy_Parameter_Search_CVPR_2022_paper.pdf
- [Laurent23] Martyushev, Bhayani, Pajdla, "Automatic Solver Generator for Systems of Laurent Polynomial Equations", arXiv 2307.00320 (2023) — https://arxiv.org/pdf/2307.00320
- [SparseR20] Bhayani, Kukelova, Heikkilä, "A sparse resultant based method for efficient minimal solvers", CVPR 2020 — https://openaccess.thecvf.com/content_CVPR_2020/papers/Bhayani_A_Sparse_Resultant_Based_Method_for_Efficient_Minimal_Solvers_CVPR_2020_paper.pdf
- [Hidden20] Bhayani, Kukelova, Heikkilä, "Computing stable resultant-based minimal solvers by hiding a variable", ICPR 2020 — https://arxiv.org/abs/2007.10100
- [SRConn23] Bhayani, Heikkilä, Kukelova, "Sparse resultant based minimal solvers in computer vision and their connection with the action matrix", arXiv 2301.06443 (v2 Sep 2023) — https://arxiv.org/pdf/2301.06443
- [GAPS20] Bo Li, Viktor Larsson, "GAPS: Generator for Automatic Polynomial Solvers", arXiv 2004.11765 — https://arxiv.org/pdf/2004.11765

Notation used in tables: template size = rows × columns of the elimination template (Macaulay/coefficient matrix that is G–J eliminated online); "d"/"#sol" = dimension of quotient ring = action-matrix size (unless symmetry reduces it). λ = radial distortion (division model), f = focal length; "f+E+f 6pt" = 6pt relative pose with shared unknown focal; "E+f 6pt" = one-sided focal; "λ1+F+λ2 9pt" = two different radial distortions; "F+λ 8pt" = one-sided radial distortion 8pt; "λ+F+λ 8pt" = shared radial distortion 8pt; "P4P+fr" = P4P with unknown focal + radial distortion.

## 1. Kukelova, Bujnak, Pajdla 2008 — the original automatic generator

### Takeaway
The ECCV 2008 generator automated Stewénius-style hand-crafted Gröbner basis solvers: it computes, on random instances with coefficients in a finite prime field Zp, a fixed "elimination trace" (which monomial multiples of the input polynomials to generate and how to eliminate them, usually with GRevLex), prunes unnecessary rows/columns, and emits MATLAB/Maple code doing one (or a few) G–J eliminations followed by an eigen-decomposition of the action (multiplication) matrix. It gave templates such as 32×48 (8pt radial distortion), 31×46 (6pt shared focal), 179×203 (9pt two distortions) and 238×290 (6pt λ+E+λ), roughly 2× faster than earlier hand-made solvers, but the templates it builds are far from minimal and the method relies on a fixed monomial order, a standard Gröbner basis, and a zero-dimensional ideal.

### Cited Findings
- Bibliographic: ECCV 2008, LNCS 5304, pp. 302–315; generator was distributed at cmp.felk.cvut.cz/minimal — [researchr / DOI](https://www.doi.org/10.1007/978-3-540-88690-7_23); [researchr bibtex](https://researchr.org/publication/KukelovaBP08/bibtex). (When I fetched https://cmp.felk.cvut.cz/minimal/ in Nov 2026 it returned "Page not found".)
- Two-phase method: an offline phase that finds an "Elimination trace" common to all instances, "We do computations in some finite prime field Zp", and an online phase that applies this fixed trace to real data — [KBP08-thesis §4.4](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- Justification of Zp: for a sufficiently large "lucky" prime the trace computed for random Zp coefficients equals, with high probability, the one over Q; exact arithmetic avoids numerical problems and speeds up computation; "In our practical situations a three- or four-digit prime number was sufficient" — [KBP08-thesis §4.4](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- Offline steps (Definition 20 in thesis): (1) fix monomial order ("graded reverse lexicographical ordering is often good"); (2) find the quotient-ring basis B that repeatedly appears for several random coefficient choices, computing in a finite prime field; (3) choose a strategy for generating polynomials from the ideal; (4) for chosen action monomial x^β, systematically generate higher-degree multiples of the input polynomials (Algorithms 8/9) until the polynomials q_i needed for the action matrix are obtained — an "admissible Elimination trace"; (5) remove "unnecessary" polynomials; (6) detect "unnecessary monomials"; (7) the remainder defines the online solver — [KBP08-thesis §4.4](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- Online solver = fill coefficient matrix with data, "one or more G-J eliminations", build multiplication matrix M_x^β from the reduced rows, eigenvectors give the solutions — [KBP08-thesis §4.4.4](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- Implementation: input is a MATLAB script with the equations; it calls Macaulay 2 (or Maple) to compute Gröbner basis / quotient basis over Zp; output is "Matlab or Maple code" — [KBP08-thesis Ch. 6](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- 8pt radial distortion (λ+F+λ, 16 solutions): previous hand-made solver had three G–J eliminations of 8×22, 11×30, 36×50 and a 16×16 eigenproblem; the generator's solver "consists only of one elimination of a 32 × 48 matrix"; MATLAB+MEX runtime ≈640 µs, C++ ≈120 µs (Intel i7 Q720) — [KBP08-thesis §7.1](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- 9pt two radial distortions (λ1+F+λ2, 24 solutions): G–J of 179×203 + 24×24 eigenproblem, ≈8.9 ms MATLAB+MEX (C++ 2.56 ms) vs 18.8 ms for the previous manual solver (LU of 393×390) — [KBP08-thesis §7.1](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- 6pt shared radial distortion (λ+E+λ, 52 solutions): G–J of 238×290 + 52×52 eigenproblem, ≈25.2 ms MATLAB (C++ 7.7 ms) vs 42.9 ms for the manual solver (LU of 320×363) — [KBP08-thesis §7.1](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)
- 6pt relative pose with shared focal (f+E+f, 15 solutions) template from [KBP08] is 31×46 — [Syzygy17 Table 1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf); 15 solutions — [Byrod09 §7.3](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Many later solvers (one-sided focal 6pt 21×30, rolling-shutter pose 48×56, generalized P4P+scale 48×56, P4P+focal+radial 136×152, L2 3-view triangulation 1,866×1,975, R6P 196×216, etc.) were built with this generator — [Syzygy17 Table 1, marked (*)](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- Limitation noted by later authors: the strategy of generating all monomial multiples up to some degree produces large templates; Syzygy17 explicitly flags a "main drawback of the approach in [33]" (sentence truncated in my text extraction; my reading is that it concerns template size/construction) — [Syzygy17 §1.1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf); the thesis itself notes that the generation strategy creates many "unnecessary" polynomials and removing them is order-dependent ("removing one 'unnecessary' polynomial may affect removing other") — [KBP08-thesis §4.4](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf)

### Inferences
- The defining design choices (Zp random instance, fixed GRevLex standard basis, monomial-multiplication expansion, then pruning) are exactly the axes later work attacked: syzygies (which multiples), basis selection (which B), saturation (zero-dimensionality), symmetry, and greedy search.
- Template sizes from KBP08 are the "baseline" in most later comparison tables: e.g. 6pt f+E+f 31×46 → 11×26/12×27 (Greedy22); 9pt λ1+F+λ2 179×203 → 76×100 (Greedy22); 6pt λ+E+λ 238×290 → 39×95 (Greedy22/SparseR).

### Gaps
- Could not access the ECCV 2008 PDF itself (Springer requires login; the CMP page is gone); numbers come from Kukelova's 2013 thesis and from Larsson 2017's comparison table. The original ECCV table of problems/sizes was not verified directly.
- A precise statement of the generator's explicit limitations from the 2008 paper itself was not available.

## 2. Larsson, Åström, Oskarsson CVPR 2017 — syzygy-based reduction ("autogen")

### Takeaway
Instead of expanding all monomial multiples, autogen computes a Gröbner basis over Zp, tracks how each polynomial needed for the action matrix is expressed in the original generators (p_i = Σ h_ij f_j), then reduces those coefficient vectors modulo a Gröbner basis of the first syzygy module, yielding representations of minimal max-degree; the multiplier monomials appearing in the reduced h_ij define the template. This non-iterative reduction gives much smaller templates in most of 40+ problems (e.g. 6pt λ+E+λ 238×290 → 154×210; 9pt λ1+F+λ2 179×203 → 165×200; rel. pose with direction 3pt 411×489 → 210×255, or 40×57 with symmetry) with median generator runtime 5.7 s.

### Cited Findings
- Core idea: "make use of the syzygies (i.e. the polynomial relations between the original equations)"; "Using these syzygies we can essentially parameterize the set of possible elimination templates" — [Syzygy17 abstract](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- Mechanism: compute GB on a Zp instance; normal set gives quotient basis b; reduce reducible monomials r_i·α to get action-matrix rows; track GB elements in terms of f_j (Macaulay2 `ChangeMatrix`) to obtain h_ij with r_i − Σ m_ij b_j = Σ h_ij f_j; compute a GB G_M of the syzygy module M = {s ∈ Zp[X]^n | Σ s_j f_j = 0} (Term-Over-Position-GRevLex) and replace h_i by its normal form h̃_i w.r.t. G_M. Proposition 1: with a degree-first module order, h̃ has minimal max-degree over all h + s, s ∈ M — [Syzygy17 §3.1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- Implementation: MATLAB generator; template generation & reduction "in just a few lines of Macaulay2"; outputs stand-alone MATLAB code; user can vary action monomial and monomial order; automatically detects/exploits variable-aligned symmetries (from ECCV 2016); median generator time over Table 1 problems 5.7 s — [Syzygy17 §3.2](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf). Code "available at http://www.maths.lth.se/~viktorl/" — [CVPR 2017 poster](https://openaccess.thecvf.com/content_cvpr_2017/poster/270-poster.pdf)
- Template sizes (Original → autogen without reduction → with reduction), [Syzygy17 Table 1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf):
  - Rel. pose 5pt (Stewénius): 10×20 → 10×20 → 10×20
  - Rel. pose 8pt one-sided rad. dist. (Kuang): 12×24 → 11×20 → 11×20
  - Rel. pose + one focal 6pt (Bujnak, KBP08-generated): 21×30 → 21×30 → 21×30
  - P3.5P + focal (Wu): 20×43 → 24×45 → 20×44
  - Rel. pose + const. focal 6pt (KBP08): 31×46 → 31×50 → 31×50
  - Rel. pose + rad. dist. 8pt (KBP08): 32×48 → 31×49 → 32×50
  - Rel. pose 6pt one-sided rad. dist. (Kuang): 48×70 → 34×60 → 34×60
  - Rolling shutter pose (Saurer): 48×56 → 50×55 → 47×55
  - Generalized P4P + scale (Ventura): 48×56 → 50×55 → 47×55
  - Stitching + const. focal + rad. dist. 3pt (Naroditsky): 54×77 → 96×108 → 48×66
  - Generalized rel. pose 6pt (Stewénius): 60×120‡ → 135×164 → 99×163
  - Optimal PnP (Hesch): 120×120 → 93×116 → 88×115
  - P4P + focal + rad. dist. (Bujnak, KBP08): 136×152 → 140×144 → 140×156
  - Rel. pose + rad. dist. 6pt (λ+E+λ, KBP08): 238×290 → 223×290 → 154×210
  - Rel. pose + 2 rad. dist. 9pt (KBP08): 179×203 → 355×298 → 165×200
  - Rel. pose 7pt one-sided focal + rad. dist. (Kuang): 200×231 → 249×214 → 185×204
  - Weak PnP: 234×276 → 568×498† → 189×232†; Weak PnP (2×2 sym): 104×90 → 83×90† → 49×59†
  - Rolling shutter R6P (Albl, KBP08): 196×216 → 222×230 → 204×224
  - Optimal pose w dir 4pt (Svärm): 280×252 → 371×351 → 203×239
  - Rel. pose w dir 3pt (Saurer, KBP08): 411×489 → 287×324 → 210×255; using symmetry: 94×111† → 40×57†
  - Abs. pose quivers (Kuang): 372×386 → 420×406 → 217×253
  - L2 3-view triangulation relaxed (KBP08): 274×305 → 399×384 → 239×290
  - Rel. pose w angle 4pt (Li): 270×290 → 280×304 → 266×329
  - Refractive P5P (Haner): 280×399 → 410×480 → 240×324
  - TDOA offset rank 3 (6,8 pts): 1359×754 → 1359×754 → 356×345
  - Optimal PnP (Zheng): 575×656 → 812×704 → 521×601
  - Optimal pose w dir 3pt (Svärm): 1,260×1,278 → 918×726 → 544×592
  - Refractive P6P + focal (Haner): 648×917 → 2,196×1,913† → 636×851†
  - Rel. pose + const. focal + rad. dist. 7pt (Jiang): 886×1,011 → 1,393×1,237 → 581×862
  - Dual-receiver TDOA 5pt: 2,625×2,352 → 850×1,167 → 455×768
  - Optimal PnP (rot. matrix, Nakano): 1,936×1,976 → 1,698×1,153 → 1,102×1,135
  - L2 3-view triangulation (KBP08): 1,866×1,975 → 2,647×2,584 → 1,759×2,013
  (†: variable-aligned symmetry auto-removed; ‡: original does not compute full GB)
- Rel. pose with known rotation direction 3pt: 411×489 → 40×57 "corresponds to a speed-up factor of 900" assuming cost quadratic in rows and linear in columns — [Syzygy17 §4](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- Numerics: all generated solvers have log10 residual mode below −4.8; median of the log10-error modes is −10.9; on stitching, Hesch PnP, Zheng PnP the accuracy is similar to originals; no template refinement (e.g. column pivoting) used — [Syzygy17 §5](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- The reduction "does not guarantee that the template will be smaller", and only one minimal-degree representation is found — [Syzygy17 §3.1, §4](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- GAPS (Bo Li & Viktor Larsson 2020) "wraps and improves autogen" with a friendlier interface, more functionality and better stability; repo https://github.com/prclibo/gaps; other generators listed: PavelTrutman/Automatic-Generator (reimplementation of KBP08-style generator), laurentkneip/polyjam (C++) — [GAPS20](https://arxiv.org/pdf/2004.11765)

### Inferences
- For the specific problems in the task: 6pt shared focal stayed at 31×50 (no gain vs. KBP08's 31×46); 9pt two radial 179×203 → 165×200 (small gain); 7pt one-sided focal+radial 200×231 → 185×204; P4Pfr 136×152 → 140×156 (no gain); 8pt radial 32×48 → 32×50 (no gain); 6pt λ+E+λ 238×290 → 154×210 (large gain). Big gains came on hard problems (TDOA, rel. pose with direction, dual-receiver TDOA).
- Absolute pose with unknown focal (P4Pf, without distortion) is not in Syzygy17 Table 1; only P3.5P+focal (20×44) and P4P+fr.

### Gaps
- Exact current GitHub URL of Larsson's original autogen could not be verified (GitHub returned 403 to this session; the poster points to the author's Lund homepage). GAPS repo URL is from the GAPS paper.

## 3. Larsson et al. ICCV 2017 — polynomial solvers for saturated ideals

### Takeaway
They extend the action matrix method to the saturation J = I : f_s^∞, removing unwanted (possibly positive-dimensional) solution components such as trivial/degenerate ones without ever computing generators of J online: basis and action polynomials are computed in J over Zp and then "lifted" back to I by multiplying with f_s^N. This made problems with infinitely many spurious solutions tractable and gave large speedups: vanishing-point estimation 2,860×3,060 / 1.8 s → 246×397 / 5 ms; TOA (4,6) 966×925 / 0.6 s → 569×692 / 22 ms; 3-view triangulation 209×265 template with 50 solutions in 3 ms.

### Cited Findings
- Goal: create solvers when "a (possibly infinite) subset of the solutions are undesirable"; "extends the standard action matrix method to saturated ideals. This allows us to add constraints that some polynomials should be non-zero on the solutions"; previously problems with infinitely many solutions could not be solved directly since action-matrix requires a zero-dimensional ideal — [Sat17 abstract](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf)
- Theory: for basis {b_k} of C[X]/J, the set {f_s^N b_k} is linearly independent in C[X]/I (Lemma 1), and S_N = span{f_s^N b_k} is stable under multiplication by α for large N (Lemma 2); so the action matrix of J can be computed with polynomials of I: they "lift both the basis {b_k} and the polynomials p_i into the original ideal I by multiplying with f_s^N" — [Sat17 §3–4](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf)
- Implementation extends the syzygy generator [24]; steps carried out in Macaulay2; restricted to saturating a single monomial; to saturate a general polynomial f, add variable x0 with x0 − f(x) = 0 and use an order in which monomials with x0 are larger so x0 is not in the basis (gives much smaller templates) — [Sat17 §4](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf)
- Results [Sat17 Table 1](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf):
  - Three-view triangulation: Stewénius — 47 sols, >1 s; Byröd (relaxed ideal) 225×209†, 154 sols, 60 ms; Byröd (basis select.) 225×209†, 58, 3 ms; Ours 571×676, 47 sols, 10 ms; Ours (new param.) 209×265, 50 sols, 3 ms.
  - Time-of-arrival: Kuang (4,6) 966×925†, 38 sols, 0.6 s; (5,5) 1,386×1,539†, 42, 1.4 s; Ours (4,6) 569×692, 38, 22 ms; Ours (5,5) 938×1,301, 42, 95 ms.
  - Optimal vanishing point: Mirzaei 2,860×3,060, 40 sols, 1.8 s; Ours 246×397, 40 sols, 5 ms.
- Vanishing point solver: average 5 ms vs 1.8 s on York Urban Dataset, "a speed-up of more than a factor of 300, which shows that ... it is much more beneficial to saturate directly compared to introducing an auxiliary variable" — [LarssonThesis §2.4](https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf)

### Inferences
- Saturation is the standard tool for removing spurious solutions introduced by modeling (e.g. zero-scale, depth=0, f=0); Martyushev's 2023 Laurent generator later handles positive-dimensional components directly and compares against saturated formulations (Rabinowitsch trick or G–J cascades).

### Gaps
- Sat17 does not evaluate classic small minimal problems (6pt focal etc.) with saturation, so no direct numbers for those.

## 4. Beyond Gröbner Bases (CVPR 2018) — basis selection and random basis sampling

### Takeaway
Standard monomial bases from different monomial orders (Gröbner fan, via Gfan) give very different templates; going further, any monomial set that spans the quotient ring and admits an action matrix can be used ("beyond Gröbner bases"). A heuristic random sampler of such bases found much smaller templates: e.g. F+λ 8pt 11×20 → 7×16; E+f 6pt 21×30 → 11×20; f+E+f 6pt 31×50 → 21×40; E+λ 6pt 34×60 → 14×40; stitching 3pt 48×66 → 18×36 (0.38 → 0.16 ms); P4Pfr 140×156 → 54×70; λ+E+λ 6pt 149×205 → 53×115; λ1+F+λ2 9pt 165×200 → 84×117; E+fλ 7pt 181×200 → 69×90.

### Cited Findings
- Abstract: "careful selection of the monomial bases"; "we describe how we can enumerate all such bases in an efficient way. We also show that going beyond Gröbner bases leads to more efficient solvers in many cases. We present a novel basis sampling scheme" — [BeyondGB18](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)
- Gröbner fan: finitely many reduced GBs per ideal (Mora–Robbiano); they computed fans with 12-hour abort; for 6pt shared radial distortion the fan computation ran a week before aborting; numbers of reduced GBs: e.g. P4Pfr 1745, 9pt λ1+F+λ2 6896, refractive P5P 8659 — [BeyondGB18 §2.2, Table 1](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)
- Sampling heuristic: random binary weight vector ω penalizing monomial degrees, iteratively sample basis monomials proportionally to weights from neighbours; 100 sampled bases per problem — [BeyondGB18 §3.1–3.3](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)
- [BeyondGB18 Table 1](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf) (Original | autogen [29] GRevLex | GFan+[29] (#GB) | Heuristic sampling+[29]):
  - Rel. pose F+λ 8pt (Kuang): 12×24 | 11×20 | 11×20 (10) | 7×16
  - Rel. pose E+f 6pt (Bujnak): 21×30 | 21×30 | 11×20 (66) | 11×20
  - Rel. pose f+E+f 6pt (Kukelova): 31×46 | 31×50 | 31×50 (218) | 21×40
  - Rel. pose E+λ 6pt (Kuang): 48×70 | 34×60 | 34×60 (846) | 14×40
  - Stitching fλ+R+fλ 3pt (Naroditsky): 54×77 | 48×66 | 48×66 (26) | 18×36
  - Abs. pose P4Pfr (Bujnak): 136×152 | 140×156 | 54×70 (1745) | 54×70
  - Rel. pose λ+E+λ 6pt (Kukelova): 238×290 | 149×205 | – | 53×115
  - Rel. pose λ1+F+λ2 9pt (Kukelova): 179×203 | 165×200 | 84×117 (6896) | 84×117
  - Rel. pose E+fλ 7pt (Kuang): 200×231 | 181×200 | 69×90 (3190) | 69×90
  - Rel. pose E+fλ 7pt (elim. λ): – | 52×71 | 37×56 (332) | 24×43
  - Rel. pose E+fλ 7pt (elim. fλ, Kukelova): 51×70 | 51×70 | 51×70 (3416) | 51×70
  - Abs. pose quivers (Kuang): 372×386 | 216×258 | – | 81×119
  - Rel. pose E angle+4pt (Li): 270×290 | 266×329 | – | 183×249
  - Abs. pose refractive P5P (Haner): 280×399 | 240×324 | 157×246 (8659) | 240×324
- Stitching 3pt timing (MATLAB+MEX): proposed 0.16 ms vs autogen 0.38 ms vs Byröd 0.89 ms — [BeyondGB18 Table 2](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)
- For E+fλ 7pt the best (69×90) is still larger than Kukelova's elimination-ideal formulation 51×70 — [BeyondGB18 §5](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)
- Note on table discrepancy: autogen GRevLex numbers differ slightly between papers (e.g. λ+E+λ 6pt 154×210 in Syzygy17 vs 149×205 in BeyondGB18; E+fλ 7pt 185×204 vs 181×200), presumably from generator updates — compare [Syzygy17 Table 1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf) and [BeyondGB18 Table 1](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf)

### Inferences
- Basis choice is a first-order lever: 2–4× fewer rows on many classic problems with no change in formulation; standard GRevLex is often not optimal.

### Gaps
- The paper reports template sizes, but stability comparisons only for selected problems (residual histograms for stitching); no systematic runtime table across all problems.

## 5. Byröd, Josephson, Åström (IJCV 2009) and related — numerical stability

### Takeaway
Byröd et al. generalized the action matrix method to arbitrary (non-Gröbner, even redundant) linear bases of C[x]/I and chose the basis at run time via QR with column pivoting or SVD, plus adaptive truncation; this improved accuracy by up to 10^7 (median/95th-percentile) on optimal 3-view triangulation, making it feasible in double precision (earlier it needed emulated 128-bit arithmetic), at ≈1.4× the cost of the standard solver (QR) or ≈5× (SVD).

### Cited Findings
- Contributions: "show how the action matrix can be computed in the general setting of an arbitrary linear basis for C[x]/I"; two basis-selection strategies — QR factorization with column pivoting and SVD; "adaptive scheme for truncation of the Gröbner basis"; "dramatically improved numerical stability" — [Byrod09 abstract](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Diagnosis: main source of numerical problems is the conditioning of a crucial elimination step; earlier 3-view triangulation needed "emulated 128 bit numerics" — [Byrod09 §1](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Redundant (overly large) basis avoids the sensitive elimination step and yields extra false solutions filtered by evaluation; QR pivots sorted descending give adaptive truncation by thresholding the ratio of largest/smallest pivot (τ in 10^6–10^10 gives similar results) — [Byrod09 §1, §5.3](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- 3-view triangulation (100,000 trials), 95th percentile / median error: STD 1.42·10^1 / 9.85·10^-5; SVD 1.19·10^-5 / 6.09·10^-9; SVD:eig 1.20·10^-6 / 1.29·10^-9; QR 1.78·10^-5 / 1.06·10^-8; QR:var+fast 3.41·10^-6 / 3.61·10^-9; "The improvement in precision is up to a factor 10^7" — [Byrod09 Table 3](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Adaptive basis size: minimal (50) in 94% (SVD) / 95% (QR) of cases — [Byrod09 Table 1](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Timing per call (Intel Core 2 2.13 GHz, MATLAB, solver part): SVD 41.69 ms, TRUNC 38.11 ms, QR:var+fast 10.94 ms, STD 8.03 ms — QR ≈4× faster than SVD but 40% slower than standard — [Byrod09 Table 4](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Localization with hybrid features (36 solutions): no standard solver could be built (C not numerically invertible even at higher degrees); with QR/SVD + truncation (basis 39+) it worked — [Byrod09 §7.2](https://www.maths.lth.se/vision/publdb/reports/pdf/byrod-josephson-etal-ijcv-09.pdf)
- Related earlier work (ICCV 2007 "Improving numerical accuracy...", ECCV 2008 column pivoting) and Bujnak's Sturm-sequence extraction are summarized in [Sat17 §1](https://openaccess.thecvf.com/content_ICCV_2017/papers/Larsson_Polynomial_Solvers_for_ICCV_2017_paper.pdf); Martyushev 2023 compares to Byröd STD and QR solvers for 3-view triangulation: 1.34 ms (ours, Laurent) vs 1.54 ms [STD] vs 2.07 ms [QR] vs 1.56 ms [Sat17]; fail (error>1) 1.11% / 3.09% / 0.29% / 8.91% — [Laurent23 Table 2](https://arxiv.org/pdf/2307.00320)

### Inferences
- Column-pivoting QR remains the most stable option for hard problems (lowest fail rate 0.29% for triangulation) at a moderate time cost; Greedy22 notes that the 581×659 autogen solver for 7pt fλ+E+fλ "is unstable and requires additional stability improving techniques, e.g., column pivoting".

### Gaps
- The IJCV paper's numbers for 6pt focal relative pose (34×50 template, 15 solutions) were seen but no error statistics extracted.

## 6. Martyushev, Vrablikova, Pajdla CVPR 2022 (greedy parameter search) and follow-ups

### Takeaway
Using a partial affine parameterization of templates derived from the syzygy-based parameterization, plus greedy heuristic search over parameters (and over many standard and non-standard bases and action variables, plus a Schur-complement reduction), they matched or beat state of the art on all 34 tested problems — up to 2.1, 2.5, 3.8, 6.6× smaller, and 20× for refractive P6P+focal (648×917 → 126×162 in std basis per their table). The 2023 Laurent-polynomial generator (Maple) further shrank templates/runtimes on 36 problems and handles positive-dimensional components and partial p-fold symmetries automatically.

### Cited Findings
- Abstract: affine parameterization + greedy heuristic search; 34 problems; "For all of them, we found the templates either of the same or smaller size"; "e.g., 2.1, 2.5, 3.8, 6.6 times smaller"; refractive absolute pose with unknown focal "20 times smaller"; new solver for relative pose with unknown common focal + radial distortion — [Greedy22](https://arxiv.org/abs/2203.14901)
- Method: "using the general syzygy-based parameterization of elimination templates from [33], we construct a partial (but still generic enough) parameterization of templates. Then, we apply a greedy heuristic optimization" — [Greedy22 §1.1](https://arxiv.org/pdf/2203.14901)
- Code: "The implementation of the new AG, as well as the Matlab solvers for all the minimal problems from Tab. 1 and Tab. 2, are available at http://github.com/martyushev/EliminationTemplates" — [Greedy22 §5](https://arxiv.org/pdf/2203.14901)
- [Greedy22 Table 1](https://arxiv.org/pdf/2203.14901) (d | Our std | Our nstd | Syzygy [33] | BeyondGB std | BeyondGB nstd | SparseR [5]):
  1. Rel. pose F+λ 8pt: d=8 | 11×19 | 7×15 | 11×19 | 11×19 | 7×15 | 7×16
  2. Rel. pose E+f 6pt: d=9 | 11×20 | 11×20 | 21×30 | 11×20 | 11×20 | 11×20
  3. Rel. pose f+E+f 6pt: d=15 | 12×27 | 11×26 | 31×46 | 31×46 | 21×40 | 12×30
  4. Rel. pose E+λ 6pt: d=26 | 34×60 | 14×40 | 34×60 | 34×60 | 14×40 | 14×40
  5. Stitching fλ+R+fλ 3pt: d=18 | 48×66 | 18×36 | 48×66 | 48×66 | 18×36 | 18×36
  6. Abs. pose P4P+fr: d=16 | 52×68 | 52×68 | 140×156 | 54×70 | 54×70 | 52×68
  7. P4P+fr (el. f): d=12 | 28×40 (all methods 28×40)
  8. Rel. pose λ+E+λ 6pt: d=52 | 73×125 | 39×95 | 149×201 | – | 53×105 | 39×95
  9. Rel. pose λ1+F+λ2 9pt: d=24 | 76×100 | 76×100 | 165×189 | 87×111 | 87×111 | 90×117
  10. Rel. pose E+fλ 7pt: d=19 | 55×74 | 56×75 | 185×204 | 69×88 | 69×88 | 61×80
  11. E+fλ 7pt (el. λ): d=19 | 37×56 | 22×41 | 52×71 | 37×56 | 24×43 | 22×41
  12. E+fλ 7pt (el. fλ): d=19 | 51×70 (all)
  13. Rolling shutter pose: d=8 | 47×55 (all)
  14. Triangulation (sat. im.): d=27 | 87×114 | 87×114 | 88×115 | 88×115 | 88×115 | 87×114
  15. Abs. pose refractive P5P: d=16 | 57×73 | 57×73 | 240×256 | 112×128 | 199×215 | 68×93
  16. Abs. pose quivers: d=20 | 65×85 | 66×86 | 169×189 | – | 68×88 | 68×92
  17. Unsynch. rel. pose: d=16 | 159×175 | 139×155 | 159×175 | – | 299×315 | 150×168
  18. Optimal PnP (Hesch): d=27 | 87×114 | 87×114 | 88×115 | 88×115 | 88×115 | 87×114
  19. Optimal PnP (Cayley): d=40 | 118×158 (all)
  20. Optimal pose 2pt v2: d=24 | 139×163* | 141×165* | 192×216 | – | 192×216 | 176×200
  21. Rel. pose E+angle 4pt: d=20 | 99×119* | 99×119* | 246×276 | – | 183×249 | –
  (*: reduced by Schur-complement method of Subsect. 4.2)
- [Greedy22 Table 2](https://arxiv.org/pdf/2203.14901) (d | Our std | Our nstd | Original | Syzygy): λ+F+λ 8pt d=16 31×47 | 31×47 | 32×48 | 32×48; P3.5P+focal d=10 18×28 | 19×29 | 20×43 | 20×30; Gen. P4P+scale d=8 47×55 | 47×55 | 48×56 | 47×55; Rolling shutter R6P d=20 120×140 | 120×140 | 196×216 | 204×224; Opt. pose w dir 4pt d=28 134×162* | 144×172* | 280×252 | 203×231; Opt. pose w dir 3pt d=48 397×445* | 385×433* | 1260×1278 | 544×592; L2 3-view triang. (relaxed) d=31 217×248 | 281×312 | 274×305 | 231×262; Refractive P6P+focal d=36 126×162 | 178×214 | 648×917 | 636×654; Rel. pose fλ+E+fλ 7pt d=68 209×277 | 255×323 | 886×1011 | 581×659; Gen. rel. pose+scale 7pt d=140 144×284 (all); Rel. pose E+angle 4pt v2 d=20 16×36*; Gen. rel. pose E+angle 5pt d=44 37×81*.
- 7pt shared focal + radial distortion (fλ+E+fλ, 68 solutions): Gröbner fan intractable (>24 h), so sampled 1,000 weighted orders + 500 heuristic bases; best 209×277 (weights w=(135,81,98,107,68), action variable λ); 8.5 ms MATLAB vs ~16 ms for 581×659 autogen solver and 400 ms reported for Jiang et al. 886×1011 solver; ~2% failures; the only publicly available SOTA solver (their ref. [48], main parts in C++) runs at 2.9 ms — [Greedy22 §5.1](https://arxiv.org/pdf/2203.14901)
- Laurent generator (2023): new algorithm checks whether a set of Laurent polynomials suffices to construct a template; "simple and fast; it applies to ideals with positive-dimensional components; it allows one to uncover partial p-fold symmetries automatically"; speed "exceeds the state-of-the-art in most cases"; Maple implementation + MATLAB solvers at github.com/martyushev/EliminationTemplates — [Laurent23](https://arxiv.org/pdf/2307.00320)
- [Laurent23 Table 1](https://arxiv.org/pdf/2307.00320), Ours (size, roots, ms) vs Greedy22 (size, roots, ms), MATLAB, i5-1155G7, times excl. coefficient-matrix construction:
  - λ+F+λ 8pt: 19×39, 20, 0.2 vs 31×47, 16, 0.2
  - P3.5P+focal: 12×26, 14, 0.1 vs 18×28, 10, 0.1
  - Stitching 3pt: 6×30, 24, 0.2 vs 18×36, 18, 0.2
  - P4P+fr: 42×60, 18, 0.2 vs 52×68, 16, 0.3
  - λ1+F+λ2 9pt: 73×97, 24, 0.3 vs 76×100, 24, 0.4
  - E+fλ 7pt: 34×56, 22, 0.2 vs 55×74, 19, 0.6
  - Refractive P5P: 38×58, 20, 0.2 vs 57×73, 16, 0.3
  - Unsynch. rel. pose: 59×79, 20, 0.3 vs 139×155, 16, 0.6
  - R6P: 66×92, 26, 0.3 vs 120×140, 20, 0.9
  - Opt. pose w dir 3pt: 297×356, 59, 3.0 vs 385×433, 48, 6.8
  - fλ+E+fλ 7pt: 138×210, 72, 2.3 vs 209×277, 68, 3.0
  - λ1+E+λ2 7pt: 204×289, 85, 4.2 vs 436×512, 76, 27.9
  - Gen. rel. pose 6pt: 78×155, 77, 2.5 vs 99×163, 64, 2.5
  - Rel. pose 9 lines (36 sols): 1,610×1,726, 116 roots, 205 ms (Greedy failed; original attempt 16,278×13,735)
  - Optimal PnP (Zheng): 120×175, 55, 1.5 vs 272×312, 40, 5.9
  - Refractive P6P+focal: 40×61, 21*, 0.2 vs 126×162, 36, 0.8
  - Vanishing point: 136×194†, 58, 2.1 vs 343×383, 40, 13.8
  - TOA (4,6): 427×475†, 48, 5.2 vs 863×901, 38, 17.1; TOA (5,5): 772×832†, 60, 17.8 (Greedy failed)
  - Optimal 3-view triangulation: 69×127‡, 58, 1.1 vs 401×448, 47, 5.7
  (Note: Laurent23 templates often have more roots than d, i.e. extra spurious eigenvalues traded for smaller templates; numerical log10 errors are mostly comparable, but for some problems worse, e.g. Gen. rel. pose+scale 7pt −3.16 vs −5.93 mean.)
- Template size alone is not an adequate efficiency measure: e.g. refractive P5P 38×58 template + 20×20 eig vs 57×73 + 16×16 eig — [Laurent23 §5](https://arxiv.org/pdf/2307.00320)

### Inferences
- Across the classic problems, the best known action-matrix templates (as of 2022–23) are roughly: 6pt one-sided focal 11×20; 6pt shared focal 11×26; 8pt one-sided radial 7×15; 6pt one-sided radial 14×40; 6pt shared radial (λ+E+λ) 39×95; 9pt two radial 76×100 (Greedy) / 73×97 (Laurent); 7pt one-sided focal+radial 55×74 / 34×56; P4Pfr 52×68 / 42×60; 3pt stitching 18×36 / 6×30.
- Greedy22's "Syzygy [33]" column sometimes differs from Larsson's own published numbers (e.g. f+E+f 6pt: 31×46 in Greedy22 vs 31×50 in Syzygy17; 9pt: 165×189 vs 165×200), likely because Greedy22 re-ran autogen and/or removed redundant columns; the report writer should cite the specific table.

### Gaps
- I did not find a verified 2024–2026 peer-reviewed venue for Laurent23 (it is an arXiv preprint as of my search), nor reliable newer (2024–2026) Gröbner-template-optimization papers; a web search for 2024–2025 work returned nothing beyond these. This may be a search limitation.
- GitHub repository content (languages, licenses) could not be verified directly (GitHub blocked for this session); language info comes from the papers (Greedy22: MATLAB solvers; Laurent23: Maple generator + MATLAB solvers).

## 7. Symmetry exploitation (Ask et al., Kuang et al., Larsson & Åström ECCV 2016)

### Takeaway
If every equation has a (partial / weighted) p-fold symmetry, the solution set is invariant under x ↦ ωx and the action matrix can be built on one remainder class of basis monomials, shrinking both template and eigenproblem. Larsson & Åström (ECCV 2016) generalized partial p-fold symmetry (Ask et al. 2012, Kuang et al. CVPR 2014) to weighted and multiple independent symmetries and to "unaligned" symmetries revealed by a linear change of variables; e.g. Weak PnP: 231×263 template / 32×32 action / 1.9 ms → 26×34 / 8×8 / 0.22 ms with 2×2-fold symmetry.

### Cited Findings
- ECCV 2016 paper: Larsson & Åström, "Uncovering Symmetries in Polynomial Systems", LNCS 9907, pp. 252–267, DOI 10.1007/978-3-319-46487-9_16 — [Lund record](https://lup.lub.lu.se/record/1506a56c-c165-4b88-89cd-9dbfc521c5d5)
- "The main contribution is a generalization of the partial p-fold symmetry from Ask et al. [13] and Kuang et al. [119]"; Ask et al. considered systems where each monomial degree has the same remainder mod p; Kuang et al. extended to symmetry in a subset of variables; Larsson & Åström handle multiple independent symmetries and weighted (c-weighted p-fold) symmetries — [LarssonThesis Ch. 3](https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf)
- Implementation in autogen: automatically finds variable-aligned symmetries by computing exponent remainders for all weight vectors c up to a maximum weight; partitions the quotient basis into remainder classes, picks an action monomial invariant to the symmetry, builds a template per class, keeps the smallest — [LarssonThesis §3.2.2](https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf)
- Weak PnP (33 solutions incl. trivial; 8 groups of 4): template / action matrix / runtime: 2×2-sym 26×34 / 8×8 / 0.22 ms; 2-sym 138×154 / 16×16 / 0.79 ms; no sym 231×263 / 32×32 / 1.9 ms; Wu's reparameterization 208×220 / 12×12 / 0.91 ms. KBP08 generator without symmetry gave 243×273 — [LarssonThesis Table 3.1](https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf)
- Syzygy17 automatically removed symmetries: Rel. pose w dir 3pt 210×255 → 40×57; Optimal PnP (Zheng) 521×601 → 302×342; Weak PnP 2×2 sym 49×59 — [Syzygy17 Table 1](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf)
- Laurent23 uncovers partial p-fold symmetries automatically "by constructing the solving set of cardinality less than the degree of the related ideal" (Weak PnP 2×2 sym: 10×20, 10 roots; Refractive P6P+focal: 40×61, 21 roots vs d=36) — [Laurent23 §5 remark 6, Table 1](https://arxiv.org/pdf/2307.00320)
- Greedy22 notes that refractive P5P (d=16) has a 2-fold symmetry and refractive P6P+focal a 4-fold symmetry that could further reduce templates — [Greedy22 §5 remarks](https://arxiv.org/pdf/2203.14901)

### Inferences
- Symmetry reduction is often the single largest saving (≈5–10× in template elements) when present, and is now automated in autogen and in Martyushev's Laurent generator.

### Gaps
- Kuang, Zheng, Åström CVPR 2014 "Partial symmetry in polynomial systems" and Ask et al. 2012 numbers were not retrieved directly; the ECCV 2016 paper PDF was not obtained (thesis chapter used).

## 8. Other notable improvements and hybrids (resultants, hidden variable, sparse resultants)

### Takeaway
Sparse-resultant generators (Bhayani, Kukelova, Heikkilä, CVPR 2020; arXiv 2023 journal version) using an extra (u-resultant-like) polynomial and Schur complement produce templates that equal or beat Gröbner/action-matrix templates on many problems (e.g. 6pt f+E+f 12×30, λ+E+λ 39×95, P4Pfr 52×68), and the 2023 paper proves conditions under which they are equivalent to action-matrix solvers; a hidden-variable resultant variant (ICPR 2020) trades a larger generalized eigenvalue problem for better stability.

### Cited Findings
- SparseR (CVPR 2020) proposes "a novel approach to convert the resultant constraint to an eigenvalue problem", leading for most problems "to solvers that are the same size as the best available" Gröbner solvers — [SparseR20](https://openaccess.thecvf.com/content_CVPR_2020/papers/Bhayani_A_Sparse_Resultant_Based_Method_for_Efficient_Minimal_Solvers_CVPR_2020_paper.pdf); its sizes are listed in the SparseR column of [Greedy22 Table 1](https://arxiv.org/pdf/2203.14901) (e.g. F+λ 8pt 7×16, f+E+f 6pt 12×30, λ+E+λ 39×95, λ1+F+λ2 9pt 90×117, E+fλ 7pt 61×80, P4Pfr 52×68, refractive P5P 68×93)
- Hidden variable (ICPR 2020): sparse resultant with one variable hidden; "more stable solvers than Gröbner basis solvers and earlier resultant-based solvers, especially near critical configurations"; larger eigenvalue problem but avoids inverting/eliminating large matrices; can be automated — [Hidden20](https://arxiv.org/abs/2007.10100)
- 2023 journal version: iterative scheme testing subsets of Newton polytopes; extra polynomial of special form + Schur complement; "for some camera geometry problems our extra polynomial-based method leads to smaller and more stable solvers than the state-of-the-art Gröbner basis-based solvers"; step-by-step equivalence conditions with action-matrix solvers; code "Automatic Generator for Sparse Resultant Solvers" at https://github.com/snehalbhayani/aut_gen... — [SRConn23](https://arxiv.org/pdf/2301.06443)
- Schur-complement reduction is also used by Greedy22 (Subsect. 4.2, marked * in its tables, e.g. Rel. pose E+angle 4pt 246×276 → 99×119) — [Greedy22](https://arxiv.org/pdf/2203.14901)

### Inferences
- The boundary between "Gröbner/action-matrix" and "sparse-resultant" generators has blurred: both produce a Macaulay-type template + (generalized) eigenproblem; the best numbers for many problems are tied between Greedy22 nstd and SparseR.

### Gaps
- Did not extract SRConn23's own tables (2023 sizes) or Hidden20 quantitative stability results.
- Did not research homotopy-continuation / learned-start-system alternatives (Hruby et al.), considered out of scope for Gröbner generators.

## 9. Open-source code availability (summary)

### Takeaway
All major lines have public code, mostly MATLAB-output generators relying on a CAS (Macaulay2 for Kukelova and Larsson; Maple for Martyushev 2023).

### Cited Findings
- KBP08: MATLAB input script; uses Macaulay 2 or Maple; outputs MATLAB or Maple solver code; distributed at cmp.felk.cvut.cz/minimal (now 404) — [KBP08-thesis Ch. 6](https://cmp.felk.cvut.cz/~kukelova/webthesis/docs/Kukelova-phd-2013.pdf); a third-party reimplementation exists at https://github.com/PavelTrutman/Automatic-Generator — [GAPS20](https://arxiv.org/pdf/2004.11765)
- Larsson autogen (Syzygy17/Sat17/BeyondGB18): MATLAB generator + Macaulay2 for template computation and reduction; "We have made the code for generating solvers publicly available" — [Syzygy17 §3.2](https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf); link http://www.maths.lth.se/~viktorl/ — [poster](https://openaccess.thecvf.com/content_cvpr_2017/poster/270-poster.pdf)
- GAPS (Li & Larsson 2020), wraps/improves autogen: https://github.com/prclibo/gaps — [GAPS20](https://arxiv.org/pdf/2004.11765)
- PolyJam (Kneip), C++ symbolic toolbox / GB solver generator: https://github.com/laurentkneip/polyjam — [GAPS20](https://arxiv.org/pdf/2004.11765); [LarssonThesis ref 108](https://lucris.lub.lu.se/ws/files/42748556/thesis.pdf)
- Martyushev et al.: http://github.com/martyushev/EliminationTemplates — Greedy22 generator + MATLAB solvers for 34 problems [Greedy22](https://arxiv.org/pdf/2203.14901); Maple implementation of Laurent generator + MATLAB solvers for 36 problems [Laurent23](https://arxiv.org/pdf/2307.00320)
- Bhayani et al. sparse resultant generator: https://github.com/snehalbhayani/aut_gen (URL truncated in PDF text extraction as "aut gen") — [SRConn23](https://arxiv.org/pdf/2301.06443)
- Gfan used for Gröbner fans in BeyondGB18 and Greedy22 — [BeyondGB18](https://openaccess.thecvf.com/content_cvpr_2018/papers/Larsson_Beyond_Grobner_Bases_CVPR_2018_paper.pdf); [Greedy22 footnote 3](https://arxiv.org/pdf/2203.14901)

### Inferences
- Current practical stack: Macaulay2 (autogen/GAPS) or Maple (Martyushev 2023) offline; MATLAB (or C++ via GAPS/PoseLib-style ports) online.

### Gaps
- Could not verify repository existence, licenses or languages directly on GitHub (access blocked); URLs are as cited in the papers. Kukelova's generator download page is no longer reachable.
