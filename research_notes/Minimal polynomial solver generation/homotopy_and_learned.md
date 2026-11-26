# Homotopy Continuation and Learning-Based Approaches to Minimal Problems

Primary sources read in full text: MINUS / TRPLP CVPR 2020 paper (CVF open access PDF), Hruby et al. "Learning to Solve Hard Minimal Problems" (arXiv 2112.03424 PDF, the CVPR 2022 version), GPU-HC CVPR 2022 paper (CVF PDF), the GPU-HC GitHub README, Simulator HC (arXiv 2411.03745 HTML) and Neural Predictor-Corrector (arXiv 2602.03086 HTML). Everything else comes from search snippets or abstracts and is marked as such.

Source URLs used below:
- MINUS CVPR 2020: https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf
- Hruby CVPR 2022: https://arxiv.org/abs/2112.03424 (CVF version: https://openaccess.thecvf.com/content/CVPR2022/papers/Hruby_Learning_To_Solve_Hard_Minimal_Problems_CVPR_2022_paper.pdf)
- GPU-HC CVPR 2022: https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf
- GPU-HC repo: https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU
- Simulator HC: https://arxiv.org/html/2411.03745
- NPC: https://arxiv.org/html/2602.03086v1

---

## Q1. MINUS: "Trifocal Relative Pose from Lines at Points" (Fabbri et al., CVPR 2020; TPAMI 2023)

### Takeaway
MINUS is a templated C++ coefficient-parameter homotopy continuation (HC) solver. It solves the "Chicago" problem (degree 312) and the "Cleveland" problem (degree 216) for trifocal calibrated relative pose. Each solve takes 660 ms on average (1.9 s worst case) on an Intel Core i7-7920HQ with 4 threads, against more than 1 minute for general-purpose HC software. The failure rate is about 1%. It was the first practical solver at this degree, but it is too slow for RANSAC.

### Cited Findings
- **Problems:**
  - "Chicago": three point triplets, two of which carry orientations (tangents, i.e. lines at points). It is minimal with up to 312 solutions; "the area code of Chicago is 312". Typically 3–4 of these are relevant. — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - "Cleveland" (the straight-line variant) has 216 solutions. The paper states "the Chicago problem is minimal with algebraic degree 312, and the Cleveland problem is minimal with algebraic degree 216". — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - A search snippet of the expanded preprint says D = 312 was obtained from a Gröbner basis computation. — [ND preprint](https://www3.nd.edu/%7Ejhauenst/preprints/fdfrdtwhgklpTrifocalRelativePose.pdf)
- **Method:**
  - Coefficient-parameter homotopy that "naturally exploits the algebraic degree of the problem". On top of it, the authors optimized polynomial evaluation and numerical linear algebra in C++. — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - The tracking "iteratively tracks solutions with a guarantee of global convergence". — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - The most reliable formulation uses 14×14 minors. Other formulations "demonstrated further potential for speedup by orders of magnitude" but may trade off reliability. — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
- **Timing:**
  - "Each solve using our software MINUS takes 660ms on average (1.9s in the worst case) as compared to over 1 minute on average for the best prototypes using general purpose software [Bertini, Macaulay2], both on an Intel core i7-7920HQ processor and four threads." — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - Because every path is independent, the solvers are suitable for parallel implementation. "More aggressive but potentially unsafe optimizations towards microseconds are feasible". — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
- **Synthetic results** (1000 runs):
  - On average 35 real solutions, 7 after the positive-depth check, and about 3–4 after verification with the unused third tangent. — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - The solver "fails in about 1% of cases, which are detectable". — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
- **Real-data motivation:** COLMAP fails on 75% of the authors' challenging triplets. The trifocal method needs only three matched features. — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
- **Code:** http://github.com/rfabbri/minus ("MInimal problem NUmerical Solver"). — [MINUS CVPR 2020](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
- **Journal version:** IEEE TPAMI 45(6):7870–7884 (from search snippets). — [NSF PAR](https://par.nsf.gov/biblio/10387070)
- **Hruby et al. on MINUS:** it is "globally convergent with probability 1, but needs about 500 miliseconds to track 312 complex solution paths". It "demonstrated that hard multiview problems involving 312 and 216 solutions can be solved in a stable way" but is not practical for RANSAC. — [Hruby et al. 2022](https://arxiv.org/abs/2112.03424)
- **GPU-HC timings for the same problem:** "Trifocal relative pose from lines at points (30x30)" takes 54.44 ms on GPU versus 234.86 ms for the best CPU solver (Julia, MINUS or the authors' own code, 8 cores). — [GPU-HC repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)

### Inferences
- The two papers report different averages for the same solver: 660 ms in MINUS itself and "about 500 ms" in Hruby et al. This is probably a different machine or code version. The README's 234.86 ms CPU figure suggests later CPU improvements, or that a faster tool than MINUS was used, since the README takes the best of Julia, MINUS and its own code.
- At 0.2–0.7 s per solve, all-path HC of degree around 300 is usable only for low-iteration RANSAC (a few hundred samples) or offline work.

### Gaps
- I did not read the TPAMI 2023 extended version. Updated timings, e.g. the "towards microseconds" variants, are not verified.
- The supplementary-material speed/reliability tradeoff numbers for the other formulations were not extracted.

---

## Q2. Hruby, Duff, Leykin, Pajdla, "Learning to Solve Hard Minimal Problems" (CVPR 2022)

### Takeaway
The method learns to select a start problem-solution pair ("anchor") with an MLP and then tracks one real HC path. On the "Scranton" relaxation of 4 points in 3 calibrated views (272 solutions) it solves 26.3% of inputs in 16.3 µs, which is 61.6 µs per correct solution in RANSAC. On the 5-point problem it solves 29.0% of inputs in 7.6 µs, or 26.1 µs per correct solution, against roughly 5 µs for Nistér's algebraic solver at about 100% success. The point is that work becomes independent of the number of spurious solutions; the cost is a lower per-sample success rate that RANSAC absorbs.

### Cited Findings
- **Abstract figures:**
  - Calibrated three-view relative pose "via a minimal relaxation using four points in each view. On average, we can solve a single problem in under 70 µs." — [arXiv 2112.03424](https://arxiv.org/abs/2112.03424)
  - The engineering choices are also studied on 5pt in two views. — [arXiv 2112.03424](https://arxiv.org/abs/2112.03424)
- **Scranton problem:**
  - A minimal relaxation of the overconstrained 4pt problem of Nistér–Schaffalitzky. Formulated in camera parameters, it has 272 complex solutions. The authors use a 3D-depth formulation that has an infinite spurious family (zero depths), but 272 potentially meaningful solutions remain. "The US telephone area code for Scranton PA is 272." — [Hruby et al.](https://arxiv.org/abs/2112.03424)
  - The original 4pt numerical solver ran in about 1 ms. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Pipeline:**
  - Offline: sample problem-solution (p-s) pairs from real SfM data (ETH3D models). Then pick anchors by greedy set cover on a graph whose edges mean "HC from (pi,si) to pj reaches sj". Then train a classifier σ(p) ∈ A ∪ {TRASH}.
  - Online: normalize the input, classify (about 8 µs), track one real path, and skip samples labelled TRASH.
  - Objective: minimize the effective time t = µt/ρ. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Model:**
  - MLP with 6 hidden layers of 100 neurons, PReLU-style activations, dropout before the last layer, trained for 80 epochs.
  - Classification takes "about 8µs for both 5pt and 4pt problems". — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Table 4** (anchors from Office+Terrains, evaluated on Delivery area+Facade; ρ = success %, µt = mean time, t = 100·µt/ρ):

  | Problem | Strategy | ρ | µt | t |
  |---|---|---|---|---|
  | 5pt | MLP with trash, A75 (best) | 29.0% | 7.6 µs | 26.1 µs |
  | 5pt | MLP without trash, A50 | 29.3% | 15.7 µs | 53.5 µs |
  | 4pt/Scranton | MLP with trash, A90 | 26.3% | 16.2 µs | 61.6 µs |
  | 4pt/Scranton | MLP with trash, A75 | 19.0% | 13.5 µs | 71.1 µs |
  | 5pt | Oracle "B1", track from all anchors in A50/A75/A90 | 47.3% / 74.2% / 87.7% | – | – |
  | 4pt | Oracle "B1", same | 44.2% / 72.0% / 87.9% | – | – |

  Nearest-anchor baselines (Euclidean or Mahalanobis distance) are much worse, about 5–14% success. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Headline claims:**
  - "It successfully solves about 26.3% of inputs in 16.3µs … about 4 samples suffice on average to obtain a valid candidate … in 61.6µs." This is "more than ten times speedup" over the 4pt solver of Nistér–Schaffalitzky (on the order of milliseconds). — [Hruby et al.](https://arxiv.org/abs/2112.03424)
  - "For the 5pt problem, we can solve 29.0% of inputs in about 7.6µs … in RANSAC … in average in 26.1µs." — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Comparison with algebraic and all-path HC solvers:**
  - "Off-the-shelf solvers are much slower (10^3 − 10^5 times) than current symbolic-numeric solvers." Solving 5pt with complex HC in Macaulay2 takes about 10^5 µs, compared with 5 µs for the Kukelova-type implementation of Nistér. The superscripts were lost in PDF extraction and are read here as 10^3–10^5 and 10^5 µs from context. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
  - "The most efficient symbolic-numeric solver of the 5pt problem solves for up to 10 essential matrices in about 5µs." — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Speedup sources against MINUS:**
  - The single path takes under 10 µs against about 500 ms for 312 paths, i.e. more than 1000× from tracking one path.
  - About 14× from real instead of complex arithmetic.
  - About 5× from optimized linear algebra, because sparse Jacobians admit closed-form solves.
  - The HC core builds on MINUS. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **RANSAC on the 2020 RANSAC Tutorial data (5pt):** the method needs about 4× more samples than Nistér to reach the same accuracy. That implies about 25% success, about 1.6× lower than on the ETH3D test data. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Limitations:**
  - It gives up the near-100% success of complex all-path HC.
  - The MLP "represents only what it is trained for", so it depends on the data distribution. — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Code:** https://github.com/petrhruby97/learning_minimal — [Hruby et al.](https://arxiv.org/abs/2112.03424)

### Inferences
- In effective time per correct solution, the learned 5pt solver (26 µs) is still about 5× slower than the algebraic 5pt (about 5 µs at about 100%). The approach pays off only for hard problems, such as Scranton with 272 solutions, where no fast algebraic solver exists.
- Success is bounded by anchor coverage: even the oracle reaches only about 88% with A90. The real-data success drop (to about 25%) shows sensitivity to distribution shift.

### Gaps
- I could not verify whether it was a best-paper candidate or award finalist; no source confirmed this.
- The CPU model used for the timings was not found in the main text (it may be in the supplementary material).

---

## Q3. Follow-ups 2022–2026: GPU-HC, Simulator HC, neural predictor-corrector, related tools

### Takeaway
GPU-HC (CVPR 2022; ISSAC 2022; IPDPS 2025) makes all-path HC run in milliseconds for problems with hundreds or thousands of solutions:
- 1784-solution trifocal pose with unknown focal length: 338 ms on an RTX 6000 in 2022, 16.7 ms in the current README.
- In 2022, GPU-HC beat elimination templates on everything except the smallest problems.

Simulator HC (Zhang, Dai, Xu, Kneip; arXiv 2024/25) replaces Hruby's anchor classifier with a regressor plus online simulation of a start pair. It reports higher success and about 4–17× speedups over Gröbner solvers on generalized-camera problems. An RL-learned step-size controller (NPC, arXiv 2026) cuts HC corrector iterations by about 45%.

### Cited Findings
- **GPU-HC (Chien, Fan, Abdelfattah, Tsigaridas, Tomov, Kimia, CVPR 2022):**
  - **Design:** one HC path per GPU warp, using kernel-fused batched linear solves (MAGMA-based). Start systems come from Macaulay2 monodromy. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
  - **Hardware:** 8-core 2.6 GHz Intel Xeon and an NVIDIA Quadro RTX 6000. Speedups range from 13× to 56× on polynomial benchmarks; for example, eco12 with 1024 solutions takes 227.54 ms on CPU and 17.06 ms on GPU. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
  - **Table 2** (each problem run 20 times; times in ms; ET = elimination template; "X" = ET impossible because it runs out of memory):

    | Problem | Unknowns | Solutions | ET | CPU-HC | GPU-HC | CPU/GPU | ET/GPU |
    |---|---|---|---|---|---|---|---|
    | Trifocal relative pose, unknown focal length | 18 | 1784 | X | 3326.07 | 338.27 | 10.21× | – |
    | 4-view triangulation | 14 | 296 | X | 156.28 | 18.60 | 8.32× | – |
    | 5pt relative pose with depth reconstruction | 16 | 160 | X | 150.94 | 26.89 | 5.61× | – |
    | 6pt rolling-shutter absolute pose, 1-linearized | 18 | 160 | X | 158.48 | 27.11 | 5.85× | – |
    | 3-view triangulation | 9 | 94 | 612.43 | 101.86 | 8.17 | 24.19× | 38.24× |
    | Optimal PnP (quaternion) | 4 | 128 | 36.33 | 80.26 | 7.18 | 11.18× | 5.06× |
    | P4P with unknown focal length and radial distortion | 5 | 192 | 9.03 | 130.79 | 7.51 | 17.42× | 1.2× |
    | 2-view triangulation with radial distortion | 5 | 28 | 5.92 | 66.22 | 3.06 | 21.64× | 1.93× |
    | Optimal P4P | 5 | 32 | 1.864 | 53.13 | 1.57 | 33.84× | 1.19× |
    | 3pt relative pose with homography constraint | 8 | 8 | 1.472 | 51.25 | 0.95 | 53.95× | 1.55× |
    | Relative pose with quiver, unknown focal length | 4 | 28 | 1.082 | 56.01 | 1.23 | 45.54× | 0.88× |
    | P3P | 3 | 8 | 0.063 | 39.64 | 0.22 | 180× | 0.29× |
    | 5pt relative pose without depths | 3 | 27 | 0.035 | 55.48 | 0.96 | 57.79× | 0.036× |

    — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
  - **Elimination template comparison:**
    - 3-view triangulation needs a 1866×1975 template for the full problem and 239×290 for the reduced one. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
    - A problem with 81 solutions needs a 575×656 template with an 81×81 action matrix. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
    - The authors say "it is not clear if [templates] can handle problems with ≥ 300 roots". — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
  - **Current README timings** (improved code; "1.2x-7.3x speedup over the old version" credited to the IPDPS 2025 paper; GPU model not stated; CPU = fastest of Julia, MINUS or the authors' own code on 8 cores):

    | Problem | GPU (ms) | CPU (ms) |
    |---|---|---|
    | Trifocal pose, unknown focal length | 16.67 | 346.29 |
    | Trifocal pose from lines at points (Chicago, 30×30) | 54.44 | 234.86 |
    | Generalized 3-view relative pose from 4 points | 4.77 | 428.26 |
    | Generalized 3-view relative pose from 6 lines | 5.11 | 1103 |
    | 4-view triangulation | 1.74 | 48.47 |
    | 3-view triangulation | 0.88 | 17.84 |
    | 5pt relative pose with depths | 7.26 | 27.46 |
    | P3P | 0.15 | 1.18 |

    — [GPU-HC repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)
  - **Limits:** system size at most 32×32; overdetermined systems not supported. — [GPU-HC repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)
  - **Related papers:**
    - Chien et al., "Parallel path tracking for homotopy continuation using GPU", ISSAC 2022. — [GPU-HC repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)
    - Chien, Abdelfattah, Kimia, "Accelerating Homotopy Continuation with GPUs: Application to Trifocal Pose Estimation", IPDPS 2025 (Milano), DOI 10.1109/IPDPS64566.2025.00110. Keywords include batch linear solver, RANSAC and SLAM. Seen in search snippets only. — [ICL UTK](https://icl.utk.edu/node/1886); [IEEE](https://ieeexplore.ieee.org/document/11078535/)
- **Simulator HC (Zhang, Dai, Xu, Kneip; arXiv 2411.03745, v1 Nov 2024, v2 Mar 2025; no venue listed):**
  - **Method:** a network regresses a rough solution. A matching start problem is then simulated online (the start pair is constructed, not retrieved), and one HC path is tracked to the target. One extra correspondence makes the regression overdetermined. — [Simulator HC](https://arxiv.org/html/2411.03745)
  - **Generalized P4P** (16 points, 2 px noise, 1000 trials):

    | Method | Success | Time |
    |---|---|---|
    | Gröbner basis | 100% | 1.254 ms |
    | CPU all-root HC | 100% | 9.796 ms |
    | GPU all-root HC | 100% | 1.002 ms |
    | Predicted start + Simulator HC | 100% | 0.316 ms |
    | Predicted start only | 32% | – |
    | Predicted start + LM refinement | about 54% | – |

    — [Simulator HC](https://arxiv.org/html/2411.03745)
  - **Generalized relative pose and scale** (140 solutions; Gröbner, Eigen-solver and all-path CPU/GPU HC baselines):
    - Success: 96.3% for 8-point Simulator HC (the text says about 95%) against 100% for the Gröbner basis solver.
    - Runtime: about 17× faster than Gröbner, about 5× faster than GPU-HC and about 40× faster than CPU-HC.
    - RANSAC on LaMAR (median success, time): CAB 68% in 2.08 s vs 60% in 2.83 s for Gröbner; HGE 84% vs 67%; LIN 86% vs 72%.
    - Hardware: 2.70 GHz Xeon and RTX 4090 (training).
    - Networks: 1D-CNNs trained on 2.4k (gP4P) and 64k (GRPS) samples. No code link found.
    - Hruby et al. is discussed but not benchmarked.

    — [Simulator HC](https://arxiv.org/html/2411.03745)
- **Neural Predictor-Corrector (Mai et al., arXiv 2602.03086, Feb 2026, cs.LG):**
  - **Method:** PPO reinforcement learning learns step sizes and corrector termination, applied to HC, GNC, Gaussian homotopy and annealed Langevin dynamics. It is trained on 4-view triangulation systems. — [NPC](https://arxiv.org/html/2602.03086v1)
  - **UPnP results:**
    - Classic HC: 100% success, 53 iterations, 8.25 per path.
    - NPC+HC: 100% success, 29 iterations, 3.86 per path. The time unit was not given in the extract.
    - Simulator HC needs 100 iterations. Its runtime is not comparable, because it is C++ while the others are Python.
    - NPC names Hruby et al. and Zhang et al. as works that "learn the starting system". — [NPC](https://arxiv.org/html/2602.03086v1)
- **Other related work** (seen only as search hits; content not verified):
  - "GPU-accelerated path tracker for polyhedral homotopy" (arXiv 2111.14317). — [arXiv](https://arxiv.org/pdf/2111.14317)
  - "Snapshot of Algebraic Vision" survey (arXiv 2210.11443) discusses offline/online HC, where the start system is effectively learned. — [arXiv](https://arxiv.org/pdf/2210.11443)
  - "Homotopy Iterators" (arXiv 2509.08084). — [arXiv](https://arxiv.org/pdf/2509.08084)
  - "HANN" homotopy auxiliary NN (arXiv 2509.26358). — [arXiv](https://arxiv.org/pdf/2509.26358)
- **Point-line minimal problem catalogues that create the hard HC targets:**
  - Duff, Kohn, Leykin, Pajdla, "PLMP" (ICCV 2019) and "PL1P" (ECCV 2020) list many problems with large degrees. Scranton comes from these. Hruby et al. also cite Duff, Korotynskiy, Pajdla, Regan, "Galois/monodromy groups for decomposing minimal problems" (arXiv 2105.04460). — [Hruby et al. references](https://arxiv.org/abs/2112.03424)

### Inferences
- Over time the field has moved from all-path complex HC on CPU (about 1 min with Bertini/M2, then 660 ms with MINUS), to all-path HC on GPU (tens of ms), to a single learned-start real path (about 10–300 µs). Each step trades guaranteed completeness for speed.
- Per the 2022 GPU-HC table, the crossover where GPU-HC beats CPU elimination templates sits at roughly 30 or more solutions and moderate template sizes. For P3P and 5pt (8 and 10 solutions), templates are 3–30× faster even than GPU-HC.

### Gaps
- I did not obtain IPDPS 2025 numbers beyond the README.
- I could not verify any Hruby-led 2023–2026 HC follow-up. Hruby's later works (e.g. point-line absolute pose, three-view focal length from homographies, arXiv 2501.07499) appear to be algebraic, not HC. This was not verified.
- I found no source for "PL1P" timings.
- I did not find "DeepMinimal" or a "Neural minimal solver" paper replacing algebraic solvers outright, and did not verify them. The closest are learned priors feeding classical solvers, e.g. Kocur et al. 2025, "Are minimal radial distortion solvers really necessary?" — [arXiv 2505.00866](https://arxiv.org/abs/2505.00866)

---

## Q4. Learned or optimized template choice and hybrid solver selection

### Takeaway
I found no peer-reviewed 2023–2026 work that uses machine learning to choose or optimize Gröbner or resultant elimination templates. Template optimization has used heuristic or greedy search, e.g. Martyushev–Vráblíková–Pajdla CVPR 2022 and IJCV 2023/24. Learning in this area has gone into start-system selection for HC (Hruby 2022, Zhang 2024/25) and step control (NPC 2026).

### Cited Findings
- **Martyushev, Vráblíková, Pajdla, "Optimizing Elimination Templates by Greedy Parameter Search" (CVPR 2022):**
  - Builds an affine parameterization of the template family, then searches it greedily.
  - Tested on 34 minimal problems. For refractive absolute pose with unknown focal length it found a template 20× smaller.
  - Code: github.com/martyushev/eliminationTemplates (MATLAB/Julia/Python). — [arXiv 2203.14901](https://arxiv.org/pdf/2203.14901); [GitHub](https://github.com/martyushev/eliminationtemplates)
- **Follow-up:** "Automatic Solver Generator for Systems of Laurent Polynomial Equations" (arXiv 2307.00320; IJCV) generates state-of-the-art templates for many minimal problems. — [arXiv 2307.00320](https://arxiv.org/pdf/2307.00320)
- **Macaulay2 package:** EliminationTemplates exists ("zero-dimensional polynomial solvers based on linear algebra"). — [M2 docs](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EliminationTemplates/html/index.html)
- **Hruby et al. on symbolic-numeric solvers:** "work well for problems with as many as 64 complex solutions" but "mostly fail to deliver practical results for hard computer vision problems with many solutions". — [Hruby et al.](https://arxiv.org/abs/2112.03424)
- **Unverified preprint:** an SSRN preprint (multi-camera minimal relative pose) claims a GNN that predicts the solver configuration among 5,953 match types, reducing compute by 73%. It is not peer-reviewed and was not verified. — [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7518220)

### Inferences
- Learned template selection appears to be an open niche. The parameterized template families of Martyushev et al. form a natural search space for learned or RL optimization, but no such paper was found.

### Gaps
- No sources were found for transformer- or symbolic-regression-based minimal solvers, or for "learned RANSAC solver selection" in a peer-reviewed venue.

---

## Q5. When HC beats templates; code ecosystem

### Takeaway
- Templates win for low-degree problems (up to a few tens of solutions): P3P takes 0.063 ms, 5pt 0.035 ms, and algebraic 5pt about 5 µs.
- All-path HC (GPU) wins for roughly 90 or more solutions, or wherever template construction runs out of memory: triangulation with 94–296 solutions, and problems with 160 or 1784 solutions.
- Learned single-path HC is the only µs-scale option for very hard problems such as Scranton (272), at about 25–30% per-sample success.

### Cited Findings
- **Templates beat GPU-HC on small problems:** ET/GPU ratio is 0.036× for 5pt without depths, 0.29× for P3P and 0.88× for the quiver problem. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
- **GPU-HC beats templates on larger problems:**
  - 3-view triangulation (94 solutions): 38× faster than the ET (612 ms against 8.17 ms).
  - Optimal PnP (128 solutions): 5× faster.
  - Problems with 160–1784 solutions: the ET is out of memory. — [GPU-HC CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)
- **CPU all-path HC versus a Gröbner solver:** on gP4P, CPU all-path HC is about 8× slower (9.8 ms against 1.25 ms). Learned-start single-path HC is 4× faster than Gröbner (0.316 ms). — [Simulator HC](https://arxiv.org/html/2411.03745)
- **Success rates:**
  - All-path complex HC: about 99% (MINUS fails about 1%) or 100% (GPU-HC and Simulator HC baselines).
  - Learned single path: 26–29% (Hruby) or 95–100% (Simulator HC).
  — [MINUS](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf); [Hruby](https://arxiv.org/abs/2112.03424); [Simulator HC](https://arxiv.org/html/2411.03745)
- **Code and tools:**
  - MINUS: github.com/rfabbri/minus — [MINUS](https://openaccess.thecvf.com/content_CVPR_2020/papers/Fabbri_TRPLP_-_Trifocal_Relative_Pose_From_Lines_at_Points_CVPR_2020_paper.pdf)
  - Hruby: github.com/petrhruby97/learning_minimal — [Hruby](https://arxiv.org/abs/2112.03424)
  - GPU-HC: github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU (Brown-LEMS) — [repo](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU)
  - The off-the-shelf HC tools cited by Hruby et al.:
    - HomotopyContinuation.jl (Breiding & Timme, ICMS 2018, LNCS 10931:458–465).
    - Bertini (Bates, Sommese, Hauenstein, Wampler, SIAM 2013).
    - Monodromy solving (Duff, Hill, Jensen, Lee, Leykin, Sommars, IMA J. Numer. Anal. 2018; the Macaulay2 MonodromySolver).

    — [Hruby references](https://arxiv.org/abs/2112.03424)
  - GPU-HC builds its start systems with the Macaulay2 monodromy module. — [GPU-HC](https://openaccess.thecvf.com/content/CVPR2022/papers/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.pdf)

### Inferences
- For a sparse-resultant or template generator project, HC is the natural fallback or baseline in two cases: problems above roughly 100–300 solutions, and problems where template construction runs out of memory. For problems below about 50 solutions, optimized templates remain 1–2 orders of magnitude faster than any HC variant.

### Gaps
- The HomotopyContinuation.jl documentation (monodromy_solve, parameter homotopy API) was not fetched directly. The URL is likely homotopycontinuation.org but was not verified here.
- I did not measure or find a head-to-head comparison of HomotopyContinuation.jl timings on 5pt or P3P.
