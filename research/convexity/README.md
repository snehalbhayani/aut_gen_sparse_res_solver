# Why Newton polytopes work: notes towards a paper

This folder studies the intuition in Section 4.2.1 of the thesis (*Sparse
resultant-based methods with their applications to generalized cameras*, Acta
Univ. Oul. C 903, 2023). That section asks why taking the monomial set B as the
lattice points of a convex set (the Minkowski sum of Newton polytopes) makes a
feasible resultant matrix likely. Page numbers below are the thesis's printed
page numbers.

Scripts (all exact, over GF(p), using `symgen`):

| file | what it does |
|---|---|
| `ex421.py` | reproduces the numbers of Ex. 4.2.1 (F+λ 8pt instance) |
| `exp_shapes.py` | feasibility and eigen size for Minkowski, random convex and non-convex sets of equal size |
| `face_test.py` | Bernstein face test: counts torus solutions of the initial systems in_w(F) on each facet of Q |

## 1. Open questions in the thesis (complete list)

| # | page | question | status / what this folder adds |
|---|---|---|---|
| T1 | Abstract, 49, 91 | "For a given polynomial system, what is the smallest/fastest solver?"; search the spectrum from hidden-variable GEP to extra-polynomial solvers | `symgen` measured part of the extra-polynomial side: the minimal eigen size was fastest there. The GEP end is untested. |
| T2 | 21 (contrib. 5), 58–59, 91 | Hypothesis: convexity of the Minkowski sum is crucial for feasibility. "No proof of this observation is presented." Prioritise monomial sets that are interiors of convex sets. | §2–3 below: partly supported, and given a precise mechanism |
| T3 | 34 | Can the sparse mixed resultant be written as det / minor (Macaulay style)? Open in mathematics. | not addressed |
| T4 | 45 | Optimal resultant matrix construction for non-generic systems is open (Cox et al., Ch. 7) | §3: the non-genericity that matters can be located face by face |
| T5 | 48–49 | Hidden-variable solvers have spurious roots from the two linearisations and from det being a non-trivial multiple of the resultant. Deflation "can not be performed symbolically (in the offline stage)." | §3: the extra roots are predicted offline by the face test, and removable offline by saturation (done for 6pt focal) |
| T6 | 49 | Stability depends on the sparsity and type of entries of the inverted matrix A12; prefer stable A12 | open |
| T7 | 62, 65 | A unified AM/SRes generator; border bases; "AM guarantees the smallest eigenvalue problem without iterative rank tests" | `symgen`'s matroid step gives the minimal eigen size for a fixed row set in one rank computation |
| T8 | 28 | Coefficient computation time is "impossible to be quantified", so speed is compared by matrix sizes | `symgen` measurements show coefficient evaluation can dominate run time; matrix sizes alone mislead |
| T9 | 70–71 | Partially calibrated semi-generalized pose (236 roots) and the calibrated case (64 roots) were infeasible for template generators | open (homotopy-continuation territory) |
| T10 | 78–79 | The degeneracy g33 ≈ 0 cannot be avoided in the partially calibrated homography case | open |
| T11 | 83, 92 | Hybrid RANSAC choosing the solver per iteration | open |

**Erratum, p. 60 (Ex. 4.2.1):** "6 shifts of F1, and 20 shifts of F1" should read
"6 shifts of F1 and 16 shifts of F2". Here NP(f2) = [0,1]×[0,2] and Q = [0,4]×[0,5],
so F2 fits 4·4 = 16 times (`ex421.py`).

## 2. What the Ex. 4.2.1 argument actually measures

The thesis compares a convex B (30 points) with a non-convex B′ (30 points):

- **C1 (rows):** 6 + 16 translates of the supports fit in B, against 2 + 13 in B′.
  This is a real difference.
- **C2 (choices of B1):** 1,081,575 against 735,471 subsets of size 8 whose x2-shift
  stays in B. These are exactly C(25,8) and C(24,8). The counts only say that 25
  against 24 points of B have their x2-shift inside B. The large numbers hide a
  one-point difference.

C2 is also not really about the number of subsets. symgen's exact formulation
says a template is valid iff `rank Top[:, G] = |G|`, where G must contain every
point of B whose shift leaves B. These are the points on the face of B
*exposed by the shift direction Δ = e_h*. What decides C2 is whether those
boundary columns can be eliminated, not how many subsets exist.

### Experiment: convex vs non-convex sets (`exp_shapes.py`, same instance, |B| = 30)

| x1 hidden | feasible | eigen sizes |
|---|---|---|
| Minkowski sum NP(f1)+NP(f2) | yes | 9 |
| 300 random convex lattice sets | 12.7 % | 9–13 |
| 2 random boundary swaps (convexity defect ≈ 1.8) | 11.0 % | 9–11 |
| 5 swaps (defect ≈ 4.7) | 1.3 % | 9–12 |
| 10 swaps (defect ≈ 8.2) | 0 % | — |
| **x2 hidden: every shape, including the Minkowski sum** | **0 %** | — |

Two lessons:

1. **Convexity helps, but shape matters more.** A random convex set is barely more
   likely to work than a slightly non-convex one (12.7 % against 11.0 %), and far less
   than the Minkowski sum. The useful property seems to be that Q's facets are
   aligned with the facets of the supports' Newton polytopes (compatible normal
   fans), not convexity alone.
2. **For x2 hidden, nothing works, whatever the shape.** That points to a cause that
   has nothing to do with convexity.

## 3. The mechanism: solutions at toric infinity (face test)

Take the x2-hidden case. The columns that must be eliminated are the top layer
of B (x2 = 5). The only rows that touch it are the leading forms of f1 and f2 in
x2:

    g1(x1) = −45x1³ + 375x1² − 660x1 − 180,   g2(x1) = 3x1 − 9,   g1(3) = 0.

They share the root x1 = 3, so that block is a singular Sylvester matrix. In other
words, the face of Q exposed by w = (0,1) carries a root: a solution at x2 = ∞.
This is also why BKK gives 9 roots but the instance has 8.

**General statement (conjecture, consistent with every test so far).**

Let Q be the lattice polytope that gives B, and let w be a facet normal of Q. Let
in_w(F) be the initial system (the terms of each f_i on its face maximising w·α).

**(a) Excess eigenvalues.** Each torus solution of in_w(F) is a solution on the
boundary of the toric variety X_Q. Such points stay in the kernel of the template
and show up as extra eigenvalues:

    minimal eigen size = #torus solutions + Σ over facets w of #torus solutions of in_w(F).

**(b) Hidden variable rule.** Near such a point x ≈ t^w·ξ with t → ∞.

- **w_h > 0:** x_h → ∞. No template with the partition B21 = I (eigenvalue u0) can exist.
- **w_h = 0:** x_h stays finite, giving a finite spurious eigenvalue.
- **w_h < 0:** x_h → 0, giving a spurious zero eigenvalue.

The alternative partition A21 = I (eigenvalue 1/u0, thesis p. 47) swaps the roles
of ∞ and 0. It should rescue the w_h > 0 case. symgen does not implement it yet.

**(c) Saturation.** Syzygies of the leading forms in_w(f_i) give new equations of
lower w-degree. These vanish on the true solutions but not on the face solutions,
so they remove the excess offline. (Thesis p. 48 says deflation cannot be done
symbolically offline; this is how it can be.)

### Validation (`face_test.py`; identical counts on three independent seed/prime instances)

| problem | facet w carrying torus solutions | predicted eigen size | observed eigen size |
|---|---|---|---|
| Ex. 4.2.1 instance (F+λ 8pt) | (0,1): 1 | 8 + 1 = 9 | 9 (x1 hidden) |
| 6pt focal relative pose | (0,0,1): 3 | 15 + 3 = 18 | 18 (Maple and symgen) |
| 9pt two radial distortions | (1,1,0,1): 3 | 24 + 3 = 27 | 27 (Maple and symgen) |
| optimal PnP (Hesch) | none | 27 | 27 |
| 7pt E+fλ (partially eliminated) | none | 19 | 19 |
| 6pt focal, saturated (`problem_relpose_6p_focal_sat`) | none | 15 | 15 |

Hidden-variable predictions:

| problem | prediction | observed |
|---|---|---|
| 6pt focal | a3 infeasible (w3 = 1); a1, a2 feasible | a3 gives no template; a1, a2 work |
| 9pt two radial | a1, a2, a4 infeasible (w = (1,1,0,1)); only a3 feasible | 0/234 feasible for a1, a2 and a4; 74/163 feasible for a3, 68 of them with eigen size 27 |
| F+λ 8pt | x2 infeasible | 0 % for every shape |

The 9pt result also explains why the repository configuration uses
`hiddenVarNum = 3`: it is the only variable that works.

## 4. A paper this could become

**Working title:** *Why Newton polytopes work: toric infinity, spurious
eigenvalues and hidden-variable feasibility in sparse resultant minimal solvers.*

**Claims, each backed by an exact, reproducible test:**

1. Template validity is exactly `rank Top[:, G] = |G|`, and the minimal eigen size
   for a row set is its corank (matroid argument). Both are certified over GF(p).
2. **Face theorem.** The excess of the eigen size over the true root count equals the
   number of solutions of the facet initial systems on X_Q. This answers T5: the
   spurious roots are predicted offline.
3. **Hidden-variable rule.** The signs of w_h decide which variable can be hidden and
   which partition works. This predicts feasibility before any template is built.
4. **Shape.** Rows (C1) come from erosion counts |B ⊖ A_i|. For the Minkowski sum
   these are at least the lattice points of the sum of the other polytopes (a
   Brunn–Minkowski / isoperimetric argument). Eliminability of the exposed boundary
   (C2) is decided by the face systems, not by subset counts. Together these refine
   the convexity hypothesis (T2) into "normal-fan compatible with the supports, and
   no solutions on the Δ-exposed face".
5. **Algorithm.** A face test, then automatic saturation along the bad facets, then
   template search. Payoff: 9pt radial from eigen size 27 to 24 (expected), 6pt focal
   from 18 to 15 (shown by hand).

**What still has to be done before submission:**

- Proofs:
  - (2) via the Hilbert function of the ideal in the Cox/toric ring of Q, for large
    dilations kQ.
  - (3) via the Puiseux / tropical limit of solutions.
  - Precise conditions where B is not a dilation of a compatible Q.
- Implement the alternative partition (eigenvalue 1/u0) in symgen, and check that
  it makes x2 (F+λ 8pt) and a3 (6pt focal) hideable.
- Automate saturation along a general facet w. Do 9pt radial (w = (1,1,0,1)) first:
  compute syzygies of in_w(f_i), lift them symbolically, and verify eigen size 24.
- Run the face test and the shape experiment across the whole benchmark (about 30
  problems from Martyushev et al. 2022 and the JMIV 2024 paper). Report the
  prediction accuracy of (2) and (3), and how many problems saturation improves.
- Position the work against Bender & Telen, *Toric eigenvalue methods for solving
  sparse polynomial systems* (Math. Comp. 2022), and Telen's truncated normal forms.
  They already use toric compactifications for eigenvalue methods, so the novelty has
  to be the computer-vision setting: non-generic, overdetermined systems, the
  explanation of the spurious eigenvalues and hidden-variable choices, and
  automatic saturation. Also cite Larsson et al. (saturated ideals, ICCV 2017) and
  Martyushev, Bhayani and Pajdla (Laurent generator, 2023).

**Venue timing (check the official dates):** the CVPR deadline is usually in
mid-November, so CVPR 2027 has probably passed. The ICCV 2027 deadline is
typically in early March 2027, which fits the work above.
