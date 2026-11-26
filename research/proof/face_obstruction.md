# Face obstruction theorem for sparse-resultant (extra polynomial) templates

The two theorems below hold over **any field K**, for **arbitrary coefficients**
(non-generic systems included) and for **any finite monomial set**, convex or not.
They prove the "necessary" half of the face test in `research/convexity`.
The key lemma and both theorems are **machine-checked in Lean 4 + Mathlib**:
`research/lean/FaceObstruction.lean` (`face_kernel`, `no_valid_template`,
`spurious_eigenpair`; no `sorry`). The corollary, the alternative-partition
remark and the orbit counting are not formalised.

## Setting

Exponents live in Z^n. For ξ ∈ (K\*)^n write ξ^b = ∏_j ξ_j^{b_j}.

**Equations.** f_i = Σ_{α∈S_i} c_{iα} x^α, where each S_i ⊂ Z^n is finite.

**Template.**
- A finite column set C ⊂ Z^n.
- Rows r = (t_r, i_r), each meaning the polynomial x^{t_r} f_{i_r}. The row's monomials
  must lie in C: t_r + S_{i_r} ⊆ C.
- The template matrix M ∈ K^{R×C} has entries M_{r,b} = c_{i_r, b−t_r} (zero if
  b − t_r ∉ S_{i_r}).

**Eigen basis and elimination block.**
- E ⊆ C with E + e_h ⊆ C (the l-rows e·(x_h − λ)).
- G = C \ E, and A2 = M[:, G].

**Valid template.** M restricted to the columns G is injective, i.e. A2 has full
column rank (symgen's criterion).

**Solver matrix.** Any linear map X on K^E with X(v|_E) = (v_{e+e_h})_{e∈E} for every
v ∈ ker M. This is the matrix the generated solver diagonalises.

**Weights and initial forms.** For a weight w ∈ Z^n, write ⟨w,α⟩ for the w-degree.
For each k ∈ Z let

    φ_{i,k}(ξ) = Σ_{α∈S_i, ⟨w,α⟩=k} c_{iα} ξ^α.

The initial form in_w(f_i) is the top component φ_{i, top_i} with
top_i = max_{α∈S_i} ⟨w,α⟩.

**Face solution.** A point ξ ∈ (K\*)^n with in_w(f_i)(ξ) = 0 for all i. Equivalently,
φ_{i,k}(ξ) = 0 for every k ≥ top_i, since for k > top_i the sum is empty.

## Key lemma (kernel vector at infinity)

Let C ≠ ∅, D = max_{b∈C} ⟨w,b⟩, and L = {b ∈ C : ⟨w,b⟩ = D} (the w-top layer of the
columns). If ξ is a face solution for w, then

    ṽ_b = ξ^b if b ∈ L,   ṽ_b = 0 otherwise

satisfies ṽ ≠ 0 and Mṽ = 0.

*Proof.* L ≠ ∅, and ξ^b ≠ 0 because ξ lies in the torus, so ṽ ≠ 0.

For a row r = (t, i), the entries M_{r,b} are non-zero only for b = t + α with α ∈ S_i:

    (Mṽ)_r = Σ_{α∈S_i, t+α∈L} c_{iα} ξ^{t+α}
           = ξ^t Σ_{α∈S_i, ⟨w,α⟩ = D−⟨w,t⟩} c_{iα} ξ^α
           = ξ^t φ_{i, D−⟨w,t⟩}(ξ).

Because t + S_i ⊆ C, we have ⟨w,t⟩ + top_i ≤ D, so D − ⟨w,t⟩ ≥ top_i. Hence
φ_{i, D−⟨w,t⟩}(ξ) = 0, and (Mṽ)_r = 0. ∎

## Theorem 1 (obstruction): w_h > 0 ⇒ no valid template

If some w with w_h > 0 has a face solution, then **no** template with hidden x_h and
eigenvalue λ = x_h is valid. This holds for every C, E and every choice of rows.

*Proof.* Suppose b ∈ L ∩ E. Then b + e_h ∈ C and ⟨w, b+e_h⟩ = D + w_h > D,
contradicting the maximality of D. So L ⊆ G. Then ṽ is non-zero, supported in G, and
lies in ker M, so M[:, G] is not injective. ∎

This explains every observed failure:
- **6pt focal, a3 hidden:** w = (0,0,1).
- **9pt two radial, a1, a2 or a4 hidden:** w = (1,1,0,1). Measured: 0 of 234 candidate
  bases feasible for each.
- **F+λ 8pt (thesis Ex. 4.2.1), x2 hidden:** w = (0,1). Measured: 0 % feasible for every
  shape.

Convexity of the monomial set cannot help, because the theorem holds for every C.

**Alternative partition** (eigenvalue 1/λ, rows e·(1 − λ⁻¹x_h), E − e_h ⊆ C). The
same proof with −e_h in place of e_h shows that face solutions with **w_h < 0** obstruct
it. So a face solution with w_h > 0 rules out the λ partition but not the 1/λ one.

## Theorem 2 (spurious eigenpairs): w_h ≤ 0 ⇒ extra eigenvalues

Let the template be valid, and let ξ be a face solution for w with w_h ≤ 0. Then:

1. ṽ|_E ≠ 0.
2. X(ṽ|_E) = μ·ṽ|_E for every solver matrix X, where μ = ξ_h if w_h = 0 and μ = 0 if
   w_h < 0.

*Proof.*

(1) Otherwise ṽ is supported in G, lies in ker M and is non-zero, which contradicts
validity.

(2) Take e ∈ E. Then e ∈ C and e + e_h ∈ C, and (X ṽ|_E)_e = ṽ_{e+e_h}.

- **w_h = 0:** ⟨w, e+e_h⟩ = ⟨w, e⟩, so e ∈ L ⇔ e + e_h ∈ L.
  - If e ∈ L: ṽ_{e+e_h} = ξ^{e+e_h} = ξ_h ξ^e = ξ_h ṽ_e.
  - Otherwise both entries are 0.
- **w_h < 0:** ⟨w, e+e_h⟩ < ⟨w, e⟩ ≤ D, so e + e_h ∉ L and ṽ_{e+e_h} = 0 = 0·ṽ_e. ∎

**w = 0** is the classical correctness statement. Then in_0(f_i) = f_i, L = C and ṽ is
the monomial vector of a true torus solution, with eigenvalue ξ_h.

**Corollary (lower bound on the eigenvalue problem).** Eigenvectors with distinct
eigenvalues are linearly independent. Hence for every valid template

    |E| ≥ #{ distinct values ξ_h over the true torus solutions
             and the face solutions with w_h = 0 },

plus one if some face solution with w_h < 0 exists and 0 is not already counted.

Points on one orbit t^w·ξ give the same vector up to scale, and the same ξ_h when
w_h = 0, so they count once. `face_test.py` counts exactly these orbits.

**Numerical confirmation** (`research/convexity/check_spurious_6p.py`). For 6pt focal
with a2 hidden, the three spurious eigenvalues of the generated solver equal the
a2-values of the three solutions of {F33(a1,a2) = 0, det F(a1,a2) = 0}. This is the face
w = (0,0,1), with w_2 = 0. The maximum relative mismatch over 20 random instances is
1.1·10⁻¹².

## What is *not* proved

**The equality |E_min| = #true + #face solutions.** It was observed on all four repo
problems and the thesis example, and it would need:

- B large enough and compatible with the supports (normal fans);
- a Hilbert-function argument on the toric variety X_Q;
- care with multiplicities and with solutions on lower-dimensional faces.

It is false for small or badly shaped B.

**Sufficiency.** The theorems do not show that a template exists once no face
obstruction is present.

**The convexity hypothesis itself** (thesis, Section 4.2.1). The theorems show that the
mechanism behind failures and spurious eigenvalues is independent of convexity.
Whether convex shapes maximise the row count is a separate, combinatorial
(Brunn–Minkowski-type) question.
