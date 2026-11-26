# Exact eigen size and feasibility of sparse-resultant templates (proof on paper)

This proves the conjecture from `face_obstruction.md`, under explicit hypotheses:

> **minimal eigen size = number of solutions of F on the toric variety X_Q
> (torus solutions + boundary solutions, with multiplicity)**, and a template
> with hidden x_h exists **iff x_h has no pole on that solution set**.

It covers non-generic coefficients, because nothing below assumes genericity.
The proof uses standard toric geometry and Serre vanishing. It is **not**
formalised in Lean: Mathlib currently has neither toric varieties nor coherent
sheaf cohomology. Theorem 1 of `face_obstruction.md`, which is formalised, gives
the "only if" direction.

References: [CLS] Cox, Little, Schenck, *Toric Varieties* (2011);
[H] Hartshorne, *Algebraic Geometry*, III.5.

## Setting and hypotheses

- M = Z^n, K an algebraically closed field. F = {f_1, …, f_m} ⊂ K[M] with
  P_i = NP(f_i).
- Q ⊂ M_R is a full-dimensional lattice polytope, Σ = Σ_Q its normal fan, and
  X = X_Q the projective toric variety. D = Σ_ρ a_ρ D_ρ is the ample divisor of Q,
  where Q = {m : ⟨m, u_ρ⟩ ≥ −a_ρ}, and O(k) := O_X(kD).
- **(A) Compatibility.** Σ refines the normal fan of every P_i. This holds for every
  Minkowski sum Q = Σ_{i∈S} P_i + (other polytopes) with S ⊇ {all i used in rows}.
  Then each P_i is the polytope of a basepoint-free Cartier divisor D_i, and
  f_i ∈ H^0(O(D_i)) [CLS Prop. 4.3.3, Thm. 6.1.7].
- **Ideal sheaf.** 𝓘 = image of ⊕_i O(−D_i) → O_X, (s_i) ↦ Σ f_i s_i. Its zero
  scheme is Z ⊂ X.
- **(B) Finiteness.** Z is finite, i.e. zero-dimensional. Each point z ∈ Z lies on a
  unique torus orbit O_τ, for a face τ of Q (τ = Q means the torus).
- **Template.** B_k = kQ ∩ M, rows x^t f_i for all t with t + supp f_i ⊆ kQ, matrix
  M_k. The eigen basis is E_k = {b ∈ B_k : b + e_h ∈ B_k}, and the boundary is
  ∂_k = B_k \ E_k.
- **corank_k** = |B_k| − rank M_k. **Feasible** means M_k[:, ∂_k] is injective
  (symgen's criterion before any reduction).

## Step 0: the rows are global sections

For a torus-invariant divisor D′ = Σ a′_ρ D_ρ,

    H^0(O(D′)) = span{ x^m : ⟨m, u_ρ⟩ ≥ −a′_ρ for all ρ }    [CLS Prop. 4.3.3].

Take D′ = kD − D_i, so a′_ρ = k a_ρ − a_{i,ρ} with a_{i,ρ} = −min_{p∈P_i}⟨p, u_ρ⟩.
The condition then says exactly t + P_i ⊆ kQ, which by convexity is t + supp f_i ⊆ kQ.
So the row multipliers of f_i are a basis of H^0(O(kD − D_i)), and

    rowspace(M_k) = R_k := Σ_i f_i · H^0(O(kD − D_i)) ⊆ V_k := H^0(O(kD)) = span{x^b : b ∈ B_k}.

## Step 1: corank_k = length(Z) for k ≫ 0

Tensor the exact sequence 0 → 𝒦 → ⊕ O(−D_i) → 𝓘 → 0 with the line bundle O(kD),
and also use 0 → 𝓘(k) → O(k) → O_Z(k) → 0. By Serre vanishing [H, III.5.2], since D
is ample:

- H^1(𝒦(k)) = 0 for k ≫ 0, so H^0(⊕ O(kD − D_i)) → H^0(𝓘(k)) is onto. Hence
  **R_k = H^0(𝓘(k))**.
- H^1(𝓘(k)) = 0 for k ≫ 0, so V_k → H^0(O_Z(k)) is onto, with kernel H^0(𝓘(k)).

Z is finite, so O_Z(k) ≅ O_Z and dim H^0(O_Z(k)) = length(Z). Therefore

    corank_k = dim V_k / R_k = length(Z) = Σ_{z ∈ Z} mult_z.     ∎

The torus points of Z are the solutions of F in (K*)^n. A point on O_τ, τ ≠ Q,
corresponds to a torus orbit of solutions of the initial system in_w(F) for w in the
relative interior of the (outer) normal cone of τ. This is exactly what `face_test.py`
counts, for facets.

## Step 2: feasibility ⇔ x_h regular on Z (k ≫ 0)

Kernel vectors φ ∈ K^{B_k} of M_k are exactly the functionals on V_k that vanish on
R_k. φ is supported on ∂_k iff it also vanishes on W_k := span{x^b : b ∈ E_k}. So

    feasible  ⇔  W_k + R_k = V_k  ⇔  W_k → H^0(O_Z(k)) is onto   (k ≫ 0, by Step 1).

**(⇐)** Assume x_h is regular at every z ∈ Z. On the normal variety X this means
⟨e_h, u_ρ⟩ ≥ 0 for every ray ρ with z ∈ D_ρ. These are the rays of the cone σ_τ of
the orbit containing z (algebraic Hartogs; div(x^m) = Σ⟨m, u_ρ⟩D_ρ).

For each such face τ, choose c ≫ 0 and a lattice point q_τ deep in relint(cτ). Then:

- For ρ ∈ σ_τ(1): ⟨q_τ + e_h, u_ρ⟩ = −c a_ρ + ⟨e_h, u_ρ⟩ ≥ −c a_ρ.
- For the other rays, q_τ is at distance ≥ |⟨e_h, u_ρ⟩| from the facet.

So q_τ, q_τ + e_h ∈ cQ, and for every m ∈ (k−c)Q both q_τ + m and q_τ + m + e_h lie in
kQ. Hence

    W_k ⊇ Σ_τ x^{q_τ} · V_{k−c}.

The section s_τ = x^{q_τ} ∈ H^0(O(cD)) does not vanish on O_τ, because q_τ lies on
the face cτ [CLS Prop. 4.3.? / 6.1]. So the s_τ have no common zero on the finite
scheme Z. H^0(O_Z(k)) is a free rank-1 module over the Artinian ring
A = H^0(O_Z) = ∏_z A_z, and the images of the s_τ generate the unit ideal (each A_z
contains some s_τ as a unit). Since V_{k−c} → H^0(O_Z(k−c)) is onto for k−c ≫ 0,
W_k → H^0(O_Z(k)) is onto. So the template is feasible.

**(⇒)** Suppose x_h has a pole or an indeterminacy at some z ∈ O_τ. Then
⟨e_h, u_ρ⟩ < 0 for some ρ ∈ σ_τ(1). Pick an integral u in the relative interior of σ_τ
with ⟨e_h, u⟩ < 0 (an open, non-empty condition), and set w = −u, so w_h > 0. The point
z lifts to a torus solution of in_w(F). By **Theorem 1 (Lean: `no_valid_template`)**,
no template is feasible, for any monomial set at all.     ∎

## Step 3: the minimal eigen size

If the template is feasible, the columns ∂_k are linearly independent. By the
exchange property of the column matroid of M_k they extend, using columns of E_k
only, to a basis of the column space. symgen's matroid step does exactly this. The
resulting eigen basis has size |B_k| − rank M_k = corank_k. Any eigen basis of any
valid template on B_k has at least this size, since kernel dimension ≤ |E|. Combined
with Steps 1–2:

> **Theorem 3.** Under (A) and (B), for all k ≫ 0: a hidden-x_h template on kQ exists
> iff x_h is regular on Z. In that case the minimal eigen size equals
> length(Z) = #torus solutions + #boundary solutions (with multiplicity). The
> boundary solutions give the spurious eigenvalues x_h(z), which are finite, and 0 on
> orbits where x_h vanishes (cf. Theorem 2).

The eigenvalues are the values x_h(z): the solver matrix is the transpose of
multiplication by x_h on A = H^0(O_Z), via ker M_k ≅ A*.

## Checks against computation

| system, Q | prediction (length Z; x_h regular?) | observed |
|---|---|---|
| Ex. 4.2.1 instance, Q = NP(f1)+NP(f2), k = 1…6 (`k_sweep.py`) | 8 + 1 = 9; x1 regular, x2 has a pole | corank 9 for every k; x1 feasible with eigen size 9; x2 infeasible for every k |
| 6pt focal, Minkowski sums of 30 to 517 points | 15 + 3 = 18; a1, a2 regular, a3 has a pole | 18 for every basis; a3 has no template |
| 9pt two radial | 24 + 3 = 27; only a3 regular | 27; a1, a2, a4 have no template (0/234 bases each) |
| 6pt focal, saturated | 15 + 0 | 15 |

## Limits

- **"k ≫ 0" is not effective.** The threshold is a Castelnuovo–Mumford-type
  regularity bound. In every example above the count was already exact at k = 1,
  but the theorem does not guarantee that. For small k, Theorem 2 still gives the
  lower bound, and Step 0 still holds.
- **(A) matters.** For a convex B that is not compatible with the supports (a random
  convex lattice set), Step 0 fails, the row space is smaller, and the coranks are
  larger and shape-dependent. This is consistent with the 12.7 % feasibility of
  random convex sets against the Minkowski sum in `exp_shapes.py`. The theorem
  therefore replaces the thesis's "convexity" by **normal-fan compatibility**.
- **(B) can fail.** Z can have positive-dimensional components at infinity. Then the
  corank grows with k, and templates get larger without bound.
- **Variable recovery** (pairs x^e, x^{e+e_i}) is not covered. It is a separate
  combinatorial condition.
- **Lean status:** see `research/lean/README.md`. Theorems 1–2, Step 3 and the
  assembly `Theorem3.toric_eigen_size` are machine-checked. The Serre-vanishing step
  (R3) is an explicit hypothesis, because Mathlib has no toric varieties or coherent
  cohomology. The independence step (R1/R2) is still a hypothesis and is the next
  item to prove.
