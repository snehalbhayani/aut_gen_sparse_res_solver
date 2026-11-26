/-
Theorem 3 of `research/proof/exact_count.md` (toric eigen size), for solution points
of multiplicity 1, assembled from the machine-checked pieces:

* `FaceObstruction.face_kernel`: every solution point z of F on the toric variety
  (a torus point, w = 0, or a boundary point, a face solution of in_w F) gives a
  kernel vector `φ z` of the template;
* `EigenSize`: for a feasible template the minimal eigen size is the corank.

Remaining inputs, stated as explicit hypotheses:

* `H_span` (**R3**, the Serre-vanishing step): for the large dilation kQ, every kernel
  vector is a combination of the point vectors `φ z`.
* `H_sep` (**R1/R2**): the point vectors restricted to the shiftable columns `E` are
  linearly independent. For kQ with k large and x_h regular on Z this is elementary;
  discharging it is the next milestone.

Conclusions: corank = |Z|, the template is feasible (the boundary columns are
eliminable), and the minimal eigen size is exactly |Z|. Together with
`FaceObstruction.no_valid_template` (the converse) and
`FaceObstruction.spurious_eigenpair` (eigenvalues = x_h at the points), this is
Theorem 3.
-/
import FaceObstruction
import EigenSize

open Finset Module

noncomputable section
open Classical

namespace Theorem3

open FaceObstruction EigenSize

variable {K : Type*} [Field K] {β : Type*} [Fintype β] [DecidableEq β] {W : Type*} [AddCommGroup W] [Module K W]

/-- Abstract assembly: point vectors spanning the kernel and independent on `E = Dᶜ`
give corank = number of points, feasibility, and minimal eigen size = number of points. -/
theorem eigen_size_of_points {ζ : Type*} [Fintype ζ]
    (M : (β → K) →ₗ[K] W) (D : Finset β) (φ : ζ → β → K)
    (hker : ∀ z, φ z ∈ LinearMap.ker M)
    (H_span : ∀ v ∈ LinearMap.ker M, v ∈ Submodule.span K (Set.range φ))
    (H_sep : LinearIndependent K (fun z => fun b : (Dᶜ : Finset β) => φ z b)) :
    finrank K (LinearMap.ker M) = Fintype.card ζ ∧ InjOnSupp M D ∧
      (∃ G, D ⊆ G ∧ InjOnSupp M G ∧ (Gᶜ).card = Fintype.card ζ) ∧
      ∀ G, InjOnSupp M G → Fintype.card ζ ≤ (Gᶜ).card := by
  -- the point vectors are independent (their restrictions already are)
  let res : (β → K) →ₗ[K] ((Dᶜ : Finset β) → K) :=
    { toFun := fun v b => v b, map_add' := fun _ _ => rfl, map_smul' := fun _ _ => rfl }
  have hind : LinearIndependent K φ := LinearIndependent.of_comp res H_sep
  -- the kernel is exactly their span
  have hspan_eq : LinearMap.ker M = Submodule.span K (Set.range φ) := by
    apply le_antisymm
    · intro v hv; exact H_span v hv
    · rw [Submodule.span_le]
      rintro _ ⟨z, rfl⟩
      exact hker z
  have hrank : finrank K (LinearMap.ker M) = Fintype.card ζ := by
    rw [hspan_eq, finrank_span_eq_card hind]
  -- feasibility: a kernel vector vanishing on `E = Dᶜ` is zero
  have hfeas : InjOnSupp M D := by
    intro v hv hvD
    obtain ⟨c, rfl⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp (H_span v hv)
    have hc : ∀ z, c z = 0 := by
      apply Fintype.linearIndependent_iff.mp H_sep c
      funext ⟨b, hb⟩
      have := hvD b (Finset.mem_compl.mp hb)
      simpa [Finset.sum_apply] using this
    simp [hc]
  obtain ⟨hopt, hlow⟩ := min_eigen_size_eq_corank M D hfeas
  refine ⟨hrank, hfeas, ?_, ?_⟩
  · obtain ⟨G, hDG, hG, hcard⟩ := hopt
    exact ⟨G, hDG, hG, hcard.trans hrank⟩
  · intro G hG
    have := hlow G hG
    rwa [hrank] at this

/-! ### Concrete templates: rows `x^(t r) * p r`, columns a finite set `C` of exponents. -/

variable {n : ℕ}

/-- The template matrix as a linear map from column vectors to row values. -/
def tmplMap {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n) (C : Finset (Exp n)) :
    (↥C → K) →ₗ[K] (ι → K) where
  toFun v := fun r => ∑ b : ↥C, entry t p r b * v b
  map_add' v v' := by
    funext r
    simp [mul_add, Finset.sum_add_distrib]
  map_smul' c v := by
    funext r
    simp [Finset.mul_sum, mul_left_comm]

/-- A solution point `(w, ξ)` (torus point for `w = 0`, boundary point otherwise)
gives a kernel vector of the concrete template (`FaceObstruction.face_kernel`). -/
theorem point_mem_ker {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n)
    (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n) (ξ : Fin n → K) (hξ : ∀ j, ξ j ≠ 0)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hinit : ∀ r, ∀ k : ℤ, (∀ α ∈ (p r).supp, wdeg w α ≤ k) →
        ∑ α ∈ (p r).supp.filter (fun α => wdeg w α = k), (p r).coeff α * mpow ξ α = 0) :
    (fun b : ↥C => faceVec C hC w ξ b) ∈ LinearMap.ker (tmplMap t p C) := by
  rw [LinearMap.mem_ker]
  funext r
  have h := face_kernel t p C hC w ξ hξ hsupp hinit r
  simp only [tmplMap, LinearMap.coe_mk, AddHom.coe_mk, Pi.zero_apply]
  rw [Finset.sum_coe_sort C (fun b => entry t p r b * faceVec C hC w ξ b)]
  exact h

/-- **Theorem 3 (simple points, R3 and R1/R2 as hypotheses).** Let `Z` be a finite
family of solution points `(w z, ξ z)` of the system on the toric variety, each a face
solution of `in_(w z) F` (with `w z = 0` for torus points). If their vectors span the
kernel (R3) and stay independent on the shiftable columns `Dᶜ` (R1/R2), then the
corank is `|Z|`, the template is feasible, and the minimal eigen size is exactly `|Z|`. -/
theorem toric_eigen_size {ι : Type*} {ζ : Type*} [Fintype ζ]
    (t : ι → Exp n) (p : ι → SPoly K n) (C : Finset (Exp n)) (hC : C.Nonempty)
    (D : Finset ↥C)
    (w : ζ → Exp n) (ξ : ζ → Fin n → K) (hξ : ∀ z j, ξ z j ≠ 0)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hinit : ∀ z r, ∀ k : ℤ, (∀ α ∈ (p r).supp, wdeg (w z) α ≤ k) →
        ∑ α ∈ (p r).supp.filter (fun α => wdeg (w z) α = k), (p r).coeff α * mpow (ξ z) α = 0)
    (H_span : ∀ v ∈ LinearMap.ker (tmplMap t p C),
        v ∈ Submodule.span K (Set.range fun z => fun b : ↥C =>
          faceVec C hC (w z) (ξ z) b))
    (H_sep : LinearIndependent K (fun z => fun b : (Dᶜ : Finset ↥C) =>
        faceVec C hC (w z) (ξ z) b)) :
    finrank K (LinearMap.ker (tmplMap t p C)) = Fintype.card ζ ∧ InjOnSupp (tmplMap t p C) D ∧
      (∃ G, D ⊆ G ∧ InjOnSupp (tmplMap t p C) G ∧ (Gᶜ).card = Fintype.card ζ) ∧
      ∀ G, InjOnSupp (tmplMap t p C) G → Fintype.card ζ ≤ (Gᶜ).card :=
  eigen_size_of_points (tmplMap t p C) D
    (fun z => fun b : ↥C => faceVec C hC (w z) (ξ z) b)
    (fun z => point_mem_ker t p C hC (w z) (ξ z) (hξ z) hsupp (hinit z)) H_span H_sep

end Theorem3
