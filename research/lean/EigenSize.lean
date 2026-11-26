/-
Step 3 of `research/proof/exact_count.md`: the eigen-basis reduction reaches
exactly the corank, and no valid template can do better.

Abstract setting: columns are a finite type `β`, and the template is any linear map
`M : (β → K) →ₗ[K] W` (column vectors ↦ row values). For a set `S` of columns,
"`M` is injective on vectors supported in `S`" is the validity of a template whose
eliminated columns are `S` (so its eigen basis is the complement `Sᶜ`).
-/
import Mathlib

open Finset Module

noncomputable section
open Classical

namespace EigenSize

variable {K : Type*} [Field K] {β : Type*} [Fintype β] {W : Type*} [AddCommGroup W] [Module K W]

/-- `M` is injective on vectors vanishing outside `S` (validity with `G = S`). -/
def InjOnSupp (M : (β → K) →ₗ[K] W) (S : Finset β) : Prop :=
  ∀ v ∈ LinearMap.ker M, (∀ b ∉ S, v b = 0) → v = 0

/-- Restriction of kernel vectors to a set of columns. -/
def restr (M : (β → K) →ₗ[K] W) (E : Finset β) : LinearMap.ker M →ₗ[K] (E → K) where
  toFun v := fun e => (v : β → K) e
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

lemma restr_injective (M : (β → K) →ₗ[K] W) (S : Finset β) (h : InjOnSupp M S) :
    Function.Injective (restr M Sᶜ) := by
  rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
  intro v hv
  have h0 : ∀ b ∉ S, (v : β → K) b = 0 := by
    intro b hb
    have := congrFun hv ⟨b, Finset.mem_compl.mpr hb⟩
    simpa [restr] using this
  exact Subtype.ext (h v v.2 h0)

/-- **Lower bound.** Every valid template (eliminated columns `S`) has an eigen
basis `Sᶜ` at least as large as the kernel of `M`, i.e. the corank. -/
theorem finrank_ker_le (M : (β → K) →ₗ[K] W) (S : Finset β) (h : InjOnSupp M S) :
    finrank K (LinearMap.ker M) ≤ (Sᶜ).card := by
  have := LinearMap.finrank_le_finrank_of_injective (restr_injective M S h)
  rwa [Module.finrank_fintype_fun_eq_card, Fintype.card_coe] at this

/-- **Achievability (symgen's matroid step).** If the boundary columns `D` are
eliminable, they extend to a set `G ⊇ D` that is still eliminable, and whose
eigen basis `Gᶜ` has exactly the size of the kernel of `M`. -/
theorem exists_optimal (M : (β → K) →ₗ[K] W) (D : Finset β) (hD : InjOnSupp M D) :
    ∃ G : Finset β, D ⊆ G ∧ InjOnSupp M G ∧ (Gᶜ).card = finrank K (LinearMap.ker M) := by
  let 𝒢 := (Finset.univ : Finset (Finset β)).filter (fun G => D ⊆ G ∧ InjOnSupp M G)
  have hne : 𝒢.Nonempty := ⟨D, by simp [𝒢, hD]⟩
  obtain ⟨G, hG𝒢, hGmax⟩ := Finset.exists_max_image 𝒢 Finset.card hne
  have hG : D ⊆ G ∧ InjOnSupp M G := by simpa [𝒢] using hG𝒢
  refine ⟨G, hG.1, hG.2, le_antisymm ?_ (finrank_ker_le M G hG.2)⟩
  -- the restriction map to `Gᶜ` is surjective
  have hsurj : Function.Surjective (restr M Gᶜ) := by
    have hsingle : ∀ b : (Gᶜ : Finset β), ∃ u : LinearMap.ker M,
        restr M Gᶜ u = Pi.single b 1 := by
      intro ⟨b, hb⟩
      have hbG : b ∉ G := Finset.mem_compl.mp hb
      -- `insert b G` is strictly larger, so it is not eliminable
      have hnot : ¬ InjOnSupp M (insert b G) := by
        intro hins
        have hmem : insert b G ∈ 𝒢 := by
          simp only [𝒢, Finset.mem_filter, Finset.mem_univ, true_and]
          exact ⟨hG.1.trans (Finset.subset_insert b G), hins⟩
        have := hGmax _ hmem
        rw [Finset.card_insert_of_not_mem hbG] at this
        omega
      simp only [InjOnSupp, not_forall] at hnot
      obtain ⟨v, hvker, hvsupp, hvne⟩ := hnot
      have hvb : v b ≠ 0 := by
        intro hvb0
        apply hvne
        refine hG.2 v hvker (fun c hc => ?_)
        by_cases hcb : c = b
        · rw [hcb]; exact hvb0
        · exact hvsupp c (by simp [hcb, hc])
      refine ⟨⟨(v b)⁻¹ • v, Submodule.smul_mem _ _ hvker⟩, ?_⟩
      funext ⟨c, hc⟩
      have hcG : c ∉ G := Finset.mem_compl.mp hc
      simp only [restr, LinearMap.coe_mk, AddHom.coe_mk, Pi.smul_apply, smul_eq_mul]
      by_cases hcb : c = b
      · subst hcb
        simp [hvb]
      · have hv0 : v c = 0 := hvsupp c (by simp [hcb, hcG])
        have hne' : (⟨c, hc⟩ : (Gᶜ : Finset β)) ≠ ⟨b, hb⟩ := fun h => hcb (congrArg Subtype.val h)
        simp [hv0, Pi.single_apply, hne']
    intro f
    have hspan : f = ∑ b : (Gᶜ : Finset β), f b • (Pi.single b 1 : (Gᶜ : Finset β) → K) := by
      funext c
      simp [Finset.sum_apply, Pi.single_apply]
    choose u hu using hsingle
    refine ⟨∑ b, f b • u b, ?_⟩
    rw [map_sum]
    simp only [map_smul, hu]
    exact hspan.symm
  have := LinearMap.finrank_range_le (restr M Gᶜ)
  rw [LinearMap.range_eq_top.mpr hsurj, finrank_top, Module.finrank_fintype_fun_eq_card,
    Fintype.card_coe] at this
  exact this

/-- **Step 3.** For a feasible template, the minimal eigen size over all valid
choices of eliminated columns containing the boundary `D` is exactly the corank
`finrank (ker M)`. -/
theorem min_eigen_size_eq_corank (M : (β → K) →ₗ[K] W) (D : Finset β) (hD : InjOnSupp M D) :
    (∃ G, D ⊆ G ∧ InjOnSupp M G ∧ (Gᶜ).card = finrank K (LinearMap.ker M)) ∧
    ∀ G, InjOnSupp M G → finrank K (LinearMap.ker M) ≤ (Gᶜ).card :=
  ⟨exists_optimal M D hD, fun G hG => finrank_ker_le M G hG⟩

end EigenSize
