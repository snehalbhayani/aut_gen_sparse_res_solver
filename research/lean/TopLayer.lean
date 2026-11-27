/-
Theorem P1a of `research/proof/lower_bounds.md` (face-certificate bound), machine-checked.

For a valid template with hidden `x_h` and a weight `w` with `w h > 0`:

* the top `w`-layer `L` of the columns contains no shiftable column, so `L ⊆ G`;
* the entries of any row on `L` come from the top face of its polynomial, i.e. from the
  initial form `in_w` (`layer_entry_initial`);
* the rows meeting `L`, restricted to `L`, are injective on `K^L`, so at least `|L|` rows
  meet `L` (`top_layer_card_le`).

Together: `M[:, L]` is a complete template of the initial system `in_w(F)` on `L`, and
`|G| ≥ |L|`. This is the input of the lower bounds P1b′ and P1e.
-/
import FaceObstruction

open Finset Module

noncomputable section
open Classical

namespace TopLayer

open FaceObstruction

variable {K : Type*} [Field K] {n : ℕ}

/-- The top `w`-layer of the column set. -/
def layer (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n) : Finset (Exp n) :=
  C.filter (fun b => wdeg w b = C.sup' hC (wdeg w))

lemma layer_subset (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n) :
    layer C hC w ⊆ C := Finset.filter_subset _ _

/-- For `w h > 0` the top layer contains no shiftable column, so it lies in `G`. -/
lemma layer_disjoint (C E : Finset (Exp n)) (hC : C.Nonempty) (h : Fin n)
    (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C) (w : Exp n) (hw : 0 < w h) :
    ∀ b ∈ layer C hC w, b ∉ E := by
  intro b hb hbE
  have h1 := Finset.le_sup' (wdeg w) (hE b hbE)
  rw [wdeg_add, wdeg_single] at h1
  have h2 := (Finset.mem_filter.mp hb).2
  linarith

/-- Entries of a row on the top layer come from the top face of its polynomial
(the initial form `in_w`). -/
lemma layer_entry_initial {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n)
    (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C) (r : ι) (b : Exp n)
    (hb : b ∈ layer C hC w) (hne : entry t p r b ≠ 0) :
    b - t r ∈ (p r).supp ∧ ∀ α ∈ (p r).supp, wdeg w α ≤ wdeg w (b - t r) := by
  have hmem : b - t r ∈ (p r).supp := by
    by_contra hn
    exact hne (by simp [entry, hn])
  refine ⟨hmem, fun α hα => ?_⟩
  have h1 := Finset.le_sup' (wdeg w) (hsupp r α hα)
  have h2 := (Finset.mem_filter.mp hb).2
  have hsplit : wdeg w b = wdeg w (t r) + wdeg w (b - t r) := by
    rw [← wdeg_add]
    congr 1
    abel
  rw [wdeg_add] at h1
  linarith

/-- **Theorem P1a (machine-checked).** In a valid template (no non-zero kernel vector
supported in `G = C \ E`) with hidden `x_h` and a weight with `w h > 0`, at least `|L|`
rows meet the top layer `L`, because the rows meeting `L`, restricted to `L`, are
injective on `K^L`. With `layer_entry_initial` these restricted rows are the rows of the
initial system `in_w(F)`, so `L` is a complete set for `in_w(F)`. -/
theorem top_layer_card_le {ι : Type*} [Fintype ι] (t : ι → Exp n) (p : ι → SPoly K n)
    (C E : Finset (Exp n)) (hC : C.Nonempty) (h : Fin n)
    (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C) (w : Exp n) (hw : 0 < w h)
    (hvalid : ∀ u : Exp n → K, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0) :
    (layer C hC w).card ≤
      (Finset.univ.filter (fun r => ∃ b ∈ layer C hC w, entry t p r b ≠ 0)).card := by
  set L := layer C hC w with hL
  set R := Finset.univ.filter (fun r => ∃ b ∈ L, entry t p r b ≠ 0) with hR
  let Φ : (↥L → K) →ₗ[K] (↥R → K) :=
    { toFun := fun u r => ∑ b : ↥L, entry t p r b * u b
      map_add' := by
        intro u v
        funext r
        simp [mul_add, Finset.sum_add_distrib]
      map_smul' := by
        intro c u
        funext r
        simp [Finset.mul_sum, mul_left_comm] }
  have hinj : Function.Injective Φ := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro u hu
    let v : Exp n → K := fun b => if hb : b ∈ L then u ⟨b, hb⟩ else 0
    have hvsum : ∀ r, ∑ b ∈ C, entry t p r b * v b = ∑ b : ↥L, entry t p r b * u b := by
      intro r
      rw [← Finset.sum_subset (layer_subset C hC w) (by
        intro b _ hb
        simp [v, hL, hb])]
      rw [← Finset.sum_coe_sort L]
      refine Finset.sum_congr rfl (fun b _ => ?_)
      simp [v, b.2]
    have hv : ∀ b, v b = 0 := by
      apply hvalid v
      · intro b hb
        have hbL : b ∈ L := by
          by_contra hn
          exact hb (by simp [v, hn])
        exact ⟨layer_subset C hC w hbL, layer_disjoint C E hC h hE w hw b hbL⟩
      · intro r
        rw [hvsum r]
        by_cases hr : r ∈ R
        · have := congrFun hu ⟨r, hr⟩
          simpa [Φ] using this
        · apply Finset.sum_eq_zero
          intro b _
          have h0 : entry t p r b = 0 := by
            by_contra hne
            exact hr (Finset.mem_filter.mpr ⟨Finset.mem_univ _, b, b.2, hne⟩)
          simp [h0]
    funext b
    have := hv b
    simpa [v, b.2] using this
  have hle := LinearMap.finrank_le_finrank_of_injective hinj
  simpa [Module.finrank_fintype_fun_eq_card, Fintype.card_coe] using hle

end TopLayer
