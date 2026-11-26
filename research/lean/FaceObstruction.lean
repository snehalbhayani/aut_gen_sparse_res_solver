/-
Face obstruction theorem for sparse-resultant (extra polynomial) templates.
See research/proof/face_obstruction.md for the informal proof.

STATUS: written against Mathlib but NOT YET COMPILED. Lean could not be
installed in the session where this was written. Expect small API fixes when
building with `lake build`.

Everything is over an arbitrary field `K`, with arbitrary coefficients (no
genericity), for an arbitrary finite column set `C` (no convexity).
-/
import Mathlib

open Finset

noncomputable section
open Classical

namespace FaceObstruction

variable {K : Type*} [Field K] {n : ℕ}

/-- Exponent vectors (Laurent, so shifts and differences stay in the type). -/
abbrev Exp (n : ℕ) := Fin n → ℤ

/-- The `w`-degree of a monomial `x^b`. -/
def wdeg (w b : Exp n) : ℤ := ∑ j, w j * b j

/-- `ξ^b` for a point `ξ` of the torus. -/
def mpow (ξ : Fin n → K) (b : Exp n) : K := ∏ j, ξ j ^ b j

/-- A sparse (Laurent) polynomial, given by a finite support and coefficients. -/
structure SPoly (K : Type*) (n : ℕ) where
  supp : Finset (Exp n)
  coeff : Exp n → K

/-- Template entry: row `r` is the polynomial `x^(t r) * p r`; column `b` is `x^b`. -/
def entry {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n) (r : ι) (b : Exp n) : K :=
  if b - t r ∈ (p r).supp then (p r).coeff (b - t r) else 0

/-- The kernel vector supported on the `w`-top layer of the columns. -/
def faceVec (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n) (ξ : Fin n → K)
    (b : Exp n) : K :=
  if b ∈ C ∧ wdeg w b = C.sup' hC (wdeg w) then mpow ξ b else 0

lemma wdeg_add (w a b : Exp n) : wdeg w (a + b) = wdeg w a + wdeg w b := by
  simp [wdeg, mul_add, Finset.sum_add_distrib]

lemma wdeg_single (w : Exp n) (h : Fin n) : wdeg w (Pi.single h 1) = w h := by
  simp [wdeg, Pi.single_apply]

lemma mpow_add {ξ : Fin n → K} (hξ : ∀ j, ξ j ≠ 0) (a b : Exp n) :
    mpow ξ (a + b) = mpow ξ a * mpow ξ b := by
  unfold mpow
  rw [← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl (fun j _ => ?_)
  simp [zpow_add₀ (hξ j)]

lemma mpow_ne_zero {ξ : Fin n → K} (hξ : ∀ j, ξ j ≠ 0) (b : Exp n) : mpow ξ b ≠ 0 := by
  unfold mpow
  exact Finset.prod_ne_zero_iff.mpr (fun j _ => zpow_ne_zero _ (hξ j))

lemma mpow_single (ξ : Fin n → K) (h : Fin n) : mpow ξ (Pi.single h 1) = ξ h := by
  unfold mpow
  rw [Finset.prod_eq_single h]
  · simp
  · intro b _ hb
    simp [Pi.single_eq_of_ne hb]
  · intro hh
    exact absurd (Finset.mem_univ h) hh

lemma faceVec_top_ne_zero (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n)
    {ξ : Fin n → K} (hξ : ∀ j, ξ j ≠ 0) :
    ∃ b ∈ C, faceVec C hC w ξ b ≠ 0 := by
  obtain ⟨b, hb, hbmax⟩ := Finset.exists_mem_eq_sup' hC (wdeg w)
  refine ⟨b, hb, ?_⟩
  have hc : b ∈ C ∧ wdeg w b = C.sup' hC (wdeg w) := ⟨hb, hbmax.symm⟩
  simp only [faceVec, if_pos hc]
  exact mpow_ne_zero hξ b

/-- **Key lemma.** If `ξ` is a face solution for `w` (every `w`-component of
degree at least the top degree vanishes at `ξ`), then `faceVec` lies in the kernel
of the template matrix. -/
theorem face_kernel {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n)
    (C : Finset (Exp n)) (hC : C.Nonempty) (w : Exp n) (ξ : Fin n → K)
    (hξ : ∀ j, ξ j ≠ 0)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hinit : ∀ r, ∀ k : ℤ, (∀ α ∈ (p r).supp, wdeg w α ≤ k) →
        ∑ α ∈ (p r).supp.filter (fun α => wdeg w α = k), (p r).coeff α * mpow ξ α = 0)
    (r : ι) :
    ∑ b ∈ C, entry t p r b * faceVec C hC w ξ b = 0 := by
  set D := C.sup' hC (wdeg w) with hD
  have himg : (p r).supp.image (fun α => t r + α) ⊆ C := by
    intro b hb
    obtain ⟨α, hα, rfl⟩ := Finset.mem_image.mp hb
    exact hsupp r α hα
  -- columns outside the translate `t r + supp` contribute nothing
  have hzero : ∀ b ∈ C, b ∉ (p r).supp.image (fun α => t r + α) →
      entry t p r b * faceVec C hC w ξ b = 0 := by
    intro b _ hbnot
    have hnot : b - t r ∉ (p r).supp := by
      intro hmem
      exact hbnot (Finset.mem_image.mpr ⟨b - t r, hmem, by abel⟩)
    simp [entry, hnot]
  rw [← Finset.sum_subset himg hzero]
  rw [Finset.sum_image (by intro a _ b _ hab; exact add_left_cancel hab)]
  have key : ∀ α ∈ (p r).supp,
      entry t p r (t r + α) * faceVec C hC w ξ (t r + α)
        = mpow ξ (t r) *
          (if wdeg w α = D - wdeg w (t r) then (p r).coeff α * mpow ξ α else 0) := by
    intro α hα
    have hmem : t r + α ∈ C := hsupp r α hα
    have e1 : entry t p r (t r + α) = (p r).coeff α := by
      simp [entry, hα]
    have e2 : faceVec C hC w ξ (t r + α)
        = if wdeg w (t r) + wdeg w α = D then mpow ξ (t r + α) else 0 := by
      simp [faceVec, hmem, wdeg_add, hD]
    rw [e1, e2]
    by_cases hdeg : wdeg w α = D - wdeg w (t r)
    · have h' : wdeg w (t r) + wdeg w α = D := by omega
      rw [if_pos h', if_pos hdeg, mpow_add hξ]
      ring
    · have h' : ¬ (wdeg w (t r) + wdeg w α = D) := by omega
      rw [if_neg h', if_neg hdeg]
      ring
  rw [Finset.sum_congr rfl key, ← Finset.mul_sum, ← Finset.sum_filter]
  rw [hinit r (D - wdeg w (t r)) ?_, mul_zero]
  intro α hα
  have hle := Finset.le_sup' (wdeg w) (hsupp r α hα)
  rw [wdeg_add] at hle
  linarith

/-- **Theorem 1 (obstruction).** A face solution for a weight with `w h > 0` makes
every template with hidden variable `x_h` (eigenvalue `λ = x_h`) invalid: there is
a nonzero kernel vector supported in `G = C \ E`. -/
theorem no_valid_template {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n)
    (C E : Finset (Exp n)) (hC : C.Nonempty) (h : Fin n)
    (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (w : Exp n) (hw : 0 < w h) (ξ : Fin n → K) (hξ : ∀ j, ξ j ≠ 0)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hinit : ∀ r, ∀ k : ℤ, (∀ α ∈ (p r).supp, wdeg w α ≤ k) →
        ∑ α ∈ (p r).supp.filter (fun α => wdeg w α = k), (p r).coeff α * mpow ξ α = 0) :
    ∃ v : Exp n → K, (∃ b, v b ≠ 0) ∧ (∀ b, v b ≠ 0 → b ∈ C ∧ b ∉ E) ∧
      ∀ r, ∑ b ∈ C, entry t p r b * v b = 0 := by
  refine ⟨faceVec C hC w ξ, ?_, ?_, fun r => face_kernel t p C hC w ξ hξ hsupp hinit r⟩
  · obtain ⟨b, _, hb⟩ := faceVec_top_ne_zero C hC w hξ
    exact ⟨b, hb⟩
  · intro b hb
    by_cases hc : b ∈ C ∧ wdeg w b = C.sup' hC (wdeg w)
    · refine ⟨hc.1, fun hbE => ?_⟩
      have h1 := Finset.le_sup' (wdeg w) (hE b hbE)
      rw [wdeg_add, wdeg_single, hc.2] at h1
      linarith
    · exact absurd (show faceVec C hC w ξ b = 0 from if_neg hc) hb

/-- **Theorem 2 (spurious eigenpairs).** For a valid template, a face solution for
a weight with `w h ≤ 0` gives a nonzero vector on `E`. Its shift satisfies
`v (e + e_h) = μ * v e` with `μ = ξ h` (if `w h = 0`) or `μ = 0` (if `w h < 0`),
i.e. it is an eigenvector of every solver matrix `X` with eigenvalue `μ`. -/
theorem spurious_eigenpair {ι : Type*} (t : ι → Exp n) (p : ι → SPoly K n)
    (C E : Finset (Exp n)) (hC : C.Nonempty) (h : Fin n)
    (hEC : ∀ e ∈ E, e ∈ C) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (w : Exp n) (hw : w h ≤ 0) (ξ : Fin n → K) (hξ : ∀ j, ξ j ≠ 0)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hinit : ∀ r, ∀ k : ℤ, (∀ α ∈ (p r).supp, wdeg w α ≤ k) →
        ∑ α ∈ (p r).supp.filter (fun α => wdeg w α = k), (p r).coeff α * mpow ξ α = 0)
    (hvalid : ∀ u : Exp n → K, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0) :
    (∃ e ∈ E, faceVec C hC w ξ e ≠ 0) ∧
      ∀ e ∈ E, faceVec C hC w ξ (e + Pi.single h 1)
        = (if w h = 0 then ξ h else 0) * faceVec C hC w ξ e := by
  constructor
  · by_contra hcon
    push_neg at hcon
    obtain ⟨b, _, hb⟩ := faceVec_top_ne_zero C hC w hξ
    apply hb
    refine hvalid (faceVec C hC w ξ) ?_ (fun r => face_kernel t p C hC w ξ hξ hsupp hinit r) b
    intro c hc
    refine ⟨?_, fun hcE => hc (hcon c hcE)⟩
    by_contra hcC
    apply hc
    have hn : ¬ (c ∈ C ∧ wdeg w c = C.sup' hC (wdeg w)) := fun hh => hcC hh.1
    exact if_neg hn
  · intro e he
    have hmem := hE e he
    have hdeg : wdeg w (e + Pi.single h 1) = wdeg w e + w h := by
      rw [wdeg_add, wdeg_single]
    by_cases hw0 : w h = 0
    · rw [if_pos hw0]
      simp only [faceVec, hmem, hEC e he, true_and, hdeg, hw0, add_zero]
      split_ifs
      · rw [mpow_add hξ, mpow_single]
        ring
      · ring
    · rw [if_neg hw0, zero_mul]
      have hlt : w h < 0 := lt_of_le_of_ne hw hw0
      have hle := Finset.le_sup' (wdeg w) (hEC e he)
      simp only [faceVec, hmem, true_and, hdeg]
      rw [if_neg]
      intro hcontra
      linarith

end FaceObstruction
