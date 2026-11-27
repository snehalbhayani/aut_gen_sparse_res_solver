/-
Theorem C2 of `research/proof/complexity.md` (NP-hardness of the minimal elimination
template), the algebraic core of the reduction from SET COVER.

Part 1 (this section): the sign-trick lemma. If the sum of the element variables lies in
the span of positive-gain gadget rows `e_y + r e_x`, then those rows span every variable
they touch, so there are at least as many rows as touched variables, and the touched
element variables are exactly all elements.
-/
import Mathlib
import FaceObstruction

open Finset Module

noncomputable section
open Classical

namespace Hardness

variable {V : Type*} [Fintype V] [DecidableEq V] {ε : Type*}

/-- Gadget row `e_{y} + r e_{x}` (edge between a set variable and an element variable). -/
def gvec (yv xv : ε → V) (r : ε → ℚ) (e : ε) : V → ℚ :=
  (Pi.single (yv e) (1 : ℚ) : V → ℚ) + r e • (Pi.single (xv e) (1 : ℚ) : V → ℚ)

/-- The functional `g ↦ ∑ b, c b * g b`. -/
def dotL (c : V → ℚ) : (V → ℚ) →ₗ[ℚ] ℚ := ∑ b, c b • LinearMap.proj b

lemma dotL_apply (c g : V → ℚ) : dotL c g = ∑ b, c b * g b := by
  simp [dotL]

lemma dotL_single (c : V → ℚ) (b : V) : dotL c (Pi.single b (1 : ℚ) : V → ℚ) = c b := by
  rw [dotL_apply]
  simp [Pi.single_apply]

lemma dotL_single_left (x : V) (g : V → ℚ) : dotL (Pi.single x (1 : ℚ) : V → ℚ) g = g x := by
  rw [dotL_apply]
  simp [Pi.single_apply]

/-- **Sign-trick lemma.** Let element variables `X` and set variables be disjoint, let every
gadget row `gvec e = e_{yv e} + r e · e_{xv e}` have `xv e ∈ X`, `yv e ∉ X` and `r e > 0`.
If `∑_{x ∈ X} e_x` lies in the span of the rows indexed by `U`, then every element is
touched, and `|X| + |set variables touched| ≤ |U|`. -/
theorem gadget_count (X : Finset V) (yv xv : ε → V) (r : ε → ℚ)
    (hr : ∀ e, 0 < r e) (hx : ∀ e, xv e ∈ X) (hy : ∀ e, yv e ∉ X) (U : Finset ε)
    (hs : (∑ x ∈ X, Pi.single x (1 : ℚ)) ∈
      Submodule.span ℚ ((U.image (gvec yv xv r) : Finset (V → ℚ)) : Set (V → ℚ))) :
    X ⊆ U.image xv ∧ X.card + (U.image yv).card ≤ U.card := by
  set S := Submodule.span ℚ ((U.image (gvec yv xv r) : Finset (V → ℚ)) : Set (V → ℚ)) with hS
  set s : V → ℚ := ∑ x ∈ X, Pi.single x (1 : ℚ) with hsdef
  have hsx : ∀ b, s b = if b ∈ X then 1 else 0 := by
    intro b
    simp [hsdef, Finset.sum_apply, Pi.single_apply]
  -- a functional killing all gadget rows kills `s`
  have kill : ∀ c : V → ℚ, (∀ e ∈ U, c (yv e) + r e * c (xv e) = 0) → dotL c s = 0 := by
    intro c hc
    have hle : S ≤ LinearMap.ker (dotL c) := by
      rw [hS, Submodule.span_le]
      intro g hg
      obtain ⟨e, he, rfl⟩ := Finset.mem_image.mp (Finset.mem_coe.mp hg)
      simp only [SetLike.mem_coe, LinearMap.mem_ker, gvec, map_add, map_smul, dotL_single,
        smul_eq_mul]
      exact hc e he
    exact LinearMap.mem_ker.mp (hle hs)
  -- (1) every element is touched
  have hcover : X ⊆ U.image xv := by
    intro x hxX
    by_contra hnot
    have h0 := kill (Pi.single x (1 : ℚ) : V → ℚ) (fun e he => by
      have hyx : yv e ≠ x := fun h => hy e (h ▸ hxX)
      have hxx : xv e ≠ x := fun h => hnot (Finset.mem_image.mpr ⟨e, he, h⟩)
      simp [Pi.single_apply, hyx, hxx])
    rw [dotL_single_left, hsx, if_pos hxX] at h0
    exact one_ne_zero h0
  refine ⟨hcover, ?_⟩
  -- (2) every touched variable's unit vector lies in the span
  set T := U.image xv ∪ U.image yv with hT
  have hsingle : ∀ b ∈ T, (Pi.single b (1 : ℚ) : V → ℚ) ∈ S := by
    intro b hb
    by_cases hmem : (Pi.single b (1 : ℚ) : V → ℚ) ∈ S
    · exact hmem
    exfalso
    obtain ⟨f, hfb, hfS⟩ := Submodule.exists_dual_map_eq_bot_of_nmem hmem inferInstance
    set π : V → ℚ := fun v => f (Pi.single v (1 : ℚ) : V → ℚ) with hπ
    have hfS' : ∀ g ∈ S, f g = 0 := by
      intro g hg
      have : f g ∈ S.map f := Submodule.mem_map_of_mem hg
      rw [hfS] at this
      exact (Submodule.mem_bot ℚ).mp this
    have hedge : ∀ e ∈ U, π (yv e) + r e * π (xv e) = 0 := by
      intro e he
      have hmem : gvec yv xv r e ∈ S :=
        Submodule.subset_span (Finset.mem_coe.mpr (Finset.mem_image.mpr ⟨e, he, rfl⟩))
      have := hfS' _ hmem
      simpa [gvec, hπ, map_add, map_smul] using this
    -- the sign-flipped functional also kills every row
    set c : V → ℚ := fun v => if v ∈ X then |π v| else -|π v| with hc
    have hc_edge : ∀ e ∈ U, c (yv e) + r e * c (xv e) = 0 := by
      intro e he
      have h1 : π (yv e) = -(r e * π (xv e)) := by linarith [hedge e he]
      have habs : |π (yv e)| = r e * |π (xv e)| := by
        rw [h1, abs_neg, abs_mul, abs_of_pos (hr e)]
      simp only [hc, if_neg (hy e), if_pos (hx e)]
      rw [habs]; ring
    have hsum := kill c hc_edge
    rw [dotL_apply] at hsum
    have hcs : ∑ b, c b * s b = ∑ x ∈ X, |π x| := by
      simp only [hsx, mul_ite, mul_one, mul_zero]
      rw [Finset.sum_ite_mem, Finset.univ_inter]
      refine Finset.sum_congr rfl (fun x hxX => ?_)
      simp [hc, hxX]
    have hsum' : ∑ x ∈ X, |π x| = 0 := hcs ▸ hsum
    have hzeroX : ∀ x ∈ X, π x = 0 := by
      intro x hxX
      have := (Finset.sum_eq_zero_iff_of_nonneg (fun i _ => abs_nonneg (π i))).mp hsum' x hxX
      exact abs_eq_zero.mp this
    have hπb : π b = 0 := by
      rcases Finset.mem_union.mp hb with hbx | hby
      · obtain ⟨e, _, rfl⟩ := Finset.mem_image.mp hbx
        exact hzeroX _ (hx e)
      · obtain ⟨e, he, rfl⟩ := Finset.mem_image.mp hby
        have := hedge e he
        rw [hzeroX _ (hx e), mul_zero, add_zero] at this
        exact this
    exact hfb hπb
  -- dimension count
  have hli : LinearIndependent ℚ (fun b : (T : Finset V) => (Pi.single b.1 (1 : ℚ) : V → ℚ)) := by
    have := (Pi.basisFun ℚ V).linearIndependent.comp (fun b : (T : Finset V) => b.1)
      Subtype.val_injective
    simpa [Function.comp_def, Pi.basisFun_apply] using this
  have hspan_le : Submodule.span ℚ (Set.range fun b : (T : Finset V) =>
      (Pi.single b.1 (1 : ℚ) : V → ℚ)) ≤ S := by
    rw [Submodule.span_le]
    rintro _ ⟨b, rfl⟩
    exact hsingle b.1 b.2
  have hT_le : T.card ≤ U.card := by
    have h1 := Submodule.finrank_mono hspan_le
    rw [finrank_span_eq_card hli, Fintype.card_coe] at h1
    have h2 := finrank_span_finset_le_card (R := ℚ) (U.image (gvec yv xv r))
    have h3 : (U.image (gvec yv xv r)).card ≤ U.card := Finset.card_image_le
    exact h1.trans (h2.trans h3)
  -- |T| = |X| + |touched set variables|
  have himg : U.image xv = X := Finset.Subset.antisymm
    (fun v hv => by obtain ⟨e, _, rfl⟩ := Finset.mem_image.mp hv; exact hx e) hcover
  have hdisj : Disjoint (U.image xv) (U.image yv) := by
    rw [Finset.disjoint_left]
    intro v hv hv'
    obtain ⟨e, _, rfl⟩ := Finset.mem_image.mp hv'
    rw [himg] at hv
    exact hy e hv
  have hTcard : T.card = X.card + (U.image yv).card := by
    rw [hT, Finset.card_union_of_disjoint hdisj, himg]
  omega

/-! ## Lemma L: the used equations pin down the hidden variable -/

namespace LemmaL

open FaceObstruction

variable {n : ℕ}

lemma wdeg_zero (b : Exp n) : wdeg (0 : Exp n) b = 0 := by simp [wdeg]

lemma faceVec_zero (C : Finset (Exp n)) (hC : C.Nonempty) (ξ : Fin n → ℚ) (b : Exp n) :
    faceVec C hC 0 ξ b = if b ∈ C then mpow ξ b else 0 := by
  have hfun : (wdeg (0 : Exp n)) = fun _ => (0 : ℤ) := funext wdeg_zero
  have hsup : C.sup' hC (wdeg (0 : Exp n)) = 0 := by rw [hfun]; exact Finset.sup'_const hC 0
  unfold faceVec
  rw [hsup, wdeg_zero]
  simp

/-- **Lemma L.** In a valid template whose eigen basis `E` is nonempty and shift-closed into
the columns, `N` torus zeros of all rows with pairwise distinct `x_h`-values force
`N ≤ |E|`. -/
theorem card_le {N : ℕ} {ι : Type*} (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n)) (hEC : E ⊆ C) (h : Fin n) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (hEne : E.Nonempty)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0)
    (ξ : Fin N → Fin n → ℚ) (hξ : ∀ k j, ξ k j ≠ 0)
    (hzero : ∀ k r, ∑ α ∈ (p r).supp, (p r).coeff α * mpow (ξ k) α = 0)
    (hinj : Function.Injective fun k => ξ k h) :
    N ≤ E.card := by
  have hC : C.Nonempty := hEne.mono hEC
  set pv : Fin N → Exp n → ℚ := fun k b => if b ∈ C then mpow (ξ k) b else 0 with hpv
  have pv_ker : ∀ k r, ∑ b ∈ C, entry t p r b * pv k b = 0 := by
    intro k r
    have hk := face_kernel t p C hC 0 (ξ k) (hξ k) hsupp (fun r' k' _ => by
      by_cases hk0 : k' = 0
      · subst hk0
        have : (p r').supp.filter (fun α => wdeg (0 : Exp n) α = 0) = (p r').supp :=
          Finset.filter_true_of_mem (fun α _ => wdeg_zero α)
        rw [this]; exact hzero k r'
      · have : (p r').supp.filter (fun α => wdeg (0 : Exp n) α = k') = ∅ :=
          Finset.filter_false_of_mem (fun α _ => by rw [wdeg_zero]; exact Ne.symm hk0)
        rw [this, Finset.sum_empty]) r
    simpa [faceVec_zero, hpv] using hk
  -- restrictions of the point vectors to E are linearly independent
  have hli : LinearIndependent ℚ (fun k : Fin N => fun e : (E : Finset (Exp n)) => mpow (ξ k) e.1) := by
    rw [Fintype.linearIndependent_iff]
    intro g hg
    have hgE : ∀ e ∈ E, ∑ k, g k * mpow (ξ k) e = 0 := by
      intro e he
      have := congrFun hg ⟨e, he⟩
      simpa [Finset.sum_apply] using this
    set u : ℕ → Exp n → ℚ := fun j b => ∑ k, g k * (ξ k h) ^ j * pv k b with hu
    have u_ker : ∀ j r, ∑ b ∈ C, entry t p r b * u j b = 0 := by
      intro j r
      simp only [hu, Finset.mul_sum]
      rw [Finset.sum_comm]
      refine Finset.sum_eq_zero (fun k _ => ?_)
      have := pv_ker k r
      calc ∑ b ∈ C, entry t p r b * (g k * ξ k h ^ j * pv k b)
          = g k * ξ k h ^ j * ∑ b ∈ C, entry t p r b * pv k b := by
            rw [Finset.mul_sum]; refine Finset.sum_congr rfl (fun b _ => by ring)
        _ = 0 := by rw [this, mul_zero]
    have u_out : ∀ j b, b ∉ C → u j b = 0 := by
      intro j b hb; simp [hu, hpv, hb]
    have u_shift : ∀ j e, e ∈ E → u (j + 1) e = u j (e + Pi.single h 1) := by
      intro j e he
      simp only [hu, hpv, if_pos (hEC he), if_pos (hE e he)]
      refine Finset.sum_congr rfl (fun k _ => ?_)
      rw [mpow_add (hξ k), mpow_single, pow_succ]; ring
    have u_zero : ∀ j b, u j b = 0 := by
      intro j
      induction j with
      | zero =>
        refine hvalid (u 0) (fun b hb => ⟨by_contra fun hbC => hb (u_out 0 b hbC), fun hbE => ?_⟩)
          (u_ker 0)
        apply hb
        simp only [hu, pow_zero, mul_one, hpv, if_pos (hEC hbE)]
        exact hgE b hbE
      | succ j ih =>
        refine hvalid (u (j + 1)) (fun b hb => ⟨by_contra fun hbC => hb (u_out _ b hbC),
          fun hbE => hb (by rw [u_shift j b hbE, ih])⟩) (u_ker (j + 1))
    -- Vandermonde at one eigen-basis monomial
    obtain ⟨e0, he0⟩ := hEne
    set a : Fin N → ℚ := fun k => g k * mpow (ξ k) e0 with ha
    have hvec : Matrix.vecMul a (Matrix.vandermonde fun k => ξ k h) = 0 := by
      funext j
      have := u_zero j e0
      simp only [hu, hpv, if_pos (hEC he0)] at this
      simp only [Matrix.vecMul, dotProduct, Matrix.vandermonde_apply, ha, Pi.zero_apply]
      rw [← this]
      refine Finset.sum_congr rfl (fun k _ => by ring)
    have hdet : (Matrix.vandermonde fun k => ξ k h).det ≠ 0 :=
      Matrix.det_vandermonde_ne_zero_iff.mpr hinj
    have ha0 := Matrix.eq_zero_of_vecMul_eq_zero hdet hvec
    intro k
    have := congrFun ha0 k
    simp only [ha, Pi.zero_apply, mul_eq_zero] at this
    rcases this with h1 | h1
    · exact h1
    · exact absurd h1 (mpow_ne_zero (hξ k) e0)
  have := hli.fintype_card_le_finrank
  rw [Fintype.card_fin, Module.finrank_fintype_fun_eq_card, Fintype.card_coe] at this
  exact this

lemma vec_expand {m : ℕ} (v : Fin m → ℚ) : v = ∑ j, v j • (Pi.single j (1 : ℚ) : Fin m → ℚ) := by
  funext i
  simp [Finset.sum_apply, Pi.single_apply]

lemma dual_expand {m : ℕ} (φ : (Fin m → ℚ) →ₗ[ℚ] ℚ) (v : Fin m → ℚ) :
    φ v = ∑ j, φ (Pi.single j (1 : ℚ) : Fin m → ℚ) * v j := by
  conv_lhs => rw [vec_expand v]
  rw [map_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [map_smul, smul_eq_mul, mul_comm]

/-- **Linear corollary of Lemma L.** For a valid template whose rows are multiples of
affine-linear equations vanishing at a torus point `v0`, the linear parts of the used
equations span `e_h`: the used equations determine `x_h`. -/
theorem determines_h {ι : Type*} (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n)) (hEC : E ⊆ C) (h : Fin n) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (hEne : E.Nonempty)
    (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0)
    (lin : ι → (Fin n → ℚ)) (cst : ι → ℚ)
    (heval : ∀ r (ξ : Fin n → ℚ), (∀ j, ξ j ≠ 0) →
        ∑ α ∈ (p r).supp, (p r).coeff α * mpow ξ α = cst r + ∑ j, lin r j * ξ j)
    (v0 : Fin n → ℚ) (hv0 : ∀ j, v0 j ≠ 0) (hv0z : ∀ r, cst r + ∑ j, lin r j * v0 j = 0) :
    (Pi.single h (1 : ℚ) : Fin n → ℚ) ∈ Submodule.span ℚ (Set.range lin) := by
  by_contra hnot
  obtain ⟨φ, hφh, hφS⟩ := Submodule.exists_dual_map_eq_bot_of_nmem hnot inferInstance
  have hφlin : ∀ r, φ (lin r) = 0 := by
    intro r
    have : φ (lin r) ∈ (Submodule.span ℚ (Set.range lin)).map φ :=
      Submodule.mem_map_of_mem (Submodule.subset_span ⟨r, rfl⟩)
    rw [hφS] at this
    exact (Submodule.mem_bot ℚ).mp this
  set eh : Fin n → ℚ := Pi.single h (1 : ℚ) with heh
  set d : Fin n → ℚ := fun j => φ (Pi.single j (1 : ℚ) : Fin n → ℚ) / φ eh with hd
  have hd_lin : ∀ r, ∑ j, lin r j * d j = 0 := by
    intro r
    have h1 := hφlin r
    rw [dual_expand] at h1
    have : ∑ j, lin r j * d j = (∑ j, φ (Pi.single j (1 : ℚ) : Fin n → ℚ) * lin r j) / φ eh := by
      rw [Finset.sum_div]
      refine Finset.sum_congr rfl (fun j _ => ?_)
      simp only [hd]; ring
    rw [this, h1, zero_div]
  have hdh : d h = 1 := by simp only [hd]; exact div_self hφh
  -- values of the line parameter that would hit a coordinate hyperplane
  set B : Finset ℚ := (Finset.univ.filter fun j => d j ≠ 0).image (fun j => -v0 j / d j) with hB
  have hBc : B.card ≤ n := by
    refine Finset.card_image_le.trans ?_
    exact (Finset.card_filter_le _ _).trans (by simp)
  set R0 : Finset ℕ := Finset.range (E.card + 1 + n) with hR0
  set S : Finset ℕ := R0.filter (fun s => (s : ℚ) ∉ B) with hS
  have hS : E.card + 1 ≤ S.card := by
    have h1 := Finset.filter_card_add_filter_neg_card_eq_card (s := R0) (fun s : ℕ => (s : ℚ) ∉ B)
    have h2 : (R0.filter (fun s : ℕ => ¬ ((s : ℚ) ∉ B))).card ≤ B.card := by
      apply Finset.card_le_card_of_injOn (fun s : ℕ => (s : ℚ))
      · intro s hs
        simp only [Finset.mem_filter, not_not] at hs
        exact hs.2
      · intro a _ b _ hab
        exact Nat.cast_injective hab
    have hR0c : R0.card = E.card + 1 + n := Finset.card_range _
    have hSdef : S.card = (R0.filter (fun s : ℕ => (s : ℚ) ∉ B)).card := rfl
    omega
  obtain ⟨S', hS'S, hS'c⟩ := Finset.exists_subset_card_eq hS
  set f := S'.orderEmbOfFin hS'c with hf
  set ξ : Fin (E.card + 1) → Fin n → ℚ := fun k j => v0 j + ((f k : ℕ) : ℚ) * d j with hξdef
  have hfS : ∀ k, (((f k : ℕ)) : ℚ) ∉ B := by
    intro k
    have hmem : f k ∈ S := hS'S (Finset.orderEmbOfFin_mem S' hS'c k)
    exact (Finset.mem_filter.mp hmem).2
  have hξ : ∀ k j, ξ k j ≠ 0 := by
    intro k j h0
    by_cases hdj : d j = 0
    · simp only [hξdef, hdj, mul_zero, add_zero] at h0
      exact hv0 j h0
    · apply hfS k
      rw [hB, Finset.mem_image]
      refine ⟨j, by simp [hdj], ?_⟩
      simp only [hξdef] at h0
      field_simp
      linarith
  have hzero : ∀ k r, ∑ α ∈ (p r).supp, (p r).coeff α * mpow (ξ k) α = 0 := by
    intro k r
    rw [heval r (ξ k) (hξ k)]
    have : cst r + ∑ j, lin r j * ξ k j
        = (cst r + ∑ j, lin r j * v0 j) + ((f k : ℕ) : ℚ) * ∑ j, lin r j * d j := by
      simp only [hξdef, mul_add, Finset.sum_add_distrib, Finset.mul_sum]
      rw [add_assoc]
      congr 1
      congr 1
      refine Finset.sum_congr rfl (fun j _ => by ring)
    rw [this, hv0z, hd_lin, mul_zero, add_zero]
  have hinj : Function.Injective fun k => ξ k h := by
    intro a b hab
    simp only [hξdef, hdh, mul_one, add_right_inj] at hab
    exact f.injective (by exact_mod_cast hab)
  have := card_le t p C E hEC h hE hEne hsupp hvalid ξ hξ hzero hinj
  omega

end LemmaL

/-! ## The gadget lower bound (Theorem C2, algebraic core) -/

section Gadget

variable {n : ℕ} {ε : Type*} {ι : Type*} [Fintype ι]

/-- Linear part of the equations: `none` is `f0 = h - ∑_{x∈X} x`, `some e` is the gadget row. -/
def linOf (h : Fin n) (X : Finset (Fin n)) (yv xv : ε → Fin n) (r : ε → ℚ) :
    Option ε → (Fin n → ℚ)
  | none => (Pi.single h (1 : ℚ) : Fin n → ℚ) - ∑ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ)
  | some e => gvec yv xv r e

/-- If the linear parts of the used equations span `e_h`, then `f0` is used, the set
variables of the used gadgets cover `X`, and there are at least `1 + |X| + |sets used|`
rows. -/
theorem gadget_rows_bound (h : Fin n) (X : Finset (Fin n)) (hhX : h ∉ X)
    (yv xv : ε → Fin n) (r : ε → ℚ) (hr : ∀ e, 0 < r e) (hx : ∀ e, xv e ∈ X)
    (hy : ∀ e, yv e ∉ X) (hyh : ∀ e, yv e ≠ h) (eq : ι → Option ε)
    (hspan : (Pi.single h (1 : ℚ) : Fin n → ℚ) ∈
      Submodule.span ℚ (Set.range fun i => linOf h X yv xv r (eq i))) :
    let Ug : Finset ε := Finset.univ.biUnion fun i => (eq i).toFinset
    X ⊆ Ug.image xv ∧ 1 + X.card + (Ug.image yv).card ≤ Fintype.card ι := by
  intro Ug
  set lin0 : Fin n → ℚ := linOf h X yv xv r none with hlin0
  set Sg : Set (Fin n → ℚ) := ((Ug.image (gvec yv xv r) : Finset (Fin n → ℚ)) : Set (Fin n → ℚ))
  have hsub : Set.range (fun i => linOf h X yv xv r (eq i)) ⊆ insert lin0 Sg := by
    rintro _ ⟨i, rfl⟩
    cases hi : eq i with
    | none => simp only [hi]; exact Set.mem_insert _ _
    | some e =>
      simp only [hi]
      apply Set.mem_insert_of_mem
      simp only [Sg, Finset.coe_image, Set.mem_image, Finset.mem_coe]
      refine ⟨e, ?_, rfl⟩
      simp only [Ug, Finset.mem_biUnion, Finset.mem_univ, true_and]
      exact ⟨i, by simp [hi]⟩
  have hspan' := Submodule.span_mono hsub hspan
  obtain ⟨a, z, hz, hza⟩ := Submodule.mem_span_insert.mp hspan'
  -- evaluation at h kills the gadget span
  have evh : ∀ w ∈ Submodule.span ℚ Sg, w h = 0 := by
    intro w hw
    have hle : Submodule.span ℚ Sg ≤ LinearMap.ker (LinearMap.proj (R := ℚ) (φ := fun _ : Fin n => ℚ) h) := by
      rw [Submodule.span_le]
      intro g hg
      simp only [Sg, Finset.coe_image, Set.mem_image, Finset.mem_coe] at hg
      obtain ⟨e, _, rfl⟩ := hg
      have hxh : xv e ≠ h := fun hh => hhX (hh ▸ hx e)
      simp [gvec, Pi.single_apply, hyh e, hxh]
    exact LinearMap.mem_ker.mp (hle hw)
  have hl0h : lin0 h = 1 := by
    have : ∀ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) h = 0 := fun x hxX => by
      have hne : h ≠ x := fun hh => hhX (hh ▸ hxX)
      simp [Pi.single_apply, hne]
    simp [hlin0, linOf, Finset.sum_apply, Finset.sum_eq_zero this]
  have ha : a = 1 := by
    have := congrFun hza h
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, hl0h, mul_one, evh z hz,
      add_zero, Pi.single_eq_same] at this
    exact this.symm
  have hzs : z = ∑ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) := by
    rw [ha, one_smul] at hza
    have : z = (Pi.single h (1 : ℚ) : Fin n → ℚ) - lin0 := by rw [hza]; abel
    rw [this, hlin0]
    simp [linOf]
  rw [hzs] at hz
  obtain ⟨hcov, hcount⟩ := gadget_count X yv xv r hr hx hy Ug hz
  refine ⟨hcov, ?_⟩
  -- f0 is used, and the used equations inject into the rows
  have hf0 : ∃ i, eq i = none := by
    by_contra hno
    push_neg at hno
    have hsub2 : Set.range (fun i => linOf h X yv xv r (eq i)) ⊆ Sg := by
      rintro _ ⟨i, rfl⟩
      obtain ⟨e, he⟩ := Option.ne_none_iff_exists'.mp (hno i)
      simp only [he, linOf, Sg, Finset.coe_image, Set.mem_image, Finset.mem_coe]
      refine ⟨e, ?_, rfl⟩
      simp only [Ug, Finset.mem_biUnion, Finset.mem_univ, true_and]
      exact ⟨i, by simp [he]⟩
    have := evh _ (Submodule.span_mono hsub2 hspan)
    simp at this
  have himg : insert none (Ug.map ⟨some, Option.some_injective _⟩) ⊆ Finset.univ.image eq := by
    intro o ho
    rcases Finset.mem_insert.mp ho with rfl | ho
    · obtain ⟨i, hi⟩ := hf0
      exact Finset.mem_image.mpr ⟨i, Finset.mem_univ _, hi⟩
    · obtain ⟨e, he, rfl⟩ := Finset.mem_map.mp ho
      simp only [Ug, Finset.mem_biUnion, Finset.mem_univ, true_and, Option.mem_toFinset,
        Option.mem_def] at he
      obtain ⟨i, hi⟩ := he
      exact Finset.mem_image.mpr ⟨i, Finset.mem_univ _, hi⟩
  have hc1 := Finset.card_le_card himg
  rw [Finset.card_insert_of_not_mem (by simp), Finset.card_map] at hc1
  have hc2 : (Finset.univ.image eq).card ≤ Fintype.card ι := Finset.card_image_le.trans (by simp)
  omega

end Gadget

/-! ## Squarification and the final lower bound -/

section Final

open FaceObstruction

universe u

variable {n : ℕ} {ε : Type*}

lemma dual_expand_ft {κ : Type*} [Fintype κ] [DecidableEq κ] (f : (κ → ℚ) →ₗ[ℚ] ℚ) (v : κ → ℚ) :
    f v = ∑ g, f (Pi.single g (1 : ℚ) : κ → ℚ) * v g := by
  have hv : v = ∑ g, v g • (Pi.single g (1 : ℚ) : κ → ℚ) := by
    funext i; simp [Finset.sum_apply, Pi.single_apply]
  conv_lhs => rw [hv]
  rw [map_sum]
  refine Finset.sum_congr rfl (fun g _ => ?_)
  rw [map_smul, smul_eq_mul, mul_comm]

/-- Sums over the columns of a vector supported in `G = C \ E`. -/
lemma sum_supp {ι : Type*} (t : ι → Exp n) (p : ι → SPoly ℚ n) (C E : Finset (Exp n))
    (U : Exp n → ℚ) (hU : ∀ b, U b ≠ 0 → b ∈ C ∧ b ∉ E) (i : ι) :
    ∑ b ∈ C, entry t p i b * U b = ∑ g : ((C \ E : Finset (Exp n))), entry t p i g * U g := by
  rw [Finset.sum_coe_sort (C \ E) (fun b => entry t p i b * U b)]
  refine (Finset.sum_subset Finset.sdiff_subset (fun b hbC hbG => ?_)).symm
  have hbE : b ∈ E := by
    by_contra hbE; exact hbG (Finset.mem_sdiff.mpr ⟨hbC, hbE⟩)
  have : U b = 0 := by
    by_contra hne; exact (hU b hne).2 hbE
  rw [this, mul_zero]

/-- **Squarification.** A valid template has a valid sub-template with exactly `|C \ E|`
rows (the rows restricted to `G` span, so a basis can be extracted). -/
theorem squarify {ι : Type u} [Fintype ι] (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n))
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0) :
    ∃ (κ : Type u) (_ : Fintype κ) (a : κ → ι), Fintype.card κ = (C \ E).card ∧
      ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ k, ∑ b ∈ C, entry (t ∘ a) (p ∘ a) k b * u b = 0) → ∀ b, u b = 0 := by
  set ρ : ι → ((C \ E : Finset (Exp n)) → ℚ) := fun i g => entry t p i g with hρ
  -- extending a vector on G by zero
  set ext : ((C \ E : Finset (Exp n)) → ℚ) → Exp n → ℚ := fun u b => if hb : b ∈ C \ E then u ⟨b, hb⟩ else 0 with hext
  have ext_supp : ∀ u b, ext u b ≠ 0 → b ∈ C ∧ b ∉ E := by
    intro u b hb
    by_cases hbG : b ∈ C \ E
    · exact Finset.mem_sdiff.mp hbG
    · exact absurd (show ext u b = 0 from dif_neg hbG) hb
  have ext_sum : ∀ (κ : Type u) (a : κ → ι) (u : (C \ E : Finset (Exp n)) → ℚ) (k : κ),
      ∑ b ∈ C, entry (t ∘ a) (p ∘ a) k b * ext u b = ∑ g, ρ (a k) g * u g := by
    intro κ a u k
    rw [sum_supp (t ∘ a) (p ∘ a) C E (ext u) (ext_supp u) k]
    refine Finset.sum_congr rfl (fun g _ => ?_)
    have hg : (g : Exp n) ∈ C \ E := g.2
    simp only [hext, hρ, dif_pos hg]
    rfl
  -- the restricted rows span everything
  have hspan : Submodule.span ℚ (Set.range ρ) = ⊤ := by
    apply Classical.byContradiction
    intro hne
    obtain ⟨f, hf0, hfle⟩ := Submodule.exists_le_ker_of_lt_top _ (lt_top_iff_ne_top.mpr hne)
    set u : (C \ E : Finset (Exp n)) → ℚ := fun g => f (Pi.single g (1 : ℚ) : (C \ E : Finset (Exp n)) → ℚ) with hu
    have hrow : ∀ i, ∑ b ∈ C, entry t p i b * ext u b = 0 := by
      intro i
      have h1 := ext_sum ι id u i
      simp only [Function.comp_id, id] at h1
      rw [h1]
      have h2 : f (ρ i) = 0 := LinearMap.mem_ker.mp (hfle (Submodule.subset_span ⟨i, rfl⟩))
      rw [dual_expand_ft] at h2
      rw [← h2]
      refine Finset.sum_congr rfl (fun g _ => by rw [hu, mul_comm])
    have hz := hvalid (ext u) (ext_supp u) hrow
    apply hf0
    refine LinearMap.ext (fun v => ?_)
    rw [dual_expand_ft, LinearMap.zero_apply]
    refine Finset.sum_eq_zero (fun g _ => ?_)
    have : u g = 0 := by
      have h1 := hz g
      have h2 : ext u g = u g := dif_pos g.2
      rwa [h2] at h1
    rw [show f (Pi.single g (1 : ℚ) : (C \ E : Finset (Exp n)) → ℚ) = u g from rfl, this, zero_mul]
  -- extract a basis of rows
  obtain ⟨κ, a, ha_inj, ha_span, ha_li⟩ := exists_linearIndependent' ℚ ρ
  haveI : Finite κ := Finite.of_injective a ha_inj
  letI : Fintype κ := Fintype.ofFinite κ
  have hcard : Fintype.card κ = (C \ E).card := by
    have h1 := finrank_span_eq_card ha_li
    rw [ha_span, hspan, finrank_top, Module.finrank_fintype_fun_eq_card, Fintype.card_coe] at h1
    exact h1.symm
  refine ⟨κ, inferInstance, a, hcard, ?_⟩
  -- validity of the sub-template
  intro U hU hker
  set φ : ((C \ E : Finset (Exp n)) → ℚ) →ₗ[ℚ] ℚ := ∑ g : (C \ E : Finset (Exp n)), U g • LinearMap.proj g with hφ
  have hφapp : ∀ v : (C \ E : Finset (Exp n)) → ℚ, φ v = ∑ g : (C \ E : Finset (Exp n)), U g * v g := by intro v; simp [hφ]
  have hle : Submodule.span ℚ (Set.range (ρ ∘ a)) ≤ LinearMap.ker φ := by
    rw [Submodule.span_le]
    rintro _ ⟨k, rfl⟩
    simp only [SetLike.mem_coe, LinearMap.mem_ker, hφapp, Function.comp]
    have := hker k
    rw [sum_supp (t ∘ a) (p ∘ a) C E U hU k] at this
    rw [← this]
    refine Finset.sum_congr rfl (fun g _ => ?_)
    rw [mul_comm]
    rfl
  rw [ha_span, hspan] at hle
  intro b
  by_cases hbG : b ∈ C \ E
  · have h1 := LinearMap.mem_ker.mp (hle (Submodule.mem_top (x := (Pi.single ⟨b, hbG⟩ (1 : ℚ) : (C \ E : Finset (Exp n)) → ℚ))))
    rw [hφapp] at h1
    simpa [Pi.single_apply] using h1
  · by_contra hne
    exact hbG (Finset.mem_sdiff.mpr (hU b hne))

/-- **Theorem C2, algebraic core (machine-checked).** For the SET COVER gadget system
(f0 = h − Σ_{x∈X} x; gadget rows y_{yv e} + r_e x_{xv e} with r_e > 0; a common torus
solution v0), every valid template satisfies
`|C \ E| ≥ 1 + |X| + |set variables used|`, and the used gadgets cover X. With the flat
template achieving `1 + |X| + OPT`, the minimum elimination size is `1 + |X| + OPT`. -/
theorem min_elim_lower_bound {ι : Type u} [Fintype ι] (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n)) (hEC : E ⊆ C) (h : Fin n) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (hEne : E.Nonempty) (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0)
    (X : Finset (Fin n)) (hhX : h ∉ X) (yv xv : ε → Fin n) (r : ε → ℚ) (hr : ∀ e, 0 < r e)
    (hx : ∀ e, xv e ∈ X) (hy : ∀ e, yv e ∉ X) (hyh : ∀ e, yv e ≠ h)
    (eq : ι → Option ε) (cst : ι → ℚ)
    (heval : ∀ i (ξ : Fin n → ℚ), (∀ j, ξ j ≠ 0) → ∑ α ∈ (p i).supp, (p i).coeff α * mpow ξ α =
        cst i + ∑ j, linOf h X yv xv r (eq i) j * ξ j)
    (v0 : Fin n → ℚ) (hv0 : ∀ j, v0 j ≠ 0)
    (hv0z : ∀ i, cst i + ∑ j, linOf h X yv xv r (eq i) j * v0 j = 0) :
    ∃ Ug : Finset ε, X ⊆ Ug.image xv ∧ 1 + X.card + (Ug.image yv).card ≤ (C \ E).card := by
  obtain ⟨κ, _, a, hcard, hvalid'⟩ := squarify t p C E hvalid
  have hdet := LemmaL.determines_h (t ∘ a) (p ∘ a) C E hEC h hE hEne
    (fun k α hα => hsupp (a k) α hα) hvalid' (fun k => linOf h X yv xv r (eq (a k)))
    (cst ∘ a) (fun k ξ hξ => heval (a k) ξ hξ) v0 hv0 (fun k => hv0z (a k))
  obtain ⟨hcov, hb⟩ := gadget_rows_bound h X hhX yv xv r hr hx hy hyh (eq ∘ a) hdet
  exact ⟨_, hcov, hcard ▸ hb⟩

end Final

/-! ## Achievability: covers give spanning flat templates -/

section Achieve

variable {n : ℕ} {ε : Type*}

/-- **Achievability (linear-algebra core of Theorem C1).** If the chosen gadgets `Ug` cover
`X`, and each chosen gadget's set variable carries a "double edge" (two chosen gadgets with
the same set and element variables but different gains), then the linear parts of `f0`
and of the chosen gadgets span the unit vector of every variable in `{h} ∪ X ∪ Y_used`.
So the flat template on these variables has full rank, i.e. it is valid. -/
theorem cover_spans (h : Fin n) (X : Finset (Fin n)) (yv xv : ε → Fin n) (r : ε → ℚ)
    (hr : ∀ e, 0 < r e) (Ug : Finset ε)
    (hcov : ∀ x ∈ X, ∃ e ∈ Ug, xv e = x)
    (hdouble : ∀ e ∈ Ug, ∃ e1 ∈ Ug, ∃ e2 ∈ Ug, yv e1 = yv e ∧ yv e2 = yv e ∧
        xv e1 = xv e2 ∧ r e1 ≠ r e2) :
    let S := Submodule.span ℚ (insert (linOf h X yv xv r none)
      ((Ug.image (gvec yv xv r) : Finset (Fin n → ℚ)) : Set (Fin n → ℚ)))
    (Pi.single h (1 : ℚ) : Fin n → ℚ) ∈ S ∧ (∀ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) ∈ S) ∧
      (∀ e ∈ Ug, (Pi.single (yv e) (1 : ℚ) : Fin n → ℚ) ∈ S) := by
  intro S
  have hg : ∀ e ∈ Ug, gvec yv xv r e ∈ S := fun e he =>
    Submodule.subset_span (Set.mem_insert_of_mem _
      (Finset.mem_coe.mpr (Finset.mem_image.mpr ⟨e, he, rfl⟩)))
  have hl0 : linOf h X yv xv r none ∈ S := Submodule.subset_span (Set.mem_insert _ _)
  -- set variables
  have hy : ∀ e ∈ Ug, (Pi.single (yv e) (1 : ℚ) : Fin n → ℚ) ∈ S := by
    intro e he
    obtain ⟨e1, he1, e2, he2, hy1, hy2, hx12, hr12⟩ := hdouble e he
    have hsub : gvec yv xv r e1 - gvec yv xv r e2 =
        (r e1 - r e2) • (Pi.single (xv e1) (1 : ℚ) : Fin n → ℚ) := by
      simp only [gvec, hy1, hy2, hx12, sub_smul]; abel
    have hxe1 : (Pi.single (xv e1) (1 : ℚ) : Fin n → ℚ) ∈ S := by
      have hm : (r e1 - r e2) • (Pi.single (xv e1) (1 : ℚ) : Fin n → ℚ) ∈ S :=
        hsub ▸ S.sub_mem (hg e1 he1) (hg e2 he2)
      have hne : r e1 - r e2 ≠ 0 := sub_ne_zero.mpr hr12
      have := S.smul_mem (r e1 - r e2)⁻¹ hm
      rwa [smul_smul, inv_mul_cancel₀ hne, one_smul] at this
    have : (Pi.single (yv e) (1 : ℚ) : Fin n → ℚ) =
        gvec yv xv r e1 - r e1 • (Pi.single (xv e1) (1 : ℚ) : Fin n → ℚ) := by
      simp only [gvec, hy1]; abel
    rw [this]
    exact S.sub_mem (hg e1 he1) (S.smul_mem _ hxe1)
  -- element variables
  have hx : ∀ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) ∈ S := by
    intro x hxX
    obtain ⟨e, he, rfl⟩ := hcov x hxX
    have : (Pi.single (xv e) (1 : ℚ) : Fin n → ℚ) =
        (r e)⁻¹ • (gvec yv xv r e - (Pi.single (yv e) (1 : ℚ) : Fin n → ℚ)) := by
      simp only [gvec, add_sub_cancel_left, smul_smul, inv_mul_cancel₀ (hr e).ne', one_smul]
    rw [this]
    exact S.smul_mem _ (S.sub_mem (hg e he) (hy e he))
  -- the hidden variable
  refine ⟨?_, hx, hy⟩
  have : (Pi.single h (1 : ℚ) : Fin n → ℚ) =
      linOf h X yv xv r none + ∑ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) := by
    simp [linOf]
  rw [this]
  exact S.add_mem hl0 (S.sum_mem (fun x hxX => hx x hxX))

end Achieve

/-! ## Theorem C3: any number N of solutions (nonlinear Lemma L step) -/

section C3

open FaceObstruction

universe u

variable {n : ℕ} {ε : Type*}

/-- `m` distinct naturals avoiding a finite set of rationals. -/
lemma good_naturals (B : Finset ℚ) (m : ℕ) :
    ∃ f : Fin m → ℕ, Function.Injective f ∧ ∀ k, ((f k : ℕ) : ℚ) ∉ B := by
  set R0 : Finset ℕ := Finset.range (m + B.card) with hR0
  set S : Finset ℕ := R0.filter (fun s => (s : ℚ) ∉ B) with hS
  have hSc : m ≤ S.card := by
    have h1 := Finset.filter_card_add_filter_neg_card_eq_card (s := R0) (fun s : ℕ => (s : ℚ) ∉ B)
    have h2 : (R0.filter (fun s : ℕ => ¬ ((s : ℚ) ∉ B))).card ≤ B.card := by
      apply Finset.card_le_card_of_injOn (fun s : ℕ => (s : ℚ))
      · intro s hs
        simp only [Finset.mem_filter, not_not] at hs
        exact hs.2
      · intro a _ b _ hab
        exact Nat.cast_injective hab
    have hR0c : R0.card = m + B.card := Finset.card_range _
    have hSdef : S.card = (R0.filter (fun s : ℕ => (s : ℚ) ∉ B)).card := rfl
    omega
  obtain ⟨S', hS'S, hS'c⟩ := Finset.exists_subset_card_eq hSc
  refine ⟨fun k => S'.orderEmbOfFin hS'c k, (S'.orderEmbOfFin hS'c).injective, fun k => ?_⟩
  exact (Finset.mem_filter.mp (hS'S (Finset.orderEmbOfFin_mem S' hS'c k))).2

/-- Evaluation of the C3 equations: `none` is `f0' = P(h) - Σ_{x∈X} x`, `some e` is a gadget. -/
def evalC3 (h : Fin n) (X : Finset (Fin n)) (yv xv : ε → Fin n) (r : ε → ℚ)
    (P : Polynomial ℚ) (cst : ε → ℚ) : Option ε → (Fin n → ℚ) → ℚ
  | none, ξ => P.eval (ξ h) - ∑ x ∈ X, ξ x
  | some e, ξ => cst e + ∑ j, gvec yv xv r e j * ξ j

/-- **Nonlinear Lemma L step (Theorem C3).** In a valid template for the C3 system
(`f0' = P(h) − Σ x` with `P` non-constant, positive-gain gadgets, torus solution `v0`),
`f0'` is used and `Σ_{x∈X} e_x` lies in the span of the used gadget rows. -/
theorem determines_sum_poly {ι : Type*} [Fintype ι] (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n)) (hEC : E ⊆ C) (h : Fin n) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (hEne : E.Nonempty) (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0)
    (X : Finset (Fin n)) (hhX : h ∉ X) (yv xv : ε → Fin n) (r : ε → ℚ)
    (hx : ∀ e, xv e ∈ X) (hyh : ∀ e, yv e ≠ h)
    (P : Polynomial ℚ) (hP : 0 < P.natDegree) (cst : ε → ℚ) (eq : ι → Option ε)
    (heval : ∀ i (ξ : Fin n → ℚ), (∀ j, ξ j ≠ 0) →
        ∑ α ∈ (p i).supp, (p i).coeff α * mpow ξ α = evalC3 h X yv xv r P cst (eq i) ξ)
    (v0 : Fin n → ℚ) (hv0 : ∀ j, v0 j ≠ 0) (hv0z : ∀ o, evalC3 h X yv xv r P cst o v0 = 0) :
    (∃ i, eq i = none) ∧
      (∑ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ)) ∈ Submodule.span ℚ
        (((Finset.univ.biUnion fun i => (eq i).toFinset).image (gvec yv xv r) : Finset (Fin n → ℚ)) :
          Set (Fin n → ℚ)) := by
  have hgh : ∀ e, gvec yv xv r e h = 0 := by
    intro e
    have hxh : xv e ≠ h := fun hh => hhX (hh ▸ hx e)
    simp [gvec, Pi.single_apply, hyh e, hxh]
  -- gadget rows ignore the h-coordinate
  have gad_h : ∀ e (ξ ξ' : Fin n → ℚ), (∀ j, j ≠ h → ξ j = ξ' j) →
      evalC3 h X yv xv r P cst (some e) ξ = evalC3 h X yv xv r P cst (some e) ξ' := by
    intro e ξ ξ' hξ
    simp only [evalC3]
    congr 1
    refine Finset.sum_congr rfl (fun j _ => ?_)
    by_cases hj : j = h
    · subst hj; rw [hgh e, zero_mul, zero_mul]
    · rw [hξ j hj]
  set Ug : Finset ε := Finset.univ.biUnion fun i => (eq i).toFinset with hUg
  have hUg_mem : ∀ i e, eq i = some e → e ∈ Ug := by
    intro i e hi
    simp only [hUg, Finset.mem_biUnion, Finset.mem_univ, true_and]
    exact ⟨i, by simp [hi]⟩
  -- a family of torus points with distinct h solving all used rows contradicts Lemma L
  have contra : ∀ (ξ : Fin (E.card + 1) → Fin n → ℚ), (∀ k j, ξ k j ≠ 0) →
      (∀ k i, evalC3 h X yv xv r P cst (eq i) (ξ k) = 0) →
      Function.Injective (fun k => ξ k h) → False := by
    intro ξ hξ hz hinj
    have := LemmaL.card_le t p C E hEC h hE hEne hsupp hvalid ξ hξ
      (fun k i => by rw [heval i (ξ k) (hξ k)]; exact hz k i) hinj
    omega
  constructor
  · -- f0' must be used
    by_contra hno
    push_neg at hno
    obtain ⟨f, hfinj, hfB⟩ := good_naturals {0} (E.card + 1)
    refine contra (fun k j => if j = h then ((f k : ℕ) : ℚ) else v0 j) ?_ ?_ ?_
    · intro k j
      by_cases hj : j = h
      · simp only [hj, if_true]; intro h0; exact hfB k (by simp [h0])
      · simp only [hj, if_false]; exact hv0 j
    · intro k i
      obtain ⟨e, he⟩ := Option.ne_none_iff_exists'.mp (hno i)
      rw [he, gad_h e _ v0 (fun j hj => by simp [hj])]
      exact hv0z (some e)
    · intro a b hab
      simp only [if_true] at hab
      exact hfinj (by exact_mod_cast hab)
  · -- Σ e_x is in the gadget span
    by_contra hnot
    obtain ⟨φ, hφs, hφS⟩ := Submodule.exists_dual_map_eq_bot_of_nmem hnot inferInstance
    set sX : Fin n → ℚ := ∑ x ∈ X, (Pi.single x (1 : ℚ) : Fin n → ℚ) with hsX
    have hφg : ∀ e ∈ Ug, φ (gvec yv xv r e) = 0 := by
      intro e he
      have hmem : gvec yv xv r e ∈ Submodule.span ℚ ((Ug.image (gvec yv xv r) : Finset _) : Set _) :=
        Submodule.subset_span (Finset.mem_coe.mpr (Finset.mem_image.mpr ⟨e, he, rfl⟩))
      have : φ (gvec yv xv r e) ∈ (Submodule.span ℚ ((Ug.image (gvec yv xv r) : Finset _) :
          Set _)).map φ := Submodule.mem_map_of_mem hmem
      rw [hφS] at this
      exact (Submodule.mem_bot ℚ).mp this
    set d : Fin n → ℚ := fun j => if j = h then 0 else φ (Pi.single j (1 : ℚ) : Fin n → ℚ) / φ sX
      with hd
    have hdot : ∀ w : Fin n → ℚ, w h = 0 → ∑ j, w j * d j = φ w / φ sX := by
      intro w hw
      rw [LemmaL.dual_expand φ w, Finset.sum_div]
      refine Finset.sum_congr rfl (fun j _ => ?_)
      by_cases hj : j = h
      · subst hj; simp [hd, hw]
      · simp only [hd, if_neg hj]; ring
    have hsXh : sX h = 0 := by
      simp only [hsX, Finset.sum_apply]
      exact Finset.sum_eq_zero (fun x hxX => by
        have : h ≠ x := fun hh => hhX (hh ▸ hxX)
        simp [Pi.single_apply, this])
    have hdX : ∑ x ∈ X, d x = 1 := by
      have h1 := hdot sX hsXh
      rw [div_self hφs] at h1
      rw [← h1, hsX]
      simp only [Finset.sum_apply, Finset.sum_mul]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl (fun x _ => ?_)
      simp [Pi.single_apply]
    have hdg : ∀ e ∈ Ug, ∑ j, gvec yv xv r e j * d j = 0 := by
      intro e he
      rw [hdot _ (hgh e), hφg e he, zero_div]
    -- bad parameter values
    set sP : ℚ → ℚ := fun τ => P.eval τ - P.eval (v0 h) with hsP
    have hroots : ∀ c : ℚ, P - Polynomial.C c ≠ 0 := by
      intro c hc
      have := congrArg Polynomial.natDegree hc
      rw [Polynomial.natDegree_sub_C, Polynomial.natDegree_zero] at this
      omega
    set B : Finset ℚ := insert 0 ((Finset.univ.filter fun j : Fin n => j ≠ h ∧ d j ≠ 0).biUnion
      fun j => (P - Polynomial.C (P.eval (v0 h) - v0 j / d j)).roots.toFinset) with hB
    obtain ⟨f, hfinj, hfB⟩ := good_naturals B (E.card + 1)
    set ξ : Fin (E.card + 1) → Fin n → ℚ := fun k j =>
      if j = h then ((f k : ℕ) : ℚ) else v0 j + sP ((f k : ℕ) : ℚ) * d j with hξ
    refine contra ξ ?_ ?_ ?_
    · intro k j hzero
      by_cases hj : j = h
      · simp only [hξ, hj, if_true] at hzero
        exact hfB k (by rw [hB, hzero]; exact Finset.mem_insert_self _ _)
      · simp only [hξ, hj, if_false] at hzero
        by_cases hdj : d j = 0
        · rw [hdj, mul_zero, add_zero] at hzero; exact hv0 j hzero
        · apply hfB k
          rw [hB]
          apply Finset.mem_insert_of_mem
          rw [Finset.mem_biUnion]
          refine ⟨j, by simp [hj, hdj], ?_⟩
          rw [Multiset.mem_toFinset, Polynomial.mem_roots (hroots _), Polynomial.IsRoot,
            Polynomial.eval_sub, Polynomial.eval_C]
          simp only [hsP] at hzero
          field_simp
          linarith
    · intro k i
      cases hi : eq i with
      | none =>
        have h0 := hv0z none
        simp only [evalC3] at h0 ⊢
        have hsum : ∑ x ∈ X, ξ k x = ∑ x ∈ X, v0 x + sP ((f k : ℕ) : ℚ) := by
          have : ∀ x ∈ X, ξ k x = v0 x + sP ((f k : ℕ) : ℚ) * d x := by
            intro x hxX
            have : x ≠ h := fun hh => hhX (hh ▸ hxX)
            simp [hξ, this]
          rw [Finset.sum_congr rfl this, Finset.sum_add_distrib, ← Finset.mul_sum, hdX, mul_one]
        rw [hsum]
        simp only [hξ, if_true, hsP]
        linarith
      | some e =>
        have he := hUg_mem i e hi
        have h0 := hv0z (some e)
        simp only [evalC3] at h0 ⊢
        have : ∑ j, gvec yv xv r e j * ξ k j =
            ∑ j, gvec yv xv r e j * v0 j + sP ((f k : ℕ) : ℚ) * ∑ j, gvec yv xv r e j * d j := by
          rw [Finset.mul_sum, ← Finset.sum_add_distrib]
          refine Finset.sum_congr rfl (fun j _ => ?_)
          by_cases hj : j = h
          · subst hj; simp [hgh e]
          · simp only [hξ, if_neg hj]; ring
        rw [this, hdg e he, mul_zero, add_zero]
        exact h0
    · intro a b hab
      simp only [hξ, if_true] at hab
      exact hfinj (by exact_mod_cast hab)

/-- **Theorem C3 lower bound (machine-checked).** For the C3 system with any non-constant
`P` (degree N, so N solutions), every valid template has
`|C \ E| ≥ 1 + |X| + |set variables used|`, and the used gadgets cover `X`. -/
theorem min_elim_lower_bound_poly {ι : Type u} [Fintype ι] (t : ι → Exp n) (p : ι → SPoly ℚ n)
    (C E : Finset (Exp n)) (hEC : E ⊆ C) (h : Fin n) (hE : ∀ e ∈ E, e + Pi.single h 1 ∈ C)
    (hEne : E.Nonempty) (hsupp : ∀ r, ∀ α ∈ (p r).supp, t r + α ∈ C)
    (hvalid : ∀ u : Exp n → ℚ, (∀ b, u b ≠ 0 → b ∈ C ∧ b ∉ E) →
        (∀ r, ∑ b ∈ C, entry t p r b * u b = 0) → ∀ b, u b = 0)
    (X : Finset (Fin n)) (hhX : h ∉ X) (yv xv : ε → Fin n) (r : ε → ℚ) (hr : ∀ e, 0 < r e)
    (hx : ∀ e, xv e ∈ X) (hy : ∀ e, yv e ∉ X) (hyh : ∀ e, yv e ≠ h)
    (P : Polynomial ℚ) (hP : 0 < P.natDegree) (cst : ε → ℚ) (eq : ι → Option ε)
    (heval : ∀ i (ξ : Fin n → ℚ), (∀ j, ξ j ≠ 0) →
        ∑ α ∈ (p i).supp, (p i).coeff α * mpow ξ α = evalC3 h X yv xv r P cst (eq i) ξ)
    (v0 : Fin n → ℚ) (hv0 : ∀ j, v0 j ≠ 0) (hv0z : ∀ o, evalC3 h X yv xv r P cst o v0 = 0) :
    ∃ Ug : Finset ε, X ⊆ Ug.image xv ∧ 1 + X.card + (Ug.image yv).card ≤ (C \ E).card := by
  obtain ⟨κ, _, a, hcard, hvalid'⟩ := squarify t p C E hvalid
  obtain ⟨hf0, hz⟩ := determines_sum_poly (t ∘ a) (p ∘ a) C E hEC h hE hEne
    (fun k α hα => hsupp (a k) α hα) hvalid' X hhX yv xv r hx hyh P hP cst (eq ∘ a)
    (fun k ξ hξ => heval (a k) ξ hξ) v0 hv0 hv0z
  set Ug : Finset ε := Finset.univ.biUnion fun k => ((eq ∘ a) k).toFinset with hUg
  obtain ⟨hcov, hcount⟩ := gadget_count X yv xv r hr hx hy Ug hz
  refine ⟨Ug, hcov, ?_⟩
  have himg : insert none (Ug.map ⟨some, Option.some_injective _⟩) ⊆ Finset.univ.image (eq ∘ a) := by
    intro o ho
    rcases Finset.mem_insert.mp ho with rfl | ho
    · obtain ⟨k, hk⟩ := hf0
      exact Finset.mem_image.mpr ⟨k, Finset.mem_univ _, hk⟩
    · obtain ⟨e, he, rfl⟩ := Finset.mem_map.mp ho
      simp only [hUg, Finset.mem_biUnion, Finset.mem_univ, true_and, Option.mem_toFinset,
        Option.mem_def] at he
      obtain ⟨k, hk⟩ := he
      exact Finset.mem_image.mpr ⟨k, Finset.mem_univ _, hk⟩
  have hc1 := Finset.card_le_card himg
  rw [Finset.card_insert_of_not_mem (by simp), Finset.card_map] at hc1
  have hc2 : (Finset.univ.image (eq ∘ a)).card ≤ Fintype.card κ := Finset.card_image_le.trans (by simp)
  omega

end C3

end Hardness
