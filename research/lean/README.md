# Lean formalisation (machine-checked)

- `FaceObstruction.lean`: the key lemma and Theorems 1–2 of `../proof/face_obstruction.md`.
- `EigenSize.lean`: Step 3 of `../proof/exact_count.md`. For a feasible template, the
  eigen-basis reduction reaches exactly the corank, and every valid template has
  an eigen basis at least that large.
- `Theorem3.lean`: assembles Theorem 3 (toric eigen size, simple points) in
  `toric_eigen_size`. Two inputs are explicit hypotheses:
  - `H_span` (**R3**, the Serre-vanishing step, deliberately left as a hypothesis);
  - `H_sep` (**R1/R2**, independence of the point vectors on the shiftable columns;
    elementary for large dilations kQ, and the next milestone).

- `Hardness.lean`: Theorem C2 of `../proof/complexity.md` (NP-hardness of the
  minimal elimination template), its algebraic core:
  - `gadget_count`: the sign-trick lemma. Positive-gain gadget rows whose span contains
    Σ x span every variable they touch.
  - `LemmaL.card_le`: Lemma L for arbitrary polynomial templates. N torus zeros of the
    used rows with distinct x_h force N ≤ |E| (Vandermonde argument).
  - `LemmaL.determines_h`: for affine-linear equations, the used equations determine x_h.
  - `squarify`: every valid template has a valid sub-template with exactly |C \ E| rows.
  - `min_elim_lower_bound`: for the SET COVER gadget system, every valid template has
    |C \ E| ≥ 1 + |X| + |set variables used|, and the used gadgets cover X.
  - `cover_spans`: achievability. A cover with double edges spans every variable of the
    flat template, so the flat template is valid.

  Not formalised: the encoding of the flat template as a `FaceObstruction` template
  (routine), and the complexity-theoretic wrapper (NP membership and the Karp reduction).
  Mathlib has no framework for that.

## Formalisation status of Theorem 3

| piece | status |
|---|---|
| Theorem 1, obstruction (`no_valid_template`) | proved |
| Theorem 2, spurious eigenpairs (`spurious_eigenpair`) | proved |
| point vectors lie in the kernel (`face_kernel`, `point_mem_ker`) | proved |
| Step 3, minimal eigen size = corank (`min_eigen_size_eq_corank`) | proved |
| assembly: corank = feasible eigen size = \|Z\| (`toric_eigen_size`) | proved, from `H_span` and `H_sep` |
| R1/R2, `H_sep` for kQ with k large and x_h regular on Z | to do (elementary) |
| R3, `H_span` (Serre vanishing) | hypothesis by design |
| points of multiplicity > 1 | not covered (simple points only) |

It builds with Lean 4.19.0 and Mathlib v4.19.0. `#print axioms` for
`face_kernel`, `no_valid_template` and `spurious_eigenpair` reports only the
standard axioms propext, Classical.choice and Quot.sound, so there is no
`sorry`. To rebuild:

```bash
cd research/lean
lake exe cache get
lake build
```

The first command downloads the prebuilt Mathlib cache (about 5.6 GB unpacked).
