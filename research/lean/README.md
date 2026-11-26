# Lean formalisation (machine-checked)

- `FaceObstruction.lean`: the key lemma and Theorems 1–2 of `../proof/face_obstruction.md`.
- `EigenSize.lean`: Step 3 of `../proof/exact_count.md`. For a feasible template, the
  eigen-basis reduction reaches exactly the corank, and every valid template has
  an eigen basis at least that large.

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
