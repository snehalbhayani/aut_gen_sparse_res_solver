# Lean formalisation (machine-checked)

`FaceObstruction.lean` formalises the key lemma and Theorems 1–2 of
`../proof/face_obstruction.md`.

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
