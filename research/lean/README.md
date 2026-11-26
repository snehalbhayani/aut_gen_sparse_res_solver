# Lean formalisation (not yet compiled)

`FaceObstruction.lean` formalises the key lemma and Theorems 1–2 of
`../proof/face_obstruction.md`.

Lean could not be installed in the session where this was written, so the file
has **not been checked by Lean yet**. To build it:

```bash
cd research/lean
lake exe cache get
lake build
```

The first command downloads the prebuilt Mathlib cache. If the pinned version
(`lean-toolchain`, `lakefile.toml`, v4.19.0) is unavailable, bump both files to
a current Mathlib release tag. Expect a few small API-name fixes.
