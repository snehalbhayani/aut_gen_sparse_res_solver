# symgen: exact, deterministic sparse resultant solver generator

`symgen` is a Python/SymPy reimplementation of the offline stage of this
repository. It needs no MATLAB or Maple licence to generate solvers. It reads
the same `problems/<name>.m` files and writes a MATLAB/Octave `solver.m` with
the same calling convention as the Maple-generated solvers:
`PEPsolutions = solver(data)`, one row per solution, with columns
`[a_h, other unknowns..., lambda]`.

```bash
pip install -r symgen/requirements.txt
python3 -m symgen problems/problem_relpose_6p_focal.m       # -> solvers_symgen/<name>/solver.m
python3 tests/test_symgen.py                                # unit + regression tests
```

Useful options: `--hidden -1` tries every hidden variable. `--sizes 2,3` searches
all polynomial combinations of those sizes instead of `cfg.polyComb`. `-v` logs
every candidate. `--seed` sets the GF(p) instance. `--test N` runs N random
accuracy trials.

## What changed compared to the Maple generator

The method is the same: a hidden variable, the extra equation `x_h - lambda`, a
monomial basis from Minkowski sums of Newton polytopes, and a Schur-complement
eigenvalue problem. How the template is *chosen* is different.

### 1. One exact validity condition replaces trial and error

The template has top rows `t*f_k` and lambda-rows `e*(x_h - lambda)` for `e`
in the eigen basis `E`. Every other column is in `G`. For each solution the
monomial vector satisfies `A1 v_E + A2 v_G = 0` and `lambda v_e = v_{e+e_h}`.
So the template is valid **if and only if `rank Top[:, G] = |G|`**, and each
unknown can be read off as a ratio of two monomials.

Every decision in symgen is an exact rank test of this condition. The tests
run over GF(p) (p = 2097143), with the data coefficients set to a seeded random
point. A full rank found mod p proves that a minor is a non-zero polynomial in
the data, so **every emitted template is valid for generic data**, not just for
the sampled instance. There are no random translations, no random restarts and
no `heurisiticTemplatesize`. The same input always gives the same solver.

### 2. Smallest eigenvalue problem without search (matroid step)

Dropping the lambda-row of `e` moves column `e` from `E` to `G`. The template
stays valid iff `Top[:, e]` is independent of the current `G` columns. This
is a matroid, so a single greedy pass reaches the **minimum eigen size for that
row set**: `|columns| - rank(Top)`. The greedy order only decides which
monomials stay in `E`; low degree ones are kept.

For `relpose_6p_focal`, every candidate basis from 30 to 517 monomials gives
exactly 18. That fits a stable count of 18 solutions for this formulation
(15 real ones plus 3 spurious ones), so no template of this kind should go
lower. An independent solution count has not been computed yet.

### 3. Exact greedy template reduction

There are two removal moves. Each is accepted only if the exact rank condition
still holds:

- **column move**: drop a `G` column together with every row that touches it;
- **row move**: drop a single row.

Columns that nothing touches any more disappear. The moves are tried in a fixed
monomial order: plain degree, and one variable-weighted order per unknown. The
smallest result is kept.

### 4. Leaner generated solver

- The lambda-rows are pure shifts, so there are no `B1`/`B2` blocks. Rows of
  `X` are either unit vectors or rows of `-A2\A1`.
- Coefficients are emitted once, with common subexpressions extracted. The
  matrices are filled with vectorised linear indexing instead of one statement
  per entry, which makes `solver.m` 5-10x shorter.
- Each unknown is recovered as `v(x_i*e) / v(e)` for the `e` with the largest
  `|v(e)|` in each solution. The Maple version took a product of eigenvector
  powers from `solForm`.

## Results (Octave 8.4, 200 random Gaussian instances, same data for both)

| problem | Maple template (elim + eigen) | symgen template | median log10 residual (old → new) | correct solutions per instance (old → new) | time per solve, Octave (old → new) |
|---|---|---|---|---|---|
| relpose_6p_focal | 12 + 18 = 30 | 12 + 18 = 30 | -12.02 → **-12.18** | 14.82 → **14.98** | 7.9 → **3.2 ms** |
| relpose_7p_fr_1s | 22 + 19 = 41 | 22 + 19 = 41 | -11.21 → **-12.44** | 18.32 → **18.96** | 5.9 → **1.0 ms** |
| 9pt2radial | 90 + 27 = 117 | 90 + 27 = 117 | -7.90 → **-8.09** | 18.13 → **19.03** | 15.7 → **1.4 ms** |
| opt_pnp_hesch | no solver in repo (cfg target 113) | 87 + 27 = 114 | → -10.67 | → 25.61 | → 0.7 ms |

A solution counts as correct when max|f_k| / ‖x‖ < 1e-6. Generation takes
1–35 s per problem.

**Template sizes:** symgen matches the hand-tuned Maple templates. It does not
beat them on these problems. What is new is that it reaches them with a fixed
sequence of steps, proves each one valid, and gives evidence that the eigen
size cannot go lower. On `9pt2radial`, keeping a larger eigen basis than the
minimum only removes one eliminated column per extra eigen column, so it never
pays off.

## Removing the spurious solutions of 6p focal

The 3 extra eigenvalues of `relpose_6p_focal` are solutions at
**a3 = infinity**. The `a3²` coefficient of the trace constraint is
`W2 = F33·(2·f3·r3 − F33·F)`, where `f3` and `r3` are the third column and
row of F. W2 = 0 together with det F = 0 means `F33 = 0` (a line) cut with a
cubic, which gives exactly 3 points. The leading coefficients also satisfy
syzygies, for i, j ∈ {1, 2}:

    F33·W2_i3 − F_i3·W2_33 = 0
    F33·W2_3j − F_3j·W2_33 = 0
    F33·W2_ij − 2·F_3j·W2_i3 + F_ij·W2_33 = 0

So the same combinations of the `E_ij` have `a3`-degree 1. They vanish at every
real solution but not at the points at infinity: a **symbolic saturation**
whose coefficients are explicit polynomials in the data, with no online
elimination step. `symgen/derive_6p_focal_sat.py` writes
`problems/problem_relpose_6p_focal_sat.m`, in which the 8 new equations are
added as eqs 11-18.

The new equations are products of existing ones, so fully expanding them
in the data doubles the coefficient code. A `Problem` can instead carry
intermediate coefficient symbols (`aux`). `make_problem()` keeps the `E_ij`
coefficients as intermediates, and the generated solver computes each new
coefficient as a short sum of products of numbers it already has.

```bash
python3 -m symgen symgen.derive_6p_focal_sat:make_problem --comb 10,11,19 --hidden 2
```

Octave, 500 instances, residuals measured on the original 10 equations:

| 6p focal solver | template (elim + eigen) | eigenvalues returned | median log10 residual | correct solutions | ms per solve |
|---|---|---|---|---|---|
| Maple (repo) | 12 + 18 | 18 (3 spurious) | -12.00 | 14.82 | 8.65 |
| symgen, original equations | 12 + 18 | 18 (3 spurious) | -12.12 | 14.99 | 3.57 |
| symgen, saturated, expanded | 25 + 15 | **15** | -11.00 | 14.80 | 7.15 |
| symgen, saturated, structured | 25 + 15 | **15** | -11.01 | 14.80 | 3.87 |

The eigenproblem now has exactly the 15 geometric solutions, and every hidden
variable works (the original equations had no valid template with `a3`
hidden). The cost: a 25×25 elimination instead of 12×12, and about one digit
less accuracy. In pure linear algebra the two are close: 77 µs for 25+15
against 70 µs for 12+18.

## Eigenvalue problem size vs. the matrix to invert

`eigen_target` (CLI `--eigen-target N`) stops the eigen basis reduction early.
Each l-row kept makes the eigenvalue problem one larger and A2 one smaller.
`--fastest` sweeps this whole family, measures the LAPACK time of each
template's linear algebra (`A2\A1` + `eig`) and keeps the fastest template
that still finds all the solutions.

Measured on every problem in the repo (Octave/LAPACK, µs for `A2\A1` + `eig`):

| problem | minimal eigen size | | larger eigen sizes | | |
|---|---|---|---|---|---|
| 9pt2radial | 27 + 90: **436** | 33 + 84: 520 | 39 + 78: 623 | 57 + 60: 1185 | 80 + 40: 2953 |
| opt_pnp_hesch | 27 + 87: **416** | 34 + 77: 660 | 41 + 70: 739 | 55 + 53: 1018 | 84 + 33: 3217 |
| 7p_fr_1s | 19 + 22: **90** | 21 + 20: 93 | 23 + 18: 115 | 26 + 15: 162 | |

**The minimal eigenvalue problem was the fastest in every case.** The
nonsymmetric eigensolver (Hessenberg QR) costs about 20–30× more per n³ than
the LU solve. With this kind of template, each extra eigen column removes only
one column of A2. Trading would pay only if one extra eigen column removed
about 20 eliminated columns. Larger eigen bases also return spurious
eigenpairs, which hurt accuracy. On `opt_pnp_hesch` they do recover a few more
of the 27 true solutions (about 26.9 against 25.4). The sweep is kept for
problems where the elimination dwarfs the eigenvalue problem.

## Limits / not yet done

- `cfg.noOfRowsToReduce > 0` (pre-elimination of the input equations) is
  not implemented. All problems in the repo use 0.
- `problems/problem_stitching.m` uses an older format that needs a missing
  `eqs/` folder, so it cannot be parsed.
- A rank deficiency seen mod p can be a false alarm with probability about
  deg/p. That can only make a template less reduced, never invalid.
- Laurent multipliers (negative exponents) are not used. All rows are
  polynomial multiples.

## Layout

| file | role |
|---|---|
| `ff.py` | exact GF(p) rank and incremental span |
| `polytope.py` | exact Newton polytopes, Minkowski sums, shifted lattice points |
| `problem.py` | `Problem` class and `problems/*.m` parser |
| `template.py` | candidate bases, validity test, eigen basis and template reduction, search |
| `numeric.py` | NumPy reference solver and accuracy statistics |
| `codegen.py` | MATLAB/Octave `solver.m` and `template.json` |
| `bench.py` | runs any `solver.m` in Octave and scores it |
