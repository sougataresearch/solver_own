# Development Rules — pyrcwa

These rules formalize conventions the codebase *already follows* (verified
against the actual source, not invented) plus a few forward-looking ones
needed now that git/CI exist. Where a rule already has a concrete example
in the codebase, it's cited — these are not aspirational.

## Coding Standards

- **Python ≥ 3.10**, `from __future__ import annotations` at the top of
  every module (already universal in `src/pyrcwa/`) so modern union-type
  hints (`X | None`) work regardless of runtime version quirks.
- **Type hints on every function signature.** Already the case throughout
  `src/pyrcwa/`; keep it that way, including for new Phase 2-9 code.
- **`dataclasses` for data containers** (`Layer`, `LayerEigenmodes`,
  `SimulationResult`, `Circle`, `Rectangle`, `Pattern`,
  `PlaneWaveExcitation`) — do not introduce a class hierarchy or
  builder pattern where a dataclass suffices.
- **No broad `except`.** Let exceptions propagate; see `design.md`'s Error
  Handling section.
- **No hidden global state, no singletons, no module-level mutable state.**
  Every module in `src/pyrcwa/` today is pure functions + dataclasses; keep
  it that way.
- **Formatting/linting**: `.flake8` and `mypy.ini` exist in the parent
  `EMpy` reference project but **not yet in `pyrcwa`** — add both before
  Phase 2 lands (a `pyproject.toml` `[tool.ruff]` section, or a standalone
  `.flake8`, plus a `mypy.ini`), so new patterned-layer code is type-checked
  from day one rather than retrofitted later.

## Naming Conventions

- Module names: lowercase, single word or `snake_case`, matching the
  physics/data concept they own (`materials.py`, `eigenmodes.py`,
  `smatrix.py`, `fourier_basis.py`) — not generic names like `utils.py` or
  `helpers.py`.
- Function names describing **what physical quantity or operation** they
  produce, not implementation detail: `solve_layer_eigenmodes_uniform`,
  `interface_smatrix`, `propagation_smatrix`, `z_poynting_flux` — continue
  this pattern for Phase 3/4's `solve_layer_eigenmodes_1d`/
  `solve_layer_eigenmodes_patterned`.
- Private/internal helpers prefixed with `_` (`_solve`, `_rotate`,
  `_select_q_branch`, `_is_trivial_interface`) — already consistent, keep it.
- Physics variable names may be short and formula-like (`kx`, `ky`, `q`,
  `phi`, `kp`, `eps`) **only** inside functions whose docstring cites the
  source formula — never introduce a short cryptic name without a citing
  docstring nearby (see Documentation Standards below).

## Documentation Standards

This is the single most load-bearing rule in the project, and it already
has a consistent house style — **every module that implements a
non-trivial physics formula opens with a docstring citing its exact
source**, e.g.:

> "Formulas verified directly against `S4/S4/rcwa.cpp` (not paraphrased):
> `MakeKPMatrix` (lines 242-293)..." — `eigenmodes.py:3-6`

> "z_poynting_flux is transcribed directly from
> `S4/S4/rcwa.cpp::GetZPoyntingFlux` (lines 1846-1897), not re-derived,
> since a from-scratch re-derivation of the sign/normalization conventions
> ... risked introducing exactly the kind of subtle error this module
> needs to avoid." — `fields.py:3-8`

**Rule**: any new function implementing a physics formula (eigenmode
solve, Fourier factorization, field reconstruction, staircase geometry,
etc.) must include, in its docstring:
1. The exact external source (file + line numbers for S4/EMpy/RCWA.jl, or
   author + year + equation number for a paper), **or**
2. An explicit statement that it's independently derived, plus how it was
   validated (e.g. `excitation.py`'s honest note that its s/p sign
   convention "is internal and self-consistent (not yet matched to
   S4/EMpy's convention)").

Do not write a physics docstring from memory and mark it as sourced. If
you can't find/verify the exact source, say so explicitly and flag it for
extra test scrutiny — this project treats "looks right" as insufficient.

Additionally:
- Module-level docstrings state *why* a design choice was made (see
  `smatrix.py`'s docstring on the fixed-block-size convention), not just
  what the module contains.
- No comments restating what code already says (per general project style)
  — comments/docstrings exist to record the *why* and the *source*, not the
  *what*.

## Testing Requirements

See `testing.md` for the full strategy. Minimum bar, restated here as a
hard rule: **no new physics capability merges without an oracle-comparison
test** (analytic formula, published benchmark, or S4/EMpy/RCWA.jl
cross-check) — a test that only checks "it runs without crashing" is not
sufficient for anything under `src/pyrcwa/`.

## Git Workflow

- `pyrcwa` is now its own git repository (root at
  `Solver_own/pyrcwa/.git`), separate from the vendored reference repos
  (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl` each keep their own
  `.git` and are not submodules of `pyrcwa`).
- **Trunk-based, single branch (`main`) is fine** at solo-research-project
  scale — no mandatory feature-branch-per-task process. Use a short-lived
  branch only for a change you want to validate/convergence-test before
  merging (e.g. "does the Phase 4 general eigensolver actually match S4"),
  not as ceremony.
- **Commit at the granularity of one verified capability**, matching how
  the phases in `phases.md` are scoped — e.g. "add Toeplitz permittivity
  construction + FFT-cross-check test" is one commit, not bundled with
  unrelated changes.
- Never commit generated/output artifacts: `.gitignore` already excludes
  `__pycache__/`, `.pytest_cache/`, `*.egg-info/`, `output_RT.csv`, `*.png`.
- Never force-push, never rewrite published history, once this repo has
  any remote — solo project today, but form the habit now.

## Commit Message Format

```
<summary line, imperative mood, ≤ 72 chars>

<optional body: why this change, not what — the diff shows what>
```

Examples matching this project's actual content:
```
Add Toeplitz epsilon_hat/epsilon_inv_hat construction from Pattern shapes

Needed before the general (patterned-layer) eigenmode solver in Phase 4
can consume anything. Validated against FFT-of-rasterized-mask for a
Circle and a Rectangle pattern (tests/test_fourier_factorization.py).
```
No fixed prefix taxonomy (no mandatory `feat:`/`fix:` conventional-commits
scheme) is imposed — solo project, low volume, `git log --oneline` staying
readable is the actual bar, not a scheme.

## Code Review Checklist

Solo project — there is no second reviewer today — so treat this checklist
as **your own pre-commit self-review**, not a delegated step:

- [ ] Does every new physics formula cite its source (Documentation
      Standards above)?
- [ ] Does every new physics capability have an oracle-comparison test,
      not just a smoke test?
- [ ] Does the change avoid introducing a broad `except`, hidden global
      state, or a framework abstraction not justified by current need?
- [ ] Are type hints present and accurate (spot-check with `mypy` once
      configured)?
- [ ] Does `pytest` pass, including any newly-added tests?
- [ ] If a `NotImplementedError` was removed, does the corresponding
      `PRD.md` functional requirement get marked done, and does
      `memory.md`/`phases.md` get updated to reflect it?

## Security Rules

See `architecture.md`'s Security Considerations. Concretely:
- Never use `eval`/`exec`/`pickle` on any input, including future
  structure-definition files.
- Any future file-ingestion code (material CSVs, structure definitions)
  must raise a clear error on malformed input, not silently substitute a
  default.
- No dependency is added to `pyproject.toml` without checking it's
  actively maintained and has no known critical CVEs — small dependency
  surface is itself a security property here.

## Performance Requirements

Per `PRD.md`'s non-functional requirements: **correctness first, no
premature optimization**. Concretely:
- Do not introduce caching, memoization, or algorithmic shortcuts for
  Fourier factorization / eigenmode solves before Phase 9, unless a
  specific correctness-validated capability is measurably too slow to use
  (e.g. a convergence study in Phase 8 that can't complete in reasonable
  time) — and even then, validate the optimized path against the
  unoptimized one before trusting it.
- `scipy.linalg.lu_factor`/`lu_solve` (not `numpy.linalg.inv`) is the house
  convention for solving linear systems in `smatrix.py::_solve` — reuse
  this helper rather than reintroducing direct matrix inversion elsewhere.
- Vectorization work (Phase 9) must not change any numerical result versus
  the unvectorized path beyond floating-point-order-of-operations
  tolerance; add a regression test comparing both paths on at least one
  existing example before considering a vectorized path "done."

## AI Coding Rules (must never be violated)

These bind any AI assistant (this one included) working in this
repository, and codify the pattern already visible throughout the existing
code (e.g. `simulation.py`'s explicit "Phase 1 scope only" docstring,
`fields.py`'s explicit refusal to re-derive from scratch):

1. **Never invent a physics formula from memory and present it as
   verified.** If a formula can't be checked against S4/EMpy/RCWA.jl
   source or a citable published result, say so explicitly in the
   docstring and flag it as unverified — do not silently write a
   plausible-looking formula.
2. **Never silently implement a feature outside the current phase's
   scope.** If asked to touch patterned layers before Phase 2/3/4 land,
   either implement the prerequisite phase first or raise
   `NotImplementedError` naming the phase, exactly as `simulation.py`
   already does — never return a plausible-but-wrong number.
3. **Never remove or weaken an existing oracle-comparison test** to make a
   change pass. If a test starts failing, the formula or the test is
   wrong — find out which before touching either.
4. **Never add a dependency, a GPU/backend requirement, or a framework**
   without it being explicitly requested and recorded in `decisions.md` —
   Phase 9's GPU/autodiff backend is opt-in and later, not something to
   pull forward unasked.
5. **Never fabricate benchmark/validation numbers.** If a cross-check
   against S4 can't actually be run (S4 not built, no Lua bindings
   available, etc.), say so plainly rather than presenting made-up
   "it matches" results.
6. **Always update `memory.md`/`decisions.md`** when a phase completes or
   a non-obvious architectural choice is made, so a future session (AI or
   human) doesn't have to re-derive context that's already been settled.
7. **Never touch the vendored reference repos (`S4`, `EMpy`,
   `RigorousCoupledWaveAnalysis.jl`, `EMTutorial`) as part of `pyrcwa`
   development** — they are read-only oracles, not code to merge from or
   modify.
