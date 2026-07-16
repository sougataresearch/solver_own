# Testing Strategy — sougata_solver

Physics-software testing differs from typical application testing: the
dominant risk is not a crash but a **plausible-looking wrong number**
(wrong sign convention, wrong Fourier-factorization rule, wrong branch
cut). Every tier below exists to catch that class of bug specifically, per
`rules.md`'s core testing rule: **no new physics capability merges without
an oracle-comparison test.**

## Current State

`tests/` (pytest, `testpaths = ["tests"]` in `pyproject.toml`):
- `tests/test_analytic_fresnel.py` — Phase 1 R/T validated against
  `tests/oracles/fresnel.py` (independent closed-form Fresnel/TMM).
- `tests/test_fourier_factorization.py` — despite the name, currently only
  tests raw shape-level Fourier transforms (`Circle`/`Rectangle` DC value,
  `contains()` geometry) at the `geometry.py` level; the actual
  Fourier-factorization (Toeplitz) tests are a Phase 2 deliverable, not yet
  written — **this file will need substantial new content, not just
  extension, once Phase 2 lands.**
- `tests/test_polarimetry.py` — Jones/Mueller validated against known
  reference matrices (identity, ideal polarizer) plus a physical invariant
  (isotropic media don't couple s/p polarizations).
- `tests/conftest.py` — a seeded `rng` fixture (`np.random.default_rng(0)`)
  for any test needing reproducible randomness.
- `pyproject.toml` already defines a `slow` pytest marker for
  convergence/benchmark studies excluded from the default run — use this
  for every convergence-vs-`num_orders`/vs-`N` study in Phases 5 and 8.

## Testing Strategy By Tier

### Unit Testing

Scope: individual functions in isolation — a shape's Fourier transform at
a specific `(kx, ky)`, `_select_q_branch`'s behavior on a hand-picked
`q_sq` value, `star_product`'s algebraic identity for trivial (identity)
inputs, etc. Every new module (`fourier_factorization.py`,
`Lattice1D`/`Slab`, staircase generator) needs unit tests at this level
covering: nominal input, a degenerate/edge case (e.g. DC term, zero-size
shape, `N=1` staircase), and — where a closed form exists — an exact
comparison, not just "doesn't crash."

### Integration Testing

Scope: a full `Simulation.solve()` call exercising multiple modules
together (eigenmode solve → S-matrix cascade → field extraction). Every
script under `structures/` is implicitly an integration test (it must run
end to end and produce a sane number), and every `postprocessing/` script
is implicitly an integration test of the derived-quantity math against
already-computed raw data; formalize the most important ones as actual
pytest tests too (not just runnable scripts) — e.g.
`structures/thin_film/sio2_on_si_thin_film.py`'s SiO2-on-Si case should have a
corresponding assertion-based test, not rely on eyeballing printed output.

### System Testing

Scope: the whole pipeline against a **named external oracle** — this is
the tier that actually validates physics correctness, not just internal
consistency:
- Phase 1: analytic Fresnel/TMM (done, `tests/oracles/fresnel.py`).
- Phase 3: a published 1D binary-grating diffraction-efficiency benchmark
  (Moharam & Gaylord 1995 or equivalent — see `references.md`).
- Phase 4: S4 itself, driven as a subprocess oracle if buildable in this
  environment (check this explicitly before assuming it — see
  `memory.md`'s Known Issues), or a published 2D benchmark otherwise.
- Phase 6: a closed-form birefringent-material benchmark (uniaxial
  waveplate at normal incidence).
- **Never substitute a fabricated "it matches" claim if the real oracle
  can't be run** — say so explicitly and mark the test `xfail`/skipped
  with a clear reason, per `rules.md`'s AI Coding Rules.

### Performance Testing

Deferred until Phase 9 in substance (per `PRD.md`'s "correctness over
speed" non-functional requirement), but two things belong here earlier:
- Any `slow`-marked convergence study is itself a (manual) performance
  data point — record roughly how `num_orders`/`N` trades off against
  wall-clock time in the study's own output/docstring, so Phase 9's
  eventual profiling has a baseline to compare against.
- Phase 9's mandatory regression test: vectorized sweep must match the
  unvectorized per-point loop numerically (see `rules.md`'s Performance
  Requirements) — this is a correctness test *of* a performance change,
  not a benchmark in the traditional sense.

### Security Testing

Minimal surface (see `architecture.md`'s Security Considerations) — the
only concrete item today is malformed-CSV handling in
`material_from_csv` (`structures/thin_film/sio2_on_si_thin_film.py`). No dedicated security
test suite is warranted at current scope; revisit if a structure-definition
file format or any network/multi-user surface is ever added.

### Regression Testing

- Every bug fix gets a test reproducing the bug first (standard practice,
  not yet needed since no post-release bug has occurred, but the rule
  applies from Phase 2 onward).
- Every phase's "reduces to a simpler already-validated case" check is a
  standing regression guard, not a one-time task — e.g. Phase 4's
  patterned-layer solver should reduce to Phase 1's uniform result when
  the shape material equals the background (already listed as a Phase 4
  task in `tasks.md`); keep that test in the suite permanently, it will
  catch future refactoring bugs too.
- Run the full (non-`slow`) suite before every commit that touches
  `src/sougata_solver/` — see `rules.md`'s Code Review Checklist.

### Acceptance Testing

Defined per `PRD.md`'s Acceptance Criteria: a functional requirement is
"done" only when it has (a) a passing oracle-comparison test and (b) a
runnable example script. This file's job is to make sure every phase's
"done" claim in `phases.md`/`memory.md` is backed by an actual test in
`tests/`, not just an example that happens to print a plausible number.

## Running Tests

```bash
pytest                    # fast suite (excludes `slow`)
pytest -m slow            # convergence/benchmark studies only
pytest -m "not slow"      # explicit form of the default
pytest tests/test_analytic_fresnel.py -v   # one file, verbose
```

## What Is Explicitly Not Required (at current scope)

- Code-coverage percentage targets — a physics library's risk is
  concentrated in a small number of formula-bearing functions, not spread
  uniformly; oracle-comparison depth matters more than line coverage here.
- Mutation testing, property-based/fuzz testing frameworks — not
  proportionate to current scope; revisit only if a specific class of bug
  (e.g. numerical edge cases in `_select_q_branch`) shows a pattern that
  would benefit from property-based generation.
