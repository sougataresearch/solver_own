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

### Physical-Invariant Testing (oracle-independent)

Scope: checks that follow from Maxwell's equations / Poynting's theorem
directly, rather than from matching an external reference implementation.
These are strictly weaker than an oracle-comparison test (they can't catch
every wrong-but-self-consistent bug) but strictly cheaper and
oracle-independent — they don't need S4 buildable, a paper table
transcribed, or a second implementation trusted. Two are required starting
Phase 3 (the first patterned-layer phase), in addition to, not instead of,
the oracle-comparison test in System Testing below:

- **Energy conservation**: for lossless materials (real, positive `n`, no
  absorption), `R + T + sum(diffraction efficiencies) = 1` to within
  solver-precision tolerance; for absorbing materials,
  `R + T + sum(diffraction efficiencies) + A = 1` where `A` is computed
  from the imaginary part of the layer permittivities (Poynting-flux
  divergence), not fit to make the identity hold. A failure here means a
  sign/normalization bug regardless of whether the R/T numbers happen to
  look plausible or even happen to match a specific oracle test case by
  coincidence.
- **Convergence-rate-vs-theory**: not just "does R/T converge as
  `num_orders` increases" (that's Phase 5/8's convergence-vs-`N`/
  `num_orders` studies) but *at what rate*, compared to the rate Li (1996)
  predicts for whichever Fourier-factorization rule (`epsilon_hat` vs.
  `epsilon_inv_hat`) is in use at a discontinuous interface (see
  `references.md`'s Li 1996 entry and `design.md`'s Fourier-factorization
  section). A convergence curve that flattens at the right *value* but the
  wrong *rate* is a real finding — the classic case is using the direct
  rule where the inverse rule is needed at a discontinuity, which still
  converges, just first-order instead of the improved rate the inverse
  rule gives.

Both checks are cheap to add once `SimulationResult` already has R/T/
diffraction-order data (no new solver machinery), so add them to each new
geometry phase's existing test file rather than deferring to Phase 8 —
Phase 8's job is to run them systematically across every geometry type in
one place, not to be the first place they're checked.

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
- Phase 4a/4b: S4 itself, driven as a subprocess oracle if buildable in
  this environment (check this explicitly before assuming it — see
  `memory.md`'s Known Issues), or a published 2D benchmark otherwise. 4a
  uses a moderate-contrast case; 4b adds a deliberately near-degenerate/
  high-contrast stress case against the same oracle.
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
  standing regression guard, not a one-time task — e.g. Phase 4a's
  patterned-layer solver should reduce to Phase 1's uniform result when
  the shape material equals the background (already listed as a Phase 4a
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

## Test Taxonomy (Category 17 target 17.1)

Every test file already carries two forms of tier identification, audited
this session (not assumed) by reading a representative sample across the
full `tests/` directory before writing this section:

1. **Filename** names the *feature under test* (`test_reciprocity.py`,
   `test_geometry_validation.py`, ...) — deliberately not the tier, since
   a single file routinely spans several tiers at once (e.g.
   `test_1d_grating.py` has unit, regression, physical-invariant, *and*
   system-tier tests in one file, and says so explicitly in its own
   docstring: `"Tiers per testing.md: unit ..., regression ..., physical-
   invariant ..., and system ..."`).
2. **Docstring** cites the `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` category/
   target the file satisfies and names its oracle/invariant in prose —
   this file's own "Validation Inventory" (target 14.1) is the
   authoritative per-feature cross-reference from capability to test file
   to tier, so tier information isn't duplicated a third time per file.

**What this session adds**: a `pytest` marker, `oracle`
(`pyproject.toml`), applied to every test file that actually `import`s
from `tests/oracles/` (confirmed by grepping for that import, not
guessed) — `tests/test_1d_grating.py`, `test_2d_pillar.py`,
`test_2d_pillar_stress.py`, `test_analytic_fresnel.py`,
`test_anisotropic_inplane.py`, `test_anisotropic_uniform.py`,
`test_optical_outputs.py`, `test_thin_film_empy_cross_check.py` (8 files,
136 tests). This makes the System Testing tier **queryable**
(`pytest -m oracle`) the same way `slow` already makes the Performance-
adjacent convergence/benchmark studies queryable
(`pytest -m slow`/`pytest -m "not slow"`) — a marker was added only to
files meeting a precise, greppable criterion (imports a named external
oracle module) rather than a subjective per-file tier judgment call, so
the marker can't silently drift from what it claims to mean.

Other tiers (unit, physical-invariant, integration, regression,
acceptance) are **not** given their own marker in this pass: unlike
"oracle" (a precise import-based criterion) or "slow" (an explicit,
already-established opt-out), those tiers overlap heavily within a single
file and a single test function often serves more than one tier at once
(e.g. a reduction-to-simpler-case test is simultaneously a unit test and
a permanent regression guard) — forcing a single marker per test would
misrepresent that overlap rather than clarify it. The filename+docstring+
Validation-Inventory combination already documented above remains the
taxonomy for those tiers, per this target's own "filenames/docstrings or
pytest markers" wording (either is sufficient, not both required
everywhere).

## Running Tests

**Note**: `pyproject.toml` sets no default marker filter, so plain
`pytest` (no `-m` flag) runs *everything*, `slow`-marked tests included —
use `-m "not slow"` explicitly for the fast-development-loop subset.

```bash
pytest                    # everything, including slow-marked studies
pytest -m "not slow"      # fast-development-loop suite (excludes `slow`)
pytest -m slow            # convergence/benchmark studies only
pytest -m oracle          # system-tier tests cross-checked against a named external oracle (136 tests, ~100s)
pytest tests/test_analytic_fresnel.py -v   # one file, verbose
```

## Validation Inventory (Category 14 target 14.1)

Every public feature, its validating oracle/invariant test, its runnable
example, and any known limitation — kept in one table so a "done" claim
in `memory.md`/`phases.md` is always traceable to an actual test file, per
this file's own Acceptance Testing rule above. Updated whenever a new
public capability ships (the same discipline `rules.md` AI Coding Rule 6
already requires for `memory.md`/`decisions.md`).

| Feature | Oracle / invariant test | Example | Known limitation |
|---|---|---|---|
| Uniform multilayer R/T (Phase 1) | `tests/oracles/fresnel.py` (from-scratch Fresnel/TMM) + `tests/oracles/empy_tmm.py` (EMpy) | `structures/thin_film/sio2_on_si_thin_film.py` | none |
| Jones/Mueller polarimetry | Known reference matrices (identity, ideal polarizer) | `postprocessing/jones_mueller_ellipsometry.py` | s/p phase convention not externally matched (see Category 10 target 10.5 entry) |
| Toeplitz Fourier factorization (Phase 2) | Rasterize-and-sum + FFT-of-raster (RCWA.jl `convmat2D.py` reproduction) | n/a (library-internal) | none |
| 1D lamellar grating (Phase 3) | `tests/oracles/rcwa_1d_gaylord.py` (Moharam/Gaylord) | `structures/trench/trench_grating.py` | TM convergence is slow at sharp interfaces (Li's-rule-sensitive), documented not hidden |
| 2D patterned layer, moderate contrast (Phase 4a) | `tests/oracles/rcwa_2djl_eigenvalues.py` (independent RCWA.jl eigenoperator) | `structures/via/pillar_array.py` | **no external R/T oracle** — see Category 14 targets 14.2-14.4, evaluated and deferred, `references.md` |
| 2D patterned layer, near-degenerate/ill-conditioned (Phase 4b) | Same RCWA.jl oracle, high-contrast stress cases | `tests/test_2d_pillar_stress.py` | same external-oracle gap as above |
| Tapered sidewalls (Phase 5) | Zero-taper reduction + convergence-vs-`num_slices` study (no external oracle exists for this technique family) | `structures/via/tapered_via.py`, `structures/trench/tapered_trench.py` | convergence can be slow for steep tapers (`tests/test_staircase.py`) |
| Anisotropic uniform/patterned layers (Phase 6 / Category 1) | Closed-form birefringence benchmark + `tests/oracles/rcwa_anisotropic_inplane_jl.py` | `tests/test_anisotropic_*.py` | longitudinal coupling (target 1.5) explicitly deferred, no citable formulation found |
| Numerical-methods robustness (Category 2) | `design.md` Failure Contract, backed 1:1 by `tests/test_failure_contract.py` | n/a (library-internal) | none |
| Fourier-factorization rule inventory (Category 3) | `tests/test_fourier_factorization_rules.py`, measured convergence fixtures | n/a | FFF/NVM deferred, `decisions.md` ADR-012 |
| Geometry engine (Category 4) | From-scratch rasterized references, PNPoly point-in-polygon | `structures/via/elliptical_pillar.py`, `triangular_pillar.py` | GDS/raster import out of scope |
| Material dispersion models (Category 5) | BK7 published Sellmeier index, Rakić et al. (1998) published Au/Ag/Al/Ti coefficients | `structures/thin_film/tio2_sio2_dbr_on_si.py` | none |
| Boundary conditions/excitation (Category 6) | Symmetry invariants (rotational at normal incidence, azimuthal), Stokes reciprocity | `tests/test_polarization_states.py`, `tests/test_bottom_incidence.py` | none |
| Layer handling (Category 7) | Equivalence to uncached/unrepeated representations | `tests/test_layer_cache.py`, `tests/test_layer_absorption.py` | interior-amplitude reconstruction can overflow for thick/highly-lossy/high-`num_orders` cases (`troubleshooting.md`) |
| Solver sweeps (Category 8) | Equivalence to manual per-point `solve()` loops | `tests/test_sweep.py`, `tests/test_harmonic_convergence.py` | none |
| Real-space field reconstruction (Category 9 / Phase 7) | Analytic plane wave, transversality, interface continuity, flux-matches-R/T | `structures/trench/trench_field_cross_section.py`, `structures/via/pillar_field_cross_section.py` | same interior-amplitude overflow limitation as Category 7 |
| Optical outputs (Category 10) | `tests/oracles/fresnel.py::multilayer_complex_rt` (both polarizations), classical grating equation | `tests/test_optical_outputs.py` | per-order s/p conversion (target 10.5) deferred, no external validation achieved |
| Semiconductor OCD features (Category 11) | Closed-form rounded-rectangle area, periodicity self-consistency (overlay) | `structures/via/tsv_ocd_sweep.py`, `structures/trench/trench_ocd_sweep.py` | stochastic LER/LWR (target 11.8) explicitly deferred |
| Linear algebra (Category 12) | Bit-for-bit equivalence to pre-refactor results, measured density | `profiling/baseline_profile.py` | none |
| Performance optimization (Category 13) | Bit-for-bit-scale equivalence (eigenmode cache, vectorized sweep) | `profiling/benchmark_suite.py` | GPU/autodiff not approved (target 13.6) |
| Reciprocity (Category 14 targets 14.5/14.6) | Snell's-law-matched transmittance symmetry | `tests/test_reciprocity.py` | uniform layers only — does not extend to patterned/diffractive layers, verified not assumed |
| Harmonic convergence matrix (Category 14 target 14.7) | `sweep.find_convergence_index` across every geometry family | `tests/test_harmonic_convergence_matrix.py` | tapered/high-contrast cases converge slowly, documented per-fixture, not hidden |

**Standing gap, not silently dropped**: no external, third-party 2D R/T
oracle exists for this project (Category 14 targets 14.2-14.4) — S4 is not
buildable in this environment (no `cmake`/Lua/C++ toolchain, re-confirmed
2026-08-05) and no published numeric benchmark table was located via a
bounded literature search (same conclusion Phase 4a/4b already reached).
Every 2D-patterned-layer capability above is instead cross-validated
against an independent, structurally-different eigenoperator formula
(`RigorousCoupledWaveAnalysis.jl`, hand-transcribed since Julia is not
installed either) — a real, if narrower, form of independent validation,
not a self-consistency check, but explicitly not the same as an external
R/T match. Revisit if S4/Julia become available in a future environment.

## Validation Report (Category 14 target 14.8)

**As of 2026-08-05** (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 14, this
category's own session).

**Environment**: Python 3.12.10, NumPy 2.5.1, SciPy 1.18.0, `numpy>=1.24`/
`scipy>=1.10` required (`pyproject.toml`). No compiled extensions, no GPU.

**Results**: 637 tests pass project-wide (627 fast + 10 `slow`-marked
convergence/benchmark studies), full fast+slow suite re-run and confirmed
green as of this category's own completion. See the Validation Inventory
above for the per-feature oracle/test breakdown; see `memory.md`'s
"Current Project Status" for the full per-category narrative history.

**Tolerances actually used** (not a single project-wide constant — the
right tolerance depends on what's being compared, per `rules.md`'s
"tolerance-scale" lesson from Category 1 target 1.4's first mistake):

| Comparison class | Typical tolerance | Example |
|---|---|---|
| Bit-for-bit refactor equivalence (same formula, different code path) | `1e-12` to `1e-16` | `tests/test_linear_algebra_audit.py`, `tests/test_vectorized_sweep.py` |
| Oracle-comparison (independent formula/source) | `1e-6` to `1e-10` | `tests/test_analytic_fresnel.py`, `tests/test_1d_grating.py` |
| Energy conservation (`R+T+sum(DE)[+A]=1`) | `1e-6` to `1e-8` | `tests/test_layer_absorption.py`, `tests/test_reciprocity.py` |
| Convergence-vs-`num_orders` (genuinely slow-converging fixtures) | `1e-2` to `5e-2`, honestly matched to measured behavior, not tightened to force a pass | `tests/test_harmonic_convergence_matrix.py`'s tapered-via/high-contrast cases |

**Known, explicitly-scoped gaps** (not silently missing): Category 1
target 1.5 (longitudinal anisotropic coupling), Category 10 target 10.5
(per-order s/p conversion), Category 11 target 11.8 (stochastic LER/LWR),
Category 13 target 13.6 (GPU backend, not approved), and Category 14
targets 14.2-14.4 (external 2D R/T oracle — S4/Julia unavailable in this
environment, no published benchmark table located). Each has its own
`decisions.md` ADR or `references.md` entry recording why, not just that.

## What Is Explicitly Not Required (at current scope)

- Code-coverage percentage targets — a physics library's risk is
  concentrated in a small number of formula-bearing functions, not spread
  uniformly; oracle-comparison depth matters more than line coverage here.
- Mutation testing, property-based/fuzz testing frameworks — not
  proportionate to current scope; revisit only if a specific class of bug
  (e.g. numerical edge cases in `_select_q_branch`) shows a pattern that
  would benefit from property-based generation.
