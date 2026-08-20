# Validation Guide — What Each Oracle Proves (and Doesn't)

Target 18.8 of `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 18. This is
deliberately a **consolidation**, not new validation work — every claim
below already exists as a docstring in `tests/oracles/*.py` or a row in
`testing.md`'s Validation Inventory; nothing here asserts a new benchmark
comparison.

## How this differs from `testing.md`

`testing.md`'s **Validation Inventory** (Category 14 target 14.1) is
already a complete table, organized **by feature/category**: "this
capability is validated by that oracle, here's the known limitation."
Read it for that view — this document doesn't repeat it.

This guide is organized **by oracle** instead, and answers a narrower
question `testing.md`'s table doesn't spell out per-entry: *if a test
passes against this oracle, what has actually been established, and what
hasn't?* Two oracles both labeled "independent RCWA reference" can mean
very different strengths of evidence (a full R/T comparison vs. an
eigenvalue-only cross-check in a different gauge) — that distinction is
the whole point of this document.

## The oracles, one by one

### `tests/oracles/fresnel.py` — from-scratch analytic Fresnel/TMM

**Source**: written from scratch from the standard Born & Wolf/Macleod
thin-film (characteristic-matrix) formulation — independent of both
`sougata_solver` and any vendored repo.

**Proves**: uniform-multilayer R/T (both polarizations, any angle,
lossless or absorbing) to `1e-6`–`1e-10` (`testing.md`'s tolerance table).
Also backs `test_optical_outputs.py`'s `complex_amplitudes()` check via
`fresnel.py::multilayer_complex_rt`, at oblique incidence for **both**
s- and p-polarization.

**Does not prove**: anything about patterned (1D/2D) layers, anisotropy,
or Fourier factorization — a uniform stack has no lateral structure to get
wrong in those ways. It's also this project's *own* derivation (not a
third-party implementation), which is exactly why `empy_tmm.py` exists as
a second, differently-sourced check on the same physics (see below) —
agreement between two independently-derived implementations is stronger
evidence than either alone.

### `tests/oracles/empy_tmm.py` — transcribed from vendored EMpy

**Source**: hand-transcribed from `EMpy/EMpy/transfer_matrix.py`'s
`IsotropicTransferMatrix.solve` (a classic Abeles dynamical-matrix
formulation) — chosen over other vendored TMM sources after evaluating
them (one mixed plotting into the physics; another had undefined free
variables in its anisotropic path). `EMpy` itself is never imported at
runtime, per `decisions.md` — only read from and transcribed.

**Proves**: the same uniform-multilayer R/T as `fresnel.py`, from a
genuinely independent source and derivation route — this is what makes
Phase 1 validated against **two** oracles, not one restated.

**Does not prove**: anything `fresnel.py` doesn't already cover (same
scope: uniform layers only). Three intentional, documented deviations
from the raw EMpy source exist (found by running the transcription
against real absorbing materials) — see the file's own docstring if a
future discrepancy needs tracing back to a specific line.

### `tests/oracles/rcwa_1d_gaylord.py` — Moharam/Gaylord 1D grating

**Source**: hand-transcribed from the vendored
`Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples` project, itself citing
Moharam, Grann, Pommet & Gaylord (1995) — the classic published RCWA
benchmark for 1D lamellar gratings.

**Proves**: full R/T diffraction efficiencies for the Phase 3 1D-grating
solver (`structures/trench/trench_grating.py` is literally this
benchmark's geometry — see [`tutorials.md`](tutorials.md)'s grating tutorial).
This is a full-pipeline check: Fourier factorization, eigenmode solve, and
S-matrix cascading all have to be correct simultaneously for this to
match.

**Does not prove**: anything about 2D-patterned layers or anisotropy
(different code paths entirely). Known limitation carried in `testing.md`:
TM-polarization convergence is genuinely slow at sharp interfaces
(Li's-rule-sensitive) — documented, not hidden.

### `tests/oracles/rcwa_1d_pyrcwa.py` — second, independent 1D oracle

**Source**: hand-transcribed from the vendored `PyRCWA` project (MIT
license), which solves the **general** 2D P/Q eigenoperator and restricts
to 1D via harmonic truncation — a structurally different derivation route
from Gaylord's reduced TE-specific operator, even though both check the
same 1D physics.

**Proves**: a second, independently-routed confirmation of the same
1D-grating R/T Gaylord's oracle checks — two oracles that agree despite
using different math is stronger evidence than one oracle alone.

**Does not prove — explicitly scoped down**: **normal incidence, TE only**.
`PyRCWA`'s own oblique-incidence angle convention (`alpha`/`theta`) doesn't
obviously map onto this project's `(theta, phi)` convention, and that
mapping was deliberately left unresolved rather than asserted without
verification (same discipline as the sign-convention notes elsewhere in
this project). Do not read a passing test against this oracle as covering
oblique incidence or TM polarization — it doesn't.

### `tests/oracles/rcwa_2djl_eigenvalues.py` — 2D eigenoperator only

**Source**: hand-transcribed from `RigorousCoupledWaveAnalysis.jl`'s
isotropic patterned-layer eigenoperator (Julia is not installed in this
environment — transcribed by reading the source, never run).

**Proves**: `solve_layer_eigenmodes_patterned`'s `q^2` eigenvalues agree to
~1e-12 with a **completely different derivation route** — RCWA.jl
eliminates directly from Maxwell's curl equations into one matrix, a
different field-basis/gauge from this project's `Epsilon2 @ kp - coupling`
construction. This is exactly the check that caught a real bug during
development (an earlier draft copied the wrong S4 branch for `Epsilon2`).

**Does not prove — the most important limitation in this whole guide**:
**this checks eigenvalues only, never a full R/T pipeline.** `q^2`
eigenvalues are basis-independent so a different gauge is comparable at
all, but eigenvectors are never compared (they live in different bases by
construction) and R/T is never computed from this oracle. **No
independently-published 2D diffraction-efficiency benchmark exists for
this project** — surveyed and not found, not silently skipped: S4 isn't
buildable in this environment (no `cmake`/Lua/C++ toolchain), and RCWA.jl
itself has no fixed reference numbers to transcribe (its own test suite
uses randomized parameters and checks only internal self-consistency). See
`rcwa_2d_pillar.py`'s own docstring and `testing.md`'s "Standing gap, not
silently dropped" note. Every other piece of 2D-patterned-layer evidence
(`ky=0` reduction to the oracle-validated 1D solver, uniform-layer
reduction, energy conservation) is real and load-bearing, but it is not
the same class of evidence as `rcwa_1d_gaylord.py`'s full R/T match above
— know which one you're relying on.

### `tests/oracles/rcwa_anisotropic_inplane_jl.py` — anisotropic eigenoperator

**Source**: hand-transcribed from `RigorousCoupledWaveAnalysis.jl`'s
*uniform* anisotropic-layer branch (a different `Common.jl` function than
the isotropic-patterned oracle above, matching target 1.4's uniform-layer
scope).

**Proves**: `solve_layer_eigenmodes_uniform_inplane`'s `q^2` eigenvalues
agree to ~1e-13 with the same kind of structurally-different derivation
route as `rcwa_2djl_eigenvalues.py` — direct Maxwell-curl elimination vs.
this project's `Epsilon2 @ kp` construction.

**Does not prove**: full R/T for the anisotropic case from this oracle
alone (same eigenvalue-only caveat as above) — but unlike the isotropic
2D case, target 1.4's uniform in-plane tensor **also** reduces exactly to
target 1.3's closed-form birefringence benchmark and to the isotropic
`fresnel.py`/`empy_tmm.py` oracles when the tensor collapses to scalar
— so the uniform anisotropic case does have a genuine closed-form R/T
check, just not from this particular oracle file.

## Two things that are not oracles, and shouldn't be mistaken for one

- **`tests/regression_fixtures/`** — frozen, previously-trusted spectra
  compared against a fresh solve on every run. This catches *drift* in an
  already-validated code path; it does not itself establish physical
  correctness (that's what the oracles above already did, once). See
  `test_regression_fixtures.py`'s docstring for that fixture's provenance.
- **Physical-invariant tests** (energy conservation `R+T+sum(A)=1`,
  reciprocity, isotropic-reduction checks) — real, load-bearing evidence,
  and required *in addition to* an oracle per `rules.md`, but they are
  necessary-not-sufficient: a solver with a consistent sign error can
  still conserve energy. They catch a different class of bug than an
  oracle comparison does, which is why this project requires both, never
  either alone (`testing.md`'s Physical-Invariant Testing section).

## One honest, standing gap (repeated here deliberately, not softened)

No external, third-party full-R/T oracle exists for **any** 2D-patterned
capability in this project — isotropic or anisotropic. This is
`testing.md`'s Category 14 targets 14.2–14.4, evaluated and explicitly
not closed (S4/Julia unavailable in this environment; no published
benchmark table located). It's real, load-bearing evidence at the
eigenoperator level, and it's genuinely not the same strength of proof as
the 1D case's two independent full-R/T oracles. Revisit if S4 or Julia
ever become runnable in this environment.
