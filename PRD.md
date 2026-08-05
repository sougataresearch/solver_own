# Product Requirements Document — sougata_solver

## Problem Statement

Simulating light interaction with periodic nanostructures (thin films,
multilayer stacks, gratings/trenches, via/pillar arrays) requires solving
Maxwell's equations under periodic boundary conditions. Commercial tools
(JCMsuite, used in the vendored `EMTutorial/` reference projects) and
existing open packages (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl`, all
vendored under `Solver_own/` as reference material) solve this, but:

- Commercial tools are closed-source and not scriptable/extensible for
  research needs the way an owned codebase is.
- The vendored reference solvers are either not Python (S4 is C++/Lua,
  RCWA.jl is Julia) or are not structured the way this project's user wants
  to reason about and extend the physics (EMpy's RCWA module and mode
  solvers are a different design).
- There is no existing tool the user directly controls, end to end, that
  can be trusted for via/trench/pillar problems with tapered sidewalls
  while being auditable line-by-line against a known-correct reference.

`sougata_solver` exists to close this gap: an owned, from-scratch, fully-understood
RCWA implementation, validated against the vendored references rather than
depending on them at runtime.

## Goals

1. Correctly compute reflectance/transmittance/diffraction efficiencies for:
   thin films, arbitrary multilayer stacks, 1D-periodic lamellar gratings
   (trenches), 2D-periodic patterned layers (vias, pillars), and tapered
   (sloped-sidewall) versions of the latter two.
2. Every physics formula is traceable to a cited source (S4 source line
   numbers, a named paper, or an independently-derived analytic check) —
   never a from-memory guess.
3. Every new geometry/physics capability ships with a validation test
   against an independent oracle before being considered "done."
4. Keep the codebase small, readable, and dependency-light (NumPy/SciPy
   only) through at least Phase 8.

## Success Criteria

- Phase 1 (uniform multilayer stacks): reflectance/transmittance match
  analytic Fresnel/TMM to numerical precision across incidence angle,
  polarization, and a dispersive-material wavelength sweep. **Met** — see
  `tests/test_analytic_fresnel.py` and `tests/oracles/fresnel.py`.
- Phase 2 (Fourier-factorization core): direct and inverse-rule Toeplitz
  permittivity matrices match two independent numerical references.
  **Met** — see `tests/test_fourier_factorization.py`.
- Phase 3 (trench): diffraction efficiencies match a published 1D
  binary-grating benchmark (e.g. Moharam & Gaylord 1995) to within
  numerical-truncation-limited agreement, for at least TE and TM
  polarization at oblique incidence; energy conservation holds and the
  measured convergence rate vs. `num_orders` matches theory (`testing.md`'s
  Physical-Invariant Testing) — required starting this phase, not deferred
  to Phase 8. **Met** — see `tests/test_1d_grating.py`.
- Phase 4a (via/pillar, well-conditioned case): reflectance/transmittance
  and/or diffraction efficiencies match an S4-driven reference simulation
  of an equivalent structure (same lattice, radius, materials, wavelength)
  to within numerical-truncation-limited agreement, and satisfy the
  energy-conservation invariant (`testing.md`). **Met**, with a caveat
  honestly recorded rather than fabricated: S4 itself was never runnable in
  this environment, so the eigenoperator is instead cross-checked against
  an independent `RigorousCoupledWaveAnalysis.jl`-derived oracle (see
  `tests/test_2d_pillar.py`, `phases.md` Phase 4a).
- Phase 4b (via/pillar, near-degenerate/ill-conditioned case): the same
  agreement holds for at least one high-contrast, small-feature-to-period,
  high-`num_orders` stress case, with any ill-conditioning explicitly
  logged rather than silently degrading the result. **Met** — see
  `tests/test_2d_pillar_stress.py`.
- Phase 5 (tapered sidewalls): R/T demonstrably converges (monotonically,
  within expected discretization error) as the number of staircase slices
  increases, for both a tapered via and a tapered trench. **Met** — see
  `tests/test_staircase.py`.
- Phase 6 (anisotropic materials), Category 1 targets 1.3-1.4 and 1.6-1.8
  (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): uniform diagonal-tensor, uniform
  in-plane-coupled, and patterned anisotropic layers reduce correctly to
  their isotropic/simpler special cases, agree with a closed-form uniaxial
  benchmark and an independent `RigorousCoupledWaveAnalysis.jl`-derived
  oracle, and satisfy energy conservation. **Met for that scope.** Target
  1.5 (longitudinal `eps_xz/eps_yz/eps_zx/eps_zy` coupling) is evaluated
  and explicitly deferred — no citable, independently-benchmarkable
  formulation was located, per `references.md`'s "Target 1.5 bounded
  literature search." See `tests/test_anisotropic_uniform.py`,
  `tests/test_anisotropic_inplane.py`, `tests/test_anisotropic_patterned.py`,
  `tests/test_anisotropic_degeneracy.py`, `tests/test_mode_classification.py`.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Categories 2-5 (numerical methods,
  Fourier factorization, geometry engine, material models): every target
  (2.1-2.5, 3.1-3.6, 4.1-4.7, 5.1-5.8) is **met** — see that file's own
  per-target status entries for the exact validation evidence (each cites
  its test file(s) directly rather than being restated here). Highlights:
  a documented failure contract and eigenvalue-diagnostics report; a
  Fourier-factorization rule inventory with measured convergence fixtures;
  `Ellipse`/`Polygon` shapes and a JSON pattern-import format; five
  analytic dispersion models validated against BK7's and Rakić et al.
  (1998)'s independently-published values.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 6 (boundary conditions and
  excitation): targets 6.1-6.6 all **met** — polarization-state and
  azimuthal-rotation symmetry invariants, a characterized grazing-incidence
  boundary, an oblique-incidence Rayleigh-threshold case, and a Stokes-
  reciprocity-verified finding that bottom illumination needs no new API
  (`decisions.md` ADR-014).
- Phase 7 (real-space field reconstruction), tracked at atomic-target grain
  as `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 9: targets 9.1-9.8 all
  **met**. Reconstructed E/H fields match the analytic plane wave (uniform
  layer), exact tangential-field continuity across a genuine material
  interface, 1D periodicity, and real-space Poynting flux matching solver
  `R`/`T` to `~1e-6` for a 2D pillar (this category's own exit criterion).
  See `tests/test_field_reconstruction.py` and Phase 7's status entry in
  `phases.md`.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7 (layer handling): targets
  7.1-7.6 all **met** — construction-time layer-thickness validation, a
  repeated-layer-identity regression guard, an instance-scoped Toeplitz-
  matrix cache gated on a measured timing case (`decisions.md` ADR-016),
  and `SimulationResult.layer_absorption()` (a flux-divergence combination
  of already-validated field-reconstruction pieces, `decisions.md`
  ADR-017) satisfying the `R+T+sum(A)=1` energy-balance identity for a
  lossy fixture, with a found-and-documented numerical-overflow limitation
  for extreme thick/highly-lossy/high-`num_orders` cases
  (`troubleshooting.md`). See `tests/test_layer_validation.py`,
  `tests/test_layer_repetition.py`, `tests/test_layer_cache.py`,
  `tests/test_layer_absorption.py`.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 8 (solver sweeps and
  convergence): targets 8.1-8.8 all **met** — a typed `SweepResult`
  container and library-level wavelength/angle/polarization/thickness
  sweep functions, each confirmed equivalent to a manual per-point
  `Simulation.solve()` loop; a harmonic-order convergence study; a
  conservative convergence criterion (`decisions.md` ADR-018) validated
  against thin-film/trench/pillar fixtures per the category's own gating
  requirement before automatic harmonic-order selection was built on top
  of it. See `tests/test_sweep.py`, `tests/test_harmonic_convergence.py`.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 10 (optical outputs):
  targets 10.1-10.4/10.6 all **met** — complex per-order Cartesian field
  coefficients validated against both s- and p-polarization Fresnel-
  oracle comparisons, diffraction angles with a clear `None` non-
  propagating representation, a one-call conservation report, and a
  frozen output schema. Target 10.5 (per-order s/p conversion) is
  evaluated and explicitly deferred — a bounded attempt to externally
  validate the polarization convention against S4's actual source found
  a plausible but numerically-unconfirmed match. See
  `tests/test_optical_outputs.py`, `references.md`'s "Target 10.5 bounded
  external-validation attempt."
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 11 (semiconductor OCD
  features): targets 11.1-11.7 all **met** — a validated CD-first OCD
  parameter object, a trapezoid trench constructor, an arc-sampled-
  `Polygon` corner-rounding geometry converging to the closed-form
  rounded-rectangle area, reproducible TSV/grating OCD example sweeps,
  and overlay confirmed already achievable with no new API
  (`decisions.md` ADR-019). Target 11.8 (stochastic LER/LWR) is evaluated
  and explicitly deferred — fundamentally in tension with RCWA's
  periodic-Fourier formulation (`decisions.md` ADR-020). See
  `tests/test_ocd.py`, `tests/test_overlay.py`.
- `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 12 (linear algebra):
  targets 12.1-12.5 all **met** — a measured baseline performance profile
  (`profiling/baseline_profile.py`), a direct-inverse audit that found
  and fixed a house-convention inconsistency (confirmed bit-for-bit
  equivalent), a factorization-reuse design note, an opt-in SVD
  diagnostic, and sparse/iterative methods evaluated and **rejected** on
  a measured 100%-dense-matrix structural finding (`decisions.md`
  ADR-021). See `tests/test_linear_algebra_audit.py`,
  `tests/test_svd_diagnostics.py`.
- No phase is marked "done" without: (a) a passing automated test against
  an oracle, and (b) a runnable example script producing physically
  plausible output.

## Functional Requirements

| ID | Requirement |
|----|-------------|
| FR-1 | Solve reflectance/transmittance for an arbitrary stack of uniform, dispersive, isotropic layers at arbitrary incidence angle/azimuth/polarization. *(done)* |
| FR-2 | Support semi-infinite incidence/transmission half-spaces of arbitrary (possibly complex/absorbing) index. *(done)* |
| FR-3 | Report Jones and Mueller-matrix polarimetric response. *(done)* |
| FR-4 | Represent 2D-periodic in-plane patterns from `Circle`, `Rectangle`, `Ellipse`, and `Polygon` (simple, analytic, no raster/GDS import) primitives, including nested/overlapping shapes with correct area subtraction. *(done — `Circle`/`Rectangle` consumed by the solver since Phase 4a; `Ellipse`/`Polygon` added `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 4 targets 4.3/4.5)* |
| FR-5 | Solve reflectance/transmittance/diffraction efficiencies for a layer patterned according to FR-4 (via/pillar). *(done — Phase 4a, hardened for near-degenerate cases in Phase 4b)* |
| FR-6 | Represent and solve 1D-periodic lamellar (line/space) patterns (trench). *(done — Phase 3)* |
| FR-7 | Represent a feature (via/trench) with linearly tapered sidewalls via staircase layer discretization, and demonstrate R/T convergence with slice count. *(done — Phase 5)* |
| FR-8 | Support anisotropic (full 3×3 tensor) materials in both uniform and patterned layers. *(partially done — Phase 6 Category 1 targets 1.3/1.4/1.6-1.8: diagonal and in-plane-coupled tensors, uniform and patterned. Longitudinal coupling, target 1.5, explicitly deferred — see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`)* |
| FR-9 | Reconstruct real-space E/H field maps at an arbitrary depth in the stack (for cross-section visualization of trench/via structures). *(done — Phase 7 / `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 9, `fields.py`, `structures/trench/trench_field_cross_section.py`, `structures/via/pillar_field_cross_section.py`)* |
| FR-10 | Ingest dispersive material data from refractiveindex.info-style CSV `n,k` exports, or build a dispersive material from an analytic model (Sellmeier, Cauchy, Lorentz, Drude, Drude-Lorentz). *(done — `structures/thin_film/sio2_on_si_thin_film.py::material_from_csv`; analytic models added `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 5 targets 5.2-5.6)* |
| FR-11 | Import a `Pattern` from a minimal, safe, non-CAD JSON format (units, isotropic-scalar materials, `Circle`/`Rectangle`/`Ellipse`/`Polygon`/`Slab`). *(done — `geometry_io.py`, Category 4 target 4.6; parser only, not yet wired into `Simulation`/`Layer` construction)* |

## Non-Functional Requirements

- **Correctness over speed.** No approximation or optimization is adopted
  if it cannot be validated against an oracle; performance work is
  explicitly deferred to Phase 9.
- **Auditability.** Every formula-bearing function must cite its source in
  its docstring (existing convention — see `eigenmodes.py`, `smatrix.py`,
  `fields.py`). This is a hard requirement, not a suggestion — see `rules.md`.
- **No hidden state / no framework magic.** Plain dataclasses and functions;
  no metaclasses, no dependency injection containers, no plugin system.
- **Single-machine, CPU, pure Python/NumPy/SciPy** through Phase 8. No
  compiled extensions, no GPU requirement.
- **Reproducibility.** Given the same inputs (materials, geometry,
  wavelength, angle, `num_orders`), output must be bit-for-bit deterministic.

## User Stories

- As the project owner, I want to define a multilayer thin-film stack and
  get R/T vs. wavelength, so I can compare against measured ellipsometry
  data. *(done)*
- As the project owner, I want to define a 1D trench grating (period, line
  width, depth, materials) and get diffraction efficiencies, so I can
  reason about a lithography/etch scatterometry target.
- As the project owner, I want to define a 2D via or pillar array (period,
  radius, depth, materials) and get R/T and diffraction efficiencies, so I
  can reason about a TSV (through-silicon-via) scatterometry target — see
  the vendored `EMTutorial/Scatterometry/ThroughSiliconVia` JCMsuite
  reference case for the kind of structure this targets.
- As the project owner, I want to specify a sidewall taper angle for a via
  or trench and see the staircase-discretized result converge as I increase
  slice count, so I can trust the tapered-sidewall approximation.
- As the project owner, I want a cross-section field-intensity plot for a
  patterned structure, so I can visually sanity-check mode confinement /
  resonances the way the vendored `EMTutorial` galleries do.

## Acceptance Criteria

For each functional requirement above, "done" means: a merged implementation
+ a passing pytest test comparing against a named oracle (analytic formula,
published table, or S4 cross-check) + a runnable example script. A phase in
`phases.md` is not considered complete until every FR it claims meets this bar.

## Constraints

- Pure Python + NumPy + SciPy only (no compiled extensions) through at
  least Phase 8 — see Non-Functional Requirements.
- Single developer, part-time / research-pace effort — no team-coordination
  process overhead (no mandatory PR review by others, no release-train
  scheduling).
- `S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl` are reference/oracle
  material only, vendored as sibling directories under `Solver_own/` — they
  are not runtime dependencies of `sougata_solver` and must not be imported by it.
- Windows development environment (PowerShell primary shell) — anything
  written into deployment/CI docs must work there, not assume a Unix-only
  toolchain.

## Risks

| Risk | Mitigation |
|------|------------|
| Subtle sign/convention bugs in eigenmode or S-matrix math (the dominant historical bug class in RCWA implementations) | Mandatory source citation + independent-oracle test for every new formula (see `rules.md`) |
| Fourier-factorization rule chosen incorrectly for patterned layers (wrong convergence rate / wrong answer at discontinuous interfaces — the classic "Li's rules" pitfall) | Explicit `epsilon_hat` vs. `epsilon_inv_hat` Toeplitz construction per Phase 2 of `phases.md`, validated against FFT-of-rasterized-mask numerically, not assumed correct by inspection |
| General (non-uniform) complex eigenproblem can have degenerate/near-degenerate eigenvalues causing numerical instability | Reuse the already-validated `_select_q_branch` outgoing-mode convention; Phase 4a scopes its own test cases away from this regime deliberately, and Phase 4b is a dedicated phase for targeted regression tests + stability handling on near-degenerate cases, rather than an implicit sub-task of Phase 4a |
| Staircase approximation for tapered sidewalls converges slowly for steep angles | Explicit convergence-vs-slice-count test/example required before Phase 5 is considered done (see PRD Success Criteria) |
| Solo-developer bus factor / knowledge loss between sessions | `memory.md` and `decisions.md` are mandatory living documents, updated at the end of every substantive session |

## Out-of-Scope Items

- Arbitrary polygon / GDS-imported / rasterized-mask geometry (explicitly
  deferred per user decision; only parametric `Circle`/`Rectangle`/1D
  `Slab` shapes are in scope — see `decisions.md`).
- GPU backend, JAX/PyTorch autodiff, batched inverse design (Meent/TORCWA-
  style) — explicitly deferred to optional Phase 9, only after Phases 2-8
  are validated.
- Non-periodic / open (aperiodic) scattering problems (isolated particles,
  FDTD-style transient simulation) — RCWA is fundamentally a periodic-BC
  method, so `sougata_solver` itself will not grow FDTD capability. This is
  **not** a statement that FDTD is abandoned project-wide: per the project
  owner (2026-07-21), a time-domain (FDTD) solver is a genuine future goal
  of the broader EM-wave-solver effort `sougata_solver` is part of — RCWA
  was deliberately chosen first for simplicity (frequency-domain, periodic
  structures are the lower-risk starting point). FDTD is expected to be a
  separate effort (its own codebase/phases, not bolted onto this RCWA
  solver's module structure) once RCWA reaches a stable state; no FDTD
  phase, module, or timeline is defined yet — see `decisions.md` for the
  recorded scope clarification. `REFERENCE/meep`, `REFERENCE/gprMax`,
  `REFERENCE/fd3d`, `REFERENCE/maxwellfdfd` (FDTD/FDFD) and
  `REFERENCE/mfem`, `REFERENCE/OpenParEM`, `REFERENCE/dolfinx`+`ufl`+`basix`+`ffcx`,
  `REFERENCE/FreeFem-sources` (FEM) are already vendored for exactly that
  future effort — not currently used by anything in `sougata_solver`, and
  not evaluated as part of any RCWA phase's reference selection (per the
  `phase-reference-picker` skill's own guidance that a different numerical
  method is not a formula source for an RCWA phase).
- Nonlinear optics, magneto-optic materials, thermal/mechanical coupling.
- A GUI or web interface — this is a Python library driven by scripts.
- Public package distribution (PyPI) / multi-user support — see
  `deployment.md` for current solo-research-tool scope.
