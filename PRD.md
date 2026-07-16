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
- Phase 3 (trench): diffraction efficiencies match a published 1D
  binary-grating benchmark (e.g. Moharam & Gaylord 1995) to within
  numerical-truncation-limited agreement, for at least TE and TM
  polarization at oblique incidence.
- Phase 4 (via/pillar): reflectance/transmittance and/or diffraction
  efficiencies match an S4-driven reference simulation of an equivalent
  structure (same lattice, radius, materials, wavelength) to within
  numerical-truncation-limited agreement.
- Phase 5 (tapered sidewalls): R/T demonstrably converges (monotonically,
  within expected discretization error) as the number of staircase slices
  increases, for both a tapered via and a tapered trench.
- No phase is marked "done" without: (a) a passing automated test against
  an oracle, and (b) a runnable example script producing physically
  plausible output.

## Functional Requirements

| ID | Requirement |
|----|-------------|
| FR-1 | Solve reflectance/transmittance for an arbitrary stack of uniform, dispersive, isotropic layers at arbitrary incidence angle/azimuth/polarization. *(done)* |
| FR-2 | Support semi-infinite incidence/transmission half-spaces of arbitrary (possibly complex/absorbing) index. *(done)* |
| FR-3 | Report Jones and Mueller-matrix polarimetric response. *(done)* |
| FR-4 | Represent 2D-periodic in-plane patterns from `Circle` and `Rectangle` primitives, including nested/overlapping shapes with correct area subtraction. *(geometry done; not yet consumed by the solver)* |
| FR-5 | Solve reflectance/transmittance/diffraction efficiencies for a layer patterned according to FR-4 (via/pillar). *(planned — Phase 4)* |
| FR-6 | Represent and solve 1D-periodic lamellar (line/space) patterns (trench). *(planned — Phase 3)* |
| FR-7 | Represent a feature (via/trench) with linearly tapered sidewalls via staircase layer discretization, and demonstrate R/T convergence with slice count. *(planned — Phase 5)* |
| FR-8 | Support anisotropic (full 3×3 tensor) materials in both uniform and patterned layers. *(planned — Phase 6)* |
| FR-9 | Reconstruct real-space E/H field maps at an arbitrary depth in the stack (for cross-section visualization of trench/via structures). *(planned — Phase 7)* |
| FR-10 | Ingest dispersive material data from refractiveindex.info-style CSV `n,k` exports. *(done — `structures/thin_film/sio2_on_si_thin_film.py::material_from_csv`)* |

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
| General (non-uniform) complex eigenproblem can have degenerate/near-degenerate eigenvalues causing numerical instability | Reuse the already-validated `_select_q_branch` outgoing-mode convention; add targeted regression tests for near-degenerate cases once Phase 4 lands |
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
  method; this is not a general-purpose EM solver.
- Nonlinear optics, magneto-optic materials, thermal/mechanical coupling.
- A GUI or web interface — this is a Python library driven by scripts.
- Public package distribution (PyPI) / multi-user support — see
  `deployment.md` for current solo-research-tool scope.
