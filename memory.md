# Project Memory — pyrcwa

Living document for future sessions (AI or human). Update this at the end
of every substantive session — see `rules.md`'s AI Coding Rules, item 6.

## Current Project Status

As of 2026-07-15:
- **Phase 1 (uniform multilayer core) is complete and validated.**
  Reflectance/transmittance for arbitrary uniform-layer stacks, arbitrary
  incidence angle/polarization, dispersive materials, and Jones/Mueller
  polarimetry all work and are tested against analytic Fresnel/TMM.
- **Phases 2-9 are planned but not started** (see `phases.md`, `tasks.md`).
  `simulation.py:98` explicitly raises `NotImplementedError` for any
  patterned layer — this is the immediate next blocker for trench/via/pillar.
- `pyrcwa` was just initialized as its own git repository (previously
  un-versioned) with a `.gitignore` covering `__pycache__/`,
  `.pytest_cache/`, `*.egg-info/`, `*.csv`, `*.png`.
- `examples/` was removed and replaced with `structures/` (build a
  lattice/layer stack/materials, run the solver) and `postprocessing/`
  (derive Jones/Mueller matrices, ellipsometric angles, and — planned —
  RI/thickness extraction, from a `structures/` script's raw output). See
  `decisions.md` ADR-009. `polarimetry.py`'s `_decompose_sp` was made
  public (`decompose_sp`) so `postprocessing/jones_mueller_ellipsometry.py`
  can reuse the solver's exact s/p convention.
- The full documentation set (this file, `README.md`, `PRD.md`,
  `architecture.md`, `design.md`, `rules.md`, `phases.md`, `tasks.md`,
  `decisions.md`, `testing.md`, `deployment.md`, `references.md`,
  `troubleshooting.md`) was just created in this session, before any
  Phase 2+ code was written — this is documentation-first, code-second by
  explicit user instruction.

## Important Decisions

Full rationale lives in `decisions.md` (ADR format). Summary:
- S-matrix (Redheffer star product), not transfer-matrix cascading —
  numerical stability for evanescent modes.
- Analytic shape Fourier transforms (S4-style), not raster+FFT
  (Meent/TORCWA-style) — accuracy for smooth boundaries.
- 1D lamellar gratings (Phase 3) come before 2D general patterned layers
  (Phase 4) — lower risk, validates Fourier-factorization machinery on the
  simpler (decoupled TE/TM) case first.
- Tapered sidewalls via staircase discretization, not new Fourier math —
  cheap, well-precedented (matches how even FEM tools like the vendored
  JCMsuite tutorials handle the same geometry).
- No arbitrary-polygon/GDS geometry support — parametric shapes
  (`Circle`, `Rectangle`, future `Slab`) only, by explicit user choice.
- GPU/autodiff backend (Meent/TORCWA-style) explicitly deferred to
  optional Phase 9, after Phases 2-8 are validated — by explicit user
  choice, so correctness work isn't chasing a moving numeric backend.
- Target audience is a solo research tool for now, not a public package —
  by explicit user choice; affects `deployment.md`'s scope (no PyPI/Docker
  yet).

## Completed Milestones

- Phase 1: uniform multilayer stacks, dispersive materials, arbitrary
  polarization/angle, Jones/Mueller polarimetry — validated against
  analytic Fresnel/TMM (`tests/test_analytic_fresnel.py`).
- Geometry primitives (`Circle`, `Rectangle`, `Pattern` with containment
  tree) implemented and unit-tested at the shape level
  (`tests/test_fourier_factorization.py`), but **not yet wired into the
  solver** — this is Phase 2/4's job, not redundant with what's tested today.
- Full project documentation suite created (this session).

## Known Issues

- No `.flake8`/`ruff` or `mypy.ini` config exists yet in `pyrcwa/` (unlike
  the sibling `EMpy` reference project, which has both) — flagged as a
  Phase 2 prerequisite task in `tasks.md`.
- `excitation.py`'s s/p polarization sign convention is explicitly
  documented as "internal and self-consistent (not yet matched to
  S4/EMpy's convention)" (`excitation.py:16-19`) — fine for Phase 1 (only
  power quantities are validated), but should be revisited if/when a
  polarization-sensitive cross-check against S4/EMpy is needed in a later
  phase.
- `scripts/` directory exists but is currently empty — no ad hoc utility
  scripts have been added yet.
- Whether S4 is actually built/runnable in this development environment
  (needed for Phase 4's cross-check oracle) has **not been checked yet** —
  first task when Phase 4 starts.

## Pending Tasks

See `tasks.md` for the full atomic checklist. Immediate next actions
(Phase 2 start):
1. Add lint/type-check config to `pyrcwa/`.
2. Create `src/pyrcwa/fourier_factorization.py` with
   `pattern_epsilon_hat`/`toeplitz_matrix`.
3. Validate against FFT-of-rasterized-mask numerically.

## Architecture Notes

- `SMatrixStack` and the rest of `smatrix.py` are **dimension-agnostic
  already** — they operate purely on `(q, phi, kp, thickness)`. Phases 3/4
  only need to *produce* `LayerEigenmodes` correctly; no changes to
  `smatrix.py` are anticipated.
- `eigenmodes.build_kp_matrix` **already accepts** a full `(n,n)`
  `epsilon_inv` matrix, not just a scalar (see the `else` branch,
  `eigenmodes.py:40-47`) — this was clearly written in anticipation of
  Phase 2/4's patterned-layer case, so Phase 2's Toeplitz output should
  target this exact interface, not a new one.
- See `architecture.md` for the full module responsibility table and data
  flow diagram.

## Technical Debt

- `excitation.py`'s polarization sign convention (see Known Issues above)
  is technical debt that's currently harmless but will need reconciling
  with S4's convention before any polarization-sensitive cross-check test
  can be written (relevant from Phase 4 onward if polarimetric, not just
  power, validation is needed against S4).
- No CI pipeline exists yet (git repo is brand new) — acceptable at current
  solo/local scope per `PRD.md`, but see `deployment.md` for the plan once
  a remote exists.

## Things Future AI Sessions Should Remember

- **Read `rules.md`'s "AI Coding Rules" section before writing any new
  physics code** — the project's entire trust model rests on never
  presenting an unverified formula as verified, and never fabricating a
  benchmark match.
- **This project's plan-mode scratch file**
  (`C:\Users\d14k4\.claude\plans\vivid-swimming-moler.md`) contains the
  originally-approved Phase 2-9 plan; `phases.md`/`tasks.md` in this repo
  are now the authoritative, living versions — keep them in sync going
  forward, don't treat the scratch plan file as the source of truth after
  this point.
- **Do not import from or modify the vendored reference repos**
  (`../S4`, `../EMpy`, `../RigorousCoupledWaveAnalysis.jl`, `../EMTutorial`)
  — they are read-only oracles.
- **The user's separate Claude-Code memory system**
  (`C:\Users\d14k4\.claude\projects\...\memory\`) is a different mechanism
  from this file — that one is cross-project and cross-session for the AI
  assistant's own use; this `memory.md` is project-scoped documentation
  living inside the `pyrcwa` repo itself, readable by any collaborator or
  future session regardless of which AI tool is used.
