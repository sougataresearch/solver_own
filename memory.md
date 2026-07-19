# Project Memory — sougata_solver

Living document for future sessions (AI or human). Update this at the end
of every substantive session — see `rules.md`'s AI Coding Rules, item 6.

## Current Project Status

As of 2026-07-18:
- **Phase 2 (Fourier-factorization core) is complete and validated.**
  `fourier_factorization.py`'s `pattern_epsilon_hat`/`toeplitz_matrix`
  build the direct and inverse-rule Toeplitz permittivity matrices from a
  `Pattern`, transcribed from `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`
  (lines 889-1029) and `S4/S4/fmm/fmm_closed.cpp::FMMGetEpsilon_ClosedForm`
  (lines 77-127). Validated against **two independent** numerical
  references (neither calling into the module under test): a from-scratch
  rasterize-and-sum, and a literal FFT-of-rasterized-mask reproduction of
  the vendored `Rigorous-Coupled-Wave-Analysis` (Python `convmat2D.py`)
  and `RigorousCoupledWaveAnalysis.jl` (`ft2d.jl::real2recip`)
  convolution-matrix algorithm — the user explicitly asked that S4 not be
  the only cross-check source, and `rules.md` names RCWA.jl as a
  sanctioned oracle. Caught one real bug during this: the first
  FFT-reference attempt used an uncentered raster grid that silently
  truncated a shape whose footprint crossed the domain edge; fixed via
  `numpy.fft.ifftshift` before the FFT (see `phases.md` Phase 2 for the
  full story). Both references agree with the analytic Toeplitz entries
  for `Circle` and `Rectangle` patterns, direct and inverse rule, at
  several nonzero G-vectors (`tests/test_fourier_factorization.py`, 12
  tests, all passing; 75 tests pass project-wide). Scalar isotropic
  materials only; anisotropic materials raise `NotImplementedError`
  naming Phase 6.
  `.flake8`/`mypy.ini` added (mypy/flake8 themselves are not installed in
  this dev environment, so linting was done by manual review, not an
  actual tool run — flag this for whoever next has the tools available).
  **Phase 3 (1D lamellar gratings) is next.**
- **Phase 1 (uniform multilayer core) is complete and validated.**
  Reflectance/transmittance for arbitrary uniform-layer stacks, arbitrary
  incidence angle/polarization, dispersive materials, and Jones/Mueller
  polarimetry all work. Validated by **two independent oracles**:
  `tests/oracles/fresnel.py` (from-scratch analytic Fresnel/TMM) and
  `tests/oracles/empy_tmm.py` (transcribed from the vendored EMpy
  reference library), cross-checked against the actual SiO2-on-Si
  structure in `tests/test_thin_film_empy_cross_check.py` across
  wavelength/angle/polarization, agreeing to `1e-8`. 68 tests pass.
- **Phases 2-9 are planned but not started** (see `phases.md`, `tasks.md`).
  `simulation.py:98` explicitly raises `NotImplementedError` for any
  patterned layer — this is the immediate next blocker for trench/via/pillar.
- The package was renamed `pyrcwa` → `sougata_solver` (all imports/docs
  updated) and is its own git repository, with a `.gitignore` covering
  `__pycache__/`, `.pytest_cache/`, `*.egg-info/`, `*.csv`, `*.png`,
  `outputs/`.
- `examples/` was removed and replaced with `structures/` (build a
  lattice/layer stack/materials, run the solver) and `postprocessing/`
  (derive Jones/Mueller matrices, ellipsometric angles, plots, and —
  planned — RI/thickness extraction, from a `structures/` script's raw
  output; never calls `Simulation.solve`). See `decisions.md` ADR-009.
  `structures/` is further grouped by category (`structures/thin_film/`
  today; `structures/trench/`, `structures/via/` etc. once those phases
  land). `polarimetry.py`'s `_decompose_sp` was made public (`decompose_sp`)
  so `postprocessing/jones_mueller_ellipsometry.py` can reuse the solver's
  exact s/p convention.
- **Every run gets its own timestamped output folder** (`output_paths.py`:
  `outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/`) containing its raw CSV/data
  *and* a `run_metadata.txt` (`write_run_metadata`) recording which script
  produced it and its key parameters — so re-running the same script with
  different settings never collides or gets mixed up (ADR-010). Plotting
  is a `postprocessing/` script (`plot_thin_film_rt.py`) that finds the
  relevant run's CSV (`find_latest_output`, or an explicit path) and saves
  its PNG back into that same run folder — plotting was briefly added
  directly to a `structures/` script and the user correctly caught that as
  a boundary violation; see ADR-010 for the fix and the reasoning.
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
- Phase 2: `fourier_factorization.py` (`pattern_epsilon_hat`,
  `toeplitz_matrix`) — validated against a from-scratch rasterize-and-sum
  reference (`tests/test_fourier_factorization.py`). **Not yet wired into
  the solver** — `simulation.py:98`'s `NotImplementedError` for patterned
  layers still stands; that's Phase 3/4's job (they *consume*
  `toeplitz_matrix`'s output as the `epsilon_inv` argument
  `eigenmodes.build_kp_matrix` already accepts).
- Full project documentation suite created (2026-07-16 session).

## Known Issues

- `.flake8`/`mypy.ini` now exist (added with Phase 2), but neither `flake8`
  nor `mypy` is actually installed in this dev environment — no run has
  verified the codebase is clean under either tool yet.
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
(Phase 3 start — 1D-periodic lamellar gratings):
1. Add `Lattice1D` and a `Slab`/`Line` shape to `geometry.py`.
2. Add `truncate_fourier_orders_1d` to `fourier_basis.py`.
3. Add `solve_layer_eigenmodes_1d(...)` to `eigenmodes.py` (decoupled
   scalar TE/TM eigenproblems), consuming Phase 2's `toeplitz_matrix`.
4. Wire a `Lattice1D` dispatch branch into `simulation.py`.
5. Validate against a published 1D binary-grating benchmark (Moharam &
   Gaylord 1995 or equivalent).
6. `structures/trench/trench_grating.py`.

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
- **`progress_log.md` (new, 2026-07-19)** is a dated, append-only log of
  discussions and their action items — check it at the start of any
  session for open `[ ]` items, verify against the actual code whether
  they've since been implemented, and add a new dated entry at the end of
  any substantive session. Distinct from this file (status snapshot) and
  `tasks.md` (phase-organized checklist).
- **The user's separate Claude-Code memory system**
  (`C:\Users\d14k4\.claude\projects\...\memory\`) is a different mechanism
  from this file — that one is cross-project and cross-session for the AI
  assistant's own use; this `memory.md` is project-scoped documentation
  living inside the `sougata_solver` repo itself, readable by any collaborator or
  future session regardless of which AI tool is used.
