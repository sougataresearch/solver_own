# Project Memory — sougata_solver

Living document for future sessions (AI or human). Update this at the end
of every substantive session — see `rules.md`'s AI Coding Rules, item 6.

## Current Project Status

As of 2026-07-24 (Phase 5), 2026-07-23 (Phase 4b), 2026-07-21 (Phase 4a) and earlier entries below:
- **Phase 5 (tapered/sloped sidewalls, staircase discretization) is
  complete.** New module `src/sougata_solver/staircase.py` with three
  generator functions (`staircase_circle_layers`, `staircase_rectangle_layers`,
  `staircase_slab_layers`), each producing `num_slices` uniform-in-z
  `Layer`s whose shape size is linearly interpolated (at each slice's
  z-midpoint) between a `top` and `bottom` value, per `decisions.md`
  ADR-004. Per the `phase-reference-picker` skill: this phase has **no
  citable source** in any vendored `REFERENCE/` repo — a grep for
  "stair"/"taper" across the whole RCWA family (`S4`, `EMpy`,
  `RigorousCoupledWaveAnalysis.jl`, `Rigorous-Coupled-Wave-Analysis`) found
  nothing; `staircase.py` is independently derived (the technique itself is
  standard/well-precedented, per ADR-004), flagged per `rules.md` AI Coding
  Rule 1, and its correctness rests on convergence-vs-`num_slices` evidence
  rather than an external oracle-comparison test — consistent with
  `phases.md` Phase 5's "no new Fourier/eigenmode math" scoping (every
  slice reuses Phase 3/4a's already-oracle-validated per-layer solve
  unchanged). Validated: (a) zero-taper (`top==bottom`) staircase
  reproduces the already-validated single-layer uniform-cross-section
  result to `1e-10` for `Circle`, `Rectangle`, and `Slab`, regardless of
  `num_slices` (regression guard); (b) energy conservation holds for
  actually-tapered cases; (c) two `slow`-marked convergence-vs-`num_slices`
  studies (`tests/test_staircase.py`), sweeping `num_slices` from 1 to 32
  (via) / 64 (trench) — both show clean monotone-shrinking successive
  differences (via: R settles to ≈0.565 by N=16-32; trench: R settles to
  ≈0.248 by N=32-64, needing more slices to converge than the via case,
  consistent with its larger relative taper). `structures/via/tapered_via.py`
  and `structures/trench/tapered_trench.py` (per ADR-010, print/save only,
  no plotting) both run end-to-end with R+T=1.0000 at every `num_slices`
  value. 130 tests pass project-wide (118 prior + 12 new — 10 fast + 2
  `slow` convergence studies — `tests/test_staircase.py`).
  **Follow-up, 2026-08-01**: at user request, renamed the tapered
  trench/via scripts' geometry constants to FDTD-style `TCD`/`BCD`/`SPACING`
  (matching an equivalent Lumerical grating structure group the user
  shared) and added `structures/via/tapered_pillar.py` (`Rectangle` case of
  `staircase_rectangle_layers`, previously the only shape type without an
  example script) — no new physics, same staircase generators; all three
  re-verified with R+T=1.0000.
- **Phase 4a (2D-periodic patterned layers, well-conditioned case) is complete
  and validated.** `solve_layer_eigenmodes_patterned` (the dense `2n x 2n`
  general eigenoperator from `S4/S4/rcwa.cpp::SolveLayerEigensystem` lines
  794-827) is implemented in `eigenmodes.py`; `simulation.py` dispatches to it
  for any 2D-patterned `Layer`. **Correction made this session** to a first
  draft (written by a different agent) that copied `solve_layer_eigenmodes_1d`'s
  `Epsilon2` construction (`epsilon_hat` top-left, `inv(epsilon_inv_hat)`
  bottom-right) into the 2D solver — this looked plausible but was never
  actually checked against the relevant S4 source. Reading
  `S4/S4/fmm/fmm_closed.cpp:109-139` in full (specifically the branch after
  line 133, which nobody had read before) shows that formula is the **1D
  special case only** (`0==Lr[2]&&Lr[3]==0`, S4's own "1D proper FFF rule"
  comment). The true-2D, no-polarization-basis path S4 actually uses by
  default is plain **Laurent's rule throughout**: `Epsilon2 =
  block_diag(epsilon_hat, epsilon_hat)` (both blocks direct-rule), and the
  matrix fed into `kp` is `inv(epsilon_hat)` (the numerical matrix-inverse
  of the direct-rule Toeplitz) — Phase 2's separately-factorized
  `epsilon_inv_hat` (`toeplitz_matrix(..., inverse=True)`) is **not
  consumed by Phase 4a at all**; it remains 1D-only (Phase 3) infrastructure.
  Fixed `solve_layer_eigenmodes_patterned`'s signature to drop the
  now-unused `epsilon_inv_hat` parameter accordingly; see the function's
  docstring for the full line-by-line citation. This is a real, known
  RCWA accuracy limitation (ordinary Laurent's rule converges slower at 2D
  discontinuities than a proper vectorial/normal-vector factorization),
  not a design choice made here — a genuine 2D Li's-rule-equivalent
  (`use_polarization_basis` in S4) is out of Phase 4a's scope and would be
  a separately-requested extension. The bug had gone undetected because the
  test suite's own "ky=0 reduces to the 1D solver" check was circular (both
  solvers used the same, wrongly-copied formula, so of course they agreed)
  — replaced with two honest tests
  (`tests/test_2d_pillar.py::test_2d_patterned_ky_zero_te_block_matches_1d`,
  which *should* agree since the TE block never depends on which
  Fourier-factorization rule is used, and
  `..._tm_block_differs_from_1d`, which now correctly asserts the two
  *diverge*, guarding against this exact regression reappearing).
  Energy conservation and the `structures/via/pillar_array.py`/`via_array.py`
  R+T=1.0000 result both still hold after the fix (re-verified). Coupling
  subtraction: `[[diag(kx^2), diag(kx*ky)], [diag(kx*ky), diag(ky^2)]]`;
  `q = branch_select(eig(op))` reuses `_select_q_branch` unmodified. S4 was
  not runnable in this environment for a subprocess cross-check oracle (per
  `tasks.md`'s ☑ item on oracle determination); per `rules.md` AI Coding
  Rule 5 (never fabricate a match), the test strategy conservatively
  relies on: (a) the `ky=0` TE-block cross-check against the
  oracle-validated 1D solver, (b) the fully-uniform reduction to the Phase 1
  result, (c) energy conservation (R+T+sum(DE)=1 to 1e-8) across
  moderate-contrast Circle+Rectangle patterns at normal and oblique
  incidence, (d) both `structures/via/pillar_array.py` and
  `structures/via/via_array.py` running end-to-end with R+T=1.0000 across a
  21-point sweep from 500 nm to 1500 nm.

  **Follow-up, same day (user asked explicitly: stop relying on S4 alone,
  use the other vendored repos too)**: surveyed `EMpy`,
  `Rigorous-Coupled-Wave-Analysis`, and `RigorousCoupledWaveAnalysis.jl`
  for anything applicable to Phase 4a. `EMpy/EMpy/RCWA.py` is 1D-only (no
  2D support) and has an author-acknowledged instability hack — ruled out.
  `Rigorous-Coupled-Wave-Analysis`'s `run_RCWA_2D` also uses plain
  Laurent's rule for 2D (third independent confirmation the fix's rule
  choice isn't S4-specific), but has no hard-coded benchmark numbers.
  `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:57-99`
  (`eigenmodes(...,l::PatternedLayer)`) does — and its formula is a
  **structurally different derivation** (direct Maxwell-curl elimination
  into one matrix, not S4's `Epsilon2 @ kp` route). Julia isn't installed
  here (`which julia` fails), so hand-transcribed into
  `tests/oracles/rcwa_2djl_eigenvalues.py`, feeding it this project's own
  already-validated `epsilon_hat` to isolate the eigenoperator-construction
  step specifically. After reconciling two conventions (both confirmed
  empirically before trusting the comparison, not assumed): RCWA.jl
  normalizes `kx,ky` by `k0`, and its eigenvalues come out an overall sign
  flip of this project's `q^2` (opposite time convention) — its `q^2`
  values match `solve_layer_eigenmodes_patterned`'s to **~1e-12** across 6
  parametrized `num_orders`/angle cases
  (`tests/test_2d_pillar.py::test_2d_patterned_eigenvalues_match_rcwa_jl_oracle`).
  This is a genuinely independent formula agreeing with the fixed
  implementation — meaningfully stronger evidence than re-reading S4 a
  second time, and exactly the kind of check that would have caught the
  original bug (the old circular "ky=0 reduces to 1D" test could not).
  **Re-audited Phase 1-3's S4 citations too**, specifically re-reading each
  cited function's *complete* body (not just the branch matching prior
  usage) given that's exactly how the Phase 4a bug slipped through —
  `MakeKPMatrix`, `SolveLayerEigensystem_uniform` (Phase 1),
  `pattern_get_fourier_transform` (Phase 2), the 1D G-vector-selection and
  `fmm_closed.cpp` branches (Phase 3): all confirmed to match the existing
  implementation exactly, no further bugs found (worth noting: S4's
  `pattern.c` has its own explicit 1D-`RECTANGLE` case that independently
  corroborates `Slab`'s formula, a nice confirmation found during the
  re-read). A full external **R/T** oracle (not just eigenvalues) is still
  missing — no independently-published 2D benchmark exists in any vendored
  repo, S4 isn't buildable here, Julia isn't installed either — explicitly
  flagged, not faked, in `tests/oracles/rcwa_2d_pillar.py`; real remaining
  risk carried to Phase 4b, not just its originally-scoped near-degenerate
  stress cases. `simulation.py`'s `NotImplementedError` for 2D patterns is
  fully removed. 107 tests pass project-wide (101 + 6 new parametrized
  eigenoperator-oracle cases).

- **Phase 4b (near-degenerate/ill-conditioned 2D cases) is complete.**
  Stress-tested `solve_layer_eigenmodes_patterned` deliberately: index
  contrast from `n=3.48` to a lossy-metal-like `-20+2j`, `num_orders` up to
  225, near-touching circular pillars (`radius=0.49*period`), a
  sub-percent-halfwidth sliver rectangle, and near-degenerate nested
  circles (`1e-4`-scale radius difference). **Honest finding**: no
  catastrophic failure — `cond(epsilon_hat)` reached ~900, `cond(phi)`
  reached ~170 in the worst cases, energy conservation and the Phase 4a
  `RigorousCoupledWaveAnalysis.jl` eigenvalue oracle both held to ~1e-10
  throughout (`tests/test_2d_pillar_stress.py`). Reported as what was
  actually observed for the closed-form isotropic `Circle`/`Rectangle`
  patterns tested — not oversold as "no pathological case can exist
  anywhere." Since no failure needed fixing, the deliverable became
  **detection**: `eigenmodes.ILL_CONDITIONED_THRESHOLD = 1e4` (~10x
  headroom above the worst case observed) triggers a `WARNING` via a new
  module-level `logger` in `eigenmodes.py` (per `design.md`'s Logging
  Strategy) whenever `cond(epsilon_hat)` or `cond(phi)` exceeds it —
  verified by dedicated tests using `caplog`/`monkeypatch` (one confirming
  silence in the ordinary case, one forcing the threshold down to confirm
  the mechanism actually fires). The "match S4 or a published benchmark"
  deliverable could **not** be met as worded — confirmed (not assumed) S4
  needs `cmake`+Lua (neither present) and Julia isn't installed either
  (`which julia` fails) — so `tasks.md` records that box honestly unchecked
  rather than silently substituting the RCWA.jl oracle and calling it done;
  the substitution is used but disclosed. 118 tests pass project-wide (107
  + 11 new: 9 stress-case oracle cross-checks + 2 logging tests).
  `troubleshooting.md` updated with the Phase 4a bug (moved to
  Already-Solved Gotchas) and this phase's stress-test findings.
  **Phase 5 (tapered sidewalls) is next**, or Phase 6 (anisotropic
  materials) — both depend only on Phase 4a, not on closing Phase 4b's
  remaining external-R/T-oracle gap, which is now a standing item (not
  tied to a specific phase number) rather than blocking anything.

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
- **FDTD is a genuine future goal (ADR-011, 2026-07-21)**, but a separate
  future effort, not a `sougata_solver` phase — `sougata_solver` stays
  RCWA-only (frequency-domain, periodic-BC); `PRD.md`'s Out-of-Scope
  wording was corrected after it was found to read as "FDTD rejected
  outright," which contradicted the actual intent. No FDTD design/timeline
  exists yet. `REFERENCE/meep`, `gprMax`, `fd3d`, `maxwellfdfd` (FDTD/FDFD)
  and `mfem`, `OpenParEM`, the FEniCS stack (FEM) are vendored for that
  future effort and remain unused/unevaluated by any current RCWA phase.

## Completed Milestones

- Phase 1: uniform multilayer stacks, dispersive materials, arbitrary
  polarization/angle, Jones/Mueller polarimetry — validated against
  analytic Fresnel/TMM (`tests/test_analytic_fresnel.py`).
- Phase 2: `fourier_factorization.py` (`pattern_epsilon_hat`,
  `toeplitz_matrix`) — validated against a from-scratch rasterize-and-sum
  reference (`tests/test_fourier_factorization.py`).
- Phase 3: `Lattice1D`/`Slab` (`geometry.py`), `truncate_fourier_orders_1d`
  (`fourier_basis.py`), `solve_layer_eigenmodes_1d` (`eigenmodes.py`),
  `Lattice1D` dispatch + `diffraction_efficiencies()` in `simulation.py` —
  validated against `tests/oracles/rcwa_1d_gaylord.py` and the
  energy-conservation/reduces-to-uniform invariants
  (`tests/test_1d_grating.py`).
- Phase 4a: `solve_layer_eigenmodes_patterned` (`eigenmodes.py`, general
  dense `2n x 2n` eigensolver, ordinary Laurent's rule per S4's true-2D
  closed-form path — corrected mid-session from an initially-wrong 1D-rule
  copy), 2D-patterned dispatch in `simulation.py`,
  `structures/via/pillar_array.py`, `structures/via/via_array.py` —
  validated by ky=0→1D TE-block reduction, fully-uniform reduction, energy
  conservation across Circle+Rectangle patterns, 21-point wavelength sweep
  producing R+T=1.0000 at each point, and an independent eigenvalue oracle
  transcribed from `RigorousCoupledWaveAnalysis.jl`
  (`tests/oracles/rcwa_2djl_eigenvalues.py`, agrees to ~1e-12)
  (`tests/test_2d_pillar.py`, 16 tests).
- Phase 4b: condition-number `WARNING` logging (`eigenmodes.ILL_CONDITIONED_THRESHOLD`)
  for `epsilon_hat`/`phi` in `solve_layer_eigenmodes_patterned` — validated
  by a deliberate near-degenerate/ill-conditioned stress sweep (no
  catastrophic failure found; energy conservation and the RCWA.jl oracle
  both held to ~1e-10) plus dedicated logging-mechanism tests
  (`tests/test_2d_pillar_stress.py`, 11 tests).
- Phase 5: `staircase.py` (`staircase_circle_layers`,
  `staircase_rectangle_layers`, `staircase_slab_layers`) — independently
  derived (no vendored-repo source exists for this discretization
  technique), validated by zero-taper regression to Phase 3/4a's
  single-layer results, energy conservation, and convergence-vs-`num_slices`
  studies (`tests/test_staircase.py`, 12 tests);
  `structures/via/tapered_via.py`, `structures/trench/tapered_trench.py`.

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
(Phase 5 is now done; Phase 6 — anisotropic materials — is the natural next
phase, and Phase 7 — field reconstruction — is an equally valid
alternative, both depending only on Phase 4a):
1. Phase 6: generalize Phase 4a's eigensolver to a full 3x3 tensor
   `Epsilon2`, source a birefringent-material closed-form benchmark
   (uniaxial waveplate at normal incidence).
2. Phase 7: real-space field reconstruction (`fields.py` extension using
   `SMatrixStack.partial_smatrix_up_to`) and cross-section plotting
   examples for trench/via.
3. Standing item, not tied to a specific phase: a full external **R/T**
   oracle for 2D patterns (S4 subprocess, a published benchmark, or
   transcribing `Rigorous-Coupled-Wave-Analysis`'s `run_RCWA_2D`) — still
   open per Phase 4a/4b's honest accounting above; revisit if/when S4 or
   Julia becomes runnable in this environment, or a literature benchmark
   is located.

## Architecture Notes

- `SMatrixStack` and the rest of `smatrix.py` (plus `fields.py`,
  `excitation.py`) were confirmed **genuinely dimension-agnostic** by
  Phase 3 landing with **zero changes** to any of the three — they operate
  purely on `(q, phi, kp, thickness)` regardless of whether a layer is
  uniform or 1D-patterned. Phase 4a should need the same: only
  `eigenmodes.py`/`simulation.py` need new code for 2D patterns.
- `eigenmodes.build_kp_matrix` **already accepts** a full `(n,n)`
  `epsilon_inv` matrix, not just a scalar (see the `else` branch,
  `eigenmodes.py:40-47`) — confirmed in Phase 3 to be exactly the right
  interface: `solve_layer_eigenmodes_1d` calls it unmodified with Phase 2's
  `epsilon_inv_hat` Toeplitz. Phase 4a's `solve_layer_eigenmodes_patterned`
  should reuse it the same way.
- `SimulationResult.diffraction_efficiencies()` (new, Phase 3) computes
  per-order R/T by calling the existing `z_poynting_flux` once per order
  with other orders' amplitudes masked to zero, rather than a new per-order
  flux formula — valid because `all_modes[0]`/`all_modes[-1]` are *always*
  the uniform incidence/transmission half-spaces (`LayerStack.__init__`)
  regardless of what's patterned in between, and `build_kp_matrix`'s
  `epsilon_inv` branch has no cross-order coupling (each order's `kappa`
  entries only touch that order's own index and its `+n` pair). Phase 4a's
  2D case can reuse this unmodified too.
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
- **Phase 4 was split into Phase 4a (well-conditioned via/pillar case) and
  Phase 4b (near-degenerate/ill-conditioned stress cases), 2026-07-19**,
  at user request for more mathematical/scientific rigor in the roadmap —
  see `phases.md`, `tasks.md`, `PRD.md`. Phase 5/6/7 now depend on 4a only
  (4b hardens the solver in place, it doesn't unlock new capability).
  `testing.md` also gained a new Physical-Invariant Testing tier (energy
  conservation, convergence-rate-vs-Li-1996-theory), required starting
  Phase 3, oracle-independent and cheaper than the S4/benchmark
  cross-check — see `rules.md` Testing Requirements.
- **A `phase-reference-picker` skill now exists**
  (`.claude/skills/phase-reference-picker/`, at the `Solver_own` workspace
  root, sibling to `sougata_solver` and `REFERENCE`), added 2026-07-19.
  Invoke it before writing physics code for any new phase — it forces a
  real per-sub-task comparison across `REFERENCE/`'s ~18 vendored repos
  (not just S4, which is easy to default to since it's the Phase 1/2
  oracle already) and an explicit transcribe-vs-derive-independently
  decision. See `references.md`'s "Choosing a Reference for a New Phase"
  and `rules.md` AI Coding Rule 8.
