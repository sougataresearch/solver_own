# Architecture Decision Record — sougata_solver

## ADR-001: S-matrix (Redheffer star product) over transfer-matrix cascading

- **Decision**: Represent every layer/interface as a scattering matrix and
  cascade with the Redheffer star product (`smatrix.py`), never as a
  transfer (ABCD-style) matrix multiplied end to end.
- **Reason**: Transfer matrices contain terms that grow exponentially for
  evanescent (decaying) modes through a thick or lossy layer; multiplying
  many such matrices together loses numerical precision catastrophically
  (the classic, well-documented instability of the "T-matrix method" in
  grating theory). S-matrices keep every intermediate quantity bounded.
- **Alternatives considered**: Transfer-matrix method (TMM) cascading —
  simpler to implement, directly matches how `EMpy`'s `transfer_matrix.py`
  reference module works for uniform stacks, but does not scale to
  patterned layers with many evanescent orders (Phase 3/4's whole point).
- **Trade-offs**: S-matrix cascading requires a matrix inversion (or LU
  solve) per interface (`_solve` in `smatrix.py`) instead of a plain
  matrix product; more arithmetic per interface, but the numerical
  stability gain is required for anything beyond Phase 1's thin, low-order
  case.
- **Impact**: Confirmed correct for Phase 1 (validated against Fresnel).
  `SMatrixStack` is dimension-agnostic and requires no changes for Phase
  3-6 — this decision made those later phases purely additive.

## ADR-002: Analytic shape Fourier transforms over raster+FFT of a pixelized mask

- **Decision**: In-plane patterns (`Circle`, `Rectangle`) expose a
  closed-form analytic Fourier transform (`jinc`/`sinc` — `geometry.py`),
  matching S4's approach, rather than rasterizing the pattern onto a pixel
  grid and taking a numerical FFT (the approach used by many RCWA
  implementations, including — per public documentation — Meent and
  TORCWA's typical workflow).
- **Reason**: Analytic transforms have no pixelization error for smooth
  boundaries (a circular via's edge is exact, not staircased at the mask
  resolution); this matters directly for the target use cases (circular
  vias, cylindrical pillars).
- **Alternatives considered**: Raster+FFT — simpler to extend to arbitrary
  shapes (any mask, including imported layouts), and is what Meent/TORCWA
  do, but was explicitly rejected for the primary shape library because it
  reintroduces a systematic error source Phase 2's Toeplitz construction
  doesn't otherwise have.
- **Trade-offs**: Every new shape type requires deriving/sourcing its own
  closed-form Fourier transform (more upfront math work per shape) instead
  of "just rasterize it" — acceptable because the shape library is
  deliberately small (`Circle`, `Rectangle`, and a future `Slab` for
  1D gratings), per the "Out-of-Scope" decision below.
- **Impact**: Directly shapes Phase 2's design (`pattern_epsilon_hat` sums
  analytic per-shape contributions, not an FFT of a rasterized array) and
  rules out arbitrary-polygon geometry without a separate, explicitly-scoped
  future decision to add raster+FFT support for that specific case.

## ADR-003: 1D lamellar gratings (Phase 3) before 2D general patterned layers (Phase 4)

- **Decision**: Implement and validate the trench (1D-periodic) case
  end-to-end before the via/pillar (2D-periodic general) case, even though
  both ultimately serve the same PRD goals.
- **Reason**: 1D gratings decouple TE and TM into independent scalar
  eigenproblems — no `Circle`-style 2D mode coupling — making them the
  lower-risk place to validate the brand-new Fourier-factorization (Phase
  2) and non-uniform-eigensolve code paths before attempting the harder,
  fully-general 2D eigenproblem (which is the single highest-risk
  remaining piece of the whole roadmap, per `phases.md` Phase 4).
- **Alternatives considered**: Go straight to 2D (via/pillar) since that's
  the geometry the vendored `Circle`/`Rectangle` shapes already support
  directly — rejected because it would conflate two new sources of risk
  (Fourier factorization + general eigendecomposition) in one
  unvalidated step, making failures harder to localize.
- **Trade-offs**: Requires building a parallel, 1D-specific code path
  (`Lattice1D`, `Slab`, `truncate_fourier_orders_1d`,
  `solve_layer_eigenmodes_1d`) that Phase 4 doesn't directly reuse (though
  it does reuse Phase 2's Fourier-factorization core and the S-matrix/field
  layers, which are dimension-agnostic).
- **Impact**: Phase ordering in `phases.md`/`tasks.md` reflects this;
  Phase 4 should be noticeably faster to validate than it would be as the
  first patterned-layer capability, because Fourier-factorization bugs
  will already have been shaken out in Phase 3.

## ADR-004: Tapered sidewalls via staircase discretization, not new Fourier math

- **Decision**: Represent a linearly-tapered via or trench sidewall as a
  stack of `N` thin layers, each with a slightly different (linearly
  interpolated) `Circle`/`Rectangle`/`Slab` size, rather than deriving a
  closed-form Fourier transform for a genuinely slanted 3D shape.
- **Reason**: This is the standard, well-precedented approach in RCWA
  (staircase/multi-slice approximation) — it requires zero new Fourier
  math, reuses Phase 3/4's per-layer solvers unchanged, and its accuracy
  is directly and cheaply verifiable via a convergence-vs-`N` study. It
  also mirrors how the vendored JCMsuite `ThroughSiliconVia` tutorial
  models the same physical geometry, just with FEM mesh refinement instead
  of layer-count refinement.
- **Alternatives considered**: A closed-form Fourier transform for a
  frustum/cone shape (would avoid discretization error entirely) — rejected
  as disproportionate effort for a solo-research-tool timeline, and
  because RCWA. This project explicitly favors the well-established,
  cheaply-validated approach over a more "elegant" one that would need its
  own from-scratch derivation and validation burden.
- **Trade-offs**: Convergence is not instantaneous — steep sidewall angles
  may require a non-trivial `N` to converge, and cost scales linearly with
  `N` (see `architecture.md`'s Scalability Considerations). Acceptable
  since `N` is a user-controlled knob, and the convergence study itself
  (Phase 5's deliverable) makes the trade-off visible rather than hidden.
- **Impact**: Phase 5 has near-zero dependency on new physics — it's
  "cheap" specifically because of this decision; do not revisit unless a
  specific structure is shown to need excessive `N` for acceptable accuracy.

## ADR-005: No arbitrary-polygon / GDS-imported geometry (parametric shapes only)

- **Decision**: The shape library stays limited to parametric primitives
  (`Circle`, `Rectangle`, and a planned `Slab` for 1D gratings) — no
  general polygon, GDS-import, or rasterized-arbitrary-mask support.
- **Reason**: Explicit user choice (asked directly, user selected "stick
  to parametric shapes"). Matches the PRD's target structures (thin film,
  multistack, trench, via, pillar), none of which need arbitrary polygons.
- **Alternatives considered**: Add polygon-vertex analytic Fourier
  transforms (still avoids raster+FFT error) or raster+FFT-of-mask support
  for GDS import — both explicitly deferred, not rejected outright; if a
  real need for imported layouts arises later, this ADR should be revisited
  with a new one, not silently overridden.
- **Trade-offs**: Cannot directly import real fab layouts; must approximate
  any non-primitive shape with the closest `Circle`/`Rectangle` combination.
- **Impact**: Keeps Phase 2's Fourier-factorization scope bounded to a
  small, closed set of analytic shape transforms — directly enables
  ADR-002's accuracy argument to hold without needing a fallback raster
  path for "everything else."

## ADR-006: GPU/autodiff backend deferred to optional Phase 9

- **Decision**: Do not pursue a torch/JAX backend, GPU batching, or
  autodiff-based inverse design (the headline features of Meent/TORCWA)
  until Phases 2-8 are complete and validated in pure NumPy/SciPy.
- **Reason**: Explicit user choice. Also a sound engineering call
  independent of preference: chasing a moving numeric backend while the
  physics itself is still being validated multiplies the surface area for
  bugs and makes it harder to tell whether a wrong answer is a physics bug
  or a backend-porting bug.
- **Alternatives considered**: Build the backend abstraction early (so
  later phases "just work" on GPU) — rejected because it's premature
  generalization against a requirement (GPU/autodiff) that isn't validated
  as needed yet, violating the project's own "no framework magic ahead of
  need" convention (see `rules.md`).
- **Trade-offs**: If/when Phase 9 is pursued, some Phase 2-8 code may need
  light refactoring to be backend-agnostic (e.g. avoiding NumPy-specific
  calls where a torch/JAX equivalent differs) — accepted as a reasonable
  future cost.
- **Impact**: `phases.md` Phase 9 explicitly requires a "decision
  checkpoint" task (re-confirm the backend is still wanted) before any
  work starts, rather than assuming it's still in scope by the time
  Phases 2-8 are done.

## ADR-007: Solo-research-tool scope for deployment/CI (no PyPI/Docker yet)

- **Decision**: `deployment.md` covers local environment setup, `pytest`,
  and (once useful) a simple GitHub Actions test-on-push workflow — not
  PyPI packaging, Docker, or production-server deployment.
- **Reason**: Explicit user choice; matches actual current usage (local
  scripts, solo developer, no external users yet).
- **Alternatives considered**: Full production-grade packaging pipeline
  from day one — rejected as unnecessary process overhead for the current
  scope; can be added later without redesigning anything, since it's purely
  additive tooling around an already-correct `pyproject.toml`-based package.
- **Trade-offs**: Revisit if/when this project is ever shared publicly or
  used by a second person.
- **Impact**: `deployment.md` stays intentionally light; `PRD.md`'s
  Out-of-Scope section records this explicitly so it isn't silently
  reintroduced as scope creep in a later phase.

## ADR-008: `sougata_solver` is its own git repository, separate from vendored reference repos

- **Decision**: `git init` was run inside `sougata_solver/` specifically (not at
  the `Solver_own/` parent level), giving `sougata_solver` its own history,
  independent of the already-git-versioned `S4`, `EMpy`, and
  `RigorousCoupledWaveAnalysis.jl` reference repos.
- **Reason**: `sougata_solver` is the user's own project; the reference repos are
  vendored, read-only oracles with their own independent upstream history
  that shouldn't be conflated with `sougata_solver`'s own commit history.
- **Alternatives considered**: A single `Solver_own`-level repo with
  everything inside — rejected because it would either require the
  reference repos to be submodules (added process overhead not currently
  justified) or would flatten their independent git history into one
  undifferentiated tree.
- **Trade-offs**: None significant at this scale — `sougata_solver/.gitignore`
  already excludes generated artifacts.
- **Impact**: `rules.md`'s Git Workflow section and `deployment.md`'s CI
  section describe `sougata_solver`'s own repo as the unit of versioning/CI, not
  `Solver_own` as a whole.

## ADR-009: Replace `examples/` with `structures/` + `postprocessing/`

- **Decision**: The generic `examples/` directory (with numeric-prefixed
  filenames like `01_fresnel_multilayer.py`) was removed entirely and
  replaced with two purpose-named directories: `structures/` (build a
  lattice/layer stack/materials and run the solver) and `postprocessing/`
  (derive Jones/Mueller matrices, ellipsometric angles, and — planned —
  RI/thickness extraction, from a `structures/` script's already-computed
  raw output). Files were renamed descriptively (e.g.
  `sio2_on_si_thin_film.py`, `custom_multistack.py`) instead of numbered.
- **Reason**: Explicit user request — `examples/` read as throwaway sample
  code rather than the actual day-to-day entry point of the project, and
  numeric filenames (`01_`, `02_`...) didn't communicate purpose. The
  deeper issue was conflating two genuinely different responsibilities in
  one file: `04_jones_mueller.py` built a stack, ran the solver, *and*
  computed a Jones/Mueller matrix all in one script, with no boundary
  between "run the physics" and "derive a quantity from the result."
- **Alternatives considered**: Keep `examples/` as one folder with
  better-named files only (no `postprocessing/` split) — rejected because
  it wouldn't address the user's specific ask that Jones/Mueller
  computation (and future RI/thickness extraction) live separately from
  structure-building/running code, and would leave `04_jones_mueller.py`'s
  build+run+analyze conflation in place.
- **Trade-offs**: Splitting the Jones/Mueller example required introducing
  a small raw-data interchange format (a CSV of per-polarization reflected
  `Ex, Ey` written by `structures/thin_film/sio2_on_si_ellipsometry_run.py` and read
  by `postprocessing/jones_mueller_ellipsometry.py`) that didn't exist
  before — more moving parts than one self-contained script, but it means
  `postprocessing/` scripts never need to call `Simulation.solve` at all,
  which is the actual property the user asked for. `polarimetry.py`'s
  internal `_decompose_sp` helper was made public (`decompose_sp`) so the
  postprocessing script reuses the solver's exact convention instead of
  duplicating that physics.
- **Impact**: Every doc that referenced `examples/NN_*.py` (`README.md`,
  `architecture.md`, `design.md`, `testing.md`, `phases.md`, `tasks.md`,
  `PRD.md`, `deployment.md`, `references.md`) was updated in the same pass
  — future phases' planned example scripts (Phase 3's `trench_grating.py`,
  Phase 4's `pillar_array.py`/`via_array.py`) now default into `structures/`
  by this same convention, and any future post-processing capability
  (RI/thickness extraction) defaults into `postprocessing/`.

## ADR-010: Plotting always lives in `postprocessing/`, never `structures/`; every run gets a `run_metadata.txt`

- **Decision**: Plotting code (matplotlib, R/T-vs-wavelength or any other
  derived view) is never added directly to a `structures/` script — it was
  briefly added to `structures/thin_film/sio2_on_si_thin_film.py` and the
  user correctly caught this as a violation of ADR-009's own boundary
  ("`structures/` never derives anything, only produces raw output").
  Plotting was moved to a new `postprocessing/plot_thin_film_rt.py`, which
  locates a `structures/` script's CSV via
  `output_paths.find_latest_output` (or an explicit path for a specific
  historical run), plots it, and saves the PNG back into that *same*
  `outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/` folder — never a new one.
  Separately, every `structures/` script that saves output now also calls
  a new `output_paths.write_run_metadata(output_dir, __file__, **params)`,
  writing a `run_metadata.txt` into that run's folder recording which
  script produced it and its key parameters (materials, thicknesses,
  angle, wavelength range, ...).
- **Reason**: Two related problems the user identified directly. First,
  plotting is a derived view of already-computed data, exactly like Jones/
  Mueller matrix construction (ADR-009) — it belongs in `postprocessing/`
  by the same logic, not bundled into the script that runs the physics.
  Second, once a `structures/` script is re-run repeatedly with different
  parameters (sweeping a thickness, trying a different material, etc.),
  its timestamped output folders become indistinguishable from each other
  by name alone — `run_metadata.txt` is the fix, making every run
  self-describing without needing to cross-reference back to whatever the
  script's `EDIT` blocks contained *at the time it was run* (which may have
  since been edited again).
- **Alternatives considered**: Logging run parameters into the CSV's
  header/filename instead of a separate metadata file — rejected as either
  cluttering the CSV (a filename encoding every parameter gets unwieldy
  past 2-3 knobs) or requiring the CSV format itself to change per script.
  A separate `run_metadata.txt` is uniform across every `structures/`
  script regardless of what it varies.
- **Trade-offs**: `output_paths.run_output_path` (single-call convenience:
  get a fresh run folder and one file path in it) is no longer sufficient
  for scripts that write both a CSV and a metadata file — those scripts
  now call `output_paths.run_output_dir` once and build both paths under
  it manually. `run_output_path` itself is kept (docstring updated to warn
  against calling it more than once per run) for any future script that
  only ever writes a single file.
- **Impact**: `structures/thin_film/sio2_on_si_thin_film.py`,
  `custom_multistack.py`, and `sio2_on_si_ellipsometry_run.py` all call
  `write_run_metadata`; `postprocessing/plot_thin_film_rt.py` is the first
  postprocessing script whose whole job is visualization rather than a new
  derived physical quantity. Any future field-visualization work (Phase 7)
  follows this same split: raw field data saved by `structures/`, plotted
  by a `postprocessing/` script into the same run folder.

## ADR-011: FDTD is a genuine future goal, but a separate effort from `sougata_solver`'s own phases

- **Decision**: `PRD.md`'s Out-of-Scope Items previously stated flatly that
  FDTD-style transient simulation is out of scope, reading as "this
  capability is rejected." Corrected: `sougata_solver` itself (this RCWA
  codebase, its module structure, its phase numbering) will not grow FDTD
  capability — RCWA is fundamentally frequency-domain/periodic-BC, so
  bolting time-domain simulation onto it doesn't fit the architecture. But
  per the project owner (2026-07-21), FDTD *is* a real future goal of the
  broader EM-wave-solver effort this project is part of; RCWA was chosen
  first specifically for simplicity (frequency-domain, periodic structures
  are the lower-risk starting point per this session's discussion). No
  FDTD phase, module, timeline, or even repository structure is decided
  yet — this ADR only records that the goal exists and clarifies
  `sougata_solver`'s own scope boundary, it does not commit to an FDTD
  design.
- **Reason**: The previous "out of scope" wording was discovered this
  session to contradict the project owner's actual intent when asked
  directly why the vendored FDTD/FEM reference repos (`meep`, `gprMax`,
  `fd3d`, `maxwellfdfd`, `mfem`, `OpenParEM`, the FEniCS stack) hadn't been
  surveyed — the honest answer combined two reasons, one still valid
  (different numerical method, not a formula source for any current RCWA
  phase, per the `phase-reference-picker` skill's own guidance) and one
  that was simply wrong (the docs saying FDTD is rejected outright, not
  "later, separately"). Docs should reflect actual intent, not create a
  false impression that a future goal was considered and declined.
- **Alternatives considered**: Leaving `PRD.md` unchanged and treating this
  as a verbal-only clarification — rejected because `memory.md`/`rules.md`
  both treat the written docs as the source of truth a future session (AI
  or human) inherits; an undocumented verbal clarification would be lost
  the moment this conversation ends, and the next session would again see
  "out of scope" and again not survey the FDTD/FEM repos, repeating the
  same gap.
- **Trade-offs**: None substantive — this is a documentation correction,
  not a scope commitment. It does not obligate any FDTD work, add a
  dependency, or change any current phase's deliverables.
- **Impact**: `PRD.md`'s Out-of-Scope Items entry rewritten to state the
  clarification and explicitly name which vendored repos are already
  sitting there for that future effort (`meep`, `gprMax`, `fd3d`,
  `maxwellfdfd` for FDTD/FDFD; `mfem`, `OpenParEM`, `dolfinx`/`ufl`/
  `basix`/`ffcx`, `FreeFem-sources` for FEM, in case that route is chosen
  instead of or alongside FDTD). No other file changed — a real FDTD
  planning session (reference survey, architecture decision on
  standalone-vs-shared-codebase, etc.) is future work, not done here.

## ADR-012: Fast Fourier Factorization (FFF) and the Normal Vector Method (NVM) — evaluated and deferred

- **Decision**: `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 3 targets 3.4
  (FFF) and 3.5 (NVM) are two separate feasibility-decision targets by
  design, but converge to the same answer here for the same underlying
  reason (below): **defer both**, do not implement either as part of
  target 3.6. Current solvers keep ordinary Laurent's-rule Toeplitz
  construction for every 2D branch (`solve_layer_eigenmodes_patterned`,
  `solve_layer_eigenmodes_patterned_inplane`) and Li's (1996) inverse-rule
  correction for the 1D branch only, exactly as already shipped and
  documented in `design.md`'s Fourier-Factorization Rule Inventory
  (Category 3 target 3.1).
- **Reason**: Both techniques exist specifically to fix a real, now
  *measured* (not hypothesized) limitation of this project's current 2D
  Fourier factorization — `tests/test_fourier_convergence.py`'s target-3.3
  high-contrast pillar fixture (`n=5`, `radius=0.2*period`) shows genuinely
  poor low-order behavior (an order-of-magnitude non-monotonic outlier at
  `num_orders=25`) before settling into ordinary, still-slow monotonic
  convergence — the exact symptom the Li/NVM/FFF literature (below) exists
  to improve. The investigation into *how* to fix it, done this session
  rather than deferred on faith:
  - **FFF** (Popov, E., & Nevière, M. (2001), "Maxwell equations in Fourier
    space: fast-converging formulation for diffraction by arbitrary
    shaped, periodic, anisotropic media," *J. Opt. Soc. Am. A* 18(11),
    2886-2894 — bibliographic details confirmed via `WebSearch` this
    session, not from memory) derives differential equations for the
    Fourier components directly, using a normal-vector field to get fast
    convergence for arbitrary-shaped anisotropic periodic media.
  - **NVM**, the specific 2D-grating formulation (Lalanne, P. (1997),
    "Improved formulation of the coupled-wave method for two-dimensional
    gratings," *J. Opt. Soc. Am. A* 14(7), 1592-1598 — also confirmed via
    `WebSearch`), is the earlier, narrower 2D-crossed-grating case of the
    same family of ideas; a related, frequently-cited follow-up is Schuster
    et al., "Normal vector method for convergence improvement using the
    RCWA for crossed gratings," *J. Opt. Soc. Am. A* 24(9), 2880 (2007).
  - Both JOSA A papers are paywalled in this environment (same situation
    Category 1 target 1.5's bounded literature search already documented
    for a different topic — see `references.md`) — the exact equations
    were **not** read or transcribed here, so nothing claiming to
    implement either paper's formula would meet `rules.md` AI Coding Rule
    1's bar. What *was* read in full is `../REFERENCE/S4`'s own
    implementation of this same family of techniques:
    `S4/S4/S4.h:49-71` (`use_polarization_basis`,
    `use_jones_vector_basis`, `use_normal_vector_basis`,
    `use_normal_vector_field` options) dispatching, per `S4.cpp:1905-1930`,
    to three separate implementation files — `fmm/fmm_PolBasisNV.cpp`
    (266 lines, the NVM path), `fmm/fmm_PolBasisJones.cpp` (378 lines),
    `fmm/fmm_PolBasisVL.cpp` (274 lines) — all three built **on top of**
    `fmm/fmm_FFT.cpp` (239 lines), i.e. S4's own closed-form analytic path
    (`fmm_closed.cpp`, already transcribed into this project) is
    *bypassed* entirely for any of these three options in favor of a
    discretized/rasterized (`use_discretized_epsilon`) representation of
    the permittivity and an FFT-generated normal-vector field over a
    `resolution`-parametrized grid (`S4.h`'s own doc comment for
    `resolution`, default 64). This is a **different Fourier-factorization
    architecture**, not an incremental correction to the existing analytic
    path — confirmed by reading the actual dispatch code
    (`S4.cpp:1905-1930`), not inferred from option names alone.
- **Alternatives considered**: Implementing a from-scratch NVM/FFF
  derivation independently (per `rules.md` AI Coding Rule 1's "derive
  independently and validate" option, used successfully for `staircase.py`,
  Phase 5) — rejected for now: unlike the staircase discretization (a
  well-precedented, low-risk geometric technique), NVM/FFF's correctness
  hinges on exactly the kind of subtle sign/normalization/field-decomposition
  convention this project's `rules.md` treats as too risky to re-derive
  without a transcribable source (the same reasoning already applied to
  the S-matrix star product and the Toeplitz subtraction rule, per
  `references.md`'s "Choosing a Reference for a New Phase" guidance) — and
  the one available from-source route (transcribing S4's `PolBasisNV`/
  `PolBasisJones`/`PolBasisVL` C++) is a genuinely large undertaking
  (~1150 lines across the four files above, needing a new FFT-based
  vector-field-generation subsystem this project doesn't have at all
  today), not a single atomic target's worth of work.
- **Trade-offs**: Deferring leaves the measured high-contrast 2D
  convergence weakness (target 3.3's fixture) unaddressed — a real,
  now-documented limitation, not swept under the rug (see
  `solve_layer_eigenmodes_patterned`'s docstring and
  `tests/test_fourier_convergence.py`). Implementing either technique
  properly would directly conflict with **ADR-002** ("Analytic shape
  Fourier transforms over raster+FFT of a pixelized mask"), since every
  one of S4's three polarization-basis paths requires the discretized/FFT
  representation ADR-002 explicitly rejected for a different reason
  (pixelization error at smooth boundaries) — revisiting that decision
  specifically for NVM/FFF, rather than for arbitrary-shape support, would
  need its own explicit re-evaluation, not a silent reversal buried inside
  a Category 3 target.
- **Impact**: `simulation.py`'s 2D solvers are unchanged; the accuracy
  limitation they already document (ordinary Laurent's rule, no Li/NVM
  correction, ordinary-convergence-rate-not-improved-rate at sharp 2D
  discontinuities) remains real and now has a measured fixture
  demonstrating it (`tests/test_fourier_convergence.py`), not just a
  docstring claim. Target 3.6 ("selected improvement") therefore has
  nothing to implement — recorded as its own explicit outcome, matching
  the register's own allowance for "explicitly decide implement/defer,"
  not silently skipped. Revisit this ADR if: (a) either paywalled paper
  becomes readable in this environment, (b) a second, structurally
  different NVM/FFF source is found for independent benchmarking (per
  target 1.5's precedent, one source alone was judged insufficient), or
  (c) a user explicitly requests reopening ADR-002 for this specific
  purpose.

## ADR-013: Revisit ADR-005 for a bounded analytic `Polygon` primitive (still no GDS/raster import)

- **Decision**: `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 4 targets 4.4
  ("Polygon design") and 4.5 ("Polygon primitive") ask for exactly the
  capability **ADR-005** deferred ("Add polygon-vertex analytic Fourier
  transforms... explicitly deferred, not rejected outright"). Per ADR-005's
  own instruction ("if a real need for imported layouts arises later, this
  ADR should be revisited with a new one, not silently overridden"), this
  is that revisit: add a `Polygon` shape using a closed-form **analytic**
  Fourier transform (no raster+FFT, no GDS import, no arbitrary-mask
  support) — narrower than what ADR-005 rejected, not a reversal of it.
  GDS import, raster masks, and self-intersecting/non-simple polygons
  remain out of scope; `Polygon` requires a simple (non-self-intersecting),
  closed, CCW-ordered vertex list, same restriction `S4/S4/pattern/pattern.h`'s
  own module docstring (lines 21-35) already imposes on every shape type.
- **Reason**: Reading `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`'s
  `POLYGON` case (lines 974-1008, alongside the already-transcribed
  `CIRCLE`/`RECTANGLE`/`ELLIPSE` cases at 951-973 in the same function)
  while investigating target 4.4 found that S4 computes a polygon's Fourier
  transform via a **closed-form boundary/edge-sum formula** — for vertices
  `v_0..v_{n-1}` (CCW, in the shape's local center-relative frame) and
  `k != 0`:
  ```
  S(k) = -i/(2*pi*|k|^2) * sum_{edges (p,q)} [(u_x k_y - u_y k_x) * sinc(k.u) * exp(-2*pi*i*k.rc)]
  ```
  where `u = v_q - v_p` (edge vector) and `rc = (v_q+v_p)/2` (edge
  midpoint) — a standard polygon-boundary Fourier-integral identity (the
  polygon's indicator function's gradient is a sum of delta functions on
  its edges; a 2D Fourier transform of a bounded region's indicator can
  always be reduced to a 1D boundary integral by Stokes' theorem, of which
  this is the closed-form evaluation for straight edges), **not** a raster
  or FFT operation at all. This directly resolves ADR-005's own
  "alternatives considered" note ("Add polygon-vertex analytic Fourier
  transforms (still avoids raster+FFT error)") in the affirmative: the
  polygon case is architecturally identical in kind to `Circle`/`Rectangle`/
  `Ellipse` (a closed-form function of `k` and the shape's parameters), not
  the raster+FFT alternative ADR-002 rejected. **Contrast with ADR-012's
  FFF/NVM deferral**: that decision was to defer specifically because
  proper vectorial 2D Fourier factorization (the `Epsilon2`/`kp` accuracy
  correction at a discontinuous interface) genuinely requires S4's
  discretized/FFT `PolBasisNV`/`PolBasisJones`/`PolBasisVL` machinery; a
  *polygon's own shape-level Fourier transform* is an unrelated, much
  simpler geometric question (the same question `Circle.fourier_transform`
  already answers) that happens to also have a closed form. Implementing
  `Polygon` does not require or imply revisiting ADR-012.
- **Accuracy contract** (target 4.4's explicit deliverable, decided before
  writing 4.5's implementation): the analytic edge-sum formula is **exact**
  for any simple (non-self-intersecting) polygon, to floating-point
  precision — there is no discretization/truncation error to characterize
  (unlike a raster+FFT approach, whose accuracy depends on grid resolution
  and would need its own convergence study). The only precondition is
  "simple polygon, correctly CCW-ordered vertices," matching `pattern.h`'s
  own stated (if S4-side unchecked) requirement; `Polygon.__post_init__`
  is not required to detect self-intersection (S4 itself does not, per
  that file's own comment, "not checked due to its complexity, so the
  input shapes must be sanitized elsewhere") — this project inherits that
  same, disclosed limitation rather than silently claiming a check that
  doesn't exist.
- **Alternatives considered**: (a) Full GDS/layout-file import — still out
  of scope, unchanged from ADR-005, no citable safe-parsing/format decision
  made here (see target 4.6 instead, which is deliberately scoped to a
  minimal *sougata_solver-native* JSON format, not a CAD/GDS format). (b)
  Raster+FFT polygon rasterization (Meent/TORCWA-style) — rejected for the
  same reason ADR-002 already gives (pixelization error at edges, and here
  it would be strictly worse than the exact closed form now available).
- **Trade-offs**: `Polygon` cannot represent curved boundaries (that's
  `Ellipse`'s/`Circle`'s job) and requires the caller to supply a correct,
  simple, CCW vertex list — no runtime self-intersection check, per the
  accuracy contract above. Still bounded scope relative to full GDS import:
  no file format, no CAD interoperability, no arbitrary raster mask.
- **Impact**: `geometry.Polygon` (target 4.5) is implemented as an analytic
  shape alongside `Circle`/`Rectangle`/`Ellipse`, sharing the same
  `Shape` ABC and `Pattern`/Fourier-factorization machinery unmodified —
  no changes needed to `fourier_factorization.py`, `eigenmodes.py`, or
  `simulation.py`. `references.md` updated with the exact `pattern.c` line
  citation. ADR-005's own text is left in place (not deleted), since GDS
  import and raster masks remain genuinely out of scope — only the narrow
  polygon-vertex sub-case it explicitly flagged as revisitable is revisited
  here.

## ADR-014: Bottom (reverse-side) illumination — already supported, no new API needed

- **Decision**: Category 6 target 6.6 asks whether reverse illumination
  (exciting from the transmission side instead of the incidence side) is
  required, and to design it separately if so. Answer, confirmed by direct
  test rather than assumed: **it is already achievable through the
  existing public `Simulation` constructor, with zero new code** — build
  `Simulation(lattice, list(reversed(layers)), num_orders, incidence=<old
  transmission material>, transmission=<old incidence material>)` and
  solve as usual. `Layer.thickness`/`Layer.pattern` carry no inherent
  z-direction, so reversing the layer list and swapping which material
  plays the incidence/transmission role is a complete, correct description
  of the mirrored problem — no new `Simulation` parameter, no
  `reverse=True` flag, no separate design needed.
- **Reason**: Verified directly (`tests/test_bottom_incidence.py`), not
  just argued: for a lossless, reciprocal, asymmetric two-film stack
  (air/1.46-quarter-wave/2.35-quarter-wave/glass), the reversed-stack
  simulation's transmittance at normal incidence matches the original
  orientation's transmittance to `~1e-15` — the Stokes transmittance-
  reciprocity relation for a lossless reciprocal medium, an independent
  physical law this project did not have to invent, confirming the
  "just reverse the list" recipe is physically correct, not merely
  "runs without crashing." Reflectance is **not** claimed equal in
  general (and isn't, in general, when tested at oblique incidence with
  mismatched incidence/transmission indices — a genuinely different R is
  the physically expected result there, not a bug) — only the
  direction-independent quantity (T at normal incidence) is used as the
  validation check, per `rules.md`'s "never fabricate a benchmark" rule.
- **Alternatives considered**: Adding a `reverse_illumination: bool`
  parameter to `Simulation`, or a `Simulation.solve_from_transmission_side(...)`
  convenience method — rejected as unnecessary API surface (per `rules.md`'s
  "don't add features not needed" guidance): the existing constructor
  already expresses this exactly, and a wrapper would only rename
  `list(reversed(layers))` plus a keyword swap without adding capability.
  Revisit only if a concrete future use case shows the manual recipe is
  genuinely error-prone in practice (e.g. a script class of bugs from
  forgetting to reverse the layer list), not preemptively.
- **Trade-offs**: None substantive — this is a documentation-and-test
  finding, not a code change. The recipe is one extra thing a user must
  know (documented in `design.md`'s API section and `CONVENTIONS.md`)
  rather than a self-explanatory constructor flag.
- **Impact**: No `simulation.py` change. `design.md`'s "API Design"
  section and `CONVENTIONS.md` both gain a short note with the recipe;
  `tests/test_bottom_incidence.py` freezes the reciprocity check as a
  permanent regression guard (a future refactor that broke this symmetry
  would be a real bug, not a style change).

## ADR-015: Interior-layer amplitude recovery via `partial_smatrix_up_to`, independently derived (not S4's `SolveInterior`)

- **Decision**: Category 9 target 9.3 (Phase 7) needs forward/backward
  modal amplitudes at an arbitrary interior interface. `phases.md`'s own
  Phase 7 deliverable already specifies the mechanism: use
  `SMatrixStack.partial_smatrix_up_to` (already implemented). This ADR
  records the specific formula built on top of it
  (`smatrix.interior_amplitudes`) and, per `rules.md` AI Coding Rule 1,
  flags it as **independently derived, not transcribed** — S4 itself
  solves this differently (`S4.cpp::SolveInterior`, called from
  `Simulation_ComputeLayerSolution`), a block-tridiagonal direct solve
  across every layer at once, a materially different algorithm from the
  star-product partial-stack approach this project's own architecture
  already committed to (`phases.md`, `architecture.md`'s dimension-agnostic
  `SMatrixStack`).
- **Reason**: Given the known full-stack incident amplitude `a0` and
  reflected amplitude `b_reflected` (both already computed by
  `Simulation.solve`), and the partial S-matrix `s_partial` from the
  incidence half-space up to interface `i`
  (`[a_i; b0] = s_partial @ [a0; b_i]`, `CONVENTIONS.md`'s S-matrix
  direction convention applied to that substack) — `b0` in this equation
  is the *same physical quantity* as `b_reflected` (both describe the
  backward amplitude in the incidence half-space; the physical field there
  cannot depend on where the stack is conceptually split), giving two
  equations solvable for the two unknowns `(a_i, b_i)`:
  ```text
  b_i = inv(S11) @ (b_reflected - S10 @ a0)
  a_i = S00 @ a0 + S01 @ b_i
  ```
  This is standard Redheffer-star-product algebra (the same block
  structure `smatrix.star_product` already uses), not a re-derivation of
  anything S4-specific — chosen over porting `SolveInterior` because this
  project's `SMatrixStack` already only builds left-cascaded partial
  products (`_partial`, one per layer boundary), not S4's full per-layer
  data structures, so this is the formula that composes with what already
  exists, per `rules.md`'s preference for building on already-validated
  blocks over introducing a second, parallel amplitude-solving mechanism.
- **Validation** (the "extra test scrutiny" AI Coding Rule 1 requires for
  independently-derived formulas, per `rules.md`): three independent
  checks, not just "runs without crashing" — (a) recovering amplitudes at
  the *last* interior interface and propagating them to the final layer's
  bottom must reproduce `a_transmitted` from the full solve exactly (an
  internal-consistency check with zero free parameters); (b) tangential
  `E`/`H` field continuity at an interface with no surface current (target
  9.5); (c) integrated real-space Poynting flux at an interior plane
  matches the already-oracle-validated per-order `z_poynting_flux`-based
  `R`/`T` (target 9.6, the category's own exit criterion). See
  `tests/test_field_reconstruction.py` for all three.
- **Alternatives considered**: Porting `SolveInterior` verbatim — rejected
  as unnecessary complexity (a full block-tridiagonal solver) for a
  capability the existing partial-stack architecture already supports with
  a few lines of algebra; revisit only if the star-product approach is
  found to be numerically worse-conditioned than direct block solving for
  some case the validation above doesn't cover.
- **Trade-offs**: One `_solve` (LU-based, per `troubleshooting.md`'s "never
  form `inv(A)` directly" rule) per interior-amplitude query, same
  numerical-stability profile as every other S-matrix operation in this
  project — no new instability class introduced.
- **Impact**: `smatrix.interior_amplitudes` (new function),
  `fields.py`'s depth-propagation ansatz (`CONVENTIONS.md`, same section)
  builds directly on its output.

## ADR-016: Instance-scoped Toeplitz-matrix cache on `Simulation`, gated on a measured timing case (Category 7 targets 7.3/7.4)

- **Decision**: cache exactly one artifact — the Toeplitz matrices built by
  `fourier_factorization.toeplitz_matrix`/`toeplitz_matrix_component` —
  keyed by `(kind, id(pattern), wavelength, ...)` on a plain dict attribute
  created fresh in `Simulation.__init__` (`self._toeplitz_cache`). Full
  design in `design.md`'s "Layer/Toeplitz Caching Design" section.
- **Reason**: `rules.md`'s Performance Requirements forbid introducing
  caching before Phase 9 "unless a specific correctness-validated
  capability is measurably too slow to use." Measured directly (not
  assumed), and corrected once after an initial measurement mistake (see
  `design.md`'s "Layer/Toeplitz Caching Design" for the full account): a
  20-point angle sweep at fixed wavelength for a single patterned layer
  (`num_orders=49`) takes ~1.44s uncached vs. ~1.01s cached, a ~30%
  reduction — real because `toeplitz_matrix`/`toeplitz_matrix_component`
  depend only on `(pattern, wavelength)`, not on incidence angle, so the
  same Toeplitz matrix is legitimately reusable across an entire angle
  sweep (Category 8 target 8.3, planned, not yet implemented). A Category 8
  sweep calling `solve()` hundreds of times would multiply that ~30%
  per-call saving into real wall-clock minutes. This satisfies the rule's
  exception clause; caching was not added speculatively just because
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` lists targets 7.3/7.4.
- **Validation**: `tests/test_layer_cache.py` — (a) equivalence: a stack
  solved with the cache populated matches the same stack solved with the
  cache forcibly cleared before every layer (i.e. every Toeplitz matrix
  recomputed from scratch), to numerical precision, per `rules.md`'s
  "validate the optimized path against the unoptimized one before
  trusting it"; (b) a call-counting monkeypatch of
  `fourier_factorization.toeplitz_matrix` confirms repeated identical
  `Pattern` objects trigger exactly one real computation, not `N`.
- **Alternatives considered**: caching the full per-layer eigenmode solve
  instead of just the Toeplitz matrix — rejected as exceeding target 7.4's
  explicit scope ("cache one safe artifact (for example, a Toeplitz
  matrix)") and as unnecessary: the measured cost above is dominated by
  Toeplitz reconstruction, not the eigensolve, for the repeated-pattern
  case that motivates this ADR. A value-based (structural-hash) cache key
  instead of `id(pattern)` — rejected because `Pattern`/`Shape` are
  unfrozen, unhashable dataclasses; building a structural hash would be new
  machinery unjustified by the one measured use case, and `id()`-based
  caching already has the right semantics for how repeated layers are
  actually constructed in this codebase (target 7.2's tests).
- **Trade-offs**: the cache never expires within a `Simulation` instance's
  lifetime; if a caller mutates a `Pattern` object in place after solving
  once and re-solves the same `Simulation` instance, they get a stale
  cached Toeplitz matrix for that `Pattern`'s `id()`. Documented as a
  caller-facing contract (`design.md`), not solver-enforced, matching the
  same construction-time-immutability assumption Category 4's shape
  validation already makes.
- **Impact**: `Simulation._cached_toeplitz`/`_cached_toeplitz_component`
  (new private methods), `Simulation._toeplitz_cache` (new instance
  attribute), no change to any public API or existing call signature.

## ADR-017: Layer-wise absorption as a flux-divergence combination of already-validated Phase 7 pieces (Category 7 targets 7.5/7.6)

- **Decision**: `SimulationResult.layer_absorption()` computes per-interior-
  layer absorbed power as `net_flux(top) - net_flux(bottom)`, normalized
  to incident power, using only already-implemented, already-validated
  Phase 7 building blocks (`smatrix.interior_amplitudes`,
  `fields.propagate_amplitudes`, `fields.z_poynting_flux`) — no new
  physics formula, no volumetric `Im(eps)*|E|^2` integral. Full derivation
  in `design.md`'s "Layer-Wise Absorption Design" section.
- **Reason**: per `rules.md` AI Coding Rule 1's preference (already applied
  in ADR-015) for composing already-validated blocks over deriving a new
  formula whenever the existing architecture already supports it — a
  volumetric integral would need a new spatial-integration formula and
  would need separately reconciling against the already-found
  `z_poynting_flux` factor-of-2 convention (Category 9 target 9.6), adding
  formula risk this flux-divergence approach avoids entirely.
- **Validation**: the energy-balance identity itself, per target 7.5's
  "define... the validation method before exposing an API" wording — `R +
  T + sum(layer_absorption()) == 1` for `tests/test_stress_regression.py`'s
  already-vetted lossy fixture (`eps=-396+80j`, sign-checked against
  `CONVENTIONS.md`'s passivity convention), finally closing the gap that
  same test file's docstring explicitly flagged ("layer-wise absorption...
  is Category 7 targets 7.5/7.6, still open" — no longer true after this
  ADR); and `layer_absorption() ~= [0, 0, ...]` for a lossless stack, to
  numerical precision.
- **Alternatives considered**: a volumetric `omega*Im(eps)*|E|^2` integral
  over each layer — rejected per the Reason above; a per-layer method
  requiring the caller to keep the original `Simulation` object around
  (matching how `tests/test_field_reconstruction.py` builds its own
  `SMatrixStack` externally) — rejected in favor of extending
  `SimulationResult` with a new `thicknesses` field so
  `layer_absorption()` is self-sufficient, consistent with
  `diffraction_efficiencies()`/`order_classification()`'s existing
  `SimulationResult`-method pattern.
- **Trade-offs**: `layer_absorption()` internally rebuilds an
  `SMatrixStack` from `self.thicknesses`/`self.all_modes` (already stored
  on `SimulationResult` after this ADR) rather than reusing one built
  during `Simulation.solve()` — a small, one-off, already-cheap
  reconstruction (S-matrix cascade is not the identified bottleneck; see
  ADR-016), not called unless a caller actually wants per-layer absorption.
- **Impact**: `SimulationResult.thicknesses` (new field, populated by
  `Simulation.solve()`), `SimulationResult.layer_absorption()` (new
  method), no change to any existing field or method.
- **Known limitation, found while validating this ADR, documented rather
  than silently worked around**: `layer_absorption()` inherits
  `interior_amplitudes`/`propagate_amplitudes`'s numerical-stability
  envelope. For a thick, highly lossy, high-`num_orders` layer, the
  deepest evanescent modes' backward-propagated amplitude can numerically
  overflow (`max(Im(q))*thickness` reaching ~38 for the `eps=-396+80j`
  fixture at `thickness=0.3`, `num_orders=25`), producing a nonsensical
  `layer_absorption()` value even though `R`/`T` themselves stay correct
  (they never reconstruct an interior amplitude). Full account and a
  regression guard on the failure symptom: `troubleshooting.md`,
  `tests/test_layer_absorption.py::test_interior_amplitude_reconstruction_can_numerically_overflow_for_thick_lossy_layers`.
  Not in scope to fix here — the validation tests for this ADR instead use
  a parameter regime (`thickness=0.05`) confirmed to stay well clear of
  this envelope.

## ADR-018: Conservative "every later point must also match" convergence criterion (Category 8 target 8.7)

- **Decision**: `sweep.find_convergence_index(values, tolerance)` returns
  the smallest index `i` such that **every later value** in the sequence
  (not just the immediate next one) stays within `tolerance` of
  `values[i]`, and requires confirmation from at least one later point (an
  index with nothing after it is never itself eligible — see the "found
  while testing" note below). Returns `None` if no such index exists.
- **Reason**: this project's own already-recorded evidence
  (`tests/test_fourier_convergence.py`, Category 3 targets 3.2/3.3) shows
  a high-contrast pattern's reflectance-vs-`num_orders` curve can be
  sharply non-monotonic at low harmonic-order counts (`num_orders=25`
  gives `R~0.214`, an order-of-magnitude outlier against both its
  immediate low-order neighbors and the eventual ~0.0236 converged value).
  A criterion that only checks "is the next point close" would misreport
  convergence at exactly this kind of transient. Requiring *every*
  remaining point to also match is the conservative choice target 8.7
  itself asks for.
- **Validation**: per target 8.8's explicit gating ("implement only after
  8.7 succeeds on thin-film, trench, and pillar fixtures"), validated
  against three structurally different cases before `auto_select_num_orders`
  (target 8.8) was implemented (`tests/test_harmonic_convergence.py`): a
  thin-film case where `num_orders` has no physical effect at all
  (converges at index 0, trivially and exactly); a moderate-contrast 1D
  TE grating that converges genuinely, at a non-trivial, non-final index;
  and the Category 3 high-contrast 2D pillar fixture's own known
  `num_orders=25` wobble, confirmed directly not to fool the criterion
  into anchoring on it.
- **Found while testing, fixed before trusting the function (not silently
  left in)**: a first version allowed the *last* index in the data to
  count as converged, since "every later value" is vacuously true when
  there are zero later values to check against. A test built to exercise
  a monotonically-diverging (never-actually-converging) sequence caught
  this immediately — the function returned the last index instead of
  `None`. Fixed by requiring at least one later point to actually confirm
  against (looping only up to `n-1`, not `n`), so returning an index is
  always a genuine confirmation, and `None` is a real, reachable outcome
  again — not merely a value the type signature claimed was possible.
- **Alternatives considered**: a relative (percentage) tolerance instead
  of an absolute one — left to the caller (pass a value computed from
  their own reference magnitude) rather than baked into the function,
  since `R`/`T` are already normalized quantities in `[0, 1]` where an
  absolute tolerance is the natural unit, unlike some other physical
  quantities this project's solvers produce.
- **Impact**: `sweep.find_convergence_index` (new function),
  `sweep.auto_select_num_orders` (new function, built directly on it, only
  added after the validation above passed).

## ADR-019: Overlay (layer-to-layer misregistration) is already achievable, no new API needed (Category 11 target 11.7)

- **Decision**: overlay between two patterned layers is already fully
  representable by the existing `Pattern`/`Layer`/`Simulation` API — give
  the lower layer's shape a `center` offset by the desired overlay error
  vector `(dx, dy)` relative to the upper layer's shape, within the same
  shared `Lattice`. No new parameter, field, or method is needed. Same
  treatment as ADR-014 (bottom illumination): "define a periodic unit-cell
  model" (this target's own wording) turns out to already be exactly what
  every multi-patterned-layer stack already does — every layer in a
  `LayerStack` shares one `Lattice`/reciprocal-vector set
  (`Simulation.__init__`), so two independently-offset shapes in two
  different layers are, by construction, already a periodic (period =
  the shared lattice) two-layer overlay model.
- **Reason**: `Simulation`/`Layer`/`Pattern` were never restricted to a
  single shape position per stack, and nothing in the existing
  construction-time validation (Category 4) assumes shapes across
  different layers share a center — this was true well before Category 11
  existed, just never previously framed as "overlay."
- **Validation**: not merely asserted — checked directly with a two-layer
  via-over-landing-pad fixture (`tests/test_overlay.py`). Three cases: (a)
  zero overlay, (b) a genuine sub-period overlay shift (confirmed to
  change R/T, as physically expected — an overlay error is a real
  structural change), and (c) a shift by *exactly one full lattice
  period*, which must reproduce case (a)'s result exactly (a periodicity
  self-consistency check specific to the "periodic unit-cell model"
  claim) — confirmed to `~1e-15`. Energy conservation (`R+T=1`) holds in
  all three cases.
- **Alternatives considered**: adding an explicit `overlay_dx`/`overlay_dy`
  parameter to `Simulation`/`LayerStack` — rejected as unnecessary
  indirection for something a caller can already do directly and
  transparently via `Shape.center`, per `rules.md`'s "don't add features
  not needed" guidance.
- **Impact**: no code change to `src/sougata_solver/` — documentation
  (this ADR, `references.md`, `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`) plus
  `tests/test_overlay.py` as the permanent regression guard for the claim.

## ADR-020: LER/LWR (stochastic edge roughness) explicitly deferred, not approximated (Category 11 target 11.8)

- **Decision**: line-edge/line-width roughness is **not** implemented,
  even as a "deterministic periodic approximation" — evaluated and
  explicitly deferred, per this target's own escape-hatch wording.
- **Reason**: genuine LER/LWR is inherently *stochastic* and, in real
  structures, has no periodicity at the base lattice's period at all — it
  is defined as a random perturbation of an edge along the (nominally
  infinite, non-periodic) line direction. RCWA is fundamentally a
  periodic-Fourier method: representing roughness at all requires either
  (a) a genuinely random, non-periodic edge (outside what any Fourier-
  modal method can represent without further approximation) or (b) a
  supercell large enough to hold many independent roughness "periods,"
  averaged over multiple random realizations — a substantial new
  capability (supercell lattice construction, per-realization random
  shape perturbation, and a statistical-averaging outer loop), not a
  "small target" scoped addition, and a materially different validation
  problem (there is no single deterministic answer to check against, only
  a statistical distribution). A *deterministic periodic* edge modulation
  (e.g. a sinusoidal wiggle at some sub-multiple of the lattice period)
  was considered as a cheaper proxy, but rejected: it approximates a
  qualitatively different physical effect (periodic edge modulation, a
  real but distinct scatterometry phenomenon in its own right) rather than
  genuine stochastic roughness, and presenting it as an "LER/LWR
  approximation" would overstate what it actually models — exactly the
  kind of overclaiming `rules.md` AI Coding Rule 1 warns against.
- **Alternatives considered**: implementing the supercell + multi-
  realization-averaging approach anyway — rejected as out of scope for a
  category of "small targets," and because averaging over how many
  realizations is "enough" for a statistically meaningful result is itself
  an unresolved design question this session did not have grounds to
  answer without further requirements from the project owner.
- **Impact**: no code change. `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`'s 11.8
  entry records this as an explicit, reasoned deferral, not a silently
  skipped target — revisit if genuine supercell/stochastic-averaging
  support is explicitly requested as its own scoped effort.

## ADR-021: Sparse/iterative linear algebra evaluated and rejected, not merely deferred (Category 12 target 12.5)

- **Decision**: sparse or iterative linear-algebra methods (sparse
  eigensolvers, Krylov-subspace iterative solvers, etc.) are **not**
  adopted for this project's dense eigenvalue/linear-solve operations —
  evaluated using `profiling/baseline_profile.py`'s measured benchmark
  dimensions (target 12.1) and rejected on structural grounds, not merely
  deferred for lack of time (unlike ADR-020's LER/LWR deferral, this is a
  "the technique doesn't apply here" finding, not a "not attempted yet"
  one).
- **Reason**: two independent, measured findings rule it out:
  1. **The coupling matrices are fully dense, not sparse.** Measured
     directly (not assumed): the direct-rule Toeplitz permittivity matrix
     `epsilon_hat` for an ordinary 2D patterned layer (`num_orders=49`,
     circular pillar) has **100% nonzero entries** (every pairwise
     Fourier-order coupling is nonzero, the expected structural
     consequence of a shape's continuous Fourier transform sampled at
     every reciprocal-lattice difference `G_i - G_j` — there is no
     banded/block-sparse pattern to exploit). Sparse linear algebra's
     benefit is proportional to how sparse the matrix actually is; here
     it is exactly zero.
  2. **Every mode is physically needed, not just a few.** Sparse
     *eigensolvers* (e.g. ARPACK/Lanczos-family) are advantageous when
     only a handful of extremal eigenvalues/eigenvectors are wanted out of
     a huge matrix — RCWA needs the *entire* mode spectrum (every Fourier
     order's forward/backward amplitude contributes to the S-matrix
     cascade and the final R/T), so there is no "few eigenvalues out of
     many" structure to exploit either, independent of sparsity.
  3. 12.1's own measurements (`design.md`'s "Linear-Algebra Baseline &
     Factorization-Reuse Design") show the eigensolve dominates runtime
     at the sizes this project actually uses (`num_orders` up to a few
     hundred at most, per Phase 4b's stress sweep) — the matrices never
     reach a size regime (tens of thousands+) where even a hypothetically
     sparse iterative method would clearly outperform a dense direct
     LAPACK solve's better constant factors and numerical robustness.
- **Validation**: the density measurement above, plus the general (well-
  established in the RCWA/Fourier-modal-method literature, not unique to
  this project) fact that Fourier-factorized permittivity matrices are
  inherently dense — cited as domain knowledge, not fabricated as if
  benchmarked against a specific paper (per `rules.md` AI Coding Rule 1,
  distinguishing a well-known structural property from a specific
  numerical claim needing citation).
- **Alternatives considered**: none pursued further, given both structural
  disqualifying factors above; a future session should re-evaluate only
  if a materially different problem regime appears (e.g. a supercell
  large enough to genuinely produce sparse coupling, which would only
  arise from work like ADR-020's deferred LER/LWR supercell approach, not
  from anything currently planned).
- **Impact**: no code change. Recorded so a future session doesn't
  re-investigate the same question without new evidence — the density
  measurement above is the reason to trust "no" rather than re-deriving it
  from first principles again.

## ADR-022: Instance-scoped per-layer eigenmode cache, implementing Category 12 target 12.3's design (Category 13 target 13.3)

- **Decision**: `Simulation._eigenmode_cache` caches each layer's full
  `LayerEigenmodes` result, keyed by `(id(layer), omega, kx.tobytes(),
  ky.tobytes())` — deliberately excluding `layer.thickness` (never
  consumed by an eigenmode solve, only by the downstream
  `propagation_smatrix` step) and the excitation's polarization
  amplitudes (only consumed downstream, by `incident_mode_amplitude`).
  Implements exactly the design `design.md`'s "Linear-Algebra Baseline &
  Factorization-Reuse Design" (Category 12 target 12.3) flagged but left
  unimplemented, since 12.3 explicitly asked only for the design.
- **Reason**: Category 12's own measurements showed the eigensolve, not
  the Toeplitz-matrix construction, dominates runtime at moderate-to-large
  `num_orders` — `sweep.sweep_polarization` (target 8.4) and
  `sweep.sweep_thickness` (target 8.5) both hold `omega`/`kx`/`ky` fixed
  across their entire sweep (polarization only affects the incident
  amplitude vector; thickness only affects propagation), so every layer's
  eigenmode solve is genuinely invariant across both — a real, not
  hypothetical, redundant-computation opportunity distinct from the
  angle-sweep scenario ADR-016's Toeplitz cache already covers.
- **Validation**: measured directly (not assumed) — a 20-point
  polarization sweep at `num_orders=49` on a single-patterned-layer pillar
  fixture: 1.52 s forced-uncached vs. 0.46 s cached (~3.3x). Equivalence
  to forced-uncached recomputation confirmed to `1e-12`
  (`tests/test_eigenmode_cache.py`), plus direct cache-entry-count checks
  for both the polarization-sweep and thickness-sweep scenarios (exactly
  one entry per layer, not one per sweep point), and a negative control
  confirming the cache correctly does *not* claim reuse across a genuine
  wavelength sweep (where `omega`/`kx`/`ky` do change every point).
- **Alternatives considered**: keying on `(id(layer.material)` or
  `id(layer.pattern))` directly instead of `id(layer)` — rejected as an
  unnecessary key-shape branch per dispatch type; `id(layer)` is uniform
  across every eigenmode-solve branch (uniform isotropic/diagonal/
  in-plane, 1D, 2D isotropic/anisotropic) since a `Layer`'s `material`/
  `pattern` reference never changes after construction, matching
  `_cached_toeplitz`'s established object-identity-keying convention
  (ADR-016).
- **Trade-offs**: same caller-facing contract as ADR-016's Toeplitz
  cache — a `Layer`'s `material`/`pattern` object must not be mutated in
  place after being wired into a `Simulation` (already an existing
  assumption, not a new one this ADR introduces).
- **Impact**: `Simulation._eigenmode_cache` (new instance attribute),
  `Simulation._cached_layer_eigenmodes`/`_solve_one_layer_eigenmodes` (new
  private methods, the latter extracted unchanged from `solve()`'s
  previous inline per-layer dispatch logic — a pure refactor, confirmed
  bit-for-bit equivalent to the pre-refactor numeric results before
  trusting it).

## ADR-023: Narrowly-scoped vectorized wavelength sweep for uniform-isotropic-only stacks (Category 13 target 13.4)

- **Decision**: `vectorized.sweep_wavelength_vectorized` batches a
  wavelength sweep across NumPy's native stacked-matrix `@`/
  `np.linalg.solve` for the one case where it's safe and simple: every
  layer uniform and isotropic, `num_orders=1` (a thin-film/multilayer
  stack, which never diffracts) — raises `ValueError` immediately for any
  patterned or anisotropic layer, or `num_orders != 1`, rather than
  silently falling back to something slower or subtly wrong.
- **Reason**: per `rules.md`'s Performance Requirements, vectorization
  work belongs to Phase 9 "unless a specific correctness-validated
  capability is measurably too slow" — the thin-film wavelength sweep is
  exactly this project's most common real use case
  (`structures/thin_film/*.py`), and every batched function is a direct,
  formula-identical re-expression of an already-cited, already-validated
  scalar function (`eigenmodes.solve_layer_eigenmodes_uniform`,
  `eigenmodes.build_kp_matrix`, `smatrix.interface_smatrix`/
  `propagation_smatrix`/`star_product`,
  `excitation.incident_mode_amplitude`) — no new physics formula, only a
  leading batch axis. `SweepResult` (Category 8) already provided the
  right output shape to slot this into without any caller-facing API
  change.
- **Validation** (`rules.md`'s "add a regression test comparing both
  paths"): `tests/test_vectorized_sweep.py` confirms bit-for-bit-scale
  (`atol=1e-12`) agreement with `sweep.sweep_wavelength`'s scalar loop
  across five polarization states, oblique/azimuthal incidence, a
  multi-layer (3 interior layer) stack, and a lossy material — every
  combination the scalar path already supports for this layer-type
  restriction.
- **Found and fixed before trusting it, exactly the discipline this rule
  requires**: a first draft of the batched eigenmode-solve helper omitted
  the `omega^2 * I` term from `build_kp_matrix`'s actual formula
  (`kp = omega^2*I - kappa`, not just `-kappa`) — caught immediately by
  the very first equivalence test run (`LinAlgError: Singular matrix`,
  not a silent wrong answer), fixed by re-reading `build_kp_matrix`'s
  exact source line-by-line rather than reconstructing the formula from
  memory a second time.
- **Measured benefit**: a 401-point wavelength sweep on a 2-interior-layer
  thin-film stack: 1.89 s scalar-loop vs. 0.060 s vectorized (~31x).
- **Alternatives considered**: extending this to patterned/anisotropic
  layers — explicitly rejected as out of scope; those paths' dense
  general eigensolve (`np.linalg.eig`, not a closed-form scalar formula)
  doesn't batch the same simple way, and attempting it now would be
  exactly the "general vectorized backend" work `rules.md` reserves for
  Phase 9, not this bounded proof.
- **Impact**: new `src/sougata_solver/vectorized.py` module. No change to
  any existing public API; `sweep_wavelength_vectorized` is purely
  additive.

## ADR-024: Parallelism evaluated and not implemented -- threading helps modestly, multiprocessing measured counterproductive (Category 13 target 13.5)

- **Decision**: no parallel-sweep API is added. Target 13.5 asks only to
  "profile and document" a parallelism decision, not to implement one —
  this ADR is that documentation.
- **Measurement** (a 20-point wavelength sweep, `num_orders=49` 2D pillar,
  the general — not vectorizable, see ADR-023 — patterned-layer path, on
  a 14-core machine):

  | Approach | Time | vs. serial |
  |---|---|---|
  | Serial | 1.644 s | 1x |
  | `ThreadPoolExecutor`, 4/8/14 workers | 1.263 / 1.197 / 1.116 s | ~1.3-1.5x |
  | `ProcessPoolExecutor`, 4/8/14 workers | 5.475 / 7.239 / 10.811 s | **0.15-0.3x (slower)** |

  A second check at `num_orders=81` (10 points, ~3s/point, a much heavier
  per-task granularity meant to rule out "process overhead just needs a
  bigger task"): serial `30.7s` vs. 8-worker process pool `44.1s` — **still
  slower**, not just at the cheap end.
- **Reason threading helps (modestly, not linearly)**: NumPy/SciPy's
  LAPACK-backed dense linear algebra (`np.linalg.eig`, `scipy.linalg.lu_factor`/
  `lu_solve`) releases the Python GIL during the underlying compiled
  calls, so genuinely concurrent execution happens for a real fraction of
  each `solve()` call — but not all of it (`Simulation.solve()`'s own
  Python-level orchestration, dict/list bookkeeping, and object
  construction stay GIL-bound), capping the achievable speedup well below
  the core count.
- **Reason multiprocessing is measured counterproductive here, at both
  task sizes tested**: the most plausible explanation (not independently
  confirmed by a separate BLAS-thread-count measurement, so stated as a
  plausible mechanism, not a proven one) is oversubscription — NumPy's
  BLAS/LAPACK backend is itself already internally multithreaded per
  process on this machine; spawning `N` additional OS processes, each
  independently trying to use multiple BLAS threads, creates far more
  concurrent threads than the 14 physical cores can usefully run,
  producing net slowdown from context-switching/cache-thrashing rather
  than the expected embarrassingly-parallel speedup. Process-spawn and
  result-pickling overhead is a secondary, compounding cost at the
  cheaper (`num_orders=49`) granularity specifically.
- **Validation**: `serial == threaded == processed` results confirmed
  identical (not just similarly-timed) across all three approaches before
  drawing any performance conclusion — correctness was never in question,
  only whether either approach is *worth doing*.
- **Alternatives considered**: pinning per-process BLAS thread counts
  (`OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1` per worker) before
  re-measuring multiprocessing — a real, standard fix for exactly this
  oversubscription pattern, and might well flip the conclusion — but not
  pursued: it adds environment-configuration complexity this "profile and
  document" target doesn't call for building, and doing so without
  re-measuring would itself violate the "never fabricate benchmark
  numbers" discipline. Left as a concrete, documented option for a future
  session if parallel sweeps become a real, requested need.
- **Impact**: no code change. If a future session wants a parallel sweep
  API, start from `ThreadPoolExecutor` (measured safe and modestly
  beneficial here) over `ProcessPoolExecutor` (measured counterproductive
  as naively applied), and re-measure — don't assume either conclusion
  transfers to a different machine's core count/BLAS configuration
  without checking.

## ADR-025: Reciprocity test scope — uniform layers only, Snell's-law-matched angles, not naive same-theta (Category 14 targets 14.5/14.6)

- **Decision**: reciprocity tests cover only **uniform (unpatterned)**
  layer stacks, comparing forward-direction transmittance at incidence
  angle `theta1` (in medium 1) against reverse-direction transmittance at
  the Snell's-law-refracted angle `theta2` (`n1*sin(theta1) =
  n2*sin(theta2)`, in medium 2) — **not** the same nominal `theta` value
  reused for both directions. Both lossless and lossy (but still
  reciprocal, i.e. ordinary absorption, not a gain/nonreciprocal medium)
  cases are covered.
- **Reason / non-obvious finding, verified numerically before writing any
  assertion (`rules.md`'s "verify before asserting" discipline)**: a
  first attempt compared T at the *same* `theta` value for both the
  forward and reversed (materials-swapped) stack — this is the naive
  reading of "reciprocity" and is **wrong**: measured directly, the
  difference between forward and reversed T grows with angle and
  reaches total mismatch (`T_reversed -> 0` from total internal
  reflection) at `theta=45 deg` for an air/glass asymmetric stack. The
  physically correct statement requires matching the transverse
  wavevector `kx` (conserved across every interface, Snell's law) — when
  angles are chosen this way instead, T reciprocity holds to `~1e-15/1e-16`
  at every angle tested (`0`/`15`/`30`/`40 deg`), including for a lossy
  reciprocal medium. This distinction between "same theta" and
  "Snell-matched theta" was not previously documented anywhere in this
  project (`tests/test_bottom_incidence.py`'s existing reciprocity check,
  Category 6 target 6.6, only tested normal incidence, where the two
  notions coincide since `theta=0` in both media — the naive/correct
  distinction was invisible until oblique incidence was actually tried).
- **Second finding, also verified numerically, not assumed**: total
  transmittance reciprocity (even at Snell-matched angles) does **not**
  hold for a **patterned** (diffractive) layer — measured directly for a
  1D lamellar grating, with mismatches up to ~0.56 at normal incidence.
  This is physically expected on reflection: reciprocity for a
  diffraction grating relates *individual diffraction orders* between
  the forward and reversed configurations (a materially more complex
  statement — which orders correspond to which — not derived or tested
  here), not the simple sum-over-all-orders total T a uniform stack's
  single-mode two-port picture gives. Scoping this category's reciprocity
  tests to uniform layers only is therefore a deliberate, verified
  boundary, not an arbitrary simplification — attempting to assert total-T
  reciprocity for a patterned layer would have been a **wrong** test, not
  a merely-incomplete one, which is exactly why it was checked first
  rather than assumed to generalize.
- **Validation**: `tests/test_reciprocity.py` — Snell-matched T
  reciprocity for a lossless stack across four angles; the same for a
  lossy (but reciprocal) stack; an explicit negative control confirming
  the *naive* same-theta comparison genuinely fails at oblique incidence
  (pinning the first finding above as a permanent regression guard, not
  just a decisions.md narrative); and a normal-incidence-only sanity
  check that a patterned layer's total T is *not* asserted reciprocal
  (documenting the scope boundary, not silently omitting it).
- **Alternatives considered**: deriving and testing the full per-order
  diffraction-grating reciprocity relation for patterned layers —
  explicitly out of scope for this pass (a materially larger derivation
  task); left as a documented possible extension, not silently dropped.
- **Impact**: `tests/test_reciprocity.py` (new file). No `src/sougata_solver/`
  code change — this category is validation-only, confirming existing
  solver behavior against a physical invariant, not adding a new
  capability.

## ADR-026: HDF5 deferred — NumPy `.npz` export is sufficient for current data shapes (Category 15 target 15.8)

- **Decision**: no HDF5 dependency or implementation is added. Target
  15.7's NumPy export (`export.py`'s `export_sweep_npz`/`load_sweep_npz`,
  a plain `.npz` archive of `parameter_values`/`reflectance`/
  `transmittance`/JSON-encoded metadata) covers every result-series shape
  this project currently produces.
- **Reason**: HDF5's actual advantages over `.npz` — hierarchical
  grouping, partial/chunked I/O into arrays too large to hold in memory,
  and cross-language structured access to deeply nested datasets — only
  matter once a single result exceeds comfortable in-memory `.npz` size
  or needs internal hierarchy `.npz`'s flat namespace can't express. This
  project's actual sweep outputs (`SweepResult`: one 1D array per swept
  parameter, at most a few thousand points; `fields.py`'s 2D field grids,
  Category 9 target 9.8, already export via existing means) are small
  flat arrays, not the scale or structure that would justify HDF5 — the
  same "evaluate before deciding" discipline this project already applied
  to Category 3's FFF/normal-vector decisions (ADR-006/007) and Category
  12/13's sparse-solver and parallelism decisions (ADR-021/024): add
  complexity only once a measured need, not a hypothetical one, exists.
- **Revisit when**: a future capability produces genuinely large (multi-GB)
  or deeply hierarchical structured output — e.g. a dense 3D field-grid
  sweep, or a multi-simulation batch run whose results are naturally
  nested by structure/wavelength/angle — that `.npz`'s flat, whole-file-
  in-memory model can no longer serve comfortably.
- **Impact**: none (`h5py` is not added to `pyproject.toml`'s dependencies).
  `export.py` (target 15.7) is the complete deliverable for this
  category's data-export scope.

## ADR-027: Test taxonomy — an `oracle` marker for a precise, greppable criterion; no per-tier marker where tiers overlap (Category 17 target 17.1)

- **Decision**: a new `pytest` marker, `oracle`, is applied to every test
  file that directly imports from `tests/oracles/` (confirmed by
  grepping for that import across all 54 pre-existing test files, not
  guessed) — 8 files, 136 tests. No marker is added for the other tiers
  `testing.md`'s "Testing Strategy By Tier" section already defines
  (unit, physical-invariant, integration, regression, acceptance).
- **Reason**: "oracle" has a precise, mechanically-checkable criterion
  (does this file import a named external-oracle module) that can't
  silently drift out of sync with reality — the marker means exactly
  what it says, by construction. The other tiers routinely overlap
  within a single test function (a reduction-to-simpler-case check is
  simultaneously a unit test and a permanent regression guard;
  `tests/test_1d_grating.py`'s own docstring already states it spans
  four tiers in one file) — forcing one marker per test would
  misrepresent that overlap, not clarify it, and risks becoming stale as
  soon as a test's actual scope changes but its marker doesn't. This
  target's own wording explicitly allows "filenames/docstrings **or**
  pytest markers" — the existing filename+docstring+Validation-Inventory
  combination (`testing.md`, target 14.1) already serves the other tiers
  adequately, audited this session for consistency via spot-checking a
  representative file sample before writing the "Test Taxonomy" section.
- **Impact**: `pyproject.toml` (`oracle` marker registration), 8 test
  files (one `pytestmark = pytest.mark.oracle` line each, no test body
  changed), `testing.md`'s new "Test Taxonomy" section. `pytest -m
  oracle` now runs exactly the system-tier suite in isolation.

## ADR-028: Performance regression guard — relative same-run ratio, not an absolute wall-clock threshold (Category 17 target 17.6)

- **Decision**: `tests/test_performance_regression.py` (new,
  `@pytest.mark.slow`) asserts that the ratio
  `eigensolve_time(num_orders=81) / eigensolve_time(num_orders=9)` for
  the same 2D-pillar fixture `profiling/baseline_profile.py` already
  uses stays below `1000` (and above `1`, guarding against a broken
  timing harness silently passing), rather than asserting either time is
  below some fixed number of milliseconds.
- **Reason**: `rules.md`'s Performance Requirements explicitly rule out
  hard-coded absolute wall-clock assertions, since timing is
  machine-dependent — a CI runner, a laptop on battery power, and a
  dedicated workstation will all report different absolute numbers for
  identical code. A ratio measured within one run on one machine cancels
  that dependence: both the numerator and denominator shift together if
  the machine is faster or slower, so the ratio itself is a genuine,
  machine-independent signal of *algorithmic* scaling behavior. The
  bound (`1000`) is set with ~6x headroom above Category 12's measured
  ~160x baseline for this exact fixture (`design.md`'s "Linear-Algebra
  Baseline & Factorization-Reuse Design" section) — generous enough to
  absorb ordinary run-to-run noise on an occasionally-run (`slow`-marked)
  guard, while still catching a genuine algorithmic regression (e.g. an
  accidentally-reintroduced `O(n^4)`-or-worse step) that would blow the
  ratio out far past any plausible noise margin.
- **Alternatives considered**: storing a cross-run baseline (e.g. a
  committed reference timing file, compared with tolerance) — rejected
  as needing periodic manual re-baselining as hardware/dependencies
  change, adding maintenance burden this target doesn't ask for; a
  same-run ratio needs no stored baseline to go stale.
- **Impact**: `tests/test_performance_regression.py` (new, 1 `slow`
  test). No `src/sougata_solver/` code change — this category adds a
  regression guard around already-existing, already-profiled behavior.

## ADR-029: 3D structure preview via matplotlib voxel rendering, reusing `plot_unit_cell`'s rasterization (no new dependency, no new analytic geometry)

- **Decision**: new `plotting.plot_structure_3d(layer_stack, lattice, *,
  resolution=40, extrusion_length=None, ax=None)` renders a full `Layer`
  stack (including a `staircase.py`-generated tapered via/pillar/trench)
  as a 3D solid, one non-cubic `Axes3D.voxels` slab per layer stacked at
  its real cumulative z-offset. In-plane cross-sections reuse a new
  `_rasterize_pattern` helper factored directly out of `plot_unit_cell`'s
  existing rasterization loop (same "later shape wins" precedence rule),
  rather than writing new analytic boundary-polygon extraction for each
  `Shape` subclass. `mpl_toolkits.mplot3d` ships inside the already-
  optional `matplotlib` dependency (`pyproject.toml`'s `dev` extra) — no
  new dependency added. `postprocessing/plot_structure_3d_preview.py` is
  the runnable demo/entry point, per ADR-009/010's `structures/`-builds
  vs. `postprocessing/`-shows split.
- **Reason**: requested directly by the project owner — "a GUI like
  Lumerical where I can see what I build based on my code," scoped down
  (via `AskUserQuestion`) to a static 3D solid preview first, with
  live/interactive parameter editing as an explicit, separate follow-up
  not covered here. Reusing `plot_unit_cell`'s rasterization instead of
  building analytic per-shape 3D boundary extraction means: zero new
  methods on `Shape` subclasses, zero new physics/geometry-formula
  citation burden (this module still only visualizes geometry the solver
  already consumes unmodified, the same exemption `plot_unit_cell`/
  `plot_layer_stack` already rely on), and a tapered structure "just
  works" with no taper-specific code — stacking each staircase `Layer`'s
  raster at its own z-offset reproduces the taper automatically, with
  step count exactly `num_slices`.
- **`Lattice1D` gap, found and fixed here, not new to this ADR**:
  `Lattice1D.b = (0, 0)` means `plot_unit_cell`'s own bounding-box logic
  already silently collapses to zero height for a 1D lattice — a latent
  gap in existing, shipped code, not something this new function
  introduced. `plot_structure_3d` makes the y-extent an explicit
  `extrusion_length` parameter (default: one period) rather than
  inheriting the silent collapse; `plot_unit_cell` itself is left
  unchanged (out of this ADR's scope — a 2D unit-cell preview of a 1D
  lattice was never a documented use case, and fixing it isn't needed for
  this target).
- **Measured performance finding, not assumed**: `Axes3D.voxels` scales
  worse than linearly in `resolution**2 * len(layer_stack)` (its hidden-
  face computation checks every voxel's neighbors) — measured directly on
  the dev machine: `resolution=20`/8 staircase slices renders in ~3s,
  `resolution=40`/16 slices takes ~2 minutes. `plot_structure_3d`'s
  docstring documents this directly (`rules.md` Performance Requirements
  — measure, don't assume), and the demo script defaults to
  `resolution=20`/`num_slices=8` for a fast first run rather than the
  function's own more-detailed default.
- **Known, honestly-documented limitation, not hidden**: this is an
  **opaque solid** render — a via (an air-filled hole in a solid
  substrate) is only visible from whichever face directly shows the
  opening; the tapered shaft itself is occluded by the surrounding
  substrate from a side view, confirmed directly by rendering the
  `postprocessing/plot_structure_3d_preview.py` tapered-via example and
  visually inspecting the output. A cutaway/transparency view would fix
  this but is out of scope for this static-preview target — flagged here
  for a future follow-up, not silently left undiscovered.
- **Alternatives considered**: an analytic per-shape 3D boundary/extrusion
  mesh (smoother-looking solid edges than a voxel raster) — rejected for
  v1 as materially more implementation work (a new boundary-sampling
  method per `Shape` subclass) for a static preview target that didn't
  need CAD-quality edges; a real 3D library (PyVista/VTK) — rejected per
  `rules.md` AI Coding Rule 4 (no new dependency without it being needed)
  since `mpl_toolkits.mplot3d` already ships with the existing optional
  `matplotlib` dependency and is sufficient for this target's static-view
  scope.
- **Impact**: `src/sougata_solver/plotting.py` (`_rasterize_pattern`
  helper factored out of `plot_unit_cell` with no behavior change,
  `plot_structure_3d` added), `postprocessing/plot_structure_3d_preview.py`
  (new), `tests/test_plotting.py` (6 new tests, structural not pixel —
  return-shape, z-extent, material-legend count, staircase voxel count,
  `Lattice1D` non-collapse regression, empty-stack error path). No
  `src/sougata_solver/` physics module touched. `phases.md` Phase 10
  records this as new, explicitly-scoped capability rather than a silent
  addition, per `rules.md` AI Coding Rule 2.
- **Correction (same day, caught by the project owner reviewing an actual
  render, not by a test)**: two real bugs found in the first shipped
  version. (1) `postprocessing/plot_structure_3d_preview.py`'s
  `_build_tapered_trench`/`_build_trench_ocd_sweep` passed only the
  finite patterned `layers` list to `plot_structure_3d`, never wrapping
  it in the actual `layer.LayerStack(layers, incidence=..., transmission=...)`
  each real `structures/*.py` script builds — so no substrate/incidence
  half-space rendered at all, making a trench look like a floating slab.
  Fixed by building a real `LayerStack` per structure, with each
  builder's incidence/transmission **confirmed by reading the actual
  `structures/*.py` file, not assumed**: `tapered_trench.py` genuinely
  has `transmission=substrate` (Si), but `tapered_via.py` and
  `trench_ocd_sweep.py` both use `incidence=air, transmission=air` — no
  substrate at all, a free-standing/suspended structure as those scripts
  are actually written. This is now documented directly in each
  `_build_*` function's docstring rather than silently assuming a
  substrate that doesn't exist in the real physics. (2) Independently,
  `plot_structure_3d`'s `math.inf`-layer end-cap formula
  (`0.5 * max(finite_thicknesses)`) used `max()` over *individual staircase
  slice* thicknesses — for a many-slice taper each slice is only
  `total_depth / num_slices`, so a substrate rendered as a sliver
  *thinner* than the trench depth, backwards from physical reality.
  Fixed to `1.5 * sum(finite_thicknesses)` (the whole finite stack's
  total depth, generously scaled up) — `plot_layer_stack`'s own
  `0.5 * max(...)` convention is deliberately **not** reused here for
  exactly this reason (see `plot_structure_3d`'s updated docstring).
  New regression test,
  `tests/test_plotting.py::test_plot_structure_3d_semi_infinite_end_caps_exceed_staircased_patterned_depth`,
  pins that a `LayerStack`'s end-caps must render thicker than the total
  patterned depth for a many-slice staircase. Both fixes re-verified by
  re-rendering all three demo structures and visually inspecting the
  output (not just re-running tests) — `tapered_trench.py`'s preview now
  shows a contiguous solid Si block (ridge and substrate are the same
  material, so they correctly merge visually) with the trench etched as
  a notch from the top, sitting above a substrate cap thicker than the
  trench depth; `trench_ocd_sweep.py`'s preview now honestly shows a
  free-standing Si ridge in air on both sides, matching its actual
  `Simulation(..., incidence=air, transmission=air)` call. 27 tests in
  `tests/test_plotting.py` (1 new), 698 tests pass project-wide, `ruff
  check .` clean.
- **Second correction, same day, again caught by the project owner
  reviewing an actual render**: fixing the first correction above (a
  single shared `1.5 * sum(finite_thicknesses)` end-cap height) traded
  one bug for another — rendering the incidence (air) half-space at the
  *same* generous thickness as the transmission (substrate) half-space
  made air look as visually dominant as the substrate, backwards from a
  real device cross-section where the substrate is the bulk and the
  incidence medium above is just context. Fixed to **asymmetric**
  end-caps, keyed by **position in the list**, not a name lookup: the
  first layer (always the incidence half-space, per `LayerStack`'s own
  `[incidence, ...finite..., transmission]` ordering) renders at
  `0.3 * sum(finite_thicknesses)`, clearly thinner than the patterned
  stack; the last layer (transmission/substrate) renders at
  `2.5 * sum(finite_thicknesses)`, clearly the dominant block.
  `tests/test_plotting.py::test_plot_structure_3d_semi_infinite_end_caps_exceed_staircased_patterned_depth`
  updated to pin both halves of this asymmetry (substrate cap exceeds the
  trench depth; incidence cap is clearly thinner than the substrate cap),
  not just a single symmetric bound. Re-verified by re-rendering all
  three demo structures again: `tapered_trench.py` now shows a thin air
  cap above a dominant Si substrate block with the trench notch near the
  top — matching the project owner's stated physical expectation
  exactly. 26 tests in `tests/test_plotting.py` (net unchanged count,
  one test rewritten not added), 698 tests pass project-wide, `ruff
  check .` clean.
- **Third correction, same day, again caught by the project owner
  reviewing an actual render**: two rounds of tuning a fabricated
  end-cap thickness (first too thin, then too visually dominant) was the
  wrong axis to iterate on — there is no size for a genuinely
  semi-infinite half-space that reads correctly next to a finite
  patterned stack. The project owner explicitly pointed back to the very
  first (pre-`LayerStack`) render of `tapered_via.py` — the finite
  patterned layer alone, no end-caps at all — as the correct reference
  style, and asked that the etched (air) region simply appear as a
  feature embedded in the substrate, extruded along z. **Reverted**
  end-cap rendering entirely: `math.inf`-thickness layers are now
  filtered out of `plot_structure_3d`'s input before rendering and never
  drawn, full stop, rather than assigned any fabricated height. This is
  simpler than either prior attempt and matches the reference exactly,
  because a via/trench's own finite patterned layer(s) already read
  correctly on their own — the background material fills the whole
  cross-section except the via/trench footprint, which already *is*
  "air etched into a solid substrate" with no separate block needed.
  `layer_stack` inputs with no finite layer now raise `ValueError`
  (`"...finite-thickness layer..."`), a new, more specific error than the
  original bare-empty-list check. Tests updated:
  `test_plot_structure_3d_semi_infinite_end_caps_exceed_staircased_patterned_depth`
  replaced with `test_plot_structure_3d_semi_infinite_layers_are_not_rendered`
  (pins that a `LayerStack`'s inf layers contribute nothing to the
  rendered z-extent) and a new `test_plot_structure_3d_all_infinite_layer_stack_raises`.
  Re-verified by re-rendering all three demo structures a final time and
  visually comparing directly against the project owner's own referenced
  screenshot — `tapered_via.py`'s new render is now visually identical to
  it, and `tapered_trench.py`/`trench_ocd_sweep.py` now show the same
  clean single-block style (majority Si, minority air at the trench,
  extruded fully along z). 27 tests in `tests/test_plotting.py`, 699
  tests pass project-wide, `ruff check .` clean. This function's
  docstring has been rewritten (not just amended) to describe the final,
  reverted behavior directly rather than layering a fourth explanation on
  top of three prior ones — the two abandoned end-cap-sizing attempts
  are preserved in this ADR's own prior two "Correction" entries above,
  not in the docstring.
- **Fourth correction, same day, again caught by the project owner
  comparing two renders side by side**: `pillar_array.py` (Si pillar in
  air) and `via_array.py` (air via in Si) — physical opposites, one
  mostly air with a solid feature, the other mostly solid with an air
  feature — rendered with visually identical-looking colors. Root cause:
  `color_map` assigned colors by **encounter order**
  (`Pattern.background` first, then `Pattern.shapes[*].material`), so
  "air" landed on whichever color slot its background/shape position
  happened to occupy in each structure — coincidentally the same slot
  both times, purely because of *where* air appeared in each pattern, not
  *what* it was. Fixed to sort `material_names` case-insensitively before
  assigning colors, so the same material name always gets the same color
  across every call — "air" now renders identically whether it's the
  background or the shape, making majority-vs-minority material
  immediately visible from color alone.
  `tests/test_plotting.py::test_plot_structure_3d_color_keyed_by_material_name_not_encounter_order`
  pins this directly: builds a pillar pattern (`background=AIR,
  shape=SI`) and a via pattern (`background=SI, shape=AIR`) and asserts
  both legends assign "air" (and "si") the identical color. Re-verified
  visually with a genuine side-by-side comparison figure (two
  `plot_structure_3d` calls composed into one figure via its existing
  `ax=` parameter) — air is now consistently blue and Si consistently
  orange in both panels, making the two structures' physical difference
  legible for the first time despite their near-identical silhouette at
  this radius/period. 28 tests in `tests/test_plotting.py`, 700 tests
  pass project-wide, `ruff check .` clean.
- **Fifth correction, same day, caught by the project owner questioning
  whether the preview's numbers actually matched the real code**: found
  correct — `postprocessing/plot_structure_3d_preview.py`'s `_build_*`
  functions hand-copied each real `structures/*.py` script's constants
  (per this ADR's own "Impact" note above) rather than reading them, and
  one had already gone stale: `_build_tapered_via` used
  `tcd=0.36e-6, bcd=0.26e-6, spacing=0.34e-6`, but
  `structures/via/tapered_via.py`'s real, current values are
  `TCD=0.48e-6, BCD=0.20e-6, SPACING=0.22e-6` — confirmed by re-reading
  that file directly, not assumed. The same builder had also silently
  swapped the shape/background material roles (`shape_material=air,
  background_material=substrate`, an etched hole) relative to the real
  script's actual call (`shape_material=via` (Si), `background_material=air`
  — a raised solid post, the same construction `pillar_array.py` uses,
  not a hole at all). **Fixed structurally, not just numerically**: every
  `_build_*` function now dynamically loads the real `structures/*.py`
  file as a module via `importlib.util.spec_from_file_location` (new
  `_load_structure_module` helper) and reads its actual current
  constants (`m.TCD`, `m.PERIOD`, ...) directly, rather than retyping a
  duplicate — this specific bug class (a hand-copied number silently
  drifting from the file it was copied from) is now impossible to
  reintroduce, since the numbers are re-read from the real file on every
  run. Loading a `structures/*.py` file this way only executes its
  module-level constant definitions (confirmed safe: no `Simulation.solve()`
  call, since every real script's `main()` is guarded by
  `if __name__ == "__main__":`, false on import); `tapered_trench.py`'s
  builder additionally reuses that script's own `_material()` helper
  function (imported from the loaded module, not reimplemented) so a
  future change to its material-name-resolution logic can't silently
  diverge either. No changes to `src/sougata_solver/plotting.py` itself —
  this was entirely a `postprocessing/` script bug; all 700 project-wide
  tests re-run and confirmed still green, `ruff check .` clean. All five
  structures re-rendered and re-verified visually after the fix —
  `tapered_via.py`'s preview now correctly shows a raised Si post (the
  same shape as `pillar_array.py`), not an etched hole.

## ADR-030: `build_geometry(**overrides)` convention across every `structures/*.py` script

- **Decision**: every one of the 18 `structures/thin_film/*.py`,
  `structures/trench/*.py`, `structures/via/*.py` scripts now exposes a
  top-level `build_geometry(**overrides) -> (layers, lattice, incidence,
  transmission)` function. `main()` calls it for the geometry-building
  portion it already had — a **pure extraction**, no behavior change; the
  same lines that already built materials/`Pattern`/`Layer`(s)/`Lattice`
  just moved into a named function. Each script's geometrically-tunable
  EDIT constants (period, radii/CDs, thickness, fill_factor, num_slices —
  *not* wavelength/angle/polarization/output-path, which aren't part of
  "the structure") became optional keyword parameters, each defaulting to
  `None` → falls back to that constant's real current value when not
  overridden.
- **Reason**: requested directly by the project owner — a generic loader
  that works on *any* `structures/*.py` file, not one hand-written
  `_build_*` wrapper per structure kind (ADR-029's fifth correction had
  already shown that duplication path breeds exactly the kind of silent
  drift bug that got caught). One function per script, reused by both
  consumers — the static preview (`postprocessing/plot_structure_3d_preview.py`,
  calls it with no overrides) and the live GUI (`postprocessing/live_structure_viewer.py`,
  calls it repeatedly with slider-driven overrides) — means there is now
  exactly one place per structure where its geometry is built, not two
  (or three, counting a hypothetical future consumer) that could drift
  apart from each other or from the real script.
- **Per-file treatment, not uniform**: 13 of the 18 files had an
  already-clean, side-effect-free geometry-building block (module level
  or the top of `main()`) — trivial extraction, sometimes keeping a
  pre-existing `ValueError` guard clause in `build_geometry()` too (e.g.
  `sio2_on_si_thin_film.py`'s substrate/film-name validation).
  `thin_film/anti_reflection_coating.py` had no `main()` at all (fully
  module-level script) — gained both a `build_geometry()` and a `main()`,
  same pattern as every other file. 3 files
  (`trench/tapered_trench.py`, `via/tapered_pillar.py`, `via/tapered_via.py`)
  rebuild their staircase geometry **inside** a `for num_slices in
  SLICE_COUNTS:` convergence-sweep loop — `build_geometry(num_slices=None, ...)`
  defaults `num_slices` to `SLICE_COUNTS[-1]` (the finest/most
  representative), and `main()`'s existing sweep loop now calls
  `build_geometry(num_slices=n, ...)` per iteration instead of inlining
  the staircase call, with identical numeric output (re-verified directly
  by re-running every sweep script, not assumed — `tapered_via.py`'s
  6-row convergence table matches the pre-refactor values exactly). 2
  files (`trench/trench_ocd_sweep.py`, `via/tsv_ocd_sweep.py`) rebuild
  geometry inside a `for bottom_cd in BOTTOM_CDS:` loop feeding
  `sweep.sweep_wavelength` — `build_geometry(bottom_cd=None, ...)`
  defaults to `BOTTOM_CDS[0]`; each `main()` still recomputes its own
  `OCDTrapezoidParams` per iteration purely to read `sidewall_angle_deg`
  for printing/metadata (a cheap derived property, not physics, so this
  small duplication is an accepted tradeoff against `build_geometry()`'s
  fixed 4-tuple return contract), and both re-run end-to-end with
  unchanged R/T/sidewall-angle output.
- **`postprocessing/plot_structure_3d_preview.py` rewrite**: every
  hand-written `_build_*` function and the `_BUILDERS` dict were deleted;
  `STRUCTURE` is now a relative path under `structures/` (e.g.
  `"via/tapered_via.py"`), loaded via the same `_load_structure_module`
  `importlib` helper ADR-029's fifth correction introduced, then simply
  `module.build_geometry(**OVERRIDES)`. Works on any conforming file,
  including ones added after this session — the generic-loading capability
  the project owner explicitly asked for.
- **Impact**: all 18 `structures/*.py` files (pattern described above, not
  enumerated per-file/per-line), `postprocessing/plot_structure_3d_preview.py`
  rewritten. `ruff check .` clean; every script directly re-run
  end-to-end and confirmed unchanged R/T (or R/T/A, or convergence-table)
  output; full 700-test project-wide suite re-run and confirmed green
  (this refactor touches no `src/sougata_solver/` physics module).

## ADR-031: Live PyVista structure viewer — new `gui` optional dependency

- **Decision**: `postprocessing/live_structure_viewer.py` (new) opens an
  interactive PyVista `Plotter` window with one slider per tunable
  `build_geometry()` parameter (auto-discovered via `inspect.signature`,
  default/range read from the real script's own constant — same
  never-hand-copy principle as ADR-030); moving a slider rebuilds the
  mesh and re-renders live. `pyvista` (pulls in `vtk`) is added as a new
  `gui` extra in `pyproject.toml`'s `[project.optional-dependencies]`,
  separate from the existing `dev` extra (not needed for the fast test
  suite/CI, ~85MB of binary wheels).
- **Reason**: requested directly by the project owner — "better option is
  to make live gui where if i change dimension, live structure will be
  changed" — with PyVista explicitly chosen (via `AskUserQuestion`) over
  matplotlib sliders specifically because matplotlib's `voxels()` is too
  slow to redraw interactively (`plot_structure_3d`'s own measured
  ~3s-2min range). Satisfies `rules.md` AI Coding Rule 4 ("never add a
  dependency... without it being explicitly requested and recorded here").
  **Verified installable and functional before building on it**, not
  assumed: `pip install pyvista` succeeded cleanly in this project's
  `.venv` (`pyvista==0.48.4`, `vtk==9.6.2`), and a basic offscreen
  `Plotter`/`RectilinearGrid`/screenshot round-trip was confirmed working
  before any GUI code was written, per this ADR's own "first
  implementation step" gate.
- **Design**: `live_structure_viewer.build_meshes(layer_stack, lattice,
  resolution)` is a pure data-prep function (no display side effects,
  same split `plotting.plot_structure_3d` already uses) — reuses
  `plotting.rasterize_pattern` (promoted from `_rasterize_pattern`, a
  one-line rename, this session) for in-plane cross-sections, skips
  `math.inf`-thickness layers entirely (matching `plot_structure_3d`'s
  own third correction), and sorts material names case-insensitively
  before assigning colors (matching that function's fourth correction) —
  so a given structure's live-GUI colors match its static-preview colors
  exactly, one convention, not two independently-drifting ones. Slider
  callbacks call `module.build_geometry(**current_slider_values)` fresh
  on every change (`num_slices` rounded to `int`, everything else left
  float) and do a full actor swap (`remove_actor` + `add_mesh`), not an
  incremental mesh patch, since changing any parameter can change the
  layer count (a taper's slice count) or the raster shape.
- **Verified, with an explicit limitation stated rather than overclaimed**:
  this sandbox cannot display a live interactive window, so verification
  here is: (1) `build_meshes()` unit-level correctness confirmed directly
  (grid count matches `num_slices`, z-extent matches total thickness,
  material names/bounds correct) against `tapered_via.py`; (2) the full
  slider-rebuild pipeline exercised programmatically in PyVista's
  `off_screen=True` mode — simulating three slider changes
  (`num_slices=32→4`, then `tcd` widened) and confirming the actor count
  and rendered screenshot update correctly each time, visually inspected
  and matching `plot_structure_3d`'s established color convention. The
  actual interactive-dragging experience (redraw latency, slider
  usability) has **not** been visually confirmed by a human in this
  session — the project owner will need to run
  `python postprocessing/live_structure_viewer.py` themselves to confirm
  that, and this is stated plainly rather than claimed as verified.
- **Impact**: `postprocessing/live_structure_viewer.py` (new),
  `pyproject.toml` (`gui` extra added), no changes to
  `src/sougata_solver/` beyond the `rasterize_pattern` rename (no
  behavior change, `tests/test_plotting.py` re-run and confirmed green).

## ADR-032: `PyRCWA` vendored as a second, structurally-different Phase 3 (1D grating) oracle

- **Decision**: `REFERENCE/PyRCWA` (github.com/vitamingcheng/PyRCWA, MIT)
  cloned as a new vendored sibling repo, added at the project owner's
  explicit request. `tests/oracles/rcwa_1d_pyrcwa.py` hand-transcribes its
  normal-incidence TE path -- a general 2D P/Q eigenoperator restricted to
  1D via harmonic truncation, structurally different from
  `rcwa_1d_gaylord.py`'s reduced TE-specific `A = KX2 - E` operator -- and
  `tests/test_1d_grating.py::test_te_matches_pyrcwa_oracle_at_normal_incidence`
  cross-checks it against `sougata_solver`'s own Phase 3 fixture.
- **Reason**: the project owner named this specific repo and, per
  `AskUserQuestion`, wanted it used as a new oracle/cross-check for an
  existing phase, not a new-phase reference. Assessed honestly before
  committing effort to it: small (v0.0.1, 6 stars, 12 commits, **no test
  suite**) by the same code-quality bar `references.md` already documents
  for why other candidates were passed over (e.g. `EMpy/EMpy/RCWA.py`'s
  "author-acknowledged hack") -- but the actual RCWA math, read directly,
  is a real, coherent implementation in the same architecture family as
  this project's own solver, so it was used anyway, on the code's own
  merits rather than its popularity signal.
- **Practical constraint, not a design choice**: `PyRCWA` cannot even be
  imported in this project's environment without a numpy-global
  monkeypatch (`np.NAN`/`np.mat`, both removed in NumPy 2.0, used
  throughout `pyrcwa/core.py`/`solver.py`) -- this reinforces, not just
  follows, the project's existing "hand-transcribe, never import a
  `REFERENCE/` repo into test infrastructure" convention (`rules.md` Rule
  7): baking a numpy-global monkeypatch into permanent test setup would
  itself be a real fragility risk.
- **Scoped to normal incidence only**: `PyRCWA`'s `(alpha, theta)`
  oblique-angle convention was not confidently mapped onto this project's
  `(theta, phi)` convention with the time available for this addition;
  rather than assert an unverified equivalence, the oracle (and its test)
  are restricted to the one regime -- normal incidence -- where both
  conventions provably coincide (confirmed by hand from `PyRCWA`'s own
  `compute_diffraction_efficiency`, `solver.py:154-194`).
- **Verified before trusting, twice, both real findings**: (1) a live run
  of the actual `PyRCWA` code (via a one-off local numpy-compatibility
  monkeypatch, never part of the permanent oracle) at increasing
  `fft_resolution` (2001/4001/8001) gave `Total R = 0.912281 / 0.912194 /
  0.912151`, monotonically converging toward `sougata_solver`'s own
  `0.912109` for the identical fixture -- confirming the formulations
  agree and the gap is `PyRCWA`'s own FFT/raster discretization error, not
  a real mismatch, *before* any transcription was trusted. (2) The first
  draft of the transcription itself had a real bug -- used a per-
  diffraction-order `Kz` matrix where `solver.py:168` actually reads a
  single scalar (the incident wave's own `kz`; `self.Kz`, though computed,
  is never read in `compute_diffraction_efficiency`) -- caught immediately
  as a loud `nan`, not a silently wrong number, fixed by re-reading the
  source line-by-line. After the fix, the transcription (using this
  project's already-validated analytic Fourier-coefficient convolution
  matrix, `rcwa_1d_gaylord._toeplitz_convolution_matrix`, reused not
  reimplemented, in place of `PyRCWA`'s own FFT-of-raster route -- the
  same substitution `rcwa_1d_gaylord.py` already made) matches
  `sougata_solver` to ~1e-10, not just the live run's ~2e-4.
- **Impact**: `REFERENCE/PyRCWA` (new vendored repo, read-only, never
  modified), `tests/oracles/rcwa_1d_pyrcwa.py` (new), one new test in
  `tests/test_1d_grating.py` (701 tests pass project-wide, up from 700),
  `references.md` updated with the new table row. No
  `src/sougata_solver/` change -- this is oracle infrastructure only.

## ADR-033: Linear-polarization `alpha` convention flipped to match a commercial RCWA tool (0=P, 90=S)

- **Decision**: `CONVENTIONS.md`'s "Worked polarization examples" table
  (Category 6 target 6.1) and its two structure-script implementations
  (`structures/thin_film/custom_multistack.py`'s `_jones_state`,
  `structures/thin_film/sio2_on_si_thin_film.py`'s
  `_polarization_amplitudes`) now use `s_amplitude=sin(alpha)`,
  `p_amplitude=cos(alpha)*exp(1j*delta)` -- the opposite trig assignment
  from before (`s=cos(alpha)`, `p=sin(alpha)*exp(1j*delta)`). Net effect:
  `alpha=0` is now pure P (was pure S/TE); `alpha=90` is now pure S (was
  pure P/TM). `custom_multistack.py`'s `POLARIZATION_STATES_DEG` entries for
  `"TE"`/`"TM"` had their `alpha_deg` swapped (`90.0`/`0.0`) to keep
  producing the correct physical amplitudes under the new formula --
  `"linear_15deg"`/`"linear_30deg"`/`RCP`/`LCP`/`elliptical_*` entries keep
  their existing `(alpha, delta)` values unchanged (verified below).
- **Reason**: the project owner is validating this solver's thin-film
  output against a commercial RCWA tool (Lumerical FDTD) and supplied its
  actual polarization-mixing script:
  `R_linear = sin(alpha)^2 * Rs_power + cos(alpha)^2 * Rp_power`, with an
  explicit in-script comment `(0=P, 90=S)`. That is the opposite reference
  axis from this project's pre-existing convention. Comparing the solver's
  `linear_15deg`/`linear_30deg` (45 deg incidence, `sio2_sio_ni_sio2_on_
  semi_infinite_si` stack) against that tool's `Linear15_Linear30.txt`
  export directly showed both a large apparent magnitude gap (peak `R`
  0.53 vs 0.35) and a reversed 15-vs-30 ordering. Back-solving the tool's
  raw `Rss`/`Rpp` from its two exported curves (a per-wavelength 2x2 linear
  solve against its own stated mixing formula) and comparing directly to
  this solver's own pure-TE/pure-TM `R` at the same angle matched to
  ~0.1% absolute (max diff 0.0013, RMS 0.0005 -- noise-level, from the
  tool's coarser non-uniform wavelength grid) -- proving the entire
  discrepancy was this labeling mismatch, not a solver or oracle physics
  error. Flipping the convention makes a solver `linear_Xdeg` state equal
  the commercial tool's `Linear X deg` state directly, with no `90-alpha`
  conversion needed for future comparisons.
- **Why RCP/LCP/elliptical entries needed no numeric change**: the phase
  term stays attached to `p_amplitude` in both the old and new formula (only
  which trig function multiplies `s`/`p` swapped) -- for the specific
  `alpha=45` used by `RCP`/`LCP`, `sin(45)=cos(45)`, so `(s_amplitude,
  p_amplitude)` is numerically identical under both formulas
  (`RCP`: `(1/sqrt2, i/sqrt2)`; `LCP`: `(1/sqrt2, -i/sqrt2)`). General
  elliptical/linear entries at other angles do change in physical meaning,
  which is the intended effect of this ADR.
- **Scope checked, not assumed**: grepped every `structures/thin_film/*.py`,
  `structures/trench/*.py`, `structures/via/*.py`, and `postprocessing/*.py`
  file for this alpha-to-amplitude formula -- only the two files above
  implement it. `tests/test_polarization_states.py`'s `_UNIT_POWER_STATES`
  hardcodes numeric `(s,p)` pairs directly (not via this formula) to test
  an unrelated invariant (normal-incidence `R`/`T` independence of
  polarization state for *any* `(s,p)` pair at fixed power) -- its
  `"linear_20deg"` label is cosmetic there and needed no change.
  `src/sougata_solver/excitation.py`, `polarimetry.py`, `sweep.py`, and
  `plotting.py` all take already-resolved `s_amplitude`/`p_amplitude`
  directly and contain no alpha-to-amplitude formula -- no core solver
  change. No prior ADR existed for this convention; this is the first.
- **Impact**: `CONVENTIONS.md`, `structures/thin_film/custom_multistack.py`,
  `structures/thin_film/sio2_on_si_thin_film.py`. No
  `src/sougata_solver/` change and no test assertions changed (the physical
  S/P reflectance computation itself was already correct and already
  matched the commercial tool to ~0.1%; only the linear/elliptical
  angle-labeling convention changed).

## ADR-034: `multistack_composite_grating.py` cross-validated against Lumerical RCWA — the mismatch was a semi-infinite-vs-finite substrate difference, not a materials or formula bug

- **Decision**: `structures/thin_film/multistack_composite_grating.py` (a
  laterally-alternating composite of two multilayer stacks -- Si/SiO2 on
  one half of a 2 um period, Ni/SiO on the other, 1D-periodic via
  `Lattice1D`/`Slab`, reproducing a structure the project owner built in
  Lumerical's RCWA solver) is now cross-validated against that Lumerical
  model to ~1% agreement: `R` max\|diff\|=0.013 (RMS 0.0045), `T`
  max\|diff\|=0.012 (RMS 0.0022), via a new
  `postprocessing/overlay_composite_grating_vs_lumerical.py` script that
  overlays this solver's `output_*_RT.csv` against Lumerical's own
  `grating_power` result (`Rs_power`/`Ts_power`, summed over all
  diffraction orders -- see that result's `n`/`m` order-index dimensions,
  confirmed via `size(Rs_power)` before trusting the sum).
- **What the mismatch actually was**: the first comparison attempt showed a
  large, structurally-shaped gap (`R` max\|diff\|=0.28, `T`
  max\|diff\|=0.69 -- absolute units, not percent) that survived a direct
  materials-data cross-check: this solver's `NK_FILE/*_KLA.txt`
  permittivity for Si/Ni/SiO2/SiO matched Lumerical's Palik-database fit
  closely at every sampled wavelength (e.g. Si `eps` 30.9/13.6 at
  400/800nm here vs. Lumerical's plotted ~30/~14; Ni `eps` -3.0/-13.1 vs.
  Lumerical's ~-2.5/~-12.5) -- ruling out materials as the cause. Root
  cause: the original Lumerical geometry drew `Si_substrate`/`Ni_substrate`
  as 5 um-deep rectangles while the RCWA computation region's own z-extent
  only reached -0.5 um. Per Lumerical's documented behavior (incidence/
  transmission media are inferred from whatever object the simulation
  region's z-boundary extends into, not from the Interfaces tab list),
  Lumerical was resolving a genuinely semi-infinite, laterally-patterned
  Si(left)/Ni(right) exit medium -- not "finite Si/Ni over air" like this
  solver's `TRANSMISSION_MATERIAL = air` model. This solver's `Simulation`
  requires `transmission=` to be one uniform `Material`; it has no way to
  represent "Si going down forever on the left, Ni going down forever on
  the right" simultaneously.
- **Confirmed, not just inferred**: a controlled thickness sweep on this
  solver's side (Si/Ni at 0.5/2/5/10 um over semi-infinite air) showed `R`
  had not converged even at 10 um for the weakly-absorbing 800nm case
  (0.405 -> 0.345 -> 0.215 -> 0.243, non-monotonic) -- genuine Fabry-Perot
  interference from the buried Si/Ni-to-air interface, matching
  `structures/README.md`'s already-documented thick-Si fringing note. This
  proves a "thick finite layer" cannot stand in for a true semi-infinite
  half-space either (at wavelengths where the material isn't strongly
  absorbing), so the fix had to be structural, not a thickness tweak.
- **Resolution**: the project owner changed the Lumerical model's
  `Si_substrate`/`Ni_substrate` objects to genuinely finite 0.5 um slabs
  (`z: -0.5` to `0`, matching what `multistack_composite_grating.py`
  already built) with real air below (RCWA region `z min` nudged to -0.6
  um so the computation actually samples that air rather than sitting
  exactly on the object boundary), then re-ran the RCWA solve. That
  produced the ~1%-agreement result above; the residual gap is consistent
  with RCWA harmonic-order truncation and Lumerical's multi-coefficient
  dispersion-fit vs. this solver's tabulated-`n,k` interpolation, not a
  remaining structural error.
- **Scope note, not a bug**: a genuinely semi-infinite, laterally-patterned
  substrate (the *original* Lumerical model) is a real capability gap in
  this solver today -- flagged here rather than worked around, per
  `rules.md` AI Coding Rule 2. Representing that exact structure would
  need new solver work (a patterned semi-infinite half-space eigenbasis as
  the exit boundary condition), not a scripting change.
- **Impact**: new `postprocessing/overlay_composite_grating_vs_lumerical.py`
  (R/T overlay against a Lumerical `grating_power` export, linear-
  interpolated onto this solver's wavelength grid for a max/RMS diff,
  robust to Lumerical's `write()` append-not-overwrite behavior via
  last-header-block detection -- see that script's
  `_load_lumerical_grating_power`). No `src/sougata_solver/` change.
  `structures/thin_film/multistack_composite_grating.py` itself needed no
  change -- it was correct from the start; the mismatch was entirely in
  the Lumerical reference model.
