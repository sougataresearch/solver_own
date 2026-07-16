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
