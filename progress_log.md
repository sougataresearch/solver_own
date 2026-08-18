# Progress Log — sougata_solver

Append-only, date-stamped log of discussions and the action items that came
out of them. This is **not** a replacement for `tasks.md` (phase-organized
checklist) or `memory.md` (living project-status snapshot) — it's the
chronological record of *what was discussed and why*, so a session that
starts cold (human or AI) can see what was raised, what was decided, and
whether it was ever actually implemented.

## How to use this file

- **At the start of any substantive session**, read the most recent entries
  first (bottom of file) and check the status of any `[ ]` open item —
  search the codebase to verify whether it's actually been implemented
  before assuming it's still pending. Update `[ ]` → `[x]` (with the date
  it was completed and where) once verified done. Never mark something done
  without checking the actual code/tests.
- **At the end of any substantive session** (discussion, decision, or code
  change), add a new dated entry below the previous one. Don't edit past
  entries except to flip a checkbox's status — history should stay visible,
  not get rewritten.
- If an item graduates into a real tracked task, add it to `tasks.md`/
  `phases.md` too and cross-reference it here (e.g. `-> tasks.md Phase 3`)
  so it isn't tracked in two disconnected places.
- If an item was discussed and explicitly resolved as "no action needed"
  (with reasoning), log it as `[x] (no action needed — see reasoning)`
  rather than leaving it open or deleting it — that prevents re-litigating
  the same question in a future session.

## Entry format

```
## YYYY-MM-DD

### Discussed
- Bullet summary of each topic raised.

### Action items
- [ ] Thing to implement, with enough context to act on later.
- [x] Thing resolved / already done — note where.
```

---

## 2026-07-19

### Discussed
- Field-convention question (`exp(-jkz)` vs `exp(+jkz)`, engineering vs
  physics sign convention) for hand derivations vs what's baked into the
  code — confirmed the codebase (`smatrix.py:108`, transcribed from S4)
  uses the physics convention (`exp(+jqz)`, `d/dt -> -jw`); hand
  derivations following the user's textbook may use the opposite
  convention. No code change — just be deliberate (`j -> -j`) when porting
  a hand-derived formula into the code.
- How polarization/transversality is generated in the solver
  (`excitation.py`) and why `Ez`/`Hz` can never appear (no such degree of
  freedom exists in the RCWA transverse-field formulation).
- How incidence direction (top vs bottom of stack) is controlled
  (`simulation.py:109-111`, drives `a_left` at `layer_stack[0]`) — confirmed
  top-down is the only implemented direction, and that's the intended
  default (user: "usually we incident from top of the sample").
- Whether `Layer` needs explicit x/y (lateral) extent for the SiO2-on-Si
  thin-film structure — confirmed uniform layers are laterally infinite by
  RCWA's formulation; x/y only becomes meaningful for patterned layers
  (`Layer.pattern`), which isn't implemented yet (`simulation.py:97-98`
  raises `NotImplementedError`).
- Investigated a real discrepancy between two `sio2_on_si_thin_film.py`
  output plots (`outputs/2026_07_16/11_06_37.../output_RT.png` vs
  `11_10_03.../output_RT.png`) — root cause confirmed via
  `run_metadata.txt`: identical structure, only difference was wavelength
  sampling (41 pts vs 401 pts). The dense oscillation in the finer-sampled
  plot is real Fabry-Perot interference from the 12 um Si layer
  (`Delta_lambda ~ lambda^2/(2 n t) ~ 6 nm` near 750 nm), not a bug; the
  41-point run was aliased/undersampled and produced a misleadingly smooth
  (wrong) curve.
- User raised three doubts about the above: (1) short-substrate
  bottom-reflection interference, (2) whether solver boundaries need
  PML-style absorbing boundaries, (3) whether surrounding medium affects
  results. Checked against actual code:
  1. Real physics, already handled exactly by the Redheffer star-product
     cascade (`smatrix.py` `star_product`/`SMatrixStack`) — no fix needed.
  2. Not applicable — RCWA has no discretized/truncated z-domain
     (`incidence`/`transmission` are literal `thickness = math.inf`
     half-spaces, `layer.py:55,57`), so there's nothing for PML to fix.
     PML is an FDTD/FEM concept for terminating a finite mesh; RCWA solves
     each layer analytically instead.
  3. Yes, surrounding medium matters and is already exposed as
     `Simulation(incidence=..., transmission=...)` — currently air/air in
     the script, matching the intended free-standing setup.
- Clarified `sougata_solver` is RCWA (Fourier Modal Method), not FEM —
  no mesh, no discretized domain, hence no PML.
- Discussed layer *slicing* (staircase discretization) as used in
  Lumerical RCWA: needed only for structures whose in-plane cross-section
  changes continuously with z (slanted sidewalls, graded index) — not
  needed for the current flat SiO2/Si stack. Confirmed this is **already
  captured** in `decisions.md` ("Tapered sidewalls via staircase
  discretization, not new Fourier math") and blocked on Phase 2+ patterned
  layers (`tasks.md` Phase 4) landing first; `SMatrixStack` already
  cascades an arbitrary number of layers today, so no change needed there
  once per-slice patterned eigenmodes exist.
- Set up this file (`progress_log.md`) itself, at user's request, as a
  dated, checkable discussion/action-item log distinct from the personal
  Claude-Code memory system (which is cross-project/cross-session for the
  AI assistant, not repo-local project documentation).
- User asked for a Claude Code skill that, before coding any phase, checks
  `REFERENCE/` and picks the best-fit vendored repo per structure type
  (thin film, multilayer, 1D trench, via, etc.) instead of one repo being
  reused everywhere. Added `phase-reference-picker`
  (`.claude/skills/phase-reference-picker/` at the workspace root). First
  draft implicitly leaned on S4 as the default answer for most phases;
  user caught this and asked for a real per-sub-task comparison plus an
  explicit choice between transcribing a source vs. deriving independently
  and using it only as a cross-check — skill rewritten accordingly, and
  `references.md`/`rules.md`/this file cross-referenced to it.

- User asked how to make the roadmap more scientifically/mathematically
  rigorous — whether to split phases further or add more complex
  requirements. Recommended, and (on approval) implemented: (1) split
  Phase 4 into 4a (well-conditioned via/pillar) and 4b (near-degenerate/
  ill-conditioned stress cases), since `design.md` already called the
  general eigendecomposition "the highest-risk remaining algorithm" but
  the original single Phase 4 let an easy passing test case stand in for
  that whole risk; (2) added a new oracle-independent Physical-Invariant
  Testing tier to `testing.md` (energy conservation, convergence-rate-vs-
  Li-1996-theory), required starting Phase 3, layered on top of (not
  replacing) the existing oracle-comparison requirement. Deliberately did
  *not* add new out-of-scope physics (magneto-optic, arbitrary polygons)
  since `PRD.md` already excludes those by explicit prior decision.
  Updated: `phases.md`, `tasks.md`, `PRD.md`, `testing.md`, `rules.md`,
  `design.md`, `memory.md`.

### Action items
- [x] Add `phase-reference-picker` skill and cross-reference it from
  `references.md` ("Choosing a Reference for a New Phase"), `rules.md`
  (AI Coding Rule 8), and `memory.md` (Things Future AI Sessions Should
  Remember) — done 2026-07-19.
- [x] Split Phase 4 into 4a/4b and add Physical-Invariant Testing tier —
  done 2026-07-19, see `phases.md`/`tasks.md`/`PRD.md`/`testing.md`.
- [x] (no action needed — see reasoning) Sign convention in code vs hand
  calculation — no mismatch to fix, just keep them separate.
- [x] (no action needed — see reasoning) Transversality of E/H fields — 
  already structural, nothing to add.
- [x] (no action needed — see reasoning) Incidence direction — top-down is
  already the implemented and intended default.
- [x] (no action needed — see reasoning) Lateral (x/y) extent for uniform
  layers — not applicable until patterned layers exist.
- [x] (no action needed — see reasoning) PML / absorbing boundaries —
  not applicable to RCWA's analytic half-space formulation.
- [x] (no action needed — see reasoning) Surrounding-medium sensitivity —
  already exposed via `incidence`/`transmission` materials, working as
  intended.
- [ ] **Layer slicing (staircase discretization) for slanted/graded
  patterned layers** — not yet implementable at all (blocked on Phase 2+
  patterned-layer support, `tasks.md` Phase 4, `simulation.py:97-98`).
  Already recorded as a design decision in `decisions.md`. Revisit once
  Phase 4 (2D-periodic patterned layers) lands — check `tasks.md` Phase 4
  checklist status first.

---

## 2026-07-21

### Discussed
- Reviewed Phase 4a completion status from scratch: `solve_layer_eigenmodes_patterned`
  (dense `2n x 2n`, `rcwa.cpp:794-827`), `simulation.py` 2D dispatch, and
  `structures/via/pillar_array.py` / `via_array.py` were already implemented.
  `tasks.md` already showed Phase 4a as ☑ complete.
- Verified all 98 non-slow tests pass (including all 9 `test_2d_pillar.py` tests:
  reduce-to-uniform, ky=0→1D, 6 energy-conservation parametrize cases, end-to-end).
- Found a signature bug in both `structures/via/pillar_array.py` and
  `via_array.py`: `write_run_metadata(out, {dict})` was passing a dict as the
  `script_path` argument instead of `__file__` — caught because
  `output_paths.write_run_metadata` takes `(output_dir, script_path: str, **params)`
  not `(output_dir, dict)`. Scripts ran without raising (Python happily wrote the
  dict's `str()` as the script path) but recorded wrong metadata.
- Both structure scripts were also not doing a wavelength sweep (only one wavelength
  point) and had no console R/T print table, inconsistent with `trench_grating.py`.
- Fixed both scripts: correct `write_run_metadata` call, 21-point wavelength sweep
  (500–1500 nm), console R/T table matching `trench_grating.py` style.
- Ran both scripts post-fix: R+T = 1.0000 at all 21 points for both pillar and via.
- Noted that `memory.md`/`progress_log.md` had not been updated to reflect Phase 4a
  completion (still said "Phase 4a is next") — updated both this session.
- The `staircase discretization for slanted/graded patterned layers` item
  (progress_log.md 2026-07-19) was blocked on Phase 4 landing — Phase 4a is now
  done, so this is unblocked; it belongs in Phase 5 (`tasks.md` Phase 5).

### Action items
- [x] Fix `write_run_metadata` call signature in `pillar_array.py` and
  `via_array.py` — done 2026-07-21.
- [x] Add wavelength sweep + console R/T table to both via structure scripts —
  done 2026-07-21 (matches `trench_grating.py` pattern).
- [x] Update `memory.md` to reflect Phase 4a completion — done 2026-07-21.
- [ ] **Phase 4b** (near-degenerate/ill-conditioned stress cases) — next phase.
  See `tasks.md` Phase 4b for the checklist. Start with: high-contrast Si pillar
  at small radius/period ratio and high `num_orders`, add condition-number
  `WARNING` logging to `eigenmodes.py`.
- [ ] **Layer slicing (staircase discretization)** — unblocked by Phase 4a;
  now trackable as Phase 5 (`tasks.md` Phase 5). Check `tasks.md` Phase 5
  checklist before starting.

---

## 2026-08-01

### Discussed
- User shared a Lumerical/FDTD structure-group script building a tapered
  1D grating parametrized by `tcd`/`bcd`/`depth`/`spacing`/`zSpan`/
  `yCompensation`/`zCompensation`/`grating_number`, and asked what each
  parameter meant, then asked to build the equivalent trench/via/pillar
  structures with this parameter style.
- Mapped FDTD params to Phase 5's already-shipped `staircase.py` generators:
  `tcd`/`bcd` -> `top_halfwidth`/`bottom_halfwidth` (trench) or
  `top_radius`/`bottom_radius` (via) or square-pillar side length;
  `depth` -> `thickness`; `spacing` -> derives `period = tcd + spacing`.
  `zSpan`, `yCompensation`/`zCompensation`, and `grating_number` have no
  sougata_solver equivalent (no finite 3D mesh/extrusion, no absolute
  z-origin offset, and `Lattice`/`Lattice1D` are infinite-periodic so no
  finite-array replication is needed).
- Confirmed with the user (via question) to keep the existing per-shape
  example-script pattern rather than add a new shared
  tcd/bcd-to-halfwidth conversion helper module — Phase 5 didn't scope a
  new shared module, and the existing scripts already work fine, so
  renaming constants in place is the smaller, correctly-scoped change.
- Renamed `structures/trench/tapered_trench.py` and
  `structures/via/tapered_via.py`'s top/bottom-size constants to
  `TCD`/`BCD`/`SPACING`, `PERIOD` derived. Added
  `structures/via/tapered_pillar.py` (the `Rectangle` case of
  `staircase_rectangle_layers`, square pillar with equal x/y halfwidths) —
  this shape type had generator support in `staircase.py` since Phase 5 but
  no example script yet. All three re-run end-to-end: R+T=1.0000 at every
  `num_slices` value.
- User then asked to commit and push to
  `https://github.com/sougataresearch/solver_own.git`. Found the repo had
  a large uncommitted backlog beyond today's change: `git log` showed the
  last commit only covers through Phase 2, while `phases.md` records
  Phases 3-8 as done — all of that source/tests/structures was sitting
  uncommitted. Also found two untracked stray files not covered by
  `.gitignore`: `demo.fsp` (25MB Lumerical binary) and `OUTPUT_RCWA/` (old
  run-output artifacts, same category `outputs/` already ignores). Asked
  the user before including either; user chose to exclude both. Added
  `OUTPUT_RCWA/` and `demo.fsp` to `.gitignore`.

### Action items
- [x] Rename tapered trench/via constants to FDTD-style tcd/bcd/spacing —
  done 2026-08-01.
- [x] Add `structures/via/tapered_pillar.py` — done 2026-08-01.
- [x] Add `OUTPUT_RCWA/`/`demo.fsp` to `.gitignore` — done 2026-08-01.
- [x] Commit and push the full backlog (Phases 3-8 plus today's changes) to
  `origin/main` — done 2026-08-01 (`f5fd81b`).

---

## 2026-08-03

### Discussed
- User asked to complete `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 1
  targets 1.3-1.8 (Phase 6 anisotropic-materials work) one by one. Read
  `rules.md`, `phases.md`, `memory.md`, `progress_log.md`, `references.md`,
  `CONVENTIONS.md`, `tasks.md`, and the current `eigenmodes.py`/
  `materials.py`/`simulation.py`/`layer.py` before writing any code, per
  `CLAUDE.md`'s workspace instructions. Planned the six targets (plan mode,
  approved) and started with 1.3.
- Implemented target 1.3 (uniform diagonal-tensor layers) — see `memory.md`'s
  new entry for the full account, including a caught-before-shipping test
  mistake (initial oracle-axis mapping was backwards; the test itself caught
  it, not a separate review pass) around the existing `Epsilon2` block-index
  convention (`CONVENTIONS.md`'s `u = [-Ey; Ex]` ordering).

### Action items
- [x] Target 1.3 (uniform diagonal tensor) — done 2026-08-03, see
  `memory.md`, `tasks.md` Phase 6, `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`.
- [x] Target 1.4 (in-plane tensor coupling: `eps_xx, eps_xy, eps_yx, eps_yy,
  eps_zz`) — done 2026-08-03; oracle hand-transcribed from
  `RigorousCoupledWaveAnalysis.jl`'s `AnisotropicLayer` path (`Common.jl:134-165`)
  into `tests/oracles/rcwa_anisotropic_inplane_jl.py`. See `memory.md` for
  the empirically-determined kx/ky-swap-and-negate convention finding and
  two test-authoring mistakes caught (tolerance scale, non-Hermitian
  energy-conservation material) before this was marked done.
- [x] Target 1.5 (longitudinal coupling) — evaluated and explicitly deferred
  2026-08-03; bounded `WebSearch`/`WebFetch` literature search found no
  source both readable in this environment and independently
  benchmarkable. See `references.md`'s "Target 1.5 bounded literature
  search" entry and `memory.md`.
- [x] Target 1.6 (patterned anisotropic layers) — done 2026-08-03. New
  citation found this session: `S4/S4/fmm/fmm_closed.cpp`'s `have_tensor`
  branch (lines 165-256), not covered by the original Phase 6 reference
  audit. See `memory.md` for the full account.
- [x] Target 1.7 (degeneracy policy) — done 2026-08-03:
  `eigenmodes._canonical_mode_order`, scoped to the three anisotropic
  dense eigensolvers only (not the pre-existing Phase 4a isotropic
  solver). See `memory.md`.
- [x] Target 1.8 (mode classification) — done 2026-08-03, last target in
  this session's approved plan. Found (and documented, not fixed — out of
  this target's scope) an exact-Rayleigh-threshold `NaN` division-by-zero
  in `smatrix.py`, tied to Category 6 target 6.4. See
  `troubleshooting.md` and `memory.md`.
- **Session summary**: all six planned targets (1.3-1.8) addressed;
  1.3/1.4/1.6/1.7/1.8 shipped, 1.5 explicitly deferred after a bounded
  literature search. 186 tests pass project-wide (123 at session start,
  65 new). Committed as two commits on top of the `f5fd81b`/`9b28627`
  backlog already pushed to `origin/main` on 2026-08-01, rebased via
  `git cherry-pick` after discovering the local checkout had independently
  (and redundantly) re-committed the same Phase 3-5 material already on
  the remote — see `memory.md`.

---

## 2026-08-04

### Discussed
- User asked to complete `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 2
  targets 2.1-2.5 (Numerical methods) one by one. Read `rules.md`,
  `phases.md`, `architecture.md`, `memory.md`, `testing.md`, `design.md`,
  `troubleshooting.md`, and the current `eigenmodes.py`/`layer.py`/
  `smatrix.py`/`simulation.py` before writing any code, per `CLAUDE.md`'s
  workspace instructions.

### Action items
- [x] Target 2.1 (failure contract) — done 2026-08-04. New "Failure
  Contract" section in `design.md`, `tests/test_failure_contract.py` (17
  tests). See `memory.md`.
- [x] Target 2.2 (eigenvalue report) — done 2026-08-04.
  `layer.EigenmodeDiagnostics`, `LayerEigenmodes.diagnostics` (new optional
  field). `tests/test_eigenvalue_diagnostics.py` (6 tests). See `memory.md`.
- [x] Target 2.3 (sweep mode matching) — done 2026-08-04, scoped to the
  three anisotropic dense solvers after an attempted extension to Phase
  4a's isotropic solver broke an existing regression test (reverted, not
  the test relaxed — see `memory.md` and
  `eigenmodes.solve_layer_eigenmodes_patterned`'s docstring).
  `tests/test_sweep_mode_matching.py` (4 tests).
- [x] Target 2.4 (degeneracy warning) — done 2026-08-04.
  `eigenmodes.DEGENERATE_GAP_THRESHOLD`/`_warn_on_small_eigenvalue_gap`,
  same three solvers as 2.3 (found during testing that the Phase 4a
  isotropic solver has routine, harmless `C4v`-symmetry near-degeneracy
  that would make the warning noisy there, so deliberately excluded).
  `tests/test_degeneracy_warning.py` (5 tests).
- [x] Target 2.5 (stress regression) — done 2026-08-04, with a
  sign-convention finding: Phase 4b's `n=-20+2j` "lossy-metal-like" index
  is actually a **gain** medium under this project's phasor convention
  (`Im(eps)<0`), only discoverable by actually calling `Simulation.solve()`
  end to end (which Phase 4b's own eigenvalue-only stress test never did)
  — fixed by using a correctly-signed lossy metal for the new fixture, not
  by touching Phase 4b's already-shipped file. `tests/test_stress_regression.py`
  (2 tests). See `memory.md` for the full account.
- **Session summary**: all five Category 2 targets (2.1-2.5) complete.
  227 tests pass project-wide (186 fast at session start, 34 new fast + 7
  unchanged `slow`). No existing oracle-comparison/regression test was
  weakened to make a new one pass (`rules.md` AI Coding Rule 3) — one
  attempted change (2.3's canonical-ordering extension) was reverted
  instead, and documented as a negative finding rather than silently
  dropped. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`,
  `tasks.md` (new cross-cutting section pointing at the atomic-targets
  file, since Category 2 isn't phase-scoped), and this file.

### Discussed (same day, continued)
- User asked to complete Category 3 (Fourier factorization) targets 3.1-3.6
  the same way. Read `design.md`'s existing Algorithm 3, `references.md`,
  `fourier_factorization.py`, and the existing 1D/2D convergence tests
  before starting, per `CLAUDE.md`'s workspace instructions.

### Action items (Category 3)
- [x] Target 3.1 (rule inventory) — done 2026-08-04. New "Fourier-
  factorization rule inventory" table in `design.md` (Algorithm 3a),
  `tests/test_fourier_factorization_rules.py` (6 tests). Finding:
  `epsilon_inv_hat` is only actually consumed as an inverse-rule Toeplitz
  in the 1D TM block; every 2D path uses a numerical inverse of the
  direct-rule matrix instead. See `memory.md`.
- [x] Targets 3.2/3.3 (1D/2D convergence fixtures) — done 2026-08-04.
  `tests/test_fourier_convergence.py` (2 `slow` tests), high-contrast
  `n=10` 1D grating and `n=5` 2D pillar, with actually-measured (not
  assumed) convergence tables recorded in the test docstrings. Notable
  finding: the 2D fixture's `num_orders=25` point is an order-of-magnitude
  non-monotonic outlier — real, measured evidence for exactly the
  weakness targets 3.4/3.5 investigate. See `memory.md`.
- [x] Targets 3.4/3.5 (FFF/NVM feasibility) — evaluated and explicitly
  deferred, 2026-08-04. Bibliographic details for Popov & Nevière (2001)
  and Lalanne (1997) confirmed via `WebSearch`; neither paper's full text
  was fetchable (paywalled). `../REFERENCE/S4` read in full instead:
  `fmm_PolBasisNV.cpp`/`fmm_PolBasisJones.cpp`/`fmm_PolBasisVL.cpp`
  (~900 lines combined), all built on a discretized/FFT permittivity
  representation that conflicts with this project's already-shipped
  ADR-002 (analytic Fourier transforms, raster+FFT rejected). See
  `decisions.md` ADR-012 and `references.md`.
- [x] Target 3.6 (selected improvement) — done 2026-08-04 (no action
  needed — both 3.4 and 3.5 concluded defer, so nothing was approved to
  implement; recorded as this target's own explicit outcome, not skipped).
- **Session summary**: all six Category 3 targets (3.1-3.6) complete.
  232 tests pass project-wide (227 at the start of this sub-session: 220
  fast + 7 slow -- 226 fast + 9 slow now). No existing test weakened.
  Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `references.md`, `decisions.md` (new ADR-012), and this file.

### Discussed (same day, continued again)
- User asked to complete Category 4 (Geometry engine) targets 4.1-4.7 the
  same way. Read `geometry.py`, `decisions.md` ADR-005 (parametric-shapes-
  only), `staircase.py`, and `test_fourier_factorization.py`'s existing
  shape-test conventions before starting, per `CLAUDE.md`'s workspace
  instructions. Noted upfront that targets 4.4/4.5/4.6 (polygon, imported
  geometry) appear to directly conflict with ADR-005 -- resolved by reading
  ADR-005's own text, which explicitly anticipates and permits exactly this
  situation ("if a real need for imported layouts arises later, revisit
  with a new ADR, not silently override") -- see ADR-013 below.

### Action items (Category 4)
- [x] Target 4.1 (geometry validation API) — done 2026-08-04.
  `geometry._require_finite`/`_require_positive`, wired into
  `Lattice`/`Lattice1D`/`Circle`/`Rectangle`/`Slab`.
  `tests/test_geometry_validation.py` (29 tests). See `memory.md`.
- [x] Target 4.2 (unit-cell bounds policy) — done 2026-08-04.
  `geometry.validate_pattern_fits_lattice`, wired into `Simulation.__init__`.
  `tests/test_unit_cell_bounds.py` (6 tests). See `memory.md`.
- [x] Target 4.3 (Ellipse primitive) — done 2026-08-04. Transcribed from
  `S4/S4/pattern/pattern.c`'s `ELLIPSE` case. `tests/test_ellipse.py`
  (19 tests), `structures/via/elliptical_pillar.py`. See `memory.md`.
- [x] Target 4.4 (Polygon design decision) — done 2026-08-04. Found S4's
  own analytic (not raster/FFT) polygon Fourier transform while
  investigating this target -- decided analytic, recorded as `decisions.md`
  ADR-013 (a narrow, explicit revisit of ADR-005, not a silent override).
- [x] Target 4.5 (Polygon primitive) — done 2026-08-04. `geometry.Polygon`,
  reduces exactly to `Rectangle` for a square; `Polygon.signed_distance_normal`
  independently derived rather than transcribed from S4's own version,
  which was found to pick the farthest (not nearest) edge -- see
  `memory.md`. `tests/test_polygon.py` (24 tests),
  `structures/via/triangular_pillar.py`.
- [x] Target 4.6 (imported geometry format) — done 2026-08-04.
  `geometry_io.py`, JSON-only, no `eval`/`exec`, parser-only (not wired
  into `Simulation`/`Layer`, per the target's own scoping).
  `tests/test_geometry_io.py` (17 tests).
- [x] Target 4.7 (profile slicing API) — done 2026-08-04.
  `staircase.slice_profile`; the three existing taper generators
  refactored into thin wrappers around it, regression-verified unchanged
  against the full pre-existing `tests/test_staircase.py` suite.
  `tests/test_profile_slicing.py` (6 tests).
- **Session summary**: all seven Category 4 targets (4.1-4.7) complete.
  336 tests pass project-wide (232 at the start of this sub-session: 226
  fast + 9 slow -- 327 fast + 9 slow now, 101 new fast tests). No existing
  test weakened; the `staircase.py` refactor was verified, not assumed,
  behavior-preserving. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`,
  `memory.md`, `tasks.md`, `references.md`, `decisions.md` (new ADR-013),
  `architecture.md` (module table + directory listing), and this file.

### Discussed (same day, continued once more)
- User pointed out that `structures/via/elliptical_pillar.py`/
  `triangular_pillar.py` reported `R+T=1.0000` rather than `R+T+A=1.0000` --
  confirmed this is correct (both scripts use real, lossless `n=3.48`/`1.0`
  materials, so `A=0` by Poynting's theorem, not by omission), but surfaced
  a real, separate gap: no `absorbance()`/`A` computation exists anywhere
  in the library yet (confirmed by grep, not assumed) -- tracked as
  Category 7 targets 7.5/7.6, still open. User then asked to do Category 5
  (Material models) targets 5.1-5.8 the same way. Read `materials.py` in
  full, and searched all vendored `REFERENCE/` repos for Sellmeier/Cauchy/
  Lorentz/Drude implementations before starting, per `CLAUDE.md`'s
  workspace instructions -- found genuinely useful, previously-unused
  sources: `EMpy/EMpy/materials.py` (Sellmeier/Cauchy),
  `Rigorous-Coupled-Wave-Analysis/TMM_examples/TMM_Drude.py` (Drude), and
  `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl` (a full,
  published Lorentz-Drude metal model with real Au/Ag/Al/Ti coefficients,
  citing Rakić et al. 1998).

### Action items (Category 5)
- [x] Target 5.1 (material validation) — done 2026-08-04. Construction-
  *and* call-time validation on `Material`. `tests/test_material_validation.py`
  (13 tests). See `memory.md`.
- [x] Targets 5.2/5.3 (Sellmeier/Cauchy) — done 2026-08-04, transcribed
  from `EMpy/EMpy/materials.py`. BK7 validated against an independently-
  published `n_d=1.5168` (confirmed via `WebSearch`). 8 tests.
- [x] Targets 5.4/5.5/5.6 (Lorentz/Drude/Drude-Lorentz) — done 2026-08-04,
  transcribed from `rakic.jl` (Rakić et al. 1998). Causality/sign-
  convention independently re-derived and confirmed correct *before*
  shipping (unlike Category 2 target 2.5's after-the-fact catch). Drude
  formula cross-checked between two independently vendored sources. A
  bounded search for a second independent metal-optics reference (Johnson
  & Christy 1972) found the citation but not the fetchable data — recorded
  honestly in `references.md`, not silently skipped. 17 tests.
- [x] Target 5.7 (tensor-material solver wiring) — done 2026-08-04. Gate
  (Category 1 targets 1.3/1.4/1.6) confirmed already met; closed the
  previously-untested dispersive-tensor-material x tensor-eigensolver
  combination. 6 tests.
- [x] Target 5.8 (material provenance) — done 2026-08-04. `Material.source`
  threaded through every classmethod and `geometry_io`'s JSON schema.
  **Found and fixed a real pre-existing bug**: `write_run_metadata` lacked
  explicit UTF-8 encoding, raising `UnicodeEncodeError` on Windows the
  first time non-ASCII citation text was actually written through it. 13
  tests (11 + 2 in `test_geometry_io.py`).
- **Session summary**: all eight Category 5 targets (5.1-5.8) complete.
  393 tests pass project-wide (336 at the start of this sub-session: 327
  fast + 9 slow -- 384 fast + 9 slow now, 57 new fast tests). No existing
  test weakened; one genuine bug found and fixed along the way
  (`output_paths.write_run_metadata`'s encoding). Updated
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `references.md`, and this file.

---

## 2026-08-04 (continued)

### Discussed
- User asked to complete Category 6 (Boundary conditions and excitation)
  targets 6.1-6.6 the same way. Read `excitation.py`, `CONVENTIONS.md`
  (already covers most of target 6.1), `smatrix.py`, and existing
  polarization/mode-classification tests before starting, per `CLAUDE.md`'s
  workspace instructions.

### Action items
- [x] Target 6.1 (convention note) — done 2026-08-04. Confirmed
  `CONVENTIONS.md` already satisfies this; added a worked-examples table.
- [x] Targets 6.2/6.3 (normal/oblique polarization regression) — done
  2026-08-04. Two symmetry invariants verified numerically then encoded:
  polarization-state independence at normal incidence, azimuthal-rotation
  invariance at oblique incidence. `tests/test_polarization_states.py`
  (89 tests). See `memory.md`.
- [x] Target 6.4 (grazing-incidence boundary) — done 2026-08-04.
  Characterized directly: any `theta<90deg` supported; exactly `90deg`
  raises `ValueError` (not `NaN`), traced to a real floating-point
  coincidence. `tests/test_grazing_incidence.py` (9 tests). Also backfilled
  several missing Category 4/5 rows into `design.md`'s Failure Contract
  while there.
- [x] Target 6.5 (Rayleigh-threshold test) — done 2026-08-04 (found already
  satisfied by Category 1 target 1.8's `tests/test_mode_classification.py`
  at normal incidence; added the oblique-incidence case that coverage
  lacked). `tests/test_oblique_rayleigh_threshold.py` (8 tests).
- [x] Target 6.6 (bottom-incidence decision) — done 2026-08-04, with a
  better answer than "defer": already achievable via the existing
  `Simulation` constructor, verified via Stokes transmittance reciprocity
  (`~1e-15` agreement), not just asserted. `decisions.md` ADR-014,
  `tests/test_bottom_incidence.py` (3 tests).
- **Session summary**: all six Category 6 targets (6.1-6.6) complete. 502
  tests pass project-wide (393 at the start of this sub-session: 384 fast
  + 9 slow -- 493 fast + 9 slow now, 109 new fast tests). No existing test
  weakened. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`,
  `tasks.md`, `decisions.md` (new ADR-014), `CONVENTIONS.md`, `design.md`,
  and this file.

## 2026-08-05

### Discussed
- Completing Phase 7 (Real-Space Field Reconstruction & Visualization),
  tracked at atomic-target grain as `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
  Category 9 (targets 9.1-9.8), then committing and pushing both this
  category and the still-uncommitted Category 6 work from the prior
  sub-session.
- How to reconstruct real-space E/H fields from per-order modal Fourier
  coefficients: transcribed S4's `GetInPlaneFieldVector`/`GetFieldAtPoint`
  formulas (transverse + longitudinal components), citing exact line
  ranges in `references.md`.
- How to recover interior-layer mode amplitudes without porting S4's more
  complex block-tridiagonal `SolveInterior` algorithm: derived a simpler
  formula from this project's own already-implemented
  `SMatrixStack.partial_smatrix_up_to`, using standard Redheffer
  star-product algebra — recorded as `decisions.md` ADR-015 (independently
  derived, not transcribed from S4).
- A genuine factor-of-2 discrepancy found while validating flux-matches-R/T:
  this project's established `fields.z_poynting_flux` modal quadratic form
  is exactly 2x the textbook real-space `Sz = 0.5*Re(Ex*conj(Hy) -
  Ey*conj(Hx))` — harmless for R/T (ratio, cancels) but must be accounted
  for when computing absolute real-space flux from raw reconstructed
  fields. Not a bug in `z_poynting_flux` itself; documented as an
  established convention.
- `structures/` vs `postprocessing/` boundary (ADR-009/010) applied to two
  new example scripts (trench (x,z) cross-section, pillar (x,y) field map)
  and one new plotting script — both example scripts run end-to-end and
  their output PNGs visually inspected for physical sensibility.

### Action items
- [x] Target 9.1 (field-component formulas) — done 2026-08-05.
  `fields.modal_field_components`, transcribed from S4's
  `GetInPlaneFieldVector`/`GetFieldAtPoint`. See `memory.md`, `CONVENTIONS.md`.
- [x] Targets 9.2/9.3 (depth propagation, interior amplitudes) — done
  2026-08-05. `fields.propagate_amplitudes`, `smatrix.interior_amplitudes`
  (independently derived, `decisions.md` ADR-015), validated by a
  zero-free-parameter self-consistency check against the already-known
  transmitted amplitude.
- [x] Targets 9.4-9.6 (reconstruction + flux-match tests) — done 2026-08-05.
  `tests/test_field_reconstruction.py` (10 tests): analytic plane-wave
  match, transversality, interface continuity, 1D periodicity, and
  flux-matches-R/T (found and documented the missing-0.5-factor
  convention along the way).
- [x] Targets 9.7/9.8 (example scripts + grid save) — done 2026-08-05.
  `fields.save_field_grid_npz`, `structures/trench/trench_field_cross_section.py`,
  `structures/via/pillar_field_cross_section.py`,
  `postprocessing/plot_field_cross_section.py`. Both examples run
  end-to-end with `R+T=1.0000`; output PNGs visually confirmed physically
  sensible.
- **Session summary**: all eight Category 9 / Phase 7 targets (9.1-9.8)
  complete. 512 tests pass project-wide (502 at the start of this
  sub-session: 493 fast + 9 slow -- 503 fast + 9 slow now, 10 new fast
  tests). No existing test weakened. Updated
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `decisions.md` (new ADR-015), `CONVENTIONS.md`, `phases.md`, `PRD.md`,
  `references.md`, `architecture.md`, `troubleshooting.md`, and this file.
  Category 6 (uncommitted from the prior sub-session) and this category's
  work committed and pushed together.

## 2026-08-05 (Category 7)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7 (Layer
  handling, targets 7.1-7.6): layer construction-time validation, a
  repeated-layer-identity regression guard, a Toeplitz-matrix cache
  (design then implementation), and layer-wise absorption (design then
  implementation).
- A real tension between the atomic-targets register (which lists a layer
  cache as a target) and `rules.md`'s Performance Requirements (which
  forbid caching before Phase 9 unless measurably too slow) -- resolved
  the same way Category 3's FFF/NVM tension was: measure first, then
  implement only what the measurement justifies, and record the gate
  explicitly in `decisions.md`.
- A measurement mistake made and caught mid-session: a first timing
  experiment (repeated identical patterned layers in one `solve()` call)
  wrongly attributed the "extra time per repeated layer" entirely to
  Toeplitz-matrix reconstruction; isolating the two costs directly showed
  the (out-of-scope) eigensolve actually dominates at high `num_orders`.
  Corrected before it was written into `design.md`/`decisions.md` as a
  justification -- the real, measured beneficiary is a fixed-wavelength
  angle sweep, not repeated layers within one call.
- Whether layer-wise absorption needs a new physics formula (a volumetric
  `Im(eps)*|E|^2` integral) or can reuse already-validated Category 9/
  Phase 7 pieces (`interior_amplitudes`, `propagate_amplitudes`,
  `z_poynting_flux`) -- chose the latter, per `rules.md` AI Coding Rule
  1's preference for composing validated blocks over deriving new formulas
  (same treatment ADR-015 gave interior-amplitude recovery).
- A genuine numerical-stability limitation found while validating the
  absorption energy-balance identity: a thick, highly lossy, high-
  `num_orders` layer numerically overflows the deepest evanescent modes'
  backward-propagated amplitude (`propagate_amplitudes`'s
  `exp(-i*q*z)` term), giving a nonsensical absorbed-power value. Not a
  formula bug -- the same instability class transfer-matrix methods are
  known for, now documented rather than silently avoided.

### Action items
- [x] Target 7.1 (layer validation audit) — done 2026-08-05.
  `layer._require_valid_thickness`. `tests/test_layer_validation.py`
  (15 tests). See `memory.md`.
- [x] Target 7.2 (repeated-layer identity) — done 2026-08-05.
  `tests/test_layer_repetition.py` (7 tests).
- [x] Target 7.3 (layer cache design) — done 2026-08-05. `design.md`'s
  "Layer/Toeplitz Caching Design" section, gated on a measured (and once
  corrected) timing case.
- [x] Target 7.4 (layer cache implementation) — done 2026-08-05.
  `Simulation._cached_toeplitz`/`_cached_toeplitz_component`,
  `decisions.md` ADR-016. `tests/test_layer_cache.py` (4 tests).
- [x] Target 7.5 (layer-wise absorption design) — done 2026-08-05.
  `design.md`'s "Layer-Wise Absorption Design" section, `decisions.md`
  ADR-017.
- [x] Target 7.6 (layer-wise absorption implementation) — done 2026-08-05.
  `SimulationResult.layer_absorption()`. `tests/test_layer_absorption.py`
  (4 tests), including a regression guard on the found numerical-overflow
  limitation. `troubleshooting.md` updated.
- **Session summary**: all six Category 7 targets (7.1-7.6) complete. 542
  tests pass project-wide (512 at the start of this category: 503 fast +
  9 slow -- 533 fast + 9 slow now, 30 new fast tests). No existing test
  weakened; one stale docstring claim (`test_stress_regression.py`, "layer-
  wise absorption isn't implemented yet") corrected without touching that
  file's actual assertions. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`,
  `memory.md`, `tasks.md`, `decisions.md` (new ADR-016, ADR-017),
  `design.md`, `references.md`, `troubleshooting.md`, and this file.

## 2026-08-05 (Category 8)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 8 (Solver
  sweeps and convergence, targets 8.1-8.8): a typed sweep-result
  container, wavelength/angle/polarization/thickness sweep APIs, a
  harmonic-order study, a conservative convergence criterion, and
  automatic harmonic-order selection built on top of it.
- Why `harmonic_study` needs a `Simulation`-*builder* callable rather than
  a single instance: `num_orders` is a construction-time parameter, and
  Category 7's Toeplitz-cache design (ADR-016) already established that
  `g` (and therefore `num_orders`) must stay fixed for a `Simulation`
  instance's whole lifetime -- resweeping `num_orders` on one live
  instance isn't supported.
- How to define a "conservative" convergence-stopping criterion (target
  8.7) given this project's own already-recorded evidence
  (`tests/test_fourier_convergence.py`) that high-contrast patterns can
  show a sharply non-monotonic low-order wobble -- settled on "every
  later point must also stay within tolerance," not just the immediate
  next one.
- A bug found and fixed by the project's own test-first discipline: the
  first version of the convergence criterion let the very last data point
  count as trivially converged (vacuously true against zero remaining
  points) -- caught immediately by a test using a never-actually-
  converging monotonic sequence, fixed before being trusted by
  `auto_select_num_orders`.

### Action items
- [x] Target 8.1 (result-series container) — done 2026-08-05.
  `sweep.SweepResult`.
- [x] Target 8.2 (wavelength sweep API) — done 2026-08-05.
  `sweep.sweep_wavelength`.
- [x] Target 8.3 (angle sweep API) — done 2026-08-05.
  `sweep.sweep_theta`/`sweep_phi`, confirmed to reuse the Category 7
  Toeplitz cache across a whole angle sweep.
- [x] Target 8.4 (polarization sweep API) — done 2026-08-05.
  `sweep.sweep_polarization`.
- [x] Target 8.5 (thickness sweep API) — done 2026-08-05.
  `sweep.sweep_thickness`, with explicit validation and thickness
  restoration.
- [x] Target 8.6 (harmonic-study API) — done 2026-08-05.
  `sweep.harmonic_study`.
- [x] Target 8.7 (convergence criterion) — done 2026-08-05.
  `sweep.find_convergence_index`, `decisions.md` ADR-018. Validated
  against thin-film/trench/pillar fixtures.
- [x] Target 8.8 (automatic harmonic selection) — done 2026-08-05, after
  8.7's validation passed. `sweep.auto_select_num_orders`.
- **Session summary**: all eight Category 8 targets (8.1-8.8) complete.
  569 tests pass project-wide (542 at the start of this category: 533
  fast + 9 slow -- 559 fast + 10 slow now, 26 new fast tests + 1 new slow
  test). No existing test weakened. Updated
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `decisions.md` (new ADR-018), `design.md`, `architecture.md`, and this
  file.

## 2026-08-05 (Category 10)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 10 (Optical
  outputs, targets 10.1-10.6): complex per-order field coefficients,
  diffraction angles, a one-call conservation report, layer-wise-
  absorption design cross-reference, per-order s/p conversion (gated on
  external validation), and a frozen output schema.
- Whether raw Cartesian `(Ex, Ey)` per-order output (target 10.1) needs
  per-order masking the way `diffraction_efficiencies()` does -- it
  doesn't, since `fields.tangential_e_field` is linear in modal
  amplitudes (unlike the bilinear Poynting-flux formula), so it can be
  evaluated once on the full amplitude vectors.
- Validating `complex_amplitudes()` against `tests/oracles/fresnel.py`
  surfaced a genuine finding: a naively hand-written textbook `r_p`
  Fresnel formula disagrees in sign with both this solver and the oracle
  (which agree with each other exactly) -- a real, pre-existing
  p-polarization sign-convention ambiguity, not a bug.
- Whether target 10.5's long-standing "not yet matched to S4/EMpy"
  polarization-convention gap could finally be closed this session --
  attempted a bounded external-validation check by reading S4's actual
  excitation-construction C++ source in full. Found S4 itself has an
  internal comment/code inconsistency there. A plausible derivation-level
  match to this project's convention was found, but S4 isn't buildable in
  this environment for a live numeric confirmation, so target 10.5
  remains explicitly deferred rather than claimed "validated" on
  derivation alone.

### Action items
- [x] Target 10.1 (complex coefficients) — done 2026-08-05.
  `SimulationResult.complex_amplitudes()`.
- [x] Target 10.2 (diffraction angles) — done 2026-08-05.
  `SimulationResult.diffraction_angles()`, new `kx`/`ky` fields.
- [x] Target 10.3 (conservation report) — done 2026-08-05.
  `SimulationResult.energy_balance()`.
- [x] Target 10.4 (loss accounting design) — already satisfied by
  Category 7 ADR-017, cross-referenced.
- [ ] Target 10.5 (polarization conversion) — evaluated and explicitly
  deferred, 2026-08-05. A bounded external-validation attempt against
  S4's actual source found a plausible but unconfirmed match; see
  `references.md`.
- [x] Target 10.6 (output schema tests) — done 2026-08-05.
  `tests/test_optical_outputs.py`'s frozen schema tests.
- **Session summary**: five of six Category 10 targets (10.1-10.4, 10.6)
  complete; target 10.5 explicitly deferred (not a gap, a documented
  decision). 570 tests pass project-wide (559 at the start of this
  category, 11 new fast tests, no new `slow` tests). No existing test
  weakened. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`,
  `tasks.md`, `design.md`, `references.md`, `CONVENTIONS.md`, and this
  file.

## 2026-08-05 (Category 11)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 11 (Semiconductor
  OCD features, targets 11.1-11.8): a validated CD-first OCD parameter
  object, a trapezoid trench constructor, corner-rounding geometry,
  TSV/grating example templates, and overlay/LER feasibility decisions.
- How to add corner rounding without a new geometry primitive or Fourier
  formula -- settled on arc-sampled `Polygon` vertices, reusing Category
  4's already-validated analytic `Polygon` Fourier transform unchanged.
- Whether overlay (layer-to-layer misregistration) needs new API --
  investigated directly (not assumed) and found it's already fully
  achievable via existing `Shape.center` offsets across layers sharing one
  `Lattice`, verified via a shift-by-one-period periodicity self-
  consistency check (same class of finding as Category 6 target 6.6's
  bottom-illumination result).
- Whether LER/LWR (stochastic edge roughness) can be approximated
  deterministically within RCWA's periodic-Fourier formulation -- concluded
  no reasonable deterministic proxy avoids overclaiming what it actually
  models; explicitly deferred rather than implemented.

### Action items
- [x] Target 11.1 (OCD parameter object) — done 2026-08-05.
  `ocd.OCDTrapezoidParams`.
- [x] Target 11.2 (trapezoid constructor) — done 2026-08-05.
  `ocd.trapezoid_trench_layers`.
- [x] Target 11.3 (corner-rounding design) — done 2026-08-05. Arc-sampled
  `Polygon`, `num_arc_points` as the convergence parameter.
- [x] Target 11.4 (corner-rounding implementation) — done 2026-08-05.
  `ocd.rounded_rectangle_polygon`.
- [x] Target 11.5 (TSV template) — done 2026-08-05.
  `structures/via/tsv_ocd_sweep.py`.
- [x] Target 11.6 (grating template) — done 2026-08-05.
  `structures/trench/trench_ocd_sweep.py`.
- [x] Target 11.7 (overlay feasibility) — done 2026-08-05, better answer
  than "define a model": already achievable, no new API. `decisions.md`
  ADR-019.
- [ ] Target 11.8 (LER/LWR feasibility) — evaluated and explicitly
  deferred, 2026-08-05. `decisions.md` ADR-020.
- **Session summary**: seven of eight Category 11 targets (11.1-11.7)
  complete; target 11.8 explicitly deferred (a documented decision, not a
  gap). 600 tests pass project-wide (570 at the start of this category,
  30 new fast tests, no new `slow` tests). No existing test weakened.
  Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `decisions.md` (new ADR-019, ADR-020), `architecture.md`, `references.md`,
  `structures/README.md`, and this file.

## 2026-08-05 (Category 12)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 12 (Linear
  algebra, targets 12.1-12.5): a baseline timing profiler, a direct-
  inverse audit, a factorization-reuse design note, an opt-in SVD
  diagnostic, and a sparse/iterative-methods feasibility decision.
- Where diagnostic/profiling scripts belong given they're neither
  `structures/` physics runs nor deterministic `pytest` assertions
  (wall-clock timing is machine-dependent) -- created a new top-level
  `profiling/` directory, explicitly documented as never asserted against
  a hard time/memory limit in any test.
- Whether any explicit matrix inverse in the codebase is "demonstrably
  unnecessary" (target 12.2's own gating language) -- audited every call
  site and found none were, but found a real, fixable house-convention
  inconsistency (three `eigenmodes.py` sites using `np.linalg.solve(A,
  eye(n))` instead of the project's documented `scipy.linalg.lu_factor`/
  `lu_solve` convention).
- Whether sparse/iterative linear algebra would help this project's
  eigenvalue problems -- measured the actual Toeplitz coupling matrix
  density directly rather than assuming, finding it 100% dense, closing
  the question structurally rather than leaving it open for later.

### Action items
- [x] Target 12.1 (baseline profiler) — done 2026-08-05.
  `profiling/baseline_profile.py`.
- [x] Target 12.2 (direct-inverse audit) — done 2026-08-05.
  `eigenmodes._dense_inverse`, three call sites migrated.
- [x] Target 12.3 (factorization reuse design) — done 2026-08-05.
  Documented in `design.md`; no further S-matrix-level reuse found beyond
  the already-shipped trivial-interface fast path.
- [x] Target 12.4 (SVD diagnostic) — done 2026-08-05.
  `eigenmodes.svd_diagnostics`.
- [x] Target 12.5 (sparse feasibility decision) — done 2026-08-05,
  evaluated and rejected (not deferred) on measured structural grounds.
  `decisions.md` ADR-021.
- **Session summary**: all five Category 12 targets complete. 612 tests
  pass project-wide (600 at the start of this category, 12 new fast
  tests, no new `slow` tests). No existing test weakened; one house-
  convention fix confirmed bit-for-bit equivalent before trusting it.
  Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `design.md`, `decisions.md` (new ADR-021), `architecture.md`,
  `references.md`, and this file.

## 2026-08-05 (Category 13)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 13 (Performance
  optimization, targets 13.1-13.6): a repeatable benchmark suite, an
  eigenmode-reuse cache (implementing Category 12 target 12.3's deferred
  design), a bounded vectorized wavelength sweep, and measured
  parallelism/GPU decisions.
- Whether the "eigenmode-cache" opportunity Category 12 flagged as a
  design note was actually worth implementing -- measured a genuine
  ~3.3x speedup on a polarization sweep before committing to it, matching
  the same "measure before optimizing" discipline `rules.md`'s
  Performance Requirements demands.
- How to safely vectorize "a single simple sweep" (target 13.4's own
  wording) without touching the general dense-eigensolve path -- scoped
  to uniform-isotropic-only (thin-film) stacks specifically, since that
  eigensolve is already closed-form and trivially batchable, unlike the
  general patterned/anisotropic case.
- Whether process- or thread-based parallelism would give a safe sweep
  speedup -- measured both directly rather than assuming either would
  help; found threading modestly helps (GIL released during LAPACK
  calls) while multiprocessing is actually counterproductive on this
  machine, likely from oversubscribing NumPy's already-multithreaded
  BLAS backend.
- Target 13.6 (GPU decision checkpoint) explicitly requires seeking the
  project owner's approval before any GPU work, once CPU targets are
  met -- asked directly; approval was not granted, so GPU/autodiff
  backend work stays deferred to Phase 9.

### Action items
- [x] Target 13.1 (benchmark suite) — done 2026-08-05.
  `profiling/benchmark_suite.py`.
- [x] Target 13.2 (Fourier-matrix reuse) — already satisfied by Category
  7's Toeplitz cache, cross-referenced.
- [x] Target 13.3 (eigenmode reuse) — done 2026-08-05.
  `Simulation._eigenmode_cache`, `decisions.md` ADR-022.
- [x] Target 13.4 (safe vectorized sweep) — done 2026-08-05.
  `vectorized.sweep_wavelength_vectorized`, `decisions.md` ADR-023.
- [x] Target 13.5 (parallelism decision) — done 2026-08-05, measured and
  documented, not implemented. `decisions.md` ADR-024.
- [x] Target 13.6 (GPU decision checkpoint) — done 2026-08-05, explicit
  approval sought and not granted; GPU work stays deferred to Phase 9.
- **Session summary**: all six Category 13 targets resolved (three
  implemented and tested, one already satisfied, two measured/decided
  without implementation). 627 tests pass project-wide (612 at the start
  of this category, 15 new fast tests, no new `slow` tests). No existing
  test weakened; one real bug (a missing `omega^2*I` term in the
  vectorized eigenmode helper) found and fixed by the equivalence test
  itself before being trusted. Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`,
  `memory.md`, `tasks.md`, `decisions.md` (new ADR-022, ADR-023,
  ADR-024), `architecture.md`, `references.md`, and this file.

## 2026-08-07

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Categories 14
  (Validation) and 15 (User interface and API) in one session, per the
  project owner's request ("complete phase 14 and 15 and push it").
- Category 14's central design question: what electromagnetic reciprocity
  actually predicts for a reversed multilayer stack at oblique incidence
  -- a naive same-`theta` comparison was tried first, found to fail badly
  (up to total mismatch at 45 deg), then the correct Snell's-law-matched
  comparison (constant transverse `kx`) was derived and verified to
  `~1e-15/1e-16` instead. A second question -- does this extend to
  patterned/diffractive layers -- was checked directly (not assumed) and
  found **not** to hold in the simple total-T sense.
- Category 14's harmonic convergence matrix: rather than guessing
  candidate `num_orders` lists and tolerances, each of the 7 geometry
  families was measured directly first. Found a new instance of the
  already-documented low-order non-monotonic wobble (this time in a
  "moderate contrast" 2D pillar, not just the previously-documented
  high-contrast case), and a genuinely slow (not `1e-2`-meeting)
  convergence rate for the tapered-via fixture, honestly reflected in a
  looser (`2e-2`) tolerance rather than tuned away.
- Category 15's configuration schema design: whether to invent a new
  material/pattern JSON sub-schema or reuse `geometry_io.py`'s existing
  one -- reused it unchanged (same material-dict shape,
  `pattern_from_dict` called directly for patterned layers), per
  `rules.md`'s "don't invent when something adequate already exists."
- Category 15's CLI exit-code design: whether "bad config" and "solver
  failure" should share one non-zero exit code or two distinct ones --
  chose two (`2` vs `1`) so a caller can tell them apart programmatically
  without parsing stderr text.
- Category 15's NumPy-export security question: how to serialize
  per-sweep metadata (a dict) into a `.npz` archive without needing
  `allow_pickle=True` on load -- JSON-encoded the metadata into a plain
  string array instead of a pickled object array, keeping the export path
  free of the untrusted-deserialization risk class `rules.md`'s Security
  Rules already flag for `eval`/`exec`/`pickle`.
- Category 15's HDF5 decision (target 15.8): evaluated against this
  project's actual result shapes (small, flat per-sweep arrays) and found
  no genuine advantage over the already-implemented `.npz` export yet --
  deferred, not implemented, per the same "evaluate before deciding"
  discipline as ADR-006/007/021/024.
- A user question mid-session ("Phase 11 and 12 is not written properly")
  was investigated in full (both sections re-read end to end) and found
  to have no actual defect; the user confirmed via a direct follow-up
  ("Category 11 and 12 look fine, keep going with 13") that no fix was
  actually needed -- no speculative edits were made to either section as
  a result.

### Action items
- [x] Target 14.1 (validation inventory) — done 2026-08-07.
  `testing.md`'s "Validation Inventory" section.
- [x] Target 14.2 (external 2D R/T oracle) — done 2026-08-07, re-evaluated
  and still blocked (S4 unbuildable in this environment, no versioned
  published dataset found); documented, not silently left unexamined.
- [x] Targets 14.3/14.4 (moderate/high-contrast 2D R/T tests) — done
  2026-08-07, documented as blocked on 14.2 rather than a false pass.
- [x] Target 14.5 (reciprocity test design) — done 2026-08-07,
  `decisions.md` ADR-025.
- [x] Target 14.6 (reciprocity tests) — done 2026-08-07,
  `tests/test_reciprocity.py` (11 tests).
- [x] Target 14.7 (harmonic convergence matrix) — done 2026-08-07,
  `tests/test_harmonic_convergence_matrix.py` (7 tests, all 4 `slow`-
  marked tests confirmed passing, 447.5s).
- [x] Target 14.8 (validation report) — done 2026-08-07, `testing.md`'s
  "Validation Report" section.
- [x] Target 15.1 (public API inventory) — done 2026-08-07, `design.md`'s
  "Public API Inventory" section; fixed a real staleness bug in
  `src/sougata_solver/__init__.py` (missing `Lattice1D`/`Ellipse`/
  `Polygon`/`Slab` exports) found while compiling it.
- [x] Target 15.2 (configuration schema) — done 2026-08-07, `config.py`.
- [x] Target 15.3 (configuration validation) — done 2026-08-07,
  `tests/test_config.py` (19 tests including malformed-input coverage).
- [x] Target 15.4 (configuration runner) — done 2026-08-07, reproduces
  `structures/thin_film/anti_reflection_coating.py` to `1e-12`.
- [x] Target 15.5 (CLI design) — done 2026-08-07, `cli.py`'s module
  docstring.
- [x] Target 15.6 (CLI implementation) — done 2026-08-07, `cli.py`,
  `sougata-solver` console-script entry, `tests/test_cli.py` (5 tests).
- [x] Target 15.7 (NumPy export) — done 2026-08-07, `export.py`,
  `tests/test_export.py` (4 tests).
- [x] Target 15.8 (HDF5 decision) — done 2026-08-07, evaluated and
  deferred, `decisions.md` ADR-026.
- **Session summary**: all 8 Category 14 targets and all 8 Category 15
  targets resolved. New modules `config.py`, `cli.py`, `export.py`; new
  test files `tests/test_reciprocity.py`, `tests/test_harmonic_convergence_matrix.py`,
  `tests/test_config.py`, `tests/test_cli.py`, `tests/test_export.py`.
  New `decisions.md` ADR-025 (reciprocity scope) and ADR-026 (HDF5
  deferred). `pyproject.toml` gained a `sougata-solver` console-script
  entry point. No existing test weakened. Updated
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`, `tasks.md`,
  `design.md`, `decisions.md`, `architecture.md`, `references.md`, and
  this file.

## 2026-08-07 (Categories 16-17)

### Discussed
- Completing `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Categories 16
  (Visualization) and 17 (Testing and quality) in the same explanatory
  style used for Categories 1-15, per the project owner's request
  ("continue with category 16,17, explain them like other category and
  then push and commit").
- Category 16's plot-data-contract question (target 16.1): whether
  plotting functions should accept a `Simulation` directly for
  convenience, or only plain arrays/already-computed result objects --
  chose the latter, matching `decisions.md` ADR-009/010's existing
  `structures/`-vs-`postprocessing/` split (plotting never triggers a
  solve) applied at the library-function level, and pinned it with a
  direct structural test rather than only a docstring claim.
- Category 16's unit-cell rendering approach (target 16.2): whether to
  special-case matplotlib patches per shape class, or reuse each
  `Shape.contains(x, y)` (already implemented by every shape) to
  rasterize a preview grid -- chose the latter for a single
  implementation covering every shape type uniformly, and confirmed it
  respects `Pattern`'s own "later shapes take precedence" rule.
- Category 17's test-taxonomy question (target 17.1): whether to add a
  pytest marker per testing tier (unit/integration/oracle/regression) to
  every one of 54 pre-existing test files, or document the
  already-largely-consistent existing filename+docstring convention --
  chose a narrower, precise addition instead: one new marker (`oracle`)
  applied only to files meeting a mechanically-checkable criterion
  (imports `tests/oracles/`), confirmed by grep across all 54 files
  rather than a per-file subjective judgment call. `decisions.md`
  ADR-027 records why the other, heavily-overlapping tiers were left to
  the existing filename/docstring/Validation-Inventory convention.
- Category 17's performance-regression-guard design (target 17.6):
  whether to assert an absolute wall-clock threshold or a relative,
  same-run ratio -- `rules.md`'s Performance Requirements explicitly
  rule out the former (machine-dependent); chose a same-run ratio
  (`eigensolve_time(81)/eigensolve_time(9)`, bounded with ~6x headroom
  above Category 12's measured ~160x baseline for the same fixture),
  `decisions.md` ADR-028.
- Category 17's regression-fixture provenance question (target 17.4):
  whether a frozen snapshot comparison could be described as a fresh
  oracle validation -- explicitly not: documented it as a regression
  guard for an already-independently-oracle-validated code path
  (`test_analytic_fresnel.py`/`test_thin_film_empy_cross_check.py`
  already validate the same uniform-multilayer solve path), not a new
  physics claim.
- Running `ruff` (target 17.5) surfaced two genuine dead-code findings
  in already-shipped code, not just test-file import cruft: unused local
  variables `modes_inc`/`modes_trans` in
  `src/sougata_solver/vectorized.py` (Category 13's vectorized sweep) and
  `epsilon_hat` in `tests/test_field_reconstruction.py` -- both verified
  genuinely unused via a direct grep for other references before
  removing, not assumed safe to delete.

### Action items
- [x] Target 16.1 (plot data contract) — done 2026-08-07, `plotting.py`.
- [x] Target 16.2 (geometry plot) — done 2026-08-07, `plot_unit_cell`/
  `plot_layer_stack`.
- [x] Target 16.3 (R/T spectrum plot) — done 2026-08-07, `plot_rt_spectrum`.
- [x] Target 16.4 (harmonic-convergence plot) — done 2026-08-07,
  `plot_harmonic_convergence`.
- [x] Target 16.5 (diffraction-order plot) — done 2026-08-07,
  `plot_diffraction_orders`.
- [x] Target 16.6 (field-intensity plot) — done 2026-08-07,
  `plot_field_intensity`.
- [x] Target 16.7 (Poynting/phase plots) — done 2026-08-07,
  `plot_field_phase`/`plot_poynting_vector`.
- [x] Target 17.1 (test taxonomy) — done 2026-08-07, `testing.md`'s "Test
  Taxonomy" section, `oracle` marker on 8 files, `decisions.md` ADR-027.
- [x] Target 17.2 (Windows CI) — done 2026-08-07, `.github/workflows/ci.yml`.
- [x] Target 17.3 (slow-test CI policy) — done 2026-08-07,
  `.github/workflows/slow-tests.yml`.
- [x] Target 17.4 (regression fixtures) — done 2026-08-07,
  `tests/regression_fixtures/thin_film_ar_coating_reference.npz`,
  `tests/test_regression_fixtures.py`.
- [x] Target 17.5 (static-analysis setup) — done 2026-08-07, `ruff`
  configured, 24 baseline issues found and fixed.
- [x] Target 17.6 (performance regression guard) — done 2026-08-07,
  `tests/test_performance_regression.py`, `decisions.md` ADR-028.
- **Session summary**: all 7 Category 16 targets and all 6 Category 17
  targets resolved. New module `src/sougata_solver/plotting.py`; new
  test files `tests/test_plotting.py`, `tests/test_regression_fixtures.py`,
  `tests/test_performance_regression.py`; new directory
  `tests/regression_fixtures/`; new CI workflows
  `.github/workflows/ci.yml`/`slow-tests.yml`; `ruff` added and
  configured, 24 baseline issues fixed (2 genuine dead-code findings, 22
  unused imports). New `decisions.md` ADR-027 (test taxonomy) and
  ADR-028 (performance regression guard design). No existing test
  weakened. 702 tests collected project-wide (683 before this session's
  work started). Updated `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, `memory.md`,
  `tasks.md`, `design.md`, `decisions.md`, `architecture.md`,
  `references.md`, `README.md`, and this file.

## 2026-08-17 (RCWA postprocessing overlay)

### Discussed
- The project owner is no longer comparing against the KLA reflectance
  calculator; comparison oracle data now comes from an external
  RCWA_module tool, exported as `OUTPUT_RCWA/**/*.txt` files with a
  `lambda(m), Y` header (comma-delimited, wavelength in **meters**,
  reflectance column literally named `Y`) rather than KLA's
  `"Wavelength (nm)"`/`"Reflectance"` (tab-delimited, nm) format.
- `postprocessing/RCWA_plot_norm.py` (new script, overlays one RCWA
  export against one solver `output_*_RT.csv`) was raising
  `ValueError: Could not find wavelength/R columns ... got fields
  ('lambdam', 'Y')` because its column-detection helper only recognized
  headers starting with `wavelength`/`r`. Even with detection fixed, the
  meters-vs-nanometers unit mismatch (RCWA: `8e-07` = 800 nm; solver
  CSV: `4.000000e+02` = 400 nm) would have plotted the two curves on
  incompatible x-axis scales.
- Fix: extended the column matcher to also accept `lambda*`/`y` as
  wavelength/R aliases, and added an auto meters->nm conversion
  triggered when the matched wavelength field starts with `lambda` and
  its max value is under 1e-3 (unambiguously meters, never a real nm
  spectrum). Applied identically to `postprocessing/RCWA_plot_norm.py`
  and to `postprocessing/plot_rcwa_reflectance.py` (renamed from
  `plot_kla_reflectance.py`, its single-file KLA plotting counterpart,
  which is now retired since KLA is no longer the comparison source).
- Verified by running both scripts end to end against
  `OUTPUT_RCWA/Thin_Film/17_08_26/Multi/0_degree.txt`: both parse the
  file and save a plot without error. The two example runs used for
  `RCWA_plot_norm.py` (`0_degree.txt` vs.
  `outputs/2026_08_06/12_50_00_.../output_multistack_RT.csv`) have
  non-matching wavelength grids, so the script correctly falls back to
  its documented "visual guide only" path rather than printing a
  per-point max-diff.
- No `rules.md`/`testing.md` doc specifies a required RCWA-export schema
  or wavelength-unit convention for this comparison; the meters-vs-nm
  handling above is inferred directly from the two files' contents, not
  cited from a doc.

### Action items
- [x] Fix `RCWA_plot_norm.py` column detection + unit conversion for
  RCWA_module's `lambda(m), Y` export format — done 2026-08-17.
- [x] Retire KLA-specific `plot_kla_reflectance.py`, replace with
  `plot_rcwa_reflectance.py` (RCWA_module format, same column-detection
  and meters->nm fix) — done 2026-08-17.
- [ ] If RCWA_module ever exports a numeric-comparison-ready wavelength
  grid matching the solver's own sweep points, revisit
  `RCWA_plot_norm.py`'s "visual guide only" fallback to also print a
  per-point max-diff for that case.

## 2026-08-18 (ADR-033: linear-polarization `alpha` convention flip)

### Discussed
- The project owner supplied the commercial RCWA tool's (Lumerical FDTD)
  own polarization-mixing script for a "grating_power" export:
  `R_linear = sin(alpha)^2 * Rs_power + cos(alpha)^2 * Rp_power`, with an
  explicit in-script comment `(0=P, 90=S)` — the opposite reference axis
  from this project's pre-existing `s=cos(alpha), p=sin(alpha)*exp(i*delta)`
  convention (`0=S, 90=P`).
- Comparing the solver's `linear_15deg`/`linear_30deg` states (45 deg
  incidence, `sio2_sio_ni_sio2_on_semi_infinite_si` stack) against that
  tool's `Linear15_Linear30.txt` export showed both a large apparent
  magnitude gap (peak `R` 0.53 vs 0.35) and a reversed 15-vs-30 ordering.
  Back-solved the tool's raw `Rss`/`Rpp` from its two exported curves (a
  per-wavelength 2x2 linear solve against its own stated mixing formula)
  and compared directly to the solver's own pure-TE/pure-TM `R` at the same
  angle — matched to ~0.1% absolute (max diff 0.0013, RMS 0.0005,
  noise-level from the tool's coarser wavelength grid), proving the entire
  discrepancy was this labeling mismatch, not a solver or oracle physics
  error.
- Flipped the convention (`s=sin(alpha), p=cos(alpha)*exp(i*delta)`) in
  `CONVENTIONS.md`'s worked-examples table,
  `structures/thin_film/custom_multistack.py::_jones_state` (TE/TM
  `alpha_deg` entries swapped to keep producing the same physical states),
  and `structures/thin_film/sio2_on_si_thin_film.py::_polarization_amplitudes`.
  Grepped every `structures/thin_film/*.py`, `structures/trench/*.py`,
  `structures/via/*.py`, and `postprocessing/*.py` file for this
  alpha-to-amplitude formula first — confirmed only these two files
  implement it, so no other structure script needed the flip.
- Confirmed `tests/test_polarization_states.py` needed no change: it
  hardcodes numeric `(s,p)` pairs directly rather than deriving them from
  this formula, so its `"linear_20deg"` label is cosmetic there.
- Full reasoning, RCP/LCP unaffected-angle math, and the file-by-file scope
  check are recorded in `decisions.md` ADR-033.

### Action items
- [x] Flip `alpha` convention to match the commercial tool (0=P, 90=S) in
  `CONVENTIONS.md`, `custom_multistack.py`, `sio2_on_si_thin_film.py` —
  done 2026-08-18, `decisions.md` ADR-033.
- [x] Update `memory.md` with the ADR-033 summary — done 2026-08-18.
- [ ] Re-run every `structures/thin_film/*_ellipsometry_run.py` /
  `postprocessing/jones_mueller_ellipsometry.py` pair that consumes a
  named linear/elliptical polarization state, to confirm downstream
  Jones/Mueller/Psi-Delta outputs still make sense under the flipped
  convention (not yet re-verified this session — only the two thin-film
  R/T scripts above were touched).

## 2026-08-18 (composite-grating structure cross-validated against Lumerical RCWA)

### Discussed
- The project owner wanted `multistack_composite_grating.py`'s new
  laterally-alternating composite structure (Si/SiO2 on one half of a 2 um
  period, Ni/SiO on the other) cross-checked against the equivalent
  structure they built in Lumerical's RCWA solver, and asked which
  Lumerical result to export: settled on `grating_power`
  (`Rs_power`/`Ts_power`/`Rp_power`/`Tp_power`), not `grating_order`
  (per-order, too granular) or `grating_characterization` (setup/
  convergence info, not a spectrum).
- `grating_power`'s attributes turned out to be per-diffraction-order
  (`size(Rs_power) = 401 x 17 x 9`, indexed by wavelength x n-order x
  m-order), not already-summed totals like the earlier ADR-033 comparison
  assumed (that case only had one order, a plain thin film, so summing was
  never visibly necessary there). Fixed by summing over both order
  dimensions in the exported Lumerical script (`pinch()` to strip
  singleton `f`/`theta`/`phi` dims, then `sum(sum(Rs_p,3),2)`).
- First overlay attempt (`postprocessing/overlay_composite_grating_vs_lumerical.py`,
  new script) showed a large mismatch: `R` max\|diff\|=0.28, `T`
  max\|diff\|=0.69. Ruled out materials data first (this solver's
  `NK_FILE/*_KLA.txt` permittivity for Si/Ni/SiO2/SiO matched Lumerical's
  Palik-fit plots closely at every sampled wavelength) before concluding
  the real cause: Lumerical's `Si_substrate`/`Ni_substrate` objects were
  drawn 5 um deep while the RCWA region's z-extent only reached -0.5 um,
  so (per Lumerical's documented incidence/transmission-inference
  behavior) it was resolving a genuinely semi-infinite, laterally-patterned
  Si(left)/Ni(right) exit -- not "finite Si/Ni over air" like this
  solver's model. Confirmed with a thickness sweep (0.5/2/5/10 um) on this
  solver's side showing `R` still hadn't converged at 10 um for the
  800nm/weakly-absorbing case (real Fabry-Perot fringing from the buried
  interface), proving a thick-finite-layer workaround wouldn't have worked
  either.
- Project owner fixed the Lumerical model to match (`Si_substrate`/
  `Ni_substrate` made genuinely finite, 0.5 um, real air below; RCWA
  region `z min` nudged to -0.6 um). Re-running the overlay against the
  corrected export gave `R` max\|diff\|=0.013 (RMS 0.0045), `T`
  max\|diff\|=0.012 (RMS 0.0022) -- full reasoning in `decisions.md`
  ADR-034.
- Also fixed along the way: the Lumerical export script's `write()`
  appends rather than overwrites (three concatenated header+data blocks
  accumulated in one file across re-runs) -- the overlay script's loader
  now always uses the last block; and an absolute-path `write()` call
  failed on the machine actually running Lumerical (different computer,
  different drive layout) -- resolved by writing next to the `.fsp` file
  and copying the result over manually.

### Action items
- [x] Build `structures/thin_film/multistack_composite_grating.py` and
  cross-validate against Lumerical RCWA — done 2026-08-18, `decisions.md`
  ADR-034, ~1% R/T agreement.
- [x] `postprocessing/overlay_composite_grating_vs_lumerical.py` — done
  2026-08-18.
- [ ] The *original* Lumerical structure (Si/Ni as a genuinely
  semi-infinite, laterally-patterned substrate, no common material
  beneath) is not representable by this solver today (`Simulation.
  transmission` must be one uniform `Material`) — logged as a real
  capability gap in `decisions.md` ADR-034, not implemented.

## 2026-08-18 (trench/via Lumerical build guidance -- structure identity still open)

### Discussed
- Project owner asked for help building "the trench structure" from
  `structures/trench/` in Lumerical RCWA, to eventually cross-validate the
  same way `multistack_composite_grating.py` was (ADR-034). Walked through
  concrete Lumerical build steps for `trench_grating.py` (period 0.7 um,
  30% fill-factor Si ridge in air, constant `n=3.48` matching the
  Moharam/Gaylord 1995 oracle `tests/test_1d_grating.py` cross-checks --
  flagged that swapping in dispersive Palik-Si data would decouple the
  comparison from that oracle test per `rules.md` AI Coding Rule 3, so it
  should be a new script, not an edit to `trench_grating.py` itself, if
  wanted later).
- Also clarified (not yet built): only one Rectangle/Slab object is needed
  for a two-material patterned layer when one side is the RCWA background
  material (air) -- a second explicit object is only needed when *neither*
  side is the background, as in `multistack_composite_grating.py`.
- Three rounds of hand-drawn sketches from the project owner progressively
  changed the actual target structure's identity: (1) an XY sketch
  initially looked like a bounded island (finite in both x and y) rather
  than a y-invariant strip: a `Lattice1D` ridge/groove strip must touch the
  frame's top/bottom edges in an XY view, since it's invariant in y -- a
  gap there would mean an unintended 2D pattern; (2) a side/cross-section
  sketch then showed the roles inverted from `trench_grating.py`'s
  convention -- Air as the minority (trench) region cut into a Si
  majority background, not a Si ridge sitting in an air groove; (3) the
  actual top view showed the "Air" region bounded in *both* x and y (a
  finite rectangle, not a full-height strip) -- meaning the real target is
  a 2D array of rectangular air holes in a Si background (`structures/
  via/`-family, `Lattice`+`Rectangle`), not a 1D trench (`Lattice1D`) at
  all. Closest existing template identified: `structures/via/via_array.py`
  (circular air via in Si, Si semi-infinite below) -- would need only a
  `Circle`->`Rectangle` shape swap, but exact dimensions (period(s), hole
  width/height, thickness, constant-vs-dispersive Si) are still unknown --
  no numbers given yet, only hand sketches.
- Project owner will share the *actual* Lumerical trench structure
  tomorrow instead of continuing to iterate on sketches; work on this
  paused until then.

### Action items
- [ ] Build the real trench/via structure once the project owner shares
  their actual Lumerical model (not yet built -- today's sketches kept
  changing the identified geometry, so nothing should be assumed from them
  alone). Likely `structures/via/rectangular_via_array.py` based on the
  last sketch, but confirm against the real model first, including the
  incidence/transmission stack (air / Si-semi-infinite, matching
  `via_array.py`? not yet confirmed).
