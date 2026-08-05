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
