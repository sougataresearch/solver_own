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
