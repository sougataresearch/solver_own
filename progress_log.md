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

### Action items
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
