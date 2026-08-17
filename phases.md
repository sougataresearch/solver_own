# Roadmap — sougata_solver

This supersedes/formalizes the plan already approved and saved at
`C:\Users\d14k4\.claude\plans\vivid-swimming-moler.md`. Phase numbering
here is the authoritative one going forward; update this file (not the
plan-mode scratch file) as phases complete.

## Phase 1 — Uniform Multilayer Core — **DONE**

- **Objectives**: solve reflectance/transmittance for arbitrary stacks of
  uniform, isotropic, dispersive layers; support arbitrary incidence
  angle/polarization; report Jones/Mueller polarimetry.
- **Deliverables**: `materials.py`, `layer.py`, `eigenmodes.py` (uniform
  case), `smatrix.py`, `excitation.py`, `fields.py` (R/T), `polarimetry.py`,
  `simulation.py` (uniform path); Fresnel-oracle test suite; the
  `structures/` scripts covering the anti-reflection coating, SiO2-on-Si,
  custom-multistack, and custom-material-from-nk-data cases, plus the
  `postprocessing/` Jones/Mueller ellipsometry split.
- **Estimated complexity**: (retrospective) Medium — closed-form eigenmodes
  avoided the hardest numerical risk (general eigendecomposition), but the
  S-matrix sign/normalization conventions required careful,
  line-by-line-verified transcription from S4.
- **Dependencies**: none (foundation phase).
- **Status**: shipped, validated against analytic Fresnel/TMM
  (`tests/test_analytic_fresnel.py`).

## Phase 2 — Fourier-Factorization Core — **DONE**

- **Objectives**: build the dimension-agnostic infrastructure that turns a
  `Pattern` (shapes + background) into the Toeplitz permittivity matrices
  every patterned-layer eigensolver (Phase 3, 4) needs.
- **Deliverables**: `src/sougata_solver/fourier_factorization.py` with
  `pattern_epsilon_hat(...)` and `toeplitz_matrix(...)`, producing both
  `epsilon_hat` (direct) and `epsilon_inv_hat` (inverse-rule) Toeplitz
  matrices; unit tests (`tests/test_fourier_factorization.py`) comparing
  analytic Toeplitz entries against **two independent** numerical
  references for `Circle` and `Rectangle`: a from-scratch rasterize-and-sum
  (direct Riemann-sum evaluation of the Fourier integral) and a literal
  FFT-of-rasterized-mask reproducing the vendored `Rigorous-Coupled-Wave-Analysis`
  (Python `convmat2D.py`) / `RigorousCoupledWaveAnalysis.jl`
  (`ft2d.jl::real2recip`) algorithm; plus an anisotropic-material
  `NotImplementedError` guard; `.flake8`/`mypy.ini` added.
- **Estimated complexity**: (retrospective) Medium, as estimated. The
  formula was transcribed from
  `S4/S4/pattern/pattern.c::pattern_get_fourier_transform` (lines
  889-1029) for the per-shape subtraction-rule accumulation, and
  `S4/S4/fmm/fmm_closed.cpp::FMMGetEpsilon_ClosedForm` (lines 77-127) for
  how direct vs. inverse-rule Toeplitz entries are assembled from
  `G_i - G_j`. Two surprises worth recording: (1) a rectangle's sharp
  edges alias more than a circle's smooth boundary at a given grid
  resolution, so the rasterize-and-sum cross-check needed a finer grid
  (900x900) than first tried (300x300) before agreeing to 5e-3; (2) the
  first attempt at the FFT cross-check used an uncentered `[0, Lx)` raster
  grid to mirror `fft2`'s domain convention directly, which silently
  truncated any shape whose footprint crossed the domain edge (caught
  because the resulting DC term didn't match the pattern's true
  area-weighted-average permittivity) — fixed by rasterizing on the same
  centered grid as the other reference and applying `numpy.fft.ifftshift`
  before the FFT instead, verified against the DC term before trusting it
  for the full matrix.
- **Dependencies**: `geometry.py` (`Shape.fourier_transform`,
  `Pattern.containment_tree`) and `fourier_basis.py`
  (`truncate_fourier_orders`) — both already implemented, no changes
  needed to either.
- **Status**: shipped, validated against two independent numerical
  references (from-scratch rasterize-and-sum, and an FFT-of-rasterized-mask
  reproduction of the vendored `Rigorous-Coupled-Wave-Analysis`/
  `RigorousCoupledWaveAnalysis.jl` convolution-matrix algorithm) for both
  `Circle` and `Rectangle` patterns, direct and inverse-rule, at several
  nonzero G-vectors (`tests/test_fourier_factorization.py`, 12 tests).
  Scalar isotropic materials only — anisotropic materials raise
  `NotImplementedError` naming Phase 6, per the AI Coding Rules' scope
  discipline.

## Phase 3 — 1D-Periodic Lamellar Gratings (Trench) — **DONE**

- **Objectives**: solve reflectance/transmittance/diffraction efficiencies
  for a 1D-periodic patterned layer (line/space, i.e. a trench), as the
  first end-to-end patterned-layer capability — chosen before 2D because
  1D gratings decouple TE/TM into independent scalar eigenproblems, making
  this the lower-risk place to validate the Phase 2 Fourier-factorization
  machinery and the general non-uniform eigenmode-solve pattern.
- **Deliverables**: `Lattice1D`, a `Slab`/`Line` shape in `geometry.py`;
  `truncate_fourier_orders_1d` in `fourier_basis.py`;
  `solve_layer_eigenmodes_1d(...)` in `eigenmodes.py`; a `Lattice1D`
  dispatch branch in `simulation.py`; a validation test against a published
  1D binary-grating benchmark (Moharam & Gaylord 1995 or equivalent);
  `structures/trench/trench_grating.py`; the energy-conservation invariant
  (`R + T + sum(diffraction efficiencies) + A = 1`) and a measured
  convergence-rate-vs-`num_orders` check compared against Li (1996)'s
  predicted rate for the chosen Fourier-factorization rule — both defined
  in `testing.md`'s Physical-Invariant Testing, and required here as the
  first patterned-layer phase precisely because these checks are
  oracle-independent and should catch a wrong-but-benchmark-coincidental
  result before the benchmark comparison is even run.
- **Estimated complexity**: Medium-High. The physics (decoupled scalar
  TE/TM eigenproblems) is simpler than Phase 4's general case, but this
  phase is also where the *first* non-uniform eigensolver gets built and
  validated, so unexpected issues surfacing here are likely to be
  Fourier-factorization bugs from Phase 2, not 1D-specific ones — budget
  time accordingly.
- **Dependencies**: Phase 2.
- **Status**: shipped. Key finding during implementation: 1D gratings'
  "decoupled scalar TE/TM eigenproblem" turned out not to be a separate
  formula at all — reading `S4/S4/rcwa.cpp::SolveLayerEigensystem`
  (lines 684-827) in full showed S4 has no dedicated 1D/TE/TM code path;
  a 1D lattice is just the *general* non-uniform eigenoperator with the
  G-vector set restricted to `ky=0`, which then reduces algebraically to
  exactly block-diagonal (verified directly by a unit test on a random
  Toeplitz input, not just claimed). So `solve_layer_eigenmodes_1d` is a
  genuine specialization of Phase 4a's eventual general solver, not an
  independently-derived scalar formula — this is the concrete form of
  "validating the general non-uniform eigenmode-solve pattern" the
  objectives above named. A second finding, from reading
  `S4/S4/fmm/fmm_closed.cpp:77-127`'s "1D proper FFF rule" branch: the
  `TM`-like block's effective permittivity is `inv(epsilon_inv_hat)` (the
  matrix-inverse of Phase 2's *inverse-rule* Toeplitz), not the direct-rule
  Toeplitz — Li's (1996) factorization rule, and the reason Phase 2 built
  both `epsilon_hat` and `epsilon_inv_hat` in the first place. Validated
  against `tests/oracles/rcwa_1d_gaylord.py` (hand-transcribed from
  `Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples/1D_Grating_Gaylord_{TE,TM}.py`,
  which cites Moharam, Grann, Pommet & Gaylord 1995): TE agrees with the
  oracle to ~1e-10 at modest `num_orders`; TM converges to the same limit
  as the oracle but only at high `num_orders` (measured directly, a real
  convergence-rate asymmetry between the direct and inverse Fourier rules,
  consistent with Li 1996) — and the oracle's own TM source file
  self-reports "STILL NOT WORKING YET" in its module docstring, so TM's
  primary correctness evidence is the energy-conservation invariant (holds
  to `1e-8`+ across TE/TM/mixed polarization, normal and oblique incidence)
  and the reduces-to-Phase-1-uniform-result regression test, not the
  oracle comparison alone — see `tests/test_1d_grating.py` for the full,
  honestly-caveated account. 91 tests pass project-wide (75 pre-existing +
  16 new).

## Phase 4a — 2D-Periodic Patterned Layers, Well-Conditioned Case (Via, Pillar) — **DONE**

- **Objectives**: solve reflectance/transmittance/diffraction efficiencies
  for a full 2D-periodic patterned layer using the existing `Circle`/
  `Rectangle` shapes on a *moderate*, well-conditioned case (single circular
  via, moderate index contrast, moderate `num_orders`) — i.e. get the
  general eigensolver and its S4/benchmark cross-check working end to end,
  removing the `NotImplementedError` at `simulation.py:98`. Split out from
  the original single Phase 4 specifically so the highest-risk numerical
  work (near-degenerate eigenvalues, ill-conditioning) has its own phase
  with its own dedicated validation, instead of being one implicit sub-task
  inside a single large phase where a passing test on an easy case could
  mask a fragile implementation.
- **Deliverables**: `solve_layer_eigenmodes_patterned(...)` in
  `eigenmodes.py` (general non-uniform eigenproblem, transcribed from
  `S4/S4/rcwa.cpp::SolveLayerEigensystem`, lines 794-827); `simulation.py`
  wiring for the 2D patterned path; an S4-cross-check validation test (or,
  if S4 isn't buildable/runnable in this environment, an explicitly-flagged
  literature benchmark instead — never a fabricated "it matches" claim per
  `rules.md`'s AI coding rules); `structures/via/pillar_array.py`,
  `structures/via/via_array.py`; the energy-conservation check
  (`R + T + sum(diffraction efficiencies) + A = 1`, see `testing.md`'s
  Physical-Invariant Testing) applied to this phase's cases as a first,
  oracle-independent correctness signal, before the S4/benchmark cross-check.
- **Estimated complexity**: High — general complex eigendecomposition is
  still new machinery here, but the chosen test cases are deliberately
  scoped to avoid the near-degenerate/ill-conditioned regime (deferred to
  Phase 4b) so this phase's risk is "get the transcription right," not
  "get the transcription right *and* numerically stable everywhere."
- **Dependencies**: Phase 2; benefits from Phase 3 having already
  shaken out Fourier-factorization bugs on the simpler 1D case.
- **Status**: shipped, with a correction worth recording as a cautionary
  finding. A first implementation (by a different agent working on this
  phase) built `Epsilon2` for the 2D case by copying
  `solve_layer_eigenmodes_1d`'s `epsilon_hat`/`inv(epsilon_inv_hat)`
  block construction — plausible-looking, consistent with the 1D
  docstring's citation of `rcwa.cpp:794-827`, and it passed every test
  including a "ky=0 reduces to the 1D solver" check. It was wrong: that
  `Epsilon2` construction is only valid inside `fmm_closed.cpp`'s
  `0==Lr[2]&&Lr[3]==0` branch (the 1D case) — nobody had actually read the
  adjacent true-2D branch (`fmm_closed.cpp:133-139`) before writing or
  reviewing this code. Reading it (this session) shows S4's actual default
  2D closed-form behavior (no polarization basis) is plain Laurent's rule
  throughout: `Epsilon2 = block_diag(epsilon_hat, epsilon_hat)`, and even
  `kp` is built from `inv(epsilon_hat)` (a numerical matrix-inverse), not
  Phase 2's separately-factorized `epsilon_inv_hat` — that quantity turns
  out to be 1D-only in S4's own source, not general-2D infrastructure as
  first assumed. The "ky=0 reduces to 1D" test that should have caught
  this didn't, because it was circular: both solvers used the identical
  (wrong) formula, so of course they agreed — a concrete instance of
  `rules.md`'s warning that a passing test on an easy/coincidental case can
  mask a fragile implementation. Fixed and re-validated (energy
  conservation and both `structures/via/` scripts re-checked after the
  fix); the "ky=0" test was replaced with two honest ones — see
  `memory.md`'s Phase 4a entry and `eigenmodes.solve_layer_eigenmodes_patterned`'s
  docstring for the full citation.

  **Follow-up (same day, per explicit user request to stop relying on S4
  alone)**: surveyed all three other vendored RCWA-family repos (`EMpy`,
  `Rigorous-Coupled-Wave-Analysis`, `RigorousCoupledWaveAnalysis.jl`) for
  material relevant to Phase 4a. `EMpy/EMpy/RCWA.py` ruled out (1D-only, no
  2D support, an author-acknowledged instability hack). Confirmed (a third
  independent source, alongside S4 and this survey) that
  `Rigorous-Coupled-Wave-Analysis`'s `run_RCWA_2D` also uses plain
  Laurent's rule for 2D with no inverse-rule correction — reassuring that
  the fixed formula's rule choice isn't an S4 idiosyncrasy. Built a real
  eigenoperator oracle from `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:57-99`
  (`tests/oracles/rcwa_2djl_eigenvalues.py`) — a **structurally different**
  derivation (direct Maxwell-curl elimination, not S4's `Epsilon2 @ kp`
  route) that, after reconciling its `k0`-normalization and an overall
  sign-convention difference (both confirmed empirically, documented in
  the oracle's docstring), matches `solve_layer_eigenmodes_patterned`'s
  `q^2` eigenvalues to ~1e-12 across several `num_orders`/angle/pattern
  cases (`tests/test_2d_pillar.py::test_2d_patterned_eigenvalues_match_rcwa_jl_oracle`,
  6 parametrized cases). This closes the eigenoperator-correctness gap with
  a genuinely independent formula — the exact class of check that would
  have caught the original bug (unlike the old circular "ky=0 reduces to
  1D" test). A full external **R/T** oracle (not just eigenvalues) is still
  an open gap — no independently-published 2D benchmark was found in any
  vendored repo, S4 isn't buildable here (no `cmake`/Lua toolchain), and
  Julia isn't installed either (`which julia` fails) — explicitly flagged
  rather than faked (`tests/oracles/rcwa_2d_pillar.py`), carried forward as
  Phase 4b work alongside that phase's originally-scoped near-degenerate
  stress cases. 107 tests pass project-wide (101 prior + 6 new
  parametrized eigenoperator-oracle cases).

## Phase 4b — 2D-Periodic Patterned Layers, Near-Degenerate / Ill-Conditioned Cases — **DONE**

- **Objectives**: stress-test Phase 4a's general eigensolver on the cases
  that actually risk numerical instability — high index contrast, small
  feature-to-period ratio, high `num_orders` — and add explicit handling
  (or an explicit, tested failure boundary) for near-degenerate eigenvalues
  with poorly-conditioned eigenvectors. This is the concrete form of the
  risk `design.md` names as "the highest-risk remaining algorithm in the
  project" and `PRD.md`'s Risks table lists — Phase 4a alone does not
  retire that risk, since a single well-conditioned test case passing says
  nothing about behavior near degeneracy.
- **Deliverables**: a documented degenerate-eigenvalue handling strategy in
  `solve_layer_eigenmodes_patterned`'s docstring (reusing/extending the
  already-validated `_select_q_branch` outgoing-mode convention, per
  `PRD.md`'s Risks mitigation); a dedicated regression test suite
  exercising high-contrast/high-`num_orders` cases against S4 (or a
  published benchmark) and against the energy-conservation invariant; an
  explicit condition-number diagnostic (logged at `WARNING`, per
  `design.md`'s Logging Strategy) when the Toeplitz or eigenvector matrix
  is ill-conditioned, rather than silently returning a degraded answer.
- **Estimated complexity**: High — this is where a from-scratch
  implementation is most likely to diverge from S4's already-solved
  stability fixes; budget real time for iterating against the S4 cross-check
  rather than treating one pass as sufficient.
- **Dependencies**: Phase 4a (reuses its eigensolver and wiring; this phase
  only adds stress cases and stability handling, not new solve machinery).
- **Status**: shipped, with an honest (not fabricated) finding: a
  deliberate stress sweep — index contrast from `n=3.48` to a
  lossy-metal-like `-20+2j`, `num_orders` up to 225, near-touching circular
  pillars (`radius=0.49*period`), a sub-percent-halfwidth sliver rectangle,
  and near-degenerate nested circles (`1e-4`-scale radius difference) —
  did **not** find a catastrophic near-degenerate/ill-conditioned failure
  for the closed-form isotropic `Circle`/`Rectangle` patterns tested:
  `cond(epsilon_hat)` reached ~900 and `cond(phi)` reached ~170 in the
  worst cases, while energy conservation and the independent
  `RigorousCoupledWaveAnalysis.jl` eigenvalue oracle (built in Phase 4a)
  both held to ~1e-10 throughout (`tests/test_2d_pillar_stress.py`). This
  is reported as what was actually observed, not oversold as "no
  pathological case can exist" — `numpy.linalg.eig` could still misbehave
  on an input this sweep didn't probe. Given no failure to fix, the
  "handle near-degenerate eigenvalues" deliverable became **detection**:
  `cond(epsilon_hat)`/`cond(phi)` now logged at `WARNING`
  (`eigenmodes.ILL_CONDITIONED_THRESHOLD = 1e4`, ~10x headroom above the
  worst observed case, module-level `logger` per `design.md`'s Logging
  Strategy), verified by dedicated logging tests (`caplog`/`monkeypatch`
  to force the threshold and confirm the mechanism fires, plus a
  no-warning-in-the-ordinary-case check). The "match S4 or a published
  benchmark" deliverable was **not met as originally worded** — S4 needs
  `cmake`/a Lua toolchain (neither present) and Julia isn't installed
  either (`which julia` fails, both confirmed not assumed) — so stress
  cases are instead cross-checked against the same independent
  `RigorousCoupledWaveAnalysis.jl` eigenoperator oracle used for Phase 4a's
  base validation; `tasks.md` records this substitution honestly (one
  unchecked box, not silently marked done) rather than claiming a match
  that didn't happen, per `rules.md` AI Coding Rule 5. 118 tests pass
  project-wide (107 prior + 11 new: 9 stress-case oracle cross-checks + 2
  logging tests — see `tests/test_2d_pillar_stress.py` for the exact
  breakdown).

## Phase 5 — Tapered / Sloped Sidewalls (Via, Trench) — **DONE**

- **Objectives**: represent a via or trench with linearly tapered
  sidewalls via staircase (z-discretized) layer approximation, and
  demonstrate R/T convergence as slice count increases.
- **Deliverables**: a small staircase-layer-stack generator (given top/
  bottom feature size, thickness, and slice count `N`, produce `N` `Layer`s
  with linearly interpolated `Circle`/`Rectangle`/`Slab` sizes); a
  convergence-vs-`N` test/example for both a tapered via and a tapered
  trench (marked `slow` per the existing pytest marker).
- **Estimated complexity**: Low — no new Fourier/eigenmode math, purely a
  layer-stack generation helper consuming Phase 3/4a's already-validated
  per-layer solvers.
- **Dependencies**: Phase 3 and Phase 4a (needs at least one working
  patterned-layer eigensolver to stack; ideally both, to cover tapered
  trench and tapered via). Does not require Phase 4b — staircase generation
  doesn't push into the near-degenerate regime any more than the base
  patterned solver already does.
- **Status**: shipped. Per the `phase-reference-picker` skill's procedure,
  every RCWA-family repo under `REFERENCE/` was grepped for "stair"/"taper"
  before writing any code; none had a matching staircase-generator
  implementation to cite (only unrelated `meep`/`gprMax` docs hits), so
  `src/sougata_solver/staircase.py` is independently derived, per
  `rules.md` AI Coding Rule 1 — the technique itself (staircase/multi-slice
  approximation of a tapered sidewall) is standard/well-precedented and was
  already decided in `decisions.md` ADR-004, this phase just implements it.
  Three generator functions (`staircase_circle_layers`,
  `staircase_rectangle_layers`, `staircase_slab_layers`) each produce
  `num_slices` uniform-in-z `Layer`s of equal thickness, with the shape
  size linearly interpolated at each slice's **z-midpoint**
  (`frac = (i+0.5)/num_slices`) between a `top` and `bottom` value — a
  specific, documented convention choice (module docstring), not the only
  possible one (edge-averaging was considered and not used, no strong
  reason to prefer it). Because this phase adds no new eigenmode/Fourier
  formula, there is no external oracle to cross-check against (unlike every
  prior phase) — correctness instead rests on: (a) a zero-taper
  (`top_size == bottom_size`) regression test showing the staircase
  reproduces Phase 3/4a's already-oracle-validated single-uniform-layer
  result to `1e-10` regardless of `num_slices`, for all three shape types;
  (b) energy conservation for genuinely tapered cases; (c) two
  `slow`-marked convergence-vs-`num_slices` studies
  (`tests/test_staircase.py`) sweeping `num_slices` from 1 up to 32 (via)
  and 64 (trench) at a fixed wavelength/angle — both show the expected
  monotone-shrinking successive-`ΔR` trend (via: `R` settles to ≈0.565 by
  `N=16-32`; trench needs more slices to settle, ≈0.248 by `N=32-64`,
  consistent with its larger top-to-bottom taper ratio). Both
  `structures/via/tapered_via.py` and `structures/trench/tapered_trench.py`
  print (not plot, per ADR-010) the same `num_slices` sweep and confirm
  R+T=1.0000 at every point. 130 tests pass project-wide (118 prior + 12
  new — 10 fast + 2 `slow` — `tests/test_staircase.py`).

  **Follow-up, 2026-08-01**: renamed `tapered_trench.py`/`tapered_via.py`'s
  top/bottom-size constants to FDTD-style `TCD`/`BCD` (top/bottom critical
  dimension) plus `SPACING`, with `PERIOD = TCD + SPACING` derived, to match
  the equivalent Lumerical grating-structure-group parametrization the user
  referenced; added `structures/via/tapered_pillar.py` (the `Rectangle`
  case of `staircase_rectangle_layers`, not previously covered by an
  example script, with equal x/y halfwidths for a tapered square pillar).
  No new physics — same staircase generators, only naming/coverage. All
  three scripts re-verified end-to-end with R+T=1.0000.

## Phase 6 — Anisotropic Materials

- **Objectives**: support full 3×3 permittivity tensors (already exposed by
  `Material.epsilon_tensor`) in both uniform and patterned layers, removing
  `simulation.py`'s anisotropic `NotImplementedError`.
- **Deliverables**: generalize Phase 4a's eigensolver to accept a full
  tensor `Epsilon2` rather than only scalar/diagonal; validation against a
  known birefringent-material benchmark (e.g. a uniaxial waveplate at
  normal incidence, checked against closed-form ordinary/extraordinary
  index behavior).
- **Estimated complexity**: Medium — the eigensolver machinery from Phase 4a
  already generalizes to tensors mathematically; the work is mostly
  correctly wiring `Material.epsilon_tensor`'s off-diagonal terms through
  and validating the coupling terms are handled right.
- **Dependencies**: Phase 4a (reuses/extends its general eigensolver).
  Does not strictly require Phase 4b, though revisiting Phase 4b's
  near-degenerate handling once anisotropy is added is worth a follow-up
  check — anisotropic coupling can shift which cases are near-degenerate.
- **Status**: this phase is tracked at finer grain in
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 1 (added after this phase
  entry was originally written, for more scientific/mathematical rigor —
  see `progress_log.md` 2026-07-19). As of 2026-08-03, targets 1.1-1.4 and
  1.6-1.8 are shipped: uniform diagonal-tensor
  (`solve_layer_eigenmodes_uniform_diagonal`, `S4.cpp:1889-1906`), uniform
  in-plane-coupled (`solve_layer_eigenmodes_uniform_inplane`, cross-checked
  against a `RigorousCoupledWaveAnalysis.jl`-derived oracle to ~1e-13,
  `tests/oracles/rcwa_anisotropic_inplane_jl.py`), patterned anisotropic
  layers (`solve_layer_eigenmodes_patterned_inplane` +
  `fourier_factorization.toeplitz_matrix_component`, transcribed from a
  previously-unread `S4/S4/fmm/fmm_closed.cpp` branch, lines 165-256), a
  deterministic mode-ordering policy (`eigenmodes._canonical_mode_order`),
  and public propagating/evanescent mode classification
  (`eigenmodes.classify_propagating`,
  `SimulationResult.order_classification()`). Target 1.5 (longitudinal
  `eps_xz/eps_yz/eps_zx/eps_zy` coupling) is evaluated and **explicitly
  deferred** — a bounded literature search found no source both readable
  in this environment and independently benchmarkable (see
  `references.md`'s "Target 1.5 bounded literature search" entry); this
  phase's original "generalize to a full tensor `Epsilon2`" deliverable is
  therefore met for the diagonal/in-plane scope, not the fully general
  9-component tensor. 186 tests pass project-wide (123 at the start of
  this session, 65 new: `tests/test_anisotropic_uniform.py`,
  `tests/test_anisotropic_inplane.py`, `tests/test_anisotropic_patterned.py`,
  `tests/test_anisotropic_degeneracy.py`, `tests/test_mode_classification.py`).
  One honest, out-of-scope finding surfaced while validating target 1.8:
  a diffraction order sitting exactly at the Rayleigh/Wood's-anomaly
  threshold produces `NaN` R/T (a pre-existing `smatrix.py` division-by-zero
  at `q=0`, not introduced this session) — see `troubleshooting.md` and
  Category 6 target 6.4.

## Phase 7 — Real-Space Field Reconstruction & Visualization — **DONE**

- **Objectives**: reconstruct E/H(x,y,z) on a grid at an arbitrary depth in
  the stack, and produce cross-section field-intensity plots for trench/via
  structures.
- **Deliverables**: extend `fields.py` to sum Fourier components onto a
  real-space grid using `SMatrixStack.partial_smatrix_up_to` (already
  implemented, `smatrix.py:174-178`) for mode amplitudes at intermediate
  depths; example scripts producing `matplotlib` cross-section plots for a
  trench and a via.
- **Estimated complexity**: Medium — the amplitude bookkeeping is already
  in place; the new work is the inverse-Fourier-sum and its own
  correctness check (e.g. field continuity across a layer interface as a
  sanity test).
- **Dependencies**: Phase 3 or 4a (need at least one working patterned-layer
  case worth visualizing) — although the *machinery* for field
  reconstruction is dimension-agnostic and could technically be built
  against Phase 1's uniform stacks first as a stepping stone.
- **Status**: shipped, tracked at atomic-target grain as
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 9 (targets 9.1-9.8, all
  done 2026-08-05). `fields.modal_field_components`/`propagate_amplitudes`/
  `reconstruct_field_at_points` (transcribed from
  `S4/S4/rcwa.cpp::GetInPlaneFieldVector`/`GetFieldAtPoint`) and
  `smatrix.interior_amplitudes` (independently derived on top of
  `partial_smatrix_up_to`, `decisions.md` ADR-015) reconstruct full
  6-component `(Ex,Ey,Ez,Hx,Hy,Hz)` fields at any real-space point/line/grid
  and any depth. Validated against the analytic plane wave (uniform layer),
  exact tangential-field continuity across a genuine material interface,
  1D periodicity, and — the category's own exit criterion — real-space
  Poynting-flux integrals matching the solver's own `R`/`T` to `~1e-6` for
  a 2D pillar. **Honest finding along the way**: `fields.z_poynting_flux`
  turns out to be missing the textbook `0.5` time-average factor (harmless
  everywhere it's used, since `reflectance()`/`transmittance()` are ratios
  that cancel it) — never noticed before this phase because nothing
  earlier computed an *absolute* flux from raw E/H fields; see
  `CONVENTIONS.md`/`troubleshooting.md`. Two end-to-end example scripts
  (`structures/trench/trench_field_cross_section.py`,
  `structures/via/pillar_field_cross_section.py`) plus
  `postprocessing/plot_field_cross_section.py` (per `decisions.md`
  ADR-009/010's `structures/`-solves/`postprocessing/`-plots split) were
  run and their output field maps visually inspected, not just checked for
  "doesn't crash" — both show physically sensible near-field patterns
  (standing-wave lobes around a scattering pillar; periodic interference
  fringes through a grating). `tests/test_field_reconstruction.py`, 10
  tests. `fields.save_field_grid_npz` covers this phase's field-export
  need (NumPy `.npz`; CSV/HDF5 deliberately deferred, no schema designed).

## Phase 8 — Expanded Validation Suite & Example Gallery

- **Objectives**: systematic convergence-vs-`num_orders` studies for every
  geometry type; a complete example gallery mirroring the vendored
  `EMTutorial` reference cases (thin film, multistack/DBR, trench, via,
  pillar, tapered via).
- **Deliverables**: `slow`-marked convergence tests per geometry type;
  `structures/` entries for each structure type in `PRD.md`'s Success
  Criteria.
- **Estimated complexity**: Low-Medium — mostly systematic application of
  patterns already established in Phases 3-7, not new algorithmic risk.
- **Dependencies**: Phases 3-7 (validates all of them). Also the phase
  where every geometry type's energy-conservation and
  convergence-rate-vs-theory checks (`testing.md`'s Physical-Invariant
  Testing, first required starting Phase 3) get collected into one
  systematic sweep rather than living only in each phase's own test file.

## Phase 9 — Performance & Optional GPU/Autodiff Backend (later, optional)

- **Objectives**: vectorize wavelength/angle sweeps in NumPy first;
  optionally, add a torch/JAX array backend behind the same function
  signatures (Meent/TORCWA-style) for GPU batching and autodiff-based
  inverse design — only after correctness is fully validated through
  Phase 8.
- **Deliverables**: a vectorized-sweep code path with a regression test
  proving numerical equivalence to the unvectorized path (see `rules.md`'s
  Performance Requirements); optionally, a backend abstraction layer.
- **Estimated complexity**: High if the optional GPU/autodiff backend is
  pursued (requires re-expressing every eigensolve/S-matrix operation in a
  backend-agnostic way); Low-Medium for the NumPy-only vectorization step
  alone.
- **Dependencies**: Phases 2-8 (explicitly deferred until correctness is
  solid — see `decisions.md`).

## Phase 10 — Structure Visualization (3D Preview) — **DONE**

- **Objectives**: let a user visually inspect the actual geometry a
  `structures/*.py` script builds (via/pillar/trench, including tapered/
  staircased ones) as a 3D solid, before or instead of running a solve —
  a lightweight analogue of a commercial RCWA/FDTD tool's structure
  viewer. Requested directly by the project owner, scoped down (via
  `AskUserQuestion`) to a **static** 3D solid preview, Python desktop
  window, no new dependency; live/interactive parameter editing is an
  explicit, separate, not-yet-scoped follow-up — a future session should
  not assume this phase covers it.
- **Deliverables**: `plotting.plot_structure_3d(layer_stack, lattice, ...)`
  (new, `src/sougata_solver/plotting.py`), reusing a `_rasterize_pattern`
  helper factored out of Category 16's `plot_unit_cell` for in-plane
  cross-sections, rendered as stacked non-cubic `Axes3D.voxels` slabs at
  each layer's real z-offset; `postprocessing/plot_structure_3d_preview.py`
  demo/entry-point script; `tests/test_plotting.py` structural tests
  (return shape, z-extent, material-legend count, staircase voxel count,
  `Lattice1D` non-collapse regression, empty-stack error path).
- **Estimated complexity**: Low-Medium — no new physics/geometry formula
  (reuses already-validated rasterization and layer/staircase data
  structures unchanged), the real work was the 3D voxel-rendering
  mechanics and a measured performance characterization.
- **Dependencies**: Category 16 (Visualization, `plot_unit_cell`) for the
  rasterization logic it reuses; Phase 5 (staircase taper generation) for
  the tapered-structure rendering case.
- **Status**: shipped. `decisions.md` ADR-029 records the full design
  rationale, including two honest findings from this session: (a)
  `Lattice1D.b == (0, 0)` means `plot_unit_cell`'s own bounding-box logic
  already silently collapses to zero height for a 1D lattice — a
  pre-existing latent gap in shipped code, not introduced here, fixed in
  the new function via an explicit `extrusion_length` parameter (left
  unfixed in `plot_unit_cell` itself, out of this phase's scope); (b) a
  measured (not assumed) performance characterization —
  `Axes3D.voxels` scales worse than linearly in
  `resolution**2 * len(layer_stack)`, `resolution=20`/8 slices ~3s vs.
  `resolution=40`/16 slices ~2 minutes on the dev machine — documented in
  the function's own docstring and reflected in the demo script's chosen
  defaults. **Known, honestly-documented limitation**: this is an opaque
  solid render, so a via's tapered shaft is occluded by the surrounding
  substrate from a side view (only the top-face opening is visible) —
  confirmed by rendering and visually inspecting the demo script's output,
  not glossed over; a cutaway/transparency view would fix this but is out
  of scope for this static-preview target.
- **Follow-up, same day, shipped**: the project owner's own review of
  actual renders (not test failures) drove three more corrections to the
  static preview's rendering conventions before it was visually correct
  (fabricated `math.inf`-layer end-caps were tried, tuned, and ultimately
  **reverted** entirely — no size for a semi-infinite half-space reads
  correctly next to a finite patterned stack; material colors are now
  keyed by name, not encounter-order position, so the same material
  always renders the same color across every structure) — full account in
  `decisions.md` ADR-029's four "Correction" addenda. **Then, requested
  directly**: (1) a **generic** loader (`decisions.md` ADR-030) — every
  `structures/*.py` script now exposes `build_geometry(**overrides)`, so
  `postprocessing/plot_structure_3d_preview.py` works on *any* structure
  file, no per-structure Python needed, eliminating the drift-prone
  hand-copied-builder pattern ADR-029's fifth correction had caught, by
  construction rather than by discipline; (2) a **live** PyVista GUI
  (`decisions.md` ADR-031, new `postprocessing/live_structure_viewer.py`,
  new `gui` optional dependency) — sliders auto-generated from each
  script's own `build_geometry()` signature, real-time rebuild on change.
  Explicitly noted limitation: the live GUI's actual interactive-dragging
  experience has not been visually confirmed by a human in this sandbox
  (no display available) — only its data-prep/rebuild pipeline was
  verified programmatically in PyVista's offscreen mode; the project
  owner needs to run it themselves to confirm the live experience.

## Phase Sequencing Summary

```
Phase 1 (done) ──► Phase 2 ──► Phase 3 ──► Phase 4a ──► Phase 4b
                                   │            │             │
                                   │            └────► Phase 6 (extends Phase 4a)
                                   │            │
                                   │            └────► Phase 5 (needs Phase 4a only)
                                   │
                       Phase 3/4a ──────────► Phase 7 ──► Phase 8 ──► Phase 9 (optional)

Phase 5 (staircase) ─────────────────────────────────────► Phase 10
```

Phase 4b is a dependency of nothing else — it hardens Phase 4a's solver
in place (stress tests + stability handling) rather than unlocking new
capability, so Phases 5-8 only need Phase 4a to proceed. It should still
land before Phase 8's systematic convergence studies are trusted at high
`num_orders`/high-contrast settings, since that's exactly the regime Phase
4b is meant to have already stress-tested.

Phase 10 depends only on Category 16 (Visualization, for
`plot_unit_cell`'s rasterization it reuses) and Phase 5 (for the
tapered-structure rendering case) — it's a visualization-only addition
with no dependency on, or effect on, the solver phases above.
