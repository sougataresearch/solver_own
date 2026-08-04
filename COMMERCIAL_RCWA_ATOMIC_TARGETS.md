# Commercial RCWA Completion Plan — Atomic Targets

## How to use this register

Each target is deliberately small enough to implement, review, and test in
one isolated change. Complete targets in order within a category unless a
dependency is explicitly recorded as satisfied. A target is not done until
its stated check passes. Update this file, `phases.md`, `tasks.md`,
`memory.md`, and `progress_log.md` when it is completed.

Correctness and validation precede performance and convenience features.
Vendored projects under `../REFERENCE/` are read-only references/oracles,
never runtime dependencies.

---

## 1. Mathematical foundation — PARTIAL

### Already present

- Frequency-domain Maxwell formulation for isotropic periodic layers.
- Floquet-Bloch periodicity, reciprocal lattice vectors, Fourier basis, and
  wave-vector decomposition.
- Uniform, 1D, and 2D patterned-layer eigenmode formulations.
- Interface continuity, forward/backward modal amplitudes, S-matrix cascade,
  reflection, transmission, diffraction efficiencies, and power checks.

**Current scope**

Isotropic uniform, 1D, and 2D patterned RCWA is implemented. The remaining
anisotropy work is intentionally split into small, reversible milestones. Do
not remove the current anisotropic guard until the matching validation gate
passes.

### Small targets

- [x] **1.1 Reference audit:** compare vendored RCWA/TMM anisotropy
  formulations and record supported tensor components and limitations.
- [x] **1.2 Convention note:** document current field, time, propagation,
  normalization, and tensor-index conventions without changing solver behavior.
- [x] **1.3 Uniform diagonal tensor:** add a closed-form-validated uniaxial
  normal-incidence path; retain the isotropic path as the regression reference.
  Done 2026-08-03: `eigenmodes.solve_layer_eigenmodes_uniform_diagonal`,
  wired into `simulation.py`'s uniform-layer dispatch (`material.is_diagonal`
  branch); validated against a closed-form normal-incidence birefringence
  formula, the already-existing Fresnel/TMM oracle applied per principal
  axis, and reduction to `solve_layer_eigenmodes_uniform` when
  `eps_xx=eps_yy=eps_zz` (`tests/test_anisotropic_uniform.py`, 20 tests).
- [x] **1.4 In-plane tensor coupling:** support epsilon_xx, epsilon_xy,
  epsilon_yx, epsilon_yy, and epsilon_zz in uniform layers with an oracle test.
  Done 2026-08-03: `eigenmodes.solve_layer_eigenmodes_uniform_inplane`,
  wired into `simulation.py` (with a longitudinal-coupling guard naming
  target 1.5); validated against an independent eigenoperator oracle
  hand-transcribed from `RigorousCoupledWaveAnalysis.jl`
  (`tests/oracles/rcwa_anisotropic_inplane_jl.py`, agrees to ~1e-13),
  reduction to target 1.3's diagonal solver, and energy conservation for a
  Hermitian (lossless) in-plane-coupled slab (`tests/test_anisotropic_inplane.py`,
  15 tests).
- [ ] **1.5 Longitudinal coupling:** support epsilon_xz, epsilon_yz,
  epsilon_zx, and epsilon_zy only after obtaining a citable formulation and
  independent benchmark. **Evaluated and explicitly deferred, 2026-08-03**:
  per `rules.md` AI Coding Rule 1, a bounded literature search (this
  session, not just the standing `references.md` vendored-repo audit) was
  done before concluding. General-anisotropic-RCWA literature exists in
  principle (e.g. Glytsis & Gaylord 1987 JOSA A on 3D coupled-wave analysis
  of anisotropic gratings; a 1997-era Li paper on the same topic; general
  gyrotropic/bi-anisotropic RCWA formulations referenced in patent/thesis
  search results) but no source was found that is both (a) actually
  fetchable/readable as full text in this environment to transcribe exact
  equations from (an arXiv candidate, 2510.01214, returned only
  binary/undecodable PDF content via `WebFetch`; the JOSA A papers are
  paywalled), and (b) independently benchmarkable per Rule 5 (no second,
  structurally-different source to cross-check against, unlike targets
  1.3/1.4's S4 + RCWA.jl pairing). `simulation.py`'s longitudinal-coupling
  guard (added alongside target 1.4) continues to raise `NotImplementedError`
  naming this target rather than a plausible-but-unverified formula being
  written. Revisit if a readable, citable source becomes available (e.g. if
  the arXiv PDF can be decoded by another tool, or library/journal access
  becomes available) — this is a "not found this session," not a permanent
  "cannot exist" conclusion.
- [x] **1.6 Patterned anisotropic layers:** extend validated tensor handling to
  Toeplitz convolution matrices and the patterned-layer eigensolver.
  Done 2026-08-03: `fourier_factorization.pattern_epsilon_hat_component`/
  `toeplitz_matrix_component` and `eigenmodes.solve_layer_eigenmodes_patterned_inplane`,
  transcribed from `S4/S4/fmm/fmm_closed.cpp`'s `have_tensor` branch (lines
  165-256) -- a citation not previously read during the Phase 6 audit
  above (that audit only covered `S4.cpp`'s uniform path). Wired into
  `simulation.py`'s patterned-layer dispatch (isotropic pattern -> existing
  Phase 4a path unchanged; diagonal/in-plane-tensor pattern -> the new
  path; any longitudinal component anywhere in the pattern -> still
  `NotImplementedError` naming target 1.5). Validated by reduction to
  Phase 4a's isotropic solver, reduction to target 1.4's uniform-tensor
  solver for a spatially-uniform "pattern", and energy conservation for a
  genuinely patterned Hermitian (lossless) anisotropic pillar
  (`tests/test_anisotropic_patterned.py`, 12 tests).
- [x] **1.7 Degeneracy policy:** define and test deterministic handling for
  exact/near modal degeneracy. Done 2026-08-03: `eigenmodes._canonical_mode_order`
  applies a documented, deterministic sort key (rounded `Re(q)`, then
  `Im(q)`, then original `eig`-output index as a tie-break) to the three
  dense anisotropic eigensolvers (targets 1.3/1.4/1.6); builds on Phase 4b's
  existing `ILL_CONDITIONED_THRESHOLD` detection (unchanged) rather than
  replacing it. Explicitly does not claim eigenvalue continuity across a
  changing input (deferred to Category 2 target 2.3, "Sweep mode
  matching"). Validated by repeated-solve determinism, sort-key unit tests,
  and energy conservation for a deliberately near-isotropic (near-degenerate)
  patterned case (`tests/test_anisotropic_degeneracy.py`, 9 tests).
- [x] **1.8 Mode classification:** expose propagating/evanescent classification
  in public results and test the Rayleigh threshold. Done 2026-08-03:
  `eigenmodes.classify_propagating` and
  `SimulationResult.order_classification()`, reusing `_select_q_branch`'s
  own branch convention rather than a separate re-derivation. Validated
  against the analytic Rayleigh-threshold wavelength
  (`lambda_threshold = n_trans*period/m`) for a diffraction order flipping
  classification at the predicted point, plus energy conservation on both
  sides (`tests/test_mode_classification.py`, 7 tests). **Honest finding,
  not fixed by this target**: the exact threshold wavelength itself
  produces `NaN` R/T (a pre-existing `q=0` division-by-zero in
  `smatrix.py`, a genuine solver limitation belonging to Category 6 target
  6.4, not this target) — see `troubleshooting.md`.

### Exit criteria

**Category gate:** isotropic tensors reduce to current results; a uniaxial
benchmark agrees with ordinary/extraordinary theory; degeneracy regressions
remain deterministic and energy conserving.

**Status as of 2026-08-03**: targets 1.1-1.4 and 1.6-1.8 are done; 1.5
(longitudinal coupling) is evaluated and explicitly deferred (no citable +
independently-benchmarkable formulation found this session — see 1.5's own
entry and `references.md`). All three exit-criteria conditions above are
met for the scope actually shipped (diagonal + in-plane tensors, uniform
and patterned): isotropic reduction, uniaxial closed-form benchmark, and
deterministic/energy-conserving degeneracy handling are each covered by a
dedicated test file (`tests/test_anisotropic_uniform.py`,
`tests/test_anisotropic_inplane.py`, `tests/test_anisotropic_patterned.py`,
`tests/test_anisotropic_degeneracy.py`, `tests/test_mode_classification.py`
— 65 new tests this session, 186 total project-wide). The category gate
does not itself require longitudinal-coupling support, so this category is
considered closed for its shipped scope; target 1.5 remains individually
open and revisitable.

---

## 2. Numerical methods — PARTIAL

### Already present

- Redheffer star-product S-matrix cascading.
- Dense complex eigensolves, outgoing branch selection, and condition warnings.

**Current scope**

The stable isotropic solve is implemented; deterministic sweep tracking and
explicit numerical-failure policy remain.

### Small targets

- [x] **2.1 Failure contract:** document which numerical conditions raise
  `ValueError`/`LinAlgError` and which emit warnings; add unit tests.
  Done 2026-08-04: audited every `raise`/`logger.warning` call site in
  `src/sougata_solver/` (grep, not memory) into a new "Failure Contract"
  section in `design.md` (four tables: `ValueError`, `NotImplementedError`,
  `LinAlgError`, `WARNING`), plus `tests/test_failure_contract.py` (17
  tests, one per documented condition, including deliberately-singular-
  matrix cases confirming `LinAlgError` propagates uncaught).
- [x] **2.2 Eigenvalue report:** expose eigenvalue/mode-conditioning diagnostics
  in an internal result object without changing solve results.
  Done 2026-08-04: `layer.EigenmodeDiagnostics` (`cond_epsilon`, `cond_phi`,
  `min_eigenvalue_gap`, `num_propagating`, `num_evanescent`), attached as
  `LayerEigenmodes.diagnostics` (new optional field, default `None`) by
  every `eigenmodes.py` solver, reusing each solver's already-computed
  condition numbers where available. Validated by
  `tests/test_eigenvalue_diagnostics.py` (fields match independent
  recomputation; attaching diagnostics changes no other field).
- [x] **2.3 Sweep mode matching:** add deterministic mode ordering for a small
  wavelength sweep; test that the order does not arbitrarily permute.
  Done 2026-08-04, with an honest scope narrowing found during
  implementation: extending target 1.7's `_canonical_mode_order` to
  `solve_layer_eigenmodes_patterned` (Phase 4a, isotropic 2D) was tried
  first and broke two existing regression tests
  (`tests/test_2d_pillar.py`'s TE/TM-block tests, which depend on that
  solver's *natural* `eig()` output keeping a block structure that
  re-sorting destroys) — reverted per `rules.md` AI Coding Rule 3 (never
  weaken an existing oracle-comparison test to make a change pass); see
  `solve_layer_eigenmodes_patterned`'s docstring for the full account.
  Target 2.3 therefore stays scoped to the three anisotropic dense solvers
  that already carry target 1.7's ordering
  (`solve_layer_eigenmodes_uniform_diagonal`/`_inplane`,
  `solve_layer_eigenmodes_patterned_inplane`), validated by
  `tests/test_sweep_mode_matching.py`: a small non-degenerate wavelength
  sweep shows each mode's canonically-ordered trajectory changing
  smoothly (no discontinuous swap), not a claim of continuity through an
  eigenvalue crossing.
- [x] **2.4 Degeneracy warning:** detect a configurable small eigenvalue gap and
  emit one actionable warning; test both warning and no-warning paths.
  Done 2026-08-04: `eigenmodes.DEGENERATE_GAP_THRESHOLD` (`1e-6`,
  monkeypatch-configurable, same precedent as `ILL_CONDITIONED_THRESHOLD`)
  and `_warn_on_small_eigenvalue_gap`, applied to the same three
  anisotropic solvers as 2.3. Deliberately **not** applied to
  `solve_layer_eigenmodes_patterned` -- found during this target's own
  testing that an ordinary, otherwise well-conditioned circular-pillar
  case has a genuinely near-zero eigenvalue gap from routine `C4v`-lattice
  symmetry, which would misfire the warning on harmless, expected
  degeneracy; see that function's docstring. `tests/test_degeneracy_warning.py`
  covers both the natural near-degenerate case (no `monkeypatch` needed)
  and the forced-threshold mechanism test.
- [x] **2.5 Stress regression:** add one lossy high-contrast stress fixture and
  assert either valid conservation or a documented numerical failure.
  Done 2026-08-04, with a genuine sign-convention finding along the way:
  the first attempt reused Phase 4b's `n = -20+2j` "lossy-metal-like"
  index verbatim through a full `Simulation.solve()` (Phase 4b itself
  never called `solve()`, only cross-checked eigenvalues) and got
  `R+T` up to ~17 -- not a solver bug: `n=-20+2j` squares to
  `Im(eps) < 0`, which is a **gain** medium under this project's documented
  `d/dt -> -i*omega` phasor convention (`CONVENTIONS.md`), not a lossy one;
  Phase 4b's own eigenvalue-only test never exercised R/T so never caught
  the mislabel. Fixed by using a correctly-signed lossy metal
  (`eps = -396+80j`) for the new fixture (Phase 4b's already-shipped file
  left untouched, per `rules.md` AI Coding Rule 3). `tests/test_stress_regression.py`
  has two full-pipeline cases (isotropic lossy metal pillar; a lossy
  in-plane-coupled anisotropic pillar, the first stress test of target
  1.6's solver with a non-Hermitian/absorbing tensor) -- both pass the
  weaker-than-full-energy-balance passivity check (`R>=0`, `T>=0`,
  `R+T<=1`), since the full `R+T+A=1` identity needs layer-wise absorption
  (Category 7 targets 7.5/7.6), not yet implemented.

### Exit criteria

**Category gate:** all diagnostics are deterministic and no existing spectrum
changes beyond numerical tolerance.

**Status as of 2026-08-04**: all five targets (2.1-2.5) are done. 227 tests
pass project-wide (186 fast tests at the start of this session; 220 fast +
7 unchanged `slow` tests now -- 34 new fast tests: 17 in
`tests/test_failure_contract.py`, 6 in `tests/test_eigenvalue_diagnostics.py`,
4 in `tests/test_sweep_mode_matching.py`, 5 in `tests/test_degeneracy_warning.py`,
2 in `tests/test_stress_regression.py`). No existing
oracle-comparison or regression test was weakened to make a new one pass
(`rules.md` AI Coding Rule 3) -- one attempted change (extending canonical
ordering to the Phase 4a solver, target 2.3) was reverted rather than the
conflicting tests relaxed, and is documented as a negative finding, not
silently dropped. The category gate's "no existing spectrum changes beyond
numerical tolerance" is met: every new diagnostics/warning addition is
purely additive (a new optional dataclass field, new module constants, new
log lines) and the full pre-existing fast+slow suite (`pytest` /
`pytest -m slow`) was re-run and confirmed passing after every change, not
just at the end.

## 3. Fourier factorization — PARTIAL

### Already present

- Direct and inverse-rule Toeplitz permittivity matrices for current patterns.
- Analytic Fourier coefficients for slabs, circles, and rectangles.

**Current scope**

Current scalar factorization is validated; high-contrast 2D improvements are
not yet selected.

### Small targets

- [x] **3.1 Rule inventory:** record the chosen direct/inverse rule for every
  existing uniform, 1D, and 2D solver branch, with citations.
  Done 2026-08-04: `design.md`'s new "Fourier-factorization rule inventory"
  section (Algorithm 3a) tables every solver branch's rule (exact/direct/
  inverse-rule/numerical-matrix-inverse) with citations already established
  in `eigenmodes.py`, plus two findings from re-verifying end to end: (a)
  `epsilon_inv_hat` (the separately-Fourier-factorized inverse-rule
  Toeplitz) is consumed as such in exactly one place project-wide (the 1D
  TM block); every 2D path uses a numerical matrix-inverse of the
  *direct*-rule Toeplitz instead, a different operation despite similar
  naming; (b) confirmed via `tests/test_fourier_factorization_rules.py`
  (6 tests) that pins each table row against actual solver behavior
  (`LayerEigenmodes.epsilon_inv` black-box checks plus one white-box
  direct-vs-inverse-rule discrimination test for the 1D TE/TM blocks).
- [x] **3.2 1D convergence fixture:** add a fixed high-contrast lamellar case
  and record convergence versus harmonic order.
  Done 2026-08-04: `tests/test_fourier_convergence.py`, `n_ridge=10`
  (`eps=100`) TM-polarization binary grating, higher contrast than the
  existing Phase 3 `n=3.48` convergence check. Recorded (not fabricated)
  reflectance at `num_ord in {5,10,20,40,80,160,320}`; honest finding: not
  monotonic from `num_ord=5` (a real pre-asymptotic transient, reproduced
  deterministically), monotonic from `num_ord=10` onward but **not fully
  converged even at `num_ord=320`** (~6% relative error remaining at
  `num_ord=160` vs. that reference) -- an open illustration motivating
  targets 3.4/3.5, not a solver bug (the existing lower-contrast Phase 3
  test already establishes the correct converged limit).
- [x] **3.3 2D convergence fixture:** add a fixed high-contrast pillar case and
  record convergence versus harmonic order.
  Done 2026-08-04: `tests/test_fourier_convergence.py`, `n=5` (`eps=25`)
  circular pillar, `radius=0.2*period`. Recorded reflectance at
  `num_orders in {9,25,49,81,121,169,225}`; honest finding, kept in the
  record rather than dropped: `num_orders=25` gives `R=0.214`, an
  order-of-magnitude non-monotonic outlier against its low-order neighbors
  and the ~0.0236 converged value -- ordinary Laurent's-rule 2D Fourier
  factorization (no Li/NVM correction, per `solve_layer_eigenmodes_patterned`'s
  own docstring) can be not just imprecise but wildly wrong at very low
  truncation counts for a genuinely 2D pattern. Clean monotonic,
  shrinking-increment convergence only starts at `num_orders=49`.
- [x] **3.4 FFF feasibility decision:** compare the required formulation and
  available references; explicitly decide implement/defer Fast Fourier
  Factorization.
  Evaluated and explicitly deferred, 2026-08-04 -- see `decisions.md`
  ADR-012 and `references.md`'s Category 3 targets 3.4/3.5 entry for the
  full account. Popov & Nevière (2001)'s bibliographic details confirmed
  via `WebSearch`, but the paper itself is paywalled in this environment
  (not read/transcribed, per `rules.md` AI Coding Rule 1). `../S4` was read
  in full for its own implementation of this technique family instead
  (`fmm/fmm_PolBasisNV.cpp`/`fmm_PolBasisJones.cpp`/`fmm_PolBasisVL.cpp`,
  ~900 combined lines, all built on a discretized/FFT permittivity
  representation, not the analytic closed-form path this project already
  uses) -- a materially different architecture, in direct tension with
  the already-decided **ADR-002** (analytic Fourier transforms, raster+FFT
  explicitly rejected). Deferred, not implemented as part of target 3.6.
- [x] **3.5 Normal-vector feasibility decision:** compare formulations and a
  sharp-interface benchmark; explicitly decide implement/defer NVM.
  Evaluated and explicitly deferred, 2026-08-04, same investigation and
  ADR-012 as target 3.4 above (Lalanne 1997's NVM paper is the narrower 2D
  case of the same technique family S4 implements via
  `fmm_PolBasisNV.cpp`) -- see that ADR for the shared reasoning. The
  "sharp-interface benchmark" this target calls for is
  `tests/test_fourier_convergence.py`'s target-3.3 fixture itself, which
  already demonstrates the sharp-2D-interface convergence weakness NVM
  exists to fix; no further benchmark was needed to conclude the technique
  would help *if implemented* -- the decision to defer rests on
  implementation cost/architecture-fit (ADR-012), not on doubting NVM's
  applicability.
- [x] **3.6 Selected improvement:** implement only the technique approved by
  3.4/3.5, with an independent comparison and no regression in current cases.
  Done 2026-08-04 (no action needed — see reasoning): targets 3.4 and 3.5
  both concluded "defer," so target 3.6 has nothing approved to implement.
  Recorded as this target's own explicit outcome (matching the register's
  "explicitly decide implement/defer" allowance), not silently skipped --
  no new solver code was written for this target, and the existing 2D
  Laurent's-rule path is unchanged (zero regression risk since nothing
  changed).

### Exit criteria

**Category gate:** every factorization rule has a documented domain of use and
a convergence result.

**Status as of 2026-08-04**: all six targets (3.1-3.6) are done. 3.1-3.3
shipped new inventory documentation and two new test files
(`tests/test_fourier_factorization_rules.py`, 6 tests;
`tests/test_fourier_convergence.py`, 2 `slow` tests); 3.4/3.5 are evaluated
and explicitly deferred after a real investigation (S4's
`fmm_PolBasisNV`/`PolBasisJones`/`PolBasisVL` read in full, two literature
citations bibliographically verified via `WebSearch`) recorded in
`decisions.md` ADR-012; 3.6 has no implementation as a direct, explicit
consequence. The category gate is met for the scope actually investigated:
every factorization rule now has a documented domain of use (3.1's table)
and both new fixtures produced a real, recorded convergence result (3.2/3.3)
— the gate does not require 3.4/3.5 to conclude "implement," only that the
domain-of-use documentation and convergence results exist, which they do.
232 tests pass project-wide (227 at the start of this session: 220 fast + 7
slow -- 226 fast + 9 slow now).

## 4. Geometry engine — PARTIAL

### Already present

- 1D slabs and 2D circles/rectangles.
- Nested-shape handling and staircase taper geometry.

**Current scope**

The geometry engine covers current grating/via examples; general primitives
and imported profiles remain pending.

### Small targets

- [ ] **4.1 Geometry validation API:** reject non-finite dimensions, invalid
  lattice vectors, and invalid shape sizes; add unit tests.
- [ ] **4.2 Unit-cell bounds policy:** define periodic wrapping/overlap behavior
  and test shapes crossing a cell edge.
- [ ] **4.3 Ellipse primitive:** add an ellipse with DC-area and nonzero Fourier
  coefficient tests.
- [ ] **4.4 Polygon design:** choose analytic versus controlled raster/FFT
  coefficients; document the accuracy contract before implementation.
- [ ] **4.5 Polygon primitive:** implement one simple polygon path and validate
  its area/Fourier coefficients against numerical integration.
- [ ] **4.6 Imported geometry format:** define a minimal, safe format with units
  and validation; add parser tests before solver integration.
- [ ] **4.7 Profile slicing API:** define and test a geometry-to-layer-slices
  interface independently of the RCWA solve.

### Exit criteria

**Category gate:** each new shape has geometry-only tests and one end-to-end
RCWA example.

## 5. Material models — PARTIAL

### Already present

- Constant complex permittivity/index and tabulated isotropic n,k data.
- A `Material` interface that stores scalar or 3x3 tensor values.

**Current scope**

Tensor values can be represented but cannot yet be solved; analytic dispersion
models are not yet implemented.

### Small targets

- [ ] **5.1 Material validation:** validate tensor shape, finite values, and
  wavelength callback output; add negative-input tests.
- [ ] **5.2 Sellmeier model:** implement one documented Sellmeier form and test
  known refractive-index values.
- [ ] **5.3 Cauchy model:** implement and test a documented Cauchy form.
- [ ] **5.4 Lorentz model:** implement a single-resonance Lorentz oscillator
  with a causality/sign-convention test.
- [ ] **5.5 Drude model:** implement and test a Drude model against a published
  or tabulated reference curve.
- [ ] **5.6 Drude-Lorentz composition:** combine existing oscillator building
  blocks and test that zero-strength terms reduce correctly.
- [ ] **5.7 Tensor-material solver wiring:** complete only after Category 1's
  corresponding tensor milestone and benchmark pass.
- [ ] **5.8 Material provenance:** add optional source/citation metadata to
  material definitions and serialized output.

### Exit criteria

**Category gate:** every analytic material model has unit tests and a
wavelength-sweep example.

## 6. Boundary conditions and excitation — PARTIAL

### Already present

- Floquet boundaries, arbitrary incidence angle/azimuth, and Jones-vector input.
- TE/TM, linear, circular, and elliptical polarization amplitudes.

**Current scope**

The existing convention is internally consistent; external polarization-phase
validation remains required.

### Small targets

- [ ] **6.1 Convention note:** publish the existing time, propagation, and
  s/p convention in one reference document.
- [ ] **6.2 Normal-incidence polarization regression:** test TE, TM, linear,
  circular, and elliptical states against known symmetry expectations.
- [ ] **6.3 Oblique-incidence regression:** test arbitrary azimuth and mixed
  polarization with energy conservation.
- [ ] **6.4 Grazing-incidence boundary test:** define a supported near-grazing
  angle range and test behavior at its boundary.
- [ ] **6.5 Rayleigh-threshold test:** add a diffraction-order opening/closing
  case and test propagating versus evanescent classification.
- [ ] **6.6 Bottom-incidence decision:** document whether reverse illumination
  is required; design it separately if approved.

### Exit criteria

**Category gate:** polarization-sensitive external comparison passes for at
least one oblique case.

## 7. Layer handling — PARTIAL

### Already present

- Uniform/patterned finite layers, semi-infinite half-spaces, and multilayer stacks.
- Multi-slice tapered via and trench stacks.

**Current scope**

Layer assembly is correct-first and uncached; layer-specific loss reporting is
not available.

### Small targets

- [ ] **7.1 Layer validation audit:** test finite-layer thickness, half-space,
  and patterned-layer invariants at construction.
- [ ] **7.2 Repeated-layer identity:** add a test that equivalent repeated-layer
  representations produce the same R/T.
- [ ] **7.3 Layer cache design:** define cache keys and invalidation conditions;
  do not implement caching yet.
- [ ] **7.4 Layer cache implementation:** cache one safe artifact (for example,
  a Toeplitz matrix) and test equivalence to uncached results.
- [ ] **7.5 Layer-wise absorption design:** define the physical quantity and
  validation method before exposing an API.
- [ ] **7.6 Layer-wise absorption implementation:** add only after the design
  has an independent energy-balance test.

### Exit criteria

**Category gate:** cached and uncached paths agree within numerical tolerance.

## 8. Solver sweeps and convergence — PARTIAL

### Already present

- Runnable examples with wavelength and staircase-slice sweeps.
- Manual harmonic-order selection through `num_orders`.

**Current scope**

Sweeps live mostly in scripts; library-level sweep and automatic-convergence
support are pending.

### Small targets

- [ ] **8.1 Result-series container:** define a typed result container for a
  one-parameter sweep, including units and run metadata.
- [ ] **8.2 Wavelength sweep API:** promote the existing script loop into a
  library-level wavelength sweep; compare to scalar calls.
- [ ] **8.3 Angle sweep API:** add theta/phi sweeps using the same container.
- [ ] **8.4 Polarization sweep API:** add a finite list of Jones states.
- [ ] **8.5 Thickness/geometry sweep API:** add one named parameter sweep with
  explicit input validation.
- [ ] **8.6 Harmonic-study API:** return R/T and conservation residual versus
  harmonic order without automatic decisions.
- [ ] **8.7 Convergence criterion:** define a conservative stopping criterion;
  validate it against manually extended studies.
- [ ] **8.8 Automatic harmonic selection:** implement only after 8.7 succeeds
  on thin-film, trench, and pillar fixtures.

### Exit criteria

**Category gate:** each sweep is equivalent to repeated scalar solves and
stores enough data to reproduce its convergence decision.

## 9. Field calculations — NOT COMPLETED

### Already present

- Modal amplitudes and z-directed Poynting-flux evaluation for R/T.
- Partial S-matrix support that can supply internal interface amplitudes.

**Current scope**

Real-space E/H reconstruction and field exports are not implemented.

### Small targets

- [ ] **9.1 Field-convention design:** document field-component ordering and
  Fourier phase convention used by the existing modal matrices.
- [ ] **9.2 Uniform-layer reconstruction:** reconstruct E/H at one depth in a
  uniform layer and test against the analytic plane wave.
- [ ] **9.3 Interface amplitudes:** recover modal amplitudes at one internal
  interface using the S-matrix partial-stack API; test a uniform stack.
- [ ] **9.4 1D inverse Fourier sum:** reconstruct a 1D patterned field grid and
  test periodicity.
- [ ] **9.5 Field continuity test:** test tangential-field continuity at an
  interface where no surface current exists.
- [ ] **9.6 Field-flux test:** integrate reconstructed flux and compare with
  existing R/T.
- [ ] **9.7 2D field grid:** extend the verified 1D reconstruction to a 2D grid.
- [ ] **9.8 Field export:** add NumPy first, then CSV/HDF5 only after schema design.

### Exit criteria

**Category gate:** field-derived flux agrees with solver R/T on uniform, 1D,
and 2D fixtures.

## 10. Optical outputs — PARTIAL

### Already present

- Total R/T, per-order diffraction efficiencies, Jones/Mueller processing,
  and ellipsometric quantities.

**Current scope**

The primary power outputs exist; standardized complex-order and loss outputs
remain pending.

### Small targets

- [ ] **10.1 Complex coefficients:** expose complex per-order reflection and
  transmission amplitudes with an explicitly documented convention.
- [ ] **10.2 Diffraction angles:** add angles for propagating orders and a clear
  non-propagating representation for evanescent orders.
- [ ] **10.3 Conservation report:** expose incident/reflected/transmitted/loss
  components and residual in one result method.
- [ ] **10.4 Loss accounting design:** define layer-wise absorption before adding
  values to the public result.
- [ ] **10.5 Polarization conversion:** expose per-order s/p conversion only
  after the polarization convention is externally validated.
- [ ] **10.6 Output schema tests:** freeze a compact output fixture for uniform,
  1D, and 2D cases.

### Exit criteria

**Category gate:** all public outputs carry units/conventions and satisfy a
conservation check where physics permits.

## 11. Semiconductor OCD features — PARTIAL

### Already present

- Trenches, vias, pillars, and linear sidewall taper via staircase layers.

**Current scope**

Basic geometry exists, but dedicated CD/rounding/overlay/roughness workflows
are not implemented.

### Small targets

- [ ] **11.1 OCD parameter object:** define validated CD, period, height,
  material, and sidewall-angle parameters.
- [ ] **11.2 Trapezoid constructor:** generate a staircase trench from top CD,
  bottom CD, height, and slice count; test zero taper.
- [ ] **11.3 Corner-rounding design:** choose a periodic geometry approximation
  and define its convergence parameter.
- [ ] **11.4 Corner-rounding implementation:** add one geometry with area and
  convergence tests.
- [ ] **11.5 TSV template:** create a reproducible via/OCD example sweep.
- [ ] **11.6 Grating template:** create a reproducible trench/OCD example sweep.
- [ ] **11.7 Overlay feasibility decision:** define a periodic unit-cell model or
  explicitly defer overlay.
- [ ] **11.8 LER/LWR feasibility decision:** define deterministic periodic
  approximations or explicitly defer stochastic roughness.

### Exit criteria

**Category gate:** parameter changes are traceable in metadata and covered by
at least one convergence study.

## 12. Linear algebra — PARTIAL

### Already present

- NumPy/SciPy dense complex eigensolves, linear solves, and condition estimates.

**Current scope**

The solver uses dense CPU linear algebra; profiling and alternative
factorizations have not been systematically evaluated.

### Small targets

- [ ] **12.1 Baseline profiler:** measure eigensolve, matrix-solve, and
  S-matrix time/memory on fixed fixtures.
- [ ] **12.2 Direct-inverse audit:** replace only demonstrably unnecessary
  explicit inverses with linear solves; add equivalence tests.
- [ ] **12.3 Factorization reuse design:** identify safe intra-solve reuse
  opportunities without adding global cache state.
- [ ] **12.4 SVD diagnostic:** add an opt-in singular-value diagnostic for
  troublesome eigenvector matrices.
- [ ] **12.5 Sparse feasibility decision:** use benchmark dimensions to decide
  whether sparse/iterative methods are worthwhile.

### Exit criteria

**Category gate:** any algebra change preserves existing oracle tests and has
a measured conditioning or runtime justification.

## 13. Performance optimization — NOT COMPLETED

### Already present

- A pure NumPy/SciPy CPU baseline suitable for correctness validation.

**Current scope**

No profiling-driven caching, vectorization, parallelism, or GPU path exists.

### Small targets

- [ ] **13.1 Benchmark suite:** establish repeatable runtime/memory benchmarks
  for thin-film, trench, pillar, and tapered structures.
- [ ] **13.2 Fourier-matrix reuse:** cache one validated sweep-invariant matrix;
  compare cached and uncached spectra.
- [ ] **13.3 Eigenmode reuse:** cache one validated sweep-invariant mode set.
- [ ] **13.4 Safe vectorized sweep:** vectorize a single simple sweep and compare
  every result with scalar solves.
- [ ] **13.5 Parallelism decision:** profile and document whether process/thread
  parallelism gives a safe benefit.
- [ ] **13.6 GPU decision checkpoint:** seek explicit approval and choose a
  backend only after CPU targets are met.

### Exit criteria

**Category gate:** every optimization has numerical-equivalence tests and a
benchmark showing why it exists.

## 14. Validation — PARTIAL

### Already present

- Analytic Fresnel/TMM tests, 1D/2D eigenvalue oracles, energy conservation,
  stress tests, and staircase convergence tests.

**Current scope**

The major gap is a fully external 2D R/T oracle and systematic convergence
coverage across every feature family.

### Small targets

- [ ] **14.1 Validation inventory:** map every public feature to its oracle,
  invariant test, example, and known limitation.
- [ ] **14.2 External 2D R/T oracle:** make S4 runnable, locate published data,
  or import versioned commercial reference results; document provenance.
- [ ] **14.3 Moderate 2D R/T test:** compare one pillar/via case to 14.2.
- [ ] **14.4 High-contrast 2D R/T test:** compare one stress case to 14.2.
- [ ] **14.5 Reciprocity test design:** select applicable reciprocal lossless
  cases and define expected symmetry.
- [ ] **14.6 Reciprocity tests:** add the selected cases without assuming they
  apply to non-reciprocal/future gain media.
- [ ] **14.7 Harmonic convergence matrix:** run documented studies across every
  supported geometry family.
- [ ] **14.8 Validation report:** publish tolerances, versions, and results.

### Exit criteria

**Category gate:** every marketed capability has a named validation source and
an automated regression test.

## 15. User interface and API — PARTIAL

### Already present

- Python API, runnable structure scripts, CSV output, and run metadata.

**Current scope**

There is no stable configuration schema, CLI, or structured field-result
export yet.

### Small targets

- [ ] **15.1 Public API inventory:** list supported public classes/functions
  and identify unstable/internal interfaces.
- [ ] **15.2 Configuration schema:** define a minimal JSON/YAML schema without
  adding parsing dependencies unless justified.
- [ ] **15.3 Configuration validation:** validate a configuration before any
  numerical calculation; add malformed-input tests.
- [ ] **15.4 Configuration runner:** reproduce one existing thin-film example
  from a configuration file.
- [ ] **15.5 CLI design:** specify commands, exit codes, and output locations.
- [ ] **15.6 CLI implementation:** add one `run` command for the validated
  configuration workflow.
- [ ] **15.7 NumPy export:** serialize a result series and metadata to NumPy.
- [ ] **15.8 HDF5 decision/implementation:** choose and add HDF5 only when
  structured field/sweep data requires it.

### Exit criteria

**Category gate:** an official example is reproducible through the documented
public API without editing source code.

## 16. Visualization — PARTIAL

### Already present

- Selected R/T and polarimetry post-processing plots.

**Current scope**

Plotting is not yet a systematic, saved-data workflow and no field plots exist.

### Small targets

- [ ] **16.1 Plot data contract:** define plotting inputs separate from solver
  calculation and output files.
- [ ] **16.2 Geometry plot:** render the unit cell and layer stack; test axes,
  units, and no-solver dependency.
- [ ] **16.3 R/T spectrum plot:** formalize an existing post-processing plot
  with metadata labels.
- [ ] **16.4 Harmonic-convergence plot:** plot the Category-8 study result.
- [ ] **16.5 Diffraction-order plot:** visualize propagating-order efficiency.
- [ ] **16.6 Field-intensity plot:** add only after Category 9 field-flux tests pass.
- [ ] **16.7 Poynting/phase plots:** add after field component conventions are
  validated.

### Exit criteria

**Category gate:** plots are reproducible from saved result data and never
alter solver numerical results.

## 17. Testing and quality — PARTIAL

### Already present

- 130 passing unit, integration-style, oracle, invariant, and slow tests.

**Current scope**

The project has a strong local suite but no CI, static-analysis baseline, or
performance regression guard.

### Small targets

- [ ] **17.1 Test taxonomy:** consistently mark unit, integration, oracle, and
  slow tests in filenames/docstrings or pytest markers.
- [ ] **17.2 Windows CI:** run the fast suite on supported Python versions.
- [ ] **17.3 Slow-test CI policy:** schedule or manually trigger convergence
  studies separately from the fast suite.
- [ ] **17.4 Regression fixtures:** store compact trusted spectra/field outputs
  with provenance and tolerance rationale.
- [ ] **17.5 Static-analysis setup:** install/configure lint/type tools and fix
  the baseline before making them required.
- [ ] **17.6 Performance regression guard:** add only after Category 13 has
  stable benchmark baselines.

### Exit criteria

**Category gate:** CI protects the fast test suite and no documented oracle
test can be silently skipped.

## 18. Documentation — PARTIAL

### Already present

- README, PRD, architecture, design, decisions, testing, references, and
  troubleshooting documents.

**Current scope**

Internal documentation is strong; a cohesive mathematical theory guide and
user tutorial/API reference remain incomplete.

### Small targets

- [ ] **18.1 Theory outline:** create a table of contents for the mathematical
  derivation, conventions, and validation notes.
- [ ] **18.2 Core derivation:** document the isotropic uniform/1D/2D equations
  already implemented, citing code and sources.
- [ ] **18.3 Anisotropy derivation:** add only after the corresponding Category-1
  milestones are validated.
- [ ] **18.4 API reference:** generate or maintain a reference for public API
  signatures, units, and exceptions.
- [ ] **18.5 Tutorial: thin film:** reproduce a validated Fresnel/TMM example.
- [ ] **18.6 Tutorial: grating:** reproduce a validated 1D diffraction example.
- [ ] **18.7 Tutorial: via/taper:** reproduce validated 2D and staircase examples.
- [ ] **18.8 Validation guide:** explain what each benchmark proves and does not prove.

### Exit criteria

**Category gate:** a new user can reproduce a validated example without
reading implementation source.

## 19. Future extensions — DEFERRED

### Already present

- A documented future-scope list; no experimental backend has been added.

**Current scope**

These are intentionally deferred until the correctness and usability targets
in Categories 1–18 are achieved.

### Small targets

- [ ] **19.1 Use-case decision:** select one research objective before starting
  any advanced extension.
- [ ] **19.2 Adjoint/autodiff feasibility study:** compare backend options and
  differentiation limitations without changing the core solver.
- [ ] **19.3 Optimization prototype:** only after a validated gradient or
  derivative-free baseline is chosen.
- [ ] **19.4 GPU feasibility study:** only after Category 13 CPU profiling.
- [ ] **19.5 Nonlinear/time-modulated scope decision:** define required physics,
  validation source, and whether it belongs in this project.
- [ ] **19.6 Hybrid RCWA/FDTD/FEM scope decision:** define coupling interface and
  independent validation before any implementation.

### Exit criteria

**Category gate:** future work starts only with an explicit use case, a source,
and a validation plan; it never bypasses Categories 1–18.
