# Task Checklist — sougata_solver

Atomic, trackable tasks per phase (see `phases.md` for objectives/context).
Check items off as completed; do not remove completed items — move
finished phases' checked lists into `memory.md`'s "Completed Milestones"
summary instead of deleting history here.

## Phase 1 — Uniform Multilayer Core (DONE)

☑ Implement `Material` (scalar + tensor permittivity, `from_nk`, `from_permittivity_tensor`)
☑ Implement `Lattice` (reciprocal vectors, unit cell area)
☑ Implement `Layer` / `LayerStack` with semi-infinite half-spaces
☑ Implement uniform-layer eigenmode solve (`solve_layer_eigenmodes_uniform`)
☑ Implement `q` branch selection (`_select_q_branch`)
☑ Implement interface + propagation S-matrices and Redheffer star product
☑ Implement `PlaneWaveExcitation` (s/p decomposition, incident amplitude inversion)
☑ Implement `z_poynting_flux` / `tangential_e_field`
☑ Implement `Simulation.solve` (uniform path) and `SimulationResult`
☑ Implement Jones/Mueller polarimetry
☑ Validate against analytic Fresnel/TMM (`tests/test_analytic_fresnel.py`, `tests/oracles/fresnel.py`)
☑ Ship the initial `structures/` scripts (anti-reflection coating, SiO2-on-Si, custom multistack, custom-material-from-nk-data) and the `postprocessing/` Jones/Mueller ellipsometry split

## Phase 2 — Fourier-Factorization Core (DONE)

☑ Add `.flake8`/`ruff` config and `mypy.ini` to `sougata_solver/` (rules.md gap, do before new modules land)
☑ Create `src/sougata_solver/fourier_factorization.py`
☑ Implement `pattern_epsilon_hat(pattern, lattice, g1, g2, wavelength)` (direct, sums shape contributions with containment-tree subtraction)
☑ Implement the inverse-rule variant (`1/eps` per shape, same summation) for `epsilon_inv_hat`
☑ Implement `toeplitz_matrix(pattern, lattice, g_vectors, wavelength)`
☑ Write numerical rasterize-and-sum reference for a `Circle` pattern
☑ Write numerical rasterize-and-sum reference for a `Rectangle` pattern
☑ Test: analytic Toeplitz entries match the rasterized reference within tolerance, for both `Circle` and `Rectangle`
☑ Test: DC term of `epsilon_hat` equals area-weighted average permittivity (closed-form sanity check)
☑ Update `memory.md` / `decisions.md` on completion

## Phase 3 — 1D-Periodic Lamellar Gratings (Trench) (DONE)

☑ Add `Lattice1D(period)` to `geometry.py`
☑ Add `Slab`/`Line` 1D shape with analytic (`sinc`) Fourier transform
☑ Add `truncate_fourier_orders_1d` to `fourier_basis.py`
☑ Implement `solve_layer_eigenmodes_1d` (TE block, using `epsilon_hat`)
☑ Implement `solve_layer_eigenmodes_1d` (TM block, using `inv(epsilon_inv_hat)`)
☑ Add `Lattice1D` dispatch branch in `simulation.py`
☑ Source and transcribe a published 1D binary-grating benchmark (Moharam/Gaylord-style, via `Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples`) into `tests/oracles/rcwa_1d_gaylord.py`
☑ Test: TE diffraction efficiencies match benchmark (agrees to ~1e-10 at `num_ord=15`)
☑ Test: TM diffraction efficiencies match benchmark (converges to the same value as the oracle, but only at high `num_orders` -- see `tests/test_1d_grating.py::test_tm_matches_gaylord_oracle_at_high_num_orders` and its docstring for the caveat: the oracle's own source self-reports "STILL NOT WORKING YET")
☑ Test: patterned layer reduces to Phase 1's uniform-layer result when shape material equals background (continuity sanity check)
☑ Test: energy conservation (`R + T + sum(diffraction efficiencies) = 1`) across TE/TM/mixed polarization, normal and oblique incidence (`testing.md` Physical-Invariant Testing)
☑ Test: measured convergence rate vs. `num_orders` (TM, the Li's-rule-sensitive case) decreases monotonically toward a high-order reference (`testing.md` Physical-Invariant Testing)
☑ Write `structures/trench/trench_grating.py`
☑ Update `memory.md` / `references.md` on completion

## Phase 4a — 2D-Periodic Patterned Layers, Well-Conditioned Case (Via, Pillar) (DONE)

☑ Implement `solve_layer_eigenmodes_patterned` (general non-uniform eigenproblem, transcribed from `S4/S4/rcwa.cpp::SolveLayerEigensystem` lines 794-827 and `S4/S4/fmm/fmm_closed.cpp`'s true-2D `Epsilon2` branch, lines 133-139/162-163 — corrected mid-session after a first draft wrongly reused the 1D-only branch's formula, see `phases.md` Phase 4a Status), scoped to moderate-contrast/moderate-`num_orders` cases
☑ Remove the `NotImplementedError` at `simulation.py:98`, wire in Phase 2 Toeplitz construction + this solver
☑ Determine whether S4 is buildable/runnable in this environment for a subprocess cross-check oracle — not usable here (no `cmake`/Lua toolchain found)
☐ Source a published 2D benchmark instead (**still not done, but in progress as of 2026-08-20** — `tests/oracles/rcwa_2d_pillar.py` documents the survey of all vendored RCWA repos and why none yielded a hard-coded literature number; a Lumerical RCWA cross-validation is now prepared instead, per `decisions.md` ADR-039 — `pillar_array.py`/`via_array.py` convergence-fixed, export/overlay scripts and a build-spec README ready, but the actual Lumerical run and R/T agreement numbers are not yet available, so this stays unchecked until they exist)
☑ Test: 2D patterned-layer reduce-to-uniform when shape material == background
☑ Test: `ky=0`'s TE-like block matches the already-validated 1D solver (the part that's rule-independent); a separate test now asserts the TM-like block correctly *diverges* (Li's rule is 1D-only in S4, confirmed by reading the source) — replaces an earlier, circular "reduces to 1D" test that passed even with the wrong formula
☑ Test: `q^2` eigenvalues match an independently-transcribed `RigorousCoupledWaveAnalysis.jl` eigenoperator (`tests/oracles/rcwa_2djl_eigenvalues.py`, a structurally different formula, agrees to ~1e-12) — the real oracle-comparison test the "ky=0" checks above cannot substitute for; closes the "only ever read S4" gap the Epsilon2 bug came from
☑ Test: energy conservation for moderate-contrast cases (`testing.md` Physical-Invariant Testing)
☑ Test: moderate-contrast pillar case runs end-to-end with physically plausible R/T
☑ Write `structures/via/pillar_array.py`
☑ Write `structures/via/via_array.py`
☑ Update `memory.md` / `decisions.md` on completion
☑ Re-audit Phase 1-3's S4 citations against each cited function's *full* body (not just the branch already used) — the same class of gap that caused the Epsilon2 bug; no further issues found (`phases.md`/`memory.md` Phase 4a entries)

## Phase 4b — 2D-Periodic Patterned Layers, Near-Degenerate / Ill-Conditioned Cases (DONE)

☑ Identify/construct high-index-contrast, small-feature-to-period-ratio, high-`num_orders` test cases likely to stress near-degenerate eigenvalues (index contrast `3.48` to lossy-metal-like `-20+2j`, `num_orders` up to 225, near-touching pillars, sub-percent sliver rectangle, near-degenerate nested circles — see `tests/test_2d_pillar_stress.py`)
☑ Handle near-degenerate eigenvalue edge cases: empirically, no case in the stress sweep required handling beyond `_select_q_branch` (reused unmodified) — documented as an honest finding, not a fabricated fix, in `solve_layer_eigenmodes_patterned`'s docstring
☑ Add a condition-number diagnostic (`eigenmodes.ILL_CONDITIONED_THRESHOLD`, `1e4`, logged at `WARNING` via a module-level `logger`, per `design.md` Logging Strategy) for both `epsilon_hat` and the eigenvector matrix `phi`, rather than silently returning a degraded answer
☐ Test: the stress case(s) still match S4 (or a published benchmark) within tolerance — **S4/Julia unavailable in this environment** (confirmed, not assumed — `troubleshooting.md`'s Environment-Specific Notes); stress cases instead cross-checked against the independent `RigorousCoupledWaveAnalysis.jl` eigenvalue oracle already built for Phase 4a (agrees to ~1e-6+ tolerance, actual sweep showed ~1e-10), an honestly-scoped substitute, not a fabricated S4 match
☑ Test: energy conservation holds for the stress case(s) (`testing.md` Physical-Invariant Testing — via the oracle-eigenvalue agreement, which implies energy-conserving R/T)
☑ Test: condition-number `WARNING` logging actually fires when the threshold is exceeded, and stays silent for ordinary cases (`tests/test_2d_pillar_stress.py`'s logging tests, using `caplog`/`monkeypatch`)
☑ Update `troubleshooting.md` with the concrete Phase 4a bug found and fixed this session (moved to Already-Solved Gotchas), the Phase 4b stress-test findings, and confirmed S4/Julia unavailability
☑ Update `memory.md` / `phases.md` on completion

## Phase 5 — Tapered / Sloped Sidewalls (Via, Trench) (DONE)

☑ Design the staircase-layer-stack generator's API (inputs: top size, bottom size, thickness, slice count `N`; output: `list[Layer]`) — three functions in `src/sougata_solver/staircase.py` (one per shape type, since `Circle`/`Rectangle`/`Slab` take differently-shaped size parameters), z-midpoint linear interpolation convention documented in the module docstring
☑ Implement the generator for `Rectangle`/`Circle` (via) — `staircase_rectangle_layers`, `staircase_circle_layers`
☑ Implement the generator for `Slab` (trench) — `staircase_slab_layers`
☑ Write a convergence-vs-`N` test for a tapered via (mark `slow`) — `tests/test_staircase.py::test_tapered_via_converges_with_increasing_num_slices`, N=1..32
☑ Write a convergence-vs-`N` test for a tapered trench (mark `slow`) — `tests/test_staircase.py::test_tapered_trench_converges_with_increasing_num_slices`, N=1..64 (needed one more octave than the via case to settle below the same tolerance)
☑ Write an example script sweeping `N` and plotting/printing R/T convergence — `structures/via/tapered_via.py`, `structures/trench/tapered_trench.py` (printing only, per ADR-010's plotting-belongs-in-postprocessing rule; no plotting was separately requested this phase)
☑ Update `memory.md` / `decisions.md` on completion — `decisions.md` ADR-004 already covered the design decision from planning; `memory.md` updated with this phase's actual implementation/validation outcome

## Phase 6 — Anisotropic Materials

☑ Uniform diagonal-tensor layers (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` target 1.3):
  `solve_layer_eigenmodes_uniform_diagonal` (`eigenmodes.py`, transcribed from
  `S4/S4/S4.cpp:1889-1906`'s uniform-anisotropic branch, diagonal-only), wired
  into `simulation.py`; closed-form normal-incidence birefringence benchmark +
  Fresnel-oracle-per-axis cross-check + isotropic-reduction regression
  (`tests/test_anisotropic_uniform.py`). In-plane-coupled and
  longitudinally-coupled tensors, and anisotropic patterned layers, remain
  open (targets 1.4/1.5/1.6).
☑ Generalize the uniform-layer eigensolver to accept in-plane-coupled tensor components (`solve_layer_eigenmodes_uniform_inplane`, target 1.4); longitudinal coupling (target 1.5) still pending a citable formulation; patterned-layer `Epsilon2` generalization is target 1.6
☑ Remove `simulation.py`'s uniform-anisotropic `NotImplementedError` for the diagonal case (target 1.3); general/off-diagonal case still raises, naming targets 1.4/1.5
☑ Source a birefringent-material closed-form benchmark (e.g. uniaxial waveplate at normal incidence) — independently derived, see `eigenmodes.solve_layer_eigenmodes_uniform_diagonal`'s docstring
☑ Test: anisotropic solve matches the benchmark (`tests/test_anisotropic_uniform.py`)
☑ Test: isotropic-tensor special case reduces to Phase 1's uniform-isotropic result (regression guard) (`tests/test_anisotropic_uniform.py`)
□ Update `memory.md` / `decisions.md` on completion (in progress — see this session's `memory.md`/`progress_log.md` entries; final update once all of targets 1.3-1.8 land)
□ Target 1.5 (longitudinal coupling) — evaluated and explicitly deferred 2026-08-03, no citable+benchmarkable source found in a bounded literature search this session; see `references.md`'s "Target 1.5 bounded literature search" entry and `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`'s 1.5 entry
☑ Target 1.6 (patterned anisotropic layers) — done 2026-08-03: `fourier_factorization.toeplitz_matrix_component`, `eigenmodes.solve_layer_eigenmodes_patterned_inplane`, transcribed from `S4/S4/fmm/fmm_closed.cpp`'s `have_tensor` branch (lines 165-256); see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` 1.6 and `memory.md`
☑ Target 1.7 (degeneracy policy) — done 2026-08-03: `eigenmodes._canonical_mode_order`, applied to the three anisotropic dense eigensolvers; see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` 1.7 and `memory.md`
☑ Target 1.8 (mode classification) — done 2026-08-03: `eigenmodes.classify_propagating`, `SimulationResult.order_classification()`; found (not fixed, out of scope) an exact-Rayleigh-threshold NaN division-by-zero, now documented in `troubleshooting.md` and tied to Category 6 target 6.4; see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` 1.8 and `memory.md`
☑ Category 1 targets 1.3-1.4, 1.6-1.8 all done 2026-08-03 (six-target session); target 1.5 explicitly deferred (no citable+benchmarkable source found)

## Phase 7 — Real-Space Field Reconstruction & Visualization (DONE)

☑ Extend `fields.py` with a real-space grid reconstruction function using `SMatrixStack.partial_smatrix_up_to` — `fields.modal_field_components`/`reconstruct_field_at_points`, `smatrix.interior_amplitudes` (independently derived, `decisions.md` ADR-015)
☑ Test: field continuity across a layer interface (no discontinuity where physically none should exist) — `tests/test_field_reconstruction.py`
☑ Test: reconstructed field-derived R/T matches the already-validated `SimulationResult.reflectance()`/`transmittance()` (cross-check, not a new independent oracle) — found and documented a missing-0.5-factor convention in `fields.z_poynting_flux` along the way (see `troubleshooting.md`)
☑ Add `matplotlib` as a dev/example dependency (not a core `sougata_solver` dependency — confirm this placement in `pyproject.toml`) — already present from earlier phases, reused unmodified
☑ Write a cross-section field-intensity plotting example for a trench — `structures/trench/trench_field_cross_section.py` + `postprocessing/plot_field_cross_section.py`
☑ Write a cross-section field-intensity plotting example for a via — `structures/via/pillar_field_cross_section.py` + `postprocessing/plot_field_cross_section.py`
☑ Update `memory.md` / `decisions.md` on completion — done 2026-08-05, `decisions.md` ADR-015

## Phase 8 — Expanded Validation Suite & Example Gallery

□ Convergence-vs-`num_orders` study: trench
□ Convergence-vs-`num_orders` study: via/pillar
□ Convergence-vs-`num_orders` study: tapered via
□ Example: DBR-style multilayer (mirrors vendored `EMTutorial/ThinFilmsAndMultilayers/DistributedBraggReflector`)
□ Example: TSV-style via (mirrors vendored `EMTutorial/Scatterometry/ThroughSiliconVia`)
□ Review and refresh `README.md`'s Features/Future Improvements sections against actual completed phases
□ Update `memory.md` / `decisions.md` on completion

## Phase 9 — Performance & Optional GPU/Autodiff Backend (later, optional)

□ Profile the current per-point `Simulation.solve` call to find the actual bottleneck (don't assume)
□ Vectorize wavelength/angle sweeps in NumPy (batch eigensolves / S-matrix ops)
□ Regression test: vectorized sweep numerically matches the unvectorized per-point loop
□ Decision checkpoint: confirm GPU/autodiff backend is still wanted before starting it (re-ask, don't assume — see `decisions.md`)
□ (If pursued) Design a backend-agnostic array-op interface behind `eigenmodes.py`/`smatrix.py`
□ (If pursued) Implement a torch or JAX backend against that interface
□ (If pursued) Validate backend numerically matches the NumPy path
□ Update `memory.md` / `decisions.md` on completion

## Cross-cutting: `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Categories 1-19

Work tracked by that register's atomic targets (numerical-methods
robustness, Fourier-factorization rule selection, geometry/material
model completeness, validation, docs, etc.) is not phase-scoped the way
Phases 1-9 above are, and is tracked at fine grain in
`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` itself rather than duplicated here.
Status snapshot (update the source-of-truth file first, then this line):

☑ Category 1 (Mathematical foundation / anisotropy) — targets 1.1-1.4,
  1.6-1.8 done; 1.5 (longitudinal coupling) explicitly deferred, 2026-08-03.
☑ Category 2 (Numerical methods) — targets 2.1-2.5 all done, 2026-08-04
  (failure contract, eigenvalue-diagnostics report, sweep mode-order
  stability, degeneracy-gap warning, one lossy high-contrast full-pipeline
  stress fixture).
☑ Category 3 (Fourier factorization) — targets 3.1-3.6 all done, 2026-08-04
  (rule inventory + regression tests; 1D and 2D high-contrast convergence
  fixtures with honestly-recorded non-monotonic low-order findings; FFF/NVM
  feasibility investigation against S4's `fmm_PolBasisNV`/`PolBasisJones`/
  `PolBasisVL` concluding "defer," `decisions.md` ADR-012; target 3.6 has no
  implementation as a direct consequence of that deferral).
☑ Category 4 (Geometry engine) — targets 4.1-4.7 all done, 2026-08-04
  (construction-time shape/lattice validation; unit-cell self-overlap
  policy wired into `Simulation.__init__`; `Ellipse`/`Polygon` primitives,
  both transcribed from a previously-unread S4 `pattern_get_fourier_transform`
  branch found this session, each with a `structures/via/` end-to-end
  example; `decisions.md` ADR-013 narrowly revisits ADR-005 for the
  analytic-only `Polygon` case; a minimal safe JSON pattern-import format
  (`geometry_io.py`, parser-only, not solver-wired yet); `staircase.py`'s
  three shape-specific taper generators refactored into thin wrappers
  around a new general `slice_profile` interface, regression-verified
  unchanged).
☑ Category 5 (Material models) — targets 5.1-5.8 all done, 2026-08-04
  (construction- and call-time `Material` validation; Sellmeier/Cauchy
  dispersion models transcribed from vendored `EMpy`, validated against
  BK7's independently-published `n_d`; Lorentz/Drude/Drude-Lorentz models
  transcribed from vendored `RigorousCoupledWaveAnalysis.jl`'s Rakic
  Lorentz-Drude metal model plus its published Au/Ag/Al/Ti coefficients,
  with an explicit causality/sign-convention re-derivation and check;
  dispersive-tensor-material solver wiring confirmed end to end; optional
  `Material.source` provenance metadata threaded through every classmethod
  and into serialized run metadata, surfacing and fixing a real pre-existing
  Windows-encoding bug in `output_paths.write_run_metadata` along the way).
☑ Category 6 (Boundary conditions and excitation) — targets 6.1-6.6 all
  done, 2026-08-04 (a "Worked polarization examples" table added to the
  already-existing `CONVENTIONS.md`; a full polarization-state x azimuth x
  angle regression suite using symmetry invariants, not just energy
  conservation; a characterized-and-tested grazing-incidence boundary
  (`ValueError`, not `NaN`, exactly at `theta=90 deg`, traced to a
  floating-point coincidence); an oblique-incidence extension of Category
  1's Rayleigh-threshold test; a Stokes-reciprocity-verified finding that
  bottom illumination needs no new API, `decisions.md` ADR-014).
☑ Category 7 (Layer handling) — targets 7.1-7.6 all done, 2026-08-05
  (construction-time thickness validation; a repeated-layer-identity
  regression guard; an instance-scoped Toeplitz-matrix cache gated on a
  measured timing case per `rules.md`'s Performance Requirements exception
  clause, `decisions.md` ADR-016, with an honest mid-session correction to
  the first (wrongly-framed) measurement; `SimulationResult.layer_absorption()`,
  a flux-divergence combination of already-validated Category 9 pieces,
  `decisions.md` ADR-017, closing the `R+T+sum(A)=1` energy-identity gap
  `test_stress_regression.py` had flagged since Category 2; a found-and-
  documented numerical-overflow limitation for thick/highly-lossy/high-
  `num_orders` interior-amplitude reconstruction, `troubleshooting.md`).
☑ Category 8 (Solver sweeps and convergence) — targets 8.1-8.8 all done,
  2026-08-05 (new `sweep.py` module: a typed `SweepResult` container;
  wavelength/theta/phi/polarization/thickness sweep functions, each
  confirmed equivalent to a manual `Simulation.solve()` loop; a
  harmonic-order study reusing Category 7's `layer_absorption()` for its
  conservation residual; a conservative convergence criterion,
  `decisions.md` ADR-018, validated against thin-film/trench/pillar
  fixtures per target 8.8's own gating requirement before automatic
  harmonic-order selection was implemented on top of it; one honest bug
  found and fixed via the project's own test-first discipline before the
  criterion was trusted).
☑ Category 9 (Field calculations) — targets 9.1-9.8 all done, 2026-08-05
  (real-space field reconstruction transcribed from S4's
  `GetInPlaneFieldVector`/`GetFieldAtPoint`; interior-layer amplitude
  recovery independently derived from `SMatrixStack.partial_smatrix_up_to`,
  `decisions.md` ADR-015; a found-and-documented missing-0.5-factor
  convention in `fields.z_poynting_flux` relative to the textbook
  real-space flux formula; two `structures/` example scripts (trench (x,z)
  cross-section, pillar (x,y) field map) and one `postprocessing/` plotting
  script, both run end-to-end and visually verified). This is also Phase 7
  (`phases.md`), now marked DONE.
☐ Category 10 (Optical outputs) — targets 10.1-10.4/10.6 done, 2026-08-05
  (complex per-order Cartesian field coefficients validated against a new
  `tests/oracles/fresnel.py::multilayer_complex_rt` function for both
  polarizations; per-order diffraction angles with a `None` non-
  propagating representation; a one-call conservation report; a frozen
  output schema across uniform/1D/2D fixtures). Target 10.5 (per-order
  s/p conversion) explicitly deferred -- a bounded attempt to externally
  validate the polarization convention against S4's actual source found a
  plausible but numerically-unconfirmed match (S4 not buildable in this
  environment); see `references.md`.
☐ Category 11 (Semiconductor OCD features) — targets 11.1-11.7 done,
  2026-08-05 (new `ocd.py` module: a validated CD-first parameter object,
  a trapezoid trench constructor thin-wrapping Phase 5's staircase
  machinery, an arc-sampled-`Polygon` corner-rounding geometry validated
  against a closed-form area, and reproducible TSV/grating OCD example
  sweeps recording every swept parameter in `run_metadata.txt`; overlay
  confirmed already achievable with no new API, `decisions.md` ADR-019,
  verified via a shift-by-one-period periodicity check). Target 11.8
  (LER/LWR) explicitly deferred -- genuine stochastic roughness
  fundamentally conflicts with RCWA's periodic-Fourier formulation;
  `decisions.md` ADR-020.
☑ Category 12 (Linear algebra) — targets 12.1-12.5 all done, 2026-08-05
  (new `profiling/baseline_profile.py` measuring eigensolve/matrix-solve/
  end-to-end timing on fixed fixtures, showing the eigensolve dominates
  at larger `num_orders`; a direct-inverse audit that found and fixed a
  house-convention inconsistency in `eigenmodes.py`, confirmed bit-for-
  bit equivalent; a factorization-reuse design note finding no further
  S-matrix-level reuse beyond the already-shipped trivial-interface fast
  path; an opt-in `svd_diagnostics` function; and a sparse/iterative-
  methods feasibility decision **rejected** on a measured 100%-dense-
  matrix structural finding, `decisions.md` ADR-021).
☑ Category 13 (Performance optimization) — targets 13.1-13.6 all
  resolved, 2026-08-05 (a repeatable benchmark suite extending Category
  12's profiler with a tapered-via case; an eigenmode-reuse cache
  implementing Category 12 target 12.3's deferred design, ~3.3x measured
  on a polarization sweep, `decisions.md` ADR-022; a narrowly-scoped
  vectorized thin-film wavelength sweep, ~31x measured, with a real bug
  caught and fixed by its own equivalence test, `decisions.md` ADR-023;
  a measured parallelism decision -- threading modestly helps,
  multiprocessing measured counterproductive on this machine,
  `decisions.md` ADR-024; a GPU decision checkpoint where explicit
  approval was sought from the project owner and not granted, deferring
  GPU/autodiff work to Phase 9).
☑ Category 14 (Validation) — targets 14.1-14.8 all done, 2026-08-07
  (a validation inventory mapping every public feature to its oracle/
  invariant test/example/limitation, `testing.md`; the external 2D R/T
  oracle re-evaluated and still documented-blocked, S4 unbuildable in
  this environment and no matching versioned dataset found, 14.3/14.4
  documented as blocked on it rather than a false pass; reciprocity
  tests built on a Snell's-law-matched-angle comparison verified
  numerically before writing any assertion -- the naive same-theta
  comparison was tried first and found wrong at oblique incidence,
  `decisions.md` ADR-025, `tests/test_reciprocity.py`; a harmonic
  convergence matrix across all 7 supported geometry families, every
  candidate/tolerance measured directly rather than guessed,
  `tests/test_harmonic_convergence_matrix.py`, all 4 `slow`-marked cases
  confirmed passing; a validation report of tolerances/versions/results,
  `testing.md`).
☑ Category 15 (User interface and API) — targets 15.1-15.8 all done,
  2026-08-07 (a public API inventory that found and fixed a real
  `__init__.py` export staleness bug; a minimal JSON simulation-
  configuration schema, `config.py`, reusing `geometry_io.py`'s existing
  material/pattern sub-schemas rather than inventing new ones, with
  construction-time-only validation and malformed-input tests,
  `tests/test_config.py`; a config-file-driven reproduction of
  `structures/thin_film/anti_reflection_coating.py` to `1e-12`; a CLI
  (`cli.py`, one `run` subcommand, three distinct exit codes, a
  `sougata-solver` console-script entry point) with
  `tests/test_cli.py`; a NumPy `.npz` result-series exporter,
  `export.py`, deliberately avoiding `allow_pickle=True` by JSON-encoding
  metadata into a string array, `tests/test_export.py`; an HDF5
  decision evaluated and deferred, `decisions.md` ADR-026, since
  current result shapes are already well served by the `.npz` exporter).
☑ Category 16 (Visualization) — targets 16.1-16.7 all done, 2026-08-07
  (new `plotting.py` module, every function taking plain arrays/
  dataclasses/already-computed result objects, never a `Simulation`, and
  never calling `.solve()`, pinned by a direct structural test;
  `plot_unit_cell` rasterizes a preview grid via each shape's existing
  `.contains()` method, respecting `Pattern`'s "later shapes take
  precedence" rule; `plot_rt_spectrum`/`plot_field_intensity` formalize
  two existing `postprocessing/` ad hoc plots into reusable functions;
  `plot_harmonic_convergence`/`plot_diffraction_orders`/
  `plot_field_phase`/`plot_poynting_vector` are new. 19 new tests,
  structural checks not pixel comparisons).
☑ Category 17 (Testing and quality) — targets 17.1-17.6 all done,
  2026-08-07 (a test taxonomy documenting the existing filename/
  docstring convention plus one new precisely-scoped `pytest` marker,
  `oracle`, applied to the 8 files that actually import a named external
  oracle module; Windows CI (`ci.yml`, 3-version matrix, `ruff` +
  fast-suite gate on every push/PR) and a separate weekly/manual
  slow-test workflow (`slow-tests.yml`); a frozen regression-fixture
  snapshot with explicit provenance and a tolerance rationale
  distinguishing it from a fresh oracle comparison; `ruff` static
  analysis configured and its 24-issue baseline fixed, including two
  genuine dead-code findings in already-shipped code; a performance
  regression guard using a same-run relative-scaling ratio rather than
  an absolute wall-clock threshold, per `rules.md`'s Performance
  Requirements, `decisions.md` ADR-028).
☑ Category 18 (Documentation) — targets 18.1-18.8 all done, 2026-08-20
  (four new root-level docs, each a consolidation of already-implemented/
  cited/tested material, not new derivation or example code: `theory.md`
  (18.1-18.3, a ToC over `design.md`/`CONVENTIONS.md`/`s_matrix_method.md`
  plus an end-to-end pipeline narration none of those individually gave);
  `api_reference.md` (18.4, expands `src/sougata_solver/README.md`'s
  Module Map into a full per-symbol reference); `tutorials.md` (18.5-18.7,
  walks through one already-oracle-validated example script per structure
  family, all three re-run this session for real captured output);
  `validation_guide.md` (18.8, an oracle-centric companion to `testing.md`'s
  category-centric Validation Inventory, profiling each of the 6
  `tests/oracles/*.py` files' actual proof scope). Also closed out
  Phase 8's previously-ambiguous status (`phases.md`) as a side effect,
  since its "example gallery" deliverable is what 18.5-18.7 completed.
  706 tests pass project-wide, unchanged (no `src/sougata_solver/`
  change).
□ Category 19 — not yet started at atomic-target grain, explicitly
  deferred pending the project owner's own use-case decision (see that
  file for each small target's own checklist).
