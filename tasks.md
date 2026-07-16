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

## Phase 2 — Fourier-Factorization Core

□ Add `.flake8`/`ruff` config and `mypy.ini` to `sougata_solver/` (rules.md gap, do before new modules land)
□ Create `src/sougata_solver/fourier_factorization.py`
□ Implement `pattern_epsilon_hat(pattern, g_vectors, lattice)` (direct, sums shape contributions with containment-tree subtraction)
□ Implement the inverse-rule variant (`1/eps` per shape, same summation) for `epsilon_inv_hat`
□ Implement `toeplitz_matrix(eps_hat_lookup, g_indices)`
□ Write numerical FFT-of-rasterized-mask reference for a `Circle` pattern
□ Write numerical FFT-of-rasterized-mask reference for a `Rectangle` pattern
□ Test: analytic Toeplitz entries match FFT reference within tolerance, for both `Circle` and `Rectangle`
□ Test: DC term of `epsilon_hat` equals area-weighted average permittivity (closed-form sanity check)
□ Update `memory.md` / `decisions.md` on completion

## Phase 3 — 1D-Periodic Lamellar Gratings (Trench)

□ Add `Lattice1D(period)` to `geometry.py`
□ Add `Slab`/`Line` 1D shape with analytic (`sinc`) Fourier transform
□ Add `truncate_fourier_orders_1d` to `fourier_basis.py`
□ Implement `solve_layer_eigenmodes_1d` (TE path, using `epsilon_hat`)
□ Implement `solve_layer_eigenmodes_1d` (TM path, using `epsilon_inv_hat`)
□ Add `Lattice1D` dispatch branch in `simulation.py`
□ Source and transcribe a published 1D binary-grating benchmark table (Moharam & Gaylord 1995 or equivalent) into a test oracle
□ Test: TE diffraction efficiencies match benchmark
□ Test: TM diffraction efficiencies match benchmark
□ Test: normal-incidence limit recovers Phase 1's uniform-layer Fresnel result when line/space contrast is set to zero (continuity sanity check)
□ Write `structures/trench/trench_grating.py`
□ Update `memory.md` / `decisions.md` on completion

## Phase 4 — 2D-Periodic Patterned Layers (Via, Pillar)

□ Implement `solve_layer_eigenmodes_patterned` (general non-uniform eigenproblem, transcribed from `S4/S4/rcwa.cpp::SolveLayerEigensystem` lines 794-827)
□ Handle near-degenerate eigenvalue edge cases (document the approach in the function's docstring)
□ Remove the `NotImplementedError` at `simulation.py:98`, wire in Phase 2 Toeplitz construction + this solver
□ Determine whether S4 is buildable/runnable in this environment for a subprocess cross-check oracle
□ If S4 is usable: write an S4-driven oracle test for a simple pillar array
□ If S4 is not usable: source a published 2D benchmark instead, and explicitly document why S4 wasn't used (per `rules.md` AI rule 5 — never fabricate a match)
□ Test: 2D patterned-layer R/T matches the chosen oracle
□ Test: patterned-layer solve reduces to the uniform-layer result when the pattern's shape material equals the background (degenerate-pattern sanity check)
□ Write `structures/via/pillar_array.py`
□ Write `structures/via/via_array.py`
□ Update `memory.md` / `decisions.md` on completion

## Phase 5 — Tapered / Sloped Sidewalls (Via, Trench)

□ Design the staircase-layer-stack generator's API (inputs: top size, bottom size, thickness, slice count `N`; output: `list[Layer]`)
□ Implement the generator for `Rectangle`/`Circle` (via)
□ Implement the generator for `Slab` (trench)
□ Write a convergence-vs-`N` test for a tapered via (mark `slow`)
□ Write a convergence-vs-`N` test for a tapered trench (mark `slow`)
□ Write an example script sweeping `N` and plotting/printing R/T convergence
□ Update `memory.md` / `decisions.md` on completion

## Phase 6 — Anisotropic Materials

□ Generalize Phase 4's eigensolver to accept a full 3×3 tensor `Epsilon2`
□ Remove `simulation.py`'s uniform-anisotropic `NotImplementedError`
□ Source a birefringent-material closed-form benchmark (e.g. uniaxial waveplate at normal incidence)
□ Test: anisotropic solve matches the benchmark
□ Test: isotropic-tensor special case reduces to Phase 1's uniform-isotropic result (regression guard)
□ Update `memory.md` / `decisions.md` on completion

## Phase 7 — Real-Space Field Reconstruction & Visualization

□ Extend `fields.py` with a real-space grid reconstruction function using `SMatrixStack.partial_smatrix_up_to`
□ Test: field continuity across a layer interface (no discontinuity where physically none should exist)
□ Test: reconstructed field-derived R/T matches the already-validated `SimulationResult.reflectance()`/`transmittance()` (cross-check, not a new independent oracle)
□ Add `matplotlib` as a dev/example dependency (not a core `sougata_solver` dependency — confirm this placement in `pyproject.toml`)
□ Write a cross-section field-intensity plotting example for a trench
□ Write a cross-section field-intensity plotting example for a via
□ Update `memory.md` / `decisions.md` on completion

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
