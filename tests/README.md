# `tests/` — Test Suite

`pytest`-based, 569 tests as of `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
Category 8 (559 fast + 10 `slow`-marked convergence/benchmark studies). Run
the fast suite with:

```bash
pytest                # fast suite (excludes slow)
pytest -m slow        # convergence/benchmark studies only (several minutes)
```

See [`testing.md`](../testing.md) for the full testing strategy and
[`rules.md`](../rules.md) for the project's validation discipline: nothing
in `src/sougata_solver/` is trusted until it agrees with an independent
oracle — never a self-consistency check against the same code path.

## Layout

| File | What it validates |
|---|---|
| [`conftest.py`](conftest.py) | Shared fixtures |
| [`test_analytic_fresnel.py`](test_analytic_fresnel.py) | Phase 1 uniform-multilayer R/T against `oracles/fresnel.py` (from-scratch analytic Fresnel/TMM), across wavelength/angle/polarization, lossless and absorbing structures |
| [`test_thin_film_empy_cross_check.py`](test_thin_film_empy_cross_check.py) | Same Phase 1 solve, cross-checked against `oracles/empy_tmm.py` (transcribed from the vendored EMpy library) — a **second, independent** oracle so Phase 1 isn't validated against only one reference |
| [`test_polarimetry.py`](test_polarimetry.py) | Jones/Mueller matrix construction |
| [`test_fourier_factorization.py`](test_fourier_factorization.py) | Phase 2 `pattern_epsilon_hat`/`toeplitz_matrix` against **two independent** references: a from-scratch rasterize-and-sum, and a from-scratch FFT-of-rasterized-mask reproduction of the vendored `RigorousCoupledWaveAnalysis.jl`/`convmat2D.py` algorithm |
| [`test_1d_grating.py`](test_1d_grating.py) | Phase 3 1D-grating solve against `oracles/rcwa_1d_gaylord.py`, plus energy conservation and convergence-rate invariants |
| [`test_2d_pillar.py`](test_2d_pillar.py) | Phase 4a 2D-patterned solve against `oracles/rcwa_2djl_eigenvalues.py` (independent RCWA.jl-derived eigenoperator), plus reduction/energy-conservation checks |
| [`test_2d_pillar_stress.py`](test_2d_pillar_stress.py) | Phase 4b near-degenerate/ill-conditioned stress cases and `ILL_CONDITIONED_THRESHOLD` logging |
| [`test_staircase.py`](test_staircase.py) | Phase 5 tapered-sidewall staircase discretization: zero-taper regression, energy conservation, convergence-vs-`num_slices` |
| [`test_anisotropic_uniform.py`](test_anisotropic_uniform.py) | Category 1 target 1.3 (uniform diagonal-tensor layers): closed-form birefringence benchmark, Fresnel-oracle per axis, isotropic reduction |
| [`test_anisotropic_inplane.py`](test_anisotropic_inplane.py) | Category 1 target 1.4 (uniform in-plane-coupled layers) against `oracles/rcwa_anisotropic_inplane_jl.py` |
| [`test_anisotropic_patterned.py`](test_anisotropic_patterned.py) | Category 1 target 1.6 (patterned anisotropic layers): reduction to isotropic/uniform-tensor cases, energy conservation |
| [`test_anisotropic_degeneracy.py`](test_anisotropic_degeneracy.py) | Category 1 target 1.7: deterministic mode ordering, repeated-solve determinism |
| [`test_mode_classification.py`](test_mode_classification.py) | Category 1 target 1.8: propagating/evanescent classification against the analytic Rayleigh threshold |
| [`test_failure_contract.py`](test_failure_contract.py) | Category 2 target 2.1: one test per documented `ValueError`/`NotImplementedError`/`LinAlgError` condition in `design.md`'s Failure Contract table |
| [`test_eigenvalue_diagnostics.py`](test_eigenvalue_diagnostics.py) | Category 2 target 2.2: `LayerEigenmodes.diagnostics` (`EigenmodeDiagnostics`) fields match independent recomputation, no change to solve results |
| [`test_sweep_mode_matching.py`](test_sweep_mode_matching.py) | Category 2 target 2.3: canonical mode ordering doesn't arbitrarily permute across a small wavelength sweep, for the three anisotropic dense solvers |
| [`test_degeneracy_warning.py`](test_degeneracy_warning.py) | Category 2 target 2.4: `eigenmodes.DEGENERATE_GAP_THRESHOLD` warning fires/doesn't fire correctly |
| [`test_stress_regression.py`](test_stress_regression.py) | Category 2 target 2.5: one lossy high-contrast fixture through the full `Simulation.solve()` pipeline (passivity check, since layer-wise absorption isn't implemented yet) |
| [`test_fourier_factorization_rules.py`](test_fourier_factorization_rules.py) | Category 3 target 3.1: pins `design.md`'s Fourier-factorization rule inventory table against actual solver behavior |
| [`test_fourier_convergence.py`](test_fourier_convergence.py) | Category 3 targets 3.2/3.3: fixed high-contrast 1D lamellar and 2D pillar convergence fixtures, `slow`-marked, with recorded (not just asserted) convergence tables |
| [`test_geometry_validation.py`](test_geometry_validation.py) | Category 4 target 4.1: construction-time validation for `Lattice`/`Lattice1D`/`Circle`/`Rectangle`/`Slab` |
| [`test_unit_cell_bounds.py`](test_unit_cell_bounds.py) | Category 4 target 4.2: edge-crossing shapes match a from-scratch periodic-tiling raster reference; self-overlap-across-periodic-images rejection |
| [`test_ellipse.py`](test_ellipse.py) | Category 4 target 4.3: `Ellipse` DC/area, a from-scratch rasterized cross-check, and reduction to `Circle` |
| [`test_polygon.py`](test_polygon.py) | Category 4 target 4.5: `Polygon` reduction to `Rectangle` (square case) and a from-scratch rasterized cross-check for a triangle and a non-convex L-shape |
| [`test_geometry_io.py`](test_geometry_io.py) | Category 4 target 4.6: `geometry_io`'s minimal JSON `Pattern`-import format — parsing and validation only |
| [`test_profile_slicing.py`](test_profile_slicing.py) | Category 4 target 4.7: `staircase.slice_profile`'s general geometry-to-layer-slices interface, independent of `Simulation.solve` |
| [`test_material_validation.py`](test_material_validation.py) | Category 5 target 5.1: `Material` construction- and call-time validation (tensor shape, finite values, callback output) |
| [`test_dispersion_models.py`](test_dispersion_models.py) | Category 5 targets 5.2-5.6: Sellmeier/Cauchy/Lorentz/Drude/Drude-Lorentz dispersion models, including a causality/sign-convention check and Rakić et al. (1998)'s published Au/Ag/Al/Ti coefficients |
| [`test_tensor_material_wiring.py`](test_tensor_material_wiring.py) | Category 5 target 5.7: a dispersive tensor material solving end to end through Category 1's uniform-diagonal and patterned-anisotropic eigensolvers |
| [`test_material_provenance.py`](test_material_provenance.py) | Category 5 target 5.8: optional `Material.source` citation metadata, forwarded by every `from_*` classmethod and threaded into serialized `run_metadata.txt` output |
| [`test_polarization_states.py`](test_polarization_states.py) | Category 6 targets 6.2/6.3: polarization-state x azimuth x angle regression suite using symmetry invariants (rotational symmetry at normal incidence, azimuthal invariance at oblique incidence) plus an energy-conservation sweep |
| [`test_grazing_incidence.py`](test_grazing_incidence.py) | Category 6 target 6.4: characterized grazing-incidence boundary — finite/energy-conserving up to `89.999deg`, `ValueError` (not `NaN`) exactly at `90deg` |
| [`test_oblique_rayleigh_threshold.py`](test_oblique_rayleigh_threshold.py) | Category 6 target 6.5: oblique-incidence Rayleigh-threshold order-classification crossing |
| [`test_bottom_incidence.py`](test_bottom_incidence.py) | Category 6 target 6.6: bottom (reverse-side) illumination via the existing `Simulation` constructor, verified by Stokes transmittance reciprocity |
| [`test_field_reconstruction.py`](test_field_reconstruction.py) | Category 9 targets 9.1-9.8 (Phase 7): real-space field reconstruction against the analytic plane wave, transversality, interior-amplitude self-consistency, interface field continuity, 1D periodicity, and flux-matches-R/T |
| [`test_layer_validation.py`](test_layer_validation.py) | Category 7 target 7.1: `Layer`/`LayerStack` construction-time invariants (finite-layer thickness, the semi-infinite half-space sentinel, patterned-layer background material) |
| [`test_layer_repetition.py`](test_layer_repetition.py) | Category 7 target 7.2: equivalent repeated-layer representations (split thickness, reused vs. re-constructed identical `Pattern`) give identical R/T |
| [`test_layer_cache.py`](test_layer_cache.py) | Category 7 target 7.4: the Toeplitz-matrix cache's equivalence to forced-uncached recomputation, cache-hit call counting, and angle-sweep cache reuse |
| [`test_layer_absorption.py`](test_layer_absorption.py) | Category 7 targets 7.5/7.6: `SimulationResult.layer_absorption()` against the `R+T+sum(A)=1` energy-balance identity, plus a regression guard on a found numerical-overflow limitation for thick/highly-lossy/high-`num_orders` layers |
| [`test_sweep.py`](test_sweep.py) | Category 8 targets 8.1-8.5: `sweep.SweepResult` and the wavelength/angle/polarization/thickness sweep functions, each confirmed equivalent to a manual per-point `solve()` loop |
| [`test_harmonic_convergence.py`](test_harmonic_convergence.py) | Category 8 targets 8.6-8.8: `harmonic_study`, `find_convergence_index` (validated against thin-film/trench/pillar fixtures, including Category 3's non-monotonic pillar wobble), and `auto_select_num_orders` |

## `oracles/`

Hand-transcribed reference implementations used *only* as independent
ground truth in tests — never imported by `src/sougata_solver/` itself.

| File | Source |
|---|---|
| [`fresnel.py`](oracles/fresnel.py) | From-scratch analytic Fresnel/TMM, derived independently (not transcribed from a vendored repo) |
| [`empy_tmm.py`](oracles/empy_tmm.py) | Transcribed by hand from the vendored `EMpy` reference library (`EMpy` itself is never imported at runtime — see `decisions.md`) |
| [`rcwa_1d_gaylord.py`](oracles/rcwa_1d_gaylord.py) | Transcribed from `Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples` (Moharam/Gaylord-style 1D binary grating) |
| [`rcwa_2djl_eigenvalues.py`](oracles/rcwa_2djl_eigenvalues.py) | Transcribed from `RigorousCoupledWaveAnalysis.jl`'s isotropic patterned-layer eigenoperator (Julia not installed here; hand-transcribed, not run) |
| [`rcwa_2d_pillar.py`](oracles/rcwa_2d_pillar.py) | Documents the still-open external 2D R/T oracle gap (Phase 4a/4b) rather than fabricating a benchmark |
| [`rcwa_anisotropic_inplane_jl.py`](oracles/rcwa_anisotropic_inplane_jl.py) | Transcribed from `RigorousCoupledWaveAnalysis.jl`'s uniform anisotropic-layer eigenoperator (Category 1 target 1.4) |

## Adding a test for new solver capability

Per `rules.md`: every new formula needs an oracle that is (a) independent
of the code path under test, and (b) not itself written by paraphrasing the
implementation. Acceptable oracles, in the order this project has used them
so far: an analytic closed-form solution derived by hand, a from-scratch
numerical reference (rasterize-and-sum, brute-force FFT, etc.), a
line-numbered transcription from a vendored reference repo (`S4`, `EMpy`,
`RigorousCoupledWaveAnalysis.jl` — never imported, only read from and
hand-transcribed), or a published benchmark table (e.g. Moharam & Gaylord
1995 for Phase 3's grating diffraction efficiencies). Never validate a
formula only against itself restated a different way.
