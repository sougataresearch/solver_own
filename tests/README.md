# `tests/` — Test Suite

`pytest`-based. Run with:

```bash
pytest
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
