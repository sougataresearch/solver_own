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

## `oracles/`

Hand-transcribed reference implementations used *only* as independent
ground truth in tests — never imported by `src/sougata_solver/` itself.

| File | Source |
|---|---|
| [`fresnel.py`](oracles/fresnel.py) | From-scratch analytic Fresnel/TMM, derived independently (not transcribed from a vendored repo) |
| [`empy_tmm.py`](oracles/empy_tmm.py) | Transcribed by hand from the vendored `EMpy` reference library (`EMpy` itself is never imported at runtime — see `decisions.md`) |

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
