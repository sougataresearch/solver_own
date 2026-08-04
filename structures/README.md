# `structures/` — Run These First

Runnable scripts: each one builds a `Lattice`/`Layer` stack/`Material`
set, runs `Simulation.solve()` over a wavelength (and/or angle) sweep, and
saves raw results (never plots — that's [`postprocessing/`](../postprocessing/)'s
job, see `decisions.md` ADR-009/010). One subfolder per geometry type.

```bash
python structures/thin_film/sio2_on_si_thin_film.py
```

Every run writes its own timestamped folder under `outputs/YYYY_MM_DD/HH_MM_SS_<script-name>/`
(gitignored) containing the raw CSV and a `run_metadata.txt` recording which
script produced it and its key parameters — see
[`src/sougata_solver/output_paths.py`](../src/sougata_solver/output_paths.py).
Re-running the same script with different settings never overwrites a
previous run.

## `thin_film/` (Phase 1, done)

Uniform, laterally-infinite multilayer stacks — no in-plane pattern, so
`num_orders`/`Lattice` are required by `Simulation` but unused (see the
comment in each script).

| Script | Purpose |
|---|---|
| [`sio2_on_si_thin_film.py`](thin_film/sio2_on_si_thin_film.py) | SiO2-on-Si, wavelength sweep, R/T/A to CSV — **copy this one to start a new stack** |
| [`custom_multistack.py`](thin_film/custom_multistack.py) | Reusable N-layer stack template |
| [`anti_reflection_coating.py`](thin_film/anti_reflection_coating.py) | Single- or multi-layer AR coating example |
| [`sio2_on_si_ellipsometry_run.py`](thin_film/sio2_on_si_ellipsometry_run.py) | Saves raw field data (not just R/T) for `postprocessing/jones_mueller_ellipsometry.py` to consume |

### Editing a script for your own structure

Each script has numbered `# EDIT (n):` comment blocks — material CSV paths
(or placeholder constants if the CSVs aren't found), layer thicknesses,
incident angle/azimuth/polarization (`s_amplitude`/`p_amplitude`, complex —
their ratio sets linear/circular/elliptical polarization), wavelength sweep,
and output path. No other part of the script should normally need touching.

## `trench/` (Phase 3 lamellar grating, Phase 5 tapered sidewall — done)

1D-periodic patterned layers (`Lattice1D`/`Slab`).

| Script | Purpose |
|---|---|
| [`trench_grating.py`](trench/trench_grating.py) | Uniform (untapered) lamellar grating, wavelength sweep |
| [`trench_grating_ellipsometry_run.py`](trench/trench_grating_ellipsometry_run.py) | Saves raw field data for ellipsometry postprocessing |
| [`tapered_trench.py`](trench/tapered_trench.py) | Linearly-tapered ridge (staircase-discretized), `num_slices` convergence sweep; FDTD-style `TCD`/`BCD`/`SPACING`/`PERIOD` naming for the ridge geometry |

## `via/` (Phase 4a/4b 2D patterned layers, Phase 5 tapered sidewall — done)

2D-periodic patterned layers (`Lattice`/`Circle`/`Rectangle`).

| Script | Purpose |
|---|---|
| [`pillar_array.py`](via/pillar_array.py) | Uniform (untapered) circular Si pillar, wavelength sweep |
| [`via_array.py`](via/via_array.py) | Uniform (untapered) circular via/hole, wavelength sweep |
| [`tapered_via.py`](via/tapered_via.py) | Linearly-tapered circular via (staircase-discretized), `num_slices` convergence sweep; FDTD-style `TCD`/`BCD`/`SPACING`/`PERIOD` naming (via diameters) |
| [`tapered_pillar.py`](via/tapered_pillar.py) | Linearly-tapered square pillar (staircase-discretized, `Rectangle` with equal x/y halfwidths), same FDTD-style naming (pillar side length) |
| [`elliptical_pillar.py`](via/elliptical_pillar.py) | Uniform (untapered) elliptical Si pillar, wavelength sweep — Category 4 target 4.3's end-to-end `Ellipse` example |
| [`triangular_pillar.py`](via/triangular_pillar.py) | Uniform (untapered) triangular Si pillar (`Polygon`), wavelength sweep — Category 4 targets 4.4/4.5's end-to-end `Polygon` example |

All three tapered scripts (`tapered_trench.py`, `tapered_via.py`,
`tapered_pillar.py`) build on `src/sougata_solver/staircase.py`'s
`staircase_slab_layers`/`staircase_circle_layers`/`staircase_rectangle_layers`
generators (Phase 5, `decisions.md` ADR-004) and use `TCD`/`BCD` (top/bottom
critical dimension) plus `SPACING` naming to match the equivalent FDTD
(Lumerical) grating-structure-group parametrization.

## Doubts already resolved for this folder (see `progress_log.md` 2026-07-19)

If you see dense, fast oscillation in an R/T-vs-wavelength plot for a
thick layer (e.g. the 12 um Si substrate in `sio2_on_si_thin_film.py`),
that's real Fabry-Perot interference (`Delta_lambda ~ lambda^2/(2 n t)`),
not a solver bug — increase the wavelength sample count if it looks
aliased/undersampled rather than assuming something's wrong with the
boundaries. RCWA has no mesh and no PML-equivalent to worry about; the
incidence/transmission media are exact semi-infinite half-spaces.
