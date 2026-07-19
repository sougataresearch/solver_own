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
| [`custom_material_from_nk_data.py`](thin_film/custom_material_from_nk_data.py) | Building a `Material` from your own `n,k` data (not refractiveindex.info CSV format) |
| [`sio2_on_si_ellipsometry_run.py`](thin_film/sio2_on_si_ellipsometry_run.py) | Saves raw field data (not just R/T) for `postprocessing/jones_mueller_ellipsometry.py` to consume |

### Editing a script for your own structure

Each script has numbered `# EDIT (n):` comment blocks — material CSV paths
(or placeholder constants if the CSVs aren't found), layer thicknesses,
incident angle/azimuth/polarization (`s_amplitude`/`p_amplitude`, complex —
their ratio sets linear/circular/elliptical polarization), wavelength sweep,
and output path. No other part of the script should normally need touching.

### `trench/`, `via/` — not yet present

Land with Phase 3 (1D lamellar gratings) and Phase 4 (2D via/pillar arrays)
respectively — see [`phases.md`](../phases.md) and [`tasks.md`](../tasks.md).
Until then, any `Layer` given a `pattern` (instead of a plain `material`)
will raise `NotImplementedError` in `simulation.py`.

## Doubts already resolved for this folder (see `progress_log.md` 2026-07-19)

If you see dense, fast oscillation in an R/T-vs-wavelength plot for a
thick layer (e.g. the 12 um Si substrate in `sio2_on_si_thin_film.py`),
that's real Fabry-Perot interference (`Delta_lambda ~ lambda^2/(2 n t)`),
not a solver bug — increase the wavelength sample count if it looks
aliased/undersampled rather than assuming something's wrong with the
boundaries. RCWA has no mesh and no PML-equivalent to worry about; the
incidence/transmission media are exact semi-infinite half-spaces.
