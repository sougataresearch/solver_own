# `structures/trench/` — 1D-Periodic Lamellar Gratings (Trenches)

Every script here builds a 1D-periodic patterned layer (`Lattice1D` — periodic
along x, uniform/invariant along y) using `sougata_solver`'s own geometry
objects (`Slab`, or the `staircase.py`/`ocd.py` generators for tapered/OCD
cases), runs `Simulation.solve()`, and saves raw results. None of them plot —
that's [`postprocessing/`](../../postprocessing/)'s job
(`decisions.md` ADR-009/010).

## What each script is for

| Script | Structure it builds | What it does |
|---|---|---|
| [`trench_grating.py`](trench_grating.py) | Uniform (untapered) binary grating — Si ridge / air groove, fixed fill factor | The Phase 3 reference case: same geometry as `tests/oracles/rcwa_1d_gaylord.py`'s published Moharam/Gaylord (1995) benchmark, so its printed R/T can be compared directly against `tests/test_1d_grating.py`'s oracle-comparison assertions. **Start here** if you want to see the simplest possible 1D-periodic case working end to end. |
| [`trench_grating_ellipsometry_run.py`](trench_grating_ellipsometry_run.py) | Same geometry as `trench_grating.py` | Solves at one or more (wavelength, angle) points for both s- and p-polarization and saves the *raw per-diffraction-order* reflected field data (not just R/T) to CSV — feeds [`postprocessing/jones_mueller_per_order.py`](../../postprocessing/jones_mueller_per_order.py), which derives the per-order Jones matrix, Mueller matrix, and ellipsometric angles (Psi, Delta). Run this, then that. |
| [`tapered_trench.py`](tapered_trench.py) | Depth-tapered air-filled trench etched into a Si slab, staircase-discretized, plus a uniform residual Si layer beneath the taper | Geometry rebuilt from a real Lumerical FDTD reference file the project owner's senior provided (`decisions.md` ADR-036) — TCD/BCD/depth/period all traceable to that file's actual dialogs, not invented. Sweeps `num_slices` (1 through 64) at a fixed wavelength/angle so the staircase-discretization convergence trend is visible directly (~32 slices is where R stabilizes to ~1e-4). **Cross-validated against the project owner's own Lumerical RCWA run** of the same structure to ~0.5% R RMS / ~0.15% T RMS agreement (`postprocessing/overlay_tapered_trench_vs_lumerical.py`) — the most real-world-validated script in this folder. |
| [`trench_field_cross_section.py`](trench_field_cross_section.py) | Same untapered ridge/groove geometry as `trench_grating.py` | Reconstructs full `(Ex,Ey,Ez,Hx,Hy,Hz)` fields over an (x,z) cross-section through the grating and saves the raw grid to `.npz` — feeds [`postprocessing/plot_field_cross_section.py`](../../postprocessing/plot_field_cross_section.py), which plots `|E|^2` over that cross-section. Run this, then that. |

**Moved out**: the OCD sweep example (`trench_ocd_sweep.py`) now lives in a
separate sibling project, `../../../ocd_library/sweeps/trench_ocd_sweep.py`
— it's inverse-modeling/library-generation work, not RCWA solving, and the
project owner asked for it to live outside this solver's own repo (`decisions.md`
ADR has the full account). It still reuses `tapered_trench.build_geometry()`
directly (dynamically loaded across the two repos), so it can't silently
drift from this file's geometry.

## Suggested order to run these in

1. **`trench_grating.py`** — confirms the basic 1D-periodic machinery works and matches a published benchmark; no dependencies, safe starting point.
2. **`tapered_trench.py`** — the depth-tapered, real-device-derived case; sweeps `num_slices` itself, so just run it directly (the finest-slice-count spectrum is what gets saved to CSV).
3. **`trench_grating_ellipsometry_run.py`** → **`postprocessing/jones_mueller_per_order.py`** — if you want per-order Jones/Mueller/ellipsometric-angle output instead of plain R/T.
4. **`trench_field_cross_section.py`** → **`postprocessing/plot_field_cross_section.py`** — if you want to see the actual reconstructed field pattern inside/around the grating, not just R/T.
5. **`../../../ocd_library/sweeps/trench_ocd_sweep.py`** (separate sibling project) — if you want to see how sidewall angle affects R/T at several bottom-CD values around the real device.

Every run writes its own timestamped folder under
`outputs/YYYY_MM_DD/HH_MM_SS_<script-name>/` (gitignored by default —
see [`../../README.md`](../../README.md)'s "Output files" section), so
re-running any script with different `EDIT` values never overwrites a
previous run.

## Comparing against a real Lumerical RCWA/FDTD structure

If you're trying to replicate a real device from Lumerical (as `tapered_trench.py`
does): the axis convention is the thing most likely to trip you up.
`sougata_solver` **always** uses z as depth (`CONVENTIONS.md`) — Lumerical's
own convention varies per file (check the RCWA region's General tab for
`propagation axis`, or an FDTD region's boundary conditions/source injection
axis, rather than assuming). `decisions.md` ADR-035/036 record the full,
honest account of getting this translation right for `tapered_trench.py`,
including two wrong turns that were caught and corrected — worth reading
before trying to replicate another real structure the same way.
