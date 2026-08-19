# `postprocessing/` — Run These Second

Every script here takes raw output a [`structures/`](../structures/) script
(or an external tool export) already wrote to disk, and derives a plot or a
comparison from it — **none of them re-solve the physics** (`decisions.md`
ADR-009/010's `structures/`-builds-and-solves vs. `postprocessing/`-derives-
and-shows split). Run the matching `structures/` script first to produce the
input file, then the `postprocessing/` script to visualize or compare it.

| Script | Input | Purpose |
|---|---|---|
| [`plot_thin_film_rt.py`](plot_thin_film_rt.py) | one `output_*_RT.csv`/`output_R.csv` from a `structures/thin_film/*.py` run | Plots R vs. wavelength, saved into the same run folder the CSV came from. No reference overlay. |
| [`RCWA_plot_norm.py`](RCWA_plot_norm.py) | two R-vs-wavelength files, mixing this solver's own CSV output with a KLA reflectance-calculator export (auto-detects which is which per file) | Overlays a solver run against a KLA export on one plot. |
| [`overlay_two_csv.py`](overlay_two_csv.py) | two of this solver's own `output_R.csv` runs | Generic sibling of `RCWA_plot_norm.py` for comparing two solver runs directly (e.g. two polarization states), reusing the same auto-detected-format loader so the two scripts never diverge on parsing. |
| [`plot_rcwa_reflectance.py`](plot_rcwa_reflectance.py) | a single `RCWA_module` export (`lambda(m)`/`Y` columns) | Plots that export alone vs. wavelength, y-axis fixed 0–1 — no solver data, no overlay. |
| [`overlay_composite_grating_vs_lumerical.py`](overlay_composite_grating_vs_lumerical.py) | `structures/thin_film/multistack_composite_grating.py`'s CSV + a Lumerical RCWA `grating_power` export (`export_grating_power.lsf`) | Overlays solver vs. Lumerical R/T on two panels, prints max/RMS difference (interpolated onto the solver's wavelength grid). The comparison behind `decisions.md` ADR-034. |
| [`overlay_tapered_trench_vs_lumerical.py`](overlay_tapered_trench_vs_lumerical.py) | `structures/trench/tapered_trench.py`'s CSV + a Lumerical RCWA `grating_power` export (`export_trench_grating_power.lsf`, project root) | Same overlay/max-RMS-diff approach as above, for the depth-tapered trench structure. The comparison behind `decisions.md` ADR-036 (~0.5% R RMS, ~0.15% T RMS agreement once materials and harmonic-order truncation were matched between the two tools). |
| [`plot_field_cross_section.py`](plot_field_cross_section.py) | a `.npz` field grid from `structures/trench/trench_field_cross_section.py` or `structures/via/pillar_field_cross_section.py` | Plots `\|E\|^2` over the reconstructed cross-section/map, handling both the trench (x,z) and pillar/via (x,y) `.npz` layouts. |
| [`jones_mueller_ellipsometry.py`](jones_mueller_ellipsometry.py) | raw reflected-field CSV from `structures/thin_film/sio2_on_si_ellipsometry_run.py` | Assembles the Jones reflection matrix, Mueller matrix, and ellipsometric angles (Psi, Delta) — zeroth-order/uniform-stack case. |
| [`jones_mueller_per_order.py`](jones_mueller_per_order.py) | raw per-order reflected-field CSV from `structures/trench/trench_grating_ellipsometry_run.py` | Same Jones/Mueller/Psi/Delta derivation as above, but per diffraction order (patterned-layer case) — reuses the same `sougata_solver.polarimetry` functions so the two scripts can't silently drift apart on convention. |

## The two Lumerical-comparison ".lsf" export scripts

`overlay_composite_grating_vs_lumerical.py` and
`overlay_tapered_trench_vs_lumerical.py` each expect a plain
`lambda_nm,Rs,Ts,Rp,Tp` text file, produced by running the matching `.lsf`
script (Lumerical Script Language) inside Lumerical itself against a solved
RCWA simulation:

- `export_grating_power.lsf` (referenced in `overlay_composite_grating_vs_lumerical.py`'s
  docstring) — for the composite-grating structure.
- [`export_trench_grating_power.lsf`](../export_trench_grating_power.lsf) (project root) —
  for the tapered-trench structure.

Both extract Lumerical's `grating_power` result (`Rs_power`/`Ts_power`/
`Rp_power`/`Tp_power`), sum over all diffraction orders, and write the total
R/T per wavelength — the same summed quantity this solver's own
`reflectance()`/`transmittance()` report, so the two are directly
comparable. Known gotcha: Lumerical's `write()` **appends**, not overwrites
— delete any existing output file before re-running an export, or the
corresponding overlay script's loader will (correctly) use only the most
recent header+data block and print a note that it did so.

## Adding a new postprocessing script

Follow `decisions.md` ADR-009/010: never call `Simulation.solve()` or import
`sougata_solver.simulation` from a `postprocessing/` script — load
already-written raw data (CSV/`.npz`/`.txt`) and derive/plot from that only.
Save output plots into the same run folder the input data came from (see
`src/sougata_solver/output_paths.py`), not a new location.
