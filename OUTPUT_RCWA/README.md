# `OUTPUT_RCWA/` — Lumerical RCWA reference data and scripts

External-tool (Lumerical RCWA) exports and the `.lsf` scripts that
produced them, kept alongside the `sougata_solver` structure each one
validates. Every `.lsf` file's own header comment states which
`structures/*.py` script it corresponds to and which `decisions.md` ADR
the comparison is recorded in -- this file is just the index.

| Folder | `.lsf` script | What it computes | Corresponds to |
|---|---|---|---|
| `Thin_Film/Complex_Multi/` | `multistack_composite_grating.lsf` | Exports `grating_power` (`Rs_power`/`Ts_power`/`Rp_power`/`Tp_power`), summed over all diffraction orders, to `composite_grating_lumerical_RCWA.txt` | `structures/thin_film/multistack_composite_grating.py` — `decisions.md` ADR-034 |
| `Thin_Film/Multi_layer/Incident45/` | `plot_R_linear_15_30deg.lsf` | Linear-polarization reflectance at 15/30 deg from the specular (n=0,m=0) order only | `structures/thin_film/custom_multistack.py` (`linear_15deg`/`linear_30deg` states) — `decisions.md` ADR-033 |
| `Thin_Film/Multi_layer/Incident45/` | `plot_R_linear_and_pure_SP_15_30deg.lsf` | Same, plus the pure-S/pure-P curves for a visual sanity check | `structures/thin_film/custom_multistack.py` — `decisions.md` ADR-033 |

**Two conventions worth knowing before writing a new one here:**
- `grating_power`'s `Rs_power`/`Rp_power`/etc. attributes are
  **per-diffraction-order** (indexed by wavelength x n-order x m-order),
  not pre-summed totals -- see ADR-034. For a plain uniform thin film
  (no lateral pattern), every non-specular order is exactly zero, so
  indexing straight to the specular order (as the two `Incident45`
  scripts do) is equivalent to summing. For a genuinely patterned
  structure, sum over all `n`/`m` first (`pinch()` then `sum(sum(...,3),2)`,
  as `multistack_composite_grating.lsf` does) -- indexing a single order
  would silently drop real diffracted power.
- `write()` appends rather than overwrites -- delete any existing output
  file before re-running an export script, or the loader on the
  `sougata_solver` side will see multiple concatenated header+data blocks
  (handled defensively in `postprocessing/overlay_composite_grating_vs_lumerical.py`,
  but cleaner to avoid).
