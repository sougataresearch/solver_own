"""Overlay `pillar_array.py`'s R/T output against a Lumerical RCWA
`grating_power` export (`OUTPUT_RCWA/Via/export_pillar_grating_power.lsf`).

Same known Lumerical-script quirk as
`overlay_tapered_trench_vs_lumerical.py`/`overlay_composite_grating_vs_lumerical.py`:
its `write()` appends rather than overwrites, so re-running the export
script without deleting the previous output leaves multiple header+data
blocks concatenated in one file -- `_load_lumerical_grating_power` below
finds every "lambda*" header line and uses only the last (most recent)
block.

Plots R and T on two panels, solver vs. Lumerical, and prints a max/RMS
difference -- Lumerical's wavelength grid is linearly interpolated onto
the solver's grid before differencing; this is a numerical comparison aid,
not a claim the two grids are identical.

Run with:  python postprocessing/overlay_pillar_vs_lumerical.py
"""

from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover
    raise ImportError("This script needs matplotlib: pip install -e '.[dev]'") from exc

# ============================================================================
# EDIT: the two files to overlay, and which Lumerical polarization column to
# use. "Rs"/"Ts" (s/TE) matches pillar_array.py's S_AMPLITUDE=1.0/
# P_AMPLITUDE=0.0 -- switch to "Rp"/"Tp" if the Lumerical excitation was
# actually p/TM.
# ============================================================================
SOLVER_CSV = (
    r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\outputs\2026_08_21\13_13_51_pillar_array\output_pillar_RT.csv"
  #  r"\16_07_30_pillar_array\output_pillar_RT.csv"
)
LUMERICAL_TXT = (
    r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\OUTPUT_RCWA\Via\pillar_lumerical_RCWA.txt"
  #  r"\pillar_lumerical_RCWA.txt"
)
LUMERICAL_R_FIELD = "Rs"
LUMERICAL_T_FIELD = "Ts"

# Saved into the same folder SOLVER_CSV came from -- this overlay is a
# derived view of that solver run, not a new run of its own.
PLOT_FILENAME = "pillar_vs_lumerical.png"


def _load_lumerical_grating_power(path: Path) -> np.ndarray:
    """Returns the LAST lambda*-headed block in `path` as a structured
    array, discarding any earlier block left over from write()'s
    append-not-overwrite behavior."""
    lines = path.read_text().strip().splitlines()
    header_indices = [i for i, line in enumerate(lines) if line.lower().startswith("lambda")]
    if not header_indices:
        raise ValueError(f"No header row starting with 'lambda' found in {path}")
    if len(header_indices) > 1:
        print(
            f"Note: {path.name} contains {len(header_indices)} concatenated export blocks "
            f"(write() appends rather than overwrites) -- using only the last one."
        )
    block = lines[header_indices[-1] :]
    return np.genfromtxt(block, delimiter=",", names=True)


def main():
    solver_path, lumerical_path = Path(SOLVER_CSV), Path(LUMERICAL_TXT)
    for p in (solver_path, lumerical_path):
        if not p.exists():
            raise FileNotFoundError(f"{p} does not exist")

    solver_data = np.genfromtxt(solver_path, delimiter=",", names=True)
    wl_solver, r_solver, t_solver = solver_data["wavelength_nm"], solver_data["R"], solver_data["T"]

    lumerical_data = _load_lumerical_grating_power(lumerical_path)
    wl_lumerical = lumerical_data["lambda_nm"]
    r_lumerical = lumerical_data[LUMERICAL_R_FIELD]
    t_lumerical = lumerical_data[LUMERICAL_T_FIELD]

    order = np.argsort(wl_lumerical)
    wl_lumerical, r_lumerical, t_lumerical = wl_lumerical[order], r_lumerical[order], t_lumerical[order]

    # Lumerical's grid isn't the same points as the solver's -- interpolate
    # onto the solver's grid, restricted to the overlapping wavelength range,
    # to get a numeric max/RMS difference rather than eyeballing the plot.
    in_range = (wl_solver >= wl_lumerical.min()) & (wl_solver <= wl_lumerical.max())
    r_lumerical_interp = np.interp(wl_solver[in_range], wl_lumerical, r_lumerical)
    t_lumerical_interp = np.interp(wl_solver[in_range], wl_lumerical, t_lumerical)
    r_diff = r_solver[in_range] - r_lumerical_interp
    t_diff = t_solver[in_range] - t_lumerical_interp

    print(f"Solver:    {len(wl_solver)} points, {wl_solver.min():.1f}-{wl_solver.max():.1f} nm")
    print(f"Lumerical: {len(wl_lumerical)} points, {wl_lumerical.min():.1f}-{wl_lumerical.max():.1f} nm "
          f"(field '{LUMERICAL_R_FIELD}'/'{LUMERICAL_T_FIELD}')")
    print(f"R: max|diff| = {np.max(np.abs(r_diff)):.4f}, RMS = {np.sqrt(np.mean(r_diff**2)):.4f}")
    print(f"T: max|diff| = {np.max(np.abs(t_diff)):.4f}, RMS = {np.sqrt(np.mean(t_diff**2)):.4f}")

    fig, (ax_r, ax_t) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    ax_r.plot(wl_solver, r_solver, label="sougata_solver", color="tab:blue")
    ax_r.plot(wl_lumerical, r_lumerical, "--", label="Lumerical RCWA", color="tab:orange")
    ax_r.set_ylabel("R")
    ax_r.legend()

    ax_t.plot(wl_solver, t_solver, label="sougata_solver", color="tab:blue")
    ax_t.plot(wl_lumerical, t_lumerical, "--", label="Lumerical RCWA", color="tab:orange")
    ax_t.set_xlabel("wavelength (nm)")
    ax_t.set_ylabel("T")
    ax_t.legend()

    fig.tight_layout()
    output_path = solver_path.parent / PLOT_FILENAME
    fig.savefig(output_path, dpi=150)
    print(f"\nSaved plot to {output_path}")


if __name__ == "__main__":
    main()
