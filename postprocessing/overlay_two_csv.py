"""Overlay any two R-vs-wavelength CSV/TXT files on one plot.

Generic sibling of `KLA_plot_norm.py` for the common case where *both*
files are this solver's own `output_R.csv` runs (e.g. comparing two
polarization states) rather than one KLA export + one solver run --
reuses `KLA_plot_norm._load_wavelength_r` (same auto-detected-format
loader) so the two scripts never diverge on parsing.

Run with:  python postprocessing/overlay_two_csv.py
"""

from pathlib import Path

import numpy as np

from sougata_solver.postprocessing.RCWA_plot_norm import _load_wavelength_r

# ============================================================================
# EDIT: the two files to overlay. Set these to two DIFFERENT files --
# pointing both at the same path just draws one curve twice.
# ============================================================================
PATH_1 = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\outputs\2026_08_12\15_27_13_sio2_sio_ni_sio2_on_semi_infinite_si\output_multistack_RT.csv"
LABEL_1 = "linear_15deg"

PATH_2 = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\outputs\2026_08_12\15_29_00_sio2_sio_ni_sio2_on_semi_infinite_si\output_multistack_RT.csv"
LABEL_2 = "linear_30deg"

# Saved into the SAME folder PATH_2 came from -- this overlay is a derived
# view of that run, not a new run of its own.
PLOT_FILENAME = "reflectance_overlay.png"


def main():
    path_1, path_2 = Path(PATH_1), Path(PATH_2)
    if path_1 == path_2:
        raise ValueError(f"PATH_1 and PATH_2 are the same file ({path_1}) -- point them at two different runs.")
    for p in (path_1, path_2):
        if not p.exists():
            raise FileNotFoundError(f"{p} does not exist")

    print(f"Plotting {LABEL_1}: {path_1}")
    wavelength_1, reflectance_1 = _load_wavelength_r(path_1)

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.plot(wavelength_1, reflectance_1, color="magenta", linewidth=1.5, label=LABEL_1)

    print(f"Overlaying {LABEL_2}: {path_2}")
    wavelength_2, reflectance_2 = _load_wavelength_r(path_2)
    plt.plot(wavelength_2, reflectance_2, "--", color="tab:blue", linewidth=1.5, label=LABEL_2)

    if len(wavelength_1) == len(wavelength_2) and np.allclose(wavelength_1, wavelength_2, atol=1.0):
        max_diff = np.max(np.abs(reflectance_1 - reflectance_2))
        mean_diff = np.mean(np.abs(reflectance_1 - reflectance_2))
        print(f"Max |R_{LABEL_1} - R_{LABEL_2}| = {max_diff:.4e}, mean = {mean_diff:.4e}")
    else:
        print("Wavelength grids differ -- overlay is a visual guide only, not a per-point comparison.")

    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Reflectance")
    plt.ylim(0, 1)
    plt.xlim(min(wavelength_1.min(), wavelength_2.min()), max(wavelength_1.max(), wavelength_2.max()))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plot_path = path_2.parent / PLOT_FILENAME
    plt.savefig(plot_path, dpi=300)
    print(f"Saved plot to {plot_path}")
    plt.show()


if __name__ == "__main__":
    main()
