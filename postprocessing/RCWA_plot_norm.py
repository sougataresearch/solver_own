"""Plot one R-vs-wavelength CSV/TXT, then overlay a second one on top of it.

Works with any mix of:
  - this solver's own output_R.csv ("wavelength_nm,R" header, comma-delimited), or
  - a KLA reflectance-calculator export ('"Wavelength (nm)"' + "Reflectance"
    columns, tab-delimited) -- format is auto-detected per file, so PATH_1/
    PATH_2 don't need to be told which kind they are.

Run with:  python postprocessing/KLA_plot_norm.py
"""

from pathlib import Path

import numpy as np

# ============================================================================
# EDIT: the two files to overlay -- KLA_TXT_FILE is drawn first (solid
# magenta, matching this script's original single-curve style),
# OUTPUT_R_CSV is overlaid on top of it (dashed). Set these to two
# DIFFERENT files -- pointing both at the same path just draws one curve
# twice.
# ============================================================================
KLA_TXT_FILE = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\OUTPUT_RCWA\Thin_Film\17_08_26\Multi\0_degree.txt"
LABEL_1 = "RCWA"

OUTPUT_R_CSV = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\outputs\2026_08_18\13_53_18_multistack_composite_grating\output_multistack_composite_grating_RT.csv"
LABEL_2 = "solver"

# Saved into the SAME folder OUTPUT_R_CSV came from -- this overlay is a
# derived view of that solver run, not a new run of its own.
PLOT_FILENAME = "reflectance_overlay.png"


def _find_field(names: tuple[str, ...], prefix: str) -> str | None:
    return next((n for n in names if n.lower().startswith(prefix.lower())), None)


def _load_wavelength_r(path: Path) -> tuple[np.ndarray, np.ndarray]:
    delimiter = "\t" if path.read_text().split("\n", 1)[0].count("\t") else ","
    data = np.genfromtxt(path, delimiter=delimiter, names=True)
    names = data.dtype.names
    # "wavelength"/"lambda" cover our own output and KLA exports; "lambda" also
    # covers RCWA_module exports whose header is "lambda(m)" -> parsed as "lambdam".
    wl_field = _find_field(names, "wavelength") or _find_field(names, "lambda")
    # "r"/"y" matches our own "R", KLA's "Reflectance", and RCWA_module's "Y".
    r_field = _find_field(names, "r") or _find_field(names, "y")
    if wl_field is None or r_field is None:
        raise ValueError(f"Could not find wavelength/R columns in {path} -- got fields {names}")

    wavelength = data[wl_field]
    if wl_field.startswith("lambda") and np.max(wavelength) < 1e-3:
        # RCWA_module reports lambda in meters; convert to nm to match this
        # solver's convention so both curves share the same x-axis scale.
        wavelength = wavelength * 1e9
    return wavelength, data[r_field]


def main():
    path_1, path_2 = Path(KLA_TXT_FILE), Path(OUTPUT_R_CSV)
    if path_1 == path_2:
        raise ValueError(
            f"KLA_TXT_FILE and OUTPUT_R_CSV are the same file ({path_1}) -- point them at two different runs."
        )
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
        print(f"Max |R_{LABEL_1} - R_{LABEL_2}| = {max_diff:.4e}")
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
