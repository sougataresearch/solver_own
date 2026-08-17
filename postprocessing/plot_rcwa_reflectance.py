"""Plot a single RCWA_module export (e.g. 'lambda(m)' + 'Y' columns,
comma-delimited) vs wavelength, y-axis fixed 0 to 1 -- no solver data, no
overlay (see RCWA_plot_norm.py for overlaying it against this solver's own
output_R.csv).

Run with:  python postprocessing/plot_rcwa_reflectance.py
"""

from pathlib import Path

import numpy as np

# ============================================================================
# EDIT (1): the RCWA_module export to plot.
# ============================================================================
RCWA_TXT_FILE = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\OUTPUT_RCWA\Thin_Film\17_08_26\Multi\0_degree.txt"

# ============================================================================
# EDIT (2): plot display/output filename (saved into sougata_solver/PLOT).
# ============================================================================
SHOW_PLOT = True
PLOT_DIR = Path(__file__).resolve().parents[1] / "PLOT"
PLOT_FILENAME = "rcwa_reflectance.png"


def _find_field(names: tuple[str, ...], prefix: str) -> str | None:
    return next((n for n in names if n.lower().startswith(prefix.lower())), None)


def main():
    input_path = Path(RCWA_TXT_FILE)
    if not input_path.exists():
        raise FileNotFoundError(
            f"{input_path} does not exist -- set RCWA_TXT_FILE at the top of this script."
        )
    print(f"Reading {input_path}")

    delimiter = "\t" if input_path.read_text().split("\n", 1)[0].count("\t") else ","
    data = np.genfromtxt(input_path, delimiter=delimiter, names=True)
    names = data.dtype.names
    # "lambda" covers RCWA_module's "lambda(m)" header (parsed as "lambdam").
    wl_field = _find_field(names, "wavelength") or _find_field(names, "lambda")
    # "y" matches RCWA_module's reflectance column, named "Y".
    r_field = _find_field(names, "r") or _find_field(names, "y")
    if wl_field is None or r_field is None:
        raise ValueError(f"Could not find wavelength/R columns in {input_path} -- got fields {names}")
    wavelengths_nm, reflectance = data[wl_field], data[r_field]

    if wl_field.startswith("lambda") and np.max(wavelengths_nm) < 1e-3:
        # RCWA_module reports lambda in meters; convert to nm for this plot's axis.
        wavelengths_nm = wavelengths_nm * 1e9

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(wavelengths_nm, reflectance, color="magenta", linewidth=1.5, label="RCWA")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_ylim(0, 1)
    ax.set_xlim(wavelengths_nm.min(), wavelengths_nm.max())
    ax.set_title("RCWA reflectance vs wavelength")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_path = PLOT_DIR / PLOT_FILENAME
    fig.savefig(plot_path, dpi=150)
    print(f"Saved plot to {plot_path}")
    if SHOW_PLOT:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
