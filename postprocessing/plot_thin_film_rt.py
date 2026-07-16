"""Post-processing: plot R and T vs wavelength from a CSV a structures/
script already wrote (e.g. structures/thin_film/sio2_on_si_thin_film.py or
custom_multistack.py), and save that plot into the SAME run folder the CSV
came from -- no re-solving of the physics, just visualizing already-computed
raw data. Also supports overlaying a reference CSV (e.g. a Lumerical export)
for a direct visual + numeric cross-check.

Run the structures/ script first to produce the input CSV, then:

Run with:  python postprocessing/plot_thin_film_rt.py
"""

from pathlib import Path

import numpy as np

from sougata_solver.output_paths import find_latest_output

# ============================================================================
# EDIT (1): which run to plot.
#
# By default this looks up the most recent CSV with this filename across
# every outputs/YYYY-MM-DD/HH-MM-SS_<run_name>/ folder (see
# find_latest_output) -- so if you just ran sio2_on_si_thin_film.py, this
# plots that run with no editing needed.
#
# To plot a SPECIFIC earlier run instead of "whatever's latest", set
# EXPLICIT_INPUT_CSV_PATH to that run's exact CSV path (check its
# run_metadata.txt, written alongside it, to confirm it's the one you want).
# ============================================================================
INPUT_CSV_FILENAME = "output_RT.csv"  # matches sio2_on_si_thin_film.py / custom_multistack.py
EXPLICIT_INPUT_CSV_PATH = None  # e.g. r"...\outputs\2026-07-16\10-06-07_sio2_on_si_thin_film\output_RT.csv"

# ============================================================================
# EDIT (2): optional reference data overlay (e.g. exported from a Lumerical
# FDTD or `stackrt` TMM run) for a direct cross-check -- same three columns
# as our own CSV: header "wavelength_nm,R,T" (wavelength in nanometers).
# ============================================================================
REFERENCE_CSV_PATH = None  # e.g. r"C:\path\to\lumerical_export.csv"

# ============================================================================
# EDIT (3): plot display/output filename (saved into the SAME folder the
# input CSV came from, not a new outputs/ subfolder).
# ============================================================================
SHOW_PLOT = True
PLOT_FILENAME = "output_RT.png"


def main():
    input_path = Path(EXPLICIT_INPUT_CSV_PATH) if EXPLICIT_INPUT_CSV_PATH else find_latest_output(INPUT_CSV_FILENAME)
    print(f"Reading {input_path}")
    metadata_path = input_path.parent / "run_metadata.txt"
    if metadata_path.exists():
        print(f"Run metadata: {metadata_path}")
        print(metadata_path.read_text())

    data = np.genfromtxt(input_path, delimiter=",", names=True)
    wavelengths_nm, reflectance, transmittance = data["wavelength_nm"], data["R"], data["T"]

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(wavelengths_nm, reflectance, color="tab:blue", label="R (this solver)")
    ax.plot(wavelengths_nm, transmittance, color="tab:orange", label="T (this solver)")

    if REFERENCE_CSV_PATH:
        ref = np.genfromtxt(REFERENCE_CSV_PATH, delimiter=",", names=True)
        ax.plot(ref["wavelength_nm"], ref["R"], "--", color="tab:blue", label="R (reference)")
        ax.plot(ref["wavelength_nm"], ref["T"], "--", color="tab:orange", label="T (reference)")
        # Only a rigorous per-point check if the reference uses the same
        # wavelength grid as our own run; otherwise this is a rough guide.
        if len(ref["wavelength_nm"]) == len(wavelengths_nm) and np.allclose(
            ref["wavelength_nm"], wavelengths_nm, atol=1.0
        ):
            r_diff = np.max(np.abs(reflectance - ref["R"]))
            t_diff = np.max(np.abs(transmittance - ref["T"]))
            print(f"\nMax |R - R_reference| = {r_diff:.4e}, max |T - T_reference| = {t_diff:.4e}")
        else:
            print("\nReference wavelength grid differs from this run's -- overlay is a visual guide only.")

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance / Transmittance")
    ax.set_ylim(0, 1)
    ax.set_title("R, T vs wavelength")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    # Save into the SAME folder the input CSV came from -- this plot is a
    # derived view of that run's data, not a new run of its own.
    plot_path = input_path.parent / PLOT_FILENAME
    fig.savefig(plot_path, dpi=150)
    print(f"\nSaved plot to {plot_path}")
    if SHOW_PLOT:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
