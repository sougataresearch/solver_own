"""Post-processing: plot R vs wavelength from a single CSV a structures/
script already wrote (e.g. structures/thin_film/sio2_on_si_thin_film.py or
custom_multistack.py), and save that plot into the SAME run folder the CSV
came from -- no re-solving of the physics, just visualizing already-computed
raw data, and no reference overlay (see KLA_plot_norm.py for that).

Run the structures/ script first to produce the input CSV, then:

Run with:  python postprocessing/plot_thin_film_rt.py
"""

from pathlib import Path

import numpy as np

# ============================================================================
# EDIT (1): which run to plot -- the exact CSV path a structures/ script
# printed after it finished (also written into that same folder's
# run_metadata.txt, if you need to find it again later). Copy the whole
# path, including the .csv extension.
# ============================================================================
INPUT_CSV_PATH = r"C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver\outputs\2026_08_18\13_53_18_multistack_composite_grating\output_multistack_composite_grating_RT.csv"

# ============================================================================
# EDIT (2): plot display/output filename (saved into the SAME folder the
# input CSV came from, not a new outputs/ subfolder).
# ============================================================================
SHOW_PLOT = True
PLOT_FILENAME = "output_R.png"


def main():
    input_path = Path(INPUT_CSV_PATH)
    if not input_path.exists():
        raise FileNotFoundError(
            f"{input_path} does not exist -- set INPUT_CSV_PATH at the top of this "
            "script to the exact CSV path a structures/ script printed after it ran."
        )
    print(f"Reading {input_path}")
    metadata_path = input_path.parent / "run_metadata.txt"
    if metadata_path.exists():
        print(f"Run metadata: {metadata_path}")
        print(metadata_path.read_text())

    data = np.genfromtxt(input_path, delimiter=",", names=True)
    wavelengths_nm, reflectance = data["wavelength_nm"], data["R"]

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(wavelengths_nm, reflectance, color="tab:blue", label="R (this solver)")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_ylim(0, 1)
    ax.set_title("R vs wavelength")
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
