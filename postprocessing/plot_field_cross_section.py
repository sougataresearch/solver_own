"""Post-processing: plot a reconstructed field cross-section/map from a
`.npz` file a `structures/` script already wrote (e.g.
`structures/trench/trench_field_cross_section.py` or
`structures/via/pillar_field_cross_section.py`), and save that plot into
the SAME run folder the `.npz` came from -- no re-solving of the physics,
just visualizing already-computed raw field data, per `decisions.md`
ADR-009/010 (same split `plot_thin_film_rt.py` already follows).

Handles both `.npz` layouts this project produces:
- trench cross-section (`x_m`, `z_m` axis keys) -- an (x, z) profile
  through the grating at fixed y=0.
- pillar/via field map (`x`, `y`, `z` keys, from `fields.save_field_grid_npz`)
  -- an (x, y) plane at one fixed depth.

Run the structures/ script first to produce the input `.npz`, then:

Run with:  python postprocessing/plot_field_cross_section.py
"""

from pathlib import Path

import numpy as np

# ============================================================================
# EDIT (1): which run to plot -- the exact .npz path a structures/ field
# script printed after it finished (also written into that same folder's
# run_metadata.txt).
# ============================================================================
INPUT_NPZ_PATH = r"C:\path\to\outputs\YYYY_MM_DD\HH_MM_SS_trench_field_cross_section\output_trench_field_xz.npz"

# ============================================================================
# EDIT (2): which field component to plot (intensity |Ex|^2+|Ey|^2+|Ez|^2,
# or one specific complex component's real/imaginary/magnitude).
# ============================================================================
COMPONENT = "intensity"  # "intensity", "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"

# ============================================================================
# EDIT (3): plot display/output filename (saved into the SAME folder the
# input .npz came from, not a new outputs/ subfolder).
# ============================================================================
SHOW_PLOT = True
PLOT_FILENAME = "output_field_cross_section.png"


def main():
    input_path = Path(INPUT_NPZ_PATH)
    if not input_path.exists():
        raise FileNotFoundError(
            f"{input_path} does not exist -- set INPUT_NPZ_PATH at the top of this "
            "script to the exact .npz path a structures/ field script printed after it ran."
        )
    print(f"Reading {input_path}")
    metadata_path = input_path.parent / "run_metadata.txt"
    if metadata_path.exists():
        print(f"Run metadata: {metadata_path}")
        print(metadata_path.read_text())

    data = np.load(input_path)
    is_cross_section = "x_m" in data.files  # trench_field_cross_section.py's layout
    if is_cross_section:
        horizontal, vertical = data["x_m"] * 1e9, data["z_m"] * 1e9
        horizontal_label, vertical_label = "x (nm)", "z (nm, depth into layer)"
    else:  # fields.save_field_grid_npz's layout (pillar_field_cross_section.py)
        horizontal, vertical = data["x"][:, 0] * 1e9, data["y"][0, :] * 1e9
        horizontal_label, vertical_label = "x (nm)", "y (nm)"

    if COMPONENT == "intensity":
        values = np.abs(data["Ex"]) ** 2 + np.abs(data["Ey"]) ** 2 + np.abs(data["Ez"]) ** 2
        title = "|E|^2"
    else:
        values = np.abs(data[COMPONENT])
        title = f"|{COMPONENT}|"

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    if is_cross_section:
        mesh = ax.pcolormesh(horizontal, vertical, values, shading="auto", cmap="inferno")
        ax.invert_yaxis()  # z increases downward (depth into the stack)
    else:
        mesh = ax.pcolormesh(horizontal, vertical, values.T, shading="auto", cmap="inferno")
        ax.set_aspect("equal")
    fig.colorbar(mesh, ax=ax, label=title)
    ax.set_xlabel(horizontal_label)
    ax.set_ylabel(vertical_label)
    ax.set_title(f"{title} field {'cross-section' if is_cross_section else 'map'}")
    fig.tight_layout()

    plot_path = input_path.parent / PLOT_FILENAME
    fig.savefig(plot_path, dpi=150)
    print(f"\nSaved plot to {plot_path}")
    if SHOW_PLOT:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
