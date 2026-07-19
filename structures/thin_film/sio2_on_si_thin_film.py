"""SiO2 (50 nm) on Si (12 um), both semi-infinite air above and below.

Stack (top to bottom):
    air   (incidence, semi-infinite)
    SiO2  (50 nm, finite layer)
    Si    (12 um, finite layer -- thick enough to show real interference
           fringes in R/T vs wavelength, not treated as infinite)
    air   (exit, semi-infinite)

Run with:  python structures/thin_film/sio2_on_si_thin_film.py
"""

import math
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation


def _parse_refractiveindex_csv(csv_path: str):
    """Parse the refractiveindex.info CSV export format: an 'n' block
    (header 'wl,n') and, optionally, a separate 'k' block (header 'wl,k')
    after a blank line -- the two blocks can even be on different
    wavelength grids. If there's no 'k' block, k is treated as zero
    (lossless) at every wavelength.
    """
    with open(csv_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        try:
            float(line.split(",")[0])
        except ValueError:
            if current:
                blocks.append(current)
            current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"No numeric data found in {csv_path!r}")

    def parse_block(block):
        arr = np.array([[float(x) for x in line.split(",")] for line in block])
        return arr[:, 0], arr[:, 1]

    wl_n, n_vals = parse_block(blocks[0])
    if len(blocks) > 1:
        wl_k, k_vals = parse_block(blocks[1])
    else:
        wl_k, k_vals = wl_n, np.zeros_like(wl_n)
    return wl_n, n_vals, wl_k, k_vals


def material_from_csv(name: str, csv_path: str, wavelength_unit: str = "um") -> Material:
    """refractiveindex.info-style CSV: an 'n' block and an optional 'k' block."""
    wl_n, n_vals, wl_k, k_vals = _parse_refractiveindex_csv(csv_path)
    scale = {"um": 1e-6, "nm": 1e-9, "m": 1.0}[wavelength_unit]

    n_interp = interp1d(wl_n * scale, n_vals, bounds_error=False, fill_value="extrapolate")
    k_interp = interp1d(wl_k * scale, k_vals, bounds_error=False, fill_value="extrapolate")
    return Material.from_nk(
        name,
        lambda wl: float(n_interp(wl)),
        lambda wl: float(k_interp(wl)),
    )


# ============================================================================
# EDIT (1): point these at your actual Si and SiO2 n,k CSV files
# ============================================================================
# Resolved relative to this script's own location (Solver_own/NK_FILE), not
# a hardcoded drive letter -- keeps working after moving/copying the whole
# Solver_own folder to another device where it might not be on G:\.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
SI_CSV_PATH = str(NK_DIR / "Si_nk.csv")       # <-- put your Si data file path here
SIO2_CSV_PATH = str(NK_DIR / "SiO2_nk.csv")   # <-- put your SiO2 data file path here
CSV_WAVELENGTH_UNIT = "um"                  # "um", "nm", or "m" -- match your file

# ============================================================================
# EDIT (2): layer thicknesses (meters)
# ============================================================================
SIO2_THICKNESS = 50e-9   # 50 nm
SI_THICKNESS = 2e-3     # 2 mm

# ============================================================================
# EDIT (3): incident light -- angle (degrees), polarization
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0     # angle from surface normal
AZIMUTHAL_ANGLE_DEG = 0.0    # usually 0 unless the sample is in-plane anisotropic

# s_amplitude / p_amplitude are COMPLEX -- their relative magnitude and phase
# set the polarization state:
S_AMPLITUDE = 1.0            # linear s-pol (TE)
P_AMPLITUDE = 0.0
# S_AMPLITUDE, P_AMPLITUDE = 1.0, 1.0                  # linear, 45 deg
# S_AMPLITUDE, P_AMPLITUDE = 1.0, 1j                   # circular (90 deg phase)
# S_AMPLITUDE, P_AMPLITUDE = 1.0, 0.5j                 # elliptical
# for unpolarized light, run twice (S_AMPLITUDE=1,P_AMPLITUDE=0 and vice
# versa) and average the resulting R and T.

# ============================================================================
# EDIT (4): wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.4e-6, 0.8e-6, 401)  # 400-800 nm, 401 points

# ============================================================================
# EDIT (5): where to save results (set to None to skip saving)
#
# Plotting is NOT done here -- this script only builds the structure, runs
# the solver, and saves raw R/T data. Run postprocessing/plot_thin_film_rt.py
# afterward to plot it (it finds this run's CSV automatically and saves the
# plot into this same output folder).
# ============================================================================
OUTPUT_CSV_PATH = "output_RT.csv"  # filename only; saved under outputs/YYYY_MM_DD/


def main():
    try:
        si = material_from_csv("Si", SI_CSV_PATH, CSV_WAVELENGTH_UNIT)
        sio2 = material_from_csv("SiO2", SIO2_CSV_PATH, CSV_WAVELENGTH_UNIT)
    except OSError:
        print(f"Could not find {SI_CSV_PATH!r} / {SIO2_CSV_PATH!r} -- using placeholder constants instead.\n")
        si = Material("Si", (3.9 + 0.02j) ** 2)      # rough placeholder, NOT real Si data
        sio2 = Material("SiO2", 1.46**2)              # rough placeholder, NOT real SiO2 data

    air = Material("air", 1.0)

    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused for uniform (unpatterned) layers
    layers = [
        Layer("SiO2", SIO2_THICKNESS, material=sio2),
        Layer("Si", SI_THICKNESS, material=si),
    ]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=air)

    # Collect results into plain arrays as we go -- this is the pattern to
    # reuse whenever you need the data beyond just printing it (plotting,
    # fitting, saving, comparing against a measured spectrum, etc.).
    reflectance = np.zeros(len(WAVELENGTHS))
    transmittance = np.zeros(len(WAVELENGTHS))

    print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'A':>8}")
    for i, wavelength in enumerate(WAVELENGTHS):
        excitation = PlaneWaveExcitation(
            wavelength=wavelength,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=math.radians(AZIMUTHAL_ANGLE_DEG),
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
        )
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {1 - reflectance[i] - transmittance[i]:8.4f}")

    if OUTPUT_CSV_PATH:
        output_dir = run_output_dir("sio2_on_si_thin_film")
        write_run_metadata(
            output_dir,
            __file__,
            si_csv_path=SI_CSV_PATH,
            sio2_csv_path=SIO2_CSV_PATH,
            sio2_thickness_m=SIO2_THICKNESS,
            si_thickness_m=SI_THICKNESS,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        absorptance = 1.0 - reflectance - transmittance
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance, absorptance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(
            output_path,
            table,
            delimiter=",",
            header="wavelength_nm,R,T,A",
            comments="",
        )
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")
        print(f"Run metadata: {output_dir / 'run_metadata.txt'}")
        print("To plot this run: python postprocessing/plot_thin_film_rt.py")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
