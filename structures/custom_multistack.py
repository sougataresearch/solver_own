"""Template: an arbitrary-length multilayer thin-film stack, any materials.

Copy this file for any new multistack (thin film / DBR / anti-reflection
coating / etc.) instead of editing sio2_on_si_thin_film.py in place.
Everything you're likely to change is in the numbered EDIT blocks below.

LIMITATION (see phases.md / troubleshooting.md): only *uniform* (unpatterned)
layers work today -- each layer has a thickness (z-direction) only, and is
treated as extending infinitely in x/y. Trench/via/pillar patterning (real
x/y dimensions: line width, radius, pitch) is not implemented yet.

Run with:  python structures/custom_multistack.py
"""

import math
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.geometry import Lattice
from pyrcwa.layer import Layer
from pyrcwa.materials import Material
from pyrcwa.simulation import Simulation

NK_DIR = Path(__file__).resolve().parent.parent.parent / "NK_FILE"


def _parse_refractiveindex_csv(csv_path: str):
    """refractiveindex.info-style CSV: an 'n' block, optionally an 'k' block."""
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
    """Dispersive material from a refractiveindex.info-style CSV (n, optional k)."""
    wl_n, n_vals, wl_k, k_vals = _parse_refractiveindex_csv(csv_path)
    scale = {"um": 1e-6, "nm": 1e-9, "m": 1.0}[wavelength_unit]
    n_interp = interp1d(wl_n * scale, n_vals, bounds_error=False, fill_value="extrapolate")
    k_interp = interp1d(wl_k * scale, k_vals, bounds_error=False, fill_value="extrapolate")
    return Material.from_nk(name, lambda wl: float(n_interp(wl)), lambda wl: float(k_interp(wl)))


# ============================================================================
# EDIT (1): define every material you need, one of three ways
# ============================================================================
air = Material("air", 1.0)                              # constant, lossless
tio2 = Material("TiO2", 2.4**2)                          # constant n, lossless (n=2.4)
sio2 = Material.from_nk("SiO2", n=1.46, k=0.0)           # constant n and k
si = material_from_csv("Si", str(NK_DIR / "si_nk.csv"))  # dispersive, from CSV

# ============================================================================
# EDIT (2): the stack itself -- as many Layer(name, thickness, material=...)
# entries as you want, top to bottom. incidence/transmission below are the
# semi-infinite half-spaces above/below this list (not part of it).
# ============================================================================
INCIDENCE_MATERIAL = air     # what light travels through before hitting the stack
TRANSMISSION_MATERIAL = si   # semi-infinite substrate below the stack
# ^ if you want the substrate to have a *finite* thickness instead of being
#   semi-infinite, add it as a normal Layer(...) in the list below and set
#   TRANSMISSION_MATERIAL back to whatever is truly underneath it (e.g. air).

layers = [
    Layer("TiO2", 60e-9, material=tio2),
    Layer("SiO2", 90e-9, material=sio2),
    Layer("TiO2", 60e-9, material=tio2),
    Layer("SiO2", 90e-9, material=sio2),
]

# ============================================================================
# EDIT (3): incident light -- angle (degrees), polarization
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
AZIMUTHAL_ANGLE_DEG = 0.0
S_AMPLITUDE = 1.0   # s_amplitude/p_amplitude set polarization state:
P_AMPLITUDE = 0.0   #   (1,0)=s-pol, (0,1)=p-pol, (1,1)=45deg, (1,1j)=circular

# ============================================================================
# EDIT (4): wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.4e-6, 0.8e-6, 401)

# ============================================================================
# EDIT (5): where to save results (set to None to skip saving)
# ============================================================================
OUTPUT_CSV_PATH = "output_multistack_RT.csv"


def main():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused -- only matters for patterned layers
    sim = Simulation(
        lattice, layers, num_orders=1,
        incidence=INCIDENCE_MATERIAL, transmission=TRANSMISSION_MATERIAL,
    )

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
        absorptance = 1 - reflectance[i] - transmittance[i]
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {absorptance:8.4f}")

    if OUTPUT_CSV_PATH:
        absorptance = 1.0 - reflectance - transmittance
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance, absorptance])
        np.savetxt(OUTPUT_CSV_PATH, table, delimiter=",", header="wavelength_nm,R,T,A", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {OUTPUT_CSV_PATH}")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
