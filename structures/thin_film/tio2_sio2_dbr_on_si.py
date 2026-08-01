"""TiO2/SiO2 DBR stack on a finite-thickness Si substrate (air below).

Copied from custom_multistack.py -- see that file's docstring for the
general template. Everything you're likely to change is in the numbered
EDIT blocks below.

Run with:  python structures/thin_film/tio2_sio2_dbr_on_si.py
"""

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
TIO2_DEVORE_O_YML = (
    NK_DIR / "refractiveindex.info-database" / "database" / "data" / "main"
    / "TiO2" / "nk" / "Devore-o.yml"
)


# ============================================================================
# EDIT (1): define every material you need -- all dispersive, from file
# ============================================================================
air = Material("air", 1.0)                                                # constant, lossless
# Devore (1951) rutile ordinary-ray Sellmeier fit, valid 0.43-1.53 um,
# k=0 (no absorption data in this source) -- see Material.from_refractiveindex_formula_file.
tio2 = Material.from_refractiveindex_formula_file("TiO2", str(TIO2_DEVORE_O_YML))
sio2 = Material.from_nk_file("SiO2", str(NK_DIR / "SiO2_nk.csv"))         # dispersive, from file
si = Material.from_nk_file("Si", str(NK_DIR / "si_nk.csv"))               # dispersive, from file

# ============================================================================
# EDIT (2): the stack itself -- as many Layer(name, thickness, material=...)
# entries as you want, top to bottom. incidence/transmission below are the
# semi-infinite half-spaces above/below this list (not part of it).
# ============================================================================
INCIDENCE_MATERIAL = air   # what light travels through before hitting the stack
TRANSMISSION_MATERIAL = air  # semi-infinite exit medium below the finite Si substrate

SI_SUBSTRATE_THICKNESS = 2e-3  # <-- EDIT: set this (meters) before running, e.g. 500e-6

layers = [
    Layer("TiO2", 50e-9, material=tio2),
    Layer("SiO2", 50e-9, material=sio2),
    Layer("TiO2", 50e-9, material=tio2),
    Layer("SiO2", 50e-9, material=sio2),
    Layer("Si", SI_SUBSTRATE_THICKNESS, material=si),  # finite substrate, not semi-infinite
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
# EDIT (5): where to save results (set OUTPUT_CSV_PATH to None to skip
# saving)
# ============================================================================
RUN_NAME = "tio2_sio2_dbr_on_si"
OUTPUT_CSV_PATH = "output_dbr_RT.csv"


def main():
    if SI_SUBSTRATE_THICKNESS is None:
        raise ValueError("Set SI_SUBSTRATE_THICKNESS (meters) before running -- see EDIT (2).")

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
        output_dir = run_output_dir(RUN_NAME)
        write_run_metadata(
            output_dir,
            __file__,
            layers=[(layer.name, layer.thickness) for layer in layers],
            incidence_material=INCIDENCE_MATERIAL.name,
            transmission_material=TRANSMISSION_MATERIAL.name,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        absorptance = 1.0 - reflectance - transmittance
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance, absorptance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(output_path, table, delimiter=",", header="wavelength_nm,R,T,A", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
