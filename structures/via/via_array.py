"""Phase 4a example: circular air via array in a Si substrate on a square lattice.

Stack (top to bottom):
    air              (incidence, semi-infinite)
    Si substrate with air-filled circular via (finite thickness, 2D-periodic)
    Si              (exit, semi-infinite)

The geometry mirrors a through-silicon-via (TSV) scatterometry target:
Si bulk with air holes etched through it, as in the vendored
`EMTutorial/Scatterometry/ThroughSiliconVia` JCMsuite reference case.

Run with:  python structures/via/via_array.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

# ============================================================================
# EDIT: via geometry and materials
# ============================================================================
PERIOD = 0.7e-6          # lattice period (meters)
VIA_RADIUS = 0.18e-6     # via (hole) radius (meters)
THICKNESS = 0.46e-6      # via depth (meters)
N_SUBSTRATE = 3.48       # substrate index (e.g. Si)
N_VIA = 1.0              # via fill index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 7            # 2D Fourier order truncation parameter
S_AMPLITUDE = 1.0         # 1.0/0.0 = s-pol; 0.0/1.0 = p-pol
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.5e-6, 1.5e-6, 21)

OUTPUT_CSV = "output_via_RT.csv"


def build_geometry(period=None, via_radius=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission)."""
    period = period if period is not None else PERIOD
    via_radius = via_radius if via_radius is not None else VIA_RADIUS
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", N_VIA**2)
    substrate = Material("substrate", N_SUBSTRATE**2)

    # Via (air hole) centered in the unit cell at (period/2, period/2)
    pattern = Pattern(
        background=substrate,
        shapes=[Circle(center=(period / 2, period / 2), radius=via_radius, material=air)],
    )
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = [Layer("via_layer", thickness, pattern=pattern)]
    return layers, lattice, air, substrate


def main() -> None:
    layers, lattice, air, substrate = build_geometry()
    sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=substrate)

    reflectance = np.zeros(len(WAVELENGTHS))
    transmittance = np.zeros(len(WAVELENGTHS))

    print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'R+T':>8}")
    for i, wavelength in enumerate(WAVELENGTHS):
        excitation = PlaneWaveExcitation(
            wavelength=wavelength,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=0.0,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
        )
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        rt = reflectance[i] + transmittance[i]
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {rt:8.4f}")

    out = run_output_dir("via_array")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_nm,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        via_radius_m=VIA_RADIUS,
        thickness_m=THICKNESS,
        n_substrate=N_SUBSTRATE,
        n_via=N_VIA,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
    )
    print(f"\nSaved {len(WAVELENGTHS)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
