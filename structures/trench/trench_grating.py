"""1D-periodic lamellar grating (line/space trench), the Phase 3 capability.

Stack (top to bottom):
    air              (incidence, semi-infinite)
    Si ridge/air groove grating (finite thickness, 1D-periodic along x)
    air              (exit, semi-infinite)

Same binary-grating geometry (period, fill factor, ridge index, thickness)
used as the Phase 3 system-test benchmark cross-checked against
`tests/oracles/rcwa_1d_gaylord.py` (hand-transcribed from the vendored
`Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples/1D_Grating_Gaylord_{TE,TM}.py`,
citing Moharam, Grann, Pommet & Gaylord 1995), so this script's printed
numbers can be compared directly against `tests/test_1d_grating.py`'s
oracle-comparison assertions.

Run with:  python structures/trench/trench_grating.py
"""

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

# ============================================================================
# EDIT: grating geometry and materials
# ============================================================================
PERIOD = 0.7e-6          # grating period (meters)
FILL_FACTOR = 0.3        # fraction of the period occupied by the ridge
THICKNESS = 0.46e-6      # grating (groove) depth (meters)
N_RIDGE = 3.48           # ridge index (e.g. Si)
N_GROOVE = 1.0           # groove index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORD = 15              # orders per side; total Fourier orders = 2*NUM_ORD+1
S_AMPLITUDE = 1.0         # 1.0/0.0 = pure TE (s); 0.0/1.0 = pure TM (p)
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.5e-6, 1.5e-6, 101)

OUTPUT_CSV_PATH = "output_trench_RT.csv"


def build_geometry(period=None, fill_factor=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission)."""
    period = period if period is not None else PERIOD
    fill_factor = fill_factor if fill_factor is not None else FILL_FACTOR
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", 1.0)
    ridge = Material("ridge", N_RIDGE**2)

    pattern = Pattern(background=air)
    pattern.add(Slab(center_x=-period * (1 - fill_factor) / 2, halfwidth=0.5 * fill_factor * period, material=ridge))

    lattice = Lattice1D(period)
    layers = [Layer("grating", thickness, pattern=pattern)]
    return layers, lattice, air, air


def main():
    layers, lattice, air, transmission = build_geometry()
    num_orders = 2 * NUM_ORD + 1
    sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=transmission)

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
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {reflectance[i] + transmittance[i]:8.4f}")

    if OUTPUT_CSV_PATH:
        output_dir = run_output_dir("trench_grating")
        write_run_metadata(
            output_dir,
            __file__,
            period_m=PERIOD,
            fill_factor=FILL_FACTOR,
            thickness_m=THICKNESS,
            n_ridge=N_RIDGE,
            n_groove=N_GROOVE,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            num_orders=num_orders,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(output_path, table, delimiter=",", header="wavelength_nm,R,T", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")
        print(f"Run metadata: {output_dir / 'run_metadata.txt'}")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
