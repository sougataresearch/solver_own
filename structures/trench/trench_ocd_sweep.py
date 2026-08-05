"""Category 11 target 11.6: reproducible trench/OCD example sweep.

Stack (top to bottom):
    air                              (incidence, semi-infinite)
    Si trapezoidal trench             (num_slices staircase layers,
      top_cd -> bottom_cd              1D-periodic, see ocd.py/staircase.py)
    air                              (exit, semi-infinite)

Sweeps a wavelength range (`sweep.sweep_wavelength`, Category 8) for each
of several trapezoid sidewall angles, built via
`ocd.OCDTrapezoidParams`/`ocd.trapezoid_trench_layers` (target 11.1/11.2) --
demonstrating the OCD parameter object driving both a via (11.5) and a
trench (here) geometry. Every sweep point's `sidewall_angle_deg` is saved
into `run_metadata.txt`, per Category 11's "parameter changes are
traceable in metadata" exit criterion.

Run with:  python structures/trench/trench_ocd_sweep.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.geometry import Lattice1D
from sougata_solver.materials import Material
from sougata_solver.ocd import OCDTrapezoidParams, trapezoid_trench_layers
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import sweep_wavelength

# ============================================================================
# EDIT: trench geometry and materials (CD-first, per Category 11 target 11.1)
# ============================================================================
PERIOD = 0.7e-6
TOP_CD = 0.3e-6
BOTTOM_CDS = [0.3e-6, 0.2e-6, 0.1e-6]  # sweep -> increasing sidewall angle
HEIGHT = 0.46e-6
N_RIDGE = 3.48
N_BG = 1.0
NUM_SLICES = 16

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation, wavelength range
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 15  # per side; total Fourier orders = 2*NUM_ORDERS+1 for Lattice1D
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0
WAVELENGTHS = np.linspace(0.8e-6, 1.2e-6, 21)

OUTPUT_CSV = "output_trench_ocd_sweep.csv"


def main() -> None:
    air = Material("air", N_BG**2)
    ridge = Material("ridge", N_RIDGE**2)
    lattice = Lattice1D(PERIOD)
    theta = math.radians(INCIDENT_ANGLE_DEG)

    all_wavelengths = []
    all_bottom_cd = []
    all_sidewall_angle = []
    all_reflectance = []
    all_transmittance = []

    for bottom_cd in BOTTOM_CDS:
        params = OCDTrapezoidParams(top_cd=TOP_CD, bottom_cd=bottom_cd, period=PERIOD, height=HEIGHT, shape_material=ridge, background_material=air)
        layers = trapezoid_trench_layers(params, num_slices=NUM_SLICES)
        sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=air)

        sweep = sweep_wavelength(sim, WAVELENGTHS, theta=theta, phi=0.0, s_amplitude=S_AMPLITUDE, p_amplitude=P_AMPLITUDE)
        r, t = sweep.reflectance(), sweep.transmittance()

        print(f"\nbottom_cd={bottom_cd * 1e9:.1f} nm, sidewall_angle={params.sidewall_angle_deg:.3f} deg")
        print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'R+T':>8}")
        for wavelength, r_i, t_i in zip(WAVELENGTHS, r, t):
            print(f"{wavelength * 1e9:16.1f}  {r_i:8.4f}  {t_i:8.4f}  {r_i + t_i:8.4f}")

        all_wavelengths.extend(WAVELENGTHS)
        all_bottom_cd.extend([bottom_cd] * len(WAVELENGTHS))
        all_sidewall_angle.extend([params.sidewall_angle_deg] * len(WAVELENGTHS))
        all_reflectance.extend(r)
        all_transmittance.extend(t)

    out = run_output_dir("trench_ocd_sweep")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([all_wavelengths, all_bottom_cd, all_sidewall_angle, all_reflectance, all_transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_m,bottom_cd_m,sidewall_angle_deg,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        top_cd_m=TOP_CD,
        bottom_cds_m=BOTTOM_CDS,
        height_m=HEIGHT,
        n_ridge=N_RIDGE,
        n_bg=N_BG,
        num_slices=NUM_SLICES,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
    )
    print(f"\nSaved {len(all_wavelengths)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
