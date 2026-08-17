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
from pathlib import Path

import numpy as np

from sougata_solver.geometry import Lattice1D
from sougata_solver.materials import Material
from sougata_solver.ocd import OCDTrapezoidParams, trapezoid_trench_layers
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import avoid_rayleigh_wood_anomalies, sweep_wavelength

# ============================================================================
# EDIT: trench geometry and materials (CD-first, per Category 11 target 11.1)
# ============================================================================
PERIOD = 2.032e-6        # 2.032 μm
TOP_CD = 1.383e-6        # 1.383 μm (TCD)
BOTTOM_CDS = [1.322e-6]  # 1.322 μm (BCD)
HEIGHT = 4.981e-6        # 4.981 μm depth
NUM_SLICES = 32          # Lock in converged slice count
# Load dispersive Silicon material from NK text file (matches Lumerical database)
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
ridge = Material.from_nk_file("Si", str(NK_DIR / "si_KLA.txt"), "nm")
air = Material("air", 1.0)


# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation, wavelength range
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 15  # per side; total Fourier orders = 2*NUM_ORDERS+1 for Lattice1D
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0
# 400-800nm, 1nm steps, exactly as requested -- any grid point landing
# exactly on a Rayleigh/Wood's-anomaly wavelength for this PERIOD/
# NUM_ORDERS/angle (troubleshooting.md's documented q==0 divide-by-zero,
# e.g. PERIOD=2.032um's exact 508nm order-4 collision) is nudged
# automatically; every other point, including the endpoints, is untouched.
WAVELENGTHS = avoid_rayleigh_wood_anomalies(
    np.linspace(0.4e-6, 0.8e-6, 401), period=PERIOD, num_orders=NUM_ORDERS, theta=math.radians(INCIDENT_ANGLE_DEG)
)

OUTPUT_CSV = "output_trench_ocd_sweep.csv"


def build_geometry(bottom_cd=None, period=None, top_cd=None, height=None, num_slices=None):
    """Returns (layers, lattice, incidence, transmission).

    `bottom_cd` defaults to `BOTTOM_CDS[0]` -- this file's own `main()`
    sweeps over every value in `BOTTOM_CDS` itself, calling this with each
    one explicitly.
    """
    bottom_cd = bottom_cd if bottom_cd is not None else BOTTOM_CDS[0]
    period = period if period is not None else PERIOD
    top_cd = top_cd if top_cd is not None else TOP_CD
    height = height if height is not None else HEIGHT
    num_slices = num_slices if num_slices is not None else NUM_SLICES

    lattice = Lattice1D(period)
    params = OCDTrapezoidParams(
        top_cd=top_cd, bottom_cd=bottom_cd, period=period, height=height,
        shape_material=ridge, background_material=air,
    )
    layers = trapezoid_trench_layers(params, num_slices=num_slices)
    return layers, lattice, air, air


def main() -> None:
    theta = math.radians(INCIDENT_ANGLE_DEG)

    all_wavelengths = []
    all_bottom_cd = []
    all_sidewall_angle = []
    all_reflectance = []
    all_transmittance = []

    for bottom_cd in BOTTOM_CDS:
        # sidewall_angle_deg is a cheap derived property (top_cd/bottom_cd/height
        # only, no physics computation) -- recomputed here for printing/metadata
        # since build_geometry() returns only (layers, lattice, incidence,
        # transmission), not the OCDTrapezoidParams object itself.
        params = OCDTrapezoidParams(top_cd=TOP_CD, bottom_cd=bottom_cd, period=PERIOD, height=HEIGHT, shape_material=ridge, background_material=air)
        layers, lattice, incidence, transmission = build_geometry(bottom_cd=bottom_cd)
        sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=incidence, transmission=transmission)

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
        ridge_material="Si (si_KLA.txt)",
        background_material="air",
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
