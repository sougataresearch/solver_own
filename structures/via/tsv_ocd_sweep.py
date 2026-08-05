"""Category 11 target 11.5: reproducible TSV (through-silicon-via)/OCD
example sweep.

Stack (top to bottom):
    air                             (incidence, semi-infinite)
    Si via, radius tapering          (num_slices staircase layers,
      top_cd/2 -> bottom_cd/2         2D-periodic, see staircase.py)
    air                             (exit, semi-infinite)

Sweeps a *sidewall angle* by varying `ocd.OCDTrapezoidParams.bottom_cd`
at fixed `top_cd`/`height` -- the same CD-first parametrization target
11.1 defines, reused here for a via (a circular cross-section, not
`ocd.trapezoid_trench_layers`'s trench-specific `Slab`) via
`staircase.staircase_circle_layers` directly (radius = CD/2). Prints each
sweep point's *derived* `sidewall_angle_deg` alongside R/T, and saves it
into `run_metadata.txt` for every point, satisfying Category 11's exit
criterion ("parameter changes are traceable in metadata").

Run with:  python structures/via/tsv_ocd_sweep.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.materials import Material
from sougata_solver.ocd import OCDTrapezoidParams
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_circle_layers

# ============================================================================
# EDIT: TSV geometry and materials (CD-first, per Category 11 target 11.1)
# ============================================================================
TOP_CD = 0.48e-6      # via diameter at the incidence-side face (meters)
BOTTOM_CDS = [0.48e-6, 0.40e-6, 0.32e-6, 0.24e-6, 0.16e-6]  # sweep -> increasing sidewall angle
SPACING = 0.22e-6     # edge-to-edge gap to the neighboring via at the top face (meters)
PERIOD = TOP_CD + SPACING
HEIGHT = 0.46e-6      # via depth (meters)
N_VIA = 3.48          # via index (e.g. Si)
N_BG = 1.0            # background index (air)
NUM_SLICES = 16        # staircase resolution (fixed here; see tapered_via.py for a slice-count study)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
WAVELENGTH = 1.0e-6
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 5
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

OUTPUT_CSV = "output_tsv_ocd_sweep.csv"


def main() -> None:
    air = Material("air", N_BG**2)
    via = Material("via", N_VIA**2)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    excitation = PlaneWaveExcitation(
        wavelength=WAVELENGTH,
        theta=math.radians(INCIDENT_ANGLE_DEG),
        phi=0.0,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
    )

    param_sets = [
        OCDTrapezoidParams(top_cd=TOP_CD, bottom_cd=bottom_cd, period=PERIOD, height=HEIGHT, shape_material=via, background_material=air)
        for bottom_cd in BOTTOM_CDS
    ]

    reflectance = np.zeros(len(param_sets))
    transmittance = np.zeros(len(param_sets))
    sidewall_angles = np.zeros(len(param_sets))

    print(f"{'bottom_cd (nm)':>15}  {'sidewall (deg)':>15}  {'R':>10}  {'T':>10}  {'R+T':>10}")
    for i, params in enumerate(param_sets):
        layers = staircase_circle_layers(
            center=(PERIOD / 2, PERIOD / 2),
            top_radius=params.top_cd / 2.0,
            bottom_radius=params.bottom_cd / 2.0,
            thickness=params.height,
            num_slices=NUM_SLICES,
            shape_material=params.shape_material,
            background_material=params.background_material,
        )
        sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=air)
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        sidewall_angles[i] = params.sidewall_angle_deg
        rt = reflectance[i] + transmittance[i]
        print(f"{params.bottom_cd * 1e9:15.1f}  {sidewall_angles[i]:15.3f}  {reflectance[i]:10.6f}  {transmittance[i]:10.6f}  {rt:10.6f}")

    out = run_output_dir("tsv_ocd_sweep")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([BOTTOM_CDS, sidewall_angles, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="bottom_cd_m,sidewall_angle_deg,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        top_cd_m=TOP_CD,
        bottom_cds_m=BOTTOM_CDS,
        spacing_m=SPACING,
        height_m=HEIGHT,
        n_via=N_VIA,
        n_bg=N_BG,
        num_slices=NUM_SLICES,
        wavelength_m=WAVELENGTH,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        sidewall_angles_deg=sidewall_angles.tolist(),
    )
    print(f"\nSaved {len(param_sets)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
