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

    air = Material("air", N_BG**2)
    via = Material("via", N_VIA**2)
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = staircase_circle_layers(
        center=(period / 2, period / 2),
        top_radius=top_cd / 2.0,
        bottom_radius=bottom_cd / 2.0,
        thickness=height,
        num_slices=num_slices,
        shape_material=via,
        background_material=air,
    )
    return layers, lattice, air, air


def main() -> None:
    excitation = PlaneWaveExcitation(
        wavelength=WAVELENGTH,
        theta=math.radians(INCIDENT_ANGLE_DEG),
        phi=0.0,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
    )

    reflectance = np.zeros(len(BOTTOM_CDS))
    transmittance = np.zeros(len(BOTTOM_CDS))
    sidewall_angles = np.zeros(len(BOTTOM_CDS))

    print(f"{'bottom_cd (nm)':>15}  {'sidewall (deg)':>15}  {'R':>10}  {'T':>10}  {'R+T':>10}")
    for i, bottom_cd in enumerate(BOTTOM_CDS):
        # sidewall_angle_deg is a cheap derived property (top_cd/bottom_cd/
        # height only, no physics computation) -- recomputed here for
        # printing/metadata since build_geometry() returns only
        # (layers, lattice, incidence, transmission), not the
        # OCDTrapezoidParams object itself.
        params = OCDTrapezoidParams(
            top_cd=TOP_CD, bottom_cd=bottom_cd, period=PERIOD, height=HEIGHT,
            shape_material=Material("via", N_VIA**2), background_material=Material("air", N_BG**2),
        )
        layers, lattice, incidence, transmission = build_geometry(bottom_cd=bottom_cd)
        sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=incidence, transmission=transmission)
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
    print(f"\nSaved {len(BOTTOM_CDS)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
