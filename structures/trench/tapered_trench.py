"""Phase 5 example: linearly-tapered Si trench (1D grating), staircase-discretized.

Stack (top to bottom):
    air                                (incidence, semi-infinite)
    Si ridge, halfwidth tapering        (num_slices staircase layers,
      top_halfwidth -> bottom_halfwidth  1D-periodic, see staircase.py)
    air                                (exit, semi-infinite)

Sweeps `num_slices` at a fixed wavelength/angle and prints R/T, mirroring
`tests/test_staircase.py`'s `test_tapered_trench_converges_with_increasing_num_slices`.
Per ADR-010, this script only prints/saves raw numbers -- no plotting here.

Run with:  python structures/trench/tapered_trench.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_slab_layers

# ============================================================================
# EDIT: tapered trench geometry and materials (FDTD-style tcd/bcd/depth/
# spacing naming -- see conversation citing phases.md Phase 5/staircase.py;
# tcd/bcd here are ridge *widths*, not the trench (groove) opening)
# ============================================================================
TCD = 0.35e-6                # top critical dimension: ridge width at incidence-side face (meters)
BCD = 0.105e-6               # bottom critical dimension: ridge width at transmission-side face (meters)
SPACING = 0.35e-6            # groove width at the incidence-side face (meters)
PERIOD = TCD + SPACING        # grating period (meters)
THICKNESS = 0.46e-6            # trench depth (meters)
N_RIDGE = 3.48                  # ridge index (e.g. Si)
N_GROOVE = 1.0                  # groove index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
WAVELENGTH = 1.0e-6
INCIDENT_ANGLE_DEG = 0.0
NUM_ORD = 11                # orders per side; total Fourier orders = 2*NUM_ORD+1
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: slice-count sweep
# ============================================================================
SLICE_COUNTS = [1, 2, 4, 8, 16, 32, 64]

OUTPUT_CSV = "output_tapered_trench_convergence.csv"


def main() -> None:
    air = Material("air", N_GROOVE**2)
    ridge = Material("ridge", N_RIDGE**2)
    lattice = Lattice1D(PERIOD)
    num_orders = 2 * NUM_ORD + 1
    excitation = PlaneWaveExcitation(
        wavelength=WAVELENGTH,
        theta=math.radians(INCIDENT_ANGLE_DEG),
        phi=0.0,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
    )

    reflectance = np.zeros(len(SLICE_COUNTS))
    transmittance = np.zeros(len(SLICE_COUNTS))

    print(f"{'num_slices':>10}  {'R':>10}  {'T':>10}  {'R+T':>10}")
    for i, num_slices in enumerate(SLICE_COUNTS):
        layers = staircase_slab_layers(
            center_x=0.0,
            top_halfwidth=0.5 * TCD,
            bottom_halfwidth=0.5 * BCD,
            thickness=THICKNESS,
            num_slices=num_slices,
            shape_material=ridge,
            background_material=air,
        )
        sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=air)
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        rt = reflectance[i] + transmittance[i]
        print(f"{num_slices:10d}  {reflectance[i]:10.6f}  {transmittance[i]:10.6f}  {rt:10.6f}")

    out = run_output_dir("tapered_trench")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([SLICE_COUNTS, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="num_slices,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        tcd_m=TCD,
        bcd_m=BCD,
        spacing_m=SPACING,
        thickness_m=THICKNESS,
        n_ridge=N_RIDGE,
        n_groove=N_GROOVE,
        wavelength_m=WAVELENGTH,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=num_orders,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        slice_counts=SLICE_COUNTS,
    )
    print(f"\nSaved {len(SLICE_COUNTS)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
