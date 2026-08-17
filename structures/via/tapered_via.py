"""Phase 5 example: linearly-tapered Si via, staircase-discretized.

Stack (top to bottom):
    air                          (incidence, semi-infinite)
    Si via, radius tapering       (num_slices staircase layers,
      top_radius -> bottom_radius  2D-periodic, see staircase.py)
    air                          (exit, semi-infinite)

Sweeps `num_slices` at a fixed wavelength/angle and prints R/T so the
convergence-vs-N trend (`staircase.py`'s docstring: this phase's
correctness evidence, no external oracle exists for tapered geometry) is
visible directly, mirroring `tests/test_staircase.py`'s
`test_tapered_via_converges_with_increasing_num_slices`. Per ADR-010, this
script only prints/saves raw numbers -- no plotting here.

Run with:  python structures/via/tapered_via.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_circle_layers

# ============================================================================
# EDIT: tapered via geometry and materials (FDTD-style tcd/bcd/depth/spacing
# naming -- see conversation citing phases.md Phase 5/staircase.py; tcd/bcd
# here are via *diameters*, top and bottom face)
# ============================================================================
TCD = 0.48e-6                 # top critical dimension: via diameter at incidence-side face (meters)
BCD = 0.20e-6                 # bottom critical dimension: via diameter at transmission-side face (meters)
SPACING = 0.22e-6             # edge-to-edge gap to the neighboring via at the top face (meters)
PERIOD = TCD + SPACING         # lattice period (meters)
THICKNESS = 0.46e-6             # via depth (meters)
N_VIA = 3.48                     # via index (e.g. Si)
N_BG = 1.0                        # background index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
WAVELENGTH = 1.0e-6
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 5
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: slice-count sweep
# ============================================================================
SLICE_COUNTS = [1, 2, 4, 8, 16, 32]

OUTPUT_CSV = "output_tapered_via_convergence.csv"


def build_geometry(num_slices=None, period=None, tcd=None, bcd=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission).

    `num_slices` defaults to `SLICE_COUNTS[-1]` (the finest/most
    representative geometry) -- this file's own `main()` sweeps over every
    value in `SLICE_COUNTS` itself, calling this with each one explicitly.
    """
    num_slices = num_slices if num_slices is not None else SLICE_COUNTS[-1]
    period = period if period is not None else PERIOD
    tcd = tcd if tcd is not None else TCD
    bcd = bcd if bcd is not None else BCD
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", N_BG**2)
    via = Material("via", N_VIA**2)
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = staircase_circle_layers(
        center=(period / 2, period / 2),
        top_radius=0.5 * tcd,
        bottom_radius=0.5 * bcd,
        thickness=thickness,
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

    reflectance = np.zeros(len(SLICE_COUNTS))
    transmittance = np.zeros(len(SLICE_COUNTS))

    print(f"{'num_slices':>10}  {'R':>10}  {'T':>10}  {'R+T':>10}")
    for i, num_slices in enumerate(SLICE_COUNTS):
        layers, lattice, air, transmission = build_geometry(num_slices=num_slices)
        sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=transmission)
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        rt = reflectance[i] + transmittance[i]
        print(f"{num_slices:10d}  {reflectance[i]:10.6f}  {transmittance[i]:10.6f}  {rt:10.6f}")

    out = run_output_dir("tapered_via")
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
        n_via=N_VIA,
        n_bg=N_BG,
        wavelength_m=WAVELENGTH,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        slice_counts=SLICE_COUNTS,
    )
    print(f"\nSaved {len(SLICE_COUNTS)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
