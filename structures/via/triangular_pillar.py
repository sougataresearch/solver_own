"""Category 4 targets 4.4/4.5 example: triangular Si pillar array on a
square lattice (air background) -- the end-to-end RCWA example required by
Category 4's exit criteria ("each new shape has geometry-only tests and one
end-to-end RCWA example").

Stack (top to bottom):
    air                          (incidence, semi-infinite)
    Si triangular pillar / air   (finite thickness, 2D-periodic)
    air                          (exit, semi-infinite)

Run with:  python structures/via/triangular_pillar.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice, Pattern, Polygon
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

# ============================================================================
# EDIT: pillar geometry and materials
# ============================================================================
PERIOD = 0.7e-6              # lattice period (meters)
TRIANGLE_HEIGHT = 0.32e-6    # apex-to-base height (meters)
TRIANGLE_BASE = 0.36e-6      # base width (meters)
THICKNESS = 0.46e-6          # pillar height (meters)
N_PILLAR = 3.48               # pillar refractive index (e.g. Si)
N_BG = 1.0                    # background index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 7
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.5e-6, 1.5e-6, 21)

OUTPUT_CSV = "output_triangular_pillar_RT.csv"


def _equilateral_triangle_vertices(base: float, height: float) -> tuple[tuple[float, float], ...]:
    """CCW, centered on the centroid so `Polygon.center` is the geometric center."""
    y_top = 2.0 / 3.0 * height
    y_bottom = -1.0 / 3.0 * height
    return ((0.0, y_top), (-base / 2, y_bottom), (base / 2, y_bottom))


def build_geometry(period=None, triangle_base=None, triangle_height=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission)."""
    period = period if period is not None else PERIOD
    triangle_base = triangle_base if triangle_base is not None else TRIANGLE_BASE
    triangle_height = triangle_height if triangle_height is not None else TRIANGLE_HEIGHT
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", N_BG**2)
    pillar = Material("pillar", N_PILLAR**2)

    pattern = Pattern(
        background=air,
        shapes=[
            Polygon(
                center=(period / 2, period / 2),
                vertices=_equilateral_triangle_vertices(triangle_base, triangle_height),
                material=pillar,
            )
        ],
    )
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = [Layer("pillar_layer", thickness, pattern=pattern)]
    return layers, lattice, air, air


def main() -> None:
    layers, lattice, air, transmission = build_geometry()
    sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=transmission)

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

    out = run_output_dir("triangular_pillar")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_nm,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        triangle_height_m=TRIANGLE_HEIGHT,
        triangle_base_m=TRIANGLE_BASE,
        thickness_m=THICKNESS,
        n_pillar=N_PILLAR,
        n_bg=N_BG,
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
