"""Category 9 target 9.7 / Phase 7 example: real-space 2D (x, y) field map
through a circular Si pillar, at one chosen depth within the pillar layer.

Stack (top to bottom, matching pillar_array.py's geometry):
    air                       (incidence, semi-infinite)
    Si circular pillar / air  (finite thickness, 2D-periodic)
    air                       (exit, semi-infinite)

Per `decisions.md` ADR-009/010: this script (`structures/`) only solves the
physics and saves the raw reconstructed field grid -- it does not plot.
Run `postprocessing/plot_field_cross_section.py` afterward to view it.

Run with:  python structures/via/pillar_field_cross_section.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import modal_field_components, propagate_amplitudes, reconstruct_field_at_points, save_field_grid_npz
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.smatrix import SMatrixStack, interior_amplitudes

# ============================================================================
# EDIT: pillar geometry and materials (matches pillar_array.py)
# ============================================================================
PERIOD = 0.7e-6
PILLAR_RADIUS = 0.18e-6
THICKNESS = 0.46e-6
N_PILLAR = 3.48
N_BG = 1.0

# ============================================================================
# EDIT: incident light and order truncation
# ============================================================================
WAVELENGTH = 0.6e-6
INCIDENT_ANGLE_DEG = 0.0
NUM_ORDERS = 25
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: reconstruction grid and depth (fraction of the pillar's thickness)
# ============================================================================
N_XY = 80
DEPTH_FRACTION = 0.5  # 0.0 = top of pillar layer, 1.0 = bottom

OUTPUT_NPZ_PATH = "output_pillar_field_xy.npz"


def main():
    air = Material("air", N_BG**2)
    pillar = Material("pillar", N_PILLAR**2)
    pattern = Pattern(background=air, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=PILLAR_RADIUS, material=pillar)])
    lattice = Lattice((PERIOD, 0.0), (0.0, PERIOD))
    layer = Layer("pillar_layer", THICKNESS, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=NUM_ORDERS, incidence=air, transmission=air)

    excitation = PlaneWaveExcitation(
        WAVELENGTH, math.radians(INCIDENT_ANGLE_DEG), 0.0, s_amplitude=S_AMPLITUDE, p_amplitude=P_AMPLITUDE
    )
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    print(f"R={r:.4f}  T={t:.4f}  R+T={r + t:.4f}")

    omega = excitation.omega()
    thicknesses = [ly.thickness for ly in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)
    n2 = 2 * result.num_orders
    a_top, b_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)

    modes = result.all_modes[1]
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (result.g[:, 0] * lk[0, 0] + result.g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (result.g[:, 0] * lk[0, 1] + result.g[:, 1] * lk[1, 1])
    epsilon_hat = toeplitz_matrix(pattern, lattice, result.g, WAVELENGTH, inverse=False)
    einv = np.linalg.inv(epsilon_hat)

    depth = DEPTH_FRACTION * THICKNESS
    a_z, b_z = propagate_amplitudes(modes.q, depth, a_top, b_top)
    ex, ey, ez, hx, hy, hz = modal_field_components(omega, kx, ky, modes.q, modes.kp, modes.phi, einv, a_z, b_z)

    xs = (np.arange(N_XY) + 0.5) / N_XY * PERIOD
    ys = (np.arange(N_XY) + 0.5) / N_XY * PERIOD
    x_grid, y_grid = np.meshgrid(xs, ys, indexing="ij")
    ex_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, ex)
    ey_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, ey)
    ez_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, ez)
    hx_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, hx)
    hy_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, hy)
    hz_g = reconstruct_field_at_points(kx, ky, x_grid, y_grid, hz)

    output_dir = run_output_dir("pillar_field_cross_section")
    path = save_field_grid_npz(
        str(output_dir / OUTPUT_NPZ_PATH),
        x_grid,
        y_grid,
        depth,
        ex_g,
        ey_g,
        ez_g,
        hx_g,
        hy_g,
        hz_g,
        wavelength_m=WAVELENGTH,
        period_m=PERIOD,
        pillar_radius_m=PILLAR_RADIUS,
        thickness_m=THICKNESS,
        depth_m=depth,
    )
    write_run_metadata(
        output_dir,
        __file__,
        period_m=PERIOD,
        pillar_radius_m=PILLAR_RADIUS,
        thickness_m=THICKNESS,
        n_pillar=N_PILLAR,
        n_bg=N_BG,
        wavelength_m=WAVELENGTH,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        depth_fraction=DEPTH_FRACTION,
        grid=(N_XY, N_XY),
    )
    print(f"Saved field map to {path}")
    print(f"Run metadata: {output_dir / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
