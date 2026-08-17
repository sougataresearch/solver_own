"""Category 9 / Phase 7 example: real-space field cross-section (x vs z)
through a 1D lamellar grating (trench).

Stack (top to bottom, matching trench_grating.py's geometry):
    air                          (incidence, semi-infinite)
    Si ridge/air groove grating  (finite thickness, 1D-periodic along x)
    air                          (exit, semi-infinite)

Per `decisions.md` ADR-009/010: this script (`structures/`) only solves the
physics and saves the raw reconstructed field grid -- it does not plot.
Run `postprocessing/plot_field_cross_section.py` afterward to view it.

Run with:  python structures/trench/trench_field_cross_section.py
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import modal_field_components, propagate_amplitudes, reconstruct_field_at_points
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.smatrix import SMatrixStack, interior_amplitudes

# ============================================================================
# EDIT: grating geometry and materials (matches trench_grating.py)
# ============================================================================
PERIOD = 0.7e-6
FILL_FACTOR = 0.3
THICKNESS = 0.46e-6
N_RIDGE = 3.48
N_GROOVE = 1.0

# ============================================================================
# EDIT: incident light and order truncation
# ============================================================================
WAVELENGTH = 0.6e-6
INCIDENT_ANGLE_DEG = 10.0
NUM_ORD = 9
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: reconstruction grid (x spans two periods to show repetition; z spans
# the grating layer only)
# ============================================================================
N_X = 200
N_Z = 80
X_SPAN_PERIODS = 2.0

OUTPUT_NPZ_PATH = "output_trench_field_xz.npz"


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
    pattern = layers[0].pattern
    num_orders = 2 * NUM_ORD + 1
    sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=transmission)

    excitation = PlaneWaveExcitation(
        WAVELENGTH, math.radians(INCIDENT_ANGLE_DEG), 0.0, s_amplitude=S_AMPLITUDE, p_amplitude=P_AMPLITUDE
    )
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    print(f"R={r:.4f}  T={t:.4f}  R+T={r + t:.4f}")

    omega = excitation.omega()
    thicknesses = [ly.thickness for ly in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)
    n2 = 2 * num_orders
    a_top, b_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)

    modes = result.all_modes[1]
    lk = lattice.reciprocal_vectors()
    kx0, _ky0 = excitation.k_parallel(1.0)
    kx = kx0 + 2 * np.pi * result.g[:, 0] * lk[0, 0]
    ky = np.zeros_like(kx)
    epsilon_inv_hat = toeplitz_matrix(pattern, lattice, result.g, WAVELENGTH, inverse=True)

    x = np.linspace(-0.5 * X_SPAN_PERIODS * PERIOD, 0.5 * X_SPAN_PERIODS * PERIOD, N_X)
    z = np.linspace(0.0, THICKNESS, N_Z)
    y = np.zeros_like(x)

    field_shape = (N_Z, N_X)
    ex_xz = np.zeros(field_shape, dtype=complex)
    ey_xz = np.zeros(field_shape, dtype=complex)
    ez_xz = np.zeros(field_shape, dtype=complex)
    hx_xz = np.zeros(field_shape, dtype=complex)
    hy_xz = np.zeros(field_shape, dtype=complex)
    hz_xz = np.zeros(field_shape, dtype=complex)
    for iz, depth in enumerate(z):
        a_z, b_z = propagate_amplitudes(modes.q, float(depth), a_top, b_top)
        ex, ey, ez, hx, hy, hz = modal_field_components(
            omega, kx, ky, modes.q, modes.kp, modes.phi, epsilon_inv_hat, a_z, b_z
        )
        ex_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, ex)
        ey_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, ey)
        ez_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, ez)
        hx_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, hx)
        hy_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, hy)
        hz_xz[iz, :] = reconstruct_field_at_points(kx, ky, x, y, hz)

    # A cross-section (x vs z, one fixed y=0 line at every depth) doesn't
    # match `fields.save_field_grid_npz`'s single-horizontal-plane
    # signature (that function is used as-is by
    # `structures/via/pillar_field_cross_section.py` below, target 9.7's
    # actual XY-grid deliverable) -- saved directly here instead, with
    # axis names that say what they actually are (x, z), not borrowed
    # (x, y) names that would misdescribe this cross-section's second axis.
    output_dir = run_output_dir("trench_field_cross_section")
    path = output_dir / OUTPUT_NPZ_PATH
    np.savez(
        path,
        x_m=x,
        z_m=z,
        Ex=ex_xz,
        Ey=ey_xz,
        Ez=ez_xz,
        Hx=hx_xz,
        Hy=hy_xz,
        Hz=hz_xz,
        wavelength_m=WAVELENGTH,
        period_m=PERIOD,
        thickness_m=THICKNESS,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
    )
    write_run_metadata(
        output_dir,
        __file__,
        period_m=PERIOD,
        fill_factor=FILL_FACTOR,
        thickness_m=THICKNESS,
        n_ridge=N_RIDGE,
        n_groove=N_GROOVE,
        wavelength_m=WAVELENGTH,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=num_orders,
        grid=(N_X, N_Z),
    )
    print(f"Saved field cross-section to {path}")
    print(f"Run metadata: {output_dir / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
