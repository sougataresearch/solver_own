"""Category 1 target 1.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): patterned
(2D-periodic) layers containing diagonal or in-plane-coupled anisotropic
materials.

Tiers enforced here, per `rules.md` Testing Requirements:
- regression: an isotropic pattern reduces to
  `solve_layer_eigenmodes_patterned`'s (Phase 4a) result.
- regression: a spatially-uniform (single-shape-equals-background) pattern
  reduces to `solve_layer_eigenmodes_uniform_inplane`'s (target 1.4) result.
- physical invariant: energy conservation for a genuinely patterned,
  Hermitian (lossless) anisotropic case.
- guard: longitudinal coupling anywhere in the pattern still raises
  `NotImplementedError` naming target 1.5.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.eigenmodes import (
    solve_layer_eigenmodes_patterned,
    solve_layer_eigenmodes_patterned_inplane,
    solve_layer_eigenmodes_uniform_inplane,
)
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix, toeplitz_matrix_component
from sougata_solver.geometry import Circle, Lattice, Pattern, Rectangle
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 0.7
WAVELENGTH = 1.0


def _pillar_lattice_and_g(num_orders=5):
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, num_orders, "circular")
    return lattice, g


def _kx_ky(lattice, g, omega, theta_deg=0.0):
    lk = lattice.reciprocal_vectors()
    kx0 = omega * np.sin(np.radians(theta_deg))
    kx = kx0 + 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    return kx, ky


# ---------------------------------------------------------------------------
# Regression: isotropic pattern reduces to solve_layer_eigenmodes_patterned
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 20.0])
def test_isotropic_pattern_reduces_to_phase4a_solver(theta_deg):
    lattice, g = _pillar_lattice_and_g()
    omega = 2 * np.pi / WAVELENGTH
    kx, ky = _kx_ky(lattice, g, omega, theta_deg)

    air = Material("air", 1.0)
    si = Material("si", 3.48**2)
    pattern = Pattern(background=air, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=si)])

    epsilon_hat = toeplitz_matrix(pattern, lattice, g, WAVELENGTH, inverse=False)
    modes_isotropic_path = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)

    exx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 0)
    exy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 1)
    eyx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 0)
    eyy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 1)
    ezz = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 2, 2)
    modes_tensor_path = solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, exx, exy, eyx, eyy, ezz)

    # Off-diagonal Toeplitz matrices should be exactly zero for an isotropic pattern.
    assert np.allclose(exy, 0.0)
    assert np.allclose(eyx, 0.0)
    assert np.allclose(exx, ezz)
    assert np.allclose(eyy, ezz)

    q_sq_iso = np.sort_complex(modes_isotropic_path.q**2)
    q_sq_tensor = np.sort_complex(modes_tensor_path.q**2)
    assert q_sq_tensor == pytest.approx(q_sq_iso, abs=1e-9)


# ---------------------------------------------------------------------------
# Regression: spatially-uniform "pattern" reduces to target 1.4's uniform solver
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 15.0])
def test_uniform_pattern_reduces_to_uniform_inplane_solver(theta_deg):
    lattice, g = _pillar_lattice_and_g(num_orders=3)
    omega = 2 * np.pi / WAVELENGTH
    kx, ky = _kx_ky(lattice, g, omega, theta_deg)

    tensor = np.array([[2.25, 0.3, 0], [0.3, 4.0, 0], [0, 0, 3.1]], dtype=complex)
    material = Material.from_permittivity_tensor("uniform_tensor", tensor)
    # A pattern whose one shape exactly equals the background is spatially
    # uniform -- every non-DC Fourier coefficient must vanish.
    pattern = Pattern(
        background=material,
        shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=material)],
    )

    exx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 0)
    exy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 1)
    eyx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 0)
    eyy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 1)
    ezz = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 2, 2)
    modes_patterned = solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, exx, exy, eyx, eyy, ezz)
    modes_uniform = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)

    q_sq_patterned = np.sort_complex(modes_patterned.q**2)
    q_sq_uniform = np.sort_complex(modes_uniform.q**2)
    assert q_sq_patterned == pytest.approx(q_sq_uniform, abs=1e-9)


# ---------------------------------------------------------------------------
# Physical invariant: energy conservation for a genuinely patterned,
# Hermitian (lossless) anisotropic case.
# ---------------------------------------------------------------------------


def _anisotropic_pillar_simulation(num_orders=5):
    air = Material("air", 1.0)
    tensor = np.array([[2.25, 0.3, 0], [0.3, 4.0, 0], [0, 0, 3.1]], dtype=complex)
    si_aniso = Material.from_permittivity_tensor("si_aniso", tensor)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(
        background=air,
        shapes=[
            Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=si_aniso),
            Rectangle(center=(0.15, 0.15), halfwidth=(0.08, 0.12), material=si_aniso),
        ],
    )
    layer = Layer("pillar_aniso", 0.46, pattern=pattern)
    return Simulation(lattice, [layer], num_orders=num_orders, incidence=air, transmission=air)


@pytest.mark.parametrize("theta_deg", [0.0, 20.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)])
def test_energy_conservation_anisotropic_patterned_layer(theta_deg, s_amp, p_amp):
    sim = _anisotropic_pillar_simulation(num_orders=5)
    excitation = PlaneWaveExcitation(WAVELENGTH, np.radians(theta_deg), 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


def test_anisotropic_pillar_runs_end_to_end():
    sim = _anisotropic_pillar_simulation(num_orders=5)
    excitation = PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert 0.0 <= result.reflectance() <= 1.0
    assert 0.0 <= result.transmittance() <= 1.0


# ---------------------------------------------------------------------------
# Guard: longitudinal coupling anywhere in the pattern still raises
# ---------------------------------------------------------------------------


def test_longitudinal_coupling_in_pattern_raises_not_implemented():
    air = Material("air", 1.0)
    tensor = np.array([[2.25, 0.0, 0.1], [0.0, 4.0, 0.0], [0.1, 0.0, 3.1]], dtype=complex)
    material = Material.from_permittivity_tensor("longitudinal", tensor)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=air, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=material)])
    layer = Layer("longitudinal_pillar", 0.46, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=5, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    with pytest.raises(NotImplementedError, match="1.5"):
        sim.solve(excitation)
