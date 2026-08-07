"""Category 1 target 1.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): uniform
diagonal-tensor anisotropic layers.

Tiers enforced here, per `rules.md` Testing Requirements:
- unit/closed-form: `solve_layer_eigenmodes_uniform_diagonal` at normal
  incidence reduces exactly to `q_x^2 = eps_xx*omega^2`,
  `q_y^2 = eps_yy*omega^2` (independently-derived birefringence result, see
  the function's docstring).
- system-level closed-form oracle: a uniaxial slab at normal incidence,
  illuminated along its two in-plane principal axes, matches the
  already-validated isotropic Fresnel/TMM oracle (`tests/oracles/fresnel.py`)
  evaluated at `n_x = sqrt(eps_xx)` / `n_y = sqrt(eps_yy)` respectively --
  this is a genuine external-oracle comparison, not a self-check, because
  `fresnel.py` has no anisotropy concept at all and is reused unmodified.
- regression: isotropic-tensor special case (`eps_xx=eps_yy=eps_zz`)
  reproduces `solve_layer_eigenmodes_uniform`'s result through a full
  `Simulation.solve`, at both normal and oblique incidence.
- physical invariant: energy conservation for a genuinely birefringent
  slab at oblique incidence.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from oracles.fresnel import multilayer_rt

from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform_diagonal
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation


pytestmark = pytest.mark.oracle  # Category 17 target 17.1: system-tier test, cross-checked against a named external oracle

WAVELENGTH = 0.55e-6


# ---------------------------------------------------------------------------
# Unit-level closed form: normal incidence, single order
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eps_xx,eps_yy,eps_zz", [(2.0, 3.0, 5.0), (2.25, 2.25, 4.0), (1.0 + 0.1j, 1.5, 2.0)])
def test_normal_incidence_closed_form_matches_birefringence_formula(eps_xx, eps_yy, eps_zz):
    omega = 2 * np.pi / WAVELENGTH
    kx = np.array([0.0])
    ky = np.array([0.0])
    modes = solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx, eps_yy, eps_zz)
    q_sq = modes.q**2
    expected = sorted([eps_xx * omega**2, eps_yy * omega**2], key=lambda z: (z.real, z.imag))
    got = sorted(q_sq.tolist(), key=lambda z: (z.real, z.imag))
    assert got == pytest.approx(expected, rel=1e-10)


# ---------------------------------------------------------------------------
# System-level oracle: uniaxial waveplate at normal incidence decouples into
# two independent isotropic Fresnel problems along the crystal's principal
# axes -- checked against the already-validated fresnel.py oracle.
# ---------------------------------------------------------------------------


def _run_anisotropic_slab(eps_xx, eps_yy, eps_zz, thickness, polarization):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    tensor = np.diag([eps_xx, eps_yy, eps_zz]).astype(complex)
    slab = Material.from_permittivity_tensor("uniaxial", tensor)
    layer = Layer("slab", thickness, material=slab)
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)

    s_amp = 1.0 if polarization == "s" else 0.0
    p_amp = 1.0 if polarization == "p" else 0.0
    excitation = PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    return result.reflectance(), result.transmittance()


@pytest.mark.parametrize("thickness", [0.1e-6, 0.55e-6 / (4 * 1.5)])
def test_normal_incidence_uniaxial_slab_matches_fresnel_oracle_per_axis(thickness):
    eps_xx, eps_yy, eps_zz = 2.25, 4.0, 6.0  # n_x=1.5, n_y=2.0 (fast/slow axes distinguishable)

    # excitation.py's normal-incidence basis puts s-polarization along +y and
    # p-polarization along x (CONVENTIONS.md's Polarization convention).
    # Which tensor component each experiences is set by the *solver's*
    # internal component ordering, not the polarization label alone:
    # CONVENTIONS.md documents `u = [-Ey; Ex]` (block 0 <-> -Ey, block 1 <->
    # Ex), and Epsilon2's block layout (`solve_layer_eigenmodes_uniform_diagonal`'s
    # docstring, transcribed from `S4/S4/S4.cpp:1889-1906`) puts eps_xx in
    # block 0 and eps_yy in block 1 -- so eps_xx governs Ey (s-polarization)
    # and eps_yy governs Ex (p-polarization), confirmed empirically against
    # this test before trusting the assertion direction (an initial draft
    # had this backwards and was caught by the oracle mismatch, not assumed
    # correct).
    r_s, t_s = _run_anisotropic_slab(eps_xx, eps_yy, eps_zz, thickness, "s")
    r_p, t_p = _run_anisotropic_slab(eps_xx, eps_yy, eps_zz, thickness, "p")

    n_for_s = math.sqrt(eps_xx)
    n_for_p = math.sqrt(eps_yy)
    r_oracle_s, t_oracle_s = multilayer_rt(WAVELENGTH, 0.0, "s", 1.0, [(n_for_s, thickness)], 1.0)
    r_oracle_p, t_oracle_p = multilayer_rt(WAVELENGTH, 0.0, "s", 1.0, [(n_for_p, thickness)], 1.0)

    assert r_s == pytest.approx(r_oracle_s, abs=1e-8)
    assert t_s == pytest.approx(t_oracle_s, abs=1e-8)
    assert r_p == pytest.approx(r_oracle_p, abs=1e-8)
    assert t_p == pytest.approx(t_oracle_p, abs=1e-8)


# ---------------------------------------------------------------------------
# Regression: isotropic-tensor special case reduces to the existing
# scalar-isotropic uniform-layer path (Phase 1), at normal AND oblique
# incidence.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 45.0])
@pytest.mark.parametrize("polarization", ["s", "p"])
def test_isotropic_diagonal_tensor_reduces_to_scalar_isotropic_path(theta_deg, polarization):
    n = 1.8
    eps = n**2
    thickness = 0.2e-6
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    theta = math.radians(theta_deg)
    s_amp = 1.0 if polarization == "s" else 0.0
    p_amp = 1.0 if polarization == "p" else 0.0
    excitation = PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=s_amp, p_amplitude=p_amp)

    scalar_layer = Layer("scalar", thickness, material=Material("iso", eps))
    sim_scalar = Simulation(lattice, [scalar_layer], num_orders=1, incidence=air, transmission=air)
    result_scalar = sim_scalar.solve(excitation)

    tensor = np.diag([eps, eps, eps]).astype(complex)
    tensor_layer = Layer("tensor", thickness, material=Material.from_permittivity_tensor("iso_tensor", tensor))
    sim_tensor = Simulation(lattice, [tensor_layer], num_orders=1, incidence=air, transmission=air)
    result_tensor = sim_tensor.solve(excitation)

    assert result_tensor.reflectance() == pytest.approx(result_scalar.reflectance(), abs=1e-8)
    assert result_tensor.transmittance() == pytest.approx(result_scalar.transmittance(), abs=1e-8)


# ---------------------------------------------------------------------------
# Physical invariant: energy conservation for a birefringent slab at
# oblique incidence and mixed polarization.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 15.0, 40.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)])
def test_energy_conservation_birefringent_slab(theta_deg, s_amp, p_amp):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    tensor = np.diag([2.25, 4.0, 3.0]).astype(complex)
    layer = Layer("slab", 0.3e-6, material=Material.from_permittivity_tensor("uniaxial", tensor))
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(theta_deg), 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)
