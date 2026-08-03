"""Category 1 target 1.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): uniform
in-plane-coupled anisotropic layers (`eps_xx, eps_xy, eps_yx, eps_yy,
eps_zz`).

Tiers enforced here, per `rules.md` Testing Requirements:
- independent oracle: `solve_layer_eigenmodes_uniform_inplane`'s `q^2`
  eigenvalues cross-checked against `tests/oracles/rcwa_anisotropic_inplane_jl.py`
  (hand-transcribed from `RigorousCoupledWaveAnalysis.jl`, a structurally
  different derivation than this project's `Epsilon2 @ kp - coupling`).
- regression: `eps_xy=eps_yx=0` reduces to
  `solve_layer_eigenmodes_uniform_diagonal`'s (target 1.3) result.
- physical invariant: energy conservation for a genuinely in-plane-coupled
  slab through a full `Simulation.solve`.
- guard: a longitudinally-coupled tensor (nonzero eps_xz/eps_yz/eps_zx/eps_zy)
  still raises `NotImplementedError` naming target 1.5, per `rules.md` AI
  Coding Rule 2.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from oracles.rcwa_anisotropic_inplane_jl import eigenoperator_eigenvalues_inplane

from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform_diagonal, solve_layer_eigenmodes_uniform_inplane
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6


# ---------------------------------------------------------------------------
# Independent oracle: RigorousCoupledWaveAnalysis.jl AnisotropicLayer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 15.0, 35.0])
def test_inplane_eigenvalues_match_rcwa_jl_oracle(theta_deg):
    omega = 2 * np.pi / WAVELENGTH
    num_ord = 3
    n = 2 * num_ord + 1
    period = 1e-6
    kx0 = omega * math.sin(math.radians(theta_deg))
    kx = kx0 + 2 * np.pi * np.arange(-num_ord, num_ord + 1) / period
    ky = 0.3 * omega * np.ones(n)  # nonzero, non-classical mounting to exercise coupling

    eps_xx, eps_xy, eps_yx, eps_yy, eps_zz = 2.25 + 0.02j, 0.35 - 0.05j, 0.22 + 0.03j, 4.0, 3.1

    modes = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, eps_xx, eps_xy, eps_yx, eps_yy, eps_zz)
    ours = np.sort_complex(modes.q**2)

    # See rcwa_anisotropic_inplane_jl.py's docstring for why kx/ky are
    # swapped and eps_xy/eps_yx negated when calling the oracle -- an
    # empirically-determined convention reconciliation, not a formula
    # difference.
    oracle = np.sort_complex(eigenoperator_eigenvalues_inplane(omega, ky, kx, eps_xx, -eps_xy, -eps_yx, eps_yy, eps_zz))

    # q^2 scales as omega^2 (~1e14 at optical wavelengths), so a relative
    # tolerance is the meaningful comparison here (matches Phase 4a's
    # oracle test's intent -- that test used abs=1e-6 at omega~O(10) only
    # because its wavelength=1.0 test scale made abs and rel coincide).
    assert ours == pytest.approx(oracle, rel=1e-8)


# ---------------------------------------------------------------------------
# Regression: eps_xy=eps_yx=0 reduces to target 1.3's diagonal solver
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 25.0])
def test_zero_off_diagonal_reduces_to_diagonal_solver(theta_deg):
    omega = 2 * np.pi / WAVELENGTH
    num_ord = 2
    n = 2 * num_ord + 1
    period = 1e-6
    kx0 = omega * math.sin(math.radians(theta_deg))
    kx = kx0 + 2 * np.pi * np.arange(-num_ord, num_ord + 1) / period
    ky = np.zeros(n)

    eps_xx, eps_yy, eps_zz = 2.25, 4.0, 3.1

    modes_inplane = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, eps_xx, 0.0, 0.0, eps_yy, eps_zz)
    modes_diag = solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx, eps_yy, eps_zz)

    q_sq_inplane = np.sort_complex(modes_inplane.q**2)
    q_sq_diag = np.sort_complex(modes_diag.q**2)
    assert q_sq_inplane == pytest.approx(q_sq_diag, abs=1e-10)


# ---------------------------------------------------------------------------
# System-level: energy conservation for a genuinely in-plane-coupled slab
# ---------------------------------------------------------------------------


def _inplane_simulation(eps_xx, eps_xy, eps_yx, eps_yy, eps_zz, thickness):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    tensor = np.array([[eps_xx, eps_xy, 0], [eps_yx, eps_yy, 0], [0, 0, eps_zz]], dtype=complex)
    layer = Layer("slab", thickness, material=Material.from_permittivity_tensor("inplane", tensor))
    return Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)


@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 40.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.7, 0.7)])
def test_energy_conservation_inplane_coupled_slab(theta_deg, s_amp, p_amp):
    # Hermitian in-plane block (eps_xy = conj(eps_yx), real here) is the
    # lossless-reciprocal case -- a non-Hermitian tensor (eps_xy != conj(eps_yx))
    # represents a gain/loss (non-reciprocal-power) medium, for which R+T=1
    # is not physically expected; that distinction is what this test's
    # material choice is deliberately respecting.
    sim = _inplane_simulation(2.25, 0.3, 0.3, 4.0, 3.1, 0.25e-6)
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(theta_deg), 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Guard: longitudinal coupling still raises, naming target 1.5
# ---------------------------------------------------------------------------


def test_longitudinal_coupling_raises_not_implemented():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    tensor = np.array([[2.25, 0.0, 0.1], [0.0, 4.0, 0.0], [0.1, 0.0, 3.1]], dtype=complex)
    layer = Layer("slab", 0.2e-6, material=Material.from_permittivity_tensor("longitudinal", tensor))
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    with pytest.raises(NotImplementedError, match="1.5"):
        sim.solve(excitation)
