"""Category 6 target 6.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): bottom
(reverse-side) illumination. See `decisions.md` ADR-014 for the full
decision -- summary: already achievable via the existing public
`Simulation` constructor (`list(reversed(layers))` + swapped
`incidence`/`transmission`), no new API needed. This file is the
regression guard for that recipe's correctness (a permanent test, not a
one-off investigation script) -- the Stokes transmittance-reciprocity
relation for a lossless reciprocal medium at normal incidence.
"""

from __future__ import annotations

import math

import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6


def _asymmetric_stack():
    l1 = Layer("L1", WAVELENGTH / (4 * 1.46), material=Material("L1", 1.46**2))
    l2 = Layer("L2", WAVELENGTH / (4 * 2.35), material=Material("L2", 2.35**2))
    return l1, l2


def test_reversed_stack_transmittance_matches_forward_at_normal_incidence():
    """Stokes reciprocity: for a lossless reciprocal medium, transmittance
    at normal incidence is the same illuminated from either side, even
    though the incidence/transmission indices differ (air vs. glass)."""
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    l1, l2 = _asymmetric_stack()

    sim_top = Simulation(lattice, [l1, l2], num_orders=1, incidence=air, transmission=glass)
    result_top = sim_top.solve(PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0))

    sim_bottom = Simulation(lattice, [l2, l1], num_orders=1, incidence=glass, transmission=air)
    result_bottom = sim_bottom.solve(PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0))

    assert result_bottom.transmittance() == pytest.approx(result_top.transmittance(), abs=1e-10)


def test_reversed_stack_is_independently_energy_conserving():
    """The reversed-stack simulation is not just numerically close to the
    forward one -- it is itself a fully well-posed, energy-conserving
    simulation (R+T=1), at both normal and oblique incidence."""
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    l1, l2 = _asymmetric_stack()
    sim_bottom = Simulation(lattice, [l2, l1], num_orders=1, incidence=glass, transmission=air)

    for theta_deg in [0.0, 30.0]:
        result = sim_bottom.solve(
            PlaneWaveExcitation(WAVELENGTH, math.radians(theta_deg), 0.0, s_amplitude=1.0, p_amplitude=0.0)
        )
        assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


def test_reversed_stack_reflectance_genuinely_differs_at_oblique_incidence():
    """Honest counter-check: reflectance is *not* claimed direction-
    independent in general -- confirms the forward/reversed R values
    genuinely differ at oblique incidence with mismatched media (the
    physically expected result for an asymmetric stack), so a future
    change that made them numerically equal would itself be suspicious,
    not a sign of extra correctness."""
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    l1, l2 = _asymmetric_stack()
    theta = math.radians(30.0)

    sim_top = Simulation(lattice, [l1, l2], num_orders=1, incidence=air, transmission=glass)
    result_top = sim_top.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0))

    sim_bottom = Simulation(lattice, [l2, l1], num_orders=1, incidence=glass, transmission=air)
    result_bottom = sim_bottom.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0))

    assert result_top.reflectance() != pytest.approx(result_bottom.reflectance(), abs=1e-4)
