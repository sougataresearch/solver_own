"""Category 6 targets 6.2/6.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
polarization-state regressions against known symmetry expectations, per
`CONVENTIONS.md`'s "Worked polarization examples" table.

**6.2 (normal incidence)**: for an isotropic (non-birefringent), laterally-
uniform stack illuminated at `theta=0`, the physics has full rotational
symmetry about the propagation axis -- there is no preferred in-plane
direction, so `R`/`T` cannot depend on which particular polarization state
(TE, TM, linear at any angle, RCP, LCP, elliptical) carries a given total
incident power. This is a strong, oracle-independent invariant (stronger
than energy conservation alone): every state in `CONVENTIONS.md`'s table,
at fixed `|s_amplitude|^2+|p_amplitude|^2=1`, must give the *same* `R`
and the *same* `T`, not just `R+T=1` each. Verified numerically before
writing this test (not assumed), see the module docstring's commit history
for the exact confirmation run.

**6.3 (oblique incidence)**: two independent checks. (a) Energy
conservation for mixed/elliptical polarization at oblique incidence and
nonzero azimuth `phi` -- weaker than 6.2's exact-state-independence
(oblique incidence *does* break the TE/TM degeneracy -- s and p see
different effective admittances -- so `R`/`T` legitimately differ between
states here). (b) Azimuthal-rotation invariance: a laterally-uniform stack
has no preferred azimuth either, so at *fixed* `theta` and a *fixed*
`(s_amplitude, p_amplitude)` pair (s/p are defined relative to the plane
of incidence, which rotates rigidly with `phi`), `R`/`T` must be
independent of `phi` itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6


def _build_sim():
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    layers = [Layer("film", 0.1e-6, material=Material("film", 2.0**2))]
    return Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)


# Every state normalized to unit total power (|s|^2 + |p|^2 == 1), per
# CONVENTIONS.md's "Worked polarization examples" table.
_UNIT_POWER_STATES = {
    "TE": (1.0, 0.0),
    "TM": (0.0, 1.0),
    "linear_45deg": (1 / np.sqrt(2), 1 / np.sqrt(2)),
    "linear_20deg": (math.cos(math.radians(20.0)), math.sin(math.radians(20.0))),
    "RCP": (1 / np.sqrt(2), 1j / np.sqrt(2)),
    "LCP": (1 / np.sqrt(2), -1j / np.sqrt(2)),
    "elliptical": (math.cos(0.3), math.sin(0.3) * np.exp(1j * 0.7)),
}


# ---------------------------------------------------------------------------
# 6.2 Normal-incidence polarization regression
# ---------------------------------------------------------------------------


def test_normal_incidence_rt_independent_of_polarization_state():
    sim = _build_sim()
    results = {
        name: sim.solve(PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=s, p_amplitude=p))
        for name, (s, p) in _UNIT_POWER_STATES.items()
    }
    reference_r = results["TE"].reflectance()
    reference_t = results["TE"].transmittance()
    for name, result in results.items():
        assert result.reflectance() == pytest.approx(reference_r, abs=1e-10), name
        assert result.transmittance() == pytest.approx(reference_t, abs=1e-10), name


def test_normal_incidence_rt_independent_of_polarization_state_multilayer():
    """Same invariant, a genuinely multilayer (not single-film) stack, so
    the check isn't accidentally trivial for a bare interface."""
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    layers = [
        Layer("AR1", WAVELENGTH / (4 * 1.46), material=Material("AR1", 1.46**2)),
        Layer("AR2", WAVELENGTH / (4 * 2.35), material=Material("AR2", 2.35**2)),
    ]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)
    results = {
        name: sim.solve(PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=s, p_amplitude=p))
        for name, (s, p) in _UNIT_POWER_STATES.items()
    }
    reference_r = results["TE"].reflectance()
    for name, result in results.items():
        assert result.reflectance() == pytest.approx(reference_r, abs=1e-10), name


# ---------------------------------------------------------------------------
# 6.3a Oblique-incidence energy conservation, mixed/elliptical polarization
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [10.0, 30.0, 60.0])
@pytest.mark.parametrize("phi_deg", [0.0, 37.0, 90.0, 200.0])
@pytest.mark.parametrize("state_name", list(_UNIT_POWER_STATES.keys()))
def test_oblique_incidence_energy_conservation_all_states_all_azimuths(theta_deg, phi_deg, state_name):
    sim = _build_sim()
    s, p = _UNIT_POWER_STATES[state_name]
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)
    result = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, phi, s_amplitude=s, p_amplitude=p))
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# 6.3b Azimuthal-rotation invariance (fixed theta, fixed s/p, varying phi)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [15.0, 30.0, 55.0])
def test_azimuthal_rotation_invariance_for_uniform_stack(theta_deg):
    sim = _build_sim()
    theta = math.radians(theta_deg)
    s, p = math.cos(0.4), math.sin(0.4) * np.exp(1j * 1.1)  # an arbitrary elliptical state

    reference = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=s, p_amplitude=p))
    for phi_deg in [37.0, 90.0, 180.0, 275.0]:
        result = sim.solve(
            PlaneWaveExcitation(WAVELENGTH, theta, math.radians(phi_deg), s_amplitude=s, p_amplitude=p)
        )
        assert result.reflectance() == pytest.approx(reference.reflectance(), abs=1e-9)
        assert result.transmittance() == pytest.approx(reference.transmittance(), abs=1e-9)
