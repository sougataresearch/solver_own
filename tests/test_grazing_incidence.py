"""Category 6 target 6.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): grazing-
incidence boundary. Defines and tests the supported near-grazing angle
range for a uniform (unpatterned) stack, characterized directly (not
assumed) this session:

- The solver stays finite, well-conditioned, and energy-conserving
  (`R+T=1` to `1e-8`+) for `theta` arbitrarily close to `90 deg` --
  confirmed here up to `theta=89.999 deg`.
- Exactly `theta=90 deg` raises a plain `ValueError` (propagated from
  `scipy.linalg.lu_factor`'s internal finiteness check inside
  `smatrix.interface_smatrix`, not caught/swallowed anywhere in this
  project, per `design.md`'s "no broad except" convention), not a silently
  wrong or `NaN` answer.

**Root cause, same class of issue as the Rayleigh-threshold `NaN`
documented in `troubleshooting.md`**: at exact grazing incidence, the
incidence half-space's own zeroth-order `q` is exactly `0` --
`interface_smatrix`'s `kp @ phi / q[None, :]` divides by zero. The reason
it happens at exactly `theta=90 deg` and not just "very close" is a
genuine floating-point coincidence, confirmed directly: `math.sin(math.pi/2)
== 1.0` exactly in float64 (the true value `1.0` is also the nearest
representable float), so for an `n=1` (air) incidence medium,
`kx0 = omega*sin(theta) == omega` exactly, giving `q_sq = omega^2 - kx0^2
== 0.0` exactly, not merely a tiny nonzero residual. For a different
incidence index this exact cancellation would not occur at exactly
`theta=90 deg` in floating point, but the same singularity is still
reached in the true-mathematical limit `theta -> 90 deg` regardless.

**Supported range, as documented here**: `theta < 90 deg` (any incidence
medium/multilayer combination) is supported; `theta = 90 deg` is not, and
raises `ValueError` rather than returning a plausible-but-wrong result.
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


@pytest.mark.parametrize("theta_deg", [80.0, 85.0, 89.0, 89.9, 89.99, 89.999])
def test_near_grazing_incidence_stays_finite_and_energy_conserving(theta_deg):
    sim = _build_sim()
    theta = math.radians(theta_deg)
    result = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0))
    r, t = result.reflectance(), result.transmittance()
    assert np.isfinite(r) and np.isfinite(t)
    assert (r + t) == pytest.approx(1.0, abs=1e-6)


def test_near_grazing_incidence_reflectance_increases_monotonically_toward_unity():
    """Physical sanity check alongside the finiteness check above: as
    `theta -> 90 deg`, `R -> 1` (all incident power reflects at true
    grazing incidence) -- a real, monotone trend, not just "doesn't
    crash"."""
    sim = _build_sim()
    angles_deg = [10.0, 30.0, 60.0, 80.0, 89.0, 89.9, 89.99]
    reflectances = [
        sim.solve(PlaneWaveExcitation(WAVELENGTH, math.radians(a), 0.0, s_amplitude=1.0, p_amplitude=0.0)).reflectance()
        for a in angles_deg
    ]
    assert all(r1 <= r2 + 1e-12 for r1, r2 in zip(reflectances, reflectances[1:]))
    assert reflectances[-1] > 0.999


def test_exact_grazing_incidence_raises_value_error_not_silently_wrong_answer():
    sim = _build_sim()
    with pytest.raises(ValueError):
        sim.solve(PlaneWaveExcitation(WAVELENGTH, math.radians(90.0), 0.0, s_amplitude=1.0, p_amplitude=0.0))


def test_grazing_incidence_boundary_root_cause_is_exact_zero_q():
    """Confirms the docstring's floating-point-coincidence claim directly,
    rather than just asserting the ValueError and leaving the "why" as
    prose: at theta=90deg exactly, math.sin is exactly 1.0 and the
    incidence-medium zeroth-order q is exactly 0.0 (not just small)."""
    assert math.sin(math.radians(90.0)) == 1.0

    from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform

    omega = 2 * math.pi / WAVELENGTH
    kx0 = omega * 1.0 * math.sin(math.radians(90.0)) * math.cos(0.0)
    ky0 = omega * 1.0 * math.sin(math.radians(90.0)) * math.sin(0.0)
    modes = solve_layer_eigenmodes_uniform(omega, np.array([kx0]), np.array([ky0]), eps=1.0)
    assert modes.q[0] == 0.0
