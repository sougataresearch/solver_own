"""Category 14 targets 14.5/14.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
reciprocity tests for uniform (unpatterned) layer stacks. See
`decisions.md` ADR-025 for the full account, including two findings
verified numerically before being asserted here (not assumed):

1. The correct reciprocity comparison matches incidence angles via
   Snell's law (same transverse wavevector `kx`), **not** the same
   nominal `theta` reused for both the forward and reversed (materials-
   swapped) stack -- the naive comparison is wrong and grows arbitrarily
   bad with angle (`test_naive_same_theta_comparison_fails_at_oblique_incidence`
   pins this as a permanent regression guard).
2. Total transmittance reciprocity does **not** extend to patterned
   (diffractive) layers -- deliberately scoped out, not silently ignored
   (`test_patterned_layer_total_transmittance_is_not_reciprocal`).
"""

from __future__ import annotations

import math

import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6
N_AIR = 1.0
N_GLASS = 1.5
AIR = Material("air", N_AIR**2)
GLASS = Material("glass", N_GLASS**2)


def _snell_matched_angle(theta1: float, n1: float, n2: float) -> float:
    """Refraction angle in medium 2 for a ray incident at `theta1` in
    medium 1, `n1*sin(theta1) = n2*sin(theta2)` -- the correct basis for
    comparing forward/reversed transmittance (see `decisions.md` ADR-025)."""
    return math.asin(n1 * math.sin(theta1) / n2)


def _asymmetric_uniform_stack() -> tuple[Layer, Layer]:
    l1 = Layer("L1", WAVELENGTH / (4 * 1.46), material=Material("L1", 1.46**2))
    l2 = Layer("L2", WAVELENGTH / (4 * 2.35), material=Material("L2", 2.35**2))
    return l1, l2


def _transmittance(layers: list[Layer], incidence: Material, transmission: Material, theta: float) -> float:
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    sim = Simulation(lattice, layers, num_orders=1, incidence=incidence, transmission=transmission)
    result = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0))
    return result.transmittance()


@pytest.mark.parametrize("theta_air_deg", [0.0, 15.0, 30.0, 40.0])
def test_snell_matched_transmittance_is_reciprocal_lossless(theta_air_deg):
    l1, l2 = _asymmetric_uniform_stack()
    theta_air = math.radians(theta_air_deg)
    theta_glass = _snell_matched_angle(theta_air, N_AIR, N_GLASS)

    t_forward = _transmittance([l1, l2], AIR, GLASS, theta_air)
    t_reversed = _transmittance([l2, l1], GLASS, AIR, theta_glass)

    assert t_reversed == pytest.approx(t_forward, abs=1e-10)


@pytest.mark.parametrize("theta_air_deg", [0.0, 20.0, 35.0])
def test_snell_matched_transmittance_is_reciprocal_lossy(theta_air_deg):
    """Ordinary absorption (a lossy but still reciprocal medium) does not
    break reciprocity -- only a nonreciprocal (e.g. magnetized) medium
    would, and this project has no such material model."""
    lossy = Material("lossy", (2.0 + 0.3j) ** 2)
    theta_air = math.radians(theta_air_deg)
    theta_glass = _snell_matched_angle(theta_air, N_AIR, N_GLASS)

    t_forward = _transmittance([Layer("l1", 0.1e-6, material=lossy)], AIR, GLASS, theta_air)
    t_reversed = _transmittance([Layer("l1", 0.1e-6, material=lossy)], GLASS, AIR, theta_glass)

    assert t_reversed == pytest.approx(t_forward, abs=1e-10)


@pytest.mark.parametrize("theta_deg", [15.0, 30.0, 45.0])
def test_naive_same_theta_comparison_fails_at_oblique_incidence(theta_deg):
    """Negative control, pinning `decisions.md` ADR-025's first finding:
    reusing the *same* `theta` for both the forward and reversed stack
    (instead of the Snell's-law-matched angle) is the naive, incorrect
    reciprocity comparison -- confirmed to genuinely diverge from the
    forward result at oblique incidence, growing with angle."""
    l1, l2 = _asymmetric_uniform_stack()
    theta = math.radians(theta_deg)

    t_forward = _transmittance([l1, l2], AIR, GLASS, theta)
    t_reversed_naive = _transmittance([l2, l1], GLASS, AIR, theta)

    assert abs(t_forward - t_reversed_naive) > 1e-3


def test_patterned_layer_total_transmittance_is_not_reciprocal():
    """Documents the scope boundary from `decisions.md` ADR-025's second
    finding: total transmittance reciprocity, even at Snell-matched
    angles, does not hold for a diffractive (patterned) layer -- checked
    directly, not assumed to generalize from the uniform-layer case."""
    period = 0.7e-6
    si = Material("si", 3.48**2)
    lattice = Lattice1D(period)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.15e-6, material=si))
    layer = Layer("grating", 0.3e-6, pattern=pattern)

    excitation_wavelength = 1.0e-6
    sim_forward = Simulation(lattice, [layer], num_orders=9, incidence=AIR, transmission=GLASS)
    t_forward = sim_forward.solve(
        PlaneWaveExcitation(excitation_wavelength, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    ).transmittance()

    sim_reversed = Simulation(lattice, [layer], num_orders=9, incidence=GLASS, transmission=AIR)
    t_reversed = sim_reversed.solve(
        PlaneWaveExcitation(excitation_wavelength, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    ).transmittance()

    assert abs(t_forward - t_reversed) > 0.1
