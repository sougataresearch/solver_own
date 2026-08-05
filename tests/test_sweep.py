"""Category 8 targets 8.1-8.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`sweep.SweepResult` and the wavelength/angle/polarization/thickness sweep
functions. The category's own exit criterion is the guiding test for every
function here: **each sweep is equivalent to repeated scalar
`Simulation.solve()` calls** -- every test below confirms this by building
the same result two ways (via the sweep function, and via a manual
per-point loop) and checking they match, not just that the sweep function
"runs".
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import (
    SweepResult,
    sweep_phi,
    sweep_polarization,
    sweep_theta,
    sweep_thickness,
    sweep_wavelength,
)

PERIOD = 0.7e-6
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


def _pillar_simulation(num_orders: int = 9) -> Simulation:
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.2 * PERIOD, material=SI)])
    layer = Layer("pillar", 0.3e-6, pattern=pattern)
    return Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR)


# ---------------------------------------------------------------------------
# 8.1 SweepResult container
# ---------------------------------------------------------------------------


def test_sweep_result_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        SweepResult("wavelength", "m", [1.0, 2.0], [None])


def test_sweep_result_rejects_empty_sweep():
    with pytest.raises(ValueError, match="at least one point"):
        SweepResult("wavelength", "m", [], [])


def test_sweep_result_reflectance_transmittance_arrays():
    sim = _pillar_simulation()
    wavelengths = [0.55e-6, 0.6e-6, 0.65e-6]
    sweep = sweep_wavelength(sim, wavelengths, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    r = sweep.reflectance()
    t = sweep.transmittance()
    assert r.shape == (3,)
    assert t.shape == (3,)
    for i, wavelength in enumerate(wavelengths):
        expected = sim.solve(PlaneWaveExcitation(wavelength, 0.0, 0.0, 1.0, 0.0))
        assert r[i] == pytest.approx(expected.reflectance())
        assert t[i] == pytest.approx(expected.transmittance())


# ---------------------------------------------------------------------------
# 8.2 Wavelength sweep
# ---------------------------------------------------------------------------


def test_sweep_wavelength_matches_manual_scalar_loop():
    sim = _pillar_simulation()
    wavelengths = np.linspace(0.5e-6, 0.9e-6, 9)
    sweep = sweep_wavelength(sim, wavelengths, theta=math.radians(10.0), phi=math.radians(5.0), s_amplitude=0.6, p_amplitude=0.8)

    manual_r = []
    manual_t = []
    for wavelength in wavelengths:
        result = sim.solve(PlaneWaveExcitation(wavelength, math.radians(10.0), math.radians(5.0), 0.6, 0.8))
        manual_r.append(result.reflectance())
        manual_t.append(result.transmittance())

    np.testing.assert_allclose(sweep.reflectance(), manual_r)
    np.testing.assert_allclose(sweep.transmittance(), manual_t)
    assert sweep.parameter_name == "wavelength"
    assert sweep.parameter_unit == "m"
    assert list(sweep.parameter_values) == list(wavelengths)
    assert sweep.metadata["num_orders"] == sim.num_orders


# ---------------------------------------------------------------------------
# 8.3 Angle sweeps
# ---------------------------------------------------------------------------


def test_sweep_theta_matches_manual_scalar_loop():
    sim = _pillar_simulation()
    thetas = [math.radians(a) for a in (0, 10, 20, 30)]
    sweep = sweep_theta(sim, wavelength=0.6e-6, thetas=thetas, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    manual_r = [sim.solve(PlaneWaveExcitation(0.6e-6, theta, 0.0, 1.0, 0.0)).reflectance() for theta in thetas]
    np.testing.assert_allclose(sweep.reflectance(), manual_r)
    assert sweep.parameter_name == "theta"
    assert sweep.parameter_unit == "rad"


def test_sweep_theta_reuses_toeplitz_cache_across_the_sweep():
    """Same real benefit `decisions.md` ADR-016 measured for the Toeplitz
    cache: a fixed-wavelength angle sweep should populate exactly one
    cache entry, not one per angle."""
    sim = _pillar_simulation()
    thetas = [math.radians(a) for a in (0, 5, 10, 15, 20)]
    sweep_theta(sim, wavelength=0.6e-6, thetas=thetas, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    assert len(sim._toeplitz_cache) == 1


def test_sweep_phi_matches_manual_scalar_loop():
    sim = _pillar_simulation()
    phis = [math.radians(a) for a in (0, 30, 60, 90)]
    sweep = sweep_phi(sim, wavelength=0.6e-6, theta=math.radians(15.0), phis=phis, s_amplitude=0.7, p_amplitude=0.5)

    manual_t = [sim.solve(PlaneWaveExcitation(0.6e-6, math.radians(15.0), phi, 0.7, 0.5)).transmittance() for phi in phis]
    np.testing.assert_allclose(sweep.transmittance(), manual_t)
    assert sweep.parameter_name == "phi"


# ---------------------------------------------------------------------------
# 8.4 Polarization sweep
# ---------------------------------------------------------------------------


def test_sweep_polarization_matches_manual_scalar_loop():
    sim = _pillar_simulation()
    # TE, TM, and 45-degree linear (CONVENTIONS.md's worked examples table).
    jones_states = [(1.0, 0.0), (0.0, 1.0), (1.0 / math.sqrt(2), 1.0 / math.sqrt(2))]
    sweep = sweep_polarization(sim, wavelength=0.6e-6, theta=math.radians(10.0), phi=0.0, jones_states=jones_states)

    manual_r = [
        sim.solve(PlaneWaveExcitation(0.6e-6, math.radians(10.0), 0.0, s, p)).reflectance()
        for s, p in jones_states
    ]
    np.testing.assert_allclose(sweep.reflectance(), manual_r)
    assert sweep.parameter_name == "jones_state"
    assert list(sweep.parameter_values) == jones_states


def test_sweep_polarization_rejects_empty_list():
    sim = _pillar_simulation()
    with pytest.raises(ValueError, match="non-empty"):
        sweep_polarization(sim, wavelength=0.6e-6, theta=0.0, phi=0.0, jones_states=[])


# ---------------------------------------------------------------------------
# 8.5 Thickness sweep
# ---------------------------------------------------------------------------


def test_sweep_thickness_matches_manual_reconstruction():
    """No library-level "rebuild with a different thickness" exists, so the
    manual reference here rebuilds a fresh `Simulation` per thickness --
    still an independent path from `sweep_thickness`'s in-place mutation,
    confirming both reach the same physical answer."""
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.2 * PERIOD, material=SI)])
    thicknesses = [0.1e-6, 0.2e-6, 0.3e-6, 0.4e-6]
    excitation = PlaneWaveExcitation(0.6e-6, math.radians(10.0), 0.0, 1.0, 0.0)

    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    sweep = sweep_thickness(sim, "pillar", thicknesses, excitation)

    manual_r = []
    for t in thicknesses:
        fresh_sim = Simulation(lattice, [Layer("pillar", t, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
        manual_r.append(fresh_sim.solve(excitation).reflectance())

    np.testing.assert_allclose(sweep.reflectance(), manual_r, atol=1e-12)
    assert sweep.parameter_name == "thickness[pillar]"


def test_sweep_thickness_restores_original_thickness_after_sweep():
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.2 * PERIOD, material=SI)])
    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, 1.0, 0.0)

    sweep_thickness(sim, "pillar", [0.1e-6, 0.2e-6], excitation)

    pillar_layer = [layer for layer in sim.layer_stack if layer.name == "pillar"][0]
    assert pillar_layer.thickness == pytest.approx(0.3e-6)


def test_sweep_thickness_rejects_unknown_layer_name():
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    sim = Simulation(lattice, [Layer("l1", 0.3e-6, material=SI)], num_orders=1, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="no layer named"):
        sweep_thickness(sim, "nonexistent", [0.1e-6], excitation)


@pytest.mark.parametrize("bad_thickness", [0.0, -1.0, float("nan"), float("inf")])
def test_sweep_thickness_rejects_invalid_values(bad_thickness):
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    sim = Simulation(lattice, [Layer("l1", 0.3e-6, material=SI)], num_orders=1, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="finite and > 0"):
        sweep_thickness(sim, "l1", [0.1e-6, bad_thickness], excitation)
