"""Category 13 target 13.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`Simulation._eigenmode_cache`, implementing the design Category 12 target
12.3 flagged (`design.md`'s "Linear-Algebra Baseline & Factorization-Reuse
Design") but deliberately left unimplemented there. Two things are
checked, per the same "validate the optimized path against the
unoptimized one before trusting it" discipline `decisions.md` ADR-016
already established for the Toeplitz-matrix cache:

1. **Equivalence**: results with the cache populated match results with
   the cache forcibly cleared before every layer (forcing recomputation).
2. **Cache-hit behavior**: the two scenarios the cache is actually
   designed for (a fixed-wavelength/angle polarization sweep, and a
   fixed-wavelength/angle thickness sweep) populate exactly one cache
   entry per layer, not one per sweep point.
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
from sougata_solver.sweep import sweep_polarization, sweep_thickness, sweep_wavelength

PERIOD = 0.7
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


def _pillar_simulation(num_orders: int = 9) -> Simulation:
    lattice = Lattice((PERIOD, 0.0), (0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.2 * PERIOD, material=SI)])
    layer = Layer("pillar", 0.3, pattern=pattern)
    return Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR)


def test_cached_result_matches_forced_uncached_recomputation():
    sim_cached = _pillar_simulation()
    excitation = PlaneWaveExcitation(0.6, math.radians(10.0), 0.0, s_amplitude=0.7, p_amplitude=0.5)
    result_cached = sim_cached.solve(excitation)

    sim_forced_uncached = _pillar_simulation()
    real_cached_layer_eigenmodes = sim_forced_uncached._cached_layer_eigenmodes

    def _always_miss(layer, g, wavelength, omega, kx, ky, is_1d):
        sim_forced_uncached._eigenmode_cache.clear()
        return real_cached_layer_eigenmodes(layer, g, wavelength, omega, kx, ky, is_1d)

    sim_forced_uncached._cached_layer_eigenmodes = _always_miss
    result_uncached = sim_forced_uncached.solve(excitation)

    assert result_cached.reflectance() == pytest.approx(result_uncached.reflectance(), abs=1e-12)
    assert result_cached.transmittance() == pytest.approx(result_uncached.transmittance(), abs=1e-12)
    np.testing.assert_allclose(result_cached.a_transmitted, result_uncached.a_transmitted, atol=1e-12)
    np.testing.assert_allclose(result_cached.b_reflected, result_uncached.b_reflected, atol=1e-12)


def test_eigenmode_cache_reused_across_polarization_sweep_at_fixed_wavelength_angle():
    """`sweep_polarization` never changes `omega`/`kx`/`ky` (only the
    incident amplitude vector `a0`, computed downstream of every
    eigenmode solve) -- exactly one eigenmode-cache entry per layer for
    the whole sweep, regardless of how many Jones states are swept."""
    sim = _pillar_simulation()
    jones_states = [(1.0, 0.0), (0.0, 1.0), (1.0 / math.sqrt(2), 1.0 / math.sqrt(2)), (1.0, 1.0j)]
    sweep_polarization(sim, wavelength=0.6, theta=math.radians(10.0), phi=0.0, jones_states=jones_states)

    # 3 layers (incidence, pillar, transmission) -> 3 eigenmode-cache entries.
    assert len(sim._eigenmode_cache) == 3


def test_eigenmode_cache_reused_across_thickness_sweep_at_fixed_wavelength_angle():
    """`sweep_thickness` mutates a `Layer`'s `thickness` in place -- the
    cache key must stay valid throughout (same object identity, and
    eigenmode solves never depend on thickness), giving exactly one
    cache entry per layer across the whole sweep."""
    sim = _pillar_simulation()
    excitation = PlaneWaveExcitation(0.6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    sweep_thickness(sim, "pillar", [0.1, 0.2, 0.3, 0.4], excitation)

    assert len(sim._eigenmode_cache) == 3


def test_eigenmode_cache_is_not_reused_across_a_wavelength_sweep():
    """A genuine wavelength sweep changes `omega` (and therefore `kx`/`ky`
    via `k_parallel`) at every point -- confirms the cache correctly does
    *not* claim reuse it cannot safely provide, unlike the Toeplitz cache
    (which is wavelength-keyed but angle-independent)."""
    sim = _pillar_simulation()
    wavelengths = [0.55, 0.6, 0.65, 0.72]
    sweep_wavelength(sim, wavelengths, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    # 3 layers x 4 distinct wavelengths -> 12 distinct cache entries, no reuse.
    assert len(sim._eigenmode_cache) == 3 * len(wavelengths)
