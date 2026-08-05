"""Category 7 target 7.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): the
`Simulation._toeplitz_cache` implementation (`decisions.md` ADR-016,
`design.md`'s "Layer/Toeplitz Caching Design"). Two things are checked,
per that design's own requirement and `rules.md`'s "validate the optimized
path against the unoptimized one before trusting it":

1. **Equivalence**: a repeated-pattern stack solved normally (cache
   populated) gives numerically identical R/T to the same stack solved
   with every Toeplitz-matrix cache entry forcibly evicted before each
   layer is processed (forcing recomputation every time, i.e. what the
   uncached code path used to do).
2. **Cache-hit behavior**: a call-counting monkeypatch of
   `fourier_factorization.toeplitz_matrix` confirms `N` layers sharing the
   *same* `Pattern` object trigger exactly one real computation, not `N` --
   the actual point of caching, not just "doesn't break anything."
"""

from __future__ import annotations

import numpy as np
import pytest

import sougata_solver.simulation as simulation_module
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 0.7
WAVELENGTH = 0.8
AIR = Material("air", 1.0)
SI = Material("si", 12.11)
NUM_ORDERS = 9
NUM_REPEATS = 5


def _lattice() -> Lattice:
    return Lattice((PERIOD, 0.0), (0.0, PERIOD))


def _excitation() -> PlaneWaveExcitation:
    return PlaneWaveExcitation(wavelength=WAVELENGTH, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)


def _pillar_pattern() -> Pattern:
    pattern = Pattern(background=AIR)
    pattern.add(Circle(center=(0.0, 0.0), radius=0.2, material=SI))
    return pattern


def _repeated_pattern_layers(pattern: Pattern) -> list[Layer]:
    return [Layer(f"pillar_{i}", 0.1, pattern=pattern) for i in range(NUM_REPEATS)]


def test_cached_result_matches_forced_uncached_recomputation():
    """Equivalence check: force every cache lookup to miss (by clearing the
    dict inside `_cached_toeplitz` via a wrapper) and compare against the
    normal, actually-cached path."""
    pattern = _pillar_pattern()

    sim_cached = Simulation(_lattice(), _repeated_pattern_layers(pattern), NUM_ORDERS, AIR, AIR)
    result_cached = sim_cached.solve(_excitation())

    sim_forced_uncached = Simulation(_lattice(), _repeated_pattern_layers(pattern), NUM_ORDERS, AIR, AIR)
    real_cached_toeplitz = sim_forced_uncached._cached_toeplitz

    def _always_miss_toeplitz(pattern_arg, g, wavelength, inverse):
        sim_forced_uncached._toeplitz_cache.clear()
        return real_cached_toeplitz(pattern_arg, g, wavelength, inverse)

    sim_forced_uncached._cached_toeplitz = _always_miss_toeplitz
    result_uncached = sim_forced_uncached.solve(_excitation())

    assert result_cached.reflectance() == pytest.approx(result_uncached.reflectance(), abs=1e-12)
    assert result_cached.transmittance() == pytest.approx(result_uncached.transmittance(), abs=1e-12)
    np.testing.assert_allclose(result_cached.a_transmitted, result_uncached.a_transmitted, atol=1e-12)
    np.testing.assert_allclose(result_cached.b_reflected, result_uncached.b_reflected, atol=1e-12)


def test_repeated_identical_pattern_triggers_exactly_one_real_toeplitz_call(monkeypatch):
    pattern = _pillar_pattern()
    sim = Simulation(_lattice(), _repeated_pattern_layers(pattern), NUM_ORDERS, AIR, AIR)

    call_count = 0
    real_toeplitz_matrix = simulation_module.toeplitz_matrix

    def _counting_toeplitz_matrix(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return real_toeplitz_matrix(*args, **kwargs)

    monkeypatch.setattr(simulation_module, "toeplitz_matrix", _counting_toeplitz_matrix)

    sim.solve(_excitation())

    # NUM_REPEATS identical-pattern layers, only `inverse=False` is used on
    # the 2D isotropic path -> exactly 1 real call, not NUM_REPEATS.
    assert call_count == 1, f"expected exactly 1 real toeplitz_matrix call, got {call_count}"


def test_toeplitz_cache_is_reused_across_an_angle_sweep_at_fixed_wavelength():
    """The measured justification for this cache (`design.md`'s "Layer/
    Toeplitz Caching Design", `decisions.md` ADR-016): `toeplitz_matrix`
    depends only on `(pattern, wavelength)`, not on incidence angle, so a
    fixed-wavelength angle sweep (Category 8 target 8.3, planned) should
    populate exactly one cache entry regardless of how many angles are
    solved, not one per `solve()` call."""
    pattern = _pillar_pattern()
    sim = Simulation(_lattice(), [Layer("l0", 0.1, pattern=pattern)], NUM_ORDERS, AIR, AIR)

    thetas = [0.0, 0.1, 0.2, 0.3, 0.4]
    results = []
    for theta in thetas:
        exc = PlaneWaveExcitation(wavelength=WAVELENGTH, theta=theta, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
        results.append(sim.solve(exc))

    assert len(sim._toeplitz_cache) == 1, f"expected 1 cache entry across the angle sweep, got {len(sim._toeplitz_cache)}"
    for result in results:
        assert 0.0 <= result.reflectance() <= 1.0 + 1e-8
        assert 0.0 <= result.transmittance() <= 1.0 + 1e-8


def test_distinct_patterns_are_cached_independently():
    """Two different (non-identical-object) patterns in the same stack must
    not collide in the cache -- confirms the `id(pattern)` key component is
    actually load-bearing, not incidentally unused."""
    pattern_a = _pillar_pattern()
    pattern_b = _pillar_pattern()
    pattern_b.shapes[0].radius = 0.3  # genuinely different physical layer

    layers = [Layer("a", 0.1, pattern=pattern_a), Layer("b", 0.1, pattern=pattern_b)]
    sim = Simulation(_lattice(), layers, NUM_ORDERS, AIR, AIR)
    result = sim.solve(_excitation())

    assert 0.0 <= result.reflectance() <= 1.0 + 1e-8
    assert 0.0 <= result.transmittance() <= 1.0 + 1e-8

    key_a = ("toeplitz", id(pattern_a), WAVELENGTH, False)
    key_b = ("toeplitz", id(pattern_b), WAVELENGTH, False)
    assert key_a in sim._toeplitz_cache
    assert key_b in sim._toeplitz_cache
    assert not np.array_equal(sim._toeplitz_cache[key_a], sim._toeplitz_cache[key_b])
