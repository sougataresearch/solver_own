"""Category 8 targets 8.6-8.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`sweep.harmonic_study`, `sweep.find_convergence_index`, and
`sweep.auto_select_num_orders`.

Per target 8.8's own wording ("implement only after 8.7 succeeds on
thin-film, trench, and pillar fixtures"), `find_convergence_index` is
validated against three structurally different fixture types before
`auto_select_num_orders` is trusted to use it automatically:

- **thin-film** (Phase 1, uniform layers): `num_orders` has *no* physical
  effect at all (no patterning to diffract), so convergence must be
  detected immediately, at the very first candidate.
- **trench** (Phase 3, 1D grating, TE polarization, `n=3.48`): genuinely
  converges, monotonically, at a moderate `num_orders` -- the "ordinary"
  case.
- **pillar** (Category 3's high-contrast 2D stress fixture,
  `tests/test_fourier_convergence.py`): deliberately non-monotonic at low
  `num_orders` (`num_orders=25` is a wild outlier vs. its neighbors) --
  the exact scenario `find_convergence_index`'s "conservative" design
  (every *later* point must also stay within tolerance, not just the next
  one) exists to not be fooled by.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import auto_select_num_orders, find_convergence_index, harmonic_study

AIR = Material("air", 1.0)


# ---------------------------------------------------------------------------
# 8.6 harmonic_study: equivalence + conservation residual
# ---------------------------------------------------------------------------

TRENCH_PERIOD = 0.7
TRENCH_FILL_FACTOR = 0.3
TRENCH_THICKNESS = 0.46
TRENCH_N_RIDGE = 3.48
TRENCH_SI = Material("si", TRENCH_N_RIDGE**2)


def _build_trench_te_simulation(num_ord: int) -> Simulation:
    lattice = Lattice1D(TRENCH_PERIOD)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=-TRENCH_PERIOD * (1 - TRENCH_FILL_FACTOR) / 2, halfwidth=0.5 * TRENCH_FILL_FACTOR * TRENCH_PERIOD, material=TRENCH_SI))
    layer = Layer("grating", TRENCH_THICKNESS, pattern=pattern)
    return Simulation(lattice, [layer], num_orders=2 * num_ord + 1, incidence=AIR, transmission=AIR)


_TRENCH_EXCITATION = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)


def test_harmonic_study_matches_manual_scalar_loop():
    orders = [2, 4, 6, 8]
    sweep = harmonic_study(_build_trench_te_simulation, orders, _TRENCH_EXCITATION)

    manual_r = [_build_trench_te_simulation(n).solve(_TRENCH_EXCITATION).reflectance() for n in orders]
    np.testing.assert_allclose(sweep.reflectance(), manual_r)
    assert sweep.parameter_name == "num_orders"
    assert list(sweep.parameter_values) == orders


def test_harmonic_study_conservation_residual_is_near_zero_for_lossless_case():
    """This structure has no absorption -- `R+T+sum(layer_absorption())`
    must equal 1 at every harmonic-order count, not just the well-
    converged ones (energy conservation doesn't depend on truncation
    accuracy for a lossless structure)."""
    orders = [2, 4, 6, 8, 10]
    sweep = harmonic_study(_build_trench_te_simulation, orders, _TRENCH_EXCITATION)
    residual = sweep.extra["conservation_residual"]
    assert residual.shape == (5,)
    assert np.all(residual < 1e-8), f"expected near-zero residual for a lossless case, got {residual}"


def test_harmonic_study_rejects_empty_orders():
    with pytest.raises(ValueError, match="non-empty"):
        harmonic_study(_build_trench_te_simulation, [], _TRENCH_EXCITATION)


# ---------------------------------------------------------------------------
# 8.7 find_convergence_index: thin-film / trench / pillar fixtures
# ---------------------------------------------------------------------------


def test_find_convergence_index_thin_film_converges_immediately():
    """A uniform (unpatterned) layer stack: num_orders has no physical
    effect, so reflectance is identical at every candidate -- convergence
    must be found at index 0."""
    lattice = Lattice(a=(1.0, 0.0), b=(0.0, 1.0))
    sio2 = Material("sio2", 1.46**2)
    si = Material("si", 3.48**2)
    values = []
    for num_orders in (1, 9, 25):
        sim = Simulation(lattice, [Layer("sio2", 0.1, material=sio2)], num_orders=num_orders, incidence=AIR, transmission=si)
        result = sim.solve(PlaneWaveExcitation(0.6, 0.0, 0.0, 1.0, 0.0))
        values.append(result.reflectance())

    assert find_convergence_index(values, tolerance=1e-10) == 0


def test_find_convergence_index_trench_te_converges_at_moderate_order():
    """TE (s-polarization) 1D grating, moderate index contrast: converges
    cleanly, not at the first candidate but well before the last one --
    the "ordinary" (not deliberately pathological) convergence case."""
    orders = [2, 4, 6, 8, 10, 15, 20]
    sweep = harmonic_study(_build_trench_te_simulation, orders, _TRENCH_EXCITATION)
    values = sweep.reflectance()

    index = find_convergence_index(values, tolerance=1e-3)
    assert index is not None
    assert 0 < index < len(orders) - 1, f"expected genuine (non-trivial, non-final) convergence, got index={index}"
    # Confirm it's actually correct: every later value really is within tolerance.
    assert np.all(np.abs(values[index:] - values[index]) <= 1e-3)


PILLAR_PERIOD = 0.7
PILLAR_RADIUS_FRAC = 0.2
PILLAR_THICKNESS = 0.3
PILLAR_N = 5.0
PILLAR_SI = Material("si_high_contrast_2d", PILLAR_N**2)


def _build_pillar_simulation(num_orders: int) -> Simulation:
    lattice = Lattice(a=(PILLAR_PERIOD, 0.0), b=(0.0, PILLAR_PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PILLAR_PERIOD / 2, PILLAR_PERIOD / 2), radius=PILLAR_RADIUS_FRAC * PILLAR_PERIOD, material=PILLAR_SI)])
    layer = Layer("pillar", PILLAR_THICKNESS, pattern=pattern)
    return Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR)


_PILLAR_EXCITATION = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)


@pytest.mark.slow
def test_find_convergence_index_pillar_is_not_fooled_by_the_low_order_wobble():
    """Reuses `tests/test_fourier_convergence.py`'s recorded high-contrast
    2D pillar fixture, whose `num_orders=25` point is a known, wild,
    non-monotonic outlier (`R~0.214` vs. neighbors and the eventual
    `~0.0236` converged value). A naive "within tolerance of the next
    point" criterion could be fooled here (e.g. `num_orders=9` and
    `num_orders=49` happen to both be small, but for unrelated reasons --
    one is pre-asymptotic, one is post-outlier); the conservative
    criterion (every *later* point must also match) must not report
    convergence until the sequence has genuinely settled."""
    orders = [9, 25, 49, 81, 121, 169, 225]
    sweep = harmonic_study(_build_pillar_simulation, orders, _PILLAR_EXCITATION)
    values = sweep.reflectance()

    # Confirm the known wobble is actually present in this run (not assumed).
    assert values[1] > 5 * values[0], "expected the known num_orders=25 outlier to be present"

    # A tight-enough tolerance that num_orders=49 (index 2) itself still
    # hasn't settled close enough to the eventual asymptote (~0.0236) --
    # the criterion must not stop there either, only once the sequence has
    # genuinely settled within tolerance for every remaining point.
    index = find_convergence_index(values, tolerance=0.005)
    assert index is not None
    assert index not in (0, 1), f"must not anchor on the num_orders=9/25 wobble, got index={index} (orders={orders})"
    assert index >= 3, f"expected num_orders=49 to still be excluded at this tolerance, got index={index} (orders={orders})"
    assert np.all(np.abs(values[index:] - values[index]) <= 0.005)


def test_find_convergence_index_returns_none_when_never_converges():
    values = [1.0, 2.0, 3.0, 4.0, 5.0]  # monotonically diverging, never settles
    assert find_convergence_index(values, tolerance=1e-6) is None


# ---------------------------------------------------------------------------
# 8.8 auto_select_num_orders
# ---------------------------------------------------------------------------


def test_auto_select_num_orders_matches_find_convergence_index_directly():
    orders = [2, 4, 6, 8, 10, 15, 20]
    selected, sweep = auto_select_num_orders(_build_trench_te_simulation, orders, _TRENCH_EXCITATION, tolerance=1e-3)

    expected_index = find_convergence_index(sweep.reflectance(), tolerance=1e-3)
    assert selected == orders[expected_index]
    assert selected in orders


def test_auto_select_num_orders_raises_when_no_candidate_converges():
    orders = [2, 4]  # far too coarse for the tight tolerance below
    with pytest.raises(ValueError, match="no num_orders"):
        auto_select_num_orders(_build_trench_te_simulation, orders, _TRENCH_EXCITATION, tolerance=1e-12)


def test_auto_select_num_orders_rejects_unsorted_candidates():
    with pytest.raises(ValueError, match="ascending"):
        auto_select_num_orders(_build_trench_te_simulation, [10, 2, 8], _TRENCH_EXCITATION, tolerance=1e-3)


def test_auto_select_num_orders_rejects_invalid_metric():
    with pytest.raises(ValueError, match="metric"):
        auto_select_num_orders(_build_trench_te_simulation, [2, 4], _TRENCH_EXCITATION, tolerance=1e-3, metric="bogus")
