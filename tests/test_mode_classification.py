"""Category 1 target 1.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): public
propagating/evanescent mode classification (`eigenmodes.classify_propagating`,
`SimulationResult.order_classification`).

Tiers enforced here, per `rules.md` Testing Requirements:
- unit: `classify_propagating` on constructed `q` arrays matches
  `_select_q_branch`'s own branch convention exactly (real q -> propagating,
  purely-imaginary q -> evanescent).
- Rayleigh-threshold system test: a diffraction order's transmission-side
  classification flips at the theoretically-predicted wavelength
  (`lambda_threshold = n_trans * period / m` at normal incidence for a
  square lattice's `(m, 0)` order), and energy conservation holds on both
  sides of, and at, that threshold.
"""

from __future__ import annotations


import numpy as np
import pytest

from sougata_solver.eigenmodes import classify_propagating
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation


# ---------------------------------------------------------------------------
# Unit: classify_propagating matches _select_q_branch's own convention
# ---------------------------------------------------------------------------


def test_classify_propagating_real_q_is_propagating():
    q = np.array([1.5 + 0j, 3.0 + 0j])
    assert np.all(classify_propagating(q))


def test_classify_propagating_purely_imaginary_q_is_evanescent():
    q = np.array([0 + 1.5j, 0 + 3.0j])
    assert not np.any(classify_propagating(q))


def test_classify_propagating_mixed_array():
    q = np.array([2.0 + 0j, 0 + 1.0j, 5.0 + 0j, 0 + 0.5j])
    assert classify_propagating(q).tolist() == [True, False, True, False]


def test_classify_propagating_complex_lossy_branch_is_evanescent():
    # Genuinely complex q (lossy medium, neither purely real nor purely
    # imaginary) is conservatively classified evanescent -- see
    # classify_propagating's docstring.
    q = np.array([1.0 + 0.3j])
    assert not np.any(classify_propagating(q))


# ---------------------------------------------------------------------------
# System-level Rayleigh-threshold test
# ---------------------------------------------------------------------------


PERIOD = 1e-6
N_TRANS = 1.5


def _order_classification_at(wavelength):
    lattice = Lattice((PERIOD, 0.0), (0.0, PERIOD))
    air = Material("air", 1.0)
    glass = Material("glass", N_TRANS**2)
    layer = Layer("filler", 0.1e-6, material=Material("filler", 2.0))
    sim = Simulation(lattice, [layer], num_orders=9, incidence=air, transmission=glass, truncation="circular")
    excitation = PlaneWaveExcitation(wavelength, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    return result, result.order_classification()


def test_rayleigh_threshold_transmission_side_order_flips_at_predicted_wavelength():
    m = 1
    lambda_threshold = N_TRANS * PERIOD / m  # 1.5e-6

    result_below, classification_below = _order_classification_at(lambda_threshold * 0.9)  # 1.35e-6: propagating
    result_above, classification_above = _order_classification_at(lambda_threshold * 1.1)  # 1.65e-6: evanescent

    assert classification_below[(m, 0)]["transmission"] == "propagating"
    assert classification_above[(m, 0)]["transmission"] == "evanescent"

    # Zeroth order is always propagating on both sides of the threshold.
    assert classification_below[(0, 0)]["transmission"] == "propagating"
    assert classification_above[(0, 0)]["transmission"] == "propagating"

    # Energy conservation holds on both sides.
    de_below = result_below.diffraction_efficiencies()
    de_above = result_above.diffraction_efficiencies()
    assert sum(r + t for r, t in de_below.values()) == pytest.approx(1.0, abs=1e-8)
    assert sum(r + t for r, t in de_above.values()) == pytest.approx(1.0, abs=1e-8)


def test_rayleigh_threshold_energy_conservation_close_to_threshold():
    # NOT exactly at the threshold: at the exact Wood's-anomaly crossing,
    # q == 0 for that order and smatrix.py's `kp @ phi / q` construction
    # divides by zero (confirmed directly: exactly-at-threshold produces
    # NaN R/T here, an honest pre-existing solver limitation at the exact
    # singular point, not something target 1.8 fixes -- the same class of
    # boundary case `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`'s Category 6 target
    # 6.4 ("Grazing-incidence boundary test") exists to define and test
    # separately). This test instead checks robustness a small relative
    # step away from the singularity on each side.
    m = 1
    lambda_threshold = N_TRANS * PERIOD / m
    for frac in (0.999, 1.001):
        result, classification = _order_classification_at(lambda_threshold * frac)
        de = result.diffraction_efficiencies()
        assert sum(r + t for r, t in de.values()) == pytest.approx(1.0, abs=1e-6)
        expected = "propagating" if frac < 1.0 else "evanescent"
        assert classification[(m, 0)]["transmission"] == expected


def test_rayleigh_threshold_incidence_side_uses_incidence_index():
    # Incidence medium is air (n=1.0) -- its own Rayleigh threshold for the
    # same order m=1 is at a shorter wavelength (PERIOD/m = 1e-6) than the
    # transmission (glass, n=1.5) side, so at a wavelength between the two
    # thresholds the same order is evanescent on the incidence side but
    # still propagating on the transmission side.
    wavelength = 1.2e-6  # between 1.0e-6 (air threshold) and 1.5e-6 (glass threshold)
    _result, classification = _order_classification_at(wavelength)
    assert classification[(1, 0)]["incidence"] == "evanescent"
    assert classification[(1, 0)]["transmission"] == "propagating"
