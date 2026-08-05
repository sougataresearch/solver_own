"""Category 6 target 6.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): diffraction-
order opening/closing (Rayleigh threshold) at **oblique** incidence.

`tests/test_mode_classification.py` (Category 1 target 1.8) already covers
this target's literal wording ("add a diffraction-order opening/closing
case and test propagating versus evanescent classification") at *normal*
incidence, against an exact analytic threshold wavelength. That target is
therefore already met by existing coverage -- see
`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`'s own status entry for the cross-
reference. This file adds the one dimension that coverage doesn't exercise
and Category 6 is specifically about (excitation-angle dependence): at
**oblique** incidence, the `+m`/`-m` order degeneracy normal incidence has
is broken (`kx0 != 0` shifts each order's effective in-plane wavevector
differently), so the two orders cross their Rayleigh threshold at
*different* wavelengths -- confirmed directly by a coarse wavelength sweep
this session (not derived from a closed form, since the oblique threshold
condition is not as simple as `lambda = n*period/m`): the `(1,0)` order
transitions between 1200-1250 nm, `(-1,0)` between 1750-1800 nm, at
`theta=15 deg`. Energy conservation holds throughout.
"""

from __future__ import annotations

import math

import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 1e-6
THETA_DEG = 15.0


def _solve_at(wavelength_nm: float):
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((PERIOD, 0.0), (0.0, PERIOD))
    layer = Layer("filler", 0.1e-6, material=Material("filler", 2.0))
    sim = Simulation(lattice, [layer], num_orders=9, incidence=air, transmission=glass, truncation="circular")
    theta = math.radians(THETA_DEG)
    excitation = PlaneWaveExcitation(wavelength_nm * 1e-9, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    return sim.solve(excitation)


def test_positive_and_negative_orders_have_different_oblique_thresholds():
    """At theta=0 the (1,0) and (-1,0) orders share one threshold
    (`test_mode_classification.py`); at oblique incidence they must not --
    confirms the classification genuinely depends on the excitation angle,
    not just the lattice/wavelength (this is Category 6's subject, not
    Category 1's)."""
    below_both = _solve_at(1200).order_classification()
    between = _solve_at(1500).order_classification()
    above_both = _solve_at(1800).order_classification()

    assert below_both[(1, 0)]["transmission"] == "propagating"
    assert below_both[(-1, 0)]["transmission"] == "propagating"

    # Strictly between the two thresholds: (1,0) has already closed,
    # (-1,0) has not yet -- the asymmetry oblique incidence predicts.
    assert between[(1, 0)]["transmission"] == "evanescent"
    assert between[(-1, 0)]["transmission"] == "propagating"

    assert above_both[(1, 0)]["transmission"] == "evanescent"
    assert above_both[(-1, 0)]["transmission"] == "evanescent"


@pytest.mark.parametrize("wavelength_nm", [1200, 1250, 1300, 1500, 1700, 1750, 1800])
def test_energy_conservation_across_oblique_order_opening_and_closing(wavelength_nm):
    result = _solve_at(wavelength_nm)
    de = result.diffraction_efficiencies()
    total = sum(r + t for r, t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)
