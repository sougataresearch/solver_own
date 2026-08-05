"""Category 11 target 11.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
layer-to-layer overlay (misregistration) is already achievable via the
existing `Pattern`/`Layer`/`Simulation` API -- no new parameter needed
(`decisions.md` ADR-019, same treatment as ADR-014's bottom-illumination
finding). This file is the permanent regression guard for that claim.
"""

from __future__ import annotations

import math

import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
PERIOD = 0.7


def _via_over_pad_result(dx: float, dy: float):
    """A via layer over a landing-pad layer, the pad shifted by `(dx, dy)`
    relative to the via -- an overlay error vector, both layers sharing
    the same `Lattice`."""
    top_pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.15, material=SI)])
    bottom_pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35 + dx, 0.35 + dy), radius=0.2, material=SI)])
    lattice = Lattice((PERIOD, 0.0), (0.0, PERIOD))
    layers = [Layer("via", 0.2, pattern=top_pattern), Layer("pad", 0.1, pattern=bottom_pattern)]
    sim = Simulation(lattice, layers, num_orders=9, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    return sim.solve(excitation)


def test_overlay_shift_changes_result_and_conserves_energy():
    """A genuine overlay error is a real structural change -- R/T must
    differ from the zero-overlay case, not be silently ignored."""
    zero = _via_over_pad_result(0.0, 0.0)
    shifted = _via_over_pad_result(0.05, -0.03)

    assert zero.reflectance() != pytest.approx(shifted.reflectance(), abs=1e-6)
    assert zero.reflectance() + zero.transmittance() == pytest.approx(1.0, abs=1e-8)
    assert shifted.reflectance() + shifted.transmittance() == pytest.approx(1.0, abs=1e-8)


def test_overlay_shift_by_one_full_period_reproduces_zero_overlay():
    """The periodicity self-consistency check specific to the "periodic
    unit-cell model" claim (ADR-019): shifting the bottom layer's shape by
    exactly one full lattice period must reproduce the zero-overlay result
    exactly, since a period-periodic pattern is unchanged by a
    whole-period translation."""
    zero = _via_over_pad_result(0.0, 0.0)
    full_period = _via_over_pad_result(PERIOD, 0.0)

    assert full_period.reflectance() == pytest.approx(zero.reflectance(), abs=1e-12)
    assert full_period.transmittance() == pytest.approx(zero.transmittance(), abs=1e-12)


def test_overlay_shift_by_one_full_period_in_y_also_reproduces_zero_overlay():
    zero = _via_over_pad_result(0.0, 0.0)
    full_period_y = _via_over_pad_result(0.0, PERIOD)

    assert full_period_y.reflectance() == pytest.approx(zero.reflectance(), abs=1e-12)
    assert full_period_y.transmittance() == pytest.approx(zero.transmittance(), abs=1e-12)
