"""Category 7 target 7.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): equivalent
repeated-layer *representations* of the same physical structure must
produce the same R/T -- this is the regression guard target 7.4's Toeplitz
cache leans on (it would silently mask a real difference as a cache-hit
coincidence if this invariant didn't already hold uncached).

Two representations are tested:
1. A single thick uniform layer vs. that same physical thickness split
   into N equal thinner layers of the same material (relies on
   `smatrix.interface_smatrix`'s `_is_trivial_interface` fast path between
   identical-material sub-layers being physically correct, not just fast).
2. A patterned layer repeated N times via the *same* `Pattern` Python
   object (object-identity reuse, the case target 7.4's cache keys on) vs.
   N *separately constructed* `Pattern` objects with equal field values
   (structurally the same physical layer, different objects) -- both must
   give the same R/T, confirming the physical result never depends on
   object identity, only on physical content.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 0.7
WAVELENGTH = 1.0
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


def _lattice() -> Lattice:
    return Lattice((PERIOD, 0.0), (0.0, PERIOD))


def _excitation() -> PlaneWaveExcitation:
    return PlaneWaveExcitation(wavelength=WAVELENGTH, theta=0.0, phi=0.0, s_amplitude=0.7, p_amplitude=0.3)


@pytest.mark.parametrize("num_slices", [1, 2, 5, 10])
def test_uniform_layer_split_into_n_equal_slices_matches_one_thick_layer(num_slices):
    total_thickness = 0.5
    one_layer = [Layer("core", total_thickness, material=SI)]
    split_layers = [
        Layer(f"core_{i}", total_thickness / num_slices, material=SI) for i in range(num_slices)
    ]

    sim_one = Simulation(_lattice(), one_layer, num_orders=1, incidence=AIR, transmission=AIR)
    sim_split = Simulation(_lattice(), split_layers, num_orders=1, incidence=AIR, transmission=AIR)

    result_one = sim_one.solve(_excitation())
    result_split = sim_split.solve(_excitation())

    assert result_one.reflectance() == pytest.approx(result_split.reflectance(), abs=1e-10)
    assert result_one.transmittance() == pytest.approx(result_split.transmittance(), abs=1e-10)


def _pillar_pattern() -> Pattern:
    pattern = Pattern(background=AIR)
    pattern.add(Circle(center=(0.0, 0.0), radius=0.2, material=SI))
    return pattern


@pytest.mark.parametrize("num_repeats", [1, 2, 4])
def test_repeated_identical_pattern_object_matches_separately_constructed_equal_patterns(num_repeats):
    shared_pattern = _pillar_pattern()
    layers_shared_object = [
        Layer(f"pillar_{i}", 0.2, pattern=shared_pattern) for i in range(num_repeats)
    ]
    layers_separate_objects = [
        Layer(f"pillar_{i}", 0.2, pattern=_pillar_pattern()) for i in range(num_repeats)
    ]
    # Confirm these really are distinct objects with equal content, not the
    # same object twice -- otherwise this test would not exercise anything
    # the "same object reused" case above doesn't already cover.
    for a, b in zip(layers_shared_object, layers_separate_objects):
        assert a.pattern is not b.pattern

    sim_shared = Simulation(_lattice(), layers_shared_object, num_orders=9, incidence=AIR, transmission=AIR)
    sim_separate = Simulation(_lattice(), layers_separate_objects, num_orders=9, incidence=AIR, transmission=AIR)

    result_shared = sim_shared.solve(_excitation())
    result_separate = sim_separate.solve(_excitation())

    assert result_shared.reflectance() == pytest.approx(result_separate.reflectance(), abs=1e-10)
    assert result_shared.transmittance() == pytest.approx(result_separate.transmittance(), abs=1e-10)
    assert 0.0 <= result_shared.reflectance() + result_shared.transmittance() <= 1.0 + 1e-8
