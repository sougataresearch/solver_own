"""Category 7 target 7.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): construction-
time invariants for `Layer`/`LayerStack` -- finite-layer thickness, the
semi-infinite half-space sentinel, and patterned-layer invariants. The
exactly-one-of-`material`/`pattern` invariant already has a dedicated test
in `tests/test_failure_contract.py::test_layer_requires_exactly_one_of_material_or_pattern`
(Category 2 target 2.1's Failure Contract audit); not duplicated here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer, LayerStack
from sougata_solver.materials import Material
from sougata_solver.smatrix import SMatrixStack

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# Finite-layer thickness invariant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("thickness", [0.0, -1.0, -0.5, float("-inf"), float("nan")])
def test_invalid_thickness_raises_value_error(thickness):
    with pytest.raises(ValueError, match="thickness"):
        Layer("bad", thickness, material=AIR)


@pytest.mark.parametrize("thickness", [1e-9, 0.5, 1.0, 1e6])
def test_positive_finite_thickness_is_accepted(thickness):
    layer = Layer("ok", thickness, material=AIR)
    assert layer.thickness == thickness


def test_positive_infinite_thickness_is_accepted_as_half_space_sentinel():
    """`math.inf` is deliberately *not* rejected -- it's the documented
    sentinel `LayerStack.__init__` itself uses for the semi-infinite
    incidence/transmission half-spaces (`layer.py`'s class docstring)."""
    layer = Layer("half_space", math.inf, material=AIR)
    assert layer.thickness == math.inf


# ---------------------------------------------------------------------------
# Half-space invariant: LayerStack always appends exactly two semi-infinite
# half-spaces, and their thickness value is never actually consumed by the
# S-matrix cascade (SMatrixStack only calls propagation_smatrix for interior
# layers, indices [1, len-2]) -- verified directly, not just asserted from
# the docstring.
# ---------------------------------------------------------------------------


def test_layer_stack_appends_incidence_and_transmission_half_spaces():
    stack = LayerStack([Layer("core", 0.5, material=SI)], incidence=AIR, transmission=AIR)
    assert len(stack) == 3
    assert stack[0].name == "incidence" and stack[0].thickness == math.inf
    assert stack[-1].name == "transmission" and stack[-1].thickness == math.inf
    assert stack[1].name == "core"


def test_half_space_thickness_value_is_never_propagated_through():
    """Build the same interior-layer physics twice, with the outer two
    (half-space) `thicknesses` entries set to two different finite values
    fed directly into `SMatrixStack` (bypassing `Layer`'s own validation,
    which is exactly the point -- this isolates whether `SMatrixStack`
    itself ever reads those two entries). If the half-space thickness were
    ever consumed, the two `full_smatrix()` results would differ."""
    omega = 2 * math.pi / 0.8
    modes_inc = solve_layer_eigenmodes_uniform(omega, np.array([0.0]), np.array([0.0]), AIR.epsilon_tensor(0.8)[0, 0])
    modes_core = solve_layer_eigenmodes_uniform(omega, np.array([0.0]), np.array([0.0]), SI.epsilon_tensor(0.8)[0, 0])
    modes_trans = solve_layer_eigenmodes_uniform(omega, np.array([0.0]), np.array([0.0]), AIR.epsilon_tensor(0.8)[0, 0])
    all_modes = [modes_inc, modes_core, modes_trans]

    stack_a = SMatrixStack([123.456, 0.3, 987.654], all_modes)
    stack_b = SMatrixStack([1.0, 0.3, 1e9], all_modes)
    np.testing.assert_array_equal(stack_a.full_smatrix(), stack_b.full_smatrix())


# ---------------------------------------------------------------------------
# Patterned-layer invariant: a patterned Layer's background material is
# reachable via background_material() regardless of how many shapes it
# holds (including zero, the "pattern with no shapes reduces to a uniform
# background" edge case already relied on elsewhere, e.g. Category 4
# target 4.2's reduction tests).
# ---------------------------------------------------------------------------


def test_patterned_layer_background_material_with_no_shapes():
    pattern = Pattern(background=SI, shapes=[])
    layer = Layer("empty_pattern", 0.5, pattern=pattern)
    assert layer.is_uniform() is False
    assert layer.background_material() is SI


def test_patterned_layer_background_material_with_shapes():
    pattern = Pattern(background=AIR)
    pattern.add(Circle(center=(0.0, 0.0), radius=0.2, material=SI))
    layer = Layer("pillar", 0.5, pattern=pattern)
    assert layer.background_material() is AIR


def test_uniform_layer_background_material_is_its_own_material():
    layer = Layer("core", 0.5, material=SI)
    assert layer.is_uniform() is True
    assert layer.background_material() is SI
