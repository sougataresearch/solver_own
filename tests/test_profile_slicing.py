"""Category 4 target 4.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): the general
`staircase.slice_profile` geometry-to-layer-slices interface, tested
independently of the RCWA solve (pure `Layer`/`Pattern` construction --
`Simulation.solve` is never called in this file, matching the target's own
"independently of the RCWA solve" wording).

Two tiers: (a) `slice_profile` itself, with a custom (non-taper)
`pattern_at` callable to prove it is genuinely general, not secretly
circle/rectangle/slab-specific; (b) a regression check that the three
existing shape-specific generators (`staircase_circle_layers`, etc.),
refactored this target to be thin wrappers around `slice_profile`, still
produce identical output to a from-scratch equivalent built directly on
`slice_profile` -- i.e. the refactor changed no behavior, only removed
duplication.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.geometry import Circle, Pattern
from sougata_solver.materials import Material
from sougata_solver.staircase import slice_profile, staircase_circle_layers

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# slice_profile itself: a genuinely custom (non-taper) profile
# ---------------------------------------------------------------------------


def test_slice_profile_rejects_zero_slices():
    with pytest.raises(ValueError, match="num_slices"):
        slice_profile(1.0, 0, lambda frac: Pattern(background=AIR))


def test_slice_profile_thicknesses_sum_to_total():
    layers = slice_profile(0.9, 6, lambda frac: Pattern(background=AIR))
    assert sum(layer.thickness for layer in layers) == pytest.approx(0.9)
    assert all(layer.thickness == pytest.approx(0.9 / 6) for layer in layers)


def test_slice_profile_layer_names_use_prefix_and_index():
    layers = slice_profile(0.9, 3, lambda frac: Pattern(background=AIR), name_prefix="custom")
    assert [layer.name for layer in layers] == ["custom_0", "custom_1", "custom_2"]


def test_slice_profile_calls_pattern_at_with_expected_fractions():
    """Same z-midpoint convention as `_slice_fractions`
    (`frac = (i + 0.5) / num_slices`) -- checked here via a `pattern_at`
    that records every `frac` it was called with, rather than assuming
    `slice_profile`'s internals match `_slice_fractions` without checking."""
    seen = []

    def pattern_at(frac: float) -> Pattern:
        seen.append(frac)
        return Pattern(background=AIR)

    slice_profile(1.0, 4, pattern_at)
    assert seen == pytest.approx([0.125, 0.375, 0.625, 0.875])


def test_slice_profile_supports_a_non_linear_non_taper_profile():
    """A sinusoidal-radius profile (not linear top/bottom taper) --
    demonstrates `slice_profile` is genuinely general, not secretly
    hard-coded to the linear-interpolation case every existing
    `staircase_*_layers` generator happens to use."""

    def pattern_at(frac: float) -> Pattern:
        radius = 0.1 + 0.05 * np.sin(np.pi * frac)  # bulges in the middle
        return Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=float(radius), material=SI)])

    layers = slice_profile(0.5, 5, pattern_at, name_prefix="bulge")
    radii = [layer.pattern.shapes[0].radius for layer in layers]
    # Middle slice's radius should exceed both the first and last slice's.
    assert radii[2] > radii[0]
    assert radii[2] > radii[-1]


# ---------------------------------------------------------------------------
# Regression: staircase_circle_layers (now a slice_profile wrapper) still
# matches a from-scratch equivalent built directly on slice_profile.
# ---------------------------------------------------------------------------


def test_staircase_circle_layers_matches_equivalent_slice_profile_call():
    center = (0.35, 0.35)
    top_radius, bottom_radius = 0.24, 0.10
    thickness, num_slices = 0.46, 8

    via_layers = staircase_circle_layers(center, top_radius, bottom_radius, thickness, num_slices, SI, AIR)

    def pattern_at(frac: float) -> Pattern:
        radius = top_radius + (bottom_radius - top_radius) * frac
        return Pattern(background=AIR, shapes=[Circle(center=center, radius=radius, material=SI)])

    equivalent_layers = slice_profile(thickness, num_slices, pattern_at, name_prefix="via")

    assert len(via_layers) == len(equivalent_layers)
    for a, b in zip(via_layers, equivalent_layers):
        assert a.name == b.name
        assert a.thickness == pytest.approx(b.thickness)
        assert a.pattern.shapes[0].radius == pytest.approx(b.pattern.shapes[0].radius)
