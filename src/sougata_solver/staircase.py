"""Staircase (z-discretized) layer-stack generator for tapered sidewalls.

Phase 5. Approximates a linearly-tapered via/trench sidewall as
`num_slices` uniform-in-z layers, each a `Circle`/`Rectangle`/`Slab`
pattern whose size is linearly interpolated between a `top` and `bottom`
value -- the standard staircase/multi-slice approximation used throughout
RCWA literature for slanted-sidewall gratings, per `decisions.md` ADR-004
("Tapered sidewalls via staircase discretization, not new Fourier math").

This module is **not** transcribed from any vendored `REFERENCE/` repo: per
the `phase-reference-picker` skill's procedure, a grep across every
RCWA-family repo (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl`,
`Rigorous-Coupled-Wave-Analysis`) for "stair"/"taper" found no matching
staircase-generator code to cite (only unrelated hits in `meep`/`gprMax`
docs, a different numerical method per that skill's own guidance) -- so
this is independently derived, per `rules.md` AI Coding Rule 1, and flagged
for the extra test scrutiny that entails. Its correctness claim rests on
the convergence-vs-`num_slices` study in `tests/test_staircase.py`, not an
external oracle: each individual slice reuses Phase 3/4a's
already-oracle-validated per-layer eigenmode solve unchanged, so the only
new risk here is the discretization itself, not a new physics formula
(matches `phases.md` Phase 5's "no new Fourier/eigenmode math" scoping).

Interpolation convention: slice `i` (0-indexed, `i=0` nearest the `top`
size, `i=num_slices-1` nearest `bottom`) uses the linearly-interpolated
size at its z-midpoint, `frac = (i + 0.5) / num_slices` -- a symmetric
choice with no a-priori reason to bias toward either slice edge. All
`num_slices` slices have equal thickness `thickness / num_slices`.

Category 4 target 4.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`) generalizes the
above into `slice_profile`, a shape-agnostic geometry-to-layer-slices
interface parametrized by an arbitrary `pattern_at(frac)` callable, of
which linear top/bottom taper is one particular choice, not the only one
the module supports. `staircase_circle_layers`/`staircase_rectangle_layers`/
`staircase_slab_layers` are now thin wrappers around `slice_profile`
(refactored, not reimplemented -- `tests/test_staircase.py`'s full existing
suite, including the Phase 5 zero-taper/energy-conservation regression
tests, passes unchanged after the refactor, confirmed rather than assumed).
See `tests/test_profile_slicing.py` for `slice_profile`'s own tests,
independent of any specific shape or of `Simulation.solve`.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from sougata_solver.geometry import Circle, Pattern, Rectangle, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material


def _slice_fractions(num_slices: int) -> np.ndarray:
    if num_slices < 1:
        raise ValueError("num_slices must be >= 1")
    return (np.arange(num_slices) + 0.5) / num_slices


def slice_profile(
    thickness: float,
    num_slices: int,
    pattern_at: Callable[[float], Pattern],
    name_prefix: str = "slice",
) -> list[Layer]:
    """Category 4 target 4.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): general
    geometry-to-layer-slices interface, independent of the RCWA solve --
    this function never imports or calls anything from `simulation.py`/
    `eigenmodes.py`, only `Layer`/`Pattern` construction, so it is testable
    (and tested, `tests/test_profile_slicing.py`) without running a solve.

    Given a callable `pattern_at(frac)` producing the cross-sectional
    `Pattern` at normalized depth `frac in (0, 1)` (same z-midpoint
    convention as the module docstring's `frac = (i + 0.5) / num_slices`),
    generate `num_slices` uniform-in-z `Layer`s of equal thickness
    `thickness / num_slices`. This is the shape-agnostic generalization of
    `staircase_circle_layers`/`staircase_rectangle_layers`/
    `staircase_slab_layers` (each now a thin wrapper around this function,
    see below) -- of which linear top/bottom taper is one particular
    `pattern_at` choice, not the only one; a caller can pass any function
    of `frac` (e.g. a non-linear taper profile, or a `Pattern` with more
    than one shape per slice), all sharing this same slicing/bookkeeping
    logic without needing a new generator function per shape type.
    """
    fracs = _slice_fractions(num_slices)
    slice_thickness = thickness / num_slices
    return [
        Layer(f"{name_prefix}_{i}", slice_thickness, pattern=pattern_at(float(frac)))
        for i, frac in enumerate(fracs)
    ]


def staircase_circle_layers(
    center: tuple[float, float],
    top_radius: float,
    bottom_radius: float,
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    name_prefix: str = "via",
) -> list[Layer]:
    """`num_slices` uniform-in-z `Circle` layers approximating a linearly
    tapered via: radius `top_radius` at z=0 tapering to `bottom_radius` at
    z=thickness. See module docstring for the interpolation convention.

    A thin wrapper around `slice_profile` (target 4.7) -- linear radius
    taper is one particular `pattern_at` choice among the general
    interface's possible callers, not a separate code path."""

    def pattern_at(frac: float) -> Pattern:
        radius = top_radius + (bottom_radius - top_radius) * frac
        return Pattern(background=background_material, shapes=[Circle(center=center, radius=radius, material=shape_material)])

    return slice_profile(thickness, num_slices, pattern_at, name_prefix)


def staircase_rectangle_layers(
    center: tuple[float, float],
    top_halfwidth: tuple[float, float],
    bottom_halfwidth: tuple[float, float],
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    angle: float = 0.0,
    name_prefix: str = "via",
) -> list[Layer]:
    """Same as `staircase_circle_layers`, for a `Rectangle` via whose
    `(hx, hy)` halfwidths taper independently between `top_halfwidth` and
    `bottom_halfwidth`. A thin wrapper around `slice_profile` (target 4.7),
    same as `staircase_circle_layers`."""
    top_hw = np.asarray(top_halfwidth, dtype=float)
    bottom_hw = np.asarray(bottom_halfwidth, dtype=float)

    def pattern_at(frac: float) -> Pattern:
        hw = top_hw + (bottom_hw - top_hw) * frac
        shape = Rectangle(center=center, halfwidth=(float(hw[0]), float(hw[1])), material=shape_material, angle=angle)
        return Pattern(background=background_material, shapes=[shape])

    return slice_profile(thickness, num_slices, pattern_at, name_prefix)


def staircase_slab_layers(
    center_x: float,
    top_halfwidth: float,
    bottom_halfwidth: float,
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    name_prefix: str = "trench",
) -> list[Layer]:
    """Same as `staircase_circle_layers`, for a `Slab` (1D trench) whose
    halfwidth tapers between `top_halfwidth` and `bottom_halfwidth`, center
    fixed at `center_x`. A thin wrapper around `slice_profile` (target 4.7),
    same as `staircase_circle_layers`."""

    def pattern_at(frac: float) -> Pattern:
        halfwidth = top_halfwidth + (bottom_halfwidth - top_halfwidth) * frac
        shape = Slab(center_x=center_x, halfwidth=halfwidth, material=shape_material)
        return Pattern(background=background_material, shapes=[shape])

    return slice_profile(thickness, num_slices, pattern_at, name_prefix)
