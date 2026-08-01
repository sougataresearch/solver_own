"""Phase 5 (tapered/sloped sidewalls) tests.

Tiers per `testing.md`: unit (slice fractions/thicknesses, degenerate
`N=1` case), regression (zero-taper staircase reduces to Phase 4a/3's
already-validated uniform-cross-section result), physical-invariant
(energy conservation for a tapered case), and Phase 5's own explicit
deliverable -- a convergence-vs-`num_slices` study (`slow`) for both a
tapered via and a tapered trench, per `phases.md` Phase 5 ("no new
Fourier/eigenmode math" means there is no external oracle for this phase;
convergence-vs-`N` is the correctness evidence, as `staircase.py`'s module
docstring explains).
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice, Lattice1D, Pattern, Slab
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import (
    _slice_fractions,
    staircase_circle_layers,
    staircase_rectangle_layers,
    staircase_slab_layers,
)

PERIOD = 0.7
THICKNESS = 0.46
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# Unit: slice fractions / thicknesses
# ---------------------------------------------------------------------------


def test_slice_fractions_midpoints():
    assert _slice_fractions(1) == pytest.approx([0.5])
    assert _slice_fractions(4) == pytest.approx([0.125, 0.375, 0.625, 0.875])


def test_slice_fractions_rejects_nonpositive_count():
    with pytest.raises(ValueError):
        _slice_fractions(0)


def test_staircase_circle_layers_thicknesses_sum_to_total():
    layers = staircase_circle_layers(
        center=(PERIOD / 2, PERIOD / 2),
        top_radius=0.20,
        bottom_radius=0.10,
        thickness=THICKNESS,
        num_slices=8,
        shape_material=SI,
        background_material=AIR,
    )
    assert len(layers) == 8
    assert sum(layer.thickness for layer in layers) == pytest.approx(THICKNESS)


def test_staircase_circle_layers_radii_interpolate_monotonically():
    layers = staircase_circle_layers(
        center=(PERIOD / 2, PERIOD / 2),
        top_radius=0.20,
        bottom_radius=0.10,
        thickness=THICKNESS,
        num_slices=5,
        shape_material=SI,
        background_material=AIR,
    )
    radii = [layer.pattern.shapes[0].radius for layer in layers]
    assert radii == sorted(radii, reverse=True)
    assert radii[0] < 0.20
    assert radii[-1] > 0.10


# ---------------------------------------------------------------------------
# Regression: zero-taper staircase reduces to a single already-validated layer
# ---------------------------------------------------------------------------


def test_staircase_circle_zero_taper_matches_single_uniform_pillar_layer():
    """`top_radius == bottom_radius` must reproduce Phase 4a's single-layer
    pillar result regardless of `num_slices` -- if the staircase generator
    or slice-thickness bookkeeping were wrong, this would drift."""
    radius = 0.18
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    from sougata_solver.geometry import Circle
    from sougata_solver.layer import Layer

    reference_pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=radius, material=SI)])
    reference_layer = Layer("pillar", THICKNESS, pattern=reference_pattern)
    reference_sim = Simulation(lattice, [reference_layer], num_orders=5, incidence=AIR, transmission=AIR)
    reference_result = reference_sim.solve(excitation)

    tapered_layers = staircase_circle_layers(
        center=(PERIOD / 2, PERIOD / 2),
        top_radius=radius,
        bottom_radius=radius,
        thickness=THICKNESS,
        num_slices=6,
        shape_material=SI,
        background_material=AIR,
    )
    tapered_sim = Simulation(lattice, tapered_layers, num_orders=5, incidence=AIR, transmission=AIR)
    tapered_result = tapered_sim.solve(excitation)

    assert tapered_result.reflectance() == pytest.approx(reference_result.reflectance(), abs=1e-10)
    assert tapered_result.transmittance() == pytest.approx(reference_result.transmittance(), abs=1e-10)


def test_staircase_slab_zero_taper_matches_single_uniform_grating_layer():
    from sougata_solver.layer import Layer

    fill_halfwidth = 0.5 * 0.3 * PERIOD
    lattice = Lattice1D(PERIOD)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    reference_pattern = Pattern(background=AIR, shapes=[Slab(center_x=0.0, halfwidth=fill_halfwidth, material=SI)])
    reference_layer = Layer("grating", THICKNESS, pattern=reference_pattern)
    reference_sim = Simulation(lattice, [reference_layer], num_orders=11, incidence=AIR, transmission=AIR)
    reference_result = reference_sim.solve(excitation)

    tapered_layers = staircase_slab_layers(
        center_x=0.0,
        top_halfwidth=fill_halfwidth,
        bottom_halfwidth=fill_halfwidth,
        thickness=THICKNESS,
        num_slices=6,
        shape_material=SI,
        background_material=AIR,
    )
    tapered_sim = Simulation(lattice, tapered_layers, num_orders=11, incidence=AIR, transmission=AIR)
    tapered_result = tapered_sim.solve(excitation)

    assert tapered_result.reflectance() == pytest.approx(reference_result.reflectance(), abs=1e-10)
    assert tapered_result.transmittance() == pytest.approx(reference_result.transmittance(), abs=1e-10)


# ---------------------------------------------------------------------------
# Physical invariant: energy conservation for an actually-tapered case
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 20.0])
def test_tapered_via_energy_conservation(theta_deg):
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    layers = staircase_circle_layers(
        center=(PERIOD / 2, PERIOD / 2),
        top_radius=0.22,
        bottom_radius=0.12,
        thickness=THICKNESS,
        num_slices=6,
        shape_material=SI,
        background_material=AIR,
    )
    sim = Simulation(lattice, layers, num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=np.radians(theta_deg), phi=0.0, s_amplitude=0.6, p_amplitude=0.8)
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)


def test_tapered_trench_energy_conservation():
    lattice = Lattice1D(PERIOD)
    layers = staircase_slab_layers(
        center_x=0.0,
        top_halfwidth=0.5 * 0.4 * PERIOD,
        bottom_halfwidth=0.5 * 0.2 * PERIOD,
        thickness=THICKNESS,
        num_slices=6,
        shape_material=SI,
        background_material=AIR,
    )
    sim = Simulation(lattice, layers, num_orders=11, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Rectangle via generator: sanity + zero-taper reduction (unit + regression)
# ---------------------------------------------------------------------------


def test_staircase_rectangle_zero_taper_matches_single_uniform_layer():
    from sougata_solver.geometry import Rectangle
    from sougata_solver.layer import Layer

    hw = (0.15, 0.10)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    reference_pattern = Pattern(background=AIR, shapes=[Rectangle(center=(PERIOD / 2, PERIOD / 2), halfwidth=hw, material=SI)])
    reference_layer = Layer("pillar", THICKNESS, pattern=reference_pattern)
    reference_sim = Simulation(lattice, [reference_layer], num_orders=5, incidence=AIR, transmission=AIR)
    reference_result = reference_sim.solve(excitation)

    tapered_layers = staircase_rectangle_layers(
        center=(PERIOD / 2, PERIOD / 2),
        top_halfwidth=hw,
        bottom_halfwidth=hw,
        thickness=THICKNESS,
        num_slices=4,
        shape_material=SI,
        background_material=AIR,
    )
    tapered_sim = Simulation(lattice, tapered_layers, num_orders=5, incidence=AIR, transmission=AIR)
    tapered_result = tapered_sim.solve(excitation)

    assert tapered_result.reflectance() == pytest.approx(reference_result.reflectance(), abs=1e-10)
    assert tapered_result.transmittance() == pytest.approx(reference_result.transmittance(), abs=1e-10)


# ---------------------------------------------------------------------------
# Phase 5's own deliverable: convergence vs. num_slices (slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_tapered_via_converges_with_increasing_num_slices():
    """Measures, not assumes, monotone-ish convergence: R at a fixed
    wavelength/angle should settle down as `num_slices` grows, with the
    successive-difference shrinking towards the finest step tried. No
    external oracle exists for this phase (see module docstring) -- this
    convergence trend *is* the correctness evidence."""
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    slice_counts = [1, 2, 4, 8, 16, 32]
    reflectances = []
    for num_slices in slice_counts:
        layers = staircase_circle_layers(
            center=(PERIOD / 2, PERIOD / 2),
            top_radius=0.24,
            bottom_radius=0.10,
            thickness=THICKNESS,
            num_slices=num_slices,
            shape_material=SI,
            background_material=AIR,
        )
        sim = Simulation(lattice, layers, num_orders=5, incidence=AIR, transmission=AIR)
        reflectances.append(sim.solve(excitation).reflectance())

    diffs = [abs(reflectances[i + 1] - reflectances[i]) for i in range(len(reflectances) - 1)]
    print(f"\ntapered via R vs. num_slices={slice_counts}: {reflectances}")
    print(f"successive |dR|: {diffs}")

    assert diffs[-1] < diffs[0]
    assert diffs[-1] < 1e-3


@pytest.mark.slow
def test_tapered_trench_converges_with_increasing_num_slices():
    lattice = Lattice1D(PERIOD)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    slice_counts = [1, 2, 4, 8, 16, 32, 64]
    reflectances = []
    for num_slices in slice_counts:
        layers = staircase_slab_layers(
            center_x=0.0,
            top_halfwidth=0.5 * 0.5 * PERIOD,
            bottom_halfwidth=0.5 * 0.15 * PERIOD,
            thickness=THICKNESS,
            num_slices=num_slices,
            shape_material=SI,
            background_material=AIR,
        )
        sim = Simulation(lattice, layers, num_orders=11, incidence=AIR, transmission=AIR)
        reflectances.append(sim.solve(excitation).reflectance())

    diffs = [abs(reflectances[i + 1] - reflectances[i]) for i in range(len(reflectances) - 1)]
    print(f"\ntapered trench R vs. num_slices={slice_counts}: {reflectances}")
    print(f"successive |dR|: {diffs}")

    assert diffs[-1] < diffs[0]
    assert diffs[-1] < 1e-3
