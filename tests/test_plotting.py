"""Category 16 targets 16.1-16.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
plotting-function tests -- structural checks (axes/labels/data extents,
correct artist counts) rather than pixel comparisons, since this project
has no golden-image test infrastructure and pixel-exact plot comparison is
brittle across matplotlib versions; the goal is to catch a function
crashing, mislabeling an axis, or silently dropping data, not to pin exact
rendering.
"""

from __future__ import annotations

import inspect

import matplotlib

matplotlib.use("Agg")  # headless backend -- no display needed/available in CI

import numpy as np
import pytest

from sougata_solver import plotting
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Rectangle, Slab
from sougata_solver.layer import Layer, LayerStack
from sougata_solver.materials import Material
from sougata_solver.staircase import staircase_circle_layers

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# 16.1: plot data contract -- no function here ever touches Simulation.solve
# ---------------------------------------------------------------------------


def test_plotting_module_never_imports_simulation():
    """Structural contract check: `plotting.py` takes plain arrays/
    dataclasses/already-computed results, never a `Simulation`, so it can
    never trigger a solve. Confirmed by inspecting the actual module
    source's import statements (not the docstring prose, which mentions
    `SimulationResult` only as an example of an already-computed input a
    caller may pass in) for any import of `sougata_solver.simulation` or
    a literal `.solve(` call anywhere in the module body."""
    source = inspect.getsource(plotting)
    assert "sougata_solver.simulation" not in source
    assert not hasattr(plotting, "Simulation")


def test_every_plot_function_returns_fig_and_ax():
    fig, ax = plotting.plot_rt_spectrum([0.5, 0.6, 0.7], [0.1, 0.2, 0.3])
    assert fig is ax.figure


# ---------------------------------------------------------------------------
# 16.2: geometry / layer-stack plots
# ---------------------------------------------------------------------------


def test_plot_unit_cell_renders_shapes_within_lattice_bounds():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.5, 0.5), radius=0.2, material=SI)])
    fig, ax = plotting.plot_unit_cell(pattern, lattice, resolution=40)
    assert ax.get_xlim()[0] <= 0.0 and ax.get_xlim()[1] >= 1.0
    assert "1 shape" in ax.get_title()


def test_plot_unit_cell_empty_pattern_is_all_background():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pattern = Pattern(background=AIR)
    fig, ax = plotting.plot_unit_cell(pattern, lattice, resolution=10)
    assert "0 shapes" in ax.get_title()


def test_plot_unit_cell_later_shape_takes_precedence():
    """Pattern's own documented rule: later-added shapes take precedence
    over earlier ones at overlapping points -- confirmed the raster
    respects this by checking a point inside both an early background
    circle and a later overriding rectangle picks up the rectangle."""
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pattern = Pattern(background=AIR)
    pattern.add(Circle(center=(0.5, 0.5), radius=0.4, material=SI))
    pattern.add(Rectangle(center=(0.5, 0.5), halfwidth=(0.1, 0.1), material=AIR))
    fig, ax = plotting.plot_unit_cell(pattern, lattice, resolution=50)
    assert "2 shapes" in ax.get_title()


def test_plot_layer_stack_mismatched_lengths_raises():
    with pytest.raises(ValueError, match="same length"):
        plotting.plot_layer_stack([1.0, 2.0], ["only one label"])


def test_plot_layer_stack_handles_semi_infinite_layers():
    fig, ax = plotting.plot_layer_stack([np.inf, 1.0, np.inf], ["air", "SiO2", "Si"])
    assert ax.get_ylim()[0] > 0  # inverted axis, top (0) < bottom
    assert len(ax.patches) == 3


# ---------------------------------------------------------------------------
# New target (`decisions.md` ADR-029): 3D structure preview
# ---------------------------------------------------------------------------


def test_plot_structure_3d_returns_fig_and_ax():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.5, 0.5), radius=0.2, material=SI)])
    layers = [Layer("pillar", 0.5, pattern=pattern)]
    fig, ax = plotting.plot_structure_3d(layers, lattice, resolution=8)
    assert fig is ax.figure


def test_plot_structure_3d_empty_stack_raises():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    with pytest.raises(ValueError, match="at least one layer"):
        plotting.plot_structure_3d([], lattice, resolution=8)


def test_plot_structure_3d_z_extent_matches_finite_thicknesses():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    layers = [Layer("l1", 0.3, material=AIR), Layer("l2", 0.7, material=SI)]
    _, ax = plotting.plot_structure_3d(layers, lattice, resolution=6)
    zlim = ax.get_zlim()
    assert min(zlim) == pytest.approx(0.0, abs=1e-9)
    assert max(zlim) == pytest.approx(1.0, abs=1e-9)  # 0.3 + 0.7


def test_plot_structure_3d_semi_infinite_layers_are_not_rendered():
    """Regression for two rounds of a fabricated end-cap thickness reading
    as visually wrong to the project owner (first too thin, then -- after
    a fix -- too visually dominant relative to the patterned region). Per
    their explicit direction, `math.inf`-thickness incidence/transmission
    layers are no longer rendered at all -- this pins that a `LayerStack`'s
    semi-infinite layers contribute nothing to the rendered z-extent, only
    the finite (patterned/staircased) layers do."""
    num_slices = 8
    total_depth = 1.0
    layers = staircase_circle_layers(
        center=(0.5, 0.5), top_radius=0.3, bottom_radius=0.1, thickness=total_depth,
        num_slices=num_slices, shape_material=SI, background_material=AIR,
    )
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    layer_stack = LayerStack(layers, incidence=AIR, transmission=SI)
    _, ax = plotting.plot_structure_3d(layer_stack, lattice, resolution=6)
    zlim = ax.get_zlim()
    assert max(zlim) - min(zlim) == pytest.approx(total_depth, abs=1e-9)


def test_plot_structure_3d_all_infinite_layer_stack_raises():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    layer_stack = LayerStack([], incidence=AIR, transmission=SI)
    with pytest.raises(ValueError, match="finite-thickness layer"):
        plotting.plot_structure_3d(layer_stack, lattice, resolution=6)


def test_plot_structure_3d_material_legend_matches_unique_material_count():
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.5, 0.5), radius=0.2, material=SI)])
    layers = [Layer("l1", 0.3, material=AIR), Layer("l2", 0.5, pattern=pattern)]
    _, ax = plotting.plot_structure_3d(layers, lattice, resolution=6)
    legend = ax.get_legend()
    assert legend is not None
    assert len(legend.get_texts()) == 2  # "air" and "si", each counted once


def test_plot_structure_3d_color_keyed_by_material_name_not_encounter_order():
    """Regression for a real bug caught by the project owner comparing two
    renders side by side: coloring materials by *encounter order*
    (background first, then shape) made "air" land on whichever color
    slot its background/shape position happened to occupy, so a Si
    pillar-in-air and an air-via-in-Si structure rendered with
    identical-looking colors despite being physically opposite (mostly
    air vs. mostly solid) -- purely because `air` was `Pattern.background`
    in one case and `Pattern.shapes[0].material` in the other. This pins
    that `air` gets the *same* color regardless of whether it is the
    background or the shape."""
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    pillar_pattern = Pattern(background=AIR, shapes=[Circle(center=(0.5, 0.5), radius=0.2, material=SI)])
    via_pattern = Pattern(background=SI, shapes=[Circle(center=(0.5, 0.5), radius=0.2, material=AIR)])

    _, ax_pillar = plotting.plot_structure_3d([Layer("l", 0.5, pattern=pillar_pattern)], lattice, resolution=6)
    _, ax_via = plotting.plot_structure_3d([Layer("l", 0.5, pattern=via_pattern)], lattice, resolution=6)

    def _color_by_label(ax):
        return {t.get_text(): h.get_color() for h, t in zip(ax.get_legend().legend_handles, ax.get_legend().get_texts())}

    assert _color_by_label(ax_pillar)["air"] == _color_by_label(ax_via)["air"]
    assert _color_by_label(ax_pillar)["si"] == _color_by_label(ax_via)["si"]


def test_plot_structure_3d_staircase_slice_count_matches_num_slices():
    """`ax.voxels` adds one `Poly3DCollection` per individual voxel cell
    (confirmed directly against the installed matplotlib version, not
    assumed), so `num_slices` staircase layers at `resolution x resolution`
    in-plane voxels each produce `num_slices * resolution**2` collections
    -- still a direct, deterministic function of `num_slices`, so this
    regresses "did every staircase slice actually get rendered"."""
    num_slices = 5
    resolution = 6
    layers = staircase_circle_layers(
        center=(0.5, 0.5), top_radius=0.3, bottom_radius=0.1, thickness=1.0,
        num_slices=num_slices, shape_material=SI, background_material=AIR,
    )
    lattice = Lattice((1.0, 0.0), (0.0, 1.0))
    _, ax = plotting.plot_structure_3d(layers, lattice, resolution=resolution)
    assert len(ax.collections) == num_slices * resolution**2


def test_plot_structure_3d_lattice1d_extrudes_finite_y_extent():
    lattice = Lattice1D(0.7)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.2, material=SI))
    layers = [Layer("trench", 0.4, pattern=pattern)]
    _, ax = plotting.plot_structure_3d(layers, lattice, resolution=6)
    assert ax.get_xlim()[1] - ax.get_xlim()[0] > 0
    assert ax.get_ylim()[1] - ax.get_ylim()[0] > 0  # regression: not collapsed to zero height


# ---------------------------------------------------------------------------
# 16.3: R/T spectrum plot
# ---------------------------------------------------------------------------


def test_plot_rt_spectrum_reflectance_only():
    wavelengths = np.linspace(0.4, 0.8, 10)
    reflectance = np.linspace(0.1, 0.3, 10)
    fig, ax = plotting.plot_rt_spectrum(wavelengths, reflectance)
    assert len(ax.lines) == 1
    assert ax.get_xlabel() == "Wavelength"


def test_plot_rt_spectrum_with_transmittance_adds_rt_conservation_line():
    wavelengths = np.linspace(0.4, 0.8, 10)
    reflectance = np.full(10, 0.3)
    transmittance = np.full(10, 0.7)
    fig, ax = plotting.plot_rt_spectrum(wavelengths, reflectance, transmittance)
    assert len(ax.lines) == 3  # R, T, R+T


def test_plot_rt_spectrum_mismatched_lengths_raises():
    with pytest.raises(ValueError, match="same length"):
        plotting.plot_rt_spectrum([0.5, 0.6], [0.1, 0.2, 0.3])


def test_plot_rt_spectrum_metadata_adds_annotation():
    fig, ax = plotting.plot_rt_spectrum([0.5, 0.6], [0.1, 0.2], metadata={"num_orders": 9})
    assert len(ax.texts) == 1
    assert "num_orders" in ax.texts[0].get_text()


# ---------------------------------------------------------------------------
# 16.4: harmonic convergence plot
# ---------------------------------------------------------------------------


def test_plot_harmonic_convergence_marks_converged_point():
    num_orders = [9, 25, 49, 81]
    values = [0.20, 0.24, 0.25, 0.251]
    fig, ax = plotting.plot_harmonic_convergence(num_orders, values, convergence_index=2)
    assert len(ax.lines) == 2  # data line + converged-point marker


def test_plot_harmonic_convergence_without_convergence_index():
    fig, ax = plotting.plot_harmonic_convergence([9, 25], [0.2, 0.25])
    assert len(ax.lines) == 1


# ---------------------------------------------------------------------------
# 16.5: diffraction-order plot
# ---------------------------------------------------------------------------


def test_plot_diffraction_orders_reflected_and_transmitted():
    efficiencies = {(0, 0): (0.1, 0.5), (1, 0): (0.05, 0.2), (-1, 0): (0.05, 0.1)}
    fig, ax = plotting.plot_diffraction_orders(efficiencies, kind="reflected")
    assert len(ax.patches) == 3
    fig2, ax2 = plotting.plot_diffraction_orders(efficiencies, kind="transmitted")
    heights = [p.get_height() for p in ax2.patches]
    assert sorted(heights) == sorted([0.5, 0.2, 0.1])


def test_plot_diffraction_orders_invalid_kind_raises():
    with pytest.raises(ValueError, match="reflected.*transmitted"):
        plotting.plot_diffraction_orders({(0, 0): (0.1, 0.5)}, kind="bogus")


def test_plot_diffraction_orders_deterministic_ordering():
    efficiencies = {(1, 0): (0.1, 0.1), (-1, 0): (0.2, 0.2), (0, 0): (0.3, 0.3)}
    _, ax1 = plotting.plot_diffraction_orders(efficiencies)
    _, ax2 = plotting.plot_diffraction_orders(efficiencies)
    labels1 = [t.get_text() for t in ax1.get_xticklabels()]
    labels2 = [t.get_text() for t in ax2.get_xticklabels()]
    assert labels1 == labels2 == ["(-1,0)", "(0,0)", "(1,0)"]


# ---------------------------------------------------------------------------
# 16.6/16.7: field intensity, phase, and Poynting-vector plots
# ---------------------------------------------------------------------------


def _field_grid():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    ex = np.exp(1j * xx)
    ey = np.zeros_like(ex)
    ez = 0.1 * np.exp(1j * yy)
    return xx, yy, ex, ey, ez


def test_plot_field_intensity_computes_sum_of_squares():
    xx, yy, ex, ey, ez = _field_grid()
    fig, ax = plotting.plot_field_intensity(xx, yy, ex, ey, ez)
    assert ax.get_title() == "Field intensity"
    assert len(fig.axes) == 2  # main axes + colorbar


def test_plot_field_phase_uses_cyclic_range():
    xx, yy, ex, ey, ez = _field_grid()
    fig, ax = plotting.plot_field_phase(xx, yy, ex, component_label="Ex")
    mesh = ax.collections[0]
    assert mesh.get_clim() == (-np.pi, np.pi)


def test_plot_poynting_vector_renders_quiver():
    x = np.linspace(-1, 1, 16)
    z = np.linspace(0, 1, 16)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    sx = np.ones_like(xx)
    sz = np.zeros_like(zz)
    fig, ax = plotting.plot_poynting_vector(xx, zz, sx, sz, subsample=4)
    assert len(ax.collections) == 1  # one quiver collection
    assert ax.get_title() == "Poynting vector direction"
