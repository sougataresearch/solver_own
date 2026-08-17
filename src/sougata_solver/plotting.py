"""Category 16 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): plotting functions.

## Plot data contract (target 16.1)

Every function here takes **plain arrays, dataclasses, or already-computed
result objects** (`geometry.Pattern`/`Lattice`, `layer.LayerStack`,
NumPy arrays, or a `SimulationResult`/`SweepResult` a caller already
solved) -- never a bare `Simulation`, and no function anywhere in this
module calls `.solve()`. This mirrors `decisions.md` ADR-009/010's
`structures/`-vs-`postprocessing/` split (plotting is strictly downstream
of a solve, never triggers one) at the library-function level instead of
only the script level: a caller who already has a `SimulationResult` can
plot it without this module ever touching the solver, and a caller who
only has raw arrays (e.g. loaded back from a `.npz`/CSV a `structures/`
script wrote) can plot those directly with no `Simulation` object at all.

Every function returns `(fig, ax)` (matplotlib `Figure`/`Axes`) rather
than saving or showing anything itself -- saving/showing is the caller's
job (`postprocessing/*.py` scripts already do this, e.g.
`plot_thin_film_rt.py`'s `fig.savefig(...)`/`plt.show()`), keeping this
module a pure "data in, figure out" library with no filesystem or display
side effects, consistent with `rules.md`'s "don't add scope beyond what
was asked" and this project's existing raw-data-only export convention
(`fields.save_field_grid_npz`, Category 9 target 9.8).

`matplotlib` is imported lazily inside each function (not at module level)
so importing `sougata_solver.plotting` doesn't force a `matplotlib`
dependency onto every user of the library -- the same lazy-import pattern
`postprocessing/plot_field_cross_section.py` already uses.
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.geometry import Pattern


def rasterize_pattern(pattern: Pattern, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Rasterize `pattern` on the `xs x ys` sample grid, one integer per
    sample point: `0` = background, `i+1` = `pattern.shapes[i]`.

    Extracted from `plot_unit_cell` (Category 16 target 16.2, the original
    caller) so `plot_structure_3d` can reuse the identical "iterate shapes
    in reverse, first `.contains(x, y)` match wins" logic that already
    implements `Pattern`'s documented "later shapes take precedence" rule
    (`geometry.Pattern`'s docstring) -- one rasterization implementation,
    not two that could silently drift apart.
    """
    grid = np.zeros((len(ys), len(xs)), dtype=int)  # 0 = background
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            for shape_index in range(len(pattern.shapes) - 1, -1, -1):
                if pattern.shapes[shape_index].contains(x, y):
                    grid[iy, ix] = shape_index + 1
                    break
    return grid


def plot_unit_cell(pattern: Pattern, lattice, *, resolution: int = 200, ax=None):
    """Category 16 target 16.2: render one periodic unit cell of `pattern`
    on `lattice` -- a rasterized preview (grid of `resolution x resolution`
    sample points, colored by whichever shape's `.contains(x, y)` matches,
    background otherwise), using `Pattern`'s own already-documented
    "later shapes take precedence" rule (`geometry.Pattern`'s docstring) so
    the preview matches the solver's actual precedence semantics rather
    than a naive first-shape-wins raster. This is a **preview raster for
    visualization only** -- it plays no role in the solve itself, which
    uses `pattern`'s analytic Fourier transforms unmodified (Category 4).

    Takes a `Pattern`/`Lattice` pair directly (not a `Simulation`), per
    target 16.1's data contract -- callable from a `structures/*.py`
    script before ever constructing a `Simulation`.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    a, b = lattice.a, lattice.b
    corners = np.array([[0, 0], a, a + b, b, [0, 0]])
    x_min, x_max = corners[:, 0].min(), corners[:, 0].max()
    y_min, y_max = corners[:, 1].min(), corners[:, 1].max()
    xs = np.linspace(x_min, x_max, resolution)
    ys = np.linspace(y_min, y_max, resolution)
    grid = rasterize_pattern(pattern, xs, ys)

    n_shapes = len(pattern.shapes)
    colors = ["#dddddd"] + [plt.cm.tab10(i % 10) for i in range(n_shapes)]
    cmap = ListedColormap(colors)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure
    ax.pcolormesh(xs, ys, grid, cmap=cmap, vmin=0, vmax=max(n_shapes, 1), shading="nearest")
    ax.plot(corners[:, 0], corners[:, 1], color="black", linewidth=1.5, linestyle="--")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Unit cell ({n_shapes} shape{'s' if n_shapes != 1 else ''})")
    return fig, ax


def plot_layer_stack(thicknesses, labels, *, ax=None):
    """Category 16 target 16.2: render the layer stack as a 1D column of
    slabs (semi-infinite incidence/transmission drawn as a fixed-height
    hatched band since `math.inf` has no natural plot height), one bar per
    finite layer sized proportionally to its thickness.

    Takes plain `thicknesses`/`labels` sequences, not a `Simulation` or
    `LayerStack` -- so a caller can plot a stack design before it's ever
    wired into a `Simulation` (or plot one reconstructed from a
    `config.py` file without importing the solver at all).
    """
    import matplotlib.pyplot as plt

    if len(thicknesses) != len(labels):
        raise ValueError(f"thicknesses ({len(thicknesses)}) and labels ({len(labels)}) must have the same length")

    finite = [t for t in thicknesses if np.isfinite(t)]
    semi_infinite_height = max(finite) * 0.5 if finite else 1.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 6))
    else:
        fig = ax.figure

    y = 0.0
    for thickness, label in zip(thicknesses, labels):
        height = thickness if np.isfinite(thickness) else semi_infinite_height
        hatch = "//" if not np.isfinite(thickness) else None
        ax.bar(0, height, bottom=y, width=1.0, edgecolor="black", hatch=hatch, label=label)
        ax.text(0.55, y + height / 2, label, va="center", ha="left", fontsize=9)
        y += height

    ax.set_xlim(-0.6, 2.2)
    ax.set_ylim(0, y)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_ylabel("depth (same units as thickness input)")
    ax.set_title("Layer stack")
    return fig, ax


def plot_structure_3d(layer_stack, lattice, *, resolution: int = 40, extrusion_length: float | None = None, ax=None):
    """New target (`decisions.md` ADR-029): render a full `Layer` stack as
    a 3D solid, one voxel slab per layer stacked at its real cumulative
    z-offset -- a static structure preview, the 3D analogue of
    `plot_unit_cell` (per-layer cross-section) x `plot_layer_stack`
    (z-ordering), so a via/pillar/trench (including a tapered one built as
    a `staircase.py` `list[Layer]`) can be visually inspected before or
    instead of running a solve.

    **Not a physics formula** -- this function computes no new solver
    quantity, only a visualization of geometry the solver already
    consumes unmodified (same class of exemption `plot_unit_cell`/
    `plot_layer_stack` already rely on, per `rules.md` Documentation
    Standards). It is independently derived (not transcribed from any
    vendored oracle), and validated by structural tests in
    `tests/test_plotting.py` (voxel z-extent/slice-count/material-count
    checks), not pixel comparison -- this project has no golden-image
    infrastructure (see `tests/test_plotting.py`'s own module docstring).

    Reuses `rasterize_pattern` (this module, factored out of
    `plot_unit_cell`) for each layer's in-plane cross-section, so a
    patterned layer's rendering uses the exact same "later shape wins"
    precedence rasterization already validated for `plot_unit_cell` --
    no new per-`Shape` geometry code. A `Slab` (1D `Lattice1D` shape,
    infinite in y) is handled automatically by its own `.contains()`
    semantics for whatever y-range is rasterized here, no special case
    needed.

    Takes `layer_stack` (a `layer.LayerStack` or a plain `list[Layer]`,
    e.g. `staircase.staircase_circle_layers(...)`'s direct return value)
    and `lattice` (`geometry.Lattice`/`Lattice1D`) directly -- never a
    bare `Simulation`, and this function never calls `.solve()`, per
    target 16.1's data contract (`plot_unit_cell`'s docstring) that every
    function in this module already follows.

    `math.inf`-thickness incidence/transmission half-spaces are **not
    rendered at all** -- reverted here after two rounds of a fabricated
    end-cap thickness reading as visually wrong (first a shared
    `max()`-based height rendered a staircased substrate *thinner* than
    the trench it should contain; fixing that to a shared `sum()`-based
    height next made the incidence half-space look as visually dominant
    as the substrate). There is no size for a genuinely semi-infinite
    half-space that reads correctly next to a finite patterned stack, so
    per the project owner's explicit direction this function only ever
    renders the **finite** layers -- the patterned/staircased region
    itself, which for a via/trench already reads correctly on its own
    (background material fills the full cross-section except the
    via/trench footprint, exactly "air etched into a solid substrate"),
    with no separate substrate block competing for visual attention. A
    `layer_stack` with no finite layer raises `ValueError`.

    `lattice.b == (0, 0)` for a `Lattice1D` grating (`geometry.Lattice1D`)
    means there is no natural y-extent to rasterize -- `plot_unit_cell`
    has this same latent gap (its bounding box collapses to zero height
    for a 1D lattice), not something new to this function.
    `extrusion_length` (default: `lattice.a`'s magnitude, i.e. one period)
    makes the y-depth an explicit, documented visualization choice for
    the 1D case rather than an implicit zero-height result.

    **Measured performance note** (not assumed, per `rules.md`
    Performance Requirements): matplotlib's `Axes3D.voxels` is slow and
    scales worse than linearly in `resolution**2 * len(layer_stack)` (its
    hidden-face computation checks every voxel's neighbors) --
    `resolution=20` with an 8-slice staircase renders in ~3s,
    `resolution=40` with a 16-slice staircase takes ~2 minutes on this
    project's dev machine. Keep `resolution` in the 15-25 range for a
    many-slice staircase preview; raise it only for a final, one-off look.
    """
    import matplotlib.pyplot as plt

    layers = list(layer_stack)
    if not layers:
        raise ValueError("layer_stack must contain at least one layer")
    layers = [ly for ly in layers if math.isfinite(ly.thickness)]
    if not layers:
        raise ValueError("layer_stack must contain at least one finite-thickness layer to render")

    a = np.asarray(lattice.a, dtype=float)
    b = np.asarray(lattice.b, dtype=float)
    is_1d = not np.any(b)
    if is_1d:
        length = extrusion_length if extrusion_length is not None else float(np.linalg.norm(a))
        x_min, x_max = -0.5 * float(np.linalg.norm(a)), 0.5 * float(np.linalg.norm(a))
        y_min, y_max = -0.5 * length, 0.5 * length
    else:
        corners = np.array([[0, 0], a, a + b, b, [0, 0]])
        x_min, x_max = corners[:, 0].min(), corners[:, 0].max()
        y_min, y_max = corners[:, 1].min(), corners[:, 1].max()
    x_edges = np.linspace(x_min, x_max, resolution + 1)
    y_edges = np.linspace(y_min, y_max, resolution + 1)
    xs = 0.5 * (x_edges[:-1] + x_edges[1:])
    ys = 0.5 * (y_edges[:-1] + y_edges[1:])

    material_names: list[str] = []
    for ly in layers:
        names = [ly.material.name] if ly.is_uniform() else [ly.pattern.background.name] + [
            shape.material.name for shape in ly.pattern.shapes
        ]
        for name in names:
            if name not in material_names:
                material_names.append(name)
    # Color by material *name*, sorted case-insensitively -- not by
    # encounter order (found and fixed this session, third finding):
    # encounter order made "air" land on whichever color slot its
    # background/shape position happened to be in, so e.g. a Si pillar-
    # in-air and an air-via-in-Si rendered with identical-looking colors
    # (air landed on the same slot both times purely by coincidence of
    # which one was `Pattern.background` first), even though physically
    # one structure is mostly air and the other mostly solid. Sorting by
    # name means the same material always gets the same color across
    # every call, so which material is majority/minority is visible from
    # color alone, not just position.
    color_map = {name: plt.cm.tab10(i % 10) for i, name in enumerate(sorted(material_names, key=str.lower))}

    if ax is None:
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(projection="3d")
    else:
        fig = ax.figure

    filled = np.ones((resolution, resolution, 1), dtype=bool)
    z = 0.0
    for ly in layers:
        height = ly.thickness
        if ly.is_uniform():
            grid = np.zeros((resolution, resolution), dtype=int)
            names_by_index = {0: ly.material.name}
        else:
            grid = rasterize_pattern(ly.pattern, xs, ys)
            names_by_index = {0: ly.pattern.background.name}
            names_by_index.update({i + 1: shape.material.name for i, shape in enumerate(ly.pattern.shapes)})
        facecolors = np.empty((resolution, resolution, 1, 4))
        for iy in range(resolution):
            for ix in range(resolution):
                facecolors[iy, ix, 0] = color_map[names_by_index[grid[iy, ix]]]
        z_edges = np.array([z, z + height])
        xg, yg, zg = np.meshgrid(x_edges, y_edges, z_edges, indexing="ij")
        ax.voxels(xg, yg, zg, filled, facecolors=facecolors.swapaxes(0, 1), edgecolor=None)
        z += height

    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=color_map[name]) for name in material_names]
    ax.legend(handles, material_names, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z (depth)")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(0.0, z)
    ax.invert_zaxis()
    ax.set_title("Structure preview")
    return fig, ax


def plot_rt_spectrum(wavelengths, reflectance, transmittance=None, *, metadata=None, ax=None):
    """Category 16 target 16.3: formalizes
    `postprocessing/plot_thin_film_rt.py`'s ad hoc R-vs-wavelength plot
    into a reusable function -- same axis labels/style, generalized to
    optionally also plot `transmittance` and an `R+T` conservation trace,
    and to attach `metadata` (a plain dict, e.g. `{"num_orders": 9,
    "materials": "air/SiO2/Si"}`) as a labeled annotation rather than
    leaving the reader to infer run parameters from the filename alone.
    """
    import matplotlib.pyplot as plt

    wavelengths = np.asarray(wavelengths, dtype=float)
    reflectance = np.asarray(reflectance, dtype=float)
    if len(wavelengths) != len(reflectance):
        raise ValueError(f"wavelengths ({len(wavelengths)}) and reflectance ({len(reflectance)}) must have the same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    ax.plot(wavelengths, reflectance, color="tab:blue", label="R")
    if transmittance is not None:
        transmittance = np.asarray(transmittance, dtype=float)
        if len(transmittance) != len(wavelengths):
            raise ValueError(f"transmittance ({len(transmittance)}) must match wavelengths ({len(wavelengths)})")
        ax.plot(wavelengths, transmittance, color="tab:orange", label="T")
        ax.plot(wavelengths, reflectance + transmittance, color="tab:gray", linestyle=":", label="R+T")
    ax.set_xlabel("Wavelength")
    ax.set_ylabel("Efficiency")
    ax.set_ylim(0, max(1.05, float(np.max(reflectance)) * 1.05))
    ax.set_title("R/T vs wavelength")
    ax.legend()
    ax.grid(True, alpha=0.3)
    if metadata:
        text = "\n".join(f"{k}: {v}" for k, v in metadata.items())
        ax.text(0.02, 0.98, text, transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    return fig, ax


def plot_harmonic_convergence(num_orders_values, values, *, convergence_index=None, value_label="Reflectance", ax=None):
    """Category 16 target 16.4: plot a Category 8
    `sweep.harmonic_study`/`find_convergence_index` result -- `value`
    (e.g. reflectance) vs. `num_orders_values`, with the converged point
    (if `convergence_index` is given, exactly the index
    `find_convergence_index` already returns) marked distinctly so the
    plot visually matches what the already-validated convergence
    criterion actually selected, not a human's eyeball guess.
    """
    import matplotlib.pyplot as plt

    num_orders_values = np.asarray(num_orders_values)
    values = np.asarray(values, dtype=float)
    if len(num_orders_values) != len(values):
        raise ValueError(f"num_orders_values ({len(num_orders_values)}) and values ({len(values)}) must have the same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    ax.plot(num_orders_values, values, "o-", color="tab:blue", label=value_label)
    if convergence_index is not None:
        ax.plot(num_orders_values[convergence_index], values[convergence_index], "o", color="tab:red",
                markersize=12, markerfacecolor="none", markeredgewidth=2, label="converged (find_convergence_index)")
    ax.set_xlabel("num_orders")
    ax.set_ylabel(value_label)
    ax.set_title("Harmonic-order convergence")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig, ax


def plot_diffraction_orders(diffraction_efficiencies, *, kind="reflected", ax=None):
    """Category 16 target 16.5: bar plot of per-order diffraction
    efficiency from `SimulationResult.diffraction_efficiencies()`'s
    already-validated `{(g1, g2): (R_order, T_order)}` dict -- `kind`
    selects `"reflected"` (index 0) or `"transmitted"` (index 1) of each
    tuple. Orders are sorted by `(g1, g2)` for a deterministic bar order
    across repeated calls (dict insertion order is otherwise incidental).
    """
    import matplotlib.pyplot as plt

    if kind not in ("reflected", "transmitted"):
        raise ValueError(f"kind must be 'reflected' or 'transmitted', got {kind!r}")
    value_index = 0 if kind == "reflected" else 1

    orders = sorted(diffraction_efficiencies.keys())
    values = [diffraction_efficiencies[order][value_index] for order in orders]
    labels = [f"({g1},{g2})" for g1, g2 in orders]

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(orders)), 5))
    else:
        fig = ax.figure
    ax.bar(range(len(orders)), values, color="tab:blue")
    ax.set_xticks(range(len(orders)))
    ax.set_xticklabels(labels, rotation=90 if len(orders) > 12 else 0)
    ax.set_xlabel("diffraction order (g1, g2)")
    ax.set_ylabel(f"{kind} efficiency")
    ax.set_title(f"Per-order {kind} efficiency")
    ax.grid(True, axis="y", alpha=0.3)
    return fig, ax


def plot_field_intensity(horizontal, vertical, ex, ey, ez, *, horizontal_label="x", vertical_label="y", ax=None):
    """Category 16 target 16.6: formalizes
    `postprocessing/plot_field_cross_section.py`'s `pcolormesh` intensity
    plot (`|E|^2 = |Ex|^2+|Ey|^2+|Ez|^2`, `fields.modal_field_components`'s
    already-validated Cartesian convention, Category 9) into a reusable
    function taking the raw field-grid arrays directly -- the same
    `.npz` layout `fields.save_field_grid_npz` (target 9.8) and the
    `structures/*_field_cross_section.py` scripts already produce.
    """
    import matplotlib.pyplot as plt

    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    mesh = ax.pcolormesh(horizontal, vertical, intensity, shading="auto", cmap="inferno")
    fig.colorbar(mesh, ax=ax, label="|E|^2")
    ax.set_xlabel(horizontal_label)
    ax.set_ylabel(vertical_label)
    ax.set_title("Field intensity")
    return fig, ax


def plot_field_phase(horizontal, vertical, component, *, component_label="Ex", horizontal_label="x", vertical_label="y", ax=None):
    """Category 16 target 16.7: phase map of one complex field component
    (`np.angle`, radians in `[-pi, pi]`) -- eligible once field-component
    conventions are validated (Category 9, `CONVENTIONS.md`), which this
    project's real-space reconstruction (`fields.reconstruct_field_at_points`)
    already is. A cyclic colormap (`twilight`) is used deliberately, not a
    sequential one, since phase wraps at `+-pi` and a sequential map would
    show a spurious discontinuity there that isn't physically real.
    """
    import matplotlib.pyplot as plt

    phase = np.angle(component)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    mesh = ax.pcolormesh(horizontal, vertical, phase, shading="auto", cmap="twilight", vmin=-np.pi, vmax=np.pi)
    fig.colorbar(mesh, ax=ax, label=f"phase({component_label}) [rad]")
    ax.set_xlabel(horizontal_label)
    ax.set_ylabel(vertical_label)
    ax.set_title(f"{component_label} phase")
    return fig, ax


def plot_poynting_vector(horizontal, vertical, sx, sz, *, horizontal_label="x", vertical_label="z", ax=None, subsample=8):
    """Category 16 target 16.7: in-plane Poynting-vector direction field
    (`quiver`) over a real-space cross-section, `sx`/`sz` already computed
    by the caller from reconstructed fields (`Sz = Re(Ex*conj(Hy) -
    Ey*conj(Hx))`, per `CONVENTIONS.md`'s and `troubleshooting.md`'s
    documented no-`0.5`-factor real-space flux convention, Category 9
    target 9.6's finding -- this function does no flux computation of its
    own, only visualizes already-computed values, per target 16.1's data
    contract). `subsample` thins the quiver grid for readability on a
    fine field grid.
    """
    import matplotlib.pyplot as plt

    horizontal, vertical = np.asarray(horizontal), np.asarray(vertical)
    sx, sz = np.asarray(sx, dtype=float), np.asarray(sz, dtype=float)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    sl = (slice(None, None, subsample), slice(None, None, subsample))
    ax.quiver(horizontal[sl], vertical[sl], sx[sl], sz[sl], color="tab:red")
    ax.set_xlabel(horizontal_label)
    ax.set_ylabel(vertical_label)
    ax.set_title("Poynting vector direction")
    ax.set_aspect("equal")
    return fig, ax
