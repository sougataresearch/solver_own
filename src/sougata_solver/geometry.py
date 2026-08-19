"""Lattice and in-plane pattern geometry.

Reciprocal-vector convention: `kx`, `ky` passed to `Shape.fourier_transform`
are in cycles per unit length (not angular frequency), i.e. a real-space
phase factor enters as ``exp(i * 2*pi * (kx*x + ky*y))``. This matches S4's
convention (`S4/S4/pattern/pattern.c`), where `Lk` (reciprocal lattice) and
`Lr` (real lattice) satisfy `Lr @ Lk.T == I` with no extra `2*pi` factor
folded into the lattice matrices themselves.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
from scipy.special import j1

from sougata_solver.materials import Material


def jinc(x: np.ndarray) -> np.ndarray:
    """`2*J1(2*pi*x) / (2*pi*x)`, with `jinc(0) = 1`.

    This is the radial analogue of `sinc` for 2D Fourier transforms of
    circularly symmetric indicator functions. Source: `pattern.c:951-953`.
    """
    x = np.asarray(x, dtype=float)
    arg = 2.0 * np.pi * x
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(np.abs(x) < 1e-14, 1.0, 2.0 * j1(arg) / arg)
    return result


class Lattice:
    """2D periodic lattice defined by real-space basis vectors `a`, `b`."""

    def __init__(self, a: tuple[float, float], b: tuple[float, float]):
        _require_finite("Lattice basis vector a", *a)
        _require_finite("Lattice basis vector b", *b)
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self._Lr = np.array([self.a, self.b])  # rows are basis vectors
        area = abs(self.a[0] * self.b[1] - self.a[1] * self.b[0])
        if not (area > 0):
            raise ValueError(
                f"Lattice basis vectors a={a}, b={b} are degenerate (zero unit-cell "
                "area) -- they must not be collinear/zero"
            )

    def reciprocal_vectors(self) -> np.ndarray:
        """Return the 2x2 reciprocal basis `Lk` such that `Lr @ Lk.T == I`
        (no `2*pi` factor; see module docstring for the convention)."""
        return np.linalg.inv(self._Lr).T

    def unit_cell_area(self) -> float:
        return abs(self.a[0] * self.b[1] - self.a[1] * self.b[0])


class Lattice1D:
    """1D periodic lattice (grating vector along `x`, invariant along `y`).

    Exposes the same `reciprocal_vectors()`/`unit_cell_area()` interface as
    `Lattice` (duck-typed) so `fourier_factorization.py`'s
    `pattern_epsilon_hat`/`toeplitz_matrix` work unmodified. Matches S4's
    1D-lattice convention (`S4/S4_internal.h:188-190`,
    `S4/S4.cpp:463-481`): the second real-space lattice vector is treated
    as zero, and so is its reciprocal counterpart, rather than computed.
    """

    def __init__(self, period: float):
        _require_finite("Lattice1D period", period)
        _require_positive("Lattice1D period", period)
        self.period = float(period)
        self.a = np.array([self.period, 0.0])
        self.b = np.array([0.0, 0.0])

    def reciprocal_vectors(self) -> np.ndarray:
        return np.array([[1.0 / self.period, 0.0], [0.0, 0.0]])

    def unit_cell_area(self) -> float:
        return self.period


def _rotate(x: float, y: float, angle: float) -> tuple[float, float]:
    c, s = np.cos(angle), np.sin(angle)
    return c * x + s * y, -s * x + c * y


def _require_finite(label: str, *values: float) -> None:
    """Category 4 target 4.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): shared
    validation helper -- fail loud, fail early, at construction (per
    `design.md`'s Error Handling conventions, already followed by
    `Layer.__post_init__`), not deep inside a solve call where a NaN/inf
    dimension would otherwise surface as a cryptic downstream `LinAlgError`
    or a silently wrong (NaN) Fourier coefficient.
    """
    for v in values:
        if not math.isfinite(v):
            raise ValueError(f"{label} must be finite, got {v!r}")


def _require_positive(label: str, *values: float) -> None:
    for v in values:
        if not (v > 0):
            raise ValueError(f"{label} must be > 0, got {v!r}")


class Shape(ABC):
    """One patterned region within a layer, tagged with a material."""

    center: tuple[float, float]
    material: Material

    @abstractmethod
    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        """2D Fourier transform of the indicator function at (kx, ky)
        (cycles per unit length), not yet normalized by unit cell area."""

    @abstractmethod
    def contains(self, x: float, y: float) -> bool:
        ...

    @abstractmethod
    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        """Outward unit normal at the boundary point nearest (x, y)."""

    @property
    @abstractmethod
    def area(self) -> float:
        ...

    @property
    @abstractmethod
    def bounding_radius(self) -> float:
        """Category 4 target 4.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        maximum distance from `center` to any point of the shape --
        used by `validate_pattern_fits_lattice` to conservatively detect
        whether a shape could overlap its own periodic images."""


@dataclass
class Circle(Shape):
    center: tuple[float, float]
    radius: float
    material: Material

    def __post_init__(self) -> None:
        _require_finite("Circle center", *self.center)
        _require_finite("Circle radius", self.radius)
        _require_positive("Circle radius", self.radius)

    @property
    def area(self) -> float:
        return np.pi * self.radius**2

    @property
    def bounding_radius(self) -> float:
        return self.radius

    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        k = np.hypot(kx, ky)
        phase = np.exp(-2j * np.pi * (kx * self.center[0] + ky * self.center[1]))
        return self.area * jinc(self.radius * k) * phase

    def contains(self, x: float, y: float) -> bool:
        dx, dy = x - self.center[0], y - self.center[1]
        return dx * dx + dy * dy <= self.radius**2

    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        dx, dy = x - self.center[0], y - self.center[1]
        r = np.hypot(dx, dy)
        if r < 1e-14:
            return np.array([1.0, 0.0])
        return np.array([dx / r, dy / r])


@dataclass
class Rectangle(Shape):
    center: tuple[float, float]
    halfwidth: tuple[float, float]
    material: Material
    angle: float = 0.0

    def __post_init__(self) -> None:
        _require_finite("Rectangle center", *self.center)
        _require_finite("Rectangle halfwidth", *self.halfwidth)
        _require_positive("Rectangle halfwidth", *self.halfwidth)
        _require_finite("Rectangle angle", self.angle)

    @property
    def area(self) -> float:
        return 4.0 * self.halfwidth[0] * self.halfwidth[1]

    @property
    def bounding_radius(self) -> float:
        return float(np.hypot(*self.halfwidth))

    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        # rotate k into the rectangle's local (unrotated) frame
        klx, kly = _rotate(kx, ky, self.angle)
        hx, hy = self.halfwidth
        phase = np.exp(-2j * np.pi * (kx * self.center[0] + ky * self.center[1]))
        return self.area * np.sinc(2 * klx * hx) * np.sinc(2 * kly * hy) * phase

    def contains(self, x: float, y: float) -> bool:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        hx, hy = self.halfwidth
        return abs(lx) <= hx and abs(ly) <= hy

    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        hx, hy = self.halfwidth
        # nearest edge in local frame: whichever axis is closer to its bound
        if abs(hx - abs(lx)) <= abs(hy - abs(ly)):
            n_local = np.array([np.sign(lx) or 1.0, 0.0])
        else:
            n_local = np.array([0.0, np.sign(ly) or 1.0])
        # rotate normal back to lab frame (inverse of _rotate)
        c, s = np.cos(self.angle), np.sin(self.angle)
        return np.array([c * n_local[0] - s * n_local[1], s * n_local[0] + c * n_local[1]])


@dataclass
class Ellipse(Shape):
    """Category 4 target 4.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): an
    axis-aligned-in-its-own-frame ellipse with semi-axes `halfwidth =
    (hx, hy)`, optionally rotated by `angle` (same CCW-radians, rotate-then-
    translate convention already used by `Rectangle`). Reduces to `Circle`
    when `hx == hy`.

    `fourier_transform` is transcribed directly from
    `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`'s `ELLIPSE`
    case (lines 955-964, read alongside `CIRCLE`'s adjacent case at
    951-954, which `geometry.jinc`'s own docstring already cites) -- an
    ellipse's Fourier transform is `Circle`'s `jinc` formula with the
    `k`-vector anisotropically rescaled by the semi-axis ratio along
    whichever axis is *not* the longer one, so the jinc argument reduces
    exactly to `hx * |k|` (or `hy * |k|`) when `hx == hy` (reduces to
    `Circle`'s formula; regression test in `tests/test_ellipse.py`). Not
    independently re-derived, per `rules.md`'s Documentation Standards
    preference for transcription over re-derivation when a source is
    directly available and the branch choice (`hx >= hy` vs. not) is easy
    to get backwards silently.

    `contains`/`signed_distance_normal` are **not** transcribed from S4's
    `shape_contains_point`/`shape_get_normal` (a numerically-equivalent but
    differently-organized focal-distance/rotation-matrix formulation,
    `pattern.c:163-179,213-236`) -- both are independently derived directly
    from the standard ellipse membership test `(x/hx)^2 + (y/hy)^2 <= 1`
    and its gradient `(x/hx^2, y/hy^2)`, elementary analytic geometry with
    no sign/normalization convention subtle enough to warrant transcription
    risk (per `rules.md` Documentation Standards option 2), validated by
    `tests/test_ellipse.py`'s reduces-to-`Circle`-when-`hx==hy` tests for
    both methods.
    """

    center: tuple[float, float]
    halfwidth: tuple[float, float]
    material: Material
    angle: float = 0.0

    def __post_init__(self) -> None:
        _require_finite("Ellipse center", *self.center)
        _require_finite("Ellipse halfwidth", *self.halfwidth)
        _require_positive("Ellipse halfwidth", *self.halfwidth)
        _require_finite("Ellipse angle", self.angle)

    @property
    def area(self) -> float:
        return np.pi * self.halfwidth[0] * self.halfwidth[1]

    @property
    def bounding_radius(self) -> float:
        return float(max(self.halfwidth))

    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        klx, kly = _rotate(kx, ky, self.angle)
        hx, hy = self.halfwidth
        if hx >= hy:
            r = (hy / hx) * kly
            z = jinc(hx * np.hypot(klx, r))
        else:
            r = (hx / hy) * klx
            z = jinc(hy * np.hypot(r, kly))
        phase = np.exp(-2j * np.pi * (kx * self.center[0] + ky * self.center[1]))
        return self.area * z * phase

    def contains(self, x: float, y: float) -> bool:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        hx, hy = self.halfwidth
        return (lx / hx) ** 2 + (ly / hy) ** 2 <= 1.0

    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        hx, hy = self.halfwidth
        n_local = np.array([lx / hx**2, ly / hy**2])
        norm = np.hypot(*n_local)
        n_local = n_local / norm if norm > 1e-14 else np.array([1.0, 0.0])
        c, s = np.cos(self.angle), np.sin(self.angle)
        return np.array([c * n_local[0] - s * n_local[1], s * n_local[0] + c * n_local[1]])


@dataclass
class Polygon(Shape):
    """Category 4 targets 4.4/4.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): a
    simple (non-self-intersecting), closed polygon with `vertices` given
    CCW in the shape's own local (center-relative, pre-rotation) frame --
    same `center`+`angle` convention as `Rectangle`/`Ellipse`. See
    `decisions.md` ADR-013 for the accuracy-contract decision (target 4.4)
    made before this implementation: analytic (exact, no discretization
    error for a simple polygon), not raster/FFT; self-intersection is not
    checked, matching `S4/S4/pattern/pattern.h`'s own disclosed limitation.

    `fourier_transform` is transcribed directly from
    `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`'s `POLYGON`
    case (lines 974-1008): for edges `(v_p -> v_q)` with `u = v_q - v_p`
    and edge midpoint `rc = (v_q+v_p)/2`,

        S(k) = -i/(2*pi*|k|^2) * sum_edges (u_x*k_y - u_y*k_x) * sinc(k.u) * exp(-2*pi*i*k.rc)

    (`sinc(x) = sin(pi*x)/(pi*x)`, matching `Rectangle`'s existing
    `np.sinc` convention -- S4's own `Sinc` is the same normalized-sinc
    function, confirmed by its `my_j0(M_PI*x)` implementation), with
    `S(0) = area` at the DC term (S4's separate `if(DC)` branch,
    `pattern.c:982-984`). Validated against a from-scratch Riemann-sum
    reference for both a triangle and a non-convex (L-shaped) polygon in
    `tests/test_polygon.py`, not just reduction to a simpler already-tested
    shape (no existing shape reduces to an arbitrary polygon).

    `contains` transcribes S4's PNPoly-based point-in-polygon test
    (`pattern.c:180-193`, itself citing
    http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html,
    a standard public-domain algorithm). `signed_distance_normal` is
    **independently derived** (nearest-edge point-to-segment distance, CCW
    outward normal `(edge_y, -edge_x)`), not transcribed from S4's own
    `shape_get_normal` `POLYGON` case (`pattern.c:256-281`) -- that
    function selects the *farthest* segment (`if(dist > maxdist)`), which
    contradicts this project's own `Shape.signed_distance_normal` contract
    ("nearest boundary point") already established by `Circle`/`Rectangle`/
    `Ellipse`; picking the nearest edge instead is the elementary,
    unambiguous correct answer for that documented contract (per `rules.md`
    Documentation Standards option 2), not a re-derivation of anything
    subtle enough to need S4's convention specifically.
    """

    center: tuple[float, float]
    vertices: tuple[tuple[float, float], ...]
    material: Material
    angle: float = 0.0

    def __post_init__(self) -> None:
        _require_finite("Polygon center", *self.center)
        _require_finite("Polygon angle", self.angle)
        if len(self.vertices) < 3:
            raise ValueError(f"Polygon requires at least 3 vertices, got {len(self.vertices)}")
        for v in self.vertices:
            _require_finite("Polygon vertex", *v)
        if not (self.area > 0):
            raise ValueError(
                f"Polygon vertices must be CCW-ordered with positive (shoelace) area, "
                f"got area={self.area!r} -- check vertex winding order"
            )

    @property
    def area(self) -> float:
        v = np.asarray(self.vertices, dtype=float)
        x, y = v[:, 0], v[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    @property
    def bounding_radius(self) -> float:
        v = np.asarray(self.vertices, dtype=float)
        return float(np.max(np.hypot(v[:, 0], v[:, 1])))

    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        klx, kly = _rotate(kx, ky, self.angle)
        v = np.asarray(self.vertices, dtype=float)
        n = len(v)

        z_re = np.zeros_like(klx, dtype=float)
        z_im = np.zeros_like(klx, dtype=float)
        for p in range(n):
            q = (p + 1) % n
            ux, uy = v[q, 0] - v[p, 0], v[q, 1] - v[p, 1]
            rcx, rcy = 0.5 * (v[q, 0] + v[p, 0]), 0.5 * (v[q, 1] + v[p, 1])
            num = (ux * kly - uy * klx) * np.sinc(klx * ux + kly * uy)
            pa = -2.0 * np.pi * (klx * rcx + kly * rcy)
            z_re = z_re + num * np.sin(pa)
            z_im = z_im - num * np.cos(pa)

        k_sq = klx**2 + kly**2
        is_dc = k_sq == 0.0
        safe_k_sq = np.where(is_dc, 1.0, k_sq)
        z = (z_re + 1j * z_im) / (2.0 * np.pi * safe_k_sq)
        z = np.where(is_dc, complex(self.area), z)

        phase = np.exp(-2j * np.pi * (kx * self.center[0] + ky * self.center[1]))
        return z * phase

    def contains(self, x: float, y: float) -> bool:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        v = self.vertices
        n = len(v)
        inside = False
        vjx, vjy = v[-1]
        for i in range(n):
            vix, viy = v[i]
            if (viy > ly) != (vjy > ly) and lx < (vjx - vix) * (ly - viy) / (vjy - viy) + vix:
                inside = not inside
            vjx, vjy = vix, viy
        return inside

    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        lx, ly = _rotate(x - self.center[0], y - self.center[1], self.angle)
        v = np.asarray(self.vertices, dtype=float)
        n = len(v)
        p_point = np.array([lx, ly])

        best_dist = np.inf
        best_normal = np.array([1.0, 0.0])
        for i in range(n):
            a, b = v[i], v[(i + 1) % n]
            edge = b - a
            edge_len_sq = float(edge @ edge)
            t = np.clip(float((p_point - a) @ edge) / edge_len_sq, 0.0, 1.0)
            closest = a + t * edge
            dist = float(np.hypot(*(p_point - closest)))
            if dist < best_dist:
                normal = np.array([edge[1], -edge[0]])
                norm = np.hypot(*normal)
                if norm > 1e-14:
                    best_dist = dist
                    best_normal = normal / norm

        c, s = np.cos(self.angle), np.sin(self.angle)
        return np.array([c * best_normal[0] - s * best_normal[1], s * best_normal[0] + c * best_normal[1]])


@dataclass
class Slab(Shape):
    """1D analogue of `Rectangle`: an interval `[center-halfwidth,
    center+halfwidth]` along `x`, infinite along `y` (invariant in the
    direction the `Lattice1D` grating is uniform along).

    `area` is a length (`2*halfwidth`), not an area, consistent with
    `Lattice1D.unit_cell_area()` also being a length — both feed the same
    `pattern_epsilon_hat` normalization (`fourier_factorization.py`)
    unmodified. `fourier_transform` ignores `ky` (always `0` for a 1D
    lattice's G-vectors) and reduces to `Rectangle`'s `x`-only sinc factor,
    per `S4/S4/pattern/pattern.c`'s 1D shape handling (same subtraction-rule
    convention already used by `Rectangle`/`Circle`).
    """

    center_x: float
    halfwidth: float
    material: Material

    def __post_init__(self) -> None:
        _require_finite("Slab center_x", self.center_x)
        _require_finite("Slab halfwidth", self.halfwidth)
        _require_positive("Slab halfwidth", self.halfwidth)

    @property
    def center(self) -> tuple[float, float]:
        """`(x, 0.0)`, matching the `Shape` ABC's tuple contract so
        `Pattern.containment_tree`'s `cx, cy = shape.center` unpacking
        works unmodified for a 1D shape."""
        return (self.center_x, 0.0)

    @property
    def area(self) -> float:
        return 2.0 * self.halfwidth

    @property
    def bounding_radius(self) -> float:
        return self.halfwidth

    def fourier_transform(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        kx = np.asarray(kx, dtype=float)
        phase = np.exp(-2j * np.pi * kx * self.center_x)
        return self.area * np.sinc(2 * kx * self.halfwidth) * phase

    def contains(self, x: float, y: float) -> bool:
        return abs(x - self.center_x) <= self.halfwidth

    def signed_distance_normal(self, x: float, y: float) -> np.ndarray:
        return np.array([np.sign(x - self.center_x) or 1.0, 0.0])


@dataclass
class Pattern:
    """Ordered list of shapes within one layer; later shapes take precedence
    over earlier ones at overlapping points (matches S4's `parent[]`
    subtraction-rule convention, `pattern.c:938`)."""

    background: Material
    shapes: list[Shape] = field(default_factory=list)
    skip_bounds_check: bool = False
    """`decisions.md` ADR-035: `validate_pattern_fits_lattice`'s
    `2*bounding_radius >= min_period` test is a conservative, cheap
    *sufficient* condition (see that function's own docstring), not an
    exact overlap test -- for an elongated, non-circular shape whose
    footprint touches (but never crosses) the cell edge along one axis, the
    circular bounding-radius bound can flag a false positive even though no
    true self-overlap occurs. Set this only after independently confirming
    (e.g. a direct real-space containment check against the shape's
    periodic images) that the shape genuinely does not overlap itself --
    this is a narrow, per-pattern escape hatch, not a general relaxation of
    the check, which stays on by default for every other pattern."""

    def add(self, shape: Shape) -> None:
        self.shapes.append(shape)

    def containment_tree(self) -> list[int | None]:
        """For each shape, return the index of the shape it is nested
        inside (the smallest-area shape added before it whose interior
        contains this shape's center), or None if it sits directly on the
        background. Used to apply the Fourier-coefficient subtraction rule
        for nested/composite shapes."""
        parents: list[int | None] = []
        for i, shape in enumerate(self.shapes):
            cx, cy = shape.center
            best_parent = None
            best_area = np.inf
            for j in range(i):
                other = self.shapes[j]
                if other.contains(cx, cy) and other.area < best_area:
                    best_parent = j
                    best_area = other.area
            parents.append(best_parent)
        return parents


def validate_pattern_fits_lattice(pattern: Pattern, lattice) -> None:
    """Category 4 target 4.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): unit-cell
    bounds policy for shapes that could overlap their own periodic images.

    **What is *not* checked here, because it does not need checking**: a
    shape whose footprint merely crosses a conceptual "cell edge" (e.g. a
    circle centered near `x=0` whose footprint pokes into `x<0`) is already
    handled correctly by the existing analytic-Fourier-transform machinery
    with **no code change and no explicit wrapping logic**. `pattern_epsilon_hat`
    evaluates each shape's continuous-plane Fourier transform at the
    lattice's discrete reciprocal-vector points; by the Poisson summation
    formula, sampling a bounded function's continuous Fourier transform at
    reciprocal-lattice frequencies is *exactly* the Fourier-series
    coefficient of that function's periodic tiling at those lattice
    vectors -- periodicity is a property of *where the Fourier transform is
    evaluated* (the reciprocal lattice points), not something that needs a
    separate real-space wrap/modulo step. Verified directly, not just
    argued, by `tests/test_unit_cell_bounds.py`'s edge-crossing test: a
    from-scratch raster reference that explicitly tiles a shape's periodic
    images and this project's unmodified analytic coefficient agree for a
    shape deliberately placed straddling a cell boundary.

    **What genuinely does need a policy, and is checked here**: the Poisson-
    summation argument above requires the shape not to overlap its own
    periodic copies (a self-intersecting periodic tiling is not the same
    physical pattern as the single shape drawn once) -- exactly the
    "must not have any shapes that intersect with each other" invariant
    `S4/S4/pattern/pattern.h`'s own module docstring (lines 21-35) states
    but explicitly does **not** check ("Currently, the no-intersection
    criterion is not checked due to its complexity, so the input shapes
    must be sanitized elsewhere"). This project chooses to check a
    conservative, cheap sufficient condition instead of leaving it entirely
    unchecked, per `rules.md`'s "fail loud, fail early" convention: reject
    if any shape's `2 * bounding_radius` is not strictly smaller than the
    shorter of the lattice's two primitive vector lengths (`|a|`, `|b|`,
    ignoring a zero vector, e.g. `Lattice1D.b`). This is a **conservative,
    not exhaustive**, check -- it only compares against the two primitive
    vectors, not the full first Brillouin-zone-adjacent shell (`a-b`,
    `a+b`, etc.), so a genuinely self-overlapping shape on a strongly
    oblique lattice could in principle slip through undetected; documented
    honestly as a known limitation rather than claimed as fully general,
    since a tighter check would need the lattice's reduced (shortest-vector)
    basis, not implemented here.
    """
    lattice_vectors = [v for v in (lattice.a, lattice.b) if float(np.hypot(*v)) > 0]
    min_period = min(float(np.hypot(*v)) for v in lattice_vectors)
    for shape in pattern.shapes:
        if 2 * shape.bounding_radius >= min_period:
            raise ValueError(
                f"shape {shape!r} has 2*bounding_radius={2 * shape.bounding_radius!r} >= "
                f"the lattice's shortest primitive vector length ({min_period!r}); it could "
                "overlap its own periodic images, which pattern_epsilon_hat's "
                "Fourier-factorization does not support (see validate_pattern_fits_lattice's "
                "docstring)"
            )
