"""Lattice and in-plane pattern geometry.

Reciprocal-vector convention: `kx`, `ky` passed to `Shape.fourier_transform`
are in cycles per unit length (not angular frequency), i.e. a real-space
phase factor enters as ``exp(i * 2*pi * (kx*x + ky*y))``. This matches S4's
convention (`S4/S4/pattern/pattern.c`), where `Lk` (reciprocal lattice) and
`Lr` (real lattice) satisfy `Lr @ Lk.T == I` with no extra `2*pi` factor
folded into the lattice matrices themselves.
"""

from __future__ import annotations

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
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self._Lr = np.array([self.a, self.b])  # rows are basis vectors

    def reciprocal_vectors(self) -> np.ndarray:
        """Return the 2x2 reciprocal basis `Lk` such that `Lr @ Lk.T == I`
        (no `2*pi` factor; see module docstring for the convention)."""
        return np.linalg.inv(self._Lr).T

    def unit_cell_area(self) -> float:
        return abs(self.a[0] * self.b[1] - self.a[1] * self.b[0])


def _rotate(x: float, y: float, angle: float) -> tuple[float, float]:
    c, s = np.cos(angle), np.sin(angle)
    return c * x + s * y, -s * x + c * y


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


@dataclass
class Circle(Shape):
    center: tuple[float, float]
    radius: float
    material: Material

    @property
    def area(self) -> float:
        return np.pi * self.radius**2

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

    @property
    def area(self) -> float:
        return 4.0 * self.halfwidth[0] * self.halfwidth[1]

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
class Pattern:
    """Ordered list of shapes within one layer; later shapes take precedence
    over earlier ones at overlapping points (matches S4's `parent[]`
    subtraction-rule convention, `pattern.c:938`)."""

    background: Material
    shapes: list[Shape] = field(default_factory=list)

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
