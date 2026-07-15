"""Layer and layer-stack data model."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from pyrcwa.geometry import Pattern
from pyrcwa.materials import Material


@dataclass
class Layer:
    """One layer in the stack.

    Use `thickness = math.inf` for the semi-infinite incidence/transmission
    half-spaces. A layer is either uniform (`material` set, `pattern` None)
    or patterned (`pattern` set, whose `Pattern.background` supplies the
    background material)."""

    name: str
    thickness: float
    material: Material | None = None
    pattern: Pattern | None = None

    def __post_init__(self):
        if (self.material is None) == (self.pattern is None):
            raise ValueError("Layer requires exactly one of `material` or `pattern`")

    def is_uniform(self) -> bool:
        return self.pattern is None

    def background_material(self) -> Material:
        return self.material if self.is_uniform() else self.pattern.background


@dataclass
class LayerEigenmodes:
    """Result of solving one layer's eigenmode problem at one wavelength."""

    q: np.ndarray            # (2n,) complex z-propagation constants
    phi: np.ndarray           # (2n,2n) complex eigenvectors
    kp: np.ndarray             # (2n,2n) complex k-parallel operator
    epsilon_inv: np.ndarray | None  # (n,n), None when is_scalar_isotropic
    is_scalar_isotropic: bool


class LayerStack:
    """Ordered list of layers, plus incidence/exit half-spaces."""

    def __init__(self, layers: list[Layer], incidence: Material, transmission: Material):
        self.layers = [
            Layer("incidence", math.inf, material=incidence),
            *layers,
            Layer("transmission", math.inf, material=transmission),
        ]

    def __len__(self) -> int:
        return len(self.layers)

    def __iter__(self):
        return iter(self.layers)

    def __getitem__(self, index: int) -> Layer:
        return self.layers[index]
