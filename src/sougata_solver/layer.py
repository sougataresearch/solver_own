"""Layer and layer-stack data model."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sougata_solver.geometry import Pattern
from sougata_solver.materials import Material


def _require_valid_thickness(thickness: float) -> None:
    """Category 7 target 7.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): fail
    loud, fail early, at construction, same "why" as `geometry._require_finite`/
    `materials._require_finite` -- a NaN or non-positive thickness would
    otherwise only surface later as a nonsensical `propagation_smatrix`
    phase (`exp(1j*q*thickness)`) or a silently-wrong stack. `math.inf` is
    deliberately **not** rejected: it is this class's own documented
    sentinel for a semi-infinite half-space (see the class docstring and
    `LayerStack.__init__`, which constructs the incidence/transmission
    layers with exactly this value); `smatrix.SMatrixStack.__init__` never
    calls `propagation_smatrix` for the first/last layer, so `+inf` is safe
    there specifically and only there.
    """
    if math.isnan(thickness):
        raise ValueError(f"Layer thickness must not be NaN, got {thickness!r}")
    if thickness != math.inf and not (thickness > 0):
        raise ValueError(
            f"Layer thickness must be > 0 (or exactly math.inf for a semi-infinite "
            f"half-space), got {thickness!r}"
        )


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
        _require_valid_thickness(self.thickness)

    def is_uniform(self) -> bool:
        return self.pattern is None

    def background_material(self) -> Material:
        return self.material if self.is_uniform() else self.pattern.background


@dataclass
class EigenmodeDiagnostics:
    """Category 2 target 2.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): eigenvalue/
    mode-conditioning diagnostics reported alongside `LayerEigenmodes.q`/
    `phi`/`kp`, never fed back into the solve itself -- purely a summary of
    already-computed quantities so a caller doing sweep diagnostics or
    investigating a `WARNING` logged by `eigenmodes.py`
    (`ILL_CONDITIONED_THRESHOLD`, `DEGENERATE_GAP_THRESHOLD`) doesn't have to
    recompute the same condition numbers/classification itself.
    """

    cond_epsilon: float          # cond() of the permittivity/Toeplitz matrix
                                  # actually inverted to build kp; exactly 1.0
                                  # for closed-form uniform-isotropic layers,
                                  # which invert nothing (phi = I).
    cond_phi: float               # cond(phi), the eigenvector matrix.
    min_eigenvalue_gap: float     # smallest pairwise |q_i - q_j| among all
                                   # distinct-index eigenvalues (inf if <2 modes).
    num_propagating: int
    num_evanescent: int


@dataclass
class LayerEigenmodes:
    """Result of solving one layer's eigenmode problem at one wavelength."""

    q: np.ndarray            # (2n,) complex z-propagation constants
    phi: np.ndarray           # (2n,2n) complex eigenvectors
    kp: np.ndarray             # (2n,2n) complex k-parallel operator
    epsilon_inv: np.ndarray | None  # (n,n), None when is_scalar_isotropic
    is_scalar_isotropic: bool
    diagnostics: "EigenmodeDiagnostics | None" = None


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
