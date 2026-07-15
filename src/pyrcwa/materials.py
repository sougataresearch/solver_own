"""Material permittivity representation.

A :class:`Material` always exposes its permittivity as a 3x3 complex tensor
via :meth:`Material.epsilon_tensor`, even for isotropic materials (`eps * I3`).
This keeps every downstream consumer (Fourier factorization, eigenmode
solver) written against a single tensor-valued interface, with `is_isotropic`
/ `is_diagonal` flags used purely as fast-path hints.
"""

from __future__ import annotations

from typing import Callable, Union

import numpy as np

_I3 = np.eye(3, dtype=complex)

ScalarOrFunc = Union[complex, float, Callable[[float], complex]]
TensorOrFunc = Union[np.ndarray, Callable[[float], np.ndarray]]


def _as_callable(value):
    if callable(value):
        return value
    return lambda wavelength: value


class Material:
    """Isotropic or anisotropic dielectric material.

    Parameters
    ----------
    name:
        Human-readable identifier.
    eps:
        Either a scalar (or callable ``wavelength -> scalar``) for an
        isotropic material, or a 3x3 array (or callable
        ``wavelength -> (3,3) array``) for an anisotropic material.
    """

    def __init__(self, name: str, eps: ScalarOrFunc | TensorOrFunc):
        self.name = name
        sample = eps(1.0) if callable(eps) else eps
        sample_arr = np.asarray(sample)
        if sample_arr.ndim == 0:
            self._kind = "isotropic"
            self._eps_fn = _as_callable(eps)
        elif sample_arr.shape == (3, 3):
            self._kind = "diagonal" if np.allclose(sample_arr, np.diag(np.diagonal(sample_arr))) else "general"
            self._eps_fn = _as_callable(eps)
        else:
            raise ValueError(f"eps must be a scalar or a 3x3 tensor, got shape {sample_arr.shape}")

    def epsilon_tensor(self, wavelength: float) -> np.ndarray:
        """Return the 3x3 complex128 permittivity tensor at `wavelength`."""
        value = self._eps_fn(wavelength)
        if self._kind == "isotropic":
            return complex(value) * _I3
        return np.asarray(value, dtype=complex)

    @property
    def is_isotropic(self) -> bool:
        return self._kind == "isotropic"

    @property
    def is_diagonal(self) -> bool:
        return self._kind in ("isotropic", "diagonal")

    @classmethod
    def from_nk(
        cls,
        name: str,
        n: ScalarOrFunc,
        k: ScalarOrFunc = 0.0,
    ) -> "Material":
        """Build an isotropic material from refractive index n (+ extinction k).

        `eps = (n + i*k)**2`. `n`/`k` may be constants or callables of
        wavelength for dispersive materials.
        """
        n_fn = _as_callable(n)
        k_fn = _as_callable(k)

        def eps_fn(wavelength: float) -> complex:
            nc = complex(n_fn(wavelength), 0) + 1j * k_fn(wavelength)
            return nc * nc

        return cls(name, eps_fn)

    @classmethod
    def from_permittivity_tensor(cls, name: str, eps3x3: TensorOrFunc) -> "Material":
        """Build a material directly from a 3x3 permittivity tensor (or a
        callable of wavelength returning one), e.g. for anisotropic crystals."""
        return cls(name, eps3x3)
