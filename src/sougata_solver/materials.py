"""Material permittivity representation.

A :class:`Material` always exposes its permittivity as a 3x3 complex tensor
via :meth:`Material.epsilon_tensor`, even for isotropic materials (`eps * I3`).
This keeps every downstream consumer (Fourier factorization, eigenmode
solver) written against a single tensor-valued interface, with `is_isotropic`
/ `is_diagonal` flags used purely as fast-path hints.
"""

from __future__ import annotations

import re
from typing import Callable, Union

import numpy as np
from scipy.interpolate import interp1d

_I3 = np.eye(3, dtype=complex)

ScalarOrFunc = Union[complex, float, Callable[[float], complex]]
TensorOrFunc = Union[np.ndarray, Callable[[float], np.ndarray]]

_WAVELENGTH_UNIT_SCALE = {"um": 1e-6, "nm": 1e-9, "m": 1.0}


def _as_callable(value):
    if callable(value):
        return value
    return lambda wavelength: value


def _split_fields(line: str) -> list[str]:
    """Split a data row on commas if present, else on whitespace/tabs --
    lets the same parser read comma-delimited .csv and whitespace/tab
    -delimited .txt files without the caller specifying a delimiter."""
    fields = line.split(",") if "," in line else line.split()
    return [f for f in fields if f]


def _read_numeric_blocks(path: str) -> list[list[list[float]]]:
    """Group a file's numeric rows into blocks, splitting wherever a
    non-numeric (header) line appears. Blank lines are ignored (not
    required as separators) -- a header row is what actually marks a new
    block, matching the refractiveindex.info export layout (an 'n' block,
    a header, then a 'k' block).
    """
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    blocks: list[list[list[float]]] = []
    current: list[list[float]] = []
    for line in lines:
        try:
            values = [float(x) for x in _split_fields(line)]
        except ValueError:
            if current:
                blocks.append(current)
                current = []
            continue
        current.append(values)
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"No numeric data found in {path!r}")
    return blocks


def _wavelength_n_k_from_blocks(blocks: list[list[list[float]]], name: str, path: str):
    """Turn parsed numeric blocks into (wl_n, n_vals, wl_k, k_vals).

    Auto-detects between two on-disk layouts:
      * one block, 3 columns/row -- (wavelength, n, k) together, the
        single-table layout used by commercial RCWA tools.
      * one block, 2 columns/row -- (wavelength, n) only, no k data in the
        file. k is treated as zero (lossless) and a warning is printed so
        this is never silently indistinguishable from a format mistake.
      * two (or more) blocks, 2 columns/row each -- refractiveindex.info's
        style: an 'n' block then a 'k' block (possibly on different
        wavelength grids). Block order is n-then-k.
    """
    if len(blocks) == 1:
        arr = np.array(blocks[0])
        ncols = arr.shape[1]
        if ncols == 3:
            return arr[:, 0], arr[:, 1], arr[:, 0], arr[:, 2]
        if ncols == 2:
            print(
                f"WARNING: no k data found for '{name}' ({path}) -- treating as "
                "lossless (k=0), please double-check your file"
            )
            return arr[:, 0], arr[:, 1], arr[:, 0], np.zeros(arr.shape[0])
        raise ValueError(f"Expected 2 or 3 columns per row in {path!r}, got {ncols}")

    n_block, k_block = np.array(blocks[0]), np.array(blocks[1])
    if n_block.shape[1] != 2 or k_block.shape[1] != 2:
        raise ValueError(
            f"Expected 2-column (wavelength, value) blocks in {path!r}, "
            f"got shapes {n_block.shape} and {k_block.shape}"
        )
    return n_block[:, 0], n_block[:, 1], k_block[:, 0], k_block[:, 1]


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
    def from_nk_file(cls, name: str, path: str, wavelength_unit: str = "um") -> "Material":
        """Build a dispersive Material from an n,k data file on disk.

        Auto-detects the on-disk layout -- see `_wavelength_n_k_from_blocks`
        for the exact rules -- so both a single-table (wavelength, n, k)
        file (comma- or whitespace/tab-delimited, the layout commercial
        RCWA tools use) and refractiveindex.info's two-block CSV export
        work through this one function without the caller having to state
        which layout they're giving it.

        `wavelength_unit` is the unit of the file's wavelength column
        ("um", "nm", or "m"); sougata_solver works in meters internally.
        """
        blocks = _read_numeric_blocks(path)
        wl_n, n_vals, wl_k, k_vals = _wavelength_n_k_from_blocks(blocks, name, path)
        scale = _WAVELENGTH_UNIT_SCALE[wavelength_unit]

        n_interp = interp1d(wl_n * scale, n_vals, bounds_error=False, fill_value="extrapolate")
        k_interp = interp1d(wl_k * scale, k_vals, bounds_error=False, fill_value="extrapolate")

        return cls.from_nk(
            name,
            lambda wavelength: float(n_interp(wavelength)),
            lambda wavelength: float(k_interp(wavelength)),
        )

    @classmethod
    def from_refractiveindex_formula_file(
        cls, name: str, path: str, wavelength_unit: str = "um"
    ) -> "Material":
        """Build a dispersive, lossless Material from a refractiveindex.info
        "formula 4" YAML entry (as vendored under
        NK_FILE/refractiveindex.info-database/database/data/.../nk/*.yml).

        Only formula type 4 is supported -- the type used by e.g. TiO2's
        Devore-o.yml -- since that's the only one needed so far; any other
        `type:` raises `NotImplementedError` naming the unsupported type
        rather than guessing at its coefficient layout.

        Formula 4 (RefractiveIndex.INFO), per
        `NK_FILE/refractiveindex.info-database/database/doc/Dispersion
        formulas.pdf` (2014-06-29), with wavelength `w` in micrometers and
        1-indexed coefficients C1..C17 (0 where absent from the file):

            n^2 = C1
                + C2 * w^C3 / (w^2 - C4^C5)      (omitted if C2 == 0)
                + C6 * w^C7 / (w^2 - C8^C9)      (omitted if C6 == 0)
                + C10*w^C11 + C12*w^C13 + C14*w^C15 + C16*w^C17

        These entries carry no absorption data, so k is always 0 (lossless)
        -- this is what the source itself provides, not a placeholder.
        Wavelengths outside the file's stated `wavelength_range` are
        extrapolated without warning, matching `from_nk_file`'s
        `interp1d(..., fill_value="extrapolate")` convention (both are
        also called once with a nonphysical probe wavelength by
        `Material.__init__`'s isotropic/anisotropic shape check).
        """
        with open(path) as f:
            text = f.read()

        type_match = re.search(r"^\s*-?\s*type:\s*(.+)$", text, re.MULTILINE)
        if type_match is None or type_match.group(1).strip() != "formula 4":
            found = type_match.group(1).strip() if type_match else "<none found>"
            raise NotImplementedError(
                f"from_refractiveindex_formula_file only supports 'formula 4', "
                f"got {found!r} in {path!r}"
            )

        coeff_match = re.search(r"coefficients:\s*(.+)", text)
        if coeff_match is None:
            raise ValueError(f"No coefficients found in {path!r}")
        c = [float(x) for x in coeff_match.group(1).split()]
        c += [0.0] * (17 - len(c))  # pad missing trailing coefficients with 0

        scale = _WAVELENGTH_UNIT_SCALE[wavelength_unit]

        def n_fn(wavelength: float) -> float:
            w = wavelength / scale  # back to micrometers, the formula's native unit
            w2 = w * w
            n2 = c[0]
            if c[1] != 0:
                n2 += c[1] * w**c[2] / (w2 - c[3] ** c[4])
            if c[5] != 0:
                n2 += c[5] * w**c[6] / (w2 - c[7] ** c[8])
            n2 += c[9] * w**c[10] + c[11] * w**c[12] + c[13] * w**c[14] + c[15] * w**c[16]
            return float(np.sqrt(n2))

        return cls.from_nk(name, n_fn, 0.0)

    @classmethod
    def from_permittivity_tensor(cls, name: str, eps3x3: TensorOrFunc) -> "Material":
        """Build a material directly from a 3x3 permittivity tensor (or a
        callable of wavelength returning one), e.g. for anisotropic crystals."""
        return cls(name, eps3x3)
