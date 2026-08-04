"""Material permittivity representation.

A :class:`Material` always exposes its permittivity as a 3x3 complex tensor
via :meth:`Material.epsilon_tensor`, even for isotropic materials (`eps * I3`).
This keeps every downstream consumer (Fourier factorization, eigenmode
solver) written against a single tensor-valued interface, with `is_isotropic`
/ `is_diagonal` flags used purely as fast-path hints.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Callable, Union

import numpy as np
from scipy.interpolate import interp1d

_I3 = np.eye(3, dtype=complex)

# Photon-energy/wavelength conversion used by from_lorentz/from_drude/
# from_drude_lorentz (Category 5 targets 5.4-5.6): E[eV] = _EV_NM / lambda[nm].
# Matches the literal constant used by the vendored
# RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl:16
# (`ω=1239.8/λ`) exactly, rather than a higher-precision h*c/e value, for
# citation fidelity with the oscillator coefficients transcribed from that
# same file (RAKIC_GOLD etc., below the Material class).
_EV_NM = 1239.8


@dataclass(frozen=True)
class LorentzOscillator:
    """One Lorentz-oscillator term in a Rakic-style Lorentz-Drude (LD)
    metal model (Category 5 target 5.6). Field names/order match
    `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl`'s
    `Oscillator` struct exactly (`f::Float64; Γ::Float64; ω::Float64`) --
    transcribed positionally, not renamed, so a reader comparing against
    that source doesn't have to remember a field reordering.
    """

    f: float
    gamma_ev: float
    omega_ev: float

ScalarOrFunc = Union[complex, float, Callable[[float], complex]]
TensorOrFunc = Union[np.ndarray, Callable[[float], np.ndarray]]

_WAVELENGTH_UNIT_SCALE = {"um": 1e-6, "nm": 1e-9, "m": 1.0}


def _require_finite(label: str, value: float) -> None:
    """Category 5 target 5.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): same
    "fail loud, fail early" convention as `geometry._require_finite` -- a
    separate copy, not a shared import, since `geometry.py` and
    `materials.py` are each meant to be usable independently (neither
    imports the other; `Material` is consumed by `geometry.Shape`, not the
    reverse) and this is a two-line function, not worth a new shared-utility
    module for.
    """
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite, got {value!r}")


def _require_finite_array(arr: np.ndarray, label: str) -> None:
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must be finite everywhere, got {arr!r}")


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
    source:
        Category 5 target 5.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        optional free-text provenance/citation string (e.g. "Rakić et al.,
        Appl. Opt. 37, 5271 (1998)", or a refractiveindex.info URL) --
        purely informational, never read by any solver code, stored as the
        public `Material.source` attribute. Every `from_*` classmethod
        below accepts and forwards it. Threading it into a run's serialized
        output is one more `output_paths.write_run_metadata(...,
        <name>_source=material.source)` keyword argument -- that function
        already accepts arbitrary `**params`, so no change to
        `output_paths.py` itself was needed; see
        `tests/test_material_provenance.py` for a worked example.
    """

    def __init__(self, name: str, eps: ScalarOrFunc | TensorOrFunc, source: str | None = None):
        """Category 5 target 5.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        validated at construction (probe wavelength `1.0`) *and* on every
        `epsilon_tensor` call thereafter -- a probe-wavelength check alone
        cannot catch a dispersion callable that only misbehaves (wrong
        shape, non-finite output) at a *different* wavelength, e.g. outside
        an interpolation table's domain, so `epsilon_tensor` re-validates
        every call rather than trusting the constructor-time sample. Fails
        loud, fails early, per `design.md`'s Error Handling conventions.
        """
        self.name = name
        self.source = source
        sample = eps(1.0) if callable(eps) else eps
        sample_arr = np.asarray(sample)
        _require_finite_array(sample_arr, f"Material {name!r}: eps")
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
        _require_finite("Material.epsilon_tensor wavelength", wavelength)
        value = self._eps_fn(wavelength)
        value_arr = np.asarray(value)
        expected_shape = () if self._kind == "isotropic" else (3, 3)
        if value_arr.shape != expected_shape:
            raise ValueError(
                f"Material {self.name!r}: eps callable returned shape {value_arr.shape} at "
                f"wavelength={wavelength!r}, expected {expected_shape or 'scalar'} "
                f"(constructed as {self._kind!r})"
            )
        _require_finite_array(value_arr, f"Material {self.name!r}: eps at wavelength={wavelength!r}")
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
        source: str | None = None,
    ) -> "Material":
        """Build an isotropic material from refractive index n (+ extinction k).

        `eps = (n + i*k)**2`. `n`/`k` may be constants or callables of
        wavelength for dispersive materials. `source`: see target 5.8's
        note on `Material.__init__`.
        """
        n_fn = _as_callable(n)
        k_fn = _as_callable(k)

        def eps_fn(wavelength: float) -> complex:
            nc = complex(n_fn(wavelength), 0) + 1j * k_fn(wavelength)
            return nc * nc

        return cls(name, eps_fn, source=source)

    @classmethod
    def from_nk_file(cls, name: str, path: str, wavelength_unit: str = "um", source: str | None = None) -> "Material":
        """Build a dispersive Material from an n,k data file on disk.

        Auto-detects the on-disk layout -- see `_wavelength_n_k_from_blocks`
        for the exact rules -- so both a single-table (wavelength, n, k)
        file (comma- or whitespace/tab-delimited, the layout commercial
        RCWA tools use) and refractiveindex.info's two-block CSV export
        work through this one function without the caller having to state
        which layout they're giving it.

        `wavelength_unit` is the unit of the file's wavelength column
        ("um", "nm", or "m"); sougata_solver works in meters internally.
        `source` defaults to `None`; passing `path` itself is a reasonable
        choice if the file has no more specific citation (target 5.8).
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
            source=source,
        )

    @classmethod
    def from_refractiveindex_formula_file(
        cls, name: str, path: str, wavelength_unit: str = "um", source: str | None = None
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

        return cls.from_nk(name, n_fn, 0.0, source=source)

    @classmethod
    def from_sellmeier(
        cls,
        name: str,
        b1: float,
        c1: float,
        b2: float,
        c2: float,
        b3: float,
        c3: float,
        wavelength_unit: str = "um",
        source: str | None = None,
    ) -> "Material":
        """Category 5 target 5.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        build a dispersive, lossless Material from the standard 3-term
        Sellmeier equation::

            n(w)^2 = 1 + B1*w^2/(w^2-C1) + B2*w^2/(w^2-C2) + B3*w^2/(w^2-C3)

        with `w` in `wavelength_unit` (default `"um"`, matching how
        Sellmeier coefficients are conventionally tabulated -- e.g. the
        SCHOTT optical glass catalog and refractiveindex.info both quote
        `C1/C2/C3` in `um^2`) and `B1..C3` dimensionless/`um^2` accordingly.

        Transcribed directly from the vendored
        `EMpy/EMpy/materials.py::RefractiveIndex.__from_sellmeier` (lines
        118-127, a 3-term Sellmeier form used by that module's own `BK7`
        preset) -- not re-derived, since the equation's well-known form
        (Sellmeier 1871) is exactly reproduced there and using a source
        already exercised by another RCWA-adjacent project reduces the risk
        of a transcription slip in the `w^2/(w^2-C_i)` term. Like
        `from_refractiveindex_formula_file`, this carries no absorption
        data, so `k` is always `0` (lossless) -- what the Sellmeier
        equation itself provides, not a placeholder.

        Validated against BK7's published Sellmeier coefficients and its
        independently-published `n_d = 1.5168` at the Fraunhofer d-line
        (587.56 nm) -- both confirmed via `WebSearch` this session (not
        transcribed from memory), see `tests/test_dispersion_models.py`.
        """
        scale = _WAVELENGTH_UNIT_SCALE[wavelength_unit]

        def n_fn(wavelength: float) -> float:
            w = wavelength / scale
            w2 = w * w
            n2 = 1.0 + b1 * w2 / (w2 - c1) + b2 * w2 / (w2 - c2) + b3 * w2 / (w2 - c3)
            return float(np.sqrt(n2))

        return cls.from_nk(name, n_fn, 0.0, source=source)

    @classmethod
    def from_cauchy(
        cls,
        name: str,
        a: float,
        b: float,
        c: float = 0.0,
        wavelength_unit: str = "um",
        source: str | None = None,
    ) -> "Material":
        """Category 5 target 5.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        build a dispersive, lossless Material from the (2- or 3-term)
        Cauchy equation::

            n(w) = A + B/w^2 + C/w^4

        with `w` in `wavelength_unit` (default `"um"`). Transcribed from
        the vendored `EMpy/EMpy/materials.py::RefractiveIndex`
        docstring's own worked Cauchy-form example (lines 65-72, a `SiN`
        dispersion function `1.887 + 0.01929/x^2 + 1.6662e-4/x^4` with `x`
        in micrometers) -- the standard Cauchy (1836) form, `C` optional
        (defaults to `0.0`, the classic 2-term form). No absorption data,
        `k` is always `0` (lossless), same as `from_sellmeier`.
        """
        scale = _WAVELENGTH_UNIT_SCALE[wavelength_unit]

        def n_fn(wavelength: float) -> float:
            w = wavelength / scale
            return a + b / w**2 + c / w**4

        return cls.from_nk(name, n_fn, 0.0, source=source)

    @classmethod
    def from_lorentz(
        cls,
        name: str,
        eps_inf: float,
        strength: float,
        omega0_ev: float,
        gamma_ev: float,
        source: str | None = None,
    ) -> "Material":
        """Category 5 target 5.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        build a Material from a single-resonance Lorentz oscillator::

            eps(w) = eps_inf + strength / (omega0_ev^2 - w^2 - i*gamma_ev*w)

        with photon energy `w = 1239.8 / wavelength_nm` [eV] (see `_EV_NM`),
        `omega0_ev`/`gamma_ev` the resonance energy/damping rate in eV, and
        `strength` in eV^2 (this is a single already-scaled parameter,
        equal to `f * omega_p^2` in the `f`/`omega_p`-split form
        `from_drude_lorentz`'s multi-oscillator composition uses below --
        a standalone single-oscillator API has no second oscillator to
        share a common `omega_p` with, so there is nothing to gain from
        splitting `strength` into two parameters here).

        Per-oscillator term transcribed from the vendored
        `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl:19`
        (`epsilon += o.f*m.ωp^2/(o.ω^2-ω^2-1im*ω*o.Γ)`), citing A. D.
        Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, "Optical
        properties of metallic films for vertical-cavity optoelectronic
        devices," Appl. Opt. 37, 5271-5283 (1998) -- `strength` here is
        that source's `f*ωp²` combined into one parameter.

        **Causality/sign-convention check, independently re-derived and
        tested (`tests/test_dispersion_models.py`), not assumed from the
        transcribed source alone**: this project's documented `d/dt ->
        -i*omega` phasor convention (`CONVENTIONS.md`) requires a passive
        (lossy, not gain) medium to have `Im(eps) > 0` -- the exact sign
        convention Category 2 target 2.5 found a naively-reused `n=-20+2j`
        "lossy metal" index actually violated (see that target's stress-
        regression fixture for the full account). For a damped classical
        oscillator driven by `exp(-i*omega*t)`, at resonance
        (`omega=omega0_ev`) the denominator reduces to `-i*gamma_ev*omega0_ev`
        exactly, so `eps - eps_inf = strength / (-i*gamma_ev*omega0_ev) =
        i*strength/(gamma_ev*omega0_ev)` -- purely positive-imaginary for
        `strength, gamma_ev, omega0_ev > 0`, confirming Rakic's Julia sign
        (`-1im*ω*o.Γ`, not `+1im*ω*o.Γ`) already matches this project's own
        convention, rather than assuming it does.
        """

        def eps_fn(wavelength: float) -> complex:
            wavelength_nm = wavelength / _WAVELENGTH_UNIT_SCALE["nm"]
            w = _EV_NM / wavelength_nm
            return eps_inf + strength / (omega0_ev**2 - w**2 - 1j * gamma_ev * w)

        return cls(name, eps_fn, source=source)

    @classmethod
    def from_drude(
        cls, name: str, eps_inf: float, omega_p_ev: float, gamma_ev: float, source: str | None = None
    ) -> "Material":
        """Category 5 target 5.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        build a Material from the free-electron Drude model::

            eps(w) = eps_inf - omega_p_ev^2 / (w * (w + i*gamma_ev))

        with `w = 1239.8 / wavelength_nm` [eV] (same convention as
        `from_lorentz`/`from_drude_lorentz`), `omega_p_ev` the plasma
        energy, `gamma_ev` the collision (damping) rate, both in eV.

        Cross-checked between **two independently vendored sources**
        agreeing on the same formula structure: `Rigorous-Coupled-Wave-Analysis/TMM_examples/TMM_Drude.py:67`
        (`drude_eps = 1 - omega_p**2/(omega**2 + 1j*omega*gamma)`, SI rad/s
        units, `eps_inf` fixed at `1`) and
        `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl:14-17`
        (`epsilon=1 .-(Ωp^2)/ω/(ω+1im*m.Γ0)`, eV units) -- algebraically
        identical (`w*(w+i*gamma) == w^2+i*w*gamma`), confirmed directly
        rather than assumed, before trusting either as the transcription
        source. This function follows `rakic.jl`'s eV convention (for
        direct compatibility with `from_drude_lorentz`'s oscillator terms)
        and generalizes `eps_inf` to a free parameter (both vendored
        sources fix it to `1`; `from_drude_lorentz` below also defaults to
        `1.0`, matching Rakic exactly, but exposes it as an override).

        Validated against Rakic's own published, tabulated Au/Ag/Al/Ti
        Drude-term coefficients (`RAKIC_GOLD` etc., below) via
        `from_drude_lorentz` (which calls this same formula for its Drude
        term) -- see `tests/test_dispersion_models.py`.
        """

        def eps_fn(wavelength: float) -> complex:
            wavelength_nm = wavelength / _WAVELENGTH_UNIT_SCALE["nm"]
            w = _EV_NM / wavelength_nm
            return eps_inf - omega_p_ev**2 / (w * (w + 1j * gamma_ev))

        return cls(name, eps_fn, source=source)

    @classmethod
    def from_drude_lorentz(
        cls,
        name: str,
        omega_p_ev: float,
        f0: float,
        gamma0_ev: float,
        oscillators: "tuple[LorentzOscillator, ...]",
        eps_inf: float = 1.0,
        source: str | None = None,
    ) -> "Material":
        """Category 5 target 5.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        Rakic's Lorentz-Drude (LD) composition -- a Drude term plus any
        number of Lorentz-oscillator terms sharing one common plasma
        energy `omega_p_ev`::

            Omega_p = sqrt(f0) * omega_p_ev
            eps(w) = eps_inf - Omega_p^2/(w*(w+i*gamma0_ev))
                     + sum_k oscillators[k].f * omega_p_ev^2
                           / (oscillators[k].omega_ev^2 - w^2 - i*w*oscillators[k].gamma_ev)

        Transcribed directly from
        `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl:14-21`
        (`LorentzDrude`) -- the Drude term reuses `from_drude`'s formula
        with an effective plasma energy `Omega_p = sqrt(f0)*omega_p_ev`
        (`rakic.jl:15`), each oscillator term matches `from_lorentz`'s
        formula with `strength = oscillators[k].f * omega_p_ev^2` (the
        split form `from_lorentz`'s standalone API intentionally avoids,
        needed here because every oscillator shares the same `omega_p_ev`).
        `eps_inf` defaults to `1.0`, matching Rakic's own convention
        exactly (`rakic.jl:17`'s literal `1 .-`) -- not extrapolated.

        `oscillators=()` (zero terms) reduces exactly to `from_drude`,
        confirmed by `tests/test_dispersion_models.py` (target 5.6's
        "zero-strength terms reduce correctly" requirement) -- more
        precisely, "zero *oscillators*"; a `LorentzOscillator` with
        `f=0.0` is algebraically equivalent (contributes exactly `0` to
        the sum for any finite `omega_ev/gamma_ev`), also tested.

        `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM`
        (module-level, below) are `rakic.jl:24-45`'s published, tabulated
        coefficients, transcribed verbatim (same numbers, same field
        order) -- the "published or tabulated reference" this target and
        target 5.5 validate against, citing A. D. Rakić, A. B. Djurišić,
        J. M. Elazar, and M. L. Majewski, "Optical properties of metallic
        films for vertical-cavity optoelectronic devices," Appl. Opt. 37,
        5271-5283 (1998).
        """
        omega_p_effective = math.sqrt(f0) * omega_p_ev

        def eps_fn(wavelength: float) -> complex:
            wavelength_nm = wavelength / _WAVELENGTH_UNIT_SCALE["nm"]
            w = _EV_NM / wavelength_nm
            eps = eps_inf - omega_p_effective**2 / (w * (w + 1j * gamma0_ev))
            for osc in oscillators:
                eps += osc.f * omega_p_ev**2 / (osc.omega_ev**2 - w**2 - 1j * w * osc.gamma_ev)
            return eps

        return cls(name, eps_fn, source=source)

    @classmethod
    def from_permittivity_tensor(cls, name: str, eps3x3: TensorOrFunc, source: str | None = None) -> "Material":
        """Build a material directly from a 3x3 permittivity tensor (or a
        callable of wavelength returning one), e.g. for anisotropic crystals."""
        return cls(name, eps3x3, source=source)


# Category 5 targets 5.5/5.6 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md): published,
# tabulated Lorentz-Drude (LD) metal-model coefficients, transcribed
# verbatim (same numbers, same field order -- (f, gamma_ev, omega_ev) per
# LorentzOscillator, (omega_p_ev, f0, gamma0_ev) per metal) from
# RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl:24-45, citing
# A. D. Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, "Optical
# properties of metallic films for vertical-cavity optoelectronic devices,"
# Appl. Opt. 37, 5271-5283 (1998) -- the same citation that file's own
# trailing comment (rakic.jl:51) gives. Pass to Material.from_drude_lorentz
# as the omega_p_ev/f0/gamma0_ev/oscillators arguments, e.g.
# Material.from_drude_lorentz("Au", *RAKIC_GOLD).
RAKIC_ALUMINUM = (
    14.98,
    0.523,
    0.047,
    (
        LorentzOscillator(0.227, 0.333, 0.162),
        LorentzOscillator(0.050, 0.312, 1.544),
        LorentzOscillator(0.166, 1.351, 1.188),
        LorentzOscillator(0.030, 3.382, 3.473),
    ),
)
RAKIC_SILVER = (
    9.010,
    0.845,
    0.048,
    (
        LorentzOscillator(0.065, 3.886, 0.816),
        LorentzOscillator(0.124, 0.452, 4.481),
        LorentzOscillator(0.011, 0.065, 8.185),
        LorentzOscillator(0.840, 0.916, 9.083),
        LorentzOscillator(5.646, 2.419, 20.290),
    ),
)
RAKIC_GOLD = (
    9.030,
    0.760,
    0.053,
    (
        LorentzOscillator(0.024, 0.241, 0.415),
        LorentzOscillator(0.010, 0.345, 0.830),
        LorentzOscillator(0.071, 0.870, 2.969),
        LorentzOscillator(0.601, 2.494, 4.304),
        LorentzOscillator(4.384, 2.214, 13.320),
    ),
)
RAKIC_TITANIUM = (
    7.290,
    0.148,
    0.082,
    (
        LorentzOscillator(0.899, 2.276, 0.777),
        LorentzOscillator(0.393, 2.518, 1.545),
        LorentzOscillator(0.187, 1.663, 2.509),
        LorentzOscillator(0.001, 1.762, 19.430),
    ),
)
