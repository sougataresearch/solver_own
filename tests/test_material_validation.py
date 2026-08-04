"""Category 5 target 5.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): `Material`
construction- and call-time validation -- tensor shape, finite values, and
wavelength-callback output, per `design.md`'s Error Handling conventions
("fail loud, fail early"). Two tiers: construction-time (probe wavelength
`1.0`) and call-time (`epsilon_tensor` re-validates every call, since a
dispersion callable can misbehave only outside the probe wavelength, e.g.
past an interpolation table's domain).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.materials import Material

# ---------------------------------------------------------------------------
# Construction-time validation
# ---------------------------------------------------------------------------


def test_material_accepts_ordinary_scalar():
    Material("m", 2.25)


def test_material_accepts_ordinary_tensor():
    Material("m", np.diag([2.25, 4.0, 3.1]).astype(complex))


@pytest.mark.parametrize("eps", [math.nan, math.inf, complex(math.nan, 0), complex(1.0, math.inf)])
def test_material_rejects_non_finite_scalar(eps):
    with pytest.raises(ValueError, match="finite"):
        Material("m", eps)


def test_material_rejects_non_finite_tensor_entry():
    tensor = np.diag([2.25, math.nan, 3.1]).astype(complex)
    with pytest.raises(ValueError, match="finite"):
        Material("m", tensor)


def test_material_rejects_wrong_tensor_shape():
    with pytest.raises(ValueError, match="scalar or a 3x3"):
        Material("m", np.eye(2))


def test_material_rejects_callable_returning_non_finite_at_probe_wavelength():
    with pytest.raises(ValueError, match="finite"):
        Material("m", lambda wavelength: float("nan"))


# ---------------------------------------------------------------------------
# Call-time validation (epsilon_tensor) -- a callable that is fine at the
# probe wavelength (1.0) but misbehaves elsewhere.
# ---------------------------------------------------------------------------


def test_epsilon_tensor_rejects_non_finite_output_away_from_probe_wavelength():
    def eps_fn(wavelength):
        return float("nan") if wavelength != 1.0 else 2.25

    material = Material("m", eps_fn)  # constructs fine (probe wavelength gives 2.25)
    material.epsilon_tensor(1.0)  # still fine
    with pytest.raises(ValueError, match="finite"):
        material.epsilon_tensor(0.5)


def test_epsilon_tensor_rejects_shape_change_away_from_probe_wavelength():
    """A callable that returns a scalar at the probe wavelength (so the
    Material is constructed as isotropic) but a tensor elsewhere -- an
    inconsistency the constructor-time check alone cannot see."""

    def eps_fn(wavelength):
        if wavelength == 1.0:
            return 2.25
        return np.eye(3, dtype=complex) * 4.0

    material = Material("m", eps_fn)
    assert material.is_isotropic
    with pytest.raises(ValueError, match="shape"):
        material.epsilon_tensor(0.5)


def test_epsilon_tensor_rejects_non_finite_wavelength():
    material = Material("m", 2.25)
    with pytest.raises(ValueError, match="finite"):
        material.epsilon_tensor(float("nan"))
    with pytest.raises(ValueError, match="finite"):
        material.epsilon_tensor(float("inf"))


def test_epsilon_tensor_accepts_ordinary_dispersive_callable():
    material = Material("m", lambda wavelength: 2.0 + 0.1 * wavelength)
    assert material.epsilon_tensor(1.0)[0, 0] == pytest.approx(2.1)
    assert material.epsilon_tensor(2.0)[0, 0] == pytest.approx(2.2)
