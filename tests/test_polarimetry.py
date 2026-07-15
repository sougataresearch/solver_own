"""Validation for Jones/Mueller reflection-matrix extraction.

`jones_to_mueller` is checked against two known reference cases (identity,
ideal polarizer). `jones_reflection_matrix` is checked against a
convention-independent physical fact (isotropic media don't couple s/p
polarizations) and against the already-validated scalar reflectance/
transmittance from Phase 1.
"""

import math

import numpy as np
import pytest

from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.geometry import Lattice
from pyrcwa.layer import Layer
from pyrcwa.materials import Material
from pyrcwa.polarimetry import jones_reflection_matrix, jones_to_mueller
from pyrcwa.simulation import Simulation

WAVELENGTH = 0.55e-6


def test_jones_to_mueller_identity():
    mueller = jones_to_mueller(np.eye(2))
    assert mueller == pytest.approx(np.eye(4), abs=1e-12)


def test_jones_to_mueller_ideal_polarizer():
    jones = np.array([[1.0, 0.0], [0.0, 0.0]])
    expected = 0.5 * np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ]
    )
    mueller = jones_to_mueller(jones)
    assert mueller == pytest.approx(expected, abs=1e-12)


def _build_sim(structure):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    layers = [Layer(f"layer{i}", d, material=Material(f"m{i}", complex(n) ** 2)) for i, (n, d) in enumerate(structure)]
    return Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)


@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 45.0])
@pytest.mark.parametrize("structure", [[], [(2.0, 0.1e-6)], [(1.46, 0.55e-6 / (4 * 1.46)), (2.35, 0.1e-6)]])
def test_isotropic_stack_has_no_cross_polarization(theta_deg, structure):
    sim = _build_sim(structure)
    theta = math.radians(theta_deg)
    jones = jones_reflection_matrix(sim, WAVELENGTH, theta, phi=0.0)
    assert jones[0, 1] == pytest.approx(0.0, abs=1e-9)  # rsp
    assert jones[1, 0] == pytest.approx(0.0, abs=1e-9)  # rps


@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 45.0, 70.0])
def test_jones_diagonal_matches_scalar_reflectance(theta_deg):
    structure = [(1.46, 0.55e-6 / (4 * 1.46)), (2.35, 0.08e-6)]
    sim = _build_sim(structure)
    theta = math.radians(theta_deg)
    jones = jones_reflection_matrix(sim, WAVELENGTH, theta, phi=0.0)

    result_s = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0))
    result_p = sim.solve(PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=0.0, p_amplitude=1.0))

    assert abs(jones[0, 0]) ** 2 == pytest.approx(result_s.reflectance(), abs=1e-8)
    assert abs(jones[1, 1]) ** 2 == pytest.approx(result_p.reflectance(), abs=1e-8)


def test_mueller_from_isotropic_stack_reflection_is_diagonal_like():
    """No cross-polarization -> the s<->p coupling terms of the Mueller
    matrix (row/col 1 vs 2) should vanish, leaving a block-diagonal-ish
    structure consistent with a simple diattenuator (no depolarization,
    no retardance-driven S2<->S3 mixing beyond what a pure phase gives)."""
    structure = [(1.46, 0.55e-6 / (4 * 1.46)), (2.35, 0.08e-6)]
    sim = _build_sim(structure)
    theta = math.radians(30.0)
    jones = jones_reflection_matrix(sim, WAVELENGTH, theta, phi=0.0)
    mueller = jones_to_mueller(jones)
    # M01 (S1 in -> S0 out coupling) should be nonzero (diattenuation),
    # but M02, M03 (S2,S3 in -> S0 out) should vanish for a diagonal Jones matrix.
    assert mueller[0, 2] == pytest.approx(0.0, abs=1e-8)
    assert mueller[0, 3] == pytest.approx(0.0, abs=1e-8)
    assert mueller[2, 0] == pytest.approx(0.0, abs=1e-8)
    assert mueller[3, 0] == pytest.approx(0.0, abs=1e-8)
