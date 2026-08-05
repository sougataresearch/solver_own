"""Category 13 target 13.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`vectorized.sweep_wavelength_vectorized`. Per `rules.md`'s Performance
Requirements ("Vectorization work... must not change any numerical result
versus the unvectorized path... add a regression test comparing both
paths on at least one existing example before considering a vectorized
path done"), every test here compares against `sweep.sweep_wavelength`'s
already-validated scalar-loop results, not against a fabricated reference.

A first draft of `_batched_uniform_layer_modes` omitted the `omega^2 * I`
term from `eigenmodes.build_kp_matrix`'s `kp = omega^2*I - kappa` formula
entirely -- caught immediately by exactly this kind of equivalence test
(a `LinAlgError: Singular matrix` on the very first real fixture tried,
not a subtle silent wrong-answer), fixed before being trusted. See
`vectorized.py`'s own comment at that line.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import sweep_wavelength
from sougata_solver.vectorized import sweep_wavelength_vectorized

AIR = Material("air", 1.0)


def _multilayer_simulation() -> Simulation:
    lattice = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    si = Material("si", 3.48**2)
    layers = [
        Layer("l1", 0.1e-6, material=Material("sio2", 1.46**2)),
        Layer("l2", 0.05e-6, material=Material("tio2", 2.4**2)),
        Layer("l3", 0.08e-6, material=Material("sio2b", 1.46**2)),
    ]
    return Simulation(lattice, layers, num_orders=1, incidence=AIR, transmission=si)


WAVELENGTHS = np.linspace(0.4e-6, 0.8e-6, 41)


@pytest.mark.parametrize(
    "s_amplitude,p_amplitude",
    [(1.0, 0.0), (0.0, 1.0), (0.7, 0.5), (0.6, 0.6j), (1.0 / math.sqrt(2), 1j / math.sqrt(2))],
)
def test_vectorized_matches_scalar_at_oblique_azimuthal_incidence(s_amplitude, p_amplitude):
    sim = _multilayer_simulation()
    theta, phi = math.radians(25.0), math.radians(35.0)

    scalar = sweep_wavelength(sim, WAVELENGTHS, theta, phi, s_amplitude=s_amplitude, p_amplitude=p_amplitude)
    vectorized = sweep_wavelength_vectorized(sim, WAVELENGTHS, theta, phi, s_amplitude=s_amplitude, p_amplitude=p_amplitude)

    np.testing.assert_allclose(vectorized.reflectance(), scalar.reflectance(), atol=1e-12)
    np.testing.assert_allclose(vectorized.transmittance(), scalar.transmittance(), atol=1e-12)


def test_vectorized_matches_scalar_at_normal_incidence():
    sim = _multilayer_simulation()
    scalar = sweep_wavelength(sim, WAVELENGTHS, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    vectorized = sweep_wavelength_vectorized(sim, WAVELENGTHS, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    np.testing.assert_allclose(vectorized.reflectance(), scalar.reflectance(), atol=1e-12)
    np.testing.assert_allclose(vectorized.transmittance(), scalar.transmittance(), atol=1e-12)


def test_vectorized_matches_scalar_for_lossy_material():
    lattice = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    lossy = Material("lossy", -20 + 2j)
    sim = Simulation(lattice, [Layer("lossy_film", 0.05e-6, material=lossy)], num_orders=1, incidence=AIR, transmission=AIR)
    theta, phi = math.radians(25.0), math.radians(35.0)

    scalar = sweep_wavelength(sim, WAVELENGTHS, theta, phi, s_amplitude=1.0, p_amplitude=0.0)
    vectorized = sweep_wavelength_vectorized(sim, WAVELENGTHS, theta, phi, s_amplitude=1.0, p_amplitude=0.0)
    np.testing.assert_allclose(vectorized.reflectance(), scalar.reflectance(), atol=1e-12)
    np.testing.assert_allclose(vectorized.transmittance(), scalar.transmittance(), atol=1e-12)


def test_vectorized_energy_conservation_holds():
    sim = _multilayer_simulation()
    result = sweep_wavelength_vectorized(sim, WAVELENGTHS, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    rt = result.reflectance() + result.transmittance()
    np.testing.assert_allclose(rt, 1.0, atol=1e-8)


def test_vectorized_rejects_patterned_layer():
    from sougata_solver.geometry import Circle, Pattern

    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    si = Material("si", 3.48**2)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=si)])
    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    with pytest.raises(ValueError, match="patterned"):
        sweep_wavelength_vectorized(sim, WAVELENGTHS, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)


def test_vectorized_rejects_anisotropic_material():
    lattice = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    aniso = Material.from_permittivity_tensor("aniso", np.diag([2.0, 2.5, 2.2]))
    sim = Simulation(lattice, [Layer("aniso_film", 0.1e-6, material=aniso)], num_orders=1, incidence=AIR, transmission=AIR)
    with pytest.raises(ValueError, match="isotropic"):
        sweep_wavelength_vectorized(sim, WAVELENGTHS, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)


def test_vectorized_rejects_num_orders_greater_than_one():
    lattice = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    sim = Simulation(lattice, [Layer("l1", 0.1e-6, material=Material("sio2", 1.46**2))], num_orders=9, incidence=AIR, transmission=AIR)
    with pytest.raises(ValueError, match="num_orders"):
        sweep_wavelength_vectorized(sim, WAVELENGTHS, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
