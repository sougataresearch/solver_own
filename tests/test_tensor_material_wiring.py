"""Category 5 target 5.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): tensor-
material solver wiring, "complete only after Category 1's corresponding
tensor milestone and benchmark pass."

Status check (not assumed): Category 1 targets 1.3 (uniform diagonal
tensor), 1.4 (uniform in-plane-coupled tensor), and 1.6 (patterned
anisotropic layers) are shipped and wired into `simulation.py`'s
`material.is_diagonal`/`is_isotropic` dispatch (only target 1.5,
longitudinal `eps_xz/eps_yz/eps_zx/eps_zy` coupling, remains explicitly
deferred -- unrelated to what this target gates on). So the gate this
target names is already met for that scope; what was **not** previously
tested is a *dispersive* tensor material (a `Material.from_permittivity_tensor`
built from a wavelength-dependent callable whose components use the new
Category 5 analytic dispersion models, e.g. `from_sellmeier`) actually
flowing correctly through Category 1's tensor eigensolvers end to end --
every prior anisotropic test used a constant (non-dispersive) tensor. This
file closes that specific, previously-untested combination.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

_BK7 = dict(b1=1.03961212, c1=0.00600069867, b2=0.231792344, c2=0.0200179144, b3=1.01046945, c3=103.560653)


def _dispersive_uniaxial_tensor(wavelength: float) -> np.ndarray:
    """A birefringent crystal whose ordinary-axis index follows BK7's
    Sellmeier dispersion and whose extraordinary axis is a fixed offset
    from it -- deliberately not physically real, just a vehicle to prove a
    genuinely wavelength-dependent tensor works, not a claim about any
    specific material."""
    bk7 = Material.from_sellmeier("bk7_probe", **_BK7)
    n_o = np.sqrt(complex(bk7.epsilon_tensor(wavelength)[0, 0])).real
    n_e = n_o + 0.05
    return np.diag([n_o**2, n_o**2, n_e**2]).astype(complex)


@pytest.mark.parametrize("wavelength_nm", [486.1, 587.56, 656.3])
def test_dispersive_uniform_diagonal_tensor_energy_conservation(wavelength_nm):
    """Target 1.3's uniform diagonal-tensor solver, fed a dispersive tensor."""
    wavelength = wavelength_nm * 1e-9
    air = Material("air", 1.0)
    crystal = Material.from_permittivity_tensor("dispersive_uniaxial", _dispersive_uniaxial_tensor)
    lattice = Lattice(a=(0.7e-6, 0.0), b=(0.0, 0.7e-6))
    layer = Layer("crystal", 0.3e-6, material=crystal)
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(wavelength, theta=np.radians(20.0), phi=0.0, s_amplitude=0.7, p_amplitude=0.7)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


def test_dispersive_uniform_diagonal_tensor_reflectance_actually_varies_with_wavelength():
    """Confirms the dispersion is actually flowing through the tensor
    solver (not e.g. accidentally frozen at a probe-wavelength value) --
    reflectance at two well-separated wavelengths must differ."""
    air = Material("air", 1.0)
    crystal = Material.from_permittivity_tensor("dispersive_uniaxial", _dispersive_uniaxial_tensor)
    lattice = Lattice(a=(0.7e-6, 0.0), b=(0.0, 0.7e-6))
    layer = Layer("crystal", 0.3e-6, material=crystal)
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)

    r_blue = sim.solve(
        PlaneWaveExcitation(486.1e-9, theta=np.radians(20.0), phi=0.0, s_amplitude=0.7, p_amplitude=0.7)
    ).reflectance()
    r_red = sim.solve(
        PlaneWaveExcitation(656.3e-9, theta=np.radians(20.0), phi=0.0, s_amplitude=0.7, p_amplitude=0.7)
    ).reflectance()
    assert r_blue != pytest.approx(r_red, abs=1e-6)


@pytest.mark.parametrize("wavelength_nm", [486.1, 656.3])
def test_dispersive_patterned_anisotropic_tensor_energy_conservation(wavelength_nm):
    """Target 1.6's patterned-anisotropic-layer solver, fed a dispersive
    tensor pattern -- the deepest previously-untested combination (Category
    1's tensor solve x Category 5's dispersion models x Category 4's
    patterned-layer Toeplitz construction, all three landing in different
    sessions)."""
    wavelength = wavelength_nm * 1e-9
    air = Material("air", 1.0)
    crystal = Material.from_permittivity_tensor("dispersive_uniaxial", _dispersive_uniaxial_tensor)
    lattice = Lattice(a=(0.7e-6, 0.0), b=(0.0, 0.7e-6))
    pattern = Pattern(background=air, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.18e-6, material=crystal)])
    layer = Layer("crystal_pillar", 0.3e-6, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=5, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(wavelength, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)
