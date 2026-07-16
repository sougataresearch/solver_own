"""Top-level simulation orchestration.

Phase 1 scope only: uniform (unpatterned), isotropic layers. Patterned
layers (Fourier factorization) and anisotropic materials are deferred to
Phase 2+ and raise `NotImplementedError` here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import z_poynting_flux
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer, LayerEigenmodes, LayerStack
from sougata_solver.materials import Material
from sougata_solver.smatrix import SMatrixStack


@dataclass
class SimulationResult:
    excitation: PlaneWaveExcitation
    num_orders: int
    zeroth_order_index: int
    all_modes: list[LayerEigenmodes]
    a0: np.ndarray
    a_transmitted: np.ndarray
    b_reflected: np.ndarray

    def reflectance(self) -> float:
        """`z_poynting_flux`'s `backward` output is a genuinely-signed z-Poynting
        flux (negative for a wave actually travelling in -z), not a
        pre-negated "reported positive power" -- verified by an independent
        direct E/H (`Sz = 0.5*Re(Ex*conj(Hy) - Ey*conj(Hx))`) computation
        against the known Fresnel reflectance for a bare interface. Hence
        the explicit negation here to report a positive reflectance."""
        omega = self.excitation.omega()
        modes_inc = self.all_modes[0]
        zeros = np.zeros_like(self.a0)
        incident_power, _ = z_poynting_flux(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, self.a0, zeros)
        _, reflected_power = z_poynting_flux(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, self.b_reflected)
        return (-reflected_power / incident_power).real

    def transmittance(self) -> float:
        omega = self.excitation.omega()
        modes_inc = self.all_modes[0]
        modes_trans = self.all_modes[-1]
        zeros_inc = np.zeros_like(self.a0)
        incident_power, _ = z_poynting_flux(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, self.a0, zeros_inc)
        zeros_trans = np.zeros_like(self.a_transmitted)
        transmitted_power, _ = z_poynting_flux(
            omega, modes_trans.q, modes_trans.kp, modes_trans.phi, self.a_transmitted, zeros_trans
        )
        return (transmitted_power / incident_power).real


class Simulation:
    """Owns a lattice, layer stack, and Fourier-order truncation; solves for
    reflected/transmitted mode amplitudes given a plane-wave excitation."""

    def __init__(
        self,
        lattice: Lattice,
        layers: list[Layer],
        num_orders: int,
        incidence: Material,
        transmission: Material,
        truncation: str = "circular",
    ):
        self.lattice = lattice
        self.layer_stack = LayerStack(layers, incidence, transmission)
        self.num_orders = num_orders
        self.truncation = truncation

    def solve(self, excitation: PlaneWaveExcitation) -> SimulationResult:
        wavelength = excitation.wavelength
        omega = excitation.omega()

        incidence_material = self.layer_stack[0].material
        eps_inc = incidence_material.epsilon_tensor(wavelength)[0, 0]
        n_inc = np.sqrt(eps_inc)
        kx0, ky0 = excitation.k_parallel(n_inc)

        g = truncate_fourier_orders(self.lattice, self.num_orders, self.truncation)
        lk = self.lattice.reciprocal_vectors()
        kx = kx0 + 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
        ky = ky0 + 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])

        zeroth_order_index = int(np.flatnonzero((g[:, 0] == 0) & (g[:, 1] == 0))[0])

        all_modes: list[LayerEigenmodes] = []
        for layer in self.layer_stack:
            if not layer.is_uniform():
                raise NotImplementedError("Patterned layers require Phase 2+ Fourier factorization")
            material = layer.material
            if not material.is_isotropic:
                raise NotImplementedError("Anisotropic layers require Phase 4")
            eps_scalar = material.epsilon_tensor(wavelength)[0, 0]
            all_modes.append(solve_layer_eigenmodes_uniform(omega, kx, ky, eps_scalar))

        thicknesses = [layer.thickness for layer in self.layer_stack]
        stack = SMatrixStack(thicknesses, all_modes)
        s_full = stack.full_smatrix()

        a0 = excitation.incident_mode_amplitude(all_modes[0], self.num_orders, zeroth_order_index)
        n2 = 2 * self.num_orders
        rhs = np.concatenate([a0, np.zeros(n2, dtype=complex)])
        out = s_full @ rhs
        a_transmitted = out[:n2]
        b_reflected = out[n2:]

        return SimulationResult(
            excitation=excitation,
            num_orders=self.num_orders,
            zeroth_order_index=zeroth_order_index,
            all_modes=all_modes,
            a0=a0,
            a_transmitted=a_transmitted,
            b_reflected=b_reflected,
        )
