"""Top-level simulation orchestration.

Phase 1: uniform (unpatterned), isotropic layers. Phase 3 adds 1D-periodic
patterned (grating) layers under a `Lattice1D`. 2D-patterned layers are
implemented in Phase 4a. Phase 6 (anisotropic materials) is in progress:
uniform diagonal-tensor (target 1.3) and uniform in-plane-coupled
(target 1.4) layers are supported; uniform longitudinally-coupled tensors
and any anisotropic patterned layer still raise `NotImplementedError`
naming the specific open target (see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
Category 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sougata_solver.eigenmodes import (
    classify_propagating,
    solve_layer_eigenmodes_1d,
    solve_layer_eigenmodes_patterned,
    solve_layer_eigenmodes_patterned_inplane,
    solve_layer_eigenmodes_uniform,
    solve_layer_eigenmodes_uniform_diagonal,
    solve_layer_eigenmodes_uniform_inplane,
)
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import z_poynting_flux
from sougata_solver.fourier_basis import truncate_fourier_orders, truncate_fourier_orders_1d
from sougata_solver.fourier_factorization import toeplitz_matrix, toeplitz_matrix_component
from sougata_solver.geometry import Lattice, Lattice1D, validate_pattern_fits_lattice
from sougata_solver.layer import Layer, LayerEigenmodes, LayerStack
from sougata_solver.materials import Material
from sougata_solver.smatrix import SMatrixStack


@dataclass
class SimulationResult:
    excitation: PlaneWaveExcitation
    num_orders: int
    zeroth_order_index: int
    g: np.ndarray
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

    def diffraction_efficiencies(self) -> dict[tuple[int, int], tuple[float, float]]:
        """Per-order reflected/transmitted diffraction efficiency, keyed by
        the full 2D reciprocal-lattice index `(g1, g2)`.

        Reuses the already-oracle-validated `z_poynting_flux` once per
        order, with every other order's mode amplitude masked to zero,
        rather than a new per-order flux formula: `all_modes[0]`/`all_modes[-1]`
        are always the uniform semi-infinite incidence/transmission half-spaces
        (`LayerStack.__init__`) regardless of what's patterned in between, and
        `build_kp_matrix`'s scalar-`epsilon_inv` branch has no cross-order
        coupling (each order's `kappa` entries only touch that order's own
        index and its `+n` pair, by inspection of `eigenmodes.build_kp_matrix`'s
        `idx`/`idx+n` indexing) -- so masking-and-re-calling isolates one
        order's flux exactly, reusing validated code instead of adding new
        formula risk. Required per `testing.md`'s Physical-Invariant Testing
        (energy conservation) and the Phase 3 system-test benchmark comparison.
        """
        omega = self.excitation.omega()
        modes_inc = self.all_modes[0]
        modes_trans = self.all_modes[-1]
        zeros = np.zeros_like(self.a0)
        incident_power, _ = z_poynting_flux(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, self.a0, zeros)

        n2 = self.a_transmitted.shape[0]
        n = self.num_orders
        block = n2 // 2
        efficiencies: dict[tuple[int, int], tuple[float, float]] = {}
        for i in range(n):
            idx = [i, i + block]
            b_masked = np.zeros_like(self.b_reflected)
            b_masked[idx] = self.b_reflected[idx]
            _, reflected_power = z_poynting_flux(
                omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, b_masked
            )
            a_masked = np.zeros_like(self.a_transmitted)
            a_masked[idx] = self.a_transmitted[idx]
            transmitted_power, _ = z_poynting_flux(
                omega, modes_trans.q, modes_trans.kp, modes_trans.phi, a_masked, zeros
            )
            de_r = float((-reflected_power / incident_power).real)
            de_t = float((transmitted_power / incident_power).real)
            key = (int(self.g[i, 0]), int(self.g[i, 1]))
            efficiencies[key] = (de_r, de_t)
        return efficiencies

    def order_classification(self) -> dict[tuple[int, int], dict[str, str]]:
        """Category 1 target 1.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        per-order propagating/evanescent classification in the incidence
        and transmission half-spaces, keyed the same way as
        `diffraction_efficiencies()`.

        Uses `eigenmodes.classify_propagating` on `all_modes[0]`/
        `all_modes[-1]`'s `q` (always the uniform semi-infinite half-spaces,
        `LayerStack.__init__`, regardless of what's patterned in between --
        the same invariant `diffraction_efficiencies()`'s docstring already
        relies on). For a uniform isotropic half-space,
        `solve_layer_eigenmodes_uniform` concatenates `q = [q_half, q_half]`
        (the same propagation constant for both polarization branches of a
        given order), so indexing `q[i]` for `i < num_orders` is sufficient
        -- classification does not depend on polarization for a uniform
        half-space.
        """
        modes_inc = self.all_modes[0]
        modes_trans = self.all_modes[-1]
        prop_inc = classify_propagating(modes_inc.q)
        prop_trans = classify_propagating(modes_trans.q)

        classification: dict[tuple[int, int], dict[str, str]] = {}
        for i in range(self.num_orders):
            key = (int(self.g[i, 0]), int(self.g[i, 1]))
            classification[key] = {
                "incidence": "propagating" if prop_inc[i] else "evanescent",
                "transmission": "propagating" if prop_trans[i] else "evanescent",
            }
        return classification


class Simulation:
    """Owns a lattice, layer stack, and Fourier-order truncation; solves for
    reflected/transmitted mode amplitudes given a plane-wave excitation."""

    def __init__(
        self,
        lattice: Lattice | Lattice1D,
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

        # Category 4 target 4.2 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md): reject
        # a patterned layer whose shapes could overlap their own periodic
        # images at construction time, not deep inside solve() -- see
        # geometry.validate_pattern_fits_lattice's docstring for the policy.
        for layer in self.layer_stack:
            if layer.pattern is not None:
                validate_pattern_fits_lattice(layer.pattern, lattice)

    def solve(self, excitation: PlaneWaveExcitation) -> SimulationResult:
        wavelength = excitation.wavelength
        omega = excitation.omega()
        is_1d = isinstance(self.lattice, Lattice1D)

        incidence_material = self.layer_stack[0].material
        eps_inc = incidence_material.epsilon_tensor(wavelength)[0, 0]
        n_inc = np.sqrt(eps_inc)
        kx0, ky0 = excitation.k_parallel(n_inc)

        if is_1d:
            g = truncate_fourier_orders_1d(self.lattice, self.num_orders)
        else:
            g = truncate_fourier_orders(self.lattice, self.num_orders, self.truncation)
        lk = self.lattice.reciprocal_vectors()
        kx = kx0 + 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
        ky = ky0 + 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])

        zeroth_order_index = int(np.flatnonzero((g[:, 0] == 0) & (g[:, 1] == 0))[0])

        all_modes: list[LayerEigenmodes] = []
        for layer in self.layer_stack:
            if layer.is_uniform():
                material = layer.material
                if material.is_isotropic:
                    eps_scalar = material.epsilon_tensor(wavelength)[0, 0]
                    all_modes.append(solve_layer_eigenmodes_uniform(omega, kx, ky, eps_scalar))
                elif material.is_diagonal:
                    # Diagonal-tensor uniform layers (Category 1 target 1.3,
                    # COMMERCIAL_RCWA_ATOMIC_TARGETS.md) -- no in-plane
                    # (eps_xy/eps_yx) or longitudinal (eps_xz/eps_yz/eps_zx/
                    # eps_zy) coupling yet; those are separate, not-yet-shipped
                    # targets 1.4/1.5.
                    eps_tensor = material.epsilon_tensor(wavelength)
                    all_modes.append(
                        solve_layer_eigenmodes_uniform_diagonal(
                            omega, kx, ky, eps_tensor[0, 0], eps_tensor[1, 1], eps_tensor[2, 2]
                        )
                    )
                else:
                    eps_tensor = material.epsilon_tensor(wavelength)
                    longitudinal = eps_tensor[[0, 1, 2, 2], [2, 2, 0, 1]]
                    if np.any(longitudinal != 0):
                        raise NotImplementedError(
                            "Longitudinally-coupled anisotropic uniform layers (nonzero "
                            "eps_xz/eps_yz/eps_zx/eps_zy) require Category 1 target 1.5 "
                            "(COMMERCIAL_RCWA_ATOMIC_TARGETS.md), not yet available"
                        )
                    # In-plane-coupled (eps_xy/eps_yx nonzero, no longitudinal
                    # coupling) uniform layers (Category 1 target 1.4).
                    all_modes.append(
                        solve_layer_eigenmodes_uniform_inplane(
                            omega,
                            kx,
                            ky,
                            eps_tensor[0, 0],
                            eps_tensor[0, 1],
                            eps_tensor[1, 0],
                            eps_tensor[1, 1],
                            eps_tensor[2, 2],
                        )
                    )
            elif is_1d:
                epsilon_hat = toeplitz_matrix(layer.pattern, self.lattice, g, wavelength, inverse=False)
                epsilon_inv_hat = toeplitz_matrix(layer.pattern, self.lattice, g, wavelength, inverse=True)
                all_modes.append(solve_layer_eigenmodes_1d(omega, kx, ky, epsilon_hat, epsilon_inv_hat))
            else:
                materials_in_pattern = [layer.pattern.background] + [s.material for s in layer.pattern.shapes]
                if all(m.is_isotropic for m in materials_in_pattern):
                    # 2D patterned layers use only the direct-rule Toeplitz --
                    # see solve_layer_eigenmodes_patterned's docstring: S4's
                    # true-2D, no-polarization-basis closed-form path (the
                    # one transcribed here) doesn't consume the inverse-rule
                    # Toeplitz at all, unlike the 1D case above.
                    epsilon_hat = toeplitz_matrix(layer.pattern, self.lattice, g, wavelength, inverse=False)
                    all_modes.append(solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat))
                else:
                    # Anisotropic patterned layer (Category 1 target 1.6,
                    # COMMERCIAL_RCWA_ATOMIC_TARGETS.md).
                    has_longitudinal = any(
                        np.any(m.epsilon_tensor(wavelength)[[0, 1, 2, 2], [2, 2, 0, 1]] != 0)
                        for m in materials_in_pattern
                    )
                    if has_longitudinal:
                        raise NotImplementedError(
                            "Longitudinally-coupled anisotropic patterned layers require "
                            "Category 1 target 1.5 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md), not yet available"
                        )
                    exx = toeplitz_matrix_component(layer.pattern, self.lattice, g, wavelength, 0, 0)
                    exy = toeplitz_matrix_component(layer.pattern, self.lattice, g, wavelength, 0, 1)
                    eyx = toeplitz_matrix_component(layer.pattern, self.lattice, g, wavelength, 1, 0)
                    eyy = toeplitz_matrix_component(layer.pattern, self.lattice, g, wavelength, 1, 1)
                    ezz = toeplitz_matrix_component(layer.pattern, self.lattice, g, wavelength, 2, 2)
                    all_modes.append(
                        solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, exx, exy, eyx, eyy, ezz)
                    )

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
            g=g,
            all_modes=all_modes,
            a0=a0,
            a_transmitted=a_transmitted,
            b_reflected=b_reflected,
        )
