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
from sougata_solver.fields import propagate_amplitudes, tangential_e_field, z_poynting_flux
from sougata_solver.fourier_basis import truncate_fourier_orders, truncate_fourier_orders_1d
from sougata_solver.fourier_factorization import toeplitz_matrix, toeplitz_matrix_component
from sougata_solver.geometry import Lattice, Lattice1D, validate_pattern_fits_lattice
from sougata_solver.layer import Layer, LayerEigenmodes, LayerStack
from sougata_solver.materials import Material
from sougata_solver.smatrix import SMatrixStack, interior_amplitudes


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
    thicknesses: list[float]
    kx: np.ndarray
    ky: np.ndarray

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

    def layer_absorption(self) -> list[float]:
        """Category 7 targets 7.5/7.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        per-interior-layer absorbed power fraction (normalized to incident
        power, same convention `reflectance()`/`transmittance()` use), one
        entry per interior (finite-thickness) layer in stack order.

        Composition of already-validated Phase 7 pieces only -- no new
        physics formula (`design.md`'s "Layer-Wise Absorption Design",
        `decisions.md` ADR-017): net z-Poynting flux entering a layer at
        its top interface (`smatrix.interior_amplitudes` at
        `partial_smatrix_up_to(i)`) minus net flux leaving at its bottom
        interface (`fields.propagate_amplitudes` through the layer's own
        thickness, then `fields.z_poynting_flux` again). `z_poynting_flux`
        already folds the forward/backward interference cross-term
        symmetrically into its `(forward, backward)` split
        (`fields.py:51-56`), so `forward + backward` is the genuine net
        total z-flux at one reference plane.
        """
        omega = self.excitation.omega()
        n2 = 2 * self.num_orders
        stack = SMatrixStack(self.thicknesses, self.all_modes)

        modes_inc = self.all_modes[0]
        zeros = np.zeros_like(self.a0)
        incident_power, _ = z_poynting_flux(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, self.a0, zeros)

        absorptions: list[float] = []
        for i in range(1, len(self.all_modes) - 1):
            modes_i = self.all_modes[i]
            a_top, b_top = interior_amplitudes(stack.partial_smatrix_up_to(i), n2, self.a0, self.b_reflected)
            a_bot, b_bot = propagate_amplitudes(modes_i.q, self.thicknesses[i], a_top, b_top)

            fwd_top, bwd_top = z_poynting_flux(omega, modes_i.q, modes_i.kp, modes_i.phi, a_top, b_top)
            fwd_bot, bwd_bot = z_poynting_flux(omega, modes_i.q, modes_i.kp, modes_i.phi, a_bot, b_bot)
            net_top = (fwd_top + bwd_top).real
            net_bot = (fwd_bot + bwd_bot).real
            absorptions.append(float((net_top - net_bot) / incident_power.real))
        return absorptions

    def complex_amplitudes(self) -> dict[tuple[int, int], dict[str, complex]]:
        """Category 10 target 10.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
        per-order complex reflected/transmitted tangential E-field
        amplitudes, keyed like `diffraction_efficiencies()`.

        **Convention, explicit per target 10.1's own wording**: these are
        raw Cartesian `(Ex, Ey)` components (`fields.tangential_e_field`,
        already cited/validated), in the same absolute scale as
        `excitation.incident_field_xy()` -- i.e. dividing a returned value
        by the incident field's corresponding component gives a genuine
        complex reflection/transmission coefficient. Deliberately **not**
        re-expressed in an s/p amplitude basis: `tangential_e_field` is
        linear in the modal amplitudes, so (unlike `diffraction_efficiencies()`,
        which is bilinear and needs per-order masking) it can be evaluated
        once on the full `b_reflected`/`a_transmitted` vectors and read off
        per order directly. Validated directly against
        `tests/oracles/fresnel.py::multilayer_complex_rt` at oblique
        incidence for **both** s- and p-polarization -- the measured ratio
        (this order's field component divided by the corresponding
        incident field component) matches that independent, from-scratch
        oracle's complex `r`/`t` to full double precision for both.

        This is a genuine independent-oracle match, but **not** by itself
        the "externally validated" bar target 10.5 requires for exposing
        an *s/p-basis* amplitude conversion: `fresnel.py` is this
        project's own from-scratch derivation (its own docstring: "not
        derived from EMpy or sougata_solver"), not a third-party source
        like S4 or EMpy. See `CONVENTIONS.md`'s Category 10 addendum for
        the full account, including a documented, non-obvious finding —
        `fresnel.py`'s admittance-based p-polarization sign convention
        already agrees with this solver's, even though a *naively hand-
        written* textbook `r_p = (n2*cos(ti) - n1*cos(tt)) / (n2*cos(ti) +
        n1*cos(tt))` formula would appear to disagree in sign (a genuine,
        pre-existing convention-choice ambiguity in p-polarization Fresnel
        formulas generally, not a bug in either this solver or the
        oracle).
        """
        omega = self.excitation.omega()
        modes_inc = self.all_modes[0]
        modes_trans = self.all_modes[-1]
        zeros_inc = np.zeros_like(self.a0)
        ex_r, ey_r = tangential_e_field(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros_inc, self.b_reflected)
        zeros_trans = np.zeros_like(self.a_transmitted)
        ex_t, ey_t = tangential_e_field(omega, modes_trans.q, modes_trans.kp, modes_trans.phi, self.a_transmitted, zeros_trans)

        amplitudes: dict[tuple[int, int], dict[str, complex]] = {}
        for i in range(self.num_orders):
            key = (int(self.g[i, 0]), int(self.g[i, 1]))
            amplitudes[key] = {
                "Ex_r": complex(ex_r[i]),
                "Ey_r": complex(ey_r[i]),
                "Ex_t": complex(ex_t[i]),
                "Ey_t": complex(ey_t[i]),
            }
        return amplitudes

    def diffraction_angles(self) -> dict[tuple[int, int], dict[str, float | None]]:
        """Category 10 target 10.2: per-order diffraction angles (radians),
        keyed like `diffraction_efficiencies()`.

        `theta` is `None` for an evanescent order on that side (reusing
        `eigenmodes.classify_propagating`'s already-validated real/
        imaginary-`q` branch classification, Category 1 target 1.8) --
        never a fabricated angle for a non-propagating order, per this
        target's own "clear non-propagating representation" wording.
        `phi` (the in-plane propagation azimuth) is well-defined
        regardless of propagating status -- it is derived purely from the
        real `(kx, ky)` grating-equation values, not from `q` -- so it is
        always reported.

        Formula (standard grating-equation geometry, not a new physics
        claim): for a propagating order, `q` *is* that order's z-wavenumber
        in the relevant half-space (real, `>=0` by `_select_q_branch`'s
        outgoing-branch convention), so `theta = atan2(sqrt(kx^2+ky^2),
        Re(q))` gives the propagation angle from `+z` directly -- no
        separate refractive-index lookup needed, since `q` already encodes
        it. `phi = atan2(ky, kx)`.
        """
        modes_inc = self.all_modes[0]
        modes_trans = self.all_modes[-1]
        prop_inc = classify_propagating(modes_inc.q)
        prop_trans = classify_propagating(modes_trans.q)

        angles: dict[tuple[int, int], dict[str, float | None]] = {}
        # kx/ky are physically real (in-plane wavevector components); they
        # carry a complex dtype only because `n_incidence` is typed complex
        # upstream (`excitation.k_parallel`) even for a lossless incidence
        # medium -- `.real` here matches how every other angle-independent
        # per-order quantity in this project already treats them.
        kx_real = self.kx.real
        ky_real = self.ky.real
        kr = np.sqrt(kx_real**2 + ky_real**2)
        phi = np.arctan2(ky_real, kx_real)
        for i in range(self.num_orders):
            key = (int(self.g[i, 0]), int(self.g[i, 1]))
            theta_r = float(np.arctan2(kr[i], modes_inc.q[i].real)) if prop_inc[i] else None
            theta_t = float(np.arctan2(kr[i], modes_trans.q[i].real)) if prop_trans[i] else None
            angles[key] = {
                "theta_r": theta_r,
                "phi_r": float(phi[i]),
                "theta_t": theta_t,
                "phi_t": float(phi[i]),
            }
        return angles

    def energy_balance(self) -> dict[str, float]:
        """Category 10 target 10.3: incident/reflected/transmitted/absorbed
        power and the residual, in one dict -- pure composition of already-
        validated methods (`reflectance()`, `transmittance()`,
        `layer_absorption()`), no new formula. `residual` is
        `1 - (reflected + transmitted + absorbed)`, which is exactly zero
        (to numerical precision) whenever energy is genuinely conserved,
        and the "conservation check where physics permits" this category's
        exit criterion requires.
        """
        r = self.reflectance()
        t = self.transmittance()
        absorbed = sum(self.layer_absorption())
        return {
            "incident": 1.0,
            "reflected": r,
            "transmitted": t,
            "absorbed": absorbed,
            "residual": 1.0 - (r + t + absorbed),
        }


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
        # Category 7 targets 7.3/7.4 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md):
        # instance-scoped Toeplitz-matrix cache, gated on a measured timing
        # case per `rules.md`'s Performance Requirements exception clause --
        # see `design.md`'s "Layer/Toeplitz Caching Design" and
        # `decisions.md` ADR-016. Never module-level (rules.md's "no hidden
        # global state" rule): a fresh `Simulation` starts with an empty
        # cache, so no entry can outlive the instance it belongs to.
        self._toeplitz_cache: dict = {}
        # Category 13 target 13.3 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md):
        # instance-scoped per-layer eigenmode cache, implementing the
        # design Category 12 target 12.3 flagged but deliberately left
        # unimplemented (`design.md`'s "Linear-Algebra Baseline &
        # Factorization-Reuse Design"). Same "no hidden global state"
        # rationale as `_toeplitz_cache` above.
        self._eigenmode_cache: dict = {}

        # Category 4 target 4.2 (COMMERCIAL_RCWA_ATOMIC_TARGETS.md): reject
        # a patterned layer whose shapes could overlap their own periodic
        # images at construction time, not deep inside solve() -- see
        # geometry.validate_pattern_fits_lattice's docstring for the policy.
        for layer in self.layer_stack:
            if layer.pattern is not None:
                validate_pattern_fits_lattice(layer.pattern, lattice)

    def _cached_toeplitz(self, pattern, g: np.ndarray, wavelength: float, inverse: bool) -> np.ndarray:
        """Category 7 target 7.4: `id(pattern)`-keyed cache -- see
        `design.md`'s "Layer/Toeplitz Caching Design" for the full key/
        invalidation rationale (`decisions.md` ADR-016)."""
        key = ("toeplitz", id(pattern), wavelength, inverse)
        if key not in self._toeplitz_cache:
            self._toeplitz_cache[key] = toeplitz_matrix(pattern, self.lattice, g, wavelength, inverse=inverse)
        return self._toeplitz_cache[key]

    def _cached_toeplitz_component(
        self, pattern, g: np.ndarray, wavelength: float, row: int, col: int
    ) -> np.ndarray:
        key = ("toeplitz_component", id(pattern), wavelength, row, col)
        if key not in self._toeplitz_cache:
            self._toeplitz_cache[key] = toeplitz_matrix_component(pattern, self.lattice, g, wavelength, row, col)
        return self._toeplitz_cache[key]

    def _cached_layer_eigenmodes(
        self, layer: Layer, g: np.ndarray, wavelength: float, omega: complex, kx: np.ndarray, ky: np.ndarray, is_1d: bool
    ) -> LayerEigenmodes:
        """Category 13 target 13.3: cache one layer's `LayerEigenmodes`
        result, keyed by `(id(layer), omega, kx, ky)` -- deliberately
        **not** including `layer.thickness` (which never affects an
        eigenmode solve; only the downstream `propagation_smatrix` step
        consumes it, per `build_kp_matrix`'s signature) and **not**
        including the excitation's polarization amplitudes (which only
        affect the incident-mode-amplitude inversion, `a0`, downstream of
        every eigenmode solve). This makes the cache automatically valid,
        with no extra invalidation logic, across a fixed-wavelength/angle
        `sweep.sweep_polarization` (target 8.4) or a fixed-wavelength/angle
        `sweep.sweep_thickness` (target 8.5) -- both leave `omega`/`kx`/`ky`
        and every `id(layer)` unchanged, even though `sweep_thickness`
        mutates a `Layer`'s `thickness` attribute in place (the same object
        identity, so the cache key is unaffected, correctly).

        `id(layer)` (not `id(layer.material)`/`id(layer.pattern)`
        directly) is used as the object-identity key component, matching
        `_cached_toeplitz`'s established convention (`decisions.md`
        ADR-016) -- since a `Layer`'s `material`/`pattern` reference never
        changes after construction, keying on the `Layer` itself is
        equivalent and keeps this cache's key shape uniform across every
        dispatch branch below (uniform isotropic/diagonal/in-plane, 1D,
        2D isotropic/anisotropic) without needing a different key shape
        per branch.
        """
        key = ("eigenmodes", id(layer), omega, kx.tobytes(), ky.tobytes())
        if key not in self._eigenmode_cache:
            self._eigenmode_cache[key] = self._solve_one_layer_eigenmodes(layer, g, wavelength, omega, kx, ky, is_1d)
        return self._eigenmode_cache[key]

    def _solve_one_layer_eigenmodes(
        self, layer: Layer, g: np.ndarray, wavelength: float, omega: complex, kx: np.ndarray, ky: np.ndarray, is_1d: bool
    ) -> LayerEigenmodes:
        if layer.is_uniform():
            material = layer.material
            if material.is_isotropic:
                eps_scalar = material.epsilon_tensor(wavelength)[0, 0]
                return solve_layer_eigenmodes_uniform(omega, kx, ky, eps_scalar)
            elif material.is_diagonal:
                # Diagonal-tensor uniform layers (Category 1 target 1.3,
                # COMMERCIAL_RCWA_ATOMIC_TARGETS.md) -- no in-plane
                # (eps_xy/eps_yx) or longitudinal (eps_xz/eps_yz/eps_zx/
                # eps_zy) coupling yet; those are separate, not-yet-shipped
                # targets 1.4/1.5.
                eps_tensor = material.epsilon_tensor(wavelength)
                return solve_layer_eigenmodes_uniform_diagonal(
                    omega, kx, ky, eps_tensor[0, 0], eps_tensor[1, 1], eps_tensor[2, 2]
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
                return solve_layer_eigenmodes_uniform_inplane(
                    omega,
                    kx,
                    ky,
                    eps_tensor[0, 0],
                    eps_tensor[0, 1],
                    eps_tensor[1, 0],
                    eps_tensor[1, 1],
                    eps_tensor[2, 2],
                )
        elif is_1d:
            epsilon_hat = self._cached_toeplitz(layer.pattern, g, wavelength, inverse=False)
            epsilon_inv_hat = self._cached_toeplitz(layer.pattern, g, wavelength, inverse=True)
            return solve_layer_eigenmodes_1d(omega, kx, ky, epsilon_hat, epsilon_inv_hat)
        else:
            materials_in_pattern = [layer.pattern.background] + [s.material for s in layer.pattern.shapes]
            if all(m.is_isotropic for m in materials_in_pattern):
                # 2D patterned layers use only the direct-rule Toeplitz --
                # see solve_layer_eigenmodes_patterned's docstring: S4's
                # true-2D, no-polarization-basis closed-form path (the
                # one transcribed here) doesn't consume the inverse-rule
                # Toeplitz at all, unlike the 1D case above.
                epsilon_hat = self._cached_toeplitz(layer.pattern, g, wavelength, inverse=False)
                return solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
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
                exx = self._cached_toeplitz_component(layer.pattern, g, wavelength, 0, 0)
                exy = self._cached_toeplitz_component(layer.pattern, g, wavelength, 0, 1)
                eyx = self._cached_toeplitz_component(layer.pattern, g, wavelength, 1, 0)
                eyy = self._cached_toeplitz_component(layer.pattern, g, wavelength, 1, 1)
                ezz = self._cached_toeplitz_component(layer.pattern, g, wavelength, 2, 2)
                return solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, exx, exy, eyx, eyy, ezz)

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

        all_modes: list[LayerEigenmodes] = [
            self._cached_layer_eigenmodes(layer, g, wavelength, omega, kx, ky, is_1d) for layer in self.layer_stack
        ]

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
            thicknesses=thicknesses,
            kx=kx,
            ky=ky,
        )
