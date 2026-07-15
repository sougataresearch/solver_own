"""Plane-wave excitation and polarization decomposition."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class PlaneWaveExcitation:
    """An incident plane wave, decomposed into s (TE) and p (TM) amplitudes.

    `theta` is the polar angle from +z (radians), `phi` the azimuthal angle
    in the xy plane (radians). Sign convention for s/p unit vectors is
    internal and self-consistent (not yet matched to S4/EMpy's convention);
    since Phase 1 validates only power quantities (`|amplitude|^2`), the
    convention choice does not affect the Fresnel validation gate.
    """

    wavelength: float
    theta: float
    phi: float
    s_amplitude: complex = 1.0
    p_amplitude: complex = 0.0

    def omega(self) -> float:
        """Vacuum angular wavenumber `2*pi/wavelength` (natural units, c=1),
        matching S4/rcwa.cpp's `omega` convention."""
        return 2.0 * math.pi / self.wavelength

    def k_parallel(self, n_incidence: complex) -> tuple[float, float]:
        """In-plane wavevector `(kx0, ky0)` of the zeroth diffraction order."""
        k0 = self.omega() * n_incidence
        kx0 = k0 * math.sin(self.theta) * math.cos(self.phi)
        ky0 = k0 * math.sin(self.theta) * math.sin(self.phi)
        return kx0, ky0

    def incident_field_xy(self) -> tuple[complex, complex]:
        """Transverse `(Ex, Ey)` of the incident wave from its s/p amplitudes.

        `s_hat = (-sin(phi), cos(phi), 0)`,
        `p_hat_xy = -cos(theta) * (cos(phi), sin(phi))` (transverse part of
        `k_hat x s_hat`).
        """
        s_hat = np.array([-math.sin(self.phi), math.cos(self.phi)])
        p_hat_xy = -math.cos(self.theta) * np.array([math.cos(self.phi), math.sin(self.phi)])
        e_xy = self.s_amplitude * s_hat + self.p_amplitude * p_hat_xy
        return complex(e_xy[0]), complex(e_xy[1])

    def incident_mode_amplitude(self, modes, num_orders: int, zeroth_order_index: int) -> np.ndarray:
        """Forward-mode amplitude vector `a0` (length `2*num_orders`) that
        produces the desired incident `(Ex, Ey)` at the zeroth order (and
        zero at all other orders) in a uniform incidence half-space, given
        a pure forward wave (`b = 0`).

        The mode-amplitude-to-field relation is *not* `E = phi @ a` — it is
        (verified directly against `rcwa.cpp::GetInPlaneFieldVector`,
        lines 1959-1995, not the commonly-paraphrased-but-wrong `E=phi@(a+b)`)::

            u = kp @ phi @ a / (omega * q)      # with b = 0
            Ex[i] = u[i + n],  Ey[i] = -u[i]    # note the index swap + sign

        so `a0` must be obtained by inverting this relation:
        `a0 = inv(kp @ phi) @ (omega * q * u_target)`.
        """
        ex, ey = self.incident_field_xy()
        ex_vec = np.zeros(num_orders, dtype=complex)
        ey_vec = np.zeros(num_orders, dtype=complex)
        ex_vec[zeroth_order_index] = ex
        ey_vec[zeroth_order_index] = ey
        u_target = np.concatenate([-ey_vec, ex_vec])

        omega = self.omega()
        rhs = omega * modes.q * u_target
        kp_phi = modes.kp @ modes.phi
        return np.linalg.solve(kp_phi, rhs)
