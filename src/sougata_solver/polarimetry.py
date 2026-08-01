"""Jones and Mueller reflection matrices for a non-depolarizing multilayer
stack, at the specular (zeroth) order or per diffraction order.

`jones_reflection_matrix`/`jones_to_mueller` give the specular-order Jones
and Mueller matrices -- the whole reflected field for a uniform
(unpatterned) stack, since it has no other orders to scatter into.
`jones_reflection_matrix_by_order` generalizes this to patterned/grating
layers, which scatter power into non-zero orders that each carry their own
Jones matrix -- pass any one order's matrix through `jones_to_mueller` the
same way.
"""

from __future__ import annotations

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import tangential_e_field


def _s_hat(phi: float) -> np.ndarray:
    return np.array([-math.sin(phi), math.cos(phi)])


def _p_dir(phi: float) -> np.ndarray:
    return np.array([math.cos(phi), math.sin(phi)])


def decompose_sp(ex: complex, ey: complex, phi: float, cos_theta_signed: float) -> tuple[complex, complex]:
    """Decompose transverse `(Ex, Ey)` into `(E_s, E_p)`.

    `cos_theta_signed` is the same signed transverse scale factor used to
    *build* the field in `PlaneWaveExcitation.incident_field_xy`
    (`-cos(theta)` for a wave travelling in `+z` i.e. incident,
    `+cos(theta)` for `-z` i.e. reflected) -- `s_hat` and `p_dir` are
    orthogonal unit vectors, but the actual p-direction vector used to
    build/read the transverse field has magnitude `|cos(theta)|`, so `E_p`
    must be recovered by dividing out that scale factor.

    Public (not `_`-prefixed) so a raw-field-saving "structures" script and
    a separate "postprocessing" script (which builds the Jones/Mueller
    matrix from that saved raw field) can both reuse this exact convention
    instead of re-deriving/duplicating it -- see `postprocessing/jones_mueller_ellipsometry.py`.
    """
    s_hat = _s_hat(phi)
    p_dir = _p_dir(phi)
    e_vec = np.array([ex, ey])
    e_s = np.dot(e_vec, s_hat)
    e_p = np.dot(e_vec, p_dir) / cos_theta_signed
    return e_s, e_p


def jones_reflection_matrix(sim, wavelength: float, theta: float, phi: float) -> np.ndarray:
    """2x2 complex Jones reflection matrix `[[rss, rsp], [rps, rpp]]` at
    the given `(wavelength, theta, phi)`, where `rXY = E_X_reflected /
    E_Y_incident` for `X, Y` in `{s, p}` (row = reflected component,
    column = incident component).

    Obtained by solving the simulation twice: once with pure s-incidence
    (recovers column `rss, rps`), once with pure p-incidence (recovers
    column `rsp, rpp`).

    For an isotropic, non-tilted (unpatterned) stack, `rsp = rps = 0`
    exactly (no cross-polarization coupling) -- a convention-independent
    physical check used to validate this function.
    """
    cos_theta = math.cos(theta)
    jones = np.zeros((2, 2), dtype=complex)
    for column, (s_amp, p_amp) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        excitation = PlaneWaveExcitation(wavelength, theta, phi, s_amplitude=s_amp, p_amplitude=p_amp)
        result = sim.solve(excitation)
        modes_inc = result.all_modes[0]
        omega = excitation.omega()
        zeros = np.zeros_like(result.a0)
        ex, ey = tangential_e_field(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, result.b_reflected)
        i = result.zeroth_order_index
        e_s, e_p = decompose_sp(ex[i], ey[i], phi, cos_theta)  # reflected: +cos(theta)
        jones[0, column] = e_s
        jones[1, column] = e_p
    return jones


def jones_reflection_matrix_by_order(sim, wavelength: float, theta: float, phi: float) -> dict[tuple[int, int], np.ndarray]:
    """Per-diffraction-order 2x2 complex Jones reflection matrix, keyed by
    the same `(g1, g2)` reciprocal-lattice index as
    `SimulationResult.diffraction_efficiencies` (`simulation.py`).

    `rXY = E_X_reflected_at_this_order / E_Y_incident_at_zeroth_order` for
    `X, Y` in `{s, p}` -- incidence is always injected at the zeroth order
    only (`PlaneWaveExcitation.incident_mode_amplitude`), so this is the
    natural per-order generalization of `jones_reflection_matrix`.

    Isolates one order's field by zeroing every other order's mode
    amplitude before calling `tangential_e_field`, the identical technique
    `SimulationResult.diffraction_efficiencies` already uses (and its
    docstring already justifies) to isolate one order's power -- reused
    here rather than re-derived. For a uniform (unpatterned) stack,
    `num_orders == 1` and this returns only the zeroth-order entry, equal
    to `jones_reflection_matrix`'s result; non-zero orders only carry
    power for patterned/grating layers (`num_orders > 1`).
    """
    cos_theta = math.cos(theta)
    jones_by_order: dict[tuple[int, int], np.ndarray] = {}
    for column, (s_amp, p_amp) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        excitation = PlaneWaveExcitation(wavelength, theta, phi, s_amplitude=s_amp, p_amplitude=p_amp)
        result = sim.solve(excitation)
        modes_inc = result.all_modes[0]
        omega = excitation.omega()
        zeros = np.zeros_like(result.a0)
        block = result.b_reflected.shape[0] // 2
        for i in range(result.num_orders):
            key = (int(result.g[i, 0]), int(result.g[i, 1]))
            if key not in jones_by_order:
                jones_by_order[key] = np.zeros((2, 2), dtype=complex)
            b_masked = np.zeros_like(result.b_reflected)
            idx = [i, i + block]
            b_masked[idx] = result.b_reflected[idx]
            ex, ey = tangential_e_field(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, b_masked)
            e_s, e_p = decompose_sp(ex[i], ey[i], phi, cos_theta)  # reflected: +cos(theta)
            jones_by_order[key][0, column] = e_s
            jones_by_order[key][1, column] = e_p
    return jones_by_order


_SIGMA = (
    np.eye(2, dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
)


def jones_to_mueller(jones: np.ndarray) -> np.ndarray:
    """Convert a 2x2 complex Jones matrix to its 4x4 real Mueller matrix.

    Stokes-parameter convention used (in the Jones matrix's own (row0,
    row1) basis, e.g. (s, p)):
        S0 = |E0|^2 + |E1|^2
        S1 = |E0|^2 - |E1|^2
        S2 = 2*Re(E0 * conj(E1))
        S3 = -2*Im(E0 * conj(E1))

    Derived from scratch via `M_ij = 0.5*Tr(sigma_i @ J @ sigma_j @ J^dagger)`
    (not copied from a remembered reference formula) and verified against
    two known cases: `J = I` gives `M = I(4)` (Mueller identity), and
    `J = diag(1, 0)` (an ideal polarizer passing only the first component)
    gives the standard `0.5*[[1,1,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]]`
    polarizer Mueller matrix.
    """
    jones = np.asarray(jones, dtype=complex)
    jones_dag = jones.conj().T
    m = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            m[i, j] = 0.5 * np.trace(_SIGMA[i] @ jones @ _SIGMA[j] @ jones_dag).real
    return m
