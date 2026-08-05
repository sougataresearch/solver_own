"""Category 13 target 13.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): a safe,
narrowly-scoped vectorized wavelength sweep.

**Scope, deliberately narrow**: only stacks where every layer is uniform
and isotropic (a thin-film/multilayer stack, `num_orders=1` -- the
standard case `structures/thin_film/*.py` already uses, since an
unpatterned layer never diffracts). This is "a single simple sweep"
(target 13.4's own wording), not a general vectorized backend for
patterned/anisotropic layers -- vectorizing the general dense eigensolve
path is explicitly Phase 9 scope (`rules.md`'s Performance Requirements),
gated on this category's own findings, not pulled forward here.

**No new physics formula.** Every function below is the *exact same*
formula already cited and validated elsewhere in this project
(`eigenmodes.solve_layer_eigenmodes_uniform`, `eigenmodes._select_q_branch`,
`eigenmodes.build_kp_matrix`, `smatrix.interface_smatrix`/
`propagation_smatrix`/`star_product`, `excitation.incident_mode_amplitude`),
re-expressed with a leading batch axis (`N` wavelengths) so NumPy's native
batched `@`/`np.linalg.solve` (which already broadcast over stacked
matrices, a standard NumPy feature, not new machinery) replace a Python
loop over wavelength for the expensive matrix operations. The final
flux/R/T step reuses `fields.z_poynting_flux` **unchanged**, in a cheap
per-wavelength loop (its cost is negligible compared to the eigensolve/
cascade this module actually vectorizes) -- keeping the new-formula
surface limited to structural re-expression, not a second independent
implementation of the physics.

Per `rules.md`'s Performance Requirements ("Vectorization work... must not
change any numerical result versus the unvectorized path... add a
regression test comparing both paths"): `tests/test_vectorized_sweep.py`
confirms `sweep_wavelength_vectorized` reproduces `sweep.sweep_wavelength`'s
scalar-loop results to numerical precision on a representative thin-film
fixture.
"""

from __future__ import annotations

import numpy as np

from sougata_solver.eigenmodes import _select_q_branch
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.layer import LayerEigenmodes
from sougata_solver.simulation import Simulation, SimulationResult
from sougata_solver.sweep import SweepResult


def _batched_solve(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """`np.linalg.solve` already broadcasts over stacked `(N, m, m)` /
    `(N, m, k)` matrices -- the batched analogue of `smatrix._solve`
    (which uses `scipy.linalg.lu_factor`/`lu_solve` per-call; SciPy's
    LU routines don't batch, so `numpy.linalg.solve`'s native stacked-
    matrix `gesv`-family LAPACK call is used here instead -- same
    underlying algorithm class, not a different one)."""
    return np.linalg.solve(a, b)


def _batched_uniform_layer_modes(omega: np.ndarray, kx: np.ndarray, ky: np.ndarray, eps: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Batched analogue of `eigenmodes.solve_layer_eigenmodes_uniform`,
    restricted to `num_orders=1` (`n=1`, `n2=2`) -- the only case a
    uniform, unpatterned layer ever needs. `omega`, `kx`, `ky`, `eps`: each
    shape `(N,)`. Returns `(q, phi, kp)` with shapes `(N,2)`, `(N,2,2)`,
    `(N,2,2)`.
    """
    n_wavelengths = omega.shape[0]
    q_sq = eps * omega**2 - kx**2 - ky**2
    q_half = _select_q_branch(q_sq)
    q = np.stack([q_half, q_half], axis=1)

    phi = np.broadcast_to(np.eye(2, dtype=complex), (n_wavelengths, 2, 2)).copy()

    # kp = omega^2 * I_2n - kappa (eigenmodes.build_kp_matrix's exact
    # formula -- the `-kappa` sign and `+= omega**2` on the diagonal are
    # both load-bearing; an earlier draft of this function omitted the
    # omega^2*I term entirely and was caught by the equivalence test this
    # module's own docstring describes, not left in).
    epsinv = 1.0 / eps
    kappa = np.zeros((n_wavelengths, 2, 2), dtype=complex)
    kappa[:, 0, 0] = ky * epsinv * ky
    kappa[:, 0, 1] = -ky * epsinv * kx
    kappa[:, 1, 0] = -kx * epsinv * ky
    kappa[:, 1, 1] = kx * epsinv * kx
    kp = -kappa
    kp[:, 0, 0] += omega**2
    kp[:, 1, 1] += omega**2
    return q, phi, kp


def _batched_interface_smatrix(
    q_l: np.ndarray, phi_l: np.ndarray, kp_l: np.ndarray, q_r: np.ndarray, phi_r: np.ndarray, kp_r: np.ndarray
) -> np.ndarray:
    """Batched analogue of `smatrix.interface_smatrix` (same `Ta`/`Tb`/
    `T2Sblocks` formula, `smatrix.py`'s own citation) -- `_is_trivial_interface`'s
    identity-matrix short-circuit is not needed here (every interface in
    this vectorized path is between two genuinely different materials by
    construction; a trivial interface would just solve to the identity
    matrix anyway, at negligible extra cost for `n2=2`)."""
    n_wavelengths, n2 = q_l.shape
    a_l, a_r = phi_l, phi_r
    b_l = kp_l @ phi_l / q_l[:, None, :]
    b_r = kp_r @ phi_r / q_r[:, None, :]

    p = _batched_solve(a_l, a_r)
    q = _batched_solve(b_l, b_r)
    ta = 0.5 * (p + q)
    tb = 0.5 * (p - q)

    identity = np.broadcast_to(np.eye(n2, dtype=complex), (n_wavelengths, n2, n2))
    s00 = _batched_solve(ta, identity)
    s10 = tb @ s00
    s01 = -s00 @ tb
    s11 = ta + tb @ s01

    s = np.zeros((n_wavelengths, 2 * n2, 2 * n2), dtype=complex)
    s[:, :n2, :n2] = s00
    s[:, :n2, n2:] = s01
    s[:, n2:, :n2] = s10
    s[:, n2:, n2:] = s11
    return s


def _batched_propagation_smatrix(q: np.ndarray, thickness: float) -> np.ndarray:
    """Batched analogue of `smatrix.propagation_smatrix`."""
    n_wavelengths, n2 = q.shape
    phase = np.exp(1j * q * thickness)
    s = np.zeros((n_wavelengths, 2 * n2, 2 * n2), dtype=complex)
    idx = np.arange(n2)
    s[:, idx, idx] = phase
    s[:, idx + n2, idx + n2] = phase
    return s


def _batched_star_product(n2: int, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Batched analogue of `smatrix.star_product`."""
    n_wavelengths = a.shape[0]
    a00, a01 = a[:, :n2, :n2], a[:, :n2, n2:]
    a10, a11 = a[:, n2:, :n2], a[:, n2:, n2:]
    b00, b01 = b[:, :n2, :n2], b[:, :n2, n2:]
    b10, b11 = b[:, n2:, :n2], b[:, n2:, n2:]
    identity = np.broadcast_to(np.eye(n2, dtype=complex), (n_wavelengths, n2, n2))

    t1 = identity - a01 @ b10
    c00 = b00 @ _batched_solve(t1, a00)
    c01 = b00 @ _batched_solve(t1, a01 @ b11) + b01

    t2 = identity - b10 @ a01
    c10 = a10 + a11 @ _batched_solve(t2, b10 @ a00)
    c11 = a11 @ _batched_solve(t2, b11)

    c = np.zeros((n_wavelengths, 2 * n2, 2 * n2), dtype=complex)
    c[:, :n2, :n2] = c00
    c[:, :n2, n2:] = c01
    c[:, n2:, :n2] = c10
    c[:, n2:, n2:] = c11
    return c


def _require_uniform_isotropic_stack(simulation: Simulation) -> None:
    for layer in simulation.layer_stack:
        if not layer.is_uniform():
            raise ValueError(
                f"sweep_wavelength_vectorized only supports uniform (unpatterned) layers, "
                f"but layer {layer.name!r} is patterned -- use sweep.sweep_wavelength instead"
            )
        if not layer.material.is_isotropic:
            raise ValueError(
                f"sweep_wavelength_vectorized only supports isotropic materials, "
                f"but layer {layer.name!r}'s material {layer.material.name!r} is anisotropic -- "
                f"use sweep.sweep_wavelength instead"
            )
    if simulation.num_orders != 1:
        raise ValueError(
            f"sweep_wavelength_vectorized only supports num_orders=1 (the natural case for an "
            f"unpatterned stack, which never diffracts), got num_orders={simulation.num_orders!r}"
        )


def sweep_wavelength_vectorized(
    simulation: Simulation,
    wavelengths: np.ndarray,
    theta: float,
    phi: float,
    s_amplitude: complex,
    p_amplitude: complex,
) -> SweepResult:
    """Category 13 target 13.4: vectorized wavelength sweep for a uniform-
    isotropic-only (thin-film) stack -- see module docstring for scope and
    the "no new physics, just batched re-expression" claim.

    Returns a `SweepResult` with the same shape/fields
    `sweep.sweep_wavelength` returns, so callers (and
    `tests/test_vectorized_sweep.py`'s equivalence check) can compare them
    directly.
    """
    _require_uniform_isotropic_stack(simulation)
    wavelengths = np.asarray(wavelengths, dtype=float)
    n_wavelengths = wavelengths.shape[0]
    omega = 2.0 * np.pi / wavelengths

    layers = list(simulation.layer_stack)
    thicknesses = [layer.thickness for layer in layers]

    incidence_material = layers[0].material
    eps_inc = np.array([complex(incidence_material.epsilon_tensor(w)[0, 0]) for w in wavelengths])
    n_inc = np.sqrt(eps_inc)
    k0 = omega * n_inc
    kx0 = k0 * np.sin(theta) * np.cos(phi)
    ky0 = k0 * np.sin(theta) * np.sin(phi)

    all_q, all_phi, all_kp = [], [], []
    for layer in layers:
        eps = np.array([complex(layer.material.epsilon_tensor(w)[0, 0]) for w in wavelengths])
        q, phi_mat, kp = _batched_uniform_layer_modes(omega, kx0, ky0, eps)
        all_q.append(q)
        all_phi.append(phi_mat)
        all_kp.append(kp)

    n2 = 2
    cumulative = _batched_interface_smatrix(all_q[0], all_phi[0], all_kp[0], all_q[1], all_phi[1], all_kp[1])
    for i in range(1, len(layers) - 1):
        prop = _batched_propagation_smatrix(all_q[i], thicknesses[i])
        cumulative = _batched_star_product(n2, cumulative, prop)
        iface = _batched_interface_smatrix(all_q[i], all_phi[i], all_kp[i], all_q[i + 1], all_phi[i + 1], all_kp[i + 1])
        cumulative = _batched_star_product(n2, cumulative, iface)
    s_full = cumulative

    # Incident field / a0 (excitation.incident_field_xy is wavelength-
    # independent for fixed theta/phi/s_amplitude/p_amplitude, since
    # s_hat/p_hat_xy depend only on the fixed angles).
    s_hat = np.array([-np.sin(phi), np.cos(phi)])
    p_hat_xy = -np.cos(theta) * np.array([np.cos(phi), np.sin(phi)])
    e_xy = s_amplitude * s_hat + p_amplitude * p_hat_xy
    ex0, ey0 = complex(e_xy[0]), complex(e_xy[1])
    u_target = np.array([-ey0, ex0], dtype=complex)
    u_target_batched = np.broadcast_to(u_target, (n_wavelengths, 2))

    rhs0 = omega[:, None] * all_q[0] * u_target_batched
    kp_phi_0 = all_kp[0] @ all_phi[0]
    a0 = _batched_solve(kp_phi_0, rhs0[..., None])[..., 0]

    rhs = np.concatenate([a0, np.zeros((n_wavelengths, n2), dtype=complex)], axis=1)
    out = (s_full @ rhs[..., None])[..., 0]
    a_transmitted = out[:, :n2]
    b_reflected = out[:, n2:]

    results: list[SimulationResult] = []
    for i in range(n_wavelengths):
        modes_inc = LayerEigenmodes(q=all_q[0][i], phi=all_phi[0][i], kp=all_kp[0][i], epsilon_inv=None, is_scalar_isotropic=True)
        modes_trans = LayerEigenmodes(q=all_q[-1][i], phi=all_phi[-1][i], kp=all_kp[-1][i], epsilon_inv=None, is_scalar_isotropic=True)
        all_modes_i = [
            LayerEigenmodes(q=all_q[j][i], phi=all_phi[j][i], kp=all_kp[j][i], epsilon_inv=None, is_scalar_isotropic=True)
            for j in range(len(layers))
        ]
        excitation_i = PlaneWaveExcitation(float(wavelengths[i]), theta, phi, s_amplitude, p_amplitude)
        results.append(
            SimulationResult(
                excitation=excitation_i,
                num_orders=1,
                zeroth_order_index=0,
                g=np.array([[0, 0]]),
                all_modes=all_modes_i,
                a0=a0[i],
                a_transmitted=a_transmitted[i],
                b_reflected=b_reflected[i],
                thicknesses=thicknesses,
                kx=np.array([kx0[i]]),
                ky=np.array([ky0[i]]),
            )
        )

    metadata = {
        "num_orders": 1,
        "truncation": simulation.truncation,
        "theta_rad": theta,
        "phi_rad": phi,
        "s_amplitude": s_amplitude,
        "p_amplitude": p_amplitude,
        "vectorized": True,
    }
    return SweepResult("wavelength", "m", list(wavelengths), results, metadata)
