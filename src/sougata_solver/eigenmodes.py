"""Per-layer eigenmode solver.

Formulas verified directly against `S4/S4/rcwa.cpp` (not paraphrased):
`MakeKPMatrix` (lines 242-293), `SolveLayerEigensystem_uniform`
(lines 422-502), and the general eigen-operator construction in
`SolveLayerEigensystem` (lines 684-827, read in full for Phase 3 — see
`solve_layer_eigenmodes_1d`'s docstring for the general-case citation and
its 1D specialization).
"""

from __future__ import annotations

import logging

import numpy as np

from sougata_solver.layer import LayerEigenmodes

logger = logging.getLogger(__name__)

# Condition-number threshold above which solve_layer_eigenmodes_patterned
# logs a WARNING (design.md's Logging Strategy). Chosen from the Phase 4b
# stress-test sweep documented in that function's docstring: cond(epsilon_hat)
# reached ~900 and cond(phi) reached ~170 in the most adversarial cases tried
# (index contrast down to -20+2j, num_orders up to 225, near-touching shapes,
# sub-percent sliver rectangles) without any loss of accuracy (energy
# conservation and the independent RCWA.jl eigenvalue oracle both held to
# ~1e-10 throughout). 1e4 gives roughly 10x headroom above the worst
# empirically-observed case, rather than being an arbitrary round number, so
# it should only fire for a genuinely worse case than anything already tested.
ILL_CONDITIONED_THRESHOLD = 1e4


def build_kp_matrix(omega: complex, kx: np.ndarray, ky: np.ndarray, epsilon_inv) -> np.ndarray:
    """Build the k-parallel operator:

        kp = omega^2 * I_2n - V @ Einv_block @ V^T

    where `V = [[diag(ky)], [-diag(kx)]]` (2n x n). `epsilon_inv` is either
    a scalar (isotropic uniform layer, `Einv_block = epsilon_inv * I_n`) or
    an `(n, n)` Fourier-space matrix (patterned/anisotropic layer).

    Source: `rcwa.cpp::MakeKPMatrix`, lines 242-293.
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    n = kx.shape[0]
    n2 = 2 * n
    idx = np.arange(n)

    kappa = np.zeros((n2, n2), dtype=complex)
    if np.ndim(epsilon_inv) == 0:
        epsinv = complex(epsilon_inv)
        kappa[idx, idx] = ky * epsinv * ky
        kappa[idx, idx + n] = -ky * epsinv * kx
        kappa[idx + n, idx] = -kx * epsinv * ky
        kappa[idx + n, idx + n] = kx * epsinv * kx
    else:
        einv = np.asarray(epsilon_inv, dtype=complex)
        ky_diag = np.diag(ky)
        kx_diag = np.diag(kx)
        kappa[:n, :n] = ky_diag @ einv @ ky_diag
        kappa[:n, n:] = -ky_diag @ einv @ kx_diag
        kappa[n:, :n] = -kx_diag @ einv @ ky_diag
        kappa[n:, n:] = kx_diag @ einv @ kx_diag

    kp = -kappa
    kp[idx, idx] += omega**2
    kp[idx + n, idx + n] += omega**2
    return kp


def _select_q_branch(q_sq: np.ndarray, tol: float = 4 * np.finfo(float).eps) -> np.ndarray:
    """Select the outgoing/decaying branch of `q = sqrt(q_sq)`, matching
    S4's real-frequency convention (`rcwa.cpp:455-467` / `847-861`):
    near-real `q_sq` is handled specially to keep purely-evanescent modes
    exactly on the positive-imaginary axis; otherwise the principal branch
    is flipped to have `Im(q) >= 0`.
    """
    q_sq = np.asarray(q_sq, dtype=complex)
    q = np.empty_like(q_sq)

    near_real = np.abs(q_sq.imag) <= tol * np.abs(q_sq.real)
    real_part = q_sq.real

    pos = near_real & (real_part >= 0)
    neg = near_real & (real_part < 0)
    q[pos] = np.sqrt(real_part[pos])
    q[neg] = 1j * np.sqrt(-real_part[neg])

    other = ~near_real
    q_other = np.sqrt(q_sq[other])
    flip = q_other.imag < 0
    q_other[flip] = -q_other[flip]
    q[other] = q_other
    return q


def solve_layer_eigenmodes_uniform(omega: complex, kx: np.ndarray, ky: np.ndarray, eps: complex) -> LayerEigenmodes:
    """Closed-form eigenmode solve for a uniform isotropic layer.

    `q[i] = q[i+n] = branch_select(eps*omega^2 - kx[i]^2 - ky[i]^2)`,
    `phi = I_2n` (the eigenbasis coincides with the plane-wave field basis
    for a homogeneous isotropic medium). Verified algebraically consistent
    with the general eigen-operator (`op = Epsilon2 @ kp - U@U^T` reduces
    to `(eps*omega^2 - kx^2 - ky^2) * I_2n` per order when `Epsilon2 = eps*I`).

    Source: `rcwa.cpp::SolveLayerEigensystem_uniform`, lines 422-502.
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    n = kx.shape[0]
    n2 = 2 * n

    q_sq = eps * omega**2 - kx**2 - ky**2
    q_half = _select_q_branch(q_sq)
    q = np.concatenate([q_half, q_half])

    phi = np.eye(n2, dtype=complex)
    kp = build_kp_matrix(omega, kx, ky, 1.0 / eps)

    return LayerEigenmodes(q=q, phi=phi, kp=kp, epsilon_inv=None, is_scalar_isotropic=True)


def solve_layer_eigenmodes_1d(
    omega: complex,
    kx: np.ndarray,
    ky: np.ndarray,
    epsilon_hat: np.ndarray,
    epsilon_inv_hat: np.ndarray,
) -> LayerEigenmodes:
    """Eigenmode solve for a 1D-periodic patterned (grating) layer under
    classical/in-plane mounting (`ky` uniformly `0` for every diffraction
    order -- conical mounting is out of scope for Phase 3 and raises
    `ValueError`).

    This is **not** a separate TE/TM formula transcribed from a
    scalar-RCWA reference -- it is S4's *general* non-uniform eigenoperator
    (`S4/S4/rcwa.cpp::SolveLayerEigensystem`, lines 684-827, read in full
    for this phase), specialized to the case `ky = 0`, where it becomes
    exactly block-diagonal (verified directly by
    `tests/test_1d_grating.py::test_op_is_block_diagonal_when_ky_zero` for a
    random complex Toeplitz input, not just claimed here). Using the general
    operator (rather than an independently-derived scalar TE/TM pair, e.g.
    `Rigorous-Coupled-Wave-Analysis`'s Rumpf-formulation matrices) is
    required because `smatrix.py`/`fields.py`/`excitation.py` already
    assume S4's specific meaning of `(q, phi, kp)`; a differently-normalized
    formulation would silently break `interface_smatrix`/`z_poynting_flux`/
    `tangential_e_field`, none of which are modified for Phase 3.

    Construction (`rcwa.cpp:770-806`, with `Epsilon2`'s per-block source
    below)::

        kp = build_kp_matrix(omega, kx, ky, epsilon_inv_hat)   # unmodified
        Epsilon2 = block_diag(epsilon_hat, inv(epsilon_inv_hat))
        op = Epsilon2 @ kp - [[diag(kx*kx), diag(kx*ky)],
                               [diag(ky*kx), diag(ky*ky)]]

    `Epsilon2`'s bottom-right block being `inv(epsilon_inv_hat)` (the
    matrix-inverse of the *inverse-rule* Toeplitz), not `epsilon_hat`
    directly, is Li's (1996) Fourier-factorization rule for the
    kx-coupled/TM-like field component at a discontinuous interface --
    verified against `S4/S4/fmm/fmm_closed.cpp:77-127`'s "1D proper FFF
    rule" branch (`0 == S->Lr[2] && 0 == S->Lr[3]`), which sets exactly this
    block to `inv(Epsilon_inv)` (see that file's `RNP::LinearSolve<'N'>`
    call at that line range).

    With `ky = 0`, `op` reduces to block-diagonal (top-left `n x n` "TE"
    block using `epsilon_hat` directly, bottom-right `n x n` "TM" block
    using `inv(epsilon_inv_hat)` via `kp`'s bottom-right block), so the two
    blocks are eigendecomposed independently (`numpy.linalg.eig`, general
    non-Hermitian since dispersive/lossy materials are supported
    project-wide, unlike the `eigh`-for-lossless shortcut used in the
    `Rigorous-Coupled-Wave-Analysis` reference) rather than as one dense
    `2n x 2n` solve -- cheaper, and free of the near-degenerate cross-block
    coupling that motivates deferring the dense 2D case to Phase 4a/4b.
    `q`'s branch selection reuses `_select_q_branch` unmodified (its
    docstring already cites `rcwa.cpp:847-861`, inside this same general
    branch).
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    if not np.allclose(ky, 0.0):
        raise ValueError(
            "solve_layer_eigenmodes_1d requires ky == 0 for every order "
            "(classical/in-plane mounting); conical mounting is out of "
            "scope for Phase 3 (1D-periodic lamellar gratings)"
        )
    n = kx.shape[0]
    epsilon_hat = np.asarray(epsilon_hat, dtype=complex)
    epsilon_inv_hat = np.asarray(epsilon_inv_hat, dtype=complex)

    kp = build_kp_matrix(omega, kx, ky, epsilon_inv_hat)

    kx_sq = np.diag(kx * kx)
    op_te = epsilon_hat @ kp[:n, :n] - kx_sq
    einv_block = np.linalg.solve(epsilon_inv_hat, np.eye(n, dtype=complex))
    op_tm = einv_block @ kp[n:, n:]

    q_sq_te, phi_te = np.linalg.eig(op_te)
    q_sq_tm, phi_tm = np.linalg.eig(op_tm)

    q = np.concatenate([_select_q_branch(q_sq_te), _select_q_branch(q_sq_tm)])
    phi = np.zeros((2 * n, 2 * n), dtype=complex)
    phi[:n, :n] = phi_te
    phi[n:, n:] = phi_tm

    return LayerEigenmodes(q=q, phi=phi, kp=kp, epsilon_inv=epsilon_inv_hat, is_scalar_isotropic=False)


def solve_layer_eigenmodes_patterned(
    omega: complex,
    kx: np.ndarray,
    ky: np.ndarray,
    epsilon_hat: np.ndarray,
) -> LayerEigenmodes:
    """General eigenmode solve for a 2D-periodic patterned (non-uniform) layer.

    This is the dense `2n x 2n` non-uniform eigenproblem transcribed from
    `S4/S4/rcwa.cpp::SolveLayerEigensystem` (lines 794-827) **and**
    `S4/S4/fmm/fmm_closed.cpp:109-139,162-163` -- the "ordinary Laurent's
    rule" (`!use_polarization_basis`), true-2D (`Lr[2]!=0 || Lr[3]!=0`)
    branch of `FMMGetEpsilon_ClosedForm`, read in full for this phase after
    an earlier draft of this function incorrectly reused
    `solve_layer_eigenmodes_1d`'s `Epsilon2` construction (that construction
    is *only* valid inside `fmm_closed.cpp`'s separate `0==Lr[2]&&0==Lr[3]`
    branch, lines 110-132 -- the 1D-only "proper FFF rule" special case, not
    the general 2D one). The true-2D, no-polarization-basis branch does
    **not** use Li's (1996) inverse-rule correction at all -- confirmed
    line-by-line: `Epsilon2` is `block_diag(epsilon_hat, epsilon_hat)`
    (`fmm_closed.cpp:135`, a plain copy of the direct-rule block, not
    `inv(epsilon_inv_hat)`), and even the matrix fed into `kp` is
    `inv(epsilon_hat)` -- the numerical matrix-inverse of the *direct*-rule
    Toeplitz (`fmm_closed.cpp:136-137`, `LinearSolve` against the identity)
    -- **not** Phase 2's separately Fourier-factorized inverse-rule Toeplitz
    (`epsilon_inv_hat` / `toeplitz_matrix(..., inverse=True)`), which
    `fmm_closed.cpp` only ever populates from `Pattern_GetFourierTransform`
    on `1/eps` inside the 1D branch (line 119-125) or the
    `use_polarization_basis` branch (not exercised here -- S4 has no such
    option wired into this project, and adding one would be a Phase 4a
    scope expansion not requested). Concretely: this is S4's known,
    documented behavior that 2D closed-form patterns without a polarization
    basis get ordinary Laurent's rule -- a real accuracy limitation at
    sharp 2D discontinuities relative to Phase 3's 1D solver, not a design
    choice made here; correcting it (a proper vectorial/normal-vector 2D
    Fourier-factorization rule) is out of Phase 4a's scope as written in
    `phases.md` and would be a separate, explicitly-requested extension.

    Construction ::

        einv = inv(epsilon_hat)                 # matrix-inverse, not Phase 2's inverse-rule Toeplitz
        kp = build_kp_matrix(omega, kx, ky, einv)
        Epsilon2 = [[epsilon_hat, 0], [0, epsilon_hat]]
        op = Epsilon2 @ kp - [[diag(kx^2),  diag(kx*ky)],
                               [diag(kx*ky), diag(ky^2)]]

    `q`'s branch selection reuses `_select_q_branch` (same convention as
    `solve_layer_eigenmodes_1d` / `solve_layer_eigenmodes_uniform`, cited
    in `_select_q_branch`'s docstring). Note this means `op`'s value at
    `ky=0` does **not** reduce to `solve_layer_eigenmodes_1d`'s result --
    that would only hold if both used the same Fourier-factorization rule,
    and they deliberately don't (Li's rule is 1D-only in S4's own source).

    **Phase 4b: degenerate/ill-conditioned eigenvalue handling.** No
    explicit degenerate-eigenvalue splitting, perturbation, or eigenvector
    re-orthogonalization is implemented here -- `_select_q_branch` already
    handles the outgoing-mode branch-selection side robustly (reused
    unmodified, as everywhere else in this module), and no additional
    handling was found to be *needed*: a deliberate stress-test sweep for
    this phase (index contrast from `n=3.48` to a lossy-metal-like
    `-20+2j`, `num_orders` up to 225, near-touching circular pillars up to
    `radius=0.49*period`, sub-percent-halfwidth sliver rectangles, and
    near-degenerate nested shapes with a `1e-4`-scale radius difference)
    found `cond(epsilon_hat)` reaching ~900 and `cond(phi)` reaching ~170
    in the worst cases tried, with energy conservation and the independent
    `RigorousCoupledWaveAnalysis.jl` eigenvalue oracle
    (`tests/oracles/rcwa_2djl_eigenvalues.py`) both holding to ~1e-10
    throughout (see `tests/test_2d_pillar_stress.py` for the frozen
    versions of these exact cases). This is an honest empirical finding,
    not proof that no pathological case exists -- `numpy.linalg.eig` on a
    genuinely near-defective matrix can still return an inaccurate
    eigenvector basis for *some* input this sweep didn't happen to probe.
    The mitigation actually shipped is **detection, not silent
    correction**: `cond(epsilon_hat)` and `cond(phi)` are both logged at
    `WARNING` (module-level `logger`, per `design.md`'s Logging Strategy)
    whenever either exceeds `ILL_CONDITIONED_THRESHOLD` (`1e4`, chosen with
    ~10x headroom above the worst case actually observed) -- so a future
    case that genuinely does hit numerical trouble surfaces as a visible
    warning to investigate, per `rules.md` AI Coding Rule 2, rather than a
    silently degraded R/T number.
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    n = kx.shape[0]
    epsilon_hat = np.asarray(epsilon_hat, dtype=complex)
    if epsilon_hat.shape != (n, n):
        raise ValueError("epsilon_hat must be an (n, n) matrix for the patterned 2D eigensolver")

    cond_epsilon_hat = np.linalg.cond(epsilon_hat)
    if cond_epsilon_hat > ILL_CONDITIONED_THRESHOLD:
        logger.warning(
            "solve_layer_eigenmodes_patterned: epsilon_hat is ill-conditioned "
            "(cond=%.3e > %.0e); the general 2D eigensolve may be numerically "
            "unreliable for this pattern/num_orders combination",
            cond_epsilon_hat,
            ILL_CONDITIONED_THRESHOLD,
        )

    einv = np.linalg.solve(epsilon_hat, np.eye(n, dtype=complex))
    kp = build_kp_matrix(omega, kx, ky, einv)

    epsilon2 = np.zeros((2 * n, 2 * n), dtype=complex)
    epsilon2[:n, :n] = epsilon_hat
    epsilon2[n:, n:] = epsilon_hat

    coupling = np.zeros((2 * n, 2 * n), dtype=complex)
    coupling[:n, :n] = np.diag(kx * kx)
    coupling[:n, n:] = np.diag(kx * ky)
    coupling[n:, :n] = np.diag(kx * ky)
    coupling[n:, n:] = np.diag(ky * ky)

    op = epsilon2 @ kp - coupling

    q_sq, phi = np.linalg.eig(op)
    q = _select_q_branch(q_sq)

    cond_phi = np.linalg.cond(phi)
    if cond_phi > ILL_CONDITIONED_THRESHOLD:
        logger.warning(
            "solve_layer_eigenmodes_patterned: eigenvector matrix phi is "
            "ill-conditioned (cond=%.3e > %.0e), suggesting near-degenerate "
            "eigenvalues; downstream S-matrix/field quantities built from "
            "phi may be numerically unreliable",
            cond_phi,
            ILL_CONDITIONED_THRESHOLD,
        )

    return LayerEigenmodes(q=q, phi=phi, kp=kp, epsilon_inv=einv, is_scalar_isotropic=False)
