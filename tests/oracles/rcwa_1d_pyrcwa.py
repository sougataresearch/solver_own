"""Second, independently-derived 1D-grating RCWA oracle for Phase 3,
supplementing `tests/oracles/rcwa_1d_gaylord.py`.

Hand-transcribed (not imported, per `rules.md` AI Coding Rule 7 -- and
practically necessary here: `REFERENCE/PyRCWA` targets a numpy version
before `np.NAN`/`np.mat` were removed in NumPy 2.0, so it cannot even be
imported in this project's environment without a numpy-global monkeypatch,
which would be a genuinely fragile thing to bake into permanent test
infrastructure) from `REFERENCE/PyRCWA` (github.com/vitamingcheng/PyRCWA,
MIT license, added to `REFERENCE/` at the project owner's explicit
request). Unlike the Gaylord oracle's reduced TE-specific `A = KX2 - E`
operator, `PyRCWA` solves the **general** 2D P/Q eigenoperator (the same
architecture family as this project's own `eigenmodes.py`/`smatrix.py`,
including an explicit "free space" gap-medium reference basis) and
restricts to 1D purely via `harmonics=(m, 0)` truncation -- a structurally
different derivation route reaching the same physics, which is what makes
it a genuinely independent second check rather than a re-verification of
the same math Gaylord's oracle already covers.

**Scoped to normal incidence, TE only** (`alpha=theta=0`, `phi=90deg` in
PyRCWA's own source-angle convention -- confirmed by hand that this
reduces to pure `Ey` excitation, i.e. TE/s-polarization in this project's
own convention, matching `PlaneWaveExcitation(theta=0, s_amplitude=1,
p_amplitude=0)`): PyRCWA's oblique-incidence angle convention (`alpha`
polar from normal, `theta` azimuthal about normal) does not obviously map
onto this project's `(theta, phi)` convention, and resolving that mapping
with confidence was out of scope for this addition -- restricting to
normal incidence (where the two conventions coincide) avoids asserting an
unverified angle-convention equivalence, the same discipline
`rcwa_1d_gaylord.py`'s own sign-convention notes already follow.

Formula source, cited by exact file:line:
- `REFERENCE/PyRCWA/pyrcwa/solver.py::compute_wave_vector` (lines 52-81):
  `Kx`/`Ky`/`Kz`/`KRef`/`KTrn` construction.
- `REFERENCE/PyRCWA/pyrcwa/core.py::RCWAFreeSpace.getPQMatrix` (lines
  24-33): vacuum reference-medium P/Q.
- `REFERENCE/PyRCWA/pyrcwa/core.py::RCWARefSide.solve`/`RCWATrnSide.solve`
  (lines 74-88, 124-138): boundary S-matrices from the incidence/
  transmission half-spaces' own P/Q.
- `REFERENCE/PyRCWA/pyrcwa/core.py::RCWASingleLoop.solve` (lines 189-207):
  the grating layer's P/Q (general, all four convolution-matrix terms)
  and its S-matrix (note the `A`/`B` inversion is on the *layer's* own
  `W`/`V`, the opposite direction from the boundary-layer formulas above
  -- transcribed exactly as split in the source, not unified).
- `REFERENCE/PyRCWA/pyrcwa/matrix.py::solve_PQMatrix`/`star_product` (lines
  5-10, 21-34): eigenmode solve and Redheffer star-product cascade.
- `REFERENCE/PyRCWA/pyrcwa/solver.py::compute_diffraction_efficiency`
  (lines 154-194): polarization vector, S-matrix application, R/T.

The material convolution matrix uses this project's own analytic
step-function Fourier coefficients (`rcwa_1d_gaylord._toeplitz_convolution_matrix`,
imported and reused, not reimplemented) rather than
`REFERENCE/PyRCWA/pyrcwa/material.py::get_convolution`'s FFT-of-rasterized-
profile route -- the same choice `rcwa_1d_gaylord.py` already made for the
same reason (avoids rasterization resolution as an extra free parameter;
both compute the same Toeplitz matrix in the fine-resolution limit).
**Verified directly, not assumed**: a live run of the actual
`REFERENCE/PyRCWA` code (via a local numpy-compatibility monkeypatch, used
only for this one-off verification, never for the permanent oracle here)
at `fft_resolution=(2001,2001)`/`(4001,4001)`/`(8001,8001)` gave
`Total R = 0.912281 / 0.912194 / 0.912151`, monotonically converging
toward `sougata_solver`'s own `0.912109` for the identical fixture --
confirming the analytic-convolution transcription below targets the same
converged answer PyRCWA's own FFT route approaches, not a different one.
`mu=1` everywhere in this fixture, so `mu_conv = mu_conv_inv = I` --
transcribed as such (not genuinely tested against a magnetic case).
"""

from __future__ import annotations

import numpy as np

from oracles.rcwa_1d_gaylord import _toeplitz_convolution_matrix


def solve_te_normal_incidence(
    wavelength: float,
    n_ridge: float,
    n_groove: float,
    fill_factor: float,
    lattice_constant: float,
    thickness: float,
    num_ord: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-order TE diffraction efficiencies `(DE_r, DE_t)` at normal
    incidence, each length `2*num_ord+1`. See module docstring for the
    full citation chain and the normal-incidence-only scoping rationale.
    """
    n = 2 * num_ord + 1
    eps_ridge = complex(n_ridge**2)
    eps_groove = complex(n_groove**2)
    identity = np.eye(n, dtype=complex)

    k0 = 2 * np.pi / wavelength
    # compute_wave_vector, solver.py:52-81 -- alpha=theta=0 => K=[0,0,k0];
    # ky=0 for every harmonic since harmonics=(num_ord, 0) truncates n to {0}.
    indices = np.arange(-num_ord, num_ord + 1)
    kx = -2 * np.pi / lattice_constant * indices  # normalized by k0 below
    kx_n = np.diag(kx / k0).astype(complex)  # "Kx" (already divided by k0, per solver.py:78)
    ky_n = np.zeros((n, n), dtype=complex)  # "Ky"
    k_ref_n = np.diag(-np.sqrt((1.0 + 0j) - np.diag(kx_n) ** 2))  # eps_ref=1
    k_trn_n = np.diag(np.sqrt((n_groove**2 + 0j) - np.diag(kx_n) ** 2))  # eps_trn = n_groove**2 (air on both sides here)
    # solver.py:168 uses the *incident wave's* scalar kz (self.K[2]/k0), not a
    # per-order matrix (self.Kz is computed by compute_wave_vector but never
    # actually read in compute_diffraction_efficiency) -- at normal incidence
    # (alpha=0) this is cos(alpha) = 1.
    kz = 1.0

    # RCWAFreeSpace.getPQMatrix, core.py:24-33
    pq11 = kx_n @ ky_n
    pq12 = identity - kx_n @ kx_n
    pq21 = ky_n @ ky_n - identity
    pq22 = -kx_n @ ky_n
    p_free = q_free = np.block([[pq11, pq12], [pq21, pq22]])
    w_free, v_free, _lambda_free = _solve_pq_matrix(p_free, q_free)

    # RCWARefSide/RCWATrnSide.getPQMatrix (isotropic half-spaces, eps=n**2, mu=1),
    # core.py:59-72/109-122 -- both incidence and transmission are air here.
    w_ref, v_ref = _isotropic_wv(kx_n, ky_n, identity, eps=1.0, mu=1.0)
    w_trn, v_trn = _isotropic_wv(kx_n, ky_n, identity, eps=n_groove**2, mu=1.0)

    # RCWARefSide.solve, core.py:74-88 (FreeSpace inverted)
    a_ref = np.linalg.inv(w_free) @ w_ref + np.linalg.inv(v_free) @ v_ref
    b_ref = np.linalg.inv(w_free) @ w_ref - np.linalg.inv(v_free) @ v_ref
    a_ref_inv = np.linalg.inv(a_ref)
    s11_ref = -a_ref_inv @ b_ref
    s12_ref = 2 * a_ref_inv
    s21_ref = 0.5 * (a_ref - b_ref @ a_ref_inv @ b_ref)
    s22_ref = b_ref @ a_ref_inv

    # RCWATrnSide.solve, core.py:124-138 (FreeSpace inverted)
    a_trn = np.linalg.inv(w_free) @ w_trn + np.linalg.inv(v_free) @ v_trn
    b_trn = np.linalg.inv(w_free) @ w_trn - np.linalg.inv(v_free) @ v_trn
    a_trn_inv = np.linalg.inv(a_trn)
    s11_trn = b_trn @ a_trn_inv
    s12_trn = 0.5 * (a_trn - b_trn @ a_trn_inv @ b_trn)
    s21_trn = 2 * a_trn_inv
    s22_trn = -a_trn_inv @ b_trn

    # RCWASingleLoop.getPQMatrix, core.py:168-187 -- mu_conv = mu_conv_inv = I (nonmagnetic).
    eps_conv = _toeplitz_convolution_matrix(num_ord, fill_factor, eps_ridge, eps_groove)
    eps_conv_inv = np.linalg.inv(eps_conv)  # material.py:75 -- numerical inverse of the direct-rule Toeplitz
    mu_conv = mu_conv_inv = identity

    p11 = kx_n @ eps_conv_inv @ ky_n
    p12 = mu_conv - kx_n @ eps_conv_inv @ kx_n
    p21 = ky_n @ eps_conv_inv @ ky_n - mu_conv
    p22 = -ky_n @ eps_conv_inv @ kx_n
    p_layer = np.block([[p11, p12], [p21, p22]])

    q11 = kx_n @ mu_conv_inv @ ky_n
    q12 = eps_conv - kx_n @ mu_conv_inv @ kx_n
    q21 = ky_n @ mu_conv_inv @ ky_n - eps_conv
    q22 = -ky_n @ mu_conv_inv @ kx_n
    q_layer = np.block([[q11, q12], [q21, q22]])

    w_layer, v_layer, lambda_layer = _solve_pq_matrix(p_layer, q_layer)

    # RCWASingleLoop.solve, core.py:189-207 (layer's own W/V inverted -- opposite
    # direction from the boundary formulas above, transcribed as-is).
    w_layer_inv = np.linalg.inv(w_layer)
    v_layer_inv = np.linalg.inv(v_layer)
    a_layer = w_layer_inv @ w_free + v_layer_inv @ v_free
    b_layer = w_layer_inv @ w_free - v_layer_inv @ v_free
    xi = np.diag(np.exp(-k0 * thickness * np.diag(lambda_layer)))

    a_layer_inv = np.linalg.inv(a_layer)
    inner = np.linalg.inv(a_layer - xi @ b_layer @ a_layer_inv @ xi @ b_layer)
    s11_layer = inner @ (xi @ b_layer @ a_layer_inv @ xi @ a_layer - b_layer)
    s12_layer = inner @ xi @ (a_layer - b_layer @ a_layer_inv @ b_layer)
    s21_layer = s12_layer
    s22_layer = s11_layer

    # matrix.py::star_product, lines 21-34.
    s_a = _star_product(
        (s11_ref, s12_ref, s21_ref, s22_ref),
        (s11_layer, s12_layer, s21_layer, s22_layer),
    )
    s_total = _star_product(s_a, (s11_trn, s12_trn, s21_trn, s22_trn))
    s11_total, s21_total = s_total[0], s_total[2]

    # compute_diffraction_efficiency, solver.py:154-194 -- alpha=theta=0,
    # phi=90deg => EP=[0,-1,0] (pure Ey, TE).
    delta = np.zeros((n, 1), dtype=complex)
    delta[num_ord, 0] = 1.0
    s_inc = np.concatenate([0.0 * delta, -1.0 * delta], axis=0)
    c_inc = np.linalg.inv(w_ref) @ s_inc
    c_ref = s11_total @ c_inc
    c_trn = s21_total @ c_inc
    r_t = w_ref @ c_ref
    t_t = w_trn @ c_trn
    r_x, r_y = r_t[:n], r_t[n:]
    t_x, t_y = t_t[:n], t_t[n:]
    r_z = -np.linalg.inv(k_ref_n) @ (kx_n @ r_x + ky_n @ r_y)
    t_z = -np.linalg.inv(k_trn_n) @ (kx_n @ t_x + ky_n @ t_y)

    r2 = (np.abs(r_x) ** 2 + np.abs(r_y) ** 2 + np.abs(r_z) ** 2).ravel()
    t2 = (np.abs(t_x) ** 2 + np.abs(t_y) ** 2 + np.abs(t_z) ** 2).ravel()
    de_r = np.abs(np.real(np.diag(k_ref_n)) / kz * r2)
    de_t = np.abs(1.0 * np.real(np.diag(k_trn_n)) / kz * t2)  # mu_ref/mu_trn = 1
    return de_r, de_t


def _isotropic_wv(kx_n: np.ndarray, ky_n: np.ndarray, identity: np.ndarray, eps: float, mu: float):
    """`RCWARefSide`/`RCWATrnSide.getPQMatrix`, `core.py:59-72`/`109-122`
    (identical formula in both -- a uniform isotropic half-space's own
    P/Q, `eps`/`mu` from that half-space's `Material`)."""
    pq11 = kx_n @ ky_n
    pq12 = eps * mu * identity - kx_n @ kx_n
    pq21 = ky_n @ ky_n - eps * mu * identity
    pq22 = -ky_n @ kx_n
    p = np.block([[pq11, pq12], [pq21, pq22]]) / eps
    q = np.block([[pq11, pq12], [pq21, pq22]]) / mu
    w, v, _lam = _solve_pq_matrix(p, q)
    return w, v


def _solve_pq_matrix(p: np.ndarray, q: np.ndarray):
    """`matrix.py::solve_PQMatrix`/`filter_matrix`, lines 5-18."""
    pq = p @ q
    real_part = np.real(pq).copy()
    imag_part = np.imag(pq).copy()
    real_part[np.abs(real_part) < 1e-16] = 0
    imag_part[np.abs(imag_part) < 1e-16] = 0
    pq_filtered = real_part + 1j * imag_part
    eigenvalues, w = np.linalg.eig(pq_filtered)
    lam = np.diag(np.sqrt(eigenvalues.astype(complex)))
    v = q @ w @ np.linalg.inv(lam)
    return w, v, lam


def _star_product(smatrix_a: tuple, smatrix_b: tuple) -> tuple:
    """`matrix.py::star_product`, lines 21-34."""
    s11_a, s12_a, s21_a, s22_a = smatrix_a
    s11_b, s12_b, s21_b, s22_b = smatrix_b
    identity = np.eye(s11_a.shape[0], dtype=complex)
    inv_1 = np.linalg.inv(identity - s11_b @ s22_a)
    inv_2 = np.linalg.inv(identity - s22_a @ s11_b)
    s11 = s11_a + s12_a @ inv_1 @ s11_b @ s21_a
    s12 = s12_a @ inv_1 @ s12_b
    s21 = s21_b @ inv_2 @ s21_a
    s22 = s22_b + s21_b @ inv_2 @ s22_a @ s12_b
    return s11, s12, s21, s22
