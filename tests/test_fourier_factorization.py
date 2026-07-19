import numpy as np

from sougata_solver.fourier_factorization import pattern_epsilon_hat, toeplitz_matrix
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.geometry import Circle, Lattice, Pattern, Rectangle
from sougata_solver.materials import Material


def test_circle_dc_value_equals_area():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.0, 0.0), radius=0.3, material=mat)
    assert np.isclose(circle.fourier_transform(0.0, 0.0), circle.area)


def test_rectangle_dc_value_equals_area():
    mat = Material("core", 4.0)
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.15), material=mat)
    assert np.isclose(rect.fourier_transform(0.0, 0.0), rect.area)


def test_circle_off_center_dc_value_still_equals_area():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.3, -0.1), radius=0.2, material=mat)
    assert np.isclose(circle.fourier_transform(0.0, 0.0), circle.area)


def test_rectangle_contains_matches_geometry():
    mat = Material("core", 4.0)
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.1), material=mat)
    assert rect.contains(0.0, 0.0)
    assert rect.contains(0.19, 0.09)
    assert not rect.contains(0.21, 0.0)
    assert not rect.contains(0.0, 0.11)


def test_circle_contains_matches_geometry():
    mat = Material("core", 4.0)
    circle = Circle(center=(0.0, 0.0), radius=0.25, material=mat)
    assert circle.contains(0.0, 0.0)
    assert circle.contains(0.24, 0.0)
    assert not circle.contains(0.26, 0.0)


# --- Toeplitz-matrix cross-checks against a rasterized-mask reference ------
#
# Reference: for a rectangular lattice a=(Lx,0), b=(0,Ly), the continuous
# Fourier coefficient hat{eps}(G) = (1/A) * integral_cell eps(r) * exp(-2j*pi*G.r) dr
# is approximated by a Riemann sum over an N x N grid covering the cell
# (centered at the origin, matching the shapes' `center` convention):
#   hat{eps}(G) ~= (1/N^2) * sum_{m,n} eps(x_m, y_n) * exp(-2j*pi*(kx*x_m + ky*y_n))
# This is an independent, from-scratch numerical evaluation (direct
# rasterize-and-sum, not calling into `fourier_factorization.py` at all),
# used only to cross-check the analytic Toeplitz entries.

_LX, _LY = 1.3, 0.9
_N_GRID = 900


def _rasterized_coefficient(pattern: Pattern, lattice: Lattice, g1: int, g2: int, inverse: bool) -> complex:
    x = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LX
    y = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")

    def value(material: Material) -> complex:
        eps = complex(material.epsilon_tensor(1.0)[0, 0])
        return 1.0 / eps if inverse else eps

    eps_grid = np.full(X.shape, value(pattern.background), dtype=complex)
    for shape in pattern.shapes:
        if isinstance(shape, Circle):
            dx, dy = X - shape.center[0], Y - shape.center[1]
            mask = dx * dx + dy * dy <= shape.radius**2
        elif isinstance(shape, Rectangle):
            assert shape.angle == 0.0, "rasterized reference only supports axis-aligned rectangles"
            hx, hy = shape.halfwidth
            mask = (np.abs(X - shape.center[0]) <= hx) & (np.abs(Y - shape.center[1]) <= hy)
        else:
            raise NotImplementedError(type(shape))
        eps_grid[mask] = value(shape.material)

    Lk = lattice.reciprocal_vectors()
    k = g1 * Lk[0] + g2 * Lk[1]
    phase = np.exp(-2j * np.pi * (k[0] * X + k[1] * Y))
    return complex(np.sum(eps_grid * phase) / (_N_GRID * _N_GRID))


def test_circle_pattern_epsilon_hat_matches_rasterized_reference():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 1.0)
    core = Material("core", 4.0)
    pattern = Pattern(background=background)
    pattern.add(Circle(center=(0.05, -0.03), radius=0.3, material=core))

    for g1, g2 in [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 2)]:
        for inverse in (False, True):
            analytic = pattern_epsilon_hat(pattern, lattice, g1, g2, wavelength=1.0, inverse=inverse)
            reference = _rasterized_coefficient(pattern, lattice, g1, g2, inverse=inverse)
            assert np.isclose(analytic, reference, atol=5e-3), (g1, g2, inverse, analytic, reference)


def test_rectangle_pattern_epsilon_hat_matches_rasterized_reference():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 2.25)
    core = Material("core", 9.0)
    pattern = Pattern(background=background)
    pattern.add(Rectangle(center=(-0.02, 0.04), halfwidth=(0.25, 0.18), material=core))

    for g1, g2 in [(0, 0), (1, 0), (0, 1), (1, 1), (2, -1)]:
        for inverse in (False, True):
            analytic = pattern_epsilon_hat(pattern, lattice, g1, g2, wavelength=1.0, inverse=inverse)
            reference = _rasterized_coefficient(pattern, lattice, g1, g2, inverse=inverse)
            assert np.isclose(analytic, reference, atol=5e-3), (g1, g2, inverse, analytic, reference)


def test_toeplitz_matrix_diagonal_equals_dc_coefficient():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 1.0)
    core = Material("core", 4.0)
    pattern = Pattern(background=background)
    pattern.add(Circle(center=(0.0, 0.0), radius=0.25, material=core))

    g_vectors = truncate_fourier_orders(lattice, num_orders=9)
    matrix = toeplitz_matrix(pattern, lattice, g_vectors, wavelength=1.0)
    dc = pattern_epsilon_hat(pattern, lattice, 0, 0, wavelength=1.0)
    assert np.allclose(np.diag(matrix), dc)


def test_toeplitz_matrix_entry_matches_g_vector_difference():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 1.0)
    core = Material("core", 4.0)
    pattern = Pattern(background=background)
    pattern.add(Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.15), material=core))

    g_vectors = truncate_fourier_orders(lattice, num_orders=9)
    matrix = toeplitz_matrix(pattern, lattice, g_vectors, wavelength=1.0, inverse=True)
    i, j = 3, 5
    dg1, dg2 = g_vectors[i] - g_vectors[j]
    expected = pattern_epsilon_hat(pattern, lattice, int(dg1), int(dg2), wavelength=1.0, inverse=True)
    assert np.isclose(matrix[i, j], expected)


def test_anisotropic_material_raises_not_implemented():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 1.0)
    aniso = Material("aniso", np.diag([2.0, 3.0, 4.0]).astype(complex))
    pattern = Pattern(background=background)
    pattern.add(Circle(center=(0.0, 0.0), radius=0.2, material=aniso))

    try:
        pattern_epsilon_hat(pattern, lattice, 1, 0, wavelength=1.0)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("expected NotImplementedError for anisotropic material in Phase 2")


# --- Cross-check against the vendored RCWA.jl / RCWA-python FFT convmat ----
#
# `rules.md`'s Testing Requirements explicitly accept an "S4/EMpy/RCWA.jl
# cross-check" as an oracle. The rasterized-reference tests above use a
# from-scratch direct Fourier sum; this test instead reproduces the actual
# algorithm used by the two vendored RCWA repos to build a Toeplitz
# permittivity ("convolution") matrix — raster the pattern on a grid, take
# its FFT, then read matrix entry `[i, j]` off the FFT array at bin
# `G_i - G_j` — see `REFERENCE/Rigorous-Coupled-Wave-Analysis/convolution_matrices/convmat2D.py`
# (`convmat2D`, lines 3-39) and `REFERENCE/RigorousCoupledWaveAnalysis.jl/src/Common/ft2d.jl`
# (`real2recip`, lines 16-27), which build the identical structure via
# `fftshift(fft2(A))/prod(size(A))` indexed at the G-vector difference.
#
# Note: `convmat2D.py` itself carries a comment flagging "NOTE: indexing
# error; N[0] actually corresponds to y and N[1] corresponds to x" — i.e.
# its own authors record an unresolved axis-order caveat in that file. To
# avoid transcribing a known-flagged bug, the FFT-bin indexing below is
# independently (re-)derived and numerically verified rather than copied
# line-for-line. `numpy.fft.fft2` assumes the sampled domain starts at
# index 0 (`x_m = m/N * Lx`); an *uncentered* grid built that way was
# tried first, but it silently truncates any shape whose footprint crosses
# the `x=0`/`y=0` domain edge (confirmed: it produced a DC term that didn't
# match the pattern's true area-weighted-average permittivity). Instead,
# `eps_grid` is rasterized on the same *centered* grid used by
# `_rasterized_coefficient` above (correct for any shape placement within
# the cell), then reordered with `numpy.fft.ifftshift` before the FFT —
# this cyclically re-indexes the centered samples into the 0-start
# convention `fft2` expects, with no truncation, and was verified against
# the DC term and several low-order entries before trusting it for the
# full matrix.

_FFT_N = 512


def _fft_convmat_reference(pattern: Pattern, lattice: Lattice, g_vectors: np.ndarray, inverse: bool) -> np.ndarray:
    x = (np.arange(_FFT_N) / _FFT_N - 0.5) * _LX
    y = (np.arange(_FFT_N) / _FFT_N - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")

    def value(material: Material) -> complex:
        eps = complex(material.epsilon_tensor(1.0)[0, 0])
        return 1.0 / eps if inverse else eps

    eps_grid = np.full(X.shape, value(pattern.background), dtype=complex)
    for shape in pattern.shapes:
        if isinstance(shape, Circle):
            dx, dy = X - shape.center[0], Y - shape.center[1]
            mask = dx * dx + dy * dy <= shape.radius**2
        elif isinstance(shape, Rectangle):
            assert shape.angle == 0.0
            hx, hy = shape.halfwidth
            mask = (np.abs(X - shape.center[0]) <= hx) & (np.abs(Y - shape.center[1]) <= hy)
        else:
            raise NotImplementedError(type(shape))
        eps_grid[mask] = value(shape.material)

    Af = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(eps_grid))) / (_FFT_N * _FFT_N)
    center = _FFT_N // 2

    n = len(g_vectors)
    matrix = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            dg1 = int(g_vectors[i, 0] - g_vectors[j, 0])
            dg2 = int(g_vectors[i, 1] - g_vectors[j, 1])
            matrix[i, j] = Af[center + dg1, center + dg2]
    return matrix


def test_toeplitz_matrix_matches_rcwa_fft_convmat_reference_circle():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 1.0)
    core = Material("core", 4.0)
    pattern = Pattern(background=background)
    pattern.add(Circle(center=(0.05, -0.03), radius=0.3, material=core))

    g_vectors = truncate_fourier_orders(lattice, num_orders=9)
    for inverse in (False, True):
        analytic = toeplitz_matrix(pattern, lattice, g_vectors, wavelength=1.0, inverse=inverse)
        reference = _fft_convmat_reference(pattern, lattice, g_vectors, inverse=inverse)
        assert np.allclose(analytic, reference, atol=5e-3), (inverse, np.max(np.abs(analytic - reference)))


def test_toeplitz_matrix_matches_rcwa_fft_convmat_reference_rectangle():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    background = Material("bg", 2.25)
    core = Material("core", 9.0)
    pattern = Pattern(background=background)
    pattern.add(Rectangle(center=(-0.02, 0.04), halfwidth=(0.25, 0.18), material=core))

    g_vectors = truncate_fourier_orders(lattice, num_orders=9)
    for inverse in (False, True):
        analytic = toeplitz_matrix(pattern, lattice, g_vectors, wavelength=1.0, inverse=inverse)
        reference = _fft_convmat_reference(pattern, lattice, g_vectors, inverse=inverse)
        assert np.allclose(analytic, reference, atol=8e-3), (inverse, np.max(np.abs(analytic - reference)))
