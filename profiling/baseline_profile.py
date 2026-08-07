"""Category 12 target 12.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): baseline
profiler. See `profiling/README.md` for what this is (and isn't) used for.

Measures, on three fixed fixtures (thin-film, 1D grating, 2D pillar), each
at a couple of `num_orders` values:

- **Eigensolve**: `eigenmodes.solve_layer_eigenmodes_*` for the patterned
  layer alone (isolated from Toeplitz construction and S-matrix work).
- **Matrix-solve**: `smatrix._solve` (the `scipy.linalg.lu_factor`/
  `lu_solve` house-convention helper) on a representative `(2n, 2n)`
  interface matrix.
- **S-matrix cascade**: `Simulation.solve()` end to end, for comparison
  against the two isolated stages above.

Uses `time.perf_counter` (wall clock, several repeats, minimum reported --
the standard way to reduce OS-scheduling noise) and `tracemalloc` (peak
Python-level allocation, not a substitute for a native profiler but
sufficient to see relative memory scaling with `num_orders`).

Run with:  python profiling/baseline_profile.py
"""

from __future__ import annotations

import time
import tracemalloc

import numpy as np

from sougata_solver.eigenmodes import solve_layer_eigenmodes_patterned
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.smatrix import _solve

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
REPEATS = 5


def _time_it(fn, *args, repeats: int = REPEATS) -> float:
    """Minimum wall-clock time over `repeats` calls (reduces OS-scheduling
    noise; taking the min, not the mean, is standard microbenchmark
    practice since noise only ever adds time, never removes it)."""
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn(*args)
        best = min(best, time.perf_counter() - t0)
    return best


def _peak_memory_kb(fn, *args) -> float:
    tracemalloc.start()
    fn(*args)
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1024.0


def profile_2d_pillar_eigensolve(num_orders: int) -> None:
    period = 0.7e-6
    lattice = Lattice((period, 0.0), (0.0, period))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(period / 2, period / 2), radius=0.2 * period, material=SI)])
    g = truncate_fourier_orders(lattice, num_orders, "circular")
    n = len(g)
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, 0.6e-6, inverse=False)
    omega = 2 * np.pi / 0.6e-6
    kx = np.zeros(n)
    ky = np.zeros(n)

    t = _time_it(solve_layer_eigenmodes_patterned, omega, kx, ky, epsilon_hat)
    mem = _peak_memory_kb(solve_layer_eigenmodes_patterned, omega, kx, ky, epsilon_hat)
    print(f"  eigensolve (n={n:4d}):        {t * 1e3:10.3f} ms   peak {mem:10.1f} KB")


def profile_matrix_solve(n2: int) -> None:
    rng = np.random.default_rng(0)
    a = rng.normal(size=(n2, n2)) + 1j * rng.normal(size=(n2, n2))
    b = rng.normal(size=(n2, n2)) + 1j * rng.normal(size=(n2, n2))

    t = _time_it(_solve, a, b)
    mem = _peak_memory_kb(_solve, a, b)
    print(f"  matrix-solve (n2={n2:4d}):     {t * 1e3:10.3f} ms   peak {mem:10.1f} KB")


def profile_full_solve(sim: Simulation, excitation: PlaneWaveExcitation, label: str) -> None:
    t = _time_it(sim.solve, excitation)
    mem = _peak_memory_kb(sim.solve, excitation)
    print(f"  Simulation.solve() [{label}]: {t * 1e3:10.3f} ms   peak {mem:10.1f} KB")


def main() -> None:
    print("=== 2D pillar: isolated eigensolve + matrix-solve vs. num_orders ===")
    for num_orders in (9, 25, 49, 81):
        profile_2d_pillar_eigensolve(num_orders)
        profile_matrix_solve(2 * num_orders)

    print("\n=== End-to-end Simulation.solve() on three fixed fixtures ===")

    # Thin film: uniform, no patterning.
    sio2 = Material("sio2", 1.46**2)
    lattice_uniform = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    sim_thin_film = Simulation(lattice_uniform, [Layer("sio2", 0.1e-6, material=sio2)], num_orders=1, incidence=AIR, transmission=SI)
    exc_thin_film = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    profile_full_solve(sim_thin_film, exc_thin_film, "thin-film")

    # 1D grating.
    lattice_1d = Lattice1D(0.7e-6)
    pattern_1d = Pattern(background=AIR)
    pattern_1d.add(Slab(center_x=0.0, halfwidth=0.15e-6, material=SI))
    sim_1d = Simulation(lattice_1d, [Layer("grating", 0.3e-6, pattern=pattern_1d)], num_orders=9, incidence=AIR, transmission=AIR)
    exc_1d = PlaneWaveExcitation(1.0e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    profile_full_solve(sim_1d, exc_1d, "1D grating, num_orders=9")

    # 2D pillar, two num_orders points.
    lattice_2d = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    pattern_2d = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=SI)])
    for num_orders in (9, 49):
        sim_2d = Simulation(lattice_2d, [Layer("pillar", 0.3e-6, pattern=pattern_2d)], num_orders=num_orders, incidence=AIR, transmission=AIR)
        exc_2d = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
        profile_full_solve(sim_2d, exc_2d, f"2D pillar, num_orders={num_orders}")


if __name__ == "__main__":
    main()
