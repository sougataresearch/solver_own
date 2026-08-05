"""Category 13 target 13.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
repeatable runtime/memory benchmarks for thin-film, trench, pillar, and
tapered structures. See `profiling/README.md` for what this is (and
isn't) used for -- extends `profiling/baseline_profile.py` (Category 12
target 12.1, which covered thin-film/1D/2D at a single depth) with the
one case it didn't: a tapered (Phase 5 staircase-discretized) structure.

"Repeatable" means: fixed geometry/material/excitation inputs every run
(no randomness), so two runs on the same machine are directly comparable,
and a run on a different machine is comparable in relative (not absolute)
terms -- the standard microbenchmark caveat, restated from
`profiling/README.md`.

Run with:  python profiling/benchmark_suite.py
"""

from __future__ import annotations

import time
import tracemalloc

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_circle_layers

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
REPEATS = 5


def _time_it(fn, *args, repeats: int = REPEATS) -> float:
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


def _report(label: str, sim: Simulation, excitation: PlaneWaveExcitation) -> None:
    t = _time_it(sim.solve, excitation)
    mem = _peak_memory_kb(sim.solve, excitation)
    print(f"  {label:<40} {t * 1e3:10.3f} ms   peak {mem:10.1f} KB")


def benchmark_thin_film() -> None:
    print("=== Thin film (uniform multilayer, no patterning) ===")
    sio2 = Material("sio2", 1.46**2)
    lattice = Lattice((1.0e-6, 0.0), (0.0, 1.0e-6))
    sim = Simulation(lattice, [Layer("sio2", 0.1e-6, material=sio2)], num_orders=1, incidence=AIR, transmission=SI)
    excitation = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    _report("thin-film", sim, excitation)


def benchmark_trench() -> None:
    print("=== Trench (1D lamellar grating) ===")
    lattice = Lattice1D(0.7e-6)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.15e-6, material=SI))
    for num_orders in (9, 25):
        sim = Simulation(lattice, [Layer("grating", 0.3e-6, pattern=pattern)], num_orders=num_orders, incidence=AIR, transmission=AIR)
        excitation = PlaneWaveExcitation(1.0e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
        _report(f"trench, num_orders={num_orders}", sim, excitation)


def benchmark_pillar() -> None:
    print("=== Pillar (2D patterned via) ===")
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=SI)])
    for num_orders in (9, 49):
        sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=num_orders, incidence=AIR, transmission=AIR)
        excitation = PlaneWaveExcitation(0.6e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
        _report(f"pillar, num_orders={num_orders}", sim, excitation)


def benchmark_tapered() -> None:
    print("=== Tapered via (Phase 5 staircase discretization) ===")
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    excitation = PlaneWaveExcitation(1.0e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    for num_slices in (4, 16):
        layers = staircase_circle_layers(
            center=(0.35e-6, 0.35e-6), top_radius=0.24e-6, bottom_radius=0.10e-6,
            thickness=0.46e-6, num_slices=num_slices, shape_material=SI, background_material=AIR,
        )
        sim = Simulation(lattice, layers, num_orders=9, incidence=AIR, transmission=AIR)
        _report(f"tapered via, num_slices={num_slices}", sim, excitation)


def main() -> None:
    benchmark_thin_film()
    benchmark_trench()
    benchmark_pillar()
    benchmark_tapered()


if __name__ == "__main__":
    main()
