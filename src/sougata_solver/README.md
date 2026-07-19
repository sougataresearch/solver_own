# `src/sougata_solver/` — Solver Core

The library itself. Everything here is imported, never run directly (run
scripts live in [`structures/`](../structures/) and
[`postprocessing/`](../postprocessing/) instead). See the repo-root
[`README.md`](../README.md) for the project overview and
[`architecture.md`](../architecture.md) for the full data-flow diagram.

## Module map

Roughly in the order data flows through them for a single `Simulation.solve()` call:

| Module | Responsibility |
|---|---|
| [`materials.py`](materials.py) | `Material`: scalar/dispersive (`from_nk`, CSV) and tensor permittivity |
| [`geometry.py`](geometry.py) | `Lattice`, `Shape` (`Circle`, `Rectangle`), `Pattern` — in-plane geometry and analytic Fourier transforms |
| [`fourier_basis.py`](fourier_basis.py) | Circular G-vector truncation (`truncate_fourier_orders`) — which Fourier orders are kept |
| [`fourier_factorization.py`](fourier_factorization.py) | `pattern_epsilon_hat`, `toeplitz_matrix` — builds the Toeplitz permittivity matrix for a patterned layer (Phase 2; not yet wired into `simulation.py`) |
| [`layer.py`](layer.py) | `Layer`, `LayerStack`, `LayerEigenmodes` — the data model for a stack, including the semi-infinite incidence/transmission half-spaces |
| [`eigenmodes.py`](eigenmodes.py) | Per-layer eigenmode solve: `q` (propagation constants), `phi` (eigenvectors), `kp` (k-parallel operator). Uniform layers only today — patterned-layer solve is Phase 4 |
| [`excitation.py`](excitation.py) | `PlaneWaveExcitation`: angle/polarization decomposition into s/p, and inversion to the incident mode-amplitude vector |
| [`smatrix.py`](smatrix.py) | Interface + propagation S-matrices, Redheffer star-product cascading (`SMatrixStack`) — dimension-agnostic, needs no changes for Phase 3/4 |
| [`fields.py`](fields.py) | `z_poynting_flux` — reflected/transmitted power from mode amplitudes |
| [`polarimetry.py`](polarimetry.py) | Jones/Mueller matrix construction from s/p amplitudes (`decompose_sp` is reused by `postprocessing/`) |
| [`simulation.py`](simulation.py) | Top-level orchestration: `Simulation.solve()` wires the above into a `SimulationResult` (`.reflectance()`, `.transmittance()`) |
| [`output_paths.py`](output_paths.py) | `outputs/YYYY_MM_DD/HH_MM_SS_<run>/` folder + `run_metadata.txt` helper, used by `structures/` scripts |

## Design rules specific to this folder

- **No formula here is original.** Every non-trivial equation is transcribed
  from a named, line-numbered source (usually `S4/S4/rcwa.cpp` or
  `S4/S4/pattern/pattern.c`) and cited in the module or function docstring —
  see [`rules.md`](../rules.md)'s AI Coding Rules before adding anything new.
- **Sign/phase convention is physics-style** (`exp(+jkz)` spatial phase,
  `d/dt -> -jw`), inherited from S4 — see `smatrix.py:108`'s propagation
  phase and `progress_log.md`'s 2026-07-19 entry for the reasoning. Hand
  derivations from a textbook using the opposite (engineering) convention
  must be sign-flipped (`j -> -j`) before porting into this module.
- **No mesh, no PML.** This is a Fourier-modal method, not FEM/FDTD — layers
  are solved analytically, and the incidence/transmission half-spaces are
  literal `thickness = math.inf` (`layer.py`), not a truncated domain. There
  is nothing analogous to an absorbing boundary condition to add here.
- **`smatrix.py` and `eigenmodes.py`'s `build_kp_matrix`** are already
  written to accept the general (patterned) case (`kp_matrix` takes a full
  `epsilon_inv` matrix, not just a scalar) — Phase 3/4 should target that
  existing interface rather than introducing a new one.
- Anisotropic materials and patterned-layer eigenmode solving currently
  raise `NotImplementedError` (`materials.py`, `simulation.py:97-101`) —
  these are Phase 4/6, not bugs.
