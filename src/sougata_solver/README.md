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
| [`materials.py`](materials.py) | `Material`: scalar/dispersive (`from_nk`, CSV) and tensor permittivity, with construction-/call-time validation (Category 5 target 5.1); analytic dispersion models `from_sellmeier`/`from_cauchy`/`from_lorentz`/`from_drude`/`from_drude_lorentz` (targets 5.2-5.6, the last with published `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM` metal presets); optional `source` citation metadata (target 5.8) |
| [`geometry.py`](geometry.py) | `Lattice`, `Shape` (`Circle`, `Rectangle`, `Ellipse`, `Polygon`, `Slab`), `Pattern` — in-plane geometry and analytic Fourier transforms; construction-time shape/lattice validation and `validate_pattern_fits_lattice` (unit-cell self-overlap policy, Category 4 targets 4.1/4.2) |
| [`geometry_io.py`](geometry_io.py) | Minimal, safe JSON `Pattern`-import format (`pattern_from_dict`/`pattern_from_json_file`) — `json` module only, never `eval`/`exec` (Category 4 target 4.6, parser-only, not wired into `Simulation`/`Layer`) |
| [`fourier_basis.py`](fourier_basis.py) | Circular G-vector truncation (`truncate_fourier_orders`) — which Fourier orders are kept |
| [`fourier_factorization.py`](fourier_factorization.py) | `pattern_epsilon_hat`/`toeplitz_matrix` (scalar isotropic) and `pattern_epsilon_hat_component`/`toeplitz_matrix_component` (per-tensor-component) — builds the Toeplitz permittivity matrix for a patterned layer |
| [`layer.py`](layer.py) | `Layer`, `LayerStack`, `LayerEigenmodes`, `EigenmodeDiagnostics` — the data model for a stack, including the semi-infinite incidence/transmission half-spaces, construction-time thickness validation (Category 7 target 7.1), and per-solve eigenvalue/conditioning diagnostics (Category 2 target 2.2) |
| [`staircase.py`](staircase.py) | `slice_profile`, a general geometry-to-layer-slices interface (Category 4 target 4.7), plus the tapered-sidewall generators built on it (`staircase_circle_layers`, `staircase_rectangle_layers`, `staircase_slab_layers`) |
| [`eigenmodes.py`](eigenmodes.py) | Per-layer eigenmode solve: `q` (propagation constants), `phi` (eigenvectors), `kp` (k-parallel operator). Covers uniform isotropic, uniform diagonal/in-plane-anisotropic, 1D-patterned, and 2D-patterned (isotropic and diagonal/in-plane-anisotropic) layers, plus `classify_propagating` (mode classification), `_canonical_mode_order` (deterministic degeneracy ordering), `ILL_CONDITIONED_THRESHOLD`/`DEGENERATE_GAP_THRESHOLD` (conditioning/near-degeneracy `WARNING`s, Category 2 targets 2.2/2.4), `_dense_inverse` (Category 12 target 12.2, house-convention-consistent full-matrix inverse) and `svd_diagnostics` (target 12.4, opt-in singular-value spectrum) |
| [`excitation.py`](excitation.py) | `PlaneWaveExcitation`: angle/polarization decomposition into s/p, and inversion to the incident mode-amplitude vector |
| [`smatrix.py`](smatrix.py) | Interface + propagation S-matrices, Redheffer star-product cascading (`SMatrixStack`) — dimension-agnostic, needs no changes for Phase 3/4; `interior_amplitudes` (Category 9 target 9.3) recovers forward/backward mode amplitudes at an interior interface, independently derived (`decisions.md` ADR-015) |
| [`fields.py`](fields.py) | `z_poynting_flux`/`tangential_e_field` — reflected/transmitted power and tangential E-field from mode amplitudes; full real-space `(Ex,Ey,Ez,Hx,Hy,Hz)` reconstruction at any point/depth (`modal_field_components`/`propagate_amplitudes`/`reconstruct_field_at_points`) and NumPy field-grid export (`save_field_grid_npz`), Category 9 targets 9.1-9.8 |
| [`polarimetry.py`](polarimetry.py) | Jones/Mueller matrix construction from s/p amplitudes (`decompose_sp` is reused by `postprocessing/`) |
| [`sweep.py`](sweep.py) | Category 8 targets 8.1-8.8: `SweepResult` (typed one-parameter-sweep container), `sweep_wavelength`/`sweep_theta`/`sweep_phi`/`sweep_polarization`/`sweep_thickness` (thin wrappers repeating `Simulation.solve()`), `harmonic_study`/`find_convergence_index`/`auto_select_num_orders` (harmonic-order convergence study and a conservative, `decisions.md` ADR-018, stopping criterion) |
| [`simulation.py`](simulation.py) | Top-level orchestration: `Simulation.solve()` wires the above into a `SimulationResult` (`.reflectance()`, `.transmittance()`, `.diffraction_efficiencies()`, `.order_classification()`, `.layer_absorption()` — Category 7 targets 7.5/7.6; `.complex_amplitudes()`/`.diffraction_angles()`/`.energy_balance()` — Category 10 targets 10.1-10.3); validates every patterned layer against `geometry.validate_pattern_fits_lattice` at construction; caches per-pattern Toeplitz matrices (Category 7 targets 7.3/7.4, `decisions.md` ADR-016) and per-layer eigenmode solves (Category 13 target 13.3, `decisions.md` ADR-022) |
| [`ocd.py`](ocd.py) | Category 11 targets 11.1/11.2/11.4: `OCDTrapezoidParams` (validated CD-first OCD parameters), `trapezoid_trench_layers` (thin wrapper around `staircase.staircase_slab_layers`), `rounded_rectangle_polygon` (arc-sampled corner rounding built on `geometry.Polygon`) — no new physics anywhere in this module |
| [`vectorized.py`](vectorized.py) | Category 13 target 13.4: `sweep_wavelength_vectorized` — a narrowly-scoped batched wavelength sweep for uniform-isotropic-only (thin-film) stacks; every batched function a formula-identical re-expression of an already-cited scalar function, no general vectorized backend (`decisions.md` ADR-023) |
| [`output_paths.py`](output_paths.py) | `outputs/YYYY_MM_DD/HH_MM_SS_<run>/` folder + `run_metadata.txt` helper (explicit UTF-8, Category 5 target 5.8), used by `structures/` scripts |

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
- **`smatrix.py` and `eigenmodes.py`'s `build_kp_matrix`** were written to
  accept the general (patterned) case (`kp_matrix` takes a full
  `epsilon_inv` matrix, not just a scalar) from the start — every later
  solver (1D, 2D patterned, anisotropic) reuses this same interface rather
  than introducing a new one.
- Longitudinally-coupled anisotropic materials (`eps_xz/eps_yz/eps_zx/eps_zy`)
  still raise `NotImplementedError` in `simulation.py`, naming
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` target 1.5 — evaluated and explicitly
  deferred (no citable + independently-benchmarkable formulation found),
  not a bug.
- **Every `ValueError`/`NotImplementedError`/`LinAlgError`/`WARNING`
  condition this package raises or logs is documented in `design.md`'s
  Failure Contract** (Category 2 target 2.1) and backed by
  `tests/test_failure_contract.py` — check that table before adding a new
  failure mode, and keep it in sync when you do.
- `geometry.Polygon`'s analytic (not raster/FFT) Fourier transform is a
  narrow, explicit revisit of `decisions.md` ADR-005 (`decisions.md`
  ADR-013) — GDS/arbitrary-mask import remains out of scope; do not treat
  `Polygon`'s existence as license to add raster-based geometry elsewhere.
- `SimulationResult.layer_absorption()` (Category 7 targets 7.5/7.6) gives
  per-layer absorbed power, validated against the `R+T+sum(A)=1`
  energy-balance identity. It inherits `interior_amplitudes`/
  `propagate_amplitudes`'s numerical-stability envelope: a thick, highly
  lossy, high-`num_orders` layer can numerically overflow the deepest
  evanescent modes' backward-propagated amplitude — see
  `troubleshooting.md` and `decisions.md` ADR-017 before trusting a result
  from a parameter regime that extreme.
