# Detailed Design — sougata_solver

`sougata_solver` has no database and no UI in the traditional sense, so those
template sections are replaced below with what actually exists: the public
Python API (in place of "API Design") and the example-script/plotting
surface (in place of "UI/UX Design") — see the "N/A sections" note at the
end.

## Algorithms

### 1. Per-layer eigenmode solve (uniform isotropic layer — done)

Source of truth: `S4/S4/rcwa.cpp::SolveLayerEigensystem_uniform`
(lines 422-502), transcribed in `eigenmodes.py::solve_layer_eigenmodes_uniform`.

For a homogeneous isotropic layer, the eigenbasis coincides with the
plane-wave basis (`phi = I`), so each diffraction order's z-propagation
constant is closed-form:

```
q[i] = branch_select(eps * omega^2 - kx[i]^2 - ky[i]^2)
```

`branch_select` (`eigenmodes.py::_select_q_branch`) chooses the
outgoing/decaying root: for real-valued `q^2`, positive values give a real
(propagating) `q`, negative values give a purely-imaginary (evanescent,
decaying-forward) `q`; for complex `q^2` (absorbing media), the principal
square root is flipped in sign if needed so `Im(q) >= 0`.

### 2. General (non-uniform) eigenmode solve (patterned layers — Phase 4a/4b)

Source of truth (not yet transcribed): `S4/S4/rcwa.cpp::SolveLayerEigensystem`,
lines 794-827. This is the general eigenproblem for a layer whose in-plane
permittivity varies periodically:

```
op = Epsilon2 @ kp - (coupling terms from anisotropy, if present)
(q^2, phi) = eig(op)
q = branch_select(q^2)          # reuse eigenmodes.py::_select_q_branch
```

where `Epsilon2` is the "direct" Toeplitz permittivity matrix and `kp` is
built by `eigenmodes.build_kp_matrix` (already implemented — it already
accepts a full `(n,n)` `epsilon_inv` matrix, not just a scalar, for exactly
this case). The two new pieces needed (Phase 2/4) are: (a) constructing
`Epsilon2`/`epsilon_inv` as Toeplitz matrices from a `Pattern`'s shapes, and
(b) the general complex eigendecomposition + degenerate-eigenvalue handling.
**This is the highest-risk remaining algorithm in the project** — general
eigendecompositions can have near-degenerate eigenvalues with
poorly-conditioned eigenvectors; the mitigation is split across two phases
(`phases.md`): Phase 4a validates the solver on moderate-contrast cases
first (cross-check against S4 itself, which has already solved this
problem correctly, rather than a from-scratch re-derivation of stability
fixes), and Phase 4b is a dedicated follow-up phase that deliberately
stress-tests near-degenerate/high-contrast cases rather than trusting
Phase 4a's easier cases to have exercised that regime.

### 3. Fourier factorization (Toeplitz permittivity construction — Phase 2)

For a patterned layer with shapes `S_1, ..., S_k` over a `background`, the
Fourier coefficient of `eps(x,y)` at reciprocal vector `G` is:

```
eps_hat(G) = [background_term(G) + sum_i shape_i_contribution(G)] / unit_cell_area
```

where each shape's contribution uses its already-implemented
`fourier_transform(kx, ky)` (`geometry.py::Circle`/`Rectangle`), and the
`Pattern.containment_tree()` (already implemented) determines which shapes
are nested inside others so overlapping-area double-counting is corrected
via the S4 subtraction-rule convention (already documented in
`geometry.py`'s `Pattern` docstring).

The **Toeplitz matrix** used in the eigenproblem is
`M[i,j] = eps_hat(G_i - G_j)` for the truncated set of `G` indices selected
by `fourier_basis.truncate_fourier_orders`. Two such matrices are needed:
`epsilon_hat` (direct) and `epsilon_inv_hat` (Toeplitz of `1/eps(x,y)`,
*not* the matrix inverse of `epsilon_hat` — this distinction, the
"Fourier factorization rule" / Li's inverse rule, is the single most
common source of wrong-but-plausible-looking RCWA results industry-wide,
and is exactly why `eigenmodes.build_kp_matrix` already takes a distinct
`epsilon_inv` argument rather than computing `inv(epsilon_hat)` internally).

### 3a. Fourier-factorization rule inventory (Category 3 target 3.1)

`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 3 target 3.1: which
direct-rule/inverse-rule/numerical-matrix-inverse choice is actually made
by every uniform, 1D, and 2D solver branch, with citations — built by
reading each solver's already-cited construction in `eigenmodes.py`, not
reconstructed from memory. "Direct rule" means a Toeplitz matrix of
`hat{eps}(G)` (Laurent's rule); "inverse rule" means a Toeplitz matrix of
`hat{1/eps}(G)` (a *separate* Fourier factorization of `1/eps(x,y)`, not
derived from the direct-rule matrix); "numerical inverse" means
`inv(direct-rule Toeplitz)`, a plain linear-algebra inverse of an
already-built direct-rule matrix — a third, distinct option this project
uses in two places below and which is easy to conflate with "inverse rule"
despite being a different Fourier-factorization choice entirely (see
`design.md`'s Algorithm 3 above and Li 1996 in `references.md`).

| Solver branch | `kp`'s `epsilon_inv` argument | `Epsilon2` | Rule | Citation |
|---|---|---|---|---|
| `solve_layer_eigenmodes_uniform` (scalar isotropic) | `1/eps` (scalar) | `eps * I` | exact (no Fourier factorization needed — spatially uniform) | `rcwa.cpp::SolveLayerEigensystem_uniform`, 422-502 |
| `solve_layer_eigenmodes_uniform_diagonal`/`_inplane` (uniform tensor) | `1/eps_zz` (scalar) | tensor components `* I_n` | exact (uniform layer, no Fourier factorization) | `S4.cpp:1889-1906` |
| `solve_layer_eigenmodes_1d`, TE-like block | `epsilon_inv_hat` (inverse-rule Toeplitz) | `epsilon_hat` (direct-rule) | direct rule | `fmm_closed.cpp:110-132`, "1D proper FFF rule" branch |
| `solve_layer_eigenmodes_1d`, TM-like block | `epsilon_inv_hat` (inverse-rule Toeplitz) | `inv(epsilon_inv_hat)` (**numerical inverse of the inverse-rule Toeplitz**) | Li's (1996) inverse rule | same citation — this is the one block in the whole project that actually uses the inverse-rule Toeplitz, and even then only after a further numerical inversion, not directly |
| `solve_layer_eigenmodes_patterned` (2D isotropic) | `inv(epsilon_hat)` (**numerical inverse of the direct-rule Toeplitz**) | `block_diag(epsilon_hat, epsilon_hat)` (direct rule, both blocks) | ordinary Laurent's rule throughout — no Li correction | `fmm_closed.cpp:133-139,162-163`, the true-2D `!use_polarization_basis` branch |
| `solve_layer_eigenmodes_patterned_inplane` (2D anisotropic) | `inv(epsilon_hat_zz)` (**numerical inverse of the direct-rule `eps_zz` Toeplitz**) | full tensor-component direct-rule Toeplitz blocks | ordinary Laurent's rule (tensor generalization of the row above) | `fmm_closed.cpp:165-256`, `have_tensor` branch |

Two findings from actually re-checking this end to end rather than
resummarizing existing docstrings from memory (per `rules.md` Documentation
Standards' "verify, don't paraphrase from memory"):

1. **Only one code path in this project uses `epsilon_inv_hat` (the
   separately-Fourier-factorized inverse-rule Toeplitz) at all**: the 1D
   TM-like block, and even there only through a further numerical
   `inv(...)`, not as `Epsilon2`/`kp` input directly. Every 2D path (both
   isotropic and anisotropic) uses `inv(direct-rule Toeplitz)` instead —
   `epsilon_inv_hat`/`toeplitz_matrix(..., inverse=True)` is genuinely
   1D-only infrastructure, exactly as already stated in `memory.md`'s
   Phase 4a entry, now confirmed as a table entry rather than prose buried
   in a phase-completion note.
2. **"Numerical inverse of a direct-rule Toeplitz" is not the same
   operation as "inverse-rule Toeplitz," even though both blocks are
   labeled `Epsilon_inv`/`einv` in code** — the 1D TM block computes
   `inv(epsilon_inv_hat)` (inverse of an *already-inverse-rule* matrix,
   Li's actual correction), while every 2D path computes `inv(epsilon_hat)`
   (inverse of the *direct*-rule matrix, ordinary Laurent's rule, no Li
   correction) — table above makes this distinction explicit rather than
   leaving a reader to infer it from four separate docstrings.

Backed by `tests/test_fourier_factorization_rules.py`, which pins each row
above against actual solver behavior (not just the table's prose) so a
future refactor that silently changes which matrix a solver inverts shows
up as a test failure.

### 3b. Ellipse and Polygon Fourier transforms (Category 4 targets 4.3-4.5)

`Ellipse` (target 4.3) and `Polygon` (targets 4.4-4.5) extend `Shape` the
same way `Circle`/`Rectangle` already do -- a closed-form `fourier_transform(kx,
ky)`, no raster/FFT, no change to `Pattern`/`fourier_factorization.py`. Both
are transcribed from `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`
(lines 889-1032, the same function `Circle`/`Rectangle`'s existing citations
already reference): `ELLIPSE` (lines 955-964) rescales `Circle`'s `jinc`
argument anisotropically by the semi-axis ratio; `POLYGON` (lines 974-1008)
is a closed-form boundary/edge-sum formula, not a raster or FFT operation --
see `decisions.md` ADR-013 for the full accuracy-contract decision (target
4.4) made before implementing `Polygon` (target 4.5), and for why this does
not revisit or depend on ADR-012's separate FFF/NVM deferral (a different,
harder problem -- correcting the *eigenoperator's* Fourier factorization at
a discontinuous interface, not a single shape's own boundary integral).

### 4. S-matrix cascading (done)

Redheffer star product, transcribed from `S4/S4r/StarProduct.hpp`
(`T2Sblocks` lines 51-65, `StarProduct` lines 83-110) —
`smatrix.py::interface_smatrix`, `propagation_smatrix`, `star_product`.
Chosen over transfer-matrix cascading specifically because evanescent modes
in a transfer matrix grow without bound through a thick/lossy layer,
destroying numerical precision; the S-matrix formulation keeps every
intermediate quantity bounded. See `architecture.md`'s Technology Choices.

### 5. Field/power extraction (done for R/T; Phase 7 for full reconstruction)

`fields.py::z_poynting_flux`, transcribed from
`S4/S4/rcwa.cpp::GetZPoyntingFlux` (lines 1846-1897) — deliberately *not*
re-derived from scratch, per the docstring, because "a from-scratch
re-derivation of the sign/normalization conventions embedded in the
`kp`/`phi` operators risked introducing exactly the kind of subtle error
this module needs to avoid." `tangential_e_field` similarly transcribes
`GetInPlaneFieldVector` (lines 1959-1995) and documents a specific, easy
mistake: `E = phi @ (a+b)` is *not* correct — that combination is actually
`H`; `E` requires `(a-b)` with an index swap and sign flip. Phase 7 extends
this from "tangential field at one interface" to "full E/H(x,y,z) on a
grid," using `SMatrixStack.partial_smatrix_up_to` to get the local mode
amplitudes at an arbitrary depth, then inverse-Fourier-summing over the
retained G-vectors.

## "UI/UX" — Script Surface (structures/ and postprocessing/)

There is no GUI. The user-facing surface is split into two directories by
responsibility (see `README.md`'s Folder Structure and `decisions.md`
ADR-009):

- **`structures/`** — build a lattice/layer stack/materials and run the
  solver; produces raw results (printed R/T, or a raw-field CSV for
  ellipsometry-style scripts). This is what you run first, and what you
  edit to change geometry, dimensions, or materials.
- **`postprocessing/`** — takes a `structures/` script's raw output and
  derives something further from it: Jones/Mueller matrices and
  ellipsometric angles today (`postprocessing/jones_mueller_ellipsometry.py`,
  reading the CSV written by `structures/thin_film/sio2_on_si_ellipsometry_run.py`),
  and — planned, not yet built — RI/thickness extraction (inverse fitting
  against measured data). This directory never calls `Simulation.solve`
  itself; it only reads already-computed data.

Concretely:

1. **Library API**, imported into small, single-purpose scripts in
   `structures/` (e.g. `sio2_on_si_thin_film.py`, `custom_multistack.py` —
   copy the latter for a new stack) and `postprocessing/`.
2. **Console output**: scripts print a table of wavelength/R/T/A to stdout
   during a sweep (see `structures/thin_film/sio2_on_si_thin_film.py::main`).
3. **CSV + metadata output**: `structures/` scripts get one output folder
   per invocation (`output_paths.run_output_dir`, `outputs/YYYY_MM_DD/
   HH_MM_SS_<run_name>/`) and write both their raw CSV and a
   `run_metadata.txt` (`output_paths.write_run_metadata`) into it — the
   metadata file records which script produced the run and its key
   parameters (materials, thicknesses, angle, wavelength range, ...) so a
   run folder is identifiable without re-reading code or guessing from the
   timestamp alone. See ADR-010 in `decisions.md`.
4. **Plotting is always `postprocessing/`, never `structures/`.**
   `postprocessing/plot_thin_film_rt.py` reads a `structures/` script's CSV
   (via `output_paths.find_latest_output` by default, or an explicit path),
   plots it, and saves the PNG back into that *same* run folder — a plot is
   a derived view of already-computed data, not a new run. This also
   applies to Phase 7's planned cross-section field-intensity plots for
   trench/via structures: the `structures/` side of that phase saves raw
   field data, and a `postprocessing/` script reads and plots it, following
   this same split.

## "API Design" — Public Python API

### Public API Inventory (Category 15 target 15.1)

**Found and fixed while compiling this inventory**: `src/sougata_solver/__init__.py`'s
top-level `__all__` re-export list had not kept pace with new geometry
primitives — `Ellipse`/`Polygon` (Category 4 targets 4.3/4.5) and
`Lattice1D`/`Slab` (Phase 3) were public, working classes for many
categories but not reachable via `from sougata_solver import ...`. Fixed
by adding all four; no behavior change, purely an export-surface fix.

**Stable public surface** (safe to depend on; changes here would be a
breaking-change event):

| Symbol | Module | Stability |
|---|---|---|
| `Material` | `materials.py` | stable — construction API (`Material(name, eps)`, `from_nk`, `from_sellmeier`, `from_cauchy`, `from_lorentz`, `from_drude`, `from_drude_lorentz`, `from_permittivity_tensor`) unchanged since introduction |
| `Lattice`, `Lattice1D` | `geometry.py` | stable |
| `Circle`, `Rectangle`, `Ellipse`, `Polygon`, `Slab` | `geometry.py` | stable (shape constructors); `Polygon.signed_distance_normal`/`contains` are lower-level, used internally |
| `Pattern` | `geometry.py` | stable |
| `Layer`, `LayerStack` | `layer.py` | stable |
| `PlaneWaveExcitation` | `excitation.py` | stable, not yet top-level re-exported (import from `sougata_solver.excitation`) |
| `Simulation`, `SimulationResult` | `simulation.py` | stable, not yet top-level re-exported (import from `sougata_solver.simulation`); every `SimulationResult` method (`.reflectance()`, `.transmittance()`, `.diffraction_efficiencies()`, `.order_classification()`, `.layer_absorption()`, `.complex_amplitudes()`, `.diffraction_angles()`, `.energy_balance()`) is public |
| `sweep.SweepResult`, `sweep_wavelength`/`sweep_theta`/`sweep_phi`/`sweep_polarization`/`sweep_thickness`/`harmonic_study`/`find_convergence_index`/`auto_select_num_orders` | `sweep.py` | stable |
| `ocd.OCDTrapezoidParams`, `trapezoid_trench_layers`, `rounded_rectangle_polygon` | `ocd.py` | stable |
| `vectorized.sweep_wavelength_vectorized` | `vectorized.py` | stable, narrow scope (see its own docstring) |
| `staircase.slice_profile`, `staircase_circle_layers`/`staircase_rectangle_layers`/`staircase_slab_layers` | `staircase.py` | stable |
| `geometry_io.pattern_from_dict`/`pattern_from_json_string`/`pattern_from_json_file` | `geometry_io.py` | stable (parser only, not solver-wired) |
| `polarimetry.py`'s Jones/Mueller functions | `polarimetry.py` | stable |
| `output_paths.py`'s run-folder helpers | `output_paths.py` | stable, intended for `structures/`/`postprocessing/` scripts |
| `eigenmodes.svd_diagnostics`/`SVDDiagnostics` | `eigenmodes.py` | stable, opt-in diagnostic |
| `layer.EigenmodeDiagnostics` | `layer.py` | stable, read-only diagnostic data |

**Internal/unstable** (implementation detail, may change without notice —
identifiable by the `_` prefix convention `rules.md`'s Naming Conventions
already establishes): every `_`-prefixed function/method across
`src/sougata_solver/` (e.g. `eigenmodes._select_q_branch`,
`eigenmodes._dense_inverse`, `smatrix._solve`,
`simulation.Simulation._cached_toeplitz`/`_cached_layer_eigenmodes`,
`vectorized._batched_*`, `geometry._require_finite`). Also unstable: the
exact internal dict-key shape of `Simulation._toeplitz_cache`/
`_eigenmode_cache` (an implementation detail of the caching design,
`decisions.md` ADR-016/ADR-022, not part of the public contract even
though the attributes themselves are inspectable).

The intended import surface (already `__all__`-exported from
`src/sougata_solver/__init__.py`):

```python
from sougata_solver import Material, Lattice, Lattice1D, Circle, Rectangle, Ellipse, Polygon, Slab, Pattern, Layer, LayerStack
```

Plus, imported directly from their submodules (not yet re-exported at
top level — worth revisiting once Phase 3/4 land and usage patterns
stabilize):

```python
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.simulation import Simulation, SimulationResult
```

**Typical call sequence** (already the pattern in every `structures/*.py` script):

```python
material = Material(...)                       # or Material.from_nk(...)
lattice = Lattice(a=..., b=...)                # or, once Phase 3 lands, Lattice1D(period)
layers = [Layer(name, thickness, material=...)] # or pattern=Pattern(...) once Phase 3/4 land
sim = Simulation(lattice, layers, num_orders, incidence, transmission)
result = sim.solve(PlaneWaveExcitation(wavelength, theta, phi, s_amplitude, p_amplitude))
result.reflectance(); result.transmittance()
```

This API is intentionally **not** builder-pattern or fluent-interface
styled — plain constructors and dataclasses only, consistent with the "no
framework magic" non-functional requirement in `PRD.md`.

**Bottom (reverse-side) illumination** (Category 6 target 6.6,
`decisions.md` ADR-014): already possible with the existing constructor,
no new parameter — reverse the layer list and swap which material plays
`incidence`/`transmission`:

```python
sim_bottom = Simulation(lattice, list(reversed(layers)), num_orders,
                          incidence=transmission_material, transmission=incidence_material)
```

Verified (not just asserted) correct via the Stokes transmittance-
reciprocity relation (`tests/test_bottom_incidence.py`): at normal
incidence, the reversed simulation's transmittance matches the original's
to `~1e-15`.

**Parameter sweeps** (Category 8 targets 8.1-8.8, `src/sougata_solver/sweep.py`):

```python
from sougata_solver.sweep import (
    SweepResult, sweep_wavelength, sweep_theta, sweep_phi,
    sweep_polarization, sweep_thickness, harmonic_study,
    find_convergence_index, auto_select_num_orders,
)

sweep = sweep_wavelength(sim, wavelengths, theta, phi, s_amplitude, p_amplitude)
sweep.reflectance(); sweep.transmittance()   # np.ndarray, one entry per point
```

Every `sweep_*` function is a thin wrapper calling `Simulation.solve()`
once per parameter value — no new solver-formula risk, per the category's
own exit criterion ("each sweep is equivalent to repeated scalar solves"),
confirmed directly by `tests/test_sweep.py` for every function (build the
same result via the sweep function and via a manual per-point loop,
compare). `sweep_theta`/`sweep_phi` (fixed-wavelength angle sweeps) are
the scenario ADR-016's Toeplitz-matrix cache was actually measured
against — reusing one `Simulation` across such a sweep gets that ~30%
reduction automatically, no extra code needed here.

`harmonic_study` takes a `Simulation`-*builder* callable
(`Callable[[int], Simulation]`), not a single instance, because
`num_orders` is fixed for a `Simulation` instance's entire lifetime (the
same invariant ADR-016's Toeplitz cache relies on) — there is no supported
way to resweep `num_orders` on one live instance. `find_convergence_index`
is a conservative data-selection rule (every *later* point must also stay
within tolerance, not just the immediate next one — see `decisions.md`
ADR-018 for why), validated against thin-film/trench/pillar fixtures
before `auto_select_num_orders` (target 8.8) was implemented, per that
target's own "implement only after 8.7 succeeds" wording.

**Public optical outputs** (Category 10 targets 10.1-10.6):
`SimulationResult.complex_amplitudes()` (raw Cartesian per-order `Ex`/`Ey`,
target 10.1), `.diffraction_angles()` (`theta=None` for evanescent orders,
target 10.2), and `.energy_balance()` (incident/reflected/transmitted/
absorbed/residual in one dict, target 10.3) round out the per-order output
surface alongside the pre-existing `diffraction_efficiencies()`/
`order_classification()`. All three are either a direct reuse of an
already-cited formula (`fields.tangential_e_field`) or a pure composition
of already-validated methods — no new physics formula introduced. Target
10.5 (per-order s/p amplitude conversion) is explicitly **not**
implemented: it requires external (S4/EMpy) validation of the
polarization convention that a bounded investigation this session could
not conclusively achieve (`references.md`'s "Target 10.5 bounded
external-validation attempt", `CONVENTIONS.md`'s Category 10 addendum) —
`complex_amplitudes()`'s Cartesian-only output was chosen specifically
because it needs no s/p convention decision at all.

## "Database Design"

Not applicable — `sougata_solver` has no database and no persistent application
state. The closest analogue is the optional CSV output described above,
which is a one-shot export, not a managed data store.

## Class Diagram (textual)

```
Material
  - epsilon_tensor(wavelength) -> (3,3) complex
  - is_isotropic / is_diagonal
  - .from_nk(...) / .from_permittivity_tensor(...)   [factory classmethods]

Lattice
  - a, b: basis vectors
  - reciprocal_vectors() -> Lk
  - unit_cell_area()

Shape (ABC)                    Pattern
  - fourier_transform(kx,ky)     - background: Material
  - contains(x,y)                 - shapes: list[Shape]
  - signed_distance_normal(x,y)   - containment_tree()
  - area
  ├── Circle(center, radius, material)
  └── Rectangle(center, halfwidth, material, angle)

Layer
  - name, thickness, material | pattern
  - is_uniform() / background_material()

LayerStack(layers, incidence, transmission)
  - wraps `layers` with two semi-infinite half-space Layers

LayerEigenmodes  (data container, produced by eigenmodes.py)
  - q, phi, kp, epsilon_inv, is_scalar_isotropic

SMatrixStack(thicknesses, all_modes)
  - full_smatrix() / partial_smatrix_up_to(i)

PlaneWaveExcitation
  - wavelength, theta, phi, s_amplitude, p_amplitude
  - omega() / k_parallel(n) / incident_field_xy() / incident_mode_amplitude(...)

Simulation(lattice, layers, num_orders, incidence, transmission)
  - solve(excitation) -> SimulationResult

SimulationResult
  - reflectance() / transmittance()
```

## Sequence Diagram (textual) — `Simulation.solve(excitation)`

```
caller -> Simulation.solve(excitation)
  Simulation -> excitation.k_parallel(n_incidence)         : kx0, ky0
  Simulation -> fourier_basis.truncate_fourier_orders(...) : g-vector list
  Simulation -> lattice.reciprocal_vectors()               : Lk
  loop for each layer in layer_stack
    Simulation -> eigenmodes.solve_layer_eigenmodes_*(...) : LayerEigenmodes
  end
  Simulation -> smatrix.SMatrixStack(thicknesses, all_modes)
  SMatrixStack -> smatrix.interface_smatrix(...) / propagation_smatrix(...) / star_product(...)
  Simulation -> excitation.incident_mode_amplitude(modes[0], ...) : a0
  Simulation -> stack.full_smatrix() @ [a0; 0]             : a_transmitted, b_reflected
  Simulation --> caller : SimulationResult
caller -> SimulationResult.reflectance()
  SimulationResult -> fields.z_poynting_flux(...) [x2: incident, reflected]
  SimulationResult --> caller : float
```

## Error Handling

Current, deliberate conventions (keep consistent in new code):

- **Fail loud, fail early, no silent fallbacks.** `Layer.__post_init__`
  raises `ValueError` immediately if neither/both of `material`/`pattern`
  are given (`layer.py:29-30`). `simulation.py` raises
  `NotImplementedError` (not a silent no-op or wrong-but-quiet result) for
  patterned or anisotropic layers it can't yet solve
  (`simulation.py:98,101`) — this is intentional: an RCWA solver that
  silently returns a plausible-looking wrong answer is far worse than one
  that crashes.
- **No broad `except` blocks anywhere in `src/sougata_solver/`.** Keep it that way —
  a caught-and-swallowed `LinAlgError` from a near-singular Toeplitz matrix
  (Phase 2/4 risk) must propagate, not be masked.
- **Validate at construction, not at use.** Prefer raising in `__init__`/
  `__post_init__` (as `Layer` already does) over deep inside a solve call,
  so configuration mistakes surface immediately rather than after an
  expensive sweep.
- **New `NotImplementedError`s for unimplemented phases must name the
  phase**, matching the existing style (`"Patterned layers require Phase
  2+ Fourier factorization"`) so a caller immediately knows whether it's a
  bug or a not-yet-built feature.

## Failure Contract (Category 2 target 2.1)

`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 2 target 2.1: an explicit,
tested inventory of which numerical/input conditions raise which exception
type, and which only emit a `WARNING`, so a caller can tell a genuine bug
from a documented, expected failure mode without reading solver source.
Built by grepping `src/sougata_solver/` for every `raise`/`logger.warning`
call site (not reconstructed from memory) — verified against the code, not
aspirational. Tested in `tests/test_failure_contract.py`.

**`ValueError` — invalid input, raised at construction or immediately on
call (fail loud, fail early, per Error Handling above):**

| Site | Condition |
|------|-----------|
| `Layer.__post_init__` | neither or both of `material`/`pattern` given |
| `smatrix.SMatrixStack.__init__` | `len(thicknesses) != len(all_modes)` |
| `staircase.*` generators | `num_slices < 1` |
| `eigenmodes.solve_layer_eigenmodes_1d` | any `ky != 0` (conical mounting, out of Phase 3 scope) |
| `eigenmodes.solve_layer_eigenmodes_patterned` | `epsilon_hat.shape != (n, n)` |
| `eigenmodes.solve_layer_eigenmodes_patterned_inplane` | any `epsilon_hat_*.shape != (n, n)` |
| `materials.Material.__init__` | `eps` sample is neither a scalar nor a `(3,3)` array |
| `materials._read_numeric_blocks` / `_wavelength_n_k_from_blocks` / `from_refractiveindex_formula_file` | malformed/empty/wrong-column-count material data files |
| `materials.Material.epsilon_tensor` | dispersion callable returns a non-finite value, or a shape mismatching the material's constructed kind, at the queried wavelength (Category 5 target 5.1 — a probe-wavelength-only check at construction can't catch this) |
| `geometry.Lattice`/`Lattice1D`/`Circle`/`Rectangle`/`Ellipse`/`Polygon`/`Slab` construction | non-finite dimension, degenerate (zero-area) lattice, or non-positive shape size (Category 4 target 4.1) |
| `geometry.validate_pattern_fits_lattice` (called from `Simulation.__init__`) | a shape's `2*bounding_radius` is not smaller than the lattice's shortest primitive vector — could overlap its own periodic image (Category 4 target 4.2) |
| `geometry_io.pattern_from_dict`/`pattern_from_json_string`/`pattern_from_json_file` | any malformed field in the JSON pattern schema — missing key, wrong type, unknown shape/unit (Category 4 target 4.6) |
| `smatrix.interface_smatrix` (via `scipy.linalg.lu_factor`'s internal `asarray_chkfinite`, not `numpy.linalg.LinAlgError`) | **exact grazing incidence** (`theta=90 deg`) — the incidence half-space's zeroth-order `q` is exactly `0.0` (confirmed a genuine floating-point coincidence for `n=1`, not merely "very small"; see `tests/test_grazing_incidence.py`, Category 6 target 6.4), so `kp @ phi / q` produces a non-finite matrix that `lu_factor` refuses outright rather than silently returning a wrong answer. Any `theta < 90 deg` is supported, confirmed finite and energy-conserving up to `89.999 deg`. |

**`NotImplementedError` — valid input, unimplemented phase/scope, always
naming the specific phase/target per Error Handling above:**

| Site | Condition |
|------|-----------|
| `simulation.Simulation.solve` (uniform-layer branch) | nonzero `eps_xz/eps_yz/eps_zx/eps_zy` (target 1.5, not yet available) |
| `simulation.Simulation.solve` (patterned-layer branch) | same, for any material in a pattern |
| `fourier_factorization._scalar_value` | anisotropic material passed to the scalar (Phase 2) Fourier-factorization path |
| `fourier_basis.truncate_fourier_orders` | unrecognized `method` string |
| `materials.Material.from_refractiveindex_formula_file` | dispersion formula type other than `"formula 4"` |

**`numpy.linalg.LinAlgError` — never caught, always propagated (per Error
Handling's "no broad `except`" rule) from a numerically singular/
non-converging linear-algebra call:**

| Site | Underlying call | Trigger |
|------|------------------|---------|
| `eigenmodes.solve_layer_eigenmodes_1d` | `np.linalg.solve(epsilon_inv_hat, I)` | exactly-singular inverse-rule Toeplitz |
| `eigenmodes.solve_layer_eigenmodes_patterned` | `np.linalg.solve(epsilon_hat, I)` | exactly-singular direct-rule Toeplitz |
| `eigenmodes.solve_layer_eigenmodes_patterned_inplane` | `np.linalg.solve(epsilon_hat_zz, I)` | exactly-singular `eps_zz` Toeplitz |
| `excitation.PlaneWaveExcitation.incident_mode_amplitude` | `np.linalg.solve(kp @ phi, rhs)` | exactly-singular `kp @ phi` |
| any `eigenmodes.solve_layer_eigenmodes_*` general solver | `np.linalg.eig(op)` | LAPACK `geev` failing to converge (rare; not observed in this project's stress testing, see Phase 4b) |
| `smatrix._solve` (used by `interface_smatrix`/`star_product`) | `scipy.linalg.lu_factor`/`lu_solve` | exactly-singular interface/star-product system |

In practice an *exactly* singular Toeplitz/interface matrix from a physical
(non-adversarially-constructed) structure is rare — the Phase 4b stress
sweep (`tests/test_2d_pillar_stress.py`) pushed condition numbers into the
hundreds without ever hitting exact singularity — but this table exists so
a future `LinAlgError` is recognized as "the documented failure mode for a
degenerate input," not a surprise, and is never silently caught.

**`logging.WARNING` (not raised, not fatal) — numerically-concerning-but-
not-necessarily-wrong, per the Logging Strategy below:**

| Site | Condition | Threshold |
|------|-----------|-----------|
| `eigenmodes.solve_layer_eigenmodes_patterned` | `cond(epsilon_hat)` or `cond(phi)` too large | `ILL_CONDITIONED_THRESHOLD` (`1e4`) |
| `eigenmodes.solve_layer_eigenmodes_patterned_inplane` | `cond(epsilon_hat_zz)` or `cond(phi)` too large | `ILL_CONDITIONED_THRESHOLD` |
| `eigenmodes.solve_layer_eigenmodes_uniform_diagonal`/`_inplane`, `solve_layer_eigenmodes_patterned_inplane` (target 2.4, new) | smallest pairwise eigenvalue gap, relative to `max|q|`, too small | `DEGENERATE_GAP_THRESHOLD` (`1e-6`) — deliberately **not** applied to `solve_layer_eigenmodes_patterned`, see that function's docstring for why (routine `C4v`-symmetry degeneracy in an ordinary case would misfire the warning) |

`materials._wavelength_n_k_from_blocks` uses a plain `print("WARNING: ...")`
for a 2-column (no-`k`-data) material file, predating the
`logging`-module convention adopted in Phase 4b — a pre-existing, minor
deviation from the Logging Strategy's "library never calls `print`" rule,
noted honestly here rather than silently left undocumented; not changed by
this target since target 2.1 is a documentation/test audit, not a
refactor, and changing it isn't otherwise motivated by this session's work.

## Logging Strategy

`src/sougata_solver/` uses the standard `logging` module in exactly one
place as of Phase 4b: `eigenmodes.py`'s module-level
`logger = logging.getLogger(__name__)`, used by
`solve_layer_eigenmodes_patterned` to emit a `WARNING` when
`cond(epsilon_hat)` or `cond(phi)` exceeds `ILL_CONDITIONED_THRESHOLD`
(`1e4`, chosen from a Phase 4b stress-test sweep — see that function's
docstring). Everywhere else, `structures/`/`postprocessing/` scripts still
use plain `print()` (`structures/thin_film/sio2_on_si_thin_film.py`). The
convention going forward:

- **`src/sougata_solver/` (the library) never calls `print` or configures
  logging.** It's a library; logging configuration belongs to the caller
  (a test uses `caplog`, e.g. `tests/test_2d_pillar_stress.py`; an
  interactive script would call `logging.basicConfig()` itself).
- **New diagnostic warnings follow the same pattern**: a module-level
  `logger = logging.getLogger(__name__)`, emitted at `WARNING` level for
  numerically-concerning-but-not-fatal conditions (near-degenerate
  eigenvalues, ill-conditioned matrices at high truncation order) — never
  at `INFO`/`DEBUG` for routine solves, to avoid noise in sweep loops that
  call `solve()` hundreds of times.
- **`structures/*.py`/`postprocessing/*.py` scripts may use `print`** freely
  — they are scripts, not library code, and this matches the existing
  convention.

## Layer/Toeplitz Caching Design (Category 7 targets 7.3/7.4)

**This section exists in tension with `rules.md`'s Performance
Requirements** ("do not introduce caching, memoization, or algorithmic
shortcuts... before Phase 9, unless a specific correctness-validated
capability is measurably too slow to use... and even then, validate the
optimized path against the unoptimized one before trusting it"). Per that
rule, this design is gated on an actual measurement, not implemented
speculatively just because `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` lists it as
a target — the same treatment Category 3 targets 3.4/3.5/3.6 got when a
similar tension came up (`decisions.md` ADR-012).

**Measurement (justifying doing 7.4 now, not deferring to Phase 9), and one
correction made honestly rather than silently along the way.** The first
measurement attempt used a repeated-pattern (DBR-like) stack of `N`
identical patterned layers *within a single `Simulation.solve()` call* and
attributed the entire "extra time per repeated layer" to redundant Toeplitz
reconstruction. Isolating the two costs directly (`toeplitz_matrix` alone
vs. `solve_layer_eigenmodes_patterned` alone, same pattern, `num_orders=49`)
shows that framing was wrong: `toeplitz_matrix` takes ~0.007 s, but the
per-layer eigensolve (dense `2n x 2n` `eig()`, **not** in this target's
cache scope — layer *order* still matters for the eigensolve's downstream
S-matrix cascade even when two layers' permittivity representation is
identical) takes ~0.059 s, roughly 8x more. Caching only the Toeplitz
matrix for that specific scenario (20 identical patterned layers in one
`solve()` call, `num_orders=49`) therefore only saves **~4%** of total
wall-clock time (0.916 s -> 0.880 s), not the much larger number a naive
reading of the first measurement would have implied — recorded here so a
future session doesn't repeat the same measurement mistake.

The actual dominant beneficiary, found by measuring the right scenario
instead: `toeplitz_matrix`/`toeplitz_matrix_component` depend only on
`(pattern, wavelength)`, **not** on `kx`/`ky` (incidence angle) — so an
angle sweep at fixed wavelength (`phases.md`/`tasks.md` Category 8 target
8.3, planned, not yet implemented) reuses the *same* cached Toeplitz
matrix across every sweep point, for even a single patterned layer (no
repeated layers needed):

| `num_orders` | 20-point angle sweep, uncached | cached | reduction |
|---|---|---|---|
| 9  | 0.0639 s | 0.0451 s | ~29% |
| 49 | 1.4372 s | 1.0067 s | ~30% |

(Measured by clearing `Simulation._toeplitz_cache` before every `solve()`
call in the sweep to force the pre-caching behavior, vs. leaving it
populated — the same `Simulation` instance, same excitation angles, same
pattern, only the cache state differs.) This is the real, measured case
that justifies implementing exactly the caching target 7.4 already scopes
("cache one safe artifact (for example, a Toeplitz matrix)"), not a
broader eigensolve or S-matrix cache — a Category 8 angle/wavelength sweep
calling `solve()` hundreds of times would multiply a ~30%-of-solve()-time
saving into a real, non-hypothetical wall-clock reduction, even though the
*within-one-call*, repeated-identical-layer case (the first, wrongly-
framed measurement) turns out to be a much smaller effect on its own.

**Cache scope.** Exactly the Toeplitz matrices built by
`fourier_factorization.toeplitz_matrix`/`toeplitz_matrix_component` —
never the eigenmode solve or S-matrix cascade, which still depend on
layer *order* (interface conditions on either side) even when two layers'
own permittivity representation is identical.

**Cache key.** `(kind, id(pattern), wavelength, ...)`, where `kind` is
`"toeplitz"` (1D/2D-isotropic path, plus an `inverse` bool) or
`"toeplitz_component"` (Category 1 target 1.6's anisotropic path, plus a
`(row, col)` tensor-component pair) — a string tag namespaces the two
call shapes so they can never collide in one dict. `id(pattern)` (Python
object identity), not a value-based structural hash, is used deliberately:
`Pattern`/`Shape` are plain mutable dataclasses (not `frozen=True`), so
they have no `__hash__` to build a value-based key from without adding new
machinery; object-identity caching also has an exactly-right false-positive
rate for this project's actual repeated-layer use case (target 7.2's test
confirms: reusing the *same* `Pattern` object across several `Layer`s is
how a DBR-like repeated unit cell is naturally constructed, and identity
comparison can never wrongly conflate two different `Pattern` objects that
merely happen to hold equal values). `wavelength` is in the key because
`toeplitz_matrix` depends on it (shape Fourier transforms are wavelength-
dependent whenever a shape's material is dispersive); the reciprocal-
lattice truncation `g` is **not** in the key because it is fixed for a
given `Simulation` instance's entire lifetime (`Simulation.__init__`
computes `num_orders`/`truncation` once; `g` itself is recomputed from
those inside `solve()` but is deterministic given them), so it can never
vary between two cache lookups on the same instance.

**Invalidation.** None needed by construction, not a design gap: the cache
is a plain instance attribute (`Simulation._toeplitz_cache: dict`, created
fresh in `__init__`, never module-level — per `rules.md`'s "no hidden
global state, no singletons, no module-level mutable state" rule), so a
new `Simulation` starts with an empty cache and no entry can outlive the
instance it belongs to. Within one instance's lifetime, the key already
captures every quantity `toeplitz_matrix`/`toeplitz_matrix_component`
depend on, so no key can go stale while its `Simulation` exists. The one
caller-facing contract this relies on, made explicit rather than solver-
enforced (matching this project's existing construction-time-validation
philosophy, Category 4 target 4.1, which similarly assumes shapes aren't
mutated after being wired into a `Layer`): **a `Pattern` object, once
passed into a `Layer` that is passed into a `Simulation`, must not have
its `background`/`shapes` mutated in place.** Mutating it and re-solving
the same `Simulation` instance would return a stale cached Toeplitz matrix
for that `Pattern`'s `id()`. This is not a new restriction — Category 4's
construction-time shape/lattice validation already assumes the same thing
implicitly (validated once, at construction, not re-checked on every
solve).

**7.4 implementation note**: see `simulation.py`'s `Simulation._cached_toeplitz`/
`_cached_toeplitz_component`, and `tests/test_layer_cache.py` for the
equivalence-to-uncached test (per this section's own requirement, and
`rules.md`'s "validate the optimized path against the unoptimized one
before trusting it").

## Layer-Wise Absorption Design (Category 7 targets 7.5/7.6)

**Physical quantity.** Per-interior-layer absorbed power, normalized to
incident power (same normalization convention `SimulationResult.reflectance()`/
`transmittance()` already use): the net time-averaged z-directed Poynting
flux entering a layer at its top interface minus the net flux leaving at
its bottom interface. For a lossless layer this is exactly zero by energy
conservation (a strong, free correctness check); for a lossy layer it is
the physical absorbed-power fraction, `A_i`, satisfying
`R + T + sum(A_i) = 1` for a stack with no diffracted orders escaping as
loss elsewhere (`testing.md`'s Physical-Invariant Testing tier already
names this identity as the target the project doesn't yet satisfy without
this capability — see `tests/test_stress_regression.py`'s docstring, which
explicitly flags Category 7 targets 7.5/7.6 as the missing piece).

**Formula — reused, not invented, per `rules.md` AI Coding Rule 1's
preference for building on already-validated blocks (the same treatment
ADR-015 gave `interior_amplitudes`).** No new physics formula is written
here at all; layer-wise absorption is a direct combination of three
already-implemented, already-oracle/consistency-validated pieces (Category
9 / Phase 7):

```text
For each interior layer i (all_modes index 1..len(all_modes)-2):
    (a_top, b_top) = interior_amplitudes(stack.partial_smatrix_up_to(i), n2, a0, b_reflected)
    (a_bot, b_bot) = propagate_amplitudes(modes_i.q, thickness_i, a_top, b_top)
    net_top = sum(z_poynting_flux(omega, modes_i.q, modes_i.kp, modes_i.phi, a_top, b_top))
    net_bot = sum(z_poynting_flux(omega, modes_i.q, modes_i.kp, modes_i.phi, a_bot, b_bot))
    A_i = (net_top - net_bot).real / incident_power.real
```

`z_poynting_flux` already returns `(forward, backward)` with the
forward/backward interference cross-term folded symmetrically into each
half (`fields.py:51-56`, the `diff`/`np.conj(diff)` split) — so
`forward + backward` is already the genuine net total z-flux at one
reference plane with both a forward and a backward wave simultaneously
present, exactly what's needed here (not an approximation or a new
derivation of that split).

**Why this reuses existing infrastructure instead of a volumetric `Im(eps)`
integral.** The textbook alternative — integrating `omega * Im(eps) *
|E|^2` over the layer's volume — would need a new spatial-integration
formula (a new source of physics-formula risk, exactly what Rule 1 warns
against) and would need reconciling against this project's already-found
`z_poynting_flux` factor-of-2 convention (`CONVENTIONS.md`, Category 9
target 9.6's finding). The flux-divergence formula above needs no new
formula at all — every piece (`interior_amplitudes`, `propagate_amplitudes`,
`z_poynting_flux`) is already independently validated
(`tests/test_field_reconstruction.py`), and Poynting-flux divergence
equaling absorbed power is a direct restatement of energy conservation,
not a new physical claim.

**Validation method (before exposing the API, per target 7.5's own
wording).** The energy-balance identity itself is the validation: `R + T +
sum(A_i) == 1` for a genuinely lossy stack (reusing
`tests/test_stress_regression.py`'s already-vetted `eps = -396+80j` lossy-
metal fixture, whose sign convention was already checked against
`CONVENTIONS.md`'s `d/dt -> -i*omega` passivity requirement), and `A_i ==
0` (to numerical precision) for every layer in a lossless stack — both
independent of any external oracle, since no published per-layer
absorption benchmark exists in any vendored `REFERENCE/` repo (confirmed
by the same `phase-reference-picker`-style audit prior categories used;
Category 7's own register entry lists no such benchmark as "already
present").

**API placement.** `SimulationResult.layer_absorption() -> list[float]`,
one entry per interior (finite-thickness) layer in stack order — mirrors
`diffraction_efficiencies()`/`order_classification()`'s existing pattern
of a `SimulationResult` method built from already-stored `all_modes`/`a0`/
`b_reflected`, rather than a free function requiring the caller to keep
the original `Simulation` object around. This needs `SimulationResult` to
also store `thicknesses` (a new field, populated by `Simulation.solve()`
from data it already computes locally) — the one small, purely-additive
API extension this target requires.

## Linear-Algebra Baseline & Factorization-Reuse Design (Category 12 targets 12.1/12.3)

**12.1 baseline measurements** (`profiling/baseline_profile.py`, run on
this development machine — absolute numbers are machine-dependent and not
asserted in any test, per `rules.md`'s Performance Requirements; the
*qualitative* scaling behavior is the load-bearing result):

| Stage | `num_orders` | measured (min of 5) |
|---|---|---|
| Isolated 2D-pillar eigensolve | 9 / 25 / 49 / 81 | 0.83 / 2.52 / 11.9 / 132 ms |
| Isolated matrix-solve (`_solve`, random `(2n,2n)`) | 9 / 25 / 49 / 81 | 0.30 / 0.25 / 0.81 / 161 ms |
| `Simulation.solve()`, thin-film | 1 | 5.0 ms |
| `Simulation.solve()`, 1D grating | 9 | 11.8 ms |
| `Simulation.solve()`, 2D pillar | 9 / 49 | 9.7 / 115 ms |

**Finding, consistent with and now quantifying ADR-016's earlier
observation**: the eigensolve is the dominant cost at moderate-to-large
`num_orders` (steep, worse-than-linear growth: 9→81 orders is a ~9x
Fourier-order increase but a ~160x eigensolve-time increase), not the
isolated matrix-solve step (which stays comparatively small until the
same large sizes, where LU-factorization cost catches up for the same
underlying reason — both operations are dense `O(n^3)`-class LAPACK calls
on the same growing matrix dimension). This directly informs 12.3 and
12.5 below.

**12.3 factorization-reuse design.** Audited every explicit linear-solve
call site in the S-matrix cascade (`smatrix.py::interface_smatrix`,
`star_product`, `interior_amplitudes`) for repeated-coefficient-matrix
reuse opportunities:

- `interface_smatrix`'s `_solve(ta, identity)`: `ta` is a fresh
  `0.5*(P+Q)` combination for every distinct layer-pair, never repeated
  within one cascade — no reuse opportunity, **except** the case already
  handled since Phase 1: `_is_trivial_interface` short-circuits an
  interface between two layers sharing identical `(q, phi, kp)` to a
  literal identity S-matrix, skipping the linear solve entirely. This is
  itself the safe, already-shipped instance of "intra-solve factorization
  reuse" this target asks about — not a new opportunity, but worth
  recording as the answer, not overlooking it because it predates
  Category 12.
- `star_product`'s two solves (`_solve(t1, ...)`, `_solve(t2, ...)`):
  `t1 = I - a01@b10` and `t2 = I - b10@a01` are different matrices in
  general (equal only if `a01`/`b10` commute, not assumed) — no reuse.
- `interior_amplitudes`'s `_solve(s11, ...)`: `s11` is a different partial
  S-matrix block per interior interface queried (`layer_absorption()`
  calls it once per interior layer) — no reuse across layers.

**Conclusion**: no further *S-matrix-level* factorization-reuse
opportunity exists without changing the cascade's fundamental per-
interface/per-layer structure (a materially larger, riskier change,
out of this target's small-target scope). The 12.1 measurements point
instead to a *different* opportunity, one layer up: since the eigensolve
(not the matrix-solve) dominates, and Category 7's Toeplitz-matrix cache
(ADR-016) already caches the Toeplitz matrix *feeding into* each
eigensolve but not the eigensolve's own result, a natural future
extension — **explicitly not implemented here**, since 12.3 asks only for
the design, matching how target 7.3 was design-only before 7.4
implemented it, and no corresponding "12.6 implementation" target exists
in this category — would be an instance-scoped `LayerEigenmodes` cache on
`Simulation`, keyed the same way ADR-016's Toeplitz cache is
(`(kind, id(pattern)-or-material-identity, wavelength)`, since `kx`/`ky`
are already fixed for a `Simulation` instance's lifetime, the same
invariant ADR-016 relies on) — directly useful for the same fixed-
wavelength angle-sweep scenario ADR-016 was measured against, since
`kx`/`ky` genuinely do change with angle and *do* affect the eigensolve
(unlike the Toeplitz matrix), so this specific extension would need
angle-independent-only reuse within one `Simulation.solve()` call for
repeated identical patterned layers, not across an angle sweep — a
narrower, different benefit than ADR-016's, worth stating precisely
rather than overclaiming. Revisit if a future session profiles a
concrete workload where this narrower reuse would actually matter.

## Configuration, CLI, and Export Design (Category 15 targets 15.2/15.5/15.7)

**15.2 configuration schema** (`config.py`). A single JSON document
describes one `Simulation` + one `PlaneWaveExcitation`: `unit`, `lattice`,
`materials` (a name → `{eps_re, eps_im, source}` dict, reusing
`geometry_io.py`'s existing material-dict shape unchanged rather than
inventing a second one), `layers` (each either `material`-referencing or
`pattern`-holding — a `pattern` value is passed straight through to
`geometry_io.pattern_from_dict`, so patterned-layer JSON has exactly one
schema project-wide, not two), `incidence`/`transmission` (material
names), `num_orders`, `truncation`, and `excitation`. Every numeric length
field is expressed in the top-level `unit` (`"m"`/`"um"`/`"nm"`) and
converted to meters before any `Simulation`/`Layer`/`Material` object is
constructed — the same unit-scale convention `geometry_io.py` already
uses for imported geometry, not a new one.

**Validation ordering (15.3)**: `simulation_from_dict` only *constructs*
objects — `Layer`, `Material`, `Lattice`/`Lattice1D`, `Simulation`,
`PlaneWaveExcitation` — and never calls `.solve()`. Every malformed-input
check (missing key, wrong type, unknown material name referenced by a
layer/incidence/transmission, non-positive `num_orders`, unrecognized
`truncation`) therefore necessarily raises before any numerical work
starts, satisfying target 15.3's ordering requirement structurally (by
what the function does, not by an extra guard) —
`tests/test_config.py::test_validation_never_reaches_a_numerical_solve`
pins this directly by monkeypatching `Simulation.solve` to fail loudly if
ever reached from an invalid config.

**15.4 configuration runner**: `simulation_from_dict`/
`simulation_from_json_string`/`simulation_from_json_file` (`config.py`)
are the runner — a caller still calls `.solve(excitation)` explicitly
afterward, keeping "build and validate a `Simulation`" and "run it"
distinct steps, mirroring `Simulation`'s own existing two-phase
construct-then-solve API rather than adding a third, config-specific
solving entry point. `tests/test_config.py::test_config_reproduces_anti_reflection_coating_example`
confirms a config-file-driven run reproduces
`structures/thin_film/anti_reflection_coating.py`'s result to `1e-12`.

**15.5/15.6 CLI design and implementation** (`cli.py`). One subcommand,
`run <config.json>`, three exit codes (`0` solved, `2` invalid
config/file, `1` any other solver failure) — kept as two distinct
non-zero codes so a caller (e.g. a shell script or CI step) can tell "your
input was wrong" apart from "the solve itself failed" without parsing
stderr text. Output reuses `output_paths.py`'s existing
`outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/` convention (via
`run_output_dir`/`write_run_metadata`, the same functions every
`structures/*.py` example already calls) rather than inventing a second
output-location convention; `--output-dir` overrides it for
scripted/CI use.

**15.7 NumPy export** (`export.py`). `export_sweep_npz`/`load_sweep_npz`
serialize a `sweep.SweepResult` to a plain `.npz` archive:
`parameter_values`/`reflectance`/`transmittance` as numeric arrays, plus
a JSON-encoded `metadata` string array (not a pickled object array) — so
`np.load` never needs `allow_pickle=True`, keeping this export path free
of the untrusted-deserialization risk class `rules.md`'s Security Rules
already flag for `eval`/`exec`/`pickle` on file-sourced data. Scope is
deliberately limited to *numeric* (1D) `parameter_values` — a discrete/
labeled sweep (e.g. Category 8's polarization-Jones-tuple sweep) raises
rather than silently truncating, since `.npz` needs a homogeneous
per-array dtype and there is no lossless numeric encoding of an arbitrary
label without guessing a convention no caller asked for.

**15.8 HDF5 — evaluated and deferred**, see ADR-026: no dependency
justifies itself against this project's actual (small, flat) result
shapes yet.

## Plotting Design (Category 16 targets 16.1-16.7)

**16.1 plot data contract**: every function in the new `plotting.py`
takes plain arrays, dataclasses, or already-computed result objects
(`geometry.Pattern`/`Lattice`, raw field-grid arrays, a
`SimulationResult.diffraction_efficiencies()` dict) — never a bare
`Simulation` — and no function calls `.solve()`. Every function returns
`(fig, ax)` rather than saving/showing, keeping the module a pure
"data in, figure out" library with no filesystem or display side
effects; saving/showing remains the caller's job, matching
`postprocessing/*.py`'s existing pattern (`plot_thin_film_rt.py`'s
`fig.savefig(...)`/`plt.show()`). `matplotlib` is imported lazily inside
each function (not at module level), so importing
`sougata_solver.plotting` doesn't force the `matplotlib` dependency onto
every library user — the same lazy-import pattern
`postprocessing/plot_field_cross_section.py` already used.

**16.2 unit-cell/layer-stack plots**: `plot_unit_cell` rasterizes a
preview grid using each `Shape.contains(x, y)` (already implemented by
every shape class) rather than special-casing matplotlib patches per
shape type — one implementation covers `Circle`/`Rectangle`/`Ellipse`/
`Polygon`/`Slab` uniformly, and it respects `Pattern`'s own documented
"later shapes take precedence" rule by iterating shapes in reverse
order and taking the first (i.e. topmost/latest) match. This is
explicitly a **preview raster for visualization only** — the solve
itself still uses `pattern`'s analytic Fourier transforms unmodified.
`plot_layer_stack` draws semi-infinite (`math.inf`-thickness) layers as
a fixed-height hatched band, since infinity has no natural bar height.

**16.3-16.5**: `plot_rt_spectrum` formalizes
`postprocessing/plot_thin_film_rt.py`'s existing ad hoc R-vs-wavelength
plot (same axis labels/style), generalized to optionally show `T` and an
`R+T` conservation trace and a `metadata` annotation.
`plot_harmonic_convergence` marks the exact index
`sweep.find_convergence_index` already selected, so the plot visually
matches the validated criterion rather than a human eyeball guess.
`plot_diffraction_orders` sorts orders by `(g1, g2)` for a deterministic
bar order across repeated calls (a `dict`'s insertion order is otherwise
incidental).

**16.6/16.7**: `plot_field_intensity` formalizes
`postprocessing/plot_field_cross_section.py`'s `pcolormesh` intensity
plot (`|E|^2`) into a reusable function. `plot_field_phase` uses a
cyclic colormap (`twilight`, not a sequential one) since phase wraps at
`+-pi` and a sequential map would show a spurious discontinuity there.
`plot_poynting_vector` visualizes already-computed `Sx`/`Sz` values (no
flux computation of its own, per the target 16.1 data contract) using
the no-`0.5`-factor real-space flux convention Category 9 target 9.6
found and documented (`CONVENTIONS.md`, `troubleshooting.md`).

## CI, Static Analysis, and Regression-Guard Design (Category 17 targets 17.2/17.3/17.5/17.6)

**17.2/17.3 CI workflows** (`.github/workflows/`): `ci.yml` runs
`pytest -m "not slow"` plus `ruff check .` on `windows-latest` across
the three Python versions `pyproject.toml`'s `requires-python` spans
(3.10-3.12), on every push/PR to `main` and on manual dispatch.
`slow-tests.yml` runs `pytest -m slow` on a weekly schedule (catches
slow drift without gating every push) plus manual dispatch, per this
target's own "schedule or manually trigger" wording. Both are ordinary
GitHub Actions YAML with no custom infrastructure.

**17.5 static analysis**: `ruff` (`[tool.ruff]`/`[tool.ruff.lint]` in
`pyproject.toml`), scoped deliberately narrow for this pass — `select =
["F", "E7"]` (pyflakes real-bug checks: unused imports/variables,
undefined names; a small pycodestyle statement-style subset) rather than
the full pycodestyle/import-sort rule families, which would flag
formatting-only nitpicks unrelated to this target's "fix the baseline
before making them required" wording. `line-length = 120` matches this
project's actual long-docstring/citation style (`rules.md`'s
Documentation Standards), not ruff's 88-char default. Running it found
and fixed 24 real baseline issues (22 unused imports, safe autofix; 2
genuinely dead local variables in `src/sougata_solver/vectorized.py` and
`tests/test_field_reconstruction.py`, verified unused via a direct grep
for other references before removing, not assumed).

**17.6 performance regression guard**: see `decisions.md` ADR-028 for
the full design rationale (a same-run relative-scaling ratio, never an
absolute wall-clock threshold, per `rules.md`'s Performance
Requirements).
