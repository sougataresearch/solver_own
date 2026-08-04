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

The intended import surface (already `__all__`-exported from
`src/sougata_solver/__init__.py`):

```python
from sougata_solver import Material, Lattice, Circle, Rectangle, Pattern, Layer, LayerStack
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
