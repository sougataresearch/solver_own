# Detailed Design — pyrcwa

`pyrcwa` has no database and no UI in the traditional sense, so those
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

### 2. General (non-uniform) eigenmode solve (patterned layers — Phase 4)

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
poorly-conditioned eigenvectors; the mitigation is the Phase 4 validation
requirement (cross-check against S4 itself, which has already solved this
problem correctly) rather than a from-scratch re-derivation of stability
fixes.

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

## "UI/UX" — Example-Script and Plotting Surface

There is no GUI. The user-facing surface is:

1. **Library API**, imported into small, single-purpose scripts (see
   `examples/`, numbered `01_...` through `04_...` today).
2. **Console output**: examples print a table of wavelength/R/T/A to
   stdout during the sweep (see `examples/03_sio2_on_si.py::main`).
3. **CSV output**: examples optionally save results
   (`OUTPUT_CSV_PATH` pattern in `examples/03_sio2_on_si.py`).
4. **(Phase 7, planned)**: `matplotlib`-based cross-section field-intensity
   plots for trench/via structures — the first genuinely visual output
   this project will produce. When implemented, follow the pattern already
   established: a `examples/NN_*.py` script that both computes and plots,
   not a separate plotting framework/module, since there's no repeated
   plotting logic yet to justify one.

## "API Design" — Public Python API

The intended import surface (already `__all__`-exported from
`src/pyrcwa/__init__.py`):

```python
from pyrcwa import Material, Lattice, Circle, Rectangle, Pattern, Layer, LayerStack
```

Plus, imported directly from their submodules (not yet re-exported at
top level — worth revisiting once Phase 3/4 land and usage patterns
stabilize):

```python
from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.simulation import Simulation, SimulationResult
```

**Typical call sequence** (already the pattern in every `examples/*.py`):

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

Not applicable — `pyrcwa` has no database and no persistent application
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
- **No broad `except` blocks anywhere in `src/pyrcwa/`.** Keep it that way —
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

## Logging Strategy

There is currently **no logging module usage anywhere in `src/pyrcwa/`** —
only `print()` in example scripts (`examples/03_sio2_on_si.py`). This is
appropriate for the library core (a numerical function should not have
side-effecting log output — it should raise or return, full stop) but
should be formalized as follows going forward:

- **`src/pyrcwa/` (the library) never calls `print` or configures
  logging.** It's a library; logging configuration belongs to the caller.
- **If/when diagnostic visibility is genuinely needed inside the library**
  (e.g. reporting Toeplitz condition number or eigenvalue-degeneracy
  warnings once Phase 4's general eigensolver lands), use the standard
  `logging` module with a module-level `logger = logging.getLogger(__name__)`,
  emitted at `WARNING` level for numerically-concerning-but-not-fatal
  conditions (e.g. near-degenerate eigenvalues, ill-conditioned Toeplitz
  matrix at high truncation order) — never at `INFO`/`DEBUG` for routine
  solves, to avoid noise in sweep loops that call `solve()` hundreds of times.
- **`examples/*.py` scripts may use `print`** freely — they are scripts, not
  library code, and this matches the existing convention.
