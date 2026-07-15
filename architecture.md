# Architecture — pyrcwa

## High-Level Architecture

`pyrcwa` is a **pipeline**, not a service: given a static description of a
lattice + layer stack + materials + one excitation, it produces mode
amplitudes, from which R/T, polarimetry, and (later) fields are derived.
There is no persistent state, no server, no I/O beyond optional file reads
(material CSV data) and file writes (example scripts saving a CSV of
results). This shape is deliberate — see `decisions.md` — a physics kernel
should be a pure function of its inputs so that oracle-comparison testing
is meaningful and reproducible.

```
Material(s)  ─┐
Lattice       ├──►  Simulation.solve(excitation)  ──►  SimulationResult
Layer(s)     ─┘            │                              │
                            │                              ├── reflectance()
              ┌─────────────┼─────────────┐                └── transmittance()
              ▼             ▼             ▼
     fourier_basis    eigenmodes      smatrix
   (truncate G-set)  (per-layer q,   (interface +
                       phi, kp)       propagation S,
                                      Redheffer star
                                      product cascade)
```

## Module Breakdown

| Module | Responsibility | Status |
|--------|-----------------|--------|
| `materials.py` | Permittivity representation (scalar or 3×3 tensor), dispersion via callables | done |
| `geometry.py` | `Lattice` (reciprocal-vector math), `Shape` (`Circle`, `Rectangle`) with analytic Fourier transforms, `Pattern` (ordered shapes + containment/subtraction tree) | geometry primitives done; not yet wired to the solver |
| `fourier_basis.py` | G-vector (Fourier order) truncation, circular selection matching S4's `gsel.c` | done |
| `layer.py` | `Layer` (uniform or patterned), `LayerStack` (with incidence/transmission half-spaces), `LayerEigenmodes` (result container) | done |
| `eigenmodes.py` | Per-layer eigenmode solve: closed-form for uniform isotropic layers today; general non-uniform solve is Phase 4 | uniform case done |
| `smatrix.py` | Interface S-matrices, propagation S-matrices, Redheffer star-product cascade (`SMatrixStack`) | done, dimension-agnostic (works unchanged once Phase 3/4 land) |
| `excitation.py` | Plane-wave s/p decomposition, incident-mode-amplitude inversion | done |
| `fields.py` | z-Poynting flux (R/T power), tangential E-field reconstruction at one interface | done for R/T; full real-space reconstruction is Phase 7 |
| `polarimetry.py` | Jones/Mueller matrix construction from simulation results | done |
| `simulation.py` | Orchestration: builds the Fourier-order set, solves every layer's eigenmodes, cascades the S-matrix stack, solves for transmitted/reflected amplitudes given an incident excitation | done for uniform layers; raises `NotImplementedError` for any patterned layer today (`simulation.py:98`) |

## Data Flow

1. **Setup** (once per structure): construct `Material`s, a `Lattice`, and a
   list of `Layer`s (uniform or — once Phase 3/4 land — patterned);
   construct a `Simulation`.
2. **Per-wavelength/angle solve** (`Simulation.solve(excitation)`):
   a. Compute the zeroth-order in-plane wavevector `(kx0, ky0)` from the
      incidence medium and excitation angle (`excitation.py::k_parallel`).
   b. Truncate the Fourier-order (G-vector) set for the lattice
      (`fourier_basis.truncate_fourier_orders`), producing the full
      `(kx, ky)` arrays for every retained diffraction order.
   c. For every layer in the stack, solve its eigenmode problem
      (`eigenmodes.solve_layer_eigenmodes_uniform` today; will dispatch to
      a patterned/1D/2D solver once later phases land), producing
      `LayerEigenmodes(q, phi, kp, ...)`.
   d. Cascade all per-layer/interface S-matrices into one full-stack
      S-matrix (`smatrix.SMatrixStack`), using the Redheffer star product
      so evanescent-mode growth is never represented as an absolute
      (ill-conditioned) transfer matrix.
   e. Convert the excitation's desired `(Ex, Ey)` at the zeroth order into
      a forward mode-amplitude vector `a0`
      (`excitation.incident_mode_amplitude`), solve
      `[a_transmitted; b_reflected] = S_full @ [a0; 0]`.
   f. Return a `SimulationResult`, from which `reflectance()`/
      `transmittance()` are derived via `fields.z_poynting_flux`.
3. **(Planned, Phase 7)** Field reconstruction at an arbitrary intermediate
   depth uses `SMatrixStack.partial_smatrix_up_to` to get local mode
   amplitudes, then inverse-Fourier-sums them onto a real-space grid.

## Component Responsibilities

- **`Simulation`** is the only object that knows how to assemble the whole
  pipeline; it owns no physics itself beyond bookkeeping (which G-vectors,
  which layer's modes go where). This keeps every individual physics
  formula testable in isolation (see `testing.md`).
- **`LayerEigenmodes`** is the contract between `eigenmodes.py` and
  everything downstream (`smatrix.py`, `fields.py`, `excitation.py`): as
  long as a per-layer solver produces `(q, phi, kp)` in the conventions
  documented in `layer.py`, the rest of the pipeline is agnostic to whether
  the layer was uniform, 1D-patterned, or 2D-patterned. This is *why* Phase
  3/4 (trench, via/pillar) are additive — they add new producers of
  `LayerEigenmodes`, not new consumers.
- **`SMatrixStack`** is dimension-agnostic and already correct for any
  future eigenmode solver, since it only operates on `(q, phi, kp)` and
  `thickness` — no change needed for Phases 3-6.

## External Services

None. `pyrcwa` has no network calls, no database, no external API
dependency at runtime. The only "external" inputs are local files: material
`n,k` CSV data (`examples/03_sio2_on_si.py::material_from_csv`) and the
vendored reference repositories (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl`)
used **only** as offline validation oracles during development/testing, never
imported by `pyrcwa` itself.

## Technology Choices

| Choice | Why |
|--------|-----|
| Pure Python + NumPy + SciPy, no compiled extension | Correctness/readability prioritized over speed through Phase 8 (see PRD non-functional requirements); `scipy.linalg.lu_factor`/`lu_solve` already used for robust interface-matrix inversion (`smatrix.py::_solve`) instead of `numpy.linalg.inv` |
| `dataclasses` for `Layer`, `LayerEigenmodes`, `SimulationResult`, shapes | Plain data containers, no framework machinery, matches the "no hidden state" non-functional requirement |
| S-matrix (Redheffer star product) over raw transfer-matrix cascading | Transfer matrices blow up numerically for evanescent (decaying) modes in thick/absorbing layers; S-matrices remain well-conditioned — this is the same reason S4 and virtually every modern RCWA implementation use S-matrices, not the classic (numerically unstable) T-matrix method |
| Analytic shape Fourier transforms (`jinc`/`sinc`) over raster+FFT of a pixelized mask | More accurate for smooth boundaries (circular vias) since there's no pixelization error; matches S4's approach; deliberately different from Meent/TORCWA's raster+FFT approach (see `decisions.md`) |
| `pytest` + a small `tests/oracles/` package of independent analytic references | Keeps "is this correct" testable without needing a live S4/Julia install for every CI run; S4-cross-check tests (Phase 4) will be an additional, optionally-skipped tier |

## Directory Structure

See `README.md`'s Folder Structure section — reproduced here for
completeness since architecture and layout are tightly coupled in a project
this size:

```
pyrcwa/src/pyrcwa/
├── materials.py
├── geometry.py
├── fourier_basis.py
├── layer.py
├── eigenmodes.py
├── smatrix.py
├── excitation.py
├── fields.py
├── polarimetry.py
└── simulation.py
```

No sub-packages yet. A `fourier_factorization.py` module is planned for
Phase 2 (see `phases.md`), and a `staircase.py` helper for Phase 5 — both
new top-level modules in `src/pyrcwa/`, not new sub-packages, to keep the
import graph flat as long as the module count stays in the low teens.
**Revisit this if the module count exceeds ~15-18** — at that point, group
into sub-packages (e.g. `pyrcwa/geometry/`, `pyrcwa/solve/`) rather than
letting `src/pyrcwa/` become an unstructured flat pile.

## Scalability Considerations

This is a single-machine, single-user numerical library, not a scaled
service — "scalability" here means:

- **Fourier-order count (`num_orders`)**: the dominant cost driver. Every
  per-layer eigenmode solve and every S-matrix operation scales with
  `num_orders` (matrix sizes are `O(num_orders)` per axis, so
  `O(num_orders^3)` per matrix operation for the general non-uniform
  eigenproblem in Phase 4). No algorithmic mitigation is planned before
  Phase 9; convergence studies (Phase 8) should report the smallest
  `num_orders` that achieves acceptable accuracy per structure type, since
  that is the real lever available today.
- **Wavelength/angle sweeps**: currently a Python `for` loop calling
  `Simulation.solve` once per point (see `examples/03_sio2_on_si.py`).
  Phase 9 explicitly proposes vectorizing this in NumPy (batched matrix
  ops) before considering any GPU/autodiff backend — see `phases.md` Phase
  9 and `decisions.md`.
- **Staircase sidewall discretization (Phase 5)**: adds a multiplicative
  factor of `N` (slice count) layers to the stack; since `SMatrixStack`
  cascades layers sequentially, cost is linear in `N`, not per-layer
  matrix-size-dependent — cheap relative to `num_orders` scaling.

## Security Considerations

`pyrcwa` has a minimal attack surface: no network, no database, no
authentication, no user-supplied code execution. The two considerations
that do apply:

- **File parsing (material CSV ingestion)**: `material_from_csv` in
  `examples/03_sio2_on_si.py` reads user-controlled file paths — acceptable
  for a local research tool run by its own author; if this ever becomes a
  shared/public tool, add path validation and handle malformed CSVs with a
  clear error rather than a raw `ValueError`/parse exception (currently
  acceptable per `PRD.md`'s solo-research-tool scope).
- **No `eval`/`exec`/`pickle` of untrusted data anywhere in the codebase** —
  keep it that way; if a future phase adds structure-definition file
  ingestion (e.g. a JSON/YAML scene format), use a schema-validated parser,
  never `eval` on structure descriptions.
