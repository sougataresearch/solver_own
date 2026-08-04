# sougata_solver

A pure-Python **Rigorous Coupled-Wave Analysis (RCWA)** solver for periodic
electromagnetic structures — thin films, multilayer stacks, 1D-periodic
lamellar gratings (trenches), 2D-periodic patterned layers (vias, pillars,
now including elliptical and simple-polygon cross-sections) including
tapered/sloped sidewalls, and (for diagonal and in-plane-coupled tensors)
anisotropic materials in both uniform and patterned layers, with dispersive
materials from analytic models (Sellmeier, Cauchy, Lorentz, Drude,
Drude-Lorentz) or tabulated `n,k` data.

This is a from-scratch, from-first-principles implementation, not a wrapper
around an existing solver. Every non-trivial formula in the codebase is
checked against a named, line-numbered source before being trusted — either
the vendored [S4](../S4) C++ reference implementation, an independent
analytic solution (Fresnel/TMM), or a classic published benchmark (e.g.
Moharam & Gaylord). This discipline is documented in [`rules.md`](rules.md)
and is the project's core engineering principle: **a physics solver is only
as trustworthy as its validation**, and silent numerical bugs (wrong sign
convention, wrong branch cut, wrong Fourier factorization rule) are the
dominant failure mode in this domain, not crashes.

## Project Overview

`sougata_solver` solves Maxwell's equations for plane-wave illumination of a stack
of periodic (or uniform) layers, returning reflected/transmitted diffraction
efficiencies, polarimetric response (Jones/Mueller), and — in a later phase
— reconstructed real-space field maps. It targets the same class of problems
as commercial tools like JCMsuite (see the vendored tutorials in
[`../EMTutorial`](../EMTutorial)) but restricted to structures that are
periodic in the lateral direction(s): thin-film stacks, distributed Bragg
reflectors, 1D gratings/trenches, and 2D via/pillar arrays.

## Objectives

1. Correctly solve uniform multilayer stacks (thin film / DBR) — **done**.
2. Correctly solve 1D-periodic lamellar gratings (trench, line/space) with
   sloped sidewalls — **done**.
3. Correctly solve 2D-periodic patterned layers (via, pillar arrays), also
   with sloped sidewalls — **done**.
4. Support dispersive, absorbing, and anisotropic materials — dispersive/
   absorbing **done** (constant/tabulated `n,k` plus five analytic dispersion
   models: Sellmeier, Cauchy, Lorentz, Drude, Drude-Lorentz); anisotropic
   **partially done** (diagonal and in-plane-coupled tensors, uniform and
   patterned; longitudinal coupling explicitly deferred — see Features
   below).
5. Validate every new capability against an independent oracle (S4, RCWA.jl,
   analytic Fresnel, or a published benchmark table) before trusting it.
6. Stay a small, readable, single-author codebase — not a framework.
7. Beyond geometry/materials/physics correctness, treat the *numerical
   robustness contract* itself as a deliverable — a documented failure
   contract (which inputs raise which exception vs. only warn), eigenvalue/
   conditioning diagnostics, and deterministic behavior across a sweep —
   **done** for the scope shipped so far (see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
   Category 2).

## Features

Current (Phase 1, shipped):
- Arbitrary-thickness multilayer stacks with semi-infinite incidence/exit
  half-spaces (`Layer`, `LayerStack` — [`src/sougata_solver/layer.py`](src/sougata_solver/layer.py))
- Dispersive materials from constant, callable, or refractiveindex.info-style
  CSV `n,k` data (`Material` — [`src/sougata_solver/materials.py`](src/sougata_solver/materials.py))
- Arbitrary incidence angle/azimuth and elliptical polarization
  (`PlaneWaveExcitation` — [`src/sougata_solver/excitation.py`](src/sougata_solver/excitation.py))
- Numerically stable Redheffer star-product S-matrix cascading
  ([`src/sougata_solver/smatrix.py`](src/sougata_solver/smatrix.py))
- Reflectance/transmittance via Poynting flux
  ([`src/sougata_solver/fields.py`](src/sougata_solver/fields.py))
- Jones/Mueller polarimetry
  ([`src/sougata_solver/polarimetry.py`](src/sougata_solver/polarimetry.py))
- Analytic in-plane Fourier transforms for `Circle` and `Rectangle` shapes
  with nested-shape subtraction
  ([`src/sougata_solver/geometry.py`](src/sougata_solver/geometry.py))
- Circular G-vector truncation for Fourier-order selection
  ([`src/sougata_solver/fourier_basis.py`](src/sougata_solver/fourier_basis.py))

Current (Phase 2, shipped):
- Toeplitz permittivity matrix construction (direct and inverse-rule) for
  patterned layers — `pattern_epsilon_hat`, `toeplitz_matrix`
  ([`src/sougata_solver/fourier_factorization.py`](src/sougata_solver/fourier_factorization.py)),
  validated against two independent numerical references (from-scratch
  rasterize-and-sum, and an FFT-of-rasterized-mask reproduction of the
  vendored `RigorousCoupledWaveAnalysis.jl`/`convmat2D.py` algorithm).

Current (Phase 3, shipped):
- 1D-periodic lamellar gratings (trench, line/space) — `Lattice1D`, `Slab`
  ([`src/sougata_solver/geometry.py`](src/sougata_solver/geometry.py)),
  `truncate_fourier_orders_1d`
  ([`src/sougata_solver/fourier_basis.py`](src/sougata_solver/fourier_basis.py)),
  `solve_layer_eigenmodes_1d`
  ([`src/sougata_solver/eigenmodes.py`](src/sougata_solver/eigenmodes.py)),
  with diffraction efficiencies via `SimulationResult.diffraction_efficiencies()`.
  Validated against `tests/oracles/rcwa_1d_gaylord.py` (Moharam & Gaylord
  1995) and the energy-conservation invariant.

Current (Phase 4a/4b, shipped):
- General (non-uniform) 2D-periodic patterned-layer eigenmode solver —
  `solve_layer_eigenmodes_patterned`
  ([`src/sougata_solver/eigenmodes.py`](src/sougata_solver/eigenmodes.py)),
  removing `simulation.py`'s prior `NotImplementedError` for any patterned
  layer. Validated against an independent eigenvalue oracle transcribed
  from `RigorousCoupledWaveAnalysis.jl` (agrees to ~1e-12), energy
  conservation, and a ky=0-reduces-to-1D cross-check.
- Condition-number `WARNING` logging (`eigenmodes.ILL_CONDITIONED_THRESHOLD`)
  for near-degenerate/ill-conditioned patterned-layer cases, after a
  deliberate high-contrast/high-order stress sweep found no catastrophic
  failure.

Current (Phase 5, shipped):
- Tapered/sloped sidewalls via staircase (z-discretized) layer
  approximation — `staircase_circle_layers`, `staircase_rectangle_layers`,
  `staircase_slab_layers`
  ([`src/sougata_solver/staircase.py`](src/sougata_solver/staircase.py)).
  Validated by a zero-taper regression to the already-oracle-validated
  single-uniform-layer result, energy conservation, and
  convergence-vs-`num_slices` studies.

Current (Phase 6, Category 1 targets 1.3-1.4/1.6-1.8, shipped; see
[`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`](COMMERCIAL_RCWA_ATOMIC_TARGETS.md)):
- Uniform diagonal-tensor and in-plane-coupled (`eps_xx, eps_xy, eps_yx,
  eps_yy, eps_zz`) anisotropic layers —
  `solve_layer_eigenmodes_uniform_diagonal`,
  `solve_layer_eigenmodes_uniform_inplane`, validated against a closed-form
  birefringence benchmark and an independent `RigorousCoupledWaveAnalysis.jl`-
  derived oracle (`tests/test_anisotropic_uniform.py`,
  `tests/test_anisotropic_inplane.py`).
- Patterned (2D-periodic) anisotropic layers —
  `solve_layer_eigenmodes_patterned_inplane`,
  `fourier_factorization.toeplitz_matrix_component`
  (`tests/test_anisotropic_patterned.py`).
- A deterministic mode-ordering policy for near-degenerate eigenvalues
  (`eigenmodes._canonical_mode_order`,
  `tests/test_anisotropic_degeneracy.py`) and public propagating/evanescent
  mode classification (`eigenmodes.classify_propagating`,
  `SimulationResult.order_classification()`,
  `tests/test_mode_classification.py`).
- **Not yet supported**: longitudinal tensor coupling
  (`eps_xz/eps_yz/eps_zx/eps_zy`) — evaluated and explicitly deferred, no
  citable + independently-benchmarkable formulation found (see
  `references.md`'s "Target 1.5 bounded literature search").

Current (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 2, Numerical methods,
targets 2.1-2.5, shipped):
- A documented failure contract (`design.md`) covering every
  `ValueError`/`NotImplementedError`/`LinAlgError`/`WARNING` condition in
  the solver, backed by `tests/test_failure_contract.py`.
- Per-solve eigenvalue/mode-conditioning diagnostics
  (`layer.EigenmodeDiagnostics`, attached as `LayerEigenmodes.diagnostics`)
  and a configurable small-eigenvalue-gap `WARNING`
  (`eigenmodes.DEGENERATE_GAP_THRESHOLD`), alongside the existing
  ill-conditioning `WARNING` from Phase 4b.
- Deterministic mode ordering across a small wavelength sweep for the
  anisotropic dense eigensolvers (reusing Category 1 target 1.7's
  canonical-ordering policy).

Current (Category 3, Fourier factorization, targets 3.1-3.6, shipped):
- A documented rule inventory (`design.md`) recording which direct/
  inverse/numerical-inverse Fourier-factorization rule every solver branch
  actually uses, backed by `tests/test_fourier_factorization_rules.py`.
- Fixed high-contrast 1D and 2D convergence fixtures with measured (not
  assumed) convergence-vs-harmonic-order data
  (`tests/test_fourier_convergence.py`).
- Fast Fourier Factorization / Normal Vector Method feasibility evaluated
  and explicitly deferred (`decisions.md` ADR-012) — the current 2D solver
  keeps ordinary Laurent's-rule Toeplitz construction.

Current (Category 4, Geometry engine, targets 4.1-4.7, shipped):
- Construction-time validation for `Lattice`/`Lattice1D`/`Circle`/
  `Rectangle`/`Slab` and a unit-cell self-overlap policy
  (`geometry.validate_pattern_fits_lattice`, wired into `Simulation`).
- `Ellipse` and `Polygon` shape primitives (`geometry.py`), both with
  closed-form analytic Fourier transforms (no raster/FFT — see
  `decisions.md` ADR-013), each with an end-to-end example
  (`structures/via/elliptical_pillar.py`, `structures/via/triangular_pillar.py`).
- A minimal, safe JSON `Pattern`-import format
  (`geometry_io.py::pattern_from_dict`/`pattern_from_json_file`).
- A general geometry-to-layer-slices interface (`staircase.slice_profile`),
  of which the existing tapered-sidewall generators are now thin wrappers.

Current (Category 5, Material models, targets 5.1-5.8, shipped):
- Construction- and call-time `Material` validation (tensor shape, finite
  values, dispersion-callback output).
- Five analytic dispersion models — `Material.from_sellmeier`/`from_cauchy`
  (validated against BK7's published Sellmeier coefficients/index) and
  `Material.from_lorentz`/`from_drude`/`from_drude_lorentz` (validated
  against Rakić et al. (1998)'s published Lorentz-Drude metal model, with
  `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM` coefficient
  presets ready to use).
- Optional `Material.source` citation metadata, threaded through every
  `from_*` classmethod and into serialized `run_metadata.txt` output.

Planned (see [`phases.md`](phases.md) for the full roadmap):
- Real-space field reconstruction and cross-section plotting (Phase 7)
- Expanded systematic validation sweep across all geometry types (Phase 8)
- Optional vectorized/GPU/autodiff backend (later; see `decisions.md`)

## Target Users

Currently a single user (project owner) doing scatterometry / thin-film /
via-trench-pillar EM simulation work, using the vendored `S4`, `EMpy`, and
`RigorousCoupledWaveAnalysis.jl` repositories as reference/validation
oracles rather than dependencies. Not (yet) intended for external/public
users — see [`PRD.md`](PRD.md) for scope.

## Tech Stack

- **Python ≥ 3.10**, pure Python + NumPy + SciPy only (no compiled
  extensions, no GPU dependency at this stage)
- `pytest` for the test suite
- `setuptools` (src-layout package, see `pyproject.toml`)
- No web framework, no database, no UI — this is a library driven by
  Python scripts, split into `structures/` (build a geometry + run the
  solver) and `postprocessing/` (derive Jones/Mueller matrices, and,
  eventually, RI/thickness extraction, from already-computed raw results)

## Folder Structure

```
sougata_solver/
├── README.md            this file
├── PRD.md                product requirements
├── architecture.md       system architecture
├── design.md              detailed design (algorithms, API, error handling)
├── rules.md               coding/testing/git/AI rules
├── phases.md              roadmap phases
├── tasks.md                atomic task checklist per phase
├── memory.md               live project status for future sessions
├── progress_log.md          dated log of discussions + action items (new 2026-07-19)
├── decisions.md            architecture decision record (ADR)
├── testing.md              testing strategy
├── deployment.md           environment/CI/release process
├── references.md            literature + reference-implementation index
├── troubleshooting.md      known numerical gotchas
├── CONVENTIONS.md           frozen field/phasor/S-matrix/tensor conventions
├── COMMERCIAL_RCWA_ATOMIC_TARGETS.md   fine-grained Phase 6+ target checklist
├── pyproject.toml
├── src/sougata_solver/        see src/sougata_solver/README.md for the module map
│   ├── materials.py         permittivity models (isotropic + tensor) + analytic dispersion
│   │                          models (Sellmeier/Cauchy/Lorentz/Drude/Drude-Lorentz)
│   ├── geometry.py           Lattice, Lattice1D, Shape (Circle/Rectangle/Ellipse/Polygon/Slab),
│   │                          Pattern, construction-time validation
│   ├── geometry_io.py         minimal safe JSON Pattern-import format (parser only)
│   ├── fourier_basis.py       G-vector truncation (2D circular + 1D)
│   ├── fourier_factorization.py  Toeplitz permittivity matrices, scalar + per-tensor-component
│   ├── layer.py                Layer, LayerStack, LayerEigenmodes, EigenmodeDiagnostics
│   ├── staircase.py             slice_profile (general layer-slicing interface) +
│   │                             tapered-sidewall staircase generators built on it
│   ├── eigenmodes.py           per-layer eigenmode solve: uniform, 1D-patterned (Phase 3),
│   │                            2D-patterned (Phase 4a/4b), and anisotropic (Phase 6);
│   │                            conditioning/degeneracy WARNING diagnostics
│   ├── smatrix.py               interface + propagation S-matrices, star product
│   ├── excitation.py            plane-wave decomposition, incident amplitude
│   ├── fields.py                  Poynting flux, tangential field reconstruction
│   ├── polarimetry.py             Jones/Mueller (reused by postprocessing/)
│   ├── simulation.py               top-level orchestration
│   └── output_paths.py             outputs/YYYY_MM_DD/HH_MM_SS_<run>/ helper
├── tests/                    pytest suite (393 tests) + `tests/oracles/` -- see tests/README.md
├── structures/                YOU RUN THESE -- see structures/README.md
│   ├── thin_film/                uniform multilayer stacks (Phase 1, done)
│   ├── trench/                    1D lamellar gratings, tapered ridges (Phase 3/5, done)
│   └── via/                        2D via/pillar arrays, tapered/elliptical/polygon
│                                    pillars (Phase 4/5, Category 4, done)
└── postprocessing/             YOU RUN THESE SECOND: take a structures/ script's raw
                                  output and derive Jones/Mueller matrices, ellipsometric
                                  angles, and (planned) RI/thickness extraction
```

Folder-level READMEs with more detail:
[`src/sougata_solver/README.md`](src/sougata_solver/README.md) ·
[`structures/README.md`](structures/README.md) ·
[`tests/README.md`](tests/README.md)

## Installation

```bash
cd sougata_solver
python -m venv .venv
.venv\Scripts\activate        # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -e ".[dev]"
```

## Usage

Run a structure end to end (builds the geometry, runs the solver, prints
R/T):

```bash
python structures/thin_film/sio2_on_si_thin_film.py
```

For a multi-layer stack of your own materials, copy
`structures/thin_film/custom_multistack.py` and edit its numbered `EDIT`
blocks. Trench (1D grating) and via/pillar (2D patterned) structures live in
`structures/trench/` and `structures/via/`, including tapered-sidewall
variants (`tapered_trench.py`, `tapered_via.py`, `tapered_pillar.py`) — see
[`structures/README.md`](structures/README.md).

For Jones/Mueller/ellipsometric-angle analysis, run the matching
`structures/thin_film/*_ellipsometry_run.py` script first (it saves raw field
data to a CSV), then the corresponding script in `postprocessing/` (it loads
that CSV and derives the Jones matrix, Mueller matrix, and Psi/Delta — no
re-solving):

```bash
python structures/thin_film/sio2_on_si_ellipsometry_run.py
python postprocessing/jones_mueller_ellipsometry.py
```

### Output files

Every script that saves a result (CSV today; plots later) writes into
`outputs/YYYY_MM_DD/HH_MM_SS_<script-name>/`, via
`src/sougata_solver/output_paths.py`: one date folder per day, and inside it
one timestamped subfolder per run, so a day's runs stay together but two
different scripts — or two runs of the same script — never overwrite each
other. `postprocessing/` scripts read the most recent matching file across
all run subfolders, so a same-day run-then-postprocess workflow needs no path
editing, and postprocessing still finds the input if run on a later day.
`outputs/` is gitignored.

Or use the library directly:

```python
from sougata_solver.materials import Material
from sougata_solver.layer import Layer
from sougata_solver.geometry import Lattice
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.simulation import Simulation
import math

air = Material("air", 1.0)
sio2 = Material("SiO2", 1.46**2)
lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused for uniform layers

sim = Simulation(lattice, [Layer("SiO2", 50e-9, material=sio2)],
                  num_orders=1, incidence=air, transmission=air)
result = sim.solve(PlaneWaveExcitation(wavelength=550e-9, theta=0.0, phi=0.0))
print(result.reflectance(), result.transmittance())
```

A dispersive metal via a published (not hand-fit) model, or a JSON-defined
pattern, both use the same `Material`/`Pattern` objects as above:

```python
from sougata_solver.materials import Material, RAKIC_GOLD
gold = Material.from_drude_lorentz("Au", *RAKIC_GOLD)   # Rakić et al. 1998

from sougata_solver.geometry_io import pattern_from_json_string
pattern = pattern_from_json_string('{"background": {"eps_re": 1.0}, '
    '"shapes": [{"type": "circle", "center": [0.35, 0.35], "radius": 0.18, '
    '"material": {"eps_re": 12.11}}]}')  # units default to meters
```

Run the test suite:

```bash
pytest                # fast suite
pytest -m slow        # convergence/benchmark studies (several minutes)
```

## Future Improvements

See [`phases.md`](phases.md) for the complete, ordered roadmap and
[`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`](COMMERCIAL_RCWA_ATOMIC_TARGETS.md) for
the fine-grained target checklist. Phases 1-5 are shipped, and
`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Categories 1-5 (mathematical foundation/
anisotropy, numerical methods, Fourier factorization, geometry engine,
material models) are all shipped except Category 1 target 1.5 (longitudinal
tensor coupling, explicitly deferred pending a citable formulation).
Remaining: Category 6 (boundary conditions/excitation) through Category 19
(future extensions) at the atomic-target level — notably Category 7's
layer-wise absorption (no `absorbance()`/`A` computation exists yet: a
lossy material's `R+T` is correctly `<1`, but nothing independently derives
the missing `A`) — plus real-space field reconstruction and cross-section
plotting (Phase 7), an expanded validation suite and example gallery
(Phase 8), and an optional vectorized/GPU/autodiff backend (Phase 9,
later).
