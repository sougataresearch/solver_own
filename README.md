# sougata_solver

A pure-Python **Rigorous Coupled-Wave Analysis (RCWA)** solver for periodic
electromagnetic structures — thin films, multilayer stacks, 1D-periodic
lamellar gratings (trenches), and 2D-periodic patterned layers (vias,
pillars), including tapered/sloped sidewalls.

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
   sloped sidewalls.
3. Correctly solve 2D-periodic patterned layers (via, pillar arrays), also
   with sloped sidewalls.
4. Support dispersive, absorbing, and (eventually) anisotropic materials.
5. Validate every new capability against an independent oracle (S4, analytic
   Fresnel, or a published benchmark table) before trusting it.
6. Stay a small, readable, single-author codebase — not a framework.

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
  vendored `RigorousCoupledWaveAnalysis.jl`/`convmat2D.py` algorithm). Not
  yet wired into `simulation.py` — the general (non-uniform) eigenmode
  solver that consumes it is Phase 4.

Planned (see [`phases.md`](phases.md) for the full roadmap):
- General (non-uniform) eigenmode solver so patterned layers actually work
  end to end (currently `simulation.py` raises `NotImplementedError` for
  any patterned layer)
- 1D lamellar gratings (trench)
- 2D patterned layers (via, pillar)
- Tapered/sloped sidewalls via staircase layer discretization
- Anisotropic materials
- Real-space field reconstruction and cross-section plotting
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
├── pyproject.toml
├── src/sougata_solver/        see src/sougata_solver/README.md for the module map
│   ├── materials.py         permittivity models (isotropic + tensor)
│   ├── geometry.py           Lattice, Shape (Circle/Rectangle), Pattern
│   ├── fourier_basis.py       G-vector truncation
│   ├── fourier_factorization.py  Toeplitz permittivity matrices (Phase 2, done)
│   ├── layer.py                Layer, LayerStack, LayerEigenmodes
│   ├── eigenmodes.py           per-layer eigenmode solve (uniform today)
│   ├── smatrix.py               interface + propagation S-matrices, star product
│   ├── excitation.py            plane-wave decomposition, incident amplitude
│   ├── fields.py                  Poynting flux, tangential field reconstruction
│   ├── polarimetry.py             Jones/Mueller (reused by postprocessing/)
│   ├── simulation.py               top-level orchestration
│   └── output_paths.py             outputs/YYYY_MM_DD/HH_MM_SS_<run>/ helper
├── tests/                    pytest suite + `tests/oracles/` -- see tests/README.md
├── structures/                YOU RUN THESE -- see structures/README.md
│   └── thin_film/                uniform multilayer stacks (Phase 1, done)
│                                (trench/ and via/ subfolders land with Phase 3/4)
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
blocks. (Trench and via/pillar structures get their own `structures/trench/`
and `structures/via/` folders once Phase 3/4 land — see `phases.md`.)

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

Run the test suite:

```bash
pytest
```

## Future Improvements

See [`phases.md`](phases.md) for the complete, ordered roadmap (Fourier
factorization, 1D trenches, 2D via/pillar, tapered sidewalls, anisotropy,
field visualization, expanded validation, optional performance backend).
