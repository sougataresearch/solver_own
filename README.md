# pyrcwa

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

`pyrcwa` solves Maxwell's equations for plane-wave illumination of a stack
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
  half-spaces (`Layer`, `LayerStack` — [`src/pyrcwa/layer.py`](src/pyrcwa/layer.py))
- Dispersive materials from constant, callable, or refractiveindex.info-style
  CSV `n,k` data (`Material` — [`src/pyrcwa/materials.py`](src/pyrcwa/materials.py))
- Arbitrary incidence angle/azimuth and elliptical polarization
  (`PlaneWaveExcitation` — [`src/pyrcwa/excitation.py`](src/pyrcwa/excitation.py))
- Numerically stable Redheffer star-product S-matrix cascading
  ([`src/pyrcwa/smatrix.py`](src/pyrcwa/smatrix.py))
- Reflectance/transmittance via Poynting flux
  ([`src/pyrcwa/fields.py`](src/pyrcwa/fields.py))
- Jones/Mueller polarimetry
  ([`src/pyrcwa/polarimetry.py`](src/pyrcwa/polarimetry.py))
- Analytic in-plane Fourier transforms for `Circle` and `Rectangle` shapes
  with nested-shape subtraction, ready to be consumed by the patterned-layer
  solver ([`src/pyrcwa/geometry.py`](src/pyrcwa/geometry.py))
- Circular G-vector truncation for Fourier-order selection
  ([`src/pyrcwa/fourier_basis.py`](src/pyrcwa/fourier_basis.py))

Planned (see [`phases.md`](phases.md) for the full roadmap):
- Fourier-factorization core + general (non-uniform) eigenmode solver so
  patterned layers actually work end to end (currently
  `simulation.py` raises `NotImplementedError` for any patterned layer)
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
pyrcwa/
├── README.md            this file
├── PRD.md                product requirements
├── architecture.md       system architecture
├── design.md              detailed design (algorithms, API, error handling)
├── rules.md               coding/testing/git/AI rules
├── phases.md              roadmap phases
├── tasks.md                atomic task checklist per phase
├── memory.md               live project status for future sessions
├── decisions.md            architecture decision record (ADR)
├── testing.md              testing strategy
├── deployment.md           environment/CI/release process
├── references.md            literature + reference-implementation index
├── troubleshooting.md      known numerical gotchas
├── pyproject.toml
├── src/pyrcwa/
│   ├── materials.py         permittivity models (isotropic + tensor)
│   ├── geometry.py           Lattice, Shape (Circle/Rectangle), Pattern
│   ├── fourier_basis.py       G-vector truncation
│   ├── layer.py                Layer, LayerStack, LayerEigenmodes
│   ├── eigenmodes.py           per-layer eigenmode solve (uniform today)
│   ├── smatrix.py               interface + propagation S-matrices, star product
│   ├── excitation.py            plane-wave decomposition, incident amplitude
│   ├── fields.py                  Poynting flux, tangential field reconstruction
│   ├── polarimetry.py             Jones/Mueller (reused by postprocessing/)
│   └── simulation.py               top-level orchestration
├── tests/                    pytest suite + `tests/oracles/` (analytic references)
├── structures/                YOU RUN THESE: define a lattice/layer stack/materials
│                                and run the solver -- e.g. sio2_on_si_thin_film.py,
│                                custom_multistack.py (copy this one for a new stack)
├── postprocessing/             YOU RUN THESE SECOND: take a structures/ script's raw
│                                 output and derive Jones/Mueller matrices, ellipsometric
│                                 angles, and (planned) RI/thickness extraction
└── scripts/                   (currently empty; ad hoc utility scripts)
```

## Installation

```bash
cd pyrcwa
python -m venv .venv
.venv\Scripts\activate        # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -e ".[dev]"
```

## Usage

Run a structure end to end (builds the geometry, runs the solver, prints
R/T):

```bash
python structures/sio2_on_si_thin_film.py
```

For a multi-layer stack of your own materials, copy
`structures/custom_multistack.py` and edit its numbered `EDIT` blocks.

For Jones/Mueller/ellipsometric-angle analysis, run the matching
`structures/*_ellipsometry_run.py` script first (it saves raw field data to
a CSV), then the corresponding script in `postprocessing/` (it loads that
CSV and derives the Jones matrix, Mueller matrix, and Psi/Delta — no
re-solving):

```bash
python structures/sio2_on_si_ellipsometry_run.py
python postprocessing/jones_mueller_ellipsometry.py
```

Or use the library directly:

```python
from pyrcwa.materials import Material
from pyrcwa.layer import Layer
from pyrcwa.geometry import Lattice
from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.simulation import Simulation
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
