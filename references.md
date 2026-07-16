# References — sougata_solver

Index of vendored reference implementations and literature this project
validates against. Update when a new phase cites a new source.

## Vendored Reference Implementations (`Solver_own/`, siblings of `sougata_solver/`)

| Repo | Language | Role |
|------|----------|------|
| [`S4`](../S4) | C++ / Lua | **Primary oracle.** Every non-trivial formula in `eigenmodes.py`, `smatrix.py`, `fields.py` is checked against a specific, cited line range in `S4/S4/rcwa.cpp` or `S4/S4r/StarProduct.hpp`. Also the source of the `geometry.py` Fourier-transform/subtraction-rule convention (`S4/S4/pattern/pattern.c`). Planned as a subprocess cross-check oracle for Phase 4 (2D patterned layers) if buildable in this environment (unverified — see `memory.md` Known Issues). |
| [`EMpy`](../EMpy) | Python | **Second Phase 1 numerical oracle**, as of the `tests/oracles/empy_tmm.py` cross-check: `EMpy/EMpy/transfer_matrix.py`'s `IsotropicTransferMatrix.solve` (lines 52-134) transcribed by hand (not imported — reference-only per this table's own rule). Chosen over `Rigorous-Coupled-Wave-Analysis/TMM_functions` (same role, but that repo mixes plotting into physics and has an anisotropic function referencing undefined free variables) as the cleanest standalone isotropic-TMM code among the vendored options. The transcription surfaced three real bugs in the EMpy source itself (`abs()` vs `Re()` transmittance prefactor, `arcsin`-branch failure for absorbing media, and a `D@P@inv(D)` sign ambiguity) — see `empy_tmm.py`'s module docstring for the full account; all three are fixed in the transcription, not reproduced. |
| [`RigorousCoupledWaveAnalysis.jl`](../RigorousCoupledWaveAnalysis.jl) | Julia | Secondary reference implementation (ETM/SRCWA submodules); available as a possible independent cross-check if S4 isn't usable for a given Phase 4/6 validation case. |
| [`Rigorous-Coupled-Wave-Analysis`](../Rigorous-Coupled-Wave-Analysis) | Python | Educational/research-grade RCWA (Rumpf-formulation) covering TMM, 1D gratings, and 2D gratings — considered for Phase 1's oracle role but passed over in favor of `EMpy` (see above) on code quality grounds. Still the best-mapped source for later phases: `RCWA_1D_examples/1D_Grating_Gaylord_TE.py`/`TM.py` directly implement the Moharam & Gaylord (1995) benchmark already cited below as the planned Phase 3 oracle (run it directly rather than hand-copying paper tables); `RCWA_2D_examples/RCWA_photonic_circle_spectra.py` is a structural analog for Phase 4's via/pillar case; `TMM_functions/anisotropic.py` is a candidate Phase 6 oracle. Not yet used as a numerical oracle — revisit when Phase 3/4/6 land. |
| [`EMTutorial`](../EMTutorial) | JCMsuite project files | Not code — FEM tutorial *geometries and setups* (thin-film DBR, through-silicon-via scatterometry, gratings, metasurfaces) used as realistic target structures for `sougata_solver`'s own `structures/` scripts (see `phases.md` Phase 8). Specifically referenced: `EMTutorial/ThinFilmsAndMultilayers/DistributedBraggReflector`, `EMTutorial/Scatterometry/ThroughSiliconVia`. |
| [`NK_FILE`](../NK_FILE) | CSV data | Si/SiO2 refractive-index data consumed by `structures/thin_film/sio2_on_si_thin_film.py::material_from_csv`. |

## External Tools Referenced (not vendored, not dependencies)

- **Meent** (KC-ML2) — open-source RCWA with NumPy/JAX/PyTorch backends,
  1D+2D lattice support, autodiff for topology optimization. Referenced in
  `decisions.md` ADR-002 (raster+FFT Fourier factorization, contrasted with
  `sougata_solver`'s analytic approach) and ADR-006 (Phase 9 GPU/autodiff backend
  precedent). Not a dependency; not imported.
- **TORCWA** — PyTorch-based, GPU-accelerated batched RCWA with autograd.
  Same referenced role as Meent above, specifically for Phase 9's optional
  vectorized/GPU backend design.
- **JCMsuite** — commercial FEM solver; source of the `EMTutorial/`
  vendored project files (tutorials only, not the solver itself).

## Literature (to be added to as phases cite specific results)

- **Moharam, M. G., & Gaylord, T. K. (1995)**, "Formulation for stable and
  efficient implementation of the rigorous coupled-wave analysis of binary
  gratings," *J. Opt. Soc. Am. A* — planned source for Phase 3's 1D
  binary-grating diffraction-efficiency validation benchmark. **Exact
  table/equation to be pinned down and cited by line/table number when
  Phase 3 actually implements the test** — do not cite this paper's
  numbers from memory; locate and transcribe the specific benchmark table
  before writing the test, per `rules.md`'s Documentation Standards.
- **Li, Lifeng (1996)**, "Use of Fourier series in the analysis of
  discontinuous periodic structures," *J. Opt. Soc. Am. A* — the
  foundational reference for the direct-vs-inverse-rule Fourier
  factorization distinction underlying Phase 2's `epsilon_hat` vs.
  `epsilon_inv_hat` Toeplitz construction (see `design.md`, Algorithm 3).
  Cite the specific rule/equation when Phase 2's docstrings are written,
  not just this general reference.

## How to Add a Reference

When a new phase's implementation cites a source (per `rules.md`'s
mandatory Documentation Standards), add it here too, with enough detail
(file + line range, or author/year/equation) that a future session can
re-locate it without re-searching from scratch.
