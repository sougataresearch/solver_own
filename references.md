# References — sougata_solver

Index of vendored reference implementations and literature this project
validates against. Update when a new phase cites a new source.

## Vendored Reference Implementations (`Solver_own/`, siblings of `sougata_solver/`)

| Repo | Language | Role |
|------|----------|------|
| [`S4`](../S4) | C++ / Lua | **Primary oracle.** Every non-trivial formula in `eigenmodes.py`, `smatrix.py`, `fields.py` is checked against a specific, cited line range in `S4/S4/rcwa.cpp` or `S4/S4r/StarProduct.hpp`. Also the source of the `geometry.py` Fourier-transform/subtraction-rule convention (`S4/S4/pattern/pattern.c`). Phase 3 additionally transcribed `S4/S4/rcwa.cpp::SolveLayerEigensystem` (lines 684-827, the *general* non-uniform eigenoperator — S4 has no separate 1D/TE/TM code path, confirmed by grepping `rcwa.cpp` for those terms) into `eigenmodes.solve_layer_eigenmodes_1d`, specialized to `ky=0`; and `S4/S4/fmm/fmm_closed.cpp:110-132`'s "1D proper FFF rule" branch (`0==Lr[2] && 0==Lr[3]`), which sets `Epsilon2`'s TM-like block to `inv(Epsilon_inv)` — confirms Li's (1996) inverse-rule correction is required there, not the direct-rule Toeplitz. **Phase 4a additionally transcribed `fmm_closed.cpp`'s adjacent true-2D branch (lines 133-139, 162-163)** into `eigenmodes.solve_layer_eigenmodes_patterned`: for a genuine 2D pattern without a polarization basis, S4 uses ordinary Laurent's rule for *both* `Epsilon2` blocks (`block_diag(epsilon_hat, epsilon_hat)`) and feeds `kp` the numerical matrix-inverse `inv(epsilon_hat)`, **not** the inverse-rule Toeplitz — Li's rule is 1D-only in S4's own source. This was caught only after a first draft wrongly reused the 1D branch's formula for the 2D case (see `phases.md` Phase 4a's Status for the full account) — recorded here as a reminder that reading only the branch matching a chosen benchmark case is not sufficient; the adjacent branch(es) need reading too before citing a source as "the general case." S4 cross-check for Phase 4a itself is deferred (not run — see `memory.md`'s Known Issues); `tests/oracles/rcwa_2d_pillar.py` explicitly flags this rather than fabricating a match. |
| [`EMpy`](../EMpy) | Python | **Second Phase 1 numerical oracle**, as of the `tests/oracles/empy_tmm.py` cross-check: `EMpy/EMpy/transfer_matrix.py`'s `IsotropicTransferMatrix.solve` (lines 52-134) transcribed by hand (not imported — reference-only per this table's own rule). Chosen over `Rigorous-Coupled-Wave-Analysis/TMM_functions` (same role, but that repo mixes plotting into physics and has an anisotropic function referencing undefined free variables) as the cleanest standalone isotropic-TMM code among the vendored options. The transcription surfaced three real bugs in the EMpy source itself (`abs()` vs `Re()` transmittance prefactor, `arcsin`-branch failure for absorbing media, and a `D@P@inv(D)` sign ambiguity) — see `empy_tmm.py`'s module docstring for the full account; all three are fixed in the transcription, not reproduced. `EMpy/EMpy/RCWA.py` (~897 lines) was surveyed for Phase 4a's oracle gap and ruled out: 1D-only (no 2D/pillar support), plain Laurent's rule, no test suite, and an author-acknowledged hack (`RCWA.py:304-316`, randomizes `kx` when `cond(A)>1e10`) — not a clean reference for anything beyond what `transfer_matrix.py` already provided. **Category 5 targets 5.2/5.3 (Sellmeier/Cauchy dispersion models)**: `EMpy/EMpy/materials.py::RefractiveIndex.__from_sellmeier` (lines 118-127) transcribed into `materials.Material.from_sellmeier`; the same module's docstring (lines 65-72) worked Cauchy-form example transcribed into `materials.Material.from_cauchy`. |
| [`RigorousCoupledWaveAnalysis.jl`](../RigorousCoupledWaveAnalysis.jl) | Julia | Secondary reference implementation (ETM/SRCWA submodules). **Phase 4a eigenoperator oracle**, as of `tests/oracles/rcwa_2djl_eigenvalues.py`: `src/Common/Common.jl:57-99` (`eigenmodes(...,l::PatternedLayer)`, isotropic branch) hand-transcribed (Julia not installed in this environment, `which julia` fails, so not run directly). A **structurally different eigenoperator derivation** from S4's `Epsilon2 @ kp - coupling` route (direct Maxwell-curl elimination into one matrix `M`, different field basis) — feeding it this project's own already-Phase-2-validated `epsilon_hat` isolates the eigenoperator-construction step specifically, and confirms `solve_layer_eigenmodes_patterned`'s `q^2` eigenvalues to ~1e-12 across several `num_orders`/angle/pattern cases, once two conventions are reconciled (both confirmed empirically, documented in the oracle's docstring): RCWA.jl normalizes `kx,ky` by `k0` (`grids.jl:78-84`), and its eigenvalues come out an overall sign-flip of this project's `q^2` (opposite time convention, same class of footnote as the other oracles' sign-convention notes). No independently-published numeric benchmark exists in this repo (`test/runtests.jl:70-111` uses `rand()` params, self-consistency only) — this closes the eigenoperator-correctness gap, not the full-R/T external-oracle gap (see `tests/oracles/rcwa_2d_pillar.py`, still open for Phase 4b). **Category 5 targets 5.4-5.6 (Lorentz/Drude/Drude-Lorentz dispersion models)**: `src/BasicMaterials/rakic.jl` (`LorentzDrude`, lines 14-21, plus the published Au/Ag/Al/Ti coefficient tables at lines 24-45) transcribed into `materials.Material.from_lorentz`/`from_drude`/`from_drude_lorentz` and the module-level `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM` constants — citing A. D. Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, "Optical properties of metallic films for vertical-cavity optoelectronic devices," *Appl. Opt.* 37, 5271-5283 (1998), the same citation `rakic.jl`'s own trailing comment gives. |
| [`Rigorous-Coupled-Wave-Analysis`](../Rigorous-Coupled-Wave-Analysis) | Python | Educational/research-grade RCWA (Rumpf-formulation) covering TMM, 1D gratings, and 2D gratings — considered for Phase 1's oracle role but passed over in favor of `EMpy` (see above) on code quality grounds. **Phase 3 oracle**, as of `tests/oracles/rcwa_1d_gaylord.py`: `RCWA_1D_examples/1D_Grating_Gaylord_TE.py` (lines 139-257) and `1D_Grating_Gaylord_TM.py` (lines 169-246) hand-transcribed (not imported — both files mix `matplotlib` calls into the physics with no `__main__` guard, and neither hard-codes a paper table number, confirmed by direct read; the scripts compute their own spectral sweep). Both cite Moharam, Grann, Pommet & Gaylord (1995) in their own header comments — see the Literature entry below. One transcription fix applied (same "fix, don't reproduce" precedent as `empy_tmm.py`): `1D_Grating_Gaylord_TM.py` hard-codes `n_groove=3.48` (same as the ridge, not a real grating) instead of `1.0` like the TE file; the oracle module uses the TE file's binary-grating geometry for both polarizations. One caveat *not* fixed, because it wasn't independently re-derivable with confidence: `1D_Grating_Gaylord_TM.py`'s own module docstring says "STILL NOT WORKING YET" — confirmed via a direct convergence sweep that both this project's `solve_layer_eigenmodes_1d` and the transcribed oracle converge to the same TM value, but only at high `num_orders` (TE agrees at `num_ord~15`; TM needs `num_ord` in the hundreds) — see `tests/test_1d_grating.py`'s TM test docstring for the full account. `RCWA_2D_examples/RCWA_photonic_circle_spectra.py` (+ `RCWA_functions/run_RCWA_simulation.py::run_RCWA_2D`, lines 13-146) was surveyed for Phase 4a's oracle gap: genuine dense-eigensolve 2D RCWA, S4-equivalent structure, plain Laurent's-rule Toeplitz only (`convmat2D.py`, no inverse-rule matrix for 2D — third independent confirmation, alongside S4 and `RigorousCoupledWaveAnalysis.jl`, that this project's Phase 4a fix's Laurent's-rule choice isn't S4-idiosyncratic). No hard-coded reference numbers (only an unresolved "compare with Fan JOSA B" provenance comment) — would need ~7 helper modules pulled together to transcribe into a full R/T oracle; not done, left as a Phase 4b option (see `tests/oracles/rcwa_2d_pillar.py`). `TMM_functions/anisotropic.py` remains a candidate Phase 6 oracle, not yet used. **Category 5 target 5.5 (Drude model)**: `TMM_examples/TMM_Drude.py:67` (`drude_eps = 1 - omega_p**2/(omega**2 + 1j*omega*gamma)`) cross-checked (algebraically identical formula structure, confirmed directly not assumed) against `RigorousCoupledWaveAnalysis.jl`'s Rakic Drude term before transcribing `materials.Material.from_drude`. |
| [`EMTutorial`](../EMTutorial) | JCMsuite project files | Not code — FEM tutorial *geometries and setups* (thin-film DBR, through-silicon-via scatterometry, gratings, metasurfaces) used as realistic target structures for `sougata_solver`'s own `structures/` scripts (see `phases.md` Phase 8). Specifically referenced: `EMTutorial/ThinFilmsAndMultilayers/DistributedBraggReflector`, `EMTutorial/Scatterometry/ThroughSiliconVia`. |
| [`NK_FILE`](../NK_FILE) | CSV data | Si/SiO2 refractive-index data consumed by `structures/thin_film/sio2_on_si_thin_film.py::material_from_csv`. |

## Phases With No Vendored-Repo Source

### Phase 6 anisotropy reference audit (2026-08-03)

This audit was performed before changing any anisotropic solver code. The
candidate implementations were compared rather than defaulting to S4:

- `S4/S4/S4.cpp` (uniform anisotropic assembly, approximately lines
  1866-1906) is a useful convention reference for the in-plane tensor block
  and epsilon_zz inverse. It is **not** a full-3x3 reference: its Lua API
  documentation (`doc/source/lua_api.rst`, material-tensor arguments)
  explicitly says that xz, yz, zx, and zy components are ignored.
- `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl` (the
  `AnisotropicLayer` eigenmode path, approximately lines 134-157) and
  `src/Common/materials.jl` support the five-component in-plane tensor
  `[epsilon_xx, epsilon_xy, epsilon_yx, epsilon_yy, epsilon_zz]`. It is the
  clearest compact RCWA cross-check for that restricted scope, not for full
  longitudinal coupling.
- `EMpy/EMpy/RCWA.py` (`AnisotropicRCWA`, beginning approximately line 455)
  contains all nine tensor-component convolution matrices, but is 1D-only
  and mixes a substantially different interface/mode formulation. Use it as
  a cross-check for a narrowly matched 1D case, not as a near-verbatim source
  for this project's 2D S-matrix implementation.
- `Rigorous-Coupled-Wave-Analysis/TMM_functions/anisotropic.py` presents a
  dense 4x4 full-tensor TMM sketch, but references undefined variables and
  has no trustworthy test coverage; it is not an implementation source or
  oracle.

Decision: derive each Phase-6 substep independently only after its exact
scope and benchmark are defined; use S4 and RCWA.jl as cross-checks for the
five-component in-plane scope. Do not claim full 3x3 longitudinal-coupling
support until a citable formulation and an independent benchmark are both
available.

**Target 1.5 bounded literature search (2026-08-03, before implementing
Category 1 target 1.5)**: searched for a citable, independently-benchmarkable
full-3x3 (longitudinal-coupling) RCWA formulation beyond the vendored-repo
audit above. General-anisotropic-RCWA literature exists (Glytsis & Gaylord,
"Rigorous three-dimensional coupled-wave diffraction analysis of single and
cascaded anisotropic gratings," JOSA A 4, 2061-2080 (1987); a gyrotropic/
bi-anisotropic RCWA formulation referenced via a University of Arizona PhD
thesis abstract, "Rigorous Coupled Wave Analysis for Gyrotropic Materials"),
but none were both fetchable as readable full text in this environment
(JOSA A is paywalled; an arXiv candidate, 2510.01214, on birefringent
holographic gratings returned only undecodable binary PDF content via
`WebFetch`) and independently benchmarkable (no second, structurally
different source located to cross-check against, unlike targets 1.3/1.4's
S4 + RCWA.jl pairing). Target 1.5 remains explicitly deferred, not
implemented — see `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`'s 1.5 entry.

**Category 3 targets 3.4/3.5 (FFF/NVM) feasibility investigation
(2026-08-04, before deciding implement/defer)**: investigated whether Fast
Fourier Factorization (Popov & Nevière 2001) or the Normal Vector Method
(Lalanne 1997) could be implemented, given `tests/test_fourier_convergence.py`
(target 3.3) measured a real high-contrast 2D convergence weakness in the
current ordinary-Laurent's-rule solver. Both papers' bibliographic details
were confirmed via `WebSearch` this session (title/author/year/journal/
volume/pages all matched independently), but both are paywalled JOSA A
articles — no full text/equations were fetchable in this environment, same
situation as target 1.5's bounded literature search. `../S4` was read in
full for its own implementation of this technique family instead:
`S4/S4.h:49-71` (`use_polarization_basis`/`use_jones_vector_basis`/
`use_normal_vector_basis`/`use_normal_vector_field` options), dispatching
via `S4.cpp:1905-1930` to `fmm/fmm_PolBasisNV.cpp` (266 lines),
`fmm/fmm_PolBasisJones.cpp` (378 lines), `fmm/fmm_PolBasisVL.cpp` (274
lines) — all three built on `fmm/fmm_FFT.cpp` (239 lines), a
discretized/FFT-based permittivity representation, not the analytic
closed-form path (`fmm_closed.cpp`) already transcribed into this project.
This is a materially different Fourier-factorization architecture, in
direct tension with **ADR-002** (analytic shape Fourier transforms,
raster+FFT explicitly rejected for a different reason). Decision: defer
both (Category 3 targets 3.4/3.5), full account and revisit conditions in
`decisions.md` ADR-012.

**Category 4 targets 4.3-4.5 (Ellipse/Polygon primitives, 2026-08-04)**:
`S4/S4/pattern/pattern.c::pattern_get_fourier_transform` (lines 889-1032,
the same function `Circle`/`Rectangle`'s existing citations already
reference) turned out to already implement `ELLIPSE` (lines 955-964) and
`POLYGON` (lines 974-1008) analytically -- both closed-form, no raster/FFT
-- discovered while investigating target 4.4's design decision, not
assumed going in (the working assumption going in was that a polygon
would need raster+FFT, per the general FFF/NVM-family reasoning in
Category 3's ADR-012; that assumption was wrong for a single shape's own
boundary Fourier transform, a different and simpler question than the
eigenoperator-level vectorial factorization ADR-012 was actually about).
`geometry.Ellipse`/`geometry.Polygon` transcribe these two cases directly;
`Polygon.contains` also transcribes S4's PNPoly point-in-polygon test
(`pattern.c:180-193`). `Polygon.signed_distance_normal` is **not**
transcribed from S4's own `shape_get_normal` `POLYGON` case
(`pattern.c:256-281`, which selects the farthest, not nearest, boundary
segment) -- independently derived instead, since it contradicts this
project's own `Shape.signed_distance_normal` contract; see
`geometry.Polygon`'s docstring for the full account. See `decisions.md`
ADR-013 for the full accuracy-contract/ADR-005-revisit decision.

**Category 5 targets 5.5/5.6 bounded literature search (2026-08-04)**: after
transcribing Rakic's Lorentz-Drude model and its published Au/Ag/Al/Ti
coefficients (`rakic.jl`, cited above), attempted to additionally
cross-check computed gold optical constants against Johnson & Christy
(1972)'s original tabulated `n`,`k` data (a second, independent published
source, per this project's usual two-source preference) via `WebSearch`/
`WebFetch`. The paper's own bibliographic details (*Phys. Rev. B* 6,
4370-4379 (1972)) were confirmed, and `refractiveindex.info`'s hosted copy
of that dataset was located, but its actual numeric table was not
fetchable in this environment (`WebFetch` on both the interactive page and
a direct CSV-export URL returned no usable tabulated values — an
interactive/JS-rendered page and an HTTP 410, respectively). Rakic's own
published, peer-reviewed coefficients (`Appl. Opt.` 37, 5271-5283 (1998))
already satisfy Category 5 target 5.5's "published or tabulated reference
curve" requirement on their own; this was an attempt at a *second*,
independent cross-check (the extra rigor this project applied for e.g.
Phase 4a's eigenoperator oracle), not a gap in the primary validation. See
`tests/test_dispersion_models.py::test_rakic_metals_are_passive_and_metallic_across_visible_nir`
for where this is recorded in code, not just here.

- **Phase 5 (tapered/sloped sidewalls, `staircase.py`)**: per the
  `phase-reference-picker` skill's procedure, every RCWA-family repo under
  `REFERENCE/` (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl`,
  `Rigorous-Coupled-Wave-Analysis`) was grepped for "stair"/"taper" before
  writing any code; no matching staircase-generator implementation exists
  in any of them (the only hits were unrelated `meep`/`gprMax` doc pages —
  a different numerical method, not a formula source per the skill's own
  guidance). `staircase.py` is therefore independently derived: the
  staircase/multi-slice discretization technique itself is standard and
  already decided in `decisions.md` ADR-004, but there is no specific
  file/line citation to record here. Correctness rests on
  convergence-vs-`num_slices` evidence (`tests/test_staircase.py`), not an
  oracle-comparison test — consistent with `phases.md` Phase 5 having no
  new Fourier/eigenmode formula to validate against a source in the first
  place.

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

- **Moharam, M. G., Grann, E. B., Pommet, D. A., & Gaylord, T. K. (1995)**,
  "Formulation for stable and efficient implementation of the rigorous
  coupled-wave analysis of binary gratings," *J. Opt. Soc. Am. A* 12(5) —
  Phase 3's source formulation, per both vendored oracle scripts'
  (`1D_Grating_Gaylord_TE.py`/`TM.py`) own header-comment citation. This
  project did **not** transcribe numbers from the paper's tables directly
  (paper not available in this environment); per `references.md`'s own
  instruction and the precedent of running the vendored oracle rather than
  hand-copying, `tests/oracles/rcwa_1d_gaylord.py` transcribes the
  *algorithm* those scripts implement and is run directly to produce
  comparison numbers — see that module's docstring for the exact
  line-range citations and the one transcription fix applied.
- **Li, Lifeng (1996)**, "Use of Fourier series in the analysis of
  discontinuous periodic structures," *J. Opt. Soc. Am. A* — the
  foundational reference for the direct-vs-inverse-rule Fourier
  factorization distinction underlying Phase 2's `epsilon_hat` vs.
  `epsilon_inv_hat` Toeplitz construction (see `design.md`, Algorithm 3).
  Cite the specific rule/equation when Phase 2's docstrings are written,
  not just this general reference.
- **Lalanne, Philippe (1997)**, "Improved formulation of the coupled-wave
  method for two-dimensional gratings," *J. Opt. Soc. Am. A* 14(7),
  1592-1598 — the foundational 2D normal-vector-method (NVM) paper, cited
  by Category 3 target 3.5's feasibility decision (`decisions.md` ADR-012).
  Bibliographic details confirmed via `WebSearch`, not read as full text in
  this environment (paywalled) — no formula transcribed from it.
- **Popov, Evgeny, & Nevière, Michel (2001)**, "Maxwell equations in
  Fourier space: fast-converging formulation for diffraction by arbitrary
  shaped, periodic, anisotropic media," *J. Opt. Soc. Am. A* 18(11),
  2886-2894 — the foundational Fast Fourier Factorization (FFF) paper,
  cited by Category 3 target 3.4's feasibility decision (`decisions.md`
  ADR-012). Same paywalled/not-transcribed status as the Lalanne (1997)
  entry above.
- **Rakić, A. D., Djurišić, A. B., Elazar, J. M., & Majewski, M. L.
  (1998)**, "Optical properties of metallic films for vertical-cavity
  optoelectronic devices," *Appl. Opt.* 37, 5271-5283 — source of the
  Lorentz-Drude (LD) metal model and its published Au/Ag/Al/Ti
  coefficients, transcribed via the vendored
  `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl` (Category 5
  targets 5.5/5.6, `materials.Material.from_drude`/`from_drude_lorentz`,
  `RAKIC_GOLD` etc.) — the paper itself was not independently fetched in
  this environment; the vendored `.jl` file's own transcription (with its
  own citation to this paper) was the source actually read, per
  `rules.md`'s "S4/EMpy/RCWA.jl cross-check" allowance.
- **Sellmeier, Wilhelm (1871)** and **Cauchy, Augustin-Louis (1836)** — the
  two foundational, standard dispersion-formula forms implemented by
  Category 5 targets 5.2/5.3 (`materials.Material.from_sellmeier`/
  `from_cauchy`); no specific paper transcribed (both are common-knowledge
  19th-century results reproduced in essentially every optics textbook and
  in the vendored `EMpy/EMpy/materials.py`, the actual transcription source
  — see the `EMpy` row above). BK7's specific Sellmeier coefficients and
  its independently-published `n_d = 1.5168` (Fraunhofer d-line, 587.56 nm)
  were both confirmed via `WebSearch` this session against SCHOTT/
  refractiveindex.info sources, not transcribed from memory.

## How to Add a Reference

When a new phase's implementation cites a source (per `rules.md`'s
mandatory Documentation Standards), add it here too, with enough detail
(file + line range, or author/year/equation) that a future session can
re-locate it without re-searching from scratch.

## Choosing a Reference for a New Phase

Before writing physics code for a new phase, use the
`phase-reference-picker` skill (`.claude/skills/phase-reference-picker/`,
at the `Solver_own` workspace root). It exists because it's easy to default
to S4 out of habit (it's the Phase 1/2 oracle already in the codebase) even
when a different vendored repo is actually the better fit for a given
sub-task — the skill forces a real per-sub-task comparison across all
plausibly-relevant repos in `REFERENCE/`, plus an explicit decision between
transcribing a source near-verbatim (when a subtle sign/normalization
convention makes re-derivation risky, as with the S-matrix star product or
the Toeplitz subtraction rule) versus deriving independently and using the
repo only as a cross-check oracle (as was done for `tests/oracles/empy_tmm.py`,
where EMpy was the reference but its own bugs were fixed rather than
copied). Whatever the skill concludes should still be recorded here, in
this file, per "How to Add a Reference" above — the skill produces the
decision, this file is where it's durably logged.
