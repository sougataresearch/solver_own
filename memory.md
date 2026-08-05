# Project Memory — sougata_solver

Living document for future sessions (AI or human). Update this at the end
of every substantive session — see `rules.md`'s AI Coding Rules, item 6.

## Current Project Status

As of 2026-08-05 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7, targets
7.1-7.6), 2026-08-05 (Phase 7 / Category 9, targets 9.1-9.8), 2026-08-04
(Category 6, targets 6.1-6.6), 2026-08-04
(Category 5, targets 5.1-5.8), 2026-08-04 (Category 4, targets 4.1-4.7),
2026-08-04 (Category 3, targets 3.1-3.6), 2026-08-04 (Category 2, targets
2.1-2.5), 2026-08-03 (Phase 6, target 1.3), 2026-07-24 (Phase 5), 2026-07-23
(Phase 4b), 2026-07-21 (Phase 4a) and earlier entries below:
- **Category 7 (Layer handling), targets 7.1-7.6, are all complete.**
  7.1: `layer._require_valid_thickness`, called from `Layer.__post_init__`
  -- rejects NaN/non-positive thickness at construction, explicitly allows
  `math.inf` (the documented semi-infinite half-space sentinel; a direct
  test confirms `SMatrixStack` never actually reads that value, since
  `propagation_smatrix` is only called for interior layers).
  `tests/test_layer_validation.py` (15 tests). 7.2: equivalent repeated-
  layer representations (a thick layer split into N thin ones; a pattern
  reused by object identity vs. N separately-constructed structurally-
  equal patterns) give identical R/T -- `tests/test_layer_repetition.py`
  (7 tests), the invariant target 7.4's cache leans on. 7.3/7.4: an
  instance-scoped Toeplitz-matrix cache on `Simulation`
  (`_cached_toeplitz`/`_cached_toeplitz_component`), deliberately gated on
  a measured timing case per `rules.md`'s Performance Requirements
  exception clause rather than implemented speculatively (`design.md`'s
  "Layer/Toeplitz Caching Design", `decisions.md` ADR-016). **Honest
  correction made mid-session**: a first measurement (repeated identical
  patterned layers within one `solve()` call) wrongly attributed the
  entire "extra time per repeat" to Toeplitz reconstruction; isolating the
  two costs directly showed the (out-of-scope, still-uncached) eigensolve
  actually dominates at high `num_orders`, so that scenario only sees a
  ~4% real speedup. The scenario that actually justifies the cache,
  measured properly: a fixed-wavelength angle sweep (Category 8 target
  8.3, planned) reuses the same Toeplitz matrix across every sweep point
  (Toeplitz depends on pattern+wavelength, not angle) -- ~30% wall-clock
  reduction over a 20-point sweep, confirmed directly.
  `tests/test_layer_cache.py` (4 tests): equivalence to forced-uncached
  recomputation, a call-counting cache-hit check, and a direct angle-sweep
  cache-reuse regression. 7.5/7.6: `SimulationResult.layer_absorption()`
  (plus a new `SimulationResult.thicknesses` field to support it) -- per-
  layer absorbed power as a z-Poynting-flux-divergence combination of
  already-validated Category 9/Phase 7 pieces (`interior_amplitudes`,
  `propagate_amplitudes`, `z_poynting_flux`), deliberately not a new
  volumetric `Im(eps)*|E|^2` formula (`design.md`'s "Layer-Wise Absorption
  Design", `decisions.md` ADR-017). Validated by the `R+T+sum(A)=1`
  energy-balance identity itself -- finally closing the gap
  `tests/test_stress_regression.py`'s docstring had flagged since Category
  2 ("layer-wise absorption isn't implemented yet"), reusing that same
  file's already-vetted `eps=-396+80j` lossy fixture.
  `tests/test_layer_absorption.py` (4 tests). **Second honest finding,
  found while validating 7.6, documented rather than silently avoided**:
  `layer_absorption()` inherits `interior_amplitudes`/`propagate_amplitudes`'s
  existing numerical-stability envelope -- a thick, highly lossy, high-
  `num_orders` case (`max(Im(q))*thickness ~= 38`) numerically overflows
  the deepest evanescent modes' backward-propagated amplitude, giving a
  nonsensical `layer_absorption() ~= 573`; the same fixture at a thinner
  layer (`~6.3`) satisfies the energy identity to `~1e-6`. Not fixed (no
  formula change, matching the same-class transfer-matrix-vs-S-matrix
  numerical-stability tradeoff this project's own S-matrix choice was
  originally made to avoid at the R/T level) -- documented in
  `troubleshooting.md` and `decisions.md` ADR-017, with a dedicated
  regression test on the failure symptom itself so a future silent
  workaround doesn't land unnoticed. 542 tests pass project-wide (512 at
  the start of this category: 503 fast + 9 slow -- 533 fast + 9 slow now,
  30 new fast tests: 15+7+4+4 across the four new test files), full
  fast+slow suite re-run and confirmed green.
- **Phase 7 (Real-Space Field Reconstruction & Visualization), tracked at
  atomic-target grain as Category 9 (targets 9.1-9.8), is complete.**
  9.1: `fields.modal_field_components` (per-order `Ex,Ey,Ez,Hx,Hy,Hz` from
  modal amplitudes) and `fields.reconstruct_field_at_points` (inverse-
  Fourier phase sum onto real-space points/grids), citing S4's
  `GetInPlaneFieldVector` (`S4.cpp:1959-1995`) and `GetFieldAtPoint`
  (`S4.cpp:1997-2074`) for the transverse+longitudinal formulas
  (`Ez = epsilon_inv @ (ky*Hx - kx*Hy) / omega`,
  `Hz = (kx*Ey - ky*Ex) / omega`). 9.2/9.3: `fields.propagate_amplitudes`
  (depth-dependence ansatz `a(z)=a_top*exp(+i*q*z)`,
  `b(z)=b_top*exp(-i*q*z)`) and `smatrix.interior_amplitudes` (interior-
  layer mode-amplitude recovery from `SMatrixStack.partial_smatrix_up_to`
  plus the already-known `a0`/`b_reflected`) — **both independently
  derived, not transcribed** (`decisions.md` ADR-015): S4 uses a
  structurally different block-tridiagonal `SolveInterior` algorithm; this
  project instead reuses its own already-implemented
  `partial_smatrix_up_to` architecture via standard Redheffer star-product
  algebra (`b_i = inv(S11) @ (b_reflected - S10 @ a0)`,
  `a_i = S00 @ a0 + S01 @ b_i`), validated by a zero-free-parameter
  self-consistency check (recovering amplitudes at the full-stack partial
  matrix exactly reproduces the already-known transmitted amplitude).
  **Real finding, not a bug**: this project's established
  `fields.z_poynting_flux` modal quadratic form is exactly **2x** the
  textbook real-space flux `Sz = 0.5*Re(Ex*conj(Hy) - Ey*conj(Hx))` —
  confirmed directly (single-order uniform-layer case: `z_poynting_flux`
  gives `1.0`, textbook formula on the same mode gives `0.5`) — harmless
  for `R`/`T` (a ratio, so the factor cancels) but must be accounted for
  (`Sz = Re(...)`, no `0.5`) when computing absolute real-space flux from
  raw reconstructed fields; documented in `CONVENTIONS.md` and
  `troubleshooting.md`, not silently absorbed into the new code. 9.4-9.6:
  validated via analytic-plane-wave match, transversality (`k.E=0`,
  `k.H=0`), field continuity across a real material interface, 1D
  periodicity, and flux-matches-`R`/`T` (using the corrected no-`0.5`
  formula, agreeing to the solver's own `T` at multiple grid resolutions)
  — `tests/test_field_reconstruction.py` (10 new tests). 9.7/9.8:
  `fields.save_field_grid_npz` (raw-data-only, per ADR-009/010) plus two
  `structures/` example scripts (`structures/trench/trench_field_cross_section.py`,
  a 1D lamellar-grating (x,z) cross-section, saved via direct `np.savez`
  since `save_field_grid_npz`'s single-z-plane signature doesn't fit a
  z-sweep; `structures/via/pillar_field_cross_section.py`, a 2D circular-
  pillar (x,y) field map at one depth, the genuine single-plane case
  `save_field_grid_npz` was designed for) and one `postprocessing/` script
  (`postprocessing/plot_field_cross_section.py`, auto-detects either
  `.npz` layout, plots `|E|^2` via `pcolormesh`) — both example scripts
  run end-to-end (`R+T=1.0000`) and their output PNGs visually inspected,
  showing physically sensible interference/near-field patterns. 512 tests
  pass project-wide (502 at the start of this category: 493 fast + 9 slow
  -- 503 fast + 9 slow now), full fast+slow suite re-run and confirmed
  green.
- **Category 6 (Boundary conditions and excitation), targets 6.1-6.6, are
  all complete.** 6.1: confirmed `CONVENTIONS.md`'s existing "Phasor and
  propagation convention"/"Polarization convention" sections already
  satisfy this target, then added a "Worked polarization examples" table
  (TE/TM/linear/RCP/LCP/elliptical amplitude pairs). 6.2/6.3: two strong,
  oracle-independent symmetry invariants, both verified numerically before
  being encoded as tests -- (a) at normal incidence, an isotropic stack's
  `R`/`T` are *identical* across every polarization state at fixed total
  power (full rotational symmetry, stronger than energy conservation
  alone); (b) at fixed `theta`/polarization, `R`/`T` are independent of
  azimuth `phi` for a laterally-uniform stack. Combined with an 84-case
  oblique-incidence energy-conservation sweep (7 states x 4 azimuths x 3
  angles), `tests/test_polarization_states.py` (89 tests). 6.4: 
  characterized the grazing-incidence boundary directly rather than
  assuming it -- any `theta<90 deg` is supported (finite, energy-
  conserving to `1e-8`+ up to `89.999 deg`); exactly `theta=90 deg` raises
  a plain `ValueError` (from `scipy.linalg.lu_factor`'s finiteness check),
  traced to a genuine floating-point coincidence
  (`math.sin(math.radians(90.0)) == 1.0` exactly in float64, making the
  incidence half-space's `q` exactly `0.0` for `n=1`, not merely small).
  Added to `design.md`'s Failure Contract, which also picked up several
  Category 4/5 rows that had never been backfilled into it.
  `tests/test_grazing_incidence.py` (9 tests). 6.5: found already satisfied
  by `tests/test_mode_classification.py` (Category 1 target 1.8, normal
  incidence) -- added the one thing that coverage didn't have and Category
  6 is specifically about: an oblique-incidence Rayleigh-threshold case,
  where the `+m`/`-m` order degeneracy breaks and the two orders cross
  threshold at different wavelengths (confirmed by a coarse sweep, not a
  closed form -- the oblique threshold condition has no simple analytic
  form the way the normal-incidence one does).
  `tests/test_oblique_rayleigh_threshold.py` (8 tests). 6.6: **found a
  better answer than "defer, not needed"** -- bottom (reverse-side)
  illumination is already achievable via the existing `Simulation`
  constructor (reverse the layer list, swap `incidence`/`transmission`;
  `Layer.thickness`/`pattern` carry no inherent z-direction), verified via
  the Stokes transmittance-reciprocity relation for a lossless reciprocal
  medium (`T` matches to `~1e-15` at normal incidence between the forward
  and reversed simulations) rather than just asserted. No new `Simulation`
  parameter added, per `rules.md`'s "don't add features not needed" --
  `decisions.md` ADR-014 records the decision, `tests/test_bottom_incidence.py`
  (3 tests) is the permanent regression guard, including an honest
  counter-check that `R` genuinely differs between directions at oblique
  incidence (not claimed direction-independent, unlike `T`). 502 tests
  pass project-wide (393 at the start of this category: 384 fast + 9 slow
  -- 493 fast + 9 slow now), full fast+slow suite re-run and confirmed
  green.
- **Category 5 (Material models), targets 5.1-5.8, are all complete.**
  5.1: `Material.__init__` and `Material.epsilon_tensor` both now validate
  (finite values, correct tensor shape) -- construction-time alone can't
  catch a dispersion callable that only misbehaves away from the probe
  wavelength (e.g. outside an interpolation table's domain), so
  `epsilon_tensor` re-validates every call. 5.2/5.3: `Material.from_sellmeier`/
  `from_cauchy`, transcribed from the vendored `EMpy/EMpy/materials.py`
  (a genuinely useful, previously-unused-for-this-purpose vendored file) --
  validated against BK7's published Sellmeier coefficients and its
  independently-published `n_d=1.5168` (both confirmed via `WebSearch`,
  computed value matched to 5 significant figures). 5.4/5.5/5.6:
  `Material.from_lorentz`/`from_drude`/`from_drude_lorentz`, all
  transcribed from `RigorousCoupledWaveAnalysis.jl/src/BasicMaterials/rakic.jl`
  -- a genuine, real published metal-dispersion model (A. D. Rakić et al.,
  Appl. Opt. 37, 5271-5283 (1998)) found in a vendored repo, including its
  published Au/Ag/Al/Ti coefficient tables, transcribed verbatim as
  `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM`. **Important
  finding, independently re-derived and tested, not assumed**: target 5.4
  explicitly required a causality/sign-convention check -- under this
  project's `d/dt->-i*omega` phasor convention (`CONVENTIONS.md`), a
  passive Lorentz oscillator must have `Im(eps)>0` at resonance; hand-
  derived the exact value (`i*strength/(gamma*omega0)`) and confirmed it
  matches Rakic's Julia sign (`-1im*ω*o.Γ`), the same class of check that
  Category 2 target 2.5 got wrong on a first attempt (a naively-reused
  `n=-20+2j` index that turned out to be a gain, not lossy, medium under
  this project's convention) -- this time verified correct *before*
  shipping, not caught after. The Drude formula was additionally cross-
  checked between **two independently vendored sources**
  (`Rigorous-Coupled-Wave-Analysis/TMM_examples/TMM_Drude.py` and
  `rakic.jl`), confirmed algebraically identical before trusting either.
  A bounded `WebSearch`/`WebFetch` attempt to also cross-check against
  Johnson & Christy (1972)'s raw tabulated Au n,k data (a second,
  independent published source) found the paper's bibliographic details
  but could not fetch the actual data table in this environment (same
  class of limitation as Category 1 target 1.5's and Category 3's bounded
  searches) -- Rakic's own published coefficients already satisfy the
  target's requirement on their own; not silently skipped, recorded in
  `references.md`. 5.7: confirmed Category 1's tensor-solver gate (targets
  1.3/1.4/1.6) is already met, then closed the one previously-untested
  combination -- a genuinely *dispersive* tensor material (built from a
  Category 5 dispersion model) flowing through Category 1's uniform-
  diagonal and patterned-anisotropic eigensolvers end to end. 5.8:
  `Material.source` (optional provenance/citation string), threaded
  through every `from_*` classmethod and `geometry_io`'s JSON schema.
  **Real bug found and fixed while validating "serialized output"**:
  `output_paths.write_run_metadata` used the platform-default text
  encoding (`cp1252` on Windows), which raised `UnicodeEncodeError` the
  moment a real non-ASCII citation ("Rakić et al.") was written through it
  for the first time -- fixed to explicit UTF-8, the first time any test
  in this project exercised that code path with non-ASCII content. 393
  tests pass project-wide (336 at the start of this category: 327 fast + 9
  slow -- 384 fast + 9 slow now), full fast+slow suite re-run and
  confirmed green.
- **Category 4 (Geometry engine), targets 4.1-4.7, are all complete.**
  4.1: `geometry._require_finite`/`_require_positive`, called from
  `Lattice`/`Lattice1D.__init__` and `Circle`/`Rectangle`/`Slab.__post_init__`
  -- non-finite dimensions, degenerate lattice vectors, and non-positive
  shape sizes now raise `ValueError` at construction instead of surfacing
  later as a NaN Fourier coefficient or a cryptic `LinAlgError`
  (`tests/test_geometry_validation.py`, 29 tests). 4.2:
  `geometry.validate_pattern_fits_lattice`, wired into `Simulation.__init__`
  -- documents (and verifies against a from-scratch periodic-tiling raster
  reference, not just argues) that a shape crossing a cell edge is already
  handled correctly by the existing reciprocal-lattice-point Fourier
  evaluation (a Poisson-summation consequence, no code change needed), and
  adds an explicit, conservative self-overlap-across-periodic-images
  rejection for the one case that genuinely isn't handled automatically
  (`tests/test_unit_cell_bounds.py`, 6 tests). 4.3/4.4/4.5: `geometry.Ellipse`
  and `geometry.Polygon`, both transcribed from
  `S4/S4/pattern/pattern.c::pattern_get_fourier_transform` (lines 889-1032)
  -- a genuine surprise found while investigating target 4.4's "analytic vs.
  raster/FFF" design decision: the working assumption going in (informed by
  Category 3's ADR-012 FFF/NVM investigation) was that a polygon would need
  raster+FFT; reading S4's actual `POLYGON` case (lines 974-1008) showed a
  closed-form boundary/edge-sum formula instead, architecturally identical
  in kind to `Circle`/`Rectangle`'s existing analytic transforms, not the
  discretized/FFT machinery ADR-012 was about (`decisions.md` ADR-013
  records the resulting decision -- analytic, exact for any simple polygon
  -- and explicitly narrows `decisions.md` ADR-005's polygon deferral rather
  than reversing it; GDS/raster import remains out of scope). `Ellipse`
  reduces exactly to `Circle` when `hx==hy`; a square `Polygon` reduces
  exactly to `Rectangle` (both regression-verified, not just plausible).
  Both validated against from-scratch rasterized references
  (`tests/test_ellipse.py`, 19 tests; `tests/test_polygon.py`, 24 tests,
  including a non-convex L-shape) and each has a `structures/via/` example
  reaching R+T=1.0000 across a 21-point wavelength sweep
  (`elliptical_pillar.py`, `triangular_pillar.py`) -- satisfying the
  Category 4 exit criteria's "geometry-only tests and one end-to-end RCWA
  example" per new shape directly. `Polygon.signed_distance_normal` is
  **not** transcribed from S4's own `POLYGON` normal formula
  (`pattern.c:256-281`), which was found to select the *farthest*, not
  nearest, boundary segment -- contradicting this project's own
  `Shape.signed_distance_normal` contract already established by
  `Circle`/`Rectangle`/`Ellipse`; implemented independently instead (an
  elementary point-to-segment-distance calculation, not risky enough to
  need transcription per `rules.md` Documentation Standards option 2). 4.6:
  `geometry_io.py`, a minimal JSON `Pattern`-import format (`unit`/
  `background`/`shapes`, isotropic-scalar materials only,
  `Circle`/`Rectangle`/`Ellipse`/`Polygon`/`Slab`) -- `json` module only,
  never `eval`/`exec` (per `rules.md` Security Rules), every malformed-input
  case raises a clear `ValueError` naming the offending key/shape index.
  Deliberately not wired into `Simulation`/`Layer` construction, matching
  the target's own "before solver integration" wording -- the returned
  `Pattern` already works with the existing public API unmodified, verified
  by one end-to-end test (`tests/test_geometry_io.py`, 17 tests). 4.7:
  `staircase.slice_profile`, a shape-agnostic `pattern_at(frac) -> Pattern`
  generalization; the three existing shape-specific taper generators
  (`staircase_circle_layers`/`_rectangle_layers`/`_slab_layers`) were
  refactored (not reimplemented) into thin wrappers around it -- the full
  pre-existing `tests/test_staircase.py` suite, including Phase 5's
  zero-taper and energy-conservation regression tests, was re-run and
  confirmed passing unchanged after the refactor, not merely assumed safe
  (`tests/test_profile_slicing.py`, 6 new tests, including a non-linear
  non-taper profile proving genuine generality). 336 tests pass
  project-wide (226 fast + 9 slow at the start of this category's work,
  327 fast + 9 slow now -- 101 new fast tests, no new `slow` tests), full
  fast+slow suite re-run and confirmed green.
- **Category 3 (Fourier factorization), targets 3.1-3.6, are all complete.**
  3.1: a "Fourier-factorization rule inventory" table added to `design.md`
  (Algorithm 3a), tabulating every uniform/1D/2D solver branch's
  direct/inverse/numerical-inverse choice with citations, backed by
  `tests/test_fourier_factorization_rules.py` (6 tests, black-box via the
  already-public `LayerEigenmodes.epsilon_inv` field plus one white-box
  1D TE/TM-block discrimination test). **Finding**: `epsilon_inv_hat` (the
  separately-Fourier-factorized inverse-rule Toeplitz) is consumed as such
  in exactly one place project-wide — the 1D TM block; every 2D path
  (isotropic and anisotropic) uses a numerical matrix-inverse of the
  *direct*-rule Toeplitz instead, a different Fourier-factorization
  operation despite both being called "epsilon_inv" in code, now made
  explicit as a table row rather than left implicit across four docstrings.
  3.2/3.3: `tests/test_fourier_convergence.py`, two new `slow` fixtures —
  a high-contrast (`n=10`) 1D lamellar grating and a high-contrast (`n=5`)
  2D pillar, each with *actually measured* (not assumed) reflectance vs.
  harmonic-order-count tables recorded in the test docstrings. **Two
  honest findings, both real, neither hidden**: the 1D fixture is not
  monotonically convergent from its very first data point (a genuine
  pre-asymptotic transient) and has not fully converged even at
  `num_ord=320` (~6% relative error remaining at `num_ord=160`); the 2D
  fixture is far more dramatic — `num_orders=25` gives `R=0.214`, an
  order-of-magnitude non-monotonic outlier against its low-order
  neighbors and the ~0.0236 eventually-converged value, a direct,
  measured illustration of ordinary Laurent's-rule 2D Fourier
  factorization's known weakness at sharp discontinuities (already
  documented in `solve_layer_eigenmodes_patterned`'s docstring, now backed
  by an actual number instead of just a citation). 3.4/3.5: Fast Fourier
  Factorization (Popov & Nevière 2001) and the Normal Vector Method
  (Lalanne 1997) feasibility decisions — both evaluated and **explicitly
  deferred** (`decisions.md` ADR-012), directly motivated by 3.2/3.3's
  measured convergence weakness (this wasn't a defer-on-faith decision;
  the technique's applicability was concretely demonstrated first). Both
  papers' bibliographic details were confirmed via `WebSearch` this
  session but neither paper's full text/equations were fetchable
  (paywalled JOSA A, same situation Category 1 target 1.5 already hit) —
  per `rules.md` AI Coding Rule 1, nothing was transcribed from either.
  `../REFERENCE/S4` was instead read in full for its own implementation
  of this technique family: `S4.h:49-71`'s `use_polarization_basis`/
  `use_normal_vector_basis`/`use_normal_vector_field` options dispatch
  (`S4.cpp:1905-1930`) to three dedicated files
  (`fmm/fmm_PolBasisNV.cpp`, `fmm_PolBasisJones.cpp`, `fmm_PolBasisVL.cpp`
  — ~900 combined lines), all built on `fmm_FFT.cpp`'s
  discretized/rasterized permittivity representation, **not** the
  analytic closed-form path this project already transcribes
  (`fmm_closed.cpp`). This is a materially different Fourier-factorization
  architecture, confirmed by reading the actual dispatch code rather than
  inferring from option names — and it directly conflicts with the
  already-shipped **ADR-002** ("analytic shape Fourier transforms over
  raster+FFT," explicitly rejecting raster+FFT for a different reason,
  pixelization error at smooth boundaries). 3.6 (selected improvement)
  therefore has nothing approved to implement — recorded as its own
  explicit outcome per the atomic-targets register's own "explicitly
  decide implement/defer" allowance, not silently skipped; zero solver
  code changed, zero regression risk. 232 tests pass project-wide (227 at
  the start of this session: 220 fast + 7 slow -- 226 fast + 9 slow now),
  full fast+slow suite re-run and confirmed green.
- **Category 2 (Numerical methods), targets 2.1-2.5, are all complete.**
  2.1: a "Failure Contract" section added to `design.md` (four tables --
  `ValueError`, `NotImplementedError`, `LinAlgError`, `WARNING` -- built by
  grepping every `raise`/`logger.warning` site in `src/sougata_solver/`,
  not from memory), backed by `tests/test_failure_contract.py` (17 tests).
  2.2: `layer.EigenmodeDiagnostics` (`cond_epsilon`, `cond_phi`,
  `min_eigenvalue_gap`, `num_propagating`, `num_evanescent`), attached as a
  new optional `LayerEigenmodes.diagnostics` field by every `eigenmodes.py`
  solver, reusing already-computed condition numbers where available --
  `tests/test_eigenvalue_diagnostics.py` (6 tests) confirms the fields
  match independent recomputation and that no other `LayerEigenmodes`
  field changed. 2.3/2.4: `eigenmodes.DEGENERATE_GAP_THRESHOLD` (`1e-6`,
  monkeypatch-configurable like `ILL_CONDITIONED_THRESHOLD`) and
  `_warn_on_small_eigenvalue_gap`, plus reuse of target 1.7's
  `_canonical_mode_order`, applied to the three anisotropic dense solvers
  (`solve_layer_eigenmodes_uniform_diagonal`/`_inplane`,
  `solve_layer_eigenmodes_patterned_inplane`) --
  `tests/test_sweep_mode_matching.py` (4 tests) and
  `tests/test_degeneracy_warning.py` (5 tests). **Two honest findings,
  both from actually running the change against the existing suite rather
  than assuming it would pass**: (a) extending `_canonical_mode_order` to
  `solve_layer_eigenmodes_patterned` (Phase 4a, isotropic 2D) was tried
  first and broke `tests/test_2d_pillar.py`'s TE/TM-block regression tests
  (that solver's *natural*, un-sorted `eig()` output happens to keep an
  exact block structure at `ky=0` that re-sorting destroys) -- reverted
  per `rules.md` AI Coding Rule 3, documented in that function's own
  docstring rather than silently dropped; (b) the same investigation found
  an ordinary circular-pillar case has a genuinely near-zero eigenvalue gap
  from routine lattice `C4v` symmetry, so the gap-warning was deliberately
  **not** applied to that solver either (it would misfire on harmless,
  expected degeneracy). 2.5: `tests/test_stress_regression.py` (2 tests),
  a lossy high-contrast fixture run through the **full**
  `Simulation.solve()` pipeline (Phase 4b's own stress sweep,
  `tests/test_2d_pillar_stress.py`, only ever cross-checked eigenvalues,
  never called `solve()`). **Third finding, a sign-convention trap, not a
  solver bug**: the first attempt reused Phase 4b's `n=-20+2j`
  "lossy-metal-like" index verbatim and got `R+T` up to ~17 through the
  full pipeline -- `n=-20+2j` squares to `Im(eps)<0`, which is a **gain**
  medium under this project's documented `d/dt -> -i*omega` phasor
  convention (`CONVENTIONS.md`), not a lossy one, so `R+T>1` was the
  numerically-correct answer for that (mislabeled) input; Phase 4b's own
  eigenvalue-only test never caught this because it never computed R/T.
  Fixed by using a correctly-signed lossy metal (`eps=-396+80j`) for the
  new fixture; Phase 4b's already-shipped, already-passing file was left
  untouched (Rule 3 again -- not a bug in that file's actual assertions,
  only an imprecise docstring label). Since layer-wise absorption
  (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7 targets 7.5/7.6) isn't
  implemented yet, the full `R+T+A=1` lossy energy identity can't be
  checked; the weaker, still-meaningful passivity check (`R>=0`, `T>=0`,
  `R+T<=1`) is used instead, honestly scoped to what's actually
  implemented. 227 tests pass project-wide (186 fast pre-session + 34 new
  fast + 7 unchanged `slow`), full fast+slow suite re-run and confirmed
  green after every change, not just once at the end.
- **Phase 6, target 1.3 (uniform diagonal-tensor anisotropic layers) is
  complete**, per `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 1's atomic
  sequencing (targets 1.4-1.8 remain open, tracked separately). New function
  `eigenmodes.solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx,
  eps_yy, eps_zz)`, transcribed from `S4/S4/S4.cpp:1889-1906` (the
  uniform-anisotropic-material branch, `0 != M->type`), specialized to the
  diagonal case (off-diagonal `abcde` in-plane components zeroed): `kp` is
  built from `1/eps_zz` (matching `Epsilon_inv` at `S4.cpp:1897`), `Epsilon2
  = block_diag(eps_xx*I_n, eps_yy*I_n)` (matching `S4.cpp:1898,1901`), and
  the general eigenoperator construction reuses `solve_layer_eigenmodes_patterned`'s
  already-cited `op = Epsilon2 @ kp - coupling` structure. Wired into
  `simulation.py`'s uniform-layer dispatch via `Material.is_diagonal`; the
  general (in-plane-coupled or longitudinally-coupled) case still raises
  `NotImplementedError` naming targets 1.4/1.5.
  **Non-obvious finding, caught by the test itself, not assumed**: the
  already-existing `Epsilon2` block-index convention (`CONVENTIONS.md`:
  `u = [-Ey; Ex]`, block 0 <-> `-Ey`, block 1 <-> `Ex`) means `eps_xx`
  (`Epsilon2`'s top-left block) governs the `Ey` field component
  (s-polarization at normal incidence), and `eps_yy` governs `Ex`
  (p-polarization) — the reverse of a naive "eps_xx acts on Ex" assumption.
  A first draft of the validation test had this backwards; caught
  immediately because the closed-form Fresnel-oracle comparison
  (`tests/test_anisotropic_uniform.py::test_normal_incidence_uniaxial_slab_matches_fresnel_oracle_per_axis`)
  failed with values that turned out to match the *other* axis exactly —
  fixed by swapping the assertion mapping, not the physics. Validated by:
  (a) closed-form normal-incidence birefringence formula
  (`q_x^2=eps_xx*omega^2`, `q_y^2=eps_yy*omega^2`) at the unit level; (b) the
  already-validated isotropic `tests/oracles/fresnel.py` oracle applied
  per-principal-axis to a full uniaxial-slab `Simulation` at normal
  incidence (a genuine external-oracle comparison, not a self-check, since
  `fresnel.py` has no anisotropy concept and was reused unmodified); (c)
  isotropic-tensor (`eps_xx=eps_yy=eps_zz`) reduction to
  `solve_layer_eigenmodes_uniform`'s result through a full `Simulation.solve`,
  at both normal and oblique incidence; (d) energy conservation for a
  genuinely birefringent slab at oblique incidence and mixed polarization.
  143 tests pass project-wide (123 prior + 20 new,
  `tests/test_anisotropic_uniform.py`). **Targets 1.4 (in-plane coupling),
  1.5 (longitudinal coupling — likely to be explicitly deferred, no citable
  formulation found by the Phase 6 reference audit in `references.md`), 1.6
  (patterned anisotropic layers), 1.7 (degeneracy policy), and 1.8 (mode
  classification) are next**, per the approved plan sequencing them one at
  a time.
- **Phase 6, target 1.4 (uniform in-plane-coupled anisotropic layers) is
  complete.** New function `eigenmodes.solve_layer_eigenmodes_uniform_inplane`
  generalizes target 1.3's diagonal solver by populating `Epsilon2`'s
  previously-zero off-diagonal quadrants with `eps_xy`/`eps_yx` (same
  `S4.cpp:1889-1906` citation, now using the full `abcde[0:8]` in-plane
  block). Wired into `simulation.py`'s uniform-layer dispatch; a
  longitudinal-coupling guard (nonzero `eps_xz/eps_yz/eps_zx/eps_zy`) still
  raises `NotImplementedError` naming target 1.5.
  **New independent oracle**: `tests/oracles/rcwa_anisotropic_inplane_jl.py`,
  hand-transcribed from `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:134-165`
  (`eigenmodes(...,l::AnisotropicLayer)`, the *uniform*-anisotropic-layer
  branch — a different `Common.jl` function from Phase 4a's patterned-layer
  oracle). **Non-obvious finding, determined empirically (not derivable
  from the S4 citation alone)**: reconciling this oracle against
  `solve_layer_eigenmodes_uniform_inplane`'s `q^2` required not just the
  already-known `k0`-normalization and overall sign flip (from Phase 4a's
  oracle), but a *third*, in-plane-coupling-specific convention: `kx`/`ky`
  swapped **and** `eps_xy`/`eps_yx` negated when calling the oracle. This
  was invisible in target 1.3's diagonal-only case (where it reduces to
  just needing `eps_xx`/`eps_yy` swapped, which is what a naive read of the
  block-index convention alone predicted) and only surfaced once
  off-diagonal terms were exercised — found by a brute-force search over
  swap/negate hypotheses (not guessed), confirmed to ~1e-13 across 20
  random-parameter trials plus the diagonal-only reduction case; see the
  oracle module's docstring for the full account. **Test-authoring mistake
  caught before shipping**: an initial oracle-comparison test used an
  absolute tolerance (`abs=1e-6`) sized for Phase 4a's oracle test, which
  runs at `wavelength=1.0` (so `omega~O(10)`, `q^2~O(100)`); this test runs
  at a realistic `wavelength=0.55e-6` (`omega~O(1e7)`, `q^2~O(1e14)`), where
  that same absolute tolerance is meaningless — caught by a spurious-looking
  failure with a genuinely tiny relative residual, fixed by switching to
  `rel=1e-8`, not by loosening the check. **Second test mistake, same
  session**: the first energy-conservation test used a non-Hermitian
  in-plane tensor (`eps_xy != conj(eps_yx)`), which is physically a
  gain/loss (non-power-reciprocal) medium — R+T=1 is not expected for that
  case, and the test failed at R+T~0.96 for a real physical reason, not a
  solver bug; fixed by using a Hermitian (lossless) tensor instead, per the
  test's updated docstring. 158 tests pass project-wide (143 prior + 15
  new, `tests/test_anisotropic_inplane.py`). **Target 1.5 (longitudinal
  coupling) is next** — per the existing Phase 6 reference audit
  (`references.md`), no vendored repo or literature source was found with a
  citable formulation for that scope; this session will do a bounded
  literature search before deciding to implement or explicitly defer, per
  the approved plan.
- **Phase 6, target 1.5 (longitudinal coupling) was evaluated and
  explicitly deferred, 2026-08-03** — not implemented. Per `rules.md` AI
  Coding Rule 1, a bounded literature search was done (`WebSearch`/`WebFetch`)
  before concluding, beyond the standing vendored-repo audit in
  `references.md`. Found that general-anisotropic-RCWA literature exists in
  principle (Glytsis & Gaylord 1987 JOSA A; a gyrotropic-RCWA PhD thesis)
  but nothing both readable in this environment (JOSA A paywalled; a
  candidate arXiv preprint, 2510.01214, returned only undecodable binary
  PDF content via `WebFetch`) and independently benchmarkable (no second
  structurally-different source found to cross-check against, unlike
  targets 1.3/1.4's S4+RCWA.jl pairing) — so no formula was written, per
  Rule 1's "say so explicitly, do not invent." `simulation.py`'s
  longitudinal-coupling guard (added with target 1.4) continues to raise
  `NotImplementedError` naming this target. See `references.md`'s "Target
  1.5 bounded literature search" entry for the full account. This is a
  "not found this session" conclusion, revisitable if a readable source
  becomes available later, not a permanent scope removal.
- **Phase 6, target 1.6 (patterned anisotropic layers) is complete.** New
  functions `fourier_factorization.pattern_epsilon_hat_component`/
  `toeplitz_matrix_component` (per-tensor-component direct-rule Toeplitz
  matrices, refactored to share the existing subtraction-rule accumulation
  with the scalar `pattern_epsilon_hat`/`toeplitz_matrix` via a new private
  `_pattern_fourier_sum`/`_toeplitz` helper, no behavior change to the
  existing scalar path — full suite re-run to confirm) and
  `eigenmodes.solve_layer_eigenmodes_patterned_inplane`. **New citation,
  not found during the original Phase 6 reference audit**: transcribed
  from `S4/S4/fmm/fmm_closed.cpp`'s `have_tensor` branch (lines 165-256),
  read in full this session — the earlier audit (`references.md`'s "Phase
  6 anisotropy reference audit") only looked at `S4.cpp`'s uniform-layer
  path and Common.jl/EMpy/Rigorous-Coupled-Wave-Analysis, never this
  branch of `fmm_closed.cpp` (whose isotropic sibling branch, lines 77-164,
  was already the Phase 4a citation). It confirms the natural
  generalization: `Epsilon2`'s four in-plane quadrants become full
  direct-rule Toeplitz matrices (same `xx/xy/yx/yy` block-index convention
  already established for the uniform case) and `kp`'s `Epsilon_inv`
  becomes the **numerical matrix inverse of the direct-rule `eps_zz`
  Toeplitz** (not a separately Fourier-factorized inverse-rule Toeplitz) —
  consistent with, and now directly citing, the same "matrix-inverse, not
  inverse-rule Toeplitz" pattern Phase 4a already established for the
  isotropic 2D case. Wired into `simulation.py`'s patterned-layer dispatch:
  isotropic pattern still uses the unmodified Phase 4a path; a
  diagonal/in-plane-tensor pattern uses the new path; any longitudinal
  component anywhere in the pattern (background or any shape) raises
  `NotImplementedError` naming target 1.5. Validated by: (a) reduction to
  Phase 4a's isotropic solver (off-diagonal Toeplitz matrices confirmed
  exactly zero, on-diagonal ones confirmed equal); (b) reduction to target
  1.4's uniform-tensor solver for a spatially-uniform "pattern" (a shape
  whose material equals the background); (c) energy conservation for a
  genuinely patterned, Hermitian (lossless) anisotropic pillar/rectangle
  case. No new test-authoring mistakes this target — the tolerance-scale
  and Hermitian-material lessons from target 1.4 were applied from the
  start. 170 tests pass project-wide (158 prior + 12 new,
  `tests/test_anisotropic_patterned.py`). **Targets 1.7 (degeneracy
  policy) and 1.8 (mode classification) are next.**
- **Phase 6, target 1.7 (degeneracy policy) is complete.** New helper
  `eigenmodes._canonical_mode_order(q, phi)` applies a documented,
  deterministic sort (rounded `Re(q)`, then `Im(q)`, then original
  `eig`-output index as the tie-break for exact/near-degenerate
  eigenvalues) to the three dense anisotropic eigensolvers introduced by
  targets 1.3/1.4/1.6 (`solve_layer_eigenmodes_uniform_diagonal`,
  `solve_layer_eigenmodes_uniform_inplane`,
  `solve_layer_eigenmodes_patterned_inplane`) — deliberately **not**
  applied to the pre-existing Phase 4a `solve_layer_eigenmodes_patterned`,
  to keep this target's blast radius scoped to the anisotropic solvers it
  actually concerns, per the approved plan. This is a policy/ordering
  layer, not a numerical fix: `numpy.linalg.eig`'s LAPACK `geev` backend
  was already confirmed deterministic for identical input (no randomness),
  verified directly rather than assumed
  (`tests/test_anisotropic_degeneracy.py::test_repeated_solve_is_deterministic`).
  Explicitly does **not** claim eigenvalue continuity across a *changing*
  input (e.g. a wavelength sweep crossing a degeneracy) — that remains
  separate, already-tracked future work (`tasks.md` Category 2 target 2.3,
  "Sweep mode matching"). Builds on, rather than replaces, Phase 4b's
  existing `ILL_CONDITIONED_THRESHOLD` `WARNING`-logging precedent
  (detection, not correction) for actual near-degeneracy. Validated by:
  sort-key unit tests (explicit tie-break-by-original-index case), repeated-
  solve bit-identical-output tests, and energy conservation for a
  deliberately near-isotropic (`eps_xx=2.2501, eps_yy=2.25`) patterned
  anisotropic case. 179 tests pass project-wide (170 prior + 9 new,
  `tests/test_anisotropic_degeneracy.py`). **Target 1.8 (mode
  classification) is next — the last of the six targets in this session's
  approved plan.**
- **Phase 6, target 1.8 (mode classification) is complete — the last of
  this session's six approved targets.** New `eigenmodes.classify_propagating(q)`
  (boolean array, reusing `_select_q_branch`'s own real/imaginary branch
  convention rather than a separate re-derivation) and
  `SimulationResult.order_classification()` (per-order
  propagating/evanescent, keyed like `diffraction_efficiencies()`, using
  `all_modes[0]`/`all_modes[-1]`'s `q`). Validated against the analytic
  Rayleigh-threshold wavelength (`lambda_threshold = n_trans*period/m` at
  normal incidence) for a diffraction order flipping classification
  exactly at the predicted point, energy conservation on both sides, and
  the already-established incidence-vs-transmission-medium-differ case
  (different Rayleigh thresholds on each side of an interface).
  **Honest finding, deliberately not fixed by this target**: evaluating
  exactly at the threshold wavelength (not just near it) produces `NaN`
  R/T — confirmed directly, not assumed: `q == 0` for that order at the
  exact crossing, and `smatrix.py::interface_smatrix`'s `kp @ phi /
  q[None, :]` divides by zero (`RuntimeWarning`s from `smatrix.py:75`
  observed). This is a genuine, pre-existing solver limitation at the
  exact Wood's-anomaly/Rayleigh singular point, not a bug this session's
  work introduced, and not something target 1.8 (which only adds
  *classification*, not threshold-crossing numerical handling) was scoped
  to fix — documented in `troubleshooting.md`'s Already-Solved-Gotchas-
  adjacent list and tied explicitly to `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
  Category 6 target 6.4 ("Grazing-incidence boundary test"), the actual
  home for defining supported near-threshold behavior. The test suite
  deliberately checks `0.999x`/`1.001x` of the threshold, not the exact
  point. 186 tests pass project-wide (179 prior + 7 new,
  `tests/test_mode_classification.py`).

  **Category 1 status (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`) as of
  2026-08-03: targets 1.1-1.4 and 1.6-1.8 are all done; target 1.5
  (longitudinal coupling) is evaluated and explicitly deferred** (bounded
  literature search found no source both readable in this environment and
  independently benchmarkable — see that target's own entry above). This
  closes the six-target plan approved at the start of this session
  (`C:\Users\sougata.bhunia\.claude\plans\polished-enchanting-hopcroft.md`).
  Session total: 65 new tests across five new test files
  (`tests/test_anisotropic_uniform.py`,
  `tests/test_anisotropic_inplane.py`,
  `tests/test_anisotropic_patterned.py`,
  `tests/test_anisotropic_degeneracy.py`,
  `tests/test_mode_classification.py`), two new oracle modules
  (`tests/oracles/rcwa_anisotropic_inplane_jl.py`, hand-transcribed from
  `RigorousCoupledWaveAnalysis.jl`), five new `eigenmodes.py` functions
  (`solve_layer_eigenmodes_uniform_diagonal`,
  `solve_layer_eigenmodes_uniform_inplane`,
  `solve_layer_eigenmodes_patterned_inplane`, `_canonical_mode_order`,
  `classify_propagating`), two new `fourier_factorization.py` functions
  (`pattern_epsilon_hat_component`, `toeplitz_matrix_component`), and one
  new `simulation.py` method (`SimulationResult.order_classification`).
  186 tests pass project-wide (123 at session start). No existing
  oracle-comparison test was weakened or removed to make a new one pass
  (`rules.md` AI Coding Rule 3) — confirmed by re-running the full
  pre-session test count after every target. **No commit was made this
  session** — the working tree already contained substantial prior
  uncommitted work (Phase 3-5 and other in-progress changes) when this
  session started, predating this session's changes, and committing was
  not requested by the user; the next session/commit should account for
  all of it together, not just this session's diff.
- **Phase 5 (tapered/sloped sidewalls, staircase discretization) is
  complete.** New module `src/sougata_solver/staircase.py` with three
  generator functions (`staircase_circle_layers`, `staircase_rectangle_layers`,
  `staircase_slab_layers`), each producing `num_slices` uniform-in-z
  `Layer`s whose shape size is linearly interpolated (at each slice's
  z-midpoint) between a `top` and `bottom` value, per `decisions.md`
  ADR-004. Per the `phase-reference-picker` skill: this phase has **no
  citable source** in any vendored `REFERENCE/` repo — a grep for
  "stair"/"taper" across the whole RCWA family (`S4`, `EMpy`,
  `RigorousCoupledWaveAnalysis.jl`, `Rigorous-Coupled-Wave-Analysis`) found
  nothing; `staircase.py` is independently derived (the technique itself is
  standard/well-precedented, per ADR-004), flagged per `rules.md` AI Coding
  Rule 1, and its correctness rests on convergence-vs-`num_slices` evidence
  rather than an external oracle-comparison test — consistent with
  `phases.md` Phase 5's "no new Fourier/eigenmode math" scoping (every
  slice reuses Phase 3/4a's already-oracle-validated per-layer solve
  unchanged). Validated: (a) zero-taper (`top==bottom`) staircase
  reproduces the already-validated single-layer uniform-cross-section
  result to `1e-10` for `Circle`, `Rectangle`, and `Slab`, regardless of
  `num_slices` (regression guard); (b) energy conservation holds for
  actually-tapered cases; (c) two `slow`-marked convergence-vs-`num_slices`
  studies (`tests/test_staircase.py`), sweeping `num_slices` from 1 to 32
  (via) / 64 (trench) — both show clean monotone-shrinking successive
  differences (via: R settles to ≈0.565 by N=16-32; trench: R settles to
  ≈0.248 by N=32-64, needing more slices to converge than the via case,
  consistent with its larger relative taper). `structures/via/tapered_via.py`
  and `structures/trench/tapered_trench.py` (per ADR-010, print/save only,
  no plotting) both run end-to-end with R+T=1.0000 at every `num_slices`
  value. 130 tests pass project-wide (118 prior + 12 new — 10 fast + 2
  `slow` convergence studies — `tests/test_staircase.py`).
  **Follow-up, 2026-08-01**: at user request, renamed the tapered
  trench/via scripts' geometry constants to FDTD-style `TCD`/`BCD`/`SPACING`
  (matching an equivalent Lumerical grating structure group the user
  shared) and added `structures/via/tapered_pillar.py` (`Rectangle` case of
  `staircase_rectangle_layers`, previously the only shape type without an
  example script) — no new physics, same staircase generators; all three
  re-verified with R+T=1.0000.
- **Phase 4a (2D-periodic patterned layers, well-conditioned case) is complete
  and validated.** `solve_layer_eigenmodes_patterned` (the dense `2n x 2n`
  general eigenoperator from `S4/S4/rcwa.cpp::SolveLayerEigensystem` lines
  794-827) is implemented in `eigenmodes.py`; `simulation.py` dispatches to it
  for any 2D-patterned `Layer`. **Correction made this session** to a first
  draft (written by a different agent) that copied `solve_layer_eigenmodes_1d`'s
  `Epsilon2` construction (`epsilon_hat` top-left, `inv(epsilon_inv_hat)`
  bottom-right) into the 2D solver — this looked plausible but was never
  actually checked against the relevant S4 source. Reading
  `S4/S4/fmm/fmm_closed.cpp:109-139` in full (specifically the branch after
  line 133, which nobody had read before) shows that formula is the **1D
  special case only** (`0==Lr[2]&&Lr[3]==0`, S4's own "1D proper FFF rule"
  comment). The true-2D, no-polarization-basis path S4 actually uses by
  default is plain **Laurent's rule throughout**: `Epsilon2 =
  block_diag(epsilon_hat, epsilon_hat)` (both blocks direct-rule), and the
  matrix fed into `kp` is `inv(epsilon_hat)` (the numerical matrix-inverse
  of the direct-rule Toeplitz) — Phase 2's separately-factorized
  `epsilon_inv_hat` (`toeplitz_matrix(..., inverse=True)`) is **not
  consumed by Phase 4a at all**; it remains 1D-only (Phase 3) infrastructure.
  Fixed `solve_layer_eigenmodes_patterned`'s signature to drop the
  now-unused `epsilon_inv_hat` parameter accordingly; see the function's
  docstring for the full line-by-line citation. This is a real, known
  RCWA accuracy limitation (ordinary Laurent's rule converges slower at 2D
  discontinuities than a proper vectorial/normal-vector factorization),
  not a design choice made here — a genuine 2D Li's-rule-equivalent
  (`use_polarization_basis` in S4) is out of Phase 4a's scope and would be
  a separately-requested extension. The bug had gone undetected because the
  test suite's own "ky=0 reduces to the 1D solver" check was circular (both
  solvers used the same, wrongly-copied formula, so of course they agreed)
  — replaced with two honest tests
  (`tests/test_2d_pillar.py::test_2d_patterned_ky_zero_te_block_matches_1d`,
  which *should* agree since the TE block never depends on which
  Fourier-factorization rule is used, and
  `..._tm_block_differs_from_1d`, which now correctly asserts the two
  *diverge*, guarding against this exact regression reappearing).
  Energy conservation and the `structures/via/pillar_array.py`/`via_array.py`
  R+T=1.0000 result both still hold after the fix (re-verified). Coupling
  subtraction: `[[diag(kx^2), diag(kx*ky)], [diag(kx*ky), diag(ky^2)]]`;
  `q = branch_select(eig(op))` reuses `_select_q_branch` unmodified. S4 was
  not runnable in this environment for a subprocess cross-check oracle (per
  `tasks.md`'s ☑ item on oracle determination); per `rules.md` AI Coding
  Rule 5 (never fabricate a match), the test strategy conservatively
  relies on: (a) the `ky=0` TE-block cross-check against the
  oracle-validated 1D solver, (b) the fully-uniform reduction to the Phase 1
  result, (c) energy conservation (R+T+sum(DE)=1 to 1e-8) across
  moderate-contrast Circle+Rectangle patterns at normal and oblique
  incidence, (d) both `structures/via/pillar_array.py` and
  `structures/via/via_array.py` running end-to-end with R+T=1.0000 across a
  21-point sweep from 500 nm to 1500 nm.

  **Follow-up, same day (user asked explicitly: stop relying on S4 alone,
  use the other vendored repos too)**: surveyed `EMpy`,
  `Rigorous-Coupled-Wave-Analysis`, and `RigorousCoupledWaveAnalysis.jl`
  for anything applicable to Phase 4a. `EMpy/EMpy/RCWA.py` is 1D-only (no
  2D support) and has an author-acknowledged instability hack — ruled out.
  `Rigorous-Coupled-Wave-Analysis`'s `run_RCWA_2D` also uses plain
  Laurent's rule for 2D (third independent confirmation the fix's rule
  choice isn't S4-specific), but has no hard-coded benchmark numbers.
  `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:57-99`
  (`eigenmodes(...,l::PatternedLayer)`) does — and its formula is a
  **structurally different derivation** (direct Maxwell-curl elimination
  into one matrix, not S4's `Epsilon2 @ kp` route). Julia isn't installed
  here (`which julia` fails), so hand-transcribed into
  `tests/oracles/rcwa_2djl_eigenvalues.py`, feeding it this project's own
  already-validated `epsilon_hat` to isolate the eigenoperator-construction
  step specifically. After reconciling two conventions (both confirmed
  empirically before trusting the comparison, not assumed): RCWA.jl
  normalizes `kx,ky` by `k0`, and its eigenvalues come out an overall sign
  flip of this project's `q^2` (opposite time convention) — its `q^2`
  values match `solve_layer_eigenmodes_patterned`'s to **~1e-12** across 6
  parametrized `num_orders`/angle cases
  (`tests/test_2d_pillar.py::test_2d_patterned_eigenvalues_match_rcwa_jl_oracle`).
  This is a genuinely independent formula agreeing with the fixed
  implementation — meaningfully stronger evidence than re-reading S4 a
  second time, and exactly the kind of check that would have caught the
  original bug (the old circular "ky=0 reduces to 1D" test could not).
  **Re-audited Phase 1-3's S4 citations too**, specifically re-reading each
  cited function's *complete* body (not just the branch matching prior
  usage) given that's exactly how the Phase 4a bug slipped through —
  `MakeKPMatrix`, `SolveLayerEigensystem_uniform` (Phase 1),
  `pattern_get_fourier_transform` (Phase 2), the 1D G-vector-selection and
  `fmm_closed.cpp` branches (Phase 3): all confirmed to match the existing
  implementation exactly, no further bugs found (worth noting: S4's
  `pattern.c` has its own explicit 1D-`RECTANGLE` case that independently
  corroborates `Slab`'s formula, a nice confirmation found during the
  re-read). A full external **R/T** oracle (not just eigenvalues) is still
  missing — no independently-published 2D benchmark exists in any vendored
  repo, S4 isn't buildable here, Julia isn't installed either — explicitly
  flagged, not faked, in `tests/oracles/rcwa_2d_pillar.py`; real remaining
  risk carried to Phase 4b, not just its originally-scoped near-degenerate
  stress cases. `simulation.py`'s `NotImplementedError` for 2D patterns is
  fully removed. 107 tests pass project-wide (101 + 6 new parametrized
  eigenoperator-oracle cases).

- **Phase 4b (near-degenerate/ill-conditioned 2D cases) is complete.**
  Stress-tested `solve_layer_eigenmodes_patterned` deliberately: index
  contrast from `n=3.48` to a lossy-metal-like `-20+2j`, `num_orders` up to
  225, near-touching circular pillars (`radius=0.49*period`), a
  sub-percent-halfwidth sliver rectangle, and near-degenerate nested
  circles (`1e-4`-scale radius difference). **Honest finding**: no
  catastrophic failure — `cond(epsilon_hat)` reached ~900, `cond(phi)`
  reached ~170 in the worst cases, energy conservation and the Phase 4a
  `RigorousCoupledWaveAnalysis.jl` eigenvalue oracle both held to ~1e-10
  throughout (`tests/test_2d_pillar_stress.py`). Reported as what was
  actually observed for the closed-form isotropic `Circle`/`Rectangle`
  patterns tested — not oversold as "no pathological case can exist
  anywhere." Since no failure needed fixing, the deliverable became
  **detection**: `eigenmodes.ILL_CONDITIONED_THRESHOLD = 1e4` (~10x
  headroom above the worst case observed) triggers a `WARNING` via a new
  module-level `logger` in `eigenmodes.py` (per `design.md`'s Logging
  Strategy) whenever `cond(epsilon_hat)` or `cond(phi)` exceeds it —
  verified by dedicated tests using `caplog`/`monkeypatch` (one confirming
  silence in the ordinary case, one forcing the threshold down to confirm
  the mechanism actually fires). The "match S4 or a published benchmark"
  deliverable could **not** be met as worded — confirmed (not assumed) S4
  needs `cmake`+Lua (neither present) and Julia isn't installed either
  (`which julia` fails) — so `tasks.md` records that box honestly unchecked
  rather than silently substituting the RCWA.jl oracle and calling it done;
  the substitution is used but disclosed. 118 tests pass project-wide (107
  + 11 new: 9 stress-case oracle cross-checks + 2 logging tests).
  `troubleshooting.md` updated with the Phase 4a bug (moved to
  Already-Solved Gotchas) and this phase's stress-test findings.
  **Phase 5 (tapered sidewalls) is next**, or Phase 6 (anisotropic
  materials) — both depend only on Phase 4a, not on closing Phase 4b's
  remaining external-R/T-oracle gap, which is now a standing item (not
  tied to a specific phase number) rather than blocking anything.

- **Phase 2 (Fourier-factorization core) is complete and validated.**
  `fourier_factorization.py`'s `pattern_epsilon_hat`/`toeplitz_matrix`
  build the direct and inverse-rule Toeplitz permittivity matrices from a
  `Pattern`, transcribed from `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`
  (lines 889-1029) and `S4/S4/fmm/fmm_closed.cpp::FMMGetEpsilon_ClosedForm`
  (lines 77-127). Validated against **two independent** numerical
  references (neither calling into the module under test): a from-scratch
  rasterize-and-sum, and a literal FFT-of-rasterized-mask reproduction of
  the vendored `Rigorous-Coupled-Wave-Analysis` (Python `convmat2D.py`)
  and `RigorousCoupledWaveAnalysis.jl` (`ft2d.jl::real2recip`)
  convolution-matrix algorithm — the user explicitly asked that S4 not be
  the only cross-check source, and `rules.md` names RCWA.jl as a
  sanctioned oracle. Caught one real bug during this: the first
  FFT-reference attempt used an uncentered raster grid that silently
  truncated a shape whose footprint crossed the domain edge; fixed via
  `numpy.fft.ifftshift` before the FFT (see `phases.md` Phase 2 for the
  full story). Both references agree with the analytic Toeplitz entries
  for `Circle` and `Rectangle` patterns, direct and inverse rule, at
  several nonzero G-vectors (`tests/test_fourier_factorization.py`, 12
  tests, all passing; 75 tests pass project-wide). Scalar isotropic
  materials only; anisotropic materials raise `NotImplementedError`
  naming Phase 6.
  `.flake8`/`mypy.ini` added (mypy/flake8 themselves are not installed in
  this dev environment, so linting was done by manual review, not an
  actual tool run — flag this for whoever next has the tools available).
  **Phase 3 (1D lamellar gratings) is next.**
- **Phase 1 (uniform multilayer core) is complete and validated.**
  Reflectance/transmittance for arbitrary uniform-layer stacks, arbitrary
  incidence angle/polarization, dispersive materials, and Jones/Mueller
  polarimetry all work. Validated by **two independent oracles**:
  `tests/oracles/fresnel.py` (from-scratch analytic Fresnel/TMM) and
  `tests/oracles/empy_tmm.py` (transcribed from the vendored EMpy
  reference library), cross-checked against the actual SiO2-on-Si
  structure in `tests/test_thin_film_empy_cross_check.py` across
  wavelength/angle/polarization, agreeing to `1e-8`. 68 tests pass.
- **Phases 2-9 are planned but not started** (see `phases.md`, `tasks.md`).
  `simulation.py:98` explicitly raises `NotImplementedError` for any
  patterned layer — this is the immediate next blocker for trench/via/pillar.
- The package was renamed `pyrcwa` → `sougata_solver` (all imports/docs
  updated) and is its own git repository, with a `.gitignore` covering
  `__pycache__/`, `.pytest_cache/`, `*.egg-info/`, `*.csv`, `*.png`,
  `outputs/`.
- `examples/` was removed and replaced with `structures/` (build a
  lattice/layer stack/materials, run the solver) and `postprocessing/`
  (derive Jones/Mueller matrices, ellipsometric angles, plots, and —
  planned — RI/thickness extraction, from a `structures/` script's raw
  output; never calls `Simulation.solve`). See `decisions.md` ADR-009.
  `structures/` is further grouped by category (`structures/thin_film/`
  today; `structures/trench/`, `structures/via/` etc. once those phases
  land). `polarimetry.py`'s `_decompose_sp` was made public (`decompose_sp`)
  so `postprocessing/jones_mueller_ellipsometry.py` can reuse the solver's
  exact s/p convention.
- **Every run gets its own timestamped output folder** (`output_paths.py`:
  `outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/`) containing its raw CSV/data
  *and* a `run_metadata.txt` (`write_run_metadata`) recording which script
  produced it and its key parameters — so re-running the same script with
  different settings never collides or gets mixed up (ADR-010). Plotting
  is a `postprocessing/` script (`plot_thin_film_rt.py`) that finds the
  relevant run's CSV (`find_latest_output`, or an explicit path) and saves
  its PNG back into that same run folder — plotting was briefly added
  directly to a `structures/` script and the user correctly caught that as
  a boundary violation; see ADR-010 for the fix and the reasoning.
- The full documentation set (this file, `README.md`, `PRD.md`,
  `architecture.md`, `design.md`, `rules.md`, `phases.md`, `tasks.md`,
  `decisions.md`, `testing.md`, `deployment.md`, `references.md`,
  `troubleshooting.md`) was just created in this session, before any
  Phase 2+ code was written — this is documentation-first, code-second by
  explicit user instruction.

## Important Decisions

Full rationale lives in `decisions.md` (ADR format). Summary:
- S-matrix (Redheffer star product), not transfer-matrix cascading —
  numerical stability for evanescent modes.
- Analytic shape Fourier transforms (S4-style), not raster+FFT
  (Meent/TORCWA-style) — accuracy for smooth boundaries.
- 1D lamellar gratings (Phase 3) come before 2D general patterned layers
  (Phase 4) — lower risk, validates Fourier-factorization machinery on the
  simpler (decoupled TE/TM) case first.
- Tapered sidewalls via staircase discretization, not new Fourier math —
  cheap, well-precedented (matches how even FEM tools like the vendored
  JCMsuite tutorials handle the same geometry).
- No arbitrary-polygon/GDS geometry support — parametric shapes
  (`Circle`, `Rectangle`, future `Slab`) only, by explicit user choice.
- GPU/autodiff backend (Meent/TORCWA-style) explicitly deferred to
  optional Phase 9, after Phases 2-8 are validated — by explicit user
  choice, so correctness work isn't chasing a moving numeric backend.
- Target audience is a solo research tool for now, not a public package —
  by explicit user choice; affects `deployment.md`'s scope (no PyPI/Docker
  yet).
- **FDTD is a genuine future goal (ADR-011, 2026-07-21)**, but a separate
  future effort, not a `sougata_solver` phase — `sougata_solver` stays
  RCWA-only (frequency-domain, periodic-BC); `PRD.md`'s Out-of-Scope
  wording was corrected after it was found to read as "FDTD rejected
  outright," which contradicted the actual intent. No FDTD design/timeline
  exists yet. `REFERENCE/meep`, `gprMax`, `fd3d`, `maxwellfdfd` (FDTD/FDFD)
  and `mfem`, `OpenParEM`, the FEniCS stack (FEM) are vendored for that
  future effort and remain unused/unevaluated by any current RCWA phase.

## Completed Milestones

- Phase 1: uniform multilayer stacks, dispersive materials, arbitrary
  polarization/angle, Jones/Mueller polarimetry — validated against
  analytic Fresnel/TMM (`tests/test_analytic_fresnel.py`).
- Phase 2: `fourier_factorization.py` (`pattern_epsilon_hat`,
  `toeplitz_matrix`) — validated against a from-scratch rasterize-and-sum
  reference (`tests/test_fourier_factorization.py`).
- Phase 3: `Lattice1D`/`Slab` (`geometry.py`), `truncate_fourier_orders_1d`
  (`fourier_basis.py`), `solve_layer_eigenmodes_1d` (`eigenmodes.py`),
  `Lattice1D` dispatch + `diffraction_efficiencies()` in `simulation.py` —
  validated against `tests/oracles/rcwa_1d_gaylord.py` and the
  energy-conservation/reduces-to-uniform invariants
  (`tests/test_1d_grating.py`).
- Phase 4a: `solve_layer_eigenmodes_patterned` (`eigenmodes.py`, general
  dense `2n x 2n` eigensolver, ordinary Laurent's rule per S4's true-2D
  closed-form path — corrected mid-session from an initially-wrong 1D-rule
  copy), 2D-patterned dispatch in `simulation.py`,
  `structures/via/pillar_array.py`, `structures/via/via_array.py` —
  validated by ky=0→1D TE-block reduction, fully-uniform reduction, energy
  conservation across Circle+Rectangle patterns, 21-point wavelength sweep
  producing R+T=1.0000 at each point, and an independent eigenvalue oracle
  transcribed from `RigorousCoupledWaveAnalysis.jl`
  (`tests/oracles/rcwa_2djl_eigenvalues.py`, agrees to ~1e-12)
  (`tests/test_2d_pillar.py`, 16 tests).
- Phase 4b: condition-number `WARNING` logging (`eigenmodes.ILL_CONDITIONED_THRESHOLD`)
  for `epsilon_hat`/`phi` in `solve_layer_eigenmodes_patterned` — validated
  by a deliberate near-degenerate/ill-conditioned stress sweep (no
  catastrophic failure found; energy conservation and the RCWA.jl oracle
  both held to ~1e-10) plus dedicated logging-mechanism tests
  (`tests/test_2d_pillar_stress.py`, 11 tests).
- Phase 5: `staircase.py` (`staircase_circle_layers`,
  `staircase_rectangle_layers`, `staircase_slab_layers`) — independently
  derived (no vendored-repo source exists for this discretization
  technique), validated by zero-taper regression to Phase 3/4a's
  single-layer results, energy conservation, and convergence-vs-`num_slices`
  studies (`tests/test_staircase.py`, 12 tests);
  `structures/via/tapered_via.py`, `structures/trench/tapered_trench.py`.

## Known Issues

- `.flake8`/`mypy.ini` now exist (added with Phase 2), but neither `flake8`
  nor `mypy` is actually installed in this dev environment — no run has
  verified the codebase is clean under either tool yet.
- `excitation.py`'s s/p polarization sign convention is explicitly
  documented as "internal and self-consistent (not yet matched to
  S4/EMpy's convention)" (`excitation.py:16-19`) — fine for Phase 1 (only
  power quantities are validated), but should be revisited if/when a
  polarization-sensitive cross-check against S4/EMpy is needed in a later
  phase.
- `scripts/` directory exists but is currently empty — no ad hoc utility
  scripts have been added yet.
- Whether S4 is actually built/runnable in this development environment
  (needed for Phase 4's cross-check oracle) has **not been checked yet** —
  first task when Phase 4 starts.

## Pending Tasks

See `tasks.md` for the full atomic checklist. Immediate next actions
(Phase 5 is now done; Phase 6 — anisotropic materials — is the natural next
phase, and Phase 7 — field reconstruction — is an equally valid
alternative, both depending only on Phase 4a):
1. Phase 6: generalize Phase 4a's eigensolver to a full 3x3 tensor
   `Epsilon2`, source a birefringent-material closed-form benchmark
   (uniaxial waveplate at normal incidence).
2. Phase 7: real-space field reconstruction (`fields.py` extension using
   `SMatrixStack.partial_smatrix_up_to`) and cross-section plotting
   examples for trench/via.
3. Standing item, not tied to a specific phase: a full external **R/T**
   oracle for 2D patterns (S4 subprocess, a published benchmark, or
   transcribing `Rigorous-Coupled-Wave-Analysis`'s `run_RCWA_2D`) — still
   open per Phase 4a/4b's honest accounting above; revisit if/when S4 or
   Julia becomes runnable in this environment, or a literature benchmark
   is located.

## Architecture Notes

- `SMatrixStack` and the rest of `smatrix.py` (plus `fields.py`,
  `excitation.py`) were confirmed **genuinely dimension-agnostic** by
  Phase 3 landing with **zero changes** to any of the three — they operate
  purely on `(q, phi, kp, thickness)` regardless of whether a layer is
  uniform or 1D-patterned. Phase 4a should need the same: only
  `eigenmodes.py`/`simulation.py` need new code for 2D patterns.
- `eigenmodes.build_kp_matrix` **already accepts** a full `(n,n)`
  `epsilon_inv` matrix, not just a scalar (see the `else` branch,
  `eigenmodes.py:40-47`) — confirmed in Phase 3 to be exactly the right
  interface: `solve_layer_eigenmodes_1d` calls it unmodified with Phase 2's
  `epsilon_inv_hat` Toeplitz. Phase 4a's `solve_layer_eigenmodes_patterned`
  should reuse it the same way.
- `SimulationResult.diffraction_efficiencies()` (new, Phase 3) computes
  per-order R/T by calling the existing `z_poynting_flux` once per order
  with other orders' amplitudes masked to zero, rather than a new per-order
  flux formula — valid because `all_modes[0]`/`all_modes[-1]` are *always*
  the uniform incidence/transmission half-spaces (`LayerStack.__init__`)
  regardless of what's patterned in between, and `build_kp_matrix`'s
  `epsilon_inv` branch has no cross-order coupling (each order's `kappa`
  entries only touch that order's own index and its `+n` pair). Phase 4a's
  2D case can reuse this unmodified too.
- See `architecture.md` for the full module responsibility table and data
  flow diagram.

## Technical Debt

- `excitation.py`'s polarization sign convention (see Known Issues above)
  is technical debt that's currently harmless but will need reconciling
  with S4's convention before any polarization-sensitive cross-check test
  can be written (relevant from Phase 4 onward if polarimetric, not just
  power, validation is needed against S4).
- No CI pipeline exists yet (git repo is brand new) — acceptable at current
  solo/local scope per `PRD.md`, but see `deployment.md` for the plan once
  a remote exists.

## Things Future AI Sessions Should Remember

- **Read `rules.md`'s "AI Coding Rules" section before writing any new
  physics code** — the project's entire trust model rests on never
  presenting an unverified formula as verified, and never fabricating a
  benchmark match.
- **This project's plan-mode scratch file**
  (`C:\Users\d14k4\.claude\plans\vivid-swimming-moler.md`) contains the
  originally-approved Phase 2-9 plan; `phases.md`/`tasks.md` in this repo
  are now the authoritative, living versions — keep them in sync going
  forward, don't treat the scratch plan file as the source of truth after
  this point.
- **Do not import from or modify the vendored reference repos**
  (`../S4`, `../EMpy`, `../RigorousCoupledWaveAnalysis.jl`, `../EMTutorial`)
  — they are read-only oracles.
- **`progress_log.md` (new, 2026-07-19)** is a dated, append-only log of
  discussions and their action items — check it at the start of any
  session for open `[ ]` items, verify against the actual code whether
  they've since been implemented, and add a new dated entry at the end of
  any substantive session. Distinct from this file (status snapshot) and
  `tasks.md` (phase-organized checklist).
- **The user's separate Claude-Code memory system**
  (`C:\Users\d14k4\.claude\projects\...\memory\`) is a different mechanism
  from this file — that one is cross-project and cross-session for the AI
  assistant's own use; this `memory.md` is project-scoped documentation
  living inside the `sougata_solver` repo itself, readable by any collaborator or
  future session regardless of which AI tool is used.
- **Phase 4 was split into Phase 4a (well-conditioned via/pillar case) and
  Phase 4b (near-degenerate/ill-conditioned stress cases), 2026-07-19**,
  at user request for more mathematical/scientific rigor in the roadmap —
  see `phases.md`, `tasks.md`, `PRD.md`. Phase 5/6/7 now depend on 4a only
  (4b hardens the solver in place, it doesn't unlock new capability).
  `testing.md` also gained a new Physical-Invariant Testing tier (energy
  conservation, convergence-rate-vs-Li-1996-theory), required starting
  Phase 3, oracle-independent and cheaper than the S4/benchmark
  cross-check — see `rules.md` Testing Requirements.
- **A `phase-reference-picker` skill now exists**
  (`.claude/skills/phase-reference-picker/`, at the `Solver_own` workspace
  root, sibling to `sougata_solver` and `REFERENCE`), added 2026-07-19.
  Invoke it before writing physics code for any new phase — it forces a
  real per-sub-task comparison across `REFERENCE/`'s ~18 vendored repos
  (not just S4, which is easy to default to since it's the Phase 1/2
  oracle already) and an explicit transcribe-vs-derive-independently
  decision. See `references.md`'s "Choosing a Reference for a New Phase"
  and `rules.md` AI Coding Rule 8.
