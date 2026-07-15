# Roadmap — pyrcwa

This supersedes/formalizes the plan already approved and saved at
`C:\Users\d14k4\.claude\plans\vivid-swimming-moler.md`. Phase numbering
here is the authoritative one going forward; update this file (not the
plan-mode scratch file) as phases complete.

## Phase 1 — Uniform Multilayer Core — **DONE**

- **Objectives**: solve reflectance/transmittance for arbitrary stacks of
  uniform, isotropic, dispersive layers; support arbitrary incidence
  angle/polarization; report Jones/Mueller polarimetry.
- **Deliverables**: `materials.py`, `layer.py`, `eigenmodes.py` (uniform
  case), `smatrix.py`, `excitation.py`, `fields.py` (R/T), `polarimetry.py`,
  `simulation.py` (uniform path); Fresnel-oracle test suite; examples
  `01`-`04`.
- **Estimated complexity**: (retrospective) Medium — closed-form eigenmodes
  avoided the hardest numerical risk (general eigendecomposition), but the
  S-matrix sign/normalization conventions required careful,
  line-by-line-verified transcription from S4.
- **Dependencies**: none (foundation phase).
- **Status**: shipped, validated against analytic Fresnel/TMM
  (`tests/test_analytic_fresnel.py`).

## Phase 2 — Fourier-Factorization Core

- **Objectives**: build the dimension-agnostic infrastructure that turns a
  `Pattern` (shapes + background) into the Toeplitz permittivity matrices
  every patterned-layer eigensolver (Phase 3, 4) needs.
- **Deliverables**: new `src/pyrcwa/fourier_factorization.py` with
  `pattern_epsilon_hat(...)` and `toeplitz_matrix(...)`, producing both
  `epsilon_hat` (direct) and `epsilon_inv_hat` (inverse-rule) Toeplitz
  matrices; unit tests comparing analytic Toeplitz entries against a
  numerical FFT-of-rasterized-mask reference for `Circle` and `Rectangle`.
- **Estimated complexity**: Medium. The math is well-defined (sum of
  already-implemented shape Fourier transforms), but getting the
  direct-vs-inverse-rule distinction right, and validating it
  independently (not just "it compiles"), is the actual work.
- **Dependencies**: `geometry.py` (`Shape.fourier_transform`,
  `Pattern.containment_tree`) and `fourier_basis.py`
  (`truncate_fourier_orders`) — both already implemented, no changes
  needed to either.

## Phase 3 — 1D-Periodic Lamellar Gratings (Trench)

- **Objectives**: solve reflectance/transmittance/diffraction efficiencies
  for a 1D-periodic patterned layer (line/space, i.e. a trench), as the
  first end-to-end patterned-layer capability — chosen before 2D because
  1D gratings decouple TE/TM into independent scalar eigenproblems, making
  this the lower-risk place to validate the Phase 2 Fourier-factorization
  machinery and the general non-uniform eigenmode-solve pattern.
- **Deliverables**: `Lattice1D`, a `Slab`/`Line` shape in `geometry.py`;
  `truncate_fourier_orders_1d` in `fourier_basis.py`;
  `solve_layer_eigenmodes_1d(...)` in `eigenmodes.py`; a `Lattice1D`
  dispatch branch in `simulation.py`; a validation test against a published
  1D binary-grating benchmark (Moharam & Gaylord 1995 or equivalent);
  `examples/05_trench_grating.py`.
- **Estimated complexity**: Medium-High. The physics (decoupled scalar
  TE/TM eigenproblems) is simpler than Phase 4's general case, but this
  phase is also where the *first* non-uniform eigensolver gets built and
  validated, so unexpected issues surfacing here are likely to be
  Fourier-factorization bugs from Phase 2, not 1D-specific ones — budget
  time accordingly.
- **Dependencies**: Phase 2.

## Phase 4 — 2D-Periodic Patterned Layers (Via, Pillar)

- **Objectives**: solve reflectance/transmittance/diffraction efficiencies
  for a full 2D-periodic patterned layer using the existing `Circle`/
  `Rectangle` shapes — i.e. make via and pillar arrays actually work,
  removing the `NotImplementedError` at `simulation.py:98`.
- **Deliverables**: `solve_layer_eigenmodes_patterned(...)` in
  `eigenmodes.py` (general non-uniform eigenproblem, transcribed from
  `S4/S4/rcwa.cpp::SolveLayerEigensystem`, lines 794-827); `simulation.py`
  wiring for the 2D patterned path; an S4-cross-check validation test (or,
  if S4 isn't buildable/runnable in this environment, an explicitly-flagged
  literature benchmark instead — never a fabricated "it matches" claim per
  `rules.md`'s AI coding rules); `examples/06_pillar_array.py`,
  `examples/07_via_array.py`.
- **Estimated complexity**: High — this is the highest-risk phase in the
  whole roadmap (general complex eigendecomposition, possible
  near-degenerate eigenvalues, the trickiest part of `rcwa.cpp` to
  transcribe correctly). See `troubleshooting.md` for known failure modes
  to watch for once this lands.
- **Dependencies**: Phase 2; benefits from Phase 3 having already
  shaken out Fourier-factorization bugs on the simpler 1D case.

## Phase 5 — Tapered / Sloped Sidewalls (Via, Trench)

- **Objectives**: represent a via or trench with linearly tapered
  sidewalls via staircase (z-discretized) layer approximation, and
  demonstrate R/T convergence as slice count increases.
- **Deliverables**: a small staircase-layer-stack generator (given top/
  bottom feature size, thickness, and slice count `N`, produce `N` `Layer`s
  with linearly interpolated `Circle`/`Rectangle`/`Slab` sizes); a
  convergence-vs-`N` test/example for both a tapered via and a tapered
  trench (marked `slow` per the existing pytest marker).
- **Estimated complexity**: Low — no new Fourier/eigenmode math, purely a
  layer-stack generation helper consuming Phase 3/4's already-validated
  per-layer solvers.
- **Dependencies**: Phase 3 and Phase 4 (needs at least one working
  patterned-layer eigensolver to stack; ideally both, to cover tapered
  trench and tapered via).

## Phase 6 — Anisotropic Materials

- **Objectives**: support full 3×3 permittivity tensors (already exposed by
  `Material.epsilon_tensor`) in both uniform and patterned layers, removing
  `simulation.py`'s anisotropic `NotImplementedError`.
- **Deliverables**: generalize Phase 4's eigensolver to accept a full
  tensor `Epsilon2` rather than only scalar/diagonal; validation against a
  known birefringent-material benchmark (e.g. a uniaxial waveplate at
  normal incidence, checked against closed-form ordinary/extraordinary
  index behavior).
- **Estimated complexity**: Medium — the eigensolver machinery from Phase 4
  already generalizes to tensors mathematically; the work is mostly
  correctly wiring `Material.epsilon_tensor`'s off-diagonal terms through
  and validating the coupling terms are handled right.
- **Dependencies**: Phase 4 (reuses/extends its general eigensolver).

## Phase 7 — Real-Space Field Reconstruction & Visualization

- **Objectives**: reconstruct E/H(x,y,z) on a grid at an arbitrary depth in
  the stack, and produce cross-section field-intensity plots for trench/via
  structures.
- **Deliverables**: extend `fields.py` to sum Fourier components onto a
  real-space grid using `SMatrixStack.partial_smatrix_up_to` (already
  implemented, `smatrix.py:174-178`) for mode amplitudes at intermediate
  depths; example scripts producing `matplotlib` cross-section plots for a
  trench and a via.
- **Estimated complexity**: Medium — the amplitude bookkeeping is already
  in place; the new work is the inverse-Fourier-sum and its own
  correctness check (e.g. field continuity across a layer interface as a
  sanity test).
- **Dependencies**: Phase 3 or 4 (need at least one working patterned-layer
  case worth visualizing) — although the *machinery* for field
  reconstruction is dimension-agnostic and could technically be built
  against Phase 1's uniform stacks first as a stepping stone.

## Phase 8 — Expanded Validation Suite & Example Gallery

- **Objectives**: systematic convergence-vs-`num_orders` studies for every
  geometry type; a complete example gallery mirroring the vendored
  `EMTutorial` reference cases (thin film, multistack/DBR, trench, via,
  pillar, tapered via).
- **Deliverables**: `slow`-marked convergence tests per geometry type;
  `examples/` entries for each structure type in `PRD.md`'s Success
  Criteria.
- **Estimated complexity**: Low-Medium — mostly systematic application of
  patterns already established in Phases 3-7, not new algorithmic risk.
- **Dependencies**: Phases 3-7 (validates all of them).

## Phase 9 — Performance & Optional GPU/Autodiff Backend (later, optional)

- **Objectives**: vectorize wavelength/angle sweeps in NumPy first;
  optionally, add a torch/JAX array backend behind the same function
  signatures (Meent/TORCWA-style) for GPU batching and autodiff-based
  inverse design — only after correctness is fully validated through
  Phase 8.
- **Deliverables**: a vectorized-sweep code path with a regression test
  proving numerical equivalence to the unvectorized path (see `rules.md`'s
  Performance Requirements); optionally, a backend abstraction layer.
- **Estimated complexity**: High if the optional GPU/autodiff backend is
  pursued (requires re-expressing every eigensolve/S-matrix operation in a
  backend-agnostic way); Low-Medium for the NumPy-only vectorization step
  alone.
- **Dependencies**: Phases 2-8 (explicitly deferred until correctness is
  solid — see `decisions.md`).

## Phase Sequencing Summary

```
Phase 1 (done) ──► Phase 2 ──► Phase 3 ──► Phase 4 ──► Phase 5
                                   │            │
                                   └────► Phase 6 (extends Phase 4)
                                                │
                       Phase 3/4 ──────────► Phase 7 ──► Phase 8 ──► Phase 9 (optional)
```
