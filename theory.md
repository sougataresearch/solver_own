# Theory — Mathematical Foundations of `sougata_solver`

Targets 18.1 ("theory outline"), 18.2 ("core derivation"), and 18.3
("anisotropy derivation") of `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category
18, closed out together in one document.

## Scope and how to use this document

This is **not** a from-scratch re-derivation of RCWA — every formula
mentioned here is already implemented, cited to a named external source,
and validated against an independent oracle elsewhere in this repo (see
`rules.md`'s Documentation Standards and the AI Coding Rules; nothing in
this file introduces a coefficient or sign not already present and tested
in code). This document's job is the one thing no single existing doc
does: **tie the pipeline together end to end** and act as a table of
contents into the detailed sources — `design.md` (algorithm-level
citations), `CONVENTIONS.md` (frozen sign/normalization conventions), and
`s_matrix_method.md` (a full worked numeric derivation of the S-matrix
stage, cross-checked against KLA). Read this first; follow the links below
for the depth you need.

For the public function signatures implementing each stage, see
[`api_reference.md`](api_reference.md). For runnable, already-validated
examples, see [`tutorials.md`](tutorials.md). For what each oracle
comparison does and does not prove, see
[`validation_guide.md`](validation_guide.md).

## The pipeline, end to end

A single `Simulation.solve()` call moves through five stages. This is the
one place that narrates all five as one pipeline; each stage's own
citation lives in `design.md`'s `## Algorithms` section (linked per stage
below).

```
 geometry (Pattern)                materials (Material)
        |                                  |
        v                                  v
 1. Fourier factorization  <---------------+
    (Toeplitz eps_hat / epsilon_inv_hat)
        |
        v
 2. Per-layer eigenmode solve
    (q, phi, kp)  -- uniform: closed-form: patterned/anisotropic: general eig()
        |
        v
 3. S-matrix cascading
    (Redheffer star product across the LayerStack)
        |
        v
 4. Field / power extraction
    (R, T, diffraction efficiencies, full E/H(x,y,z) reconstruction)
```

Stage 1 only applies to *patterned* layers — a uniform layer's
permittivity is already spatially constant, so there is nothing to
Fourier-factorize (`design.md` Algorithm 1 below). Stages 2–4 apply to
every layer type.

### Stage 1 — Fourier factorization

Source: `design.md` **Algorithm 3** ("Fourier factorization", Toeplitz
permittivity construction) and **Algorithm 3a** (the rule inventory —
which of direct-rule / inverse-rule / numerical-matrix-inverse each solver
branch actually uses, with a citation table pinned by
`tests/test_fourier_factorization_rules.py`).

The one-sentence summary: for a patterned layer, `eps(x,y)`'s Fourier
coefficients at each retained reciprocal-lattice vector `G` are assembled
into a Toeplitz matrix `M[i,j] = eps_hat(G_i - G_j)`. **Which** Toeplitz
matrix (`eps_hat` directly, or a numerical inverse of it, or a *separately
factorized* `1/eps` Toeplitz — Li's 1996 inverse rule) depends on the
solver branch; getting this distinction wrong is, per `design.md`, "the
single most common source of wrong-but-plausible-looking RCWA results
industry-wide." Do not assume; consult the Algorithm 3a table for the
exact rule any given code path uses.

### Stage 2 — Per-layer eigenmode solve

Source: `design.md` **Algorithm 1** (uniform isotropic layer — closed-form
`q[i] = branch_select(eps*omega^2 - kx[i]^2 - ky[i]^2)`) and **Algorithm 2**
(general/patterned layer — the full eigenproblem
`(q^2, phi) = eig(Epsilon2 @ kp)`, this project's highest-risk algorithm,
validated in two phases: Phase 4a on moderate-contrast cases, Phase 4b
deliberately stress-testing near-degenerate/high-contrast cases).

Every layer, regardless of type, produces the same three objects — `q`
(propagation constants), `phi` (eigenvector basis), `kp` (k-parallel
operator) — which is exactly what lets Stage 3 (S-matrix cascading) treat
every layer type identically.

### Stage 3 — S-matrix cascading

Source: `design.md` **Algorithm 4**, and the full worked derivation in
[`s_matrix_method.md`](s_matrix_method.md) — physical motivation
(`s_matrix_method.md` §1–3: why transfer matrices blow up numerically for
evanescent modes in a thick/lossy layer, and how the S-matrix formulation
keeps every intermediate quantity bounded instead), the Redheffer star
product construction (§4–5), and three fully worked numeric examples
(§6–8: one-layer SiO2/Si, two-layer TiO2+SiO2/Si, three-layer SiO2+SiO+Ni/Si)
each cross-checked against `smatrix.py`'s actual output and, for the
first example, against an external reference (KLA) — see
`s_matrix_method.md` §6 Step 6.

If you want to understand *why* the S-matrix method is used at all (not
just what it computes), start at `s_matrix_method.md` §1–3 before Stage 4
below — that motivation isn't repeated here.

### Stage 4 — Field / power extraction

Source: `design.md` **Algorithm 5**. Reflectance/transmittance come from
`fields.z_poynting_flux` (transcribed, not re-derived — see `design.md`'s
note on why a from-scratch re-derivation of the sign/normalization
conventions embedded in `kp`/`phi` was deliberately avoided). Full
real-space `(Ex,Ey,Ez,Hx,Hy,Hz)` field reconstruction (Phase 7) extends
this via `fields.modal_field_components`/`reconstruct_field_at_points`,
using `smatrix.interior_amplitudes` to get local mode amplitudes at an
arbitrary depth, then summing the retained Fourier series in real space.

## Conventions you need before reading any equation above

`CONVENTIONS.md` is the frozen reference for every sign/normalization
choice the pipeline above depends on — read it alongside, not instead of,
the stage descriptions:

- **Phasor/propagation convention** (`CONVENTIONS.md` "Phasor and
  propagation convention"): `d/dt -> -i*omega`, forward propagation factor
  `exp(+i*q*z)`. Not interchangeable with a textbook convention using the
  opposite time sign without also flipping the associated spatial signs.
- **Modal vectors** (`CONVENTIONS.md` "Modal vectors and normalization"):
  `phi` is an internal eigenvector basis, not individually
  power-normalized; the tangential-field reconstruction is
  `u = kp @ phi @ ((a-b)/(omega*q))` with `Ex = u[n:2n]`, `Ey = -u[0:n]` —
  i.e. the internal transverse ordering is `u = [-Ey; Ex]`, not `[Ex; Ey]`.
- **S-matrix direction convention** (`CONVENTIONS.md` "S-matrix direction
  convention"): `[a_right; b_left] = S @ [a_left; b_right]`, fixed across
  every layer and Fourier order.
- **Polarization convention** (`CONVENTIONS.md` "Polarization convention"):
  `s_hat`/`p_hat_xy` basis; a worked table of `(s_amplitude, p_amplitude)`
  pairs for TE/TM/linear/circular/elliptical states, matched to a
  commercial RCWA tool's own angle convention (`decisions.md` ADR-033).
- **Permittivity tensor ordering** (`CONVENTIONS.md` "Permittivity tensor
  ordering"): `Material.epsilon_tensor` returns Cartesian
  `[[exx,exy,exz],[eyx,eyy,eyz],[ezx,ezy,ezz]]` — documented ahead of full
  anisotropic support specifically so later milestones could be tested
  without ambiguity.

## Anisotropic materials

`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 1 tracks anisotropy at finer
grain than the original `phases.md` Phase 6 entry. As of 2026-08-03,
targets 1.1–1.4 and 1.6–1.8 are shipped:

| Target | Capability | Solver entry point | Validated against |
|---|---|---|---|
| 1.3 | Uniform diagonal tensor (uniaxial, normal incidence) | `eigenmodes.solve_layer_eigenmodes_uniform_diagonal` | Closed-form birefringence formula + Fresnel/TMM oracle per principal axis (`tests/test_anisotropic_uniform.py`) |
| 1.4 | Uniform in-plane-coupled tensor (`exx,exy,eyx,eyy,ezz`) | `eigenmodes.solve_layer_eigenmodes_uniform_inplane` | Independent eigenoperator oracle from `RigorousCoupledWaveAnalysis.jl` (`tests/oracles/rcwa_anisotropic_inplane_jl.py`, agrees to ~1e-13) |
| 1.6 | Patterned (2D-periodic) anisotropic layers | `eigenmodes.solve_layer_eigenmodes_patterned_inplane`, `fourier_factorization.toeplitz_matrix_component` | Reduction to Phase 4a's isotropic solver, reduction to target 1.4's uniform solver, energy conservation for a genuinely patterned lossless case (`tests/test_anisotropic_patterned.py`) |
| 1.7 | Deterministic mode-ordering for near-degenerate eigenvalues | `eigenmodes._canonical_mode_order` | Repeated-solve determinism, sort-key unit tests, energy conservation on a near-isotropic patterned case (`tests/test_anisotropic_degeneracy.py`) |
| 1.8 | Propagating/evanescent mode classification | `eigenmodes.classify_propagating`, `SimulationResult.order_classification()` | Analytic Rayleigh-threshold wavelength cross-check (`tests/test_mode_classification.py`) |

**Target 1.5 — longitudinal tensor coupling (`eps_xz`, `eps_yz`, `eps_zx`,
`eps_zy`) is explicitly deferred, not implemented.** `simulation.py`
raises `NotImplementedError` naming this target for any layer with a
non-zero longitudinal component. This is a "not found this session," not
a permanent "cannot exist" conclusion: a bounded literature search (per
`rules.md` AI Coding Rule 1) found general-anisotropic-RCWA literature
exists in principle, but no source that is both (a) actually
fetchable/readable as full text in this environment, and (b)
independently benchmarkable per Rule 5 (a second, structurally-different
source to cross-check against, the way targets 1.3/1.4 had S4 + RCWA.jl).
See `references.md`'s "Target 1.5 bounded literature search" entry and
`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` target 1.5's own entry for the full
account, including which specific candidate sources were checked and why
each was rejected. **Do not write code that assumes a full 9-component
tensor is supported** — `Material.epsilon_tensor` can store one
(`CONVENTIONS.md`'s tensor-ordering note), but the solver only accepts the
diagonal/in-plane subset above.
