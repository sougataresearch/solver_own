# Troubleshooting — sougata_solver

Known numerical gotchas, organized by where they'll bite. This file exists
because in RCWA, the failure mode is almost never a crash — it's a
plausible-looking wrong number. Add to this file the moment you find a new
one; don't rely on memory across sessions.

## Already-Solved Gotchas (documented in code, worth knowing)

- **`E = phi @ (a+b)` looks right but is wrong — that's actually `H`.**
  The correct tangential E-field relation is `u = kp @ phi @ (a-b) / (omega*q)`
  with an index swap and sign flip (`Ex = u[n:]`, `Ey = -u[:n]`). See
  `fields.py::tangential_e_field`'s docstring — this exact mistake is
  "commonly-paraphrased-but-wrong" per that docstring, i.e. it's a known
  trap in RCWA writeups generally, not just an internal note.
- **Branch selection for `q = sqrt(q_sq)` needs special handling near the
  real axis.** A naive `np.sqrt` on complex `q_sq` can put a purely-real
  (propagating, lossless) mode's `q` on the wrong side of the branch cut
  due to floating-point noise in the imaginary part. `eigenmodes.py::_select_q_branch`
  handles near-real `q_sq` specially (exact real/imaginary split) before
  falling back to "flip principal root if `Im(q) < 0`" for the general
  case. Any new eigensolver (Phase 3, 4, 6) **must reuse this function**,
  not reimplement branch selection from scratch.
- **Never form `inv(A)` directly.** `smatrix.py::_solve` uses
  `scipy.linalg.lu_factor`/`lu_solve` instead of `numpy.linalg.inv` —
  explicit matrix inversion is both slower and less numerically stable
  than an LU solve for the interface-matrix systems here. Reuse `_solve`
  for any new linear-system solve; don't reintroduce `np.linalg.inv`.
- **Transfer matrices blow up; that's why S-matrices exist.** If you ever
  find yourself tempted to multiply per-layer transfer matrices directly
  (e.g. "just to check something quickly"), don't — see ADR-001 in
  `decisions.md`. Even a quick diagnostic script should use the existing
  `SMatrixStack`/star-product machinery.
- **A source-cited formula can still be wrong if only the matching branch
  was read (Phase 4a).** `solve_layer_eigenmodes_patterned`'s `Epsilon2`
  construction was first written by copying `solve_layer_eigenmodes_1d`'s
  `epsilon_hat`/`inv(epsilon_inv_hat)` formula — plausible, cited the right
  file, and passed every test, including a "ky=0 reduces to the 1D solver"
  check that seemed like a real cross-check. It was wrong: that formula is
  S4's **1D-only** special case (`fmm_closed.cpp:110-132`,
  `0==Lr[2]&&Lr[3]==0`); the adjacent true-2D branch
  (`fmm_closed.cpp:133-139`) uses plain Laurent's rule for *both* blocks
  and was never read. The "ky=0" test couldn't catch this because both
  solvers shared the identical wrong formula — a **circular test that
  passes for the wrong reason looks identical to a real regression guard**
  until you ask what would make it fail. **Symptom to watch for**: any
  time two of this project's own functions are cross-checked against each
  other (not an external oracle) and always agree suspiciously exactly —
  ask whether they could share a bug, not just a correct formula. **Fixed
  by**: actually reading the full source function (not just the branch
  matching prior usage) whenever citing it as "the general case," and by
  building a genuinely independent oracle
  (`tests/oracles/rcwa_2djl_eigenvalues.py`, a structurally different
  formula from `RigorousCoupledWaveAnalysis.jl`) instead of only
  cross-checking this project's own code against itself. See
  `eigenmodes.solve_layer_eigenmodes_patterned`'s docstring and `phases.md`
  Phase 4a's Status for the full account.
- **Degenerate/near-degenerate eigenvalues and Toeplitz ill-conditioning
  (Phase 4b) — stress-tested, not found to be catastrophic in practice.**
  A deliberate sweep (index contrast from `3.48` to a lossy-metal-like
  `-20+2j`, `num_orders` up to 225, near-touching circular pillars, a
  sub-percent-halfwidth sliver rectangle, near-degenerate nested circles
  with a `1e-4`-scale radius difference) found `cond(epsilon_hat)` up to
  ~900 and `cond(phi)` up to ~170 in the worst cases tried — energy
  conservation and the independent `RigorousCoupledWaveAnalysis.jl`
  eigenvalue oracle both held to ~1e-10 throughout, no case actually broke.
  This is an honest empirical finding for the *closed-form isotropic
  Circle/Rectangle* patterns tested, not a proof that no pathological case
  exists — `numpy.linalg.eig` can still misbehave on an input this sweep
  didn't probe (e.g. a shape with a true measure-zero degeneracy engineered
  on purpose, or a much larger `num_orders`). **Mitigation shipped**:
  `cond(epsilon_hat)`/`cond(phi)` logged at `WARNING`
  (`eigenmodes.ILL_CONDITIONED_THRESHOLD`, `1e4`, ~10x headroom above the
  worst observed case) rather than silently returning a possibly-degraded
  result — detection, not silent correction. See
  `solve_layer_eigenmodes_patterned`'s docstring and
  `tests/test_2d_pillar_stress.py` for the frozen stress cases and the
  logging tests.
- **A diffraction order sitting exactly at the Rayleigh/Wood's-anomaly
  threshold produces `NaN` R/T, not a degraded-but-finite answer (found
  while writing Category 1 target 1.8's mode-classification test).** At
  the exact threshold, `q == 0` for that order, and
  `smatrix.py::interface_smatrix`'s `kp @ phi / q[None, :]` construction
  divides by zero. Confirmed directly (not assumed): evaluating
  `SimulationResult.diffraction_efficiencies()` at the exact threshold
  wavelength returns `NaN` for every order, with
  `RuntimeWarning: divide by zero encountered in divide` /
  `invalid value encountered in divide` from `smatrix.py:75`. This is a
  genuine, pre-existing solver limitation at the exact singular point (a
  physically infinite-length evanescent decay / infinitely-slow group
  velocity at a grazing/Wood's-anomaly order), not a bug introduced by
  target 1.8's mode-classification work, and not something that work
  attempted to fix. `tests/test_mode_classification.py`'s Rayleigh test
  deliberately checks a small relative step away from the threshold on
  each side (`0.999x`/`1.001x`) rather than the exact point. Defining a
  supported near-grazing/near-threshold behavior was
  `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 6 target 6.4's
  ("Grazing-incidence boundary test") job — see the next entry for what
  that target found for the closely related exact-grazing-*incidence*
  case (as opposed to this entry's exact-Rayleigh-*threshold* case; same
  `q==0` mechanism, different trigger).
- **Exact grazing incidence (`theta=90 deg`) raises a plain `ValueError`,
  not a `NaN` (Category 6 target 6.4, `tests/test_grazing_incidence.py`).**
  Same root mechanism as the Rayleigh-threshold entry above (`q==0` in
  `interface_smatrix`'s `kp @ phi / q` division) but a different failure
  mode: here the *entire* incidence-medium mode set shares `q=0` (not just
  one order among many), so the resulting non-finite `b_l` matrix is
  caught by `scipy.linalg.lu_factor`'s internal `asarray_chkfinite` check
  and raised as `ValueError` before ever reaching a silent `NaN`. Confirmed
  this happens at *exactly* `theta=90 deg` in float64, not merely "very
  close to it," because `math.sin(math.radians(90.0)) == 1.0` exactly (a
  genuine floating-point coincidence for the nearest representable double
  to `pi/2`), which for an `n=1` (air) incidence medium makes
  `kx0 = omega*sin(theta)` exactly equal to `omega`, giving `q_sq =
  omega^2 - kx0^2 == 0.0` exactly. **Supported range**: any `theta < 90
  deg` — confirmed finite and energy-conserving (`R+T=1` to `1e-8`+) up to
  `theta=89.999 deg`; only the exact endpoint fails, and fails loud, not
  silently.
- **A "lossy metal" index copied from a source using the opposite time
  convention is actually a gain medium here — `R+T>1`, not a bug (found
  twice: Category 2 target 2.5, then guarded against in Category 5 target
  5.4).** This project's convention is `d/dt -> -i*omega`
  (`CONVENTIONS.md`), which requires `Im(eps) > 0` for a passive/absorbing
  medium. `n = -20+2j` (a value copied verbatim from Phase 4b's own stress
  sweep, itself picked without checking sign convention) squares to
  `eps = 396-80j`, `Im(eps) < 0` — a gain medium, giving `R+T` up to ~17
  through the full `Simulation.solve()` pipeline. **Symptom to watch for**:
  `R+T` noticeably greater than 1 for a material you intended to be
  absorptive. **Fix**: check the sign of `Im(eps)` (or `Im(n*conj(n))`)
  before trusting a "lossy" material value taken from a different source
  library/paper — don't assume every published metal index already matches
  this project's phasor convention. Category 5's `Material.from_lorentz`/
  `from_drude`/`from_drude_lorentz` docstrings each independently re-derive
  and test the correct sign (`tests/test_dispersion_models.py`) rather than
  assuming the transcribed vendored formula's sign already matches.
- **`Path.write_text(...)` without an explicit `encoding` uses the platform
  default (`cp1252` on Windows), which raises `UnicodeEncodeError` on
  non-ASCII text (found when a citation containing "Rakić" was written to
  `run_metadata.txt` for the first time, Category 5 target 5.8).** Fixed
  in `output_paths.write_run_metadata` by passing `encoding="utf-8"`
  explicitly; any new code that reads/writes text files in this project
  should do the same rather than relying on the platform default, per
  `PRD.md`'s Windows-primary-shell constraint.
- **`fields.z_poynting_flux` is missing the textbook `0.5` time-average
  factor relative to `Sz = 0.5*Re(Ex*conj(Hy) - Ey*conj(Hx))` -- confirmed
  directly, not a bug (found while building Category 9's real-space field
  reconstruction, target 9.6).** `z_poynting_flux`'s own modal quadratic
  form gives exactly `2x` the textbook per-order flux for the same mode --
  harmless everywhere it's actually used (`reflectance()`/`transmittance()`
  are *ratios* of two `z_poynting_flux` outputs, so the factor of 2
  cancels), but a real-space Poynting-flux integral built independently
  from `modal_field_components`' `Ex/Ey/Hx/Hy` must use the **no-`0.5`**
  form (`Sz = Re(Ex*conj(Hy) - Ey*conj(Hx))`) to match this project's
  established convention, confirmed to reproduce `R`/`T` to full double
  precision once accounted for. See `CONVENTIONS.md`'s "Real-space field
  reconstruction" section and `tests/test_field_reconstruction.py`.
- **Interior field/amplitude reconstruction (`smatrix.interior_amplitudes`
  + `fields.propagate_amplitudes`) can numerically overflow for thick,
  highly lossy, high-`num_orders` layers -- found while implementing
  Category 7 target 7.6's layer-wise absorption, a real numerical-
  conditioning limitation, not a formula bug.** `propagate_amplitudes`'s
  backward-wave exponential, `b(z) = b_top * exp(-i*q*z)`, grows as
  `exp(Im(q)*z)` for `Im(q) >= 0` (the correct, already-validated
  convention for a *stationary reference point deep inside a decaying
  medium* -- confirmed exactly by Category 9's field-continuity tests at
  moderate `Im(q)*thickness`). For a high-loss material (`eps=-396+80j`,
  the same fixture `tests/test_stress_regression.py` uses) at high
  `num_orders`, the deepest evanescent Fourier orders can have `Im(q)`
  large enough that `Im(q)*thickness` reaches the 30s-and-up range for a
  layer that isn't even especially thick in absolute terms
  (`thickness=0.3`, `num_orders=25`: `max(Im(q))*thickness ~= 38`,
  `exp(38) ~= 3e16`) -- one such term then dominates the
  `z_poynting_flux` sum and produces a nonsensical result (a measured
  `layer_absorption()` of `~573`, versus the correct `R+T+A=1` identity
  that same material/pattern satisfies exactly at `thickness=0.05`,
  `max(Im(q))*thickness ~= 6.3`). This is the same general instability
  class the transfer-matrix method is famous for (`decisions.md`'s
  original choice of S-matrix cascading over transfer matrices was made
  for exactly this reason) — it doesn't affect `R`/`T` themselves (the
  S-matrix cascade used for those never reconstructs an interior
  amplitude), only interior-point queries (Category 9's field
  reconstruction, Category 7's per-layer absorption) for genuinely
  extreme thickness/loss/order combinations. Not fixed (no formula
  change) — `layer_absorption()`/field reconstruction should be used with
  parameters where `max(Im(q))*thickness` stays a modest double-digit
  exponent or less; `tests/test_layer_absorption.py::test_interior_amplitude_reconstruction_can_numerically_overflow_for_thick_lossy_layers`
  is a regression guard on the failure symptom itself, not a fix, so a
  future silent workaround (e.g. clipping/renormalizing) doesn't land
  unnoticed. See `decisions.md` ADR-017.

## Anticipated Gotchas (not yet encountered — flagged ahead of Phase 2-6)

- **Direct vs. inverse-rule Fourier factorization (Phase 2).** Using
  `inv(epsilon_hat_toeplitz)` where `epsilon_inv_hat_toeplitz` (built from
  `1/eps(x,y)` directly, per shape) is required will produce a
  plausible-converging-but-wrong answer, especially for TM-like
  polarization at a material interface with high index contrast. This is
  the single most common historical RCWA bug (see Li 1996 in
  `references.md`). **Symptom to watch for**: convergence with `num_orders`
  that's suspiciously slow, or a result that doesn't match a known limit
  (e.g. the Fresnel limit as pattern contrast goes to zero).
- **1D-vs-2D lattice convention mismatch (Phase 3).** Do not implement
  `Lattice1D` by reusing 2D `Lattice` with one basis vector set to a very
  large period — this was explicitly rejected during planning (see
  `phases.md` Phase 3's rationale) because it introduces spurious weak
  coupling along the "infinite" direction rather than a true decoupled 1D
  formulation. If a future session is tempted to take this shortcut for
  convenience, don't — implement the genuine 1D TE/TM formulation instead.
- **Staircase convergence for steep sidewall angles (Phase 5).** A very
  steep taper angle may require a surprisingly large `N` to converge.
  **Mitigation**: the Phase 5 convergence-vs-`N` test/example is mandatory
  specifically to catch this per-structure rather than assuming a fixed
  `N` is always sufficient (see `tasks.md` Phase 5).
- **Polarization sign-convention mismatch when cross-checking against S4/EMpy
  (Phase 4/6).** `excitation.py`'s s/p convention is explicitly documented
  as not yet matched to S4/EMpy's (see `memory.md` Known Issues). A
  polarization-resolved cross-check (as opposed to a scalar-power-only
  check like Phase 1's) may show a sign or axis-swap mismatch that is a
  *convention* difference, not a physics bug — reconcile the convention
  explicitly before concluding a discrepancy is a real bug.

## Environment-Specific Notes

- Development is on **Windows with PowerShell** as the primary shell (per
  `PRD.md` Constraints) — any new tooling (linting, CI YAML, scripts) must
  be verified to work there, not assumed to work only in a Unix shell.
- **Neither S4 nor Julia is runnable in this environment** — confirmed
  (not just assumed) during Phase 4a/4b: `S4` needs `cmake`+a Lua
  toolchain, neither found; `which julia` fails (needed to run
  `RigorousCoupledWaveAnalysis.jl` directly, e.g. to freeze a reference
  number). This is why Phase 4a/4b's oracle strategy is hand-transcription
  of vendored source into `tests/oracles/` rather than a live subprocess
  cross-check — re-check this if the environment ever changes (a new
  toolchain installed) rather than assuming it's still true.

## When You Hit Something New

Add it here immediately, in the same format: what the gotcha is, why it's
subtle (looks-right-but-isn't), and what symptom would reveal it. Future
sessions (AI or human) should be able to search this file before spending
time rediscovering a bug class that's already been mapped.
