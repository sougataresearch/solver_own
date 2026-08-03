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
  supported near-grazing/near-threshold behavior (interpolation, a
  documented exclusion zone, or similar) is `COMMERCIAL_RCWA_ATOMIC_TARGETS.md`
  Category 6 target 6.4's ("Grazing-incidence boundary test") job, not
  yet done.

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
