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
- **Degenerate/near-degenerate eigenvalues in the general non-uniform
  eigensolver (Phase 4).** `scipy.linalg.eig` on a matrix with two
  eigenvalues that are equal or very close can return a poorly-conditioned
  eigenvector basis, which then corrupts `phi`/`kp`-based quantities
  downstream (S-matrix, fields) even though `q` itself looks fine.
  **Symptom to watch for**: R/T that's sensitive to tiny perturbations in
  wavelength/angle where physically it shouldn't be, or that fails the
  "reduces to uniform case" regression test (`tasks.md` Phase 4) at low
  pattern contrast specifically (where near-degeneracy is most likely).
- **Toeplitz matrix ill-conditioning at high truncation order (Phase 2/4).**
  A very high `num_orders` with a high-index-contrast pattern can produce
  a badly-conditioned Toeplitz matrix. **Mitigation**: this is exactly the
  kind of thing the planned `WARNING`-level logging (see `design.md`'s
  Logging Strategy) should surface once implemented — condition-number
  checks are cheap (`numpy.linalg.cond`) and worth adding as a diagnostic
  in Phase 2, not just reacted to after a bad result.
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
- Whether **S4 is actually built/runnable** in this environment (needed
  for a live Phase 4 cross-check) has not been verified as of this
  writing — check this explicitly (e.g. can `S4/` be compiled, is its Lua
  interface accessible) before Phase 4 work assumes it's available; if it
  isn't, fall back to a literature benchmark per `rules.md`'s AI Coding
  Rules (never fabricate a match).

## When You Hit Something New

Add it here immediately, in the same format: what the gotcha is, why it's
subtle (looks-right-but-isn't), and what symptom would reveal it. Future
sessions (AI or human) should be able to search this file before spending
time rediscovering a bug class that's already been mapped.
