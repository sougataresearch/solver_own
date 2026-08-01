"""Status of Phase 4a's 2D oracle strategy (see `tests/test_2d_pillar.py`).

**Eigenoperator-level cross-check: closed.** `tests/oracles/rcwa_2djl_eigenvalues.py`
(hand-transcribed from `RigorousCoupledWaveAnalysis.jl`, a structurally
different eigenoperator derivation than S4's) confirms
`solve_layer_eigenmodes_patterned`'s `q^2` eigenvalues to ~1e-12 across
several `num_orders`/angle/pattern combinations. This is exactly the check
that would have caught the `Epsilon2` bug found and fixed this session (an
earlier draft copied the wrong S4 branch; both formulas agree now that the
fix matches S4's actual true-2D, no-polarization-basis behavior).

**Full R/T external-oracle comparison: still open, Phase 4b work.** No
independently-published 2D diffraction-efficiency benchmark was found
during this session's survey of the vendored RCWA-family repos:
- `Rigorous-Coupled-Wave-Analysis/RCWA_2D_examples/RCWA_photonic_circle_spectra.py`
  (+ `RCWA_functions/run_RCWA_simulation.py::run_RCWA_2D`, lines 13-146) is
  a genuine dense-eigensolve 2D implementation and could be transcribed
  into a full second-implementation R/T oracle (moderate effort -- ~7
  helper modules to pull together), but has no hard-coded reference
  numbers itself (only a "compare with Fan JOSA B" provenance comment, no
  actual transcribed values) -- it would give a second *implementation* to
  compare against, not independently-published ground truth.
- `RigorousCoupledWaveAnalysis.jl`'s `test/runtests.jl` (lines 70-111) uses
  `rand()` parameters and checks only internal self-consistency (its own
  ETM vs. SRCWA engines agree, R+T=1) -- no fixed 2D benchmark values
  either, and Julia is not installed in this environment to run it
  directly and freeze a reference number that way.
- S4 itself was not buildable/runnable in this environment (no `cmake`/Lua
  toolchain found -- see `memory.md`'s Known Issues) for a subprocess
  cross-check.

Per `rules.md` AI Coding Rule 5 (never fabricate a benchmark match), Phase
4a's correctness case rests on: the eigenoperator cross-check above, the
`ky=0` TE-block reduction to the oracle-validated 1D solver
(`test_2d_patterned_ky_zero_te_block_matches_1d`), the uniform-layer
reduction (`test_2d_patterned_layer_reduces_to_uniform_when_shapes_match_background`),
and energy conservation across moderate-contrast cases -- not a full
external R/T oracle. Closing that gap (transcribing `run_RCWA_2D`, sourcing
a literature benchmark, or getting S4/Julia runnable) is carried forward
explicitly as Phase 4b work, not silently dropped.
"""
