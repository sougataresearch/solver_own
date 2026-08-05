# `profiling/` — Diagnostic Timing Scripts

Category 12 target 12.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`). These
scripts measure wall-clock time (and, where noted, peak memory) for the
solver's three big linear-algebra stages — per-layer eigenmode solve,
matrix-solve/inverse operations, and S-matrix cascading — on a small set
of fixed, representative fixtures (thin-film, 1D grating, 2D pillar).

**Not tests.** Nothing here is asserted against a hard time/memory limit:
wall-clock timing is machine-dependent, and per `rules.md`'s Performance
Requirements ("correctness first, no premature optimization"), this
project does not gate correctness on absolute speed. These scripts exist
so that a later optimization decision (Category 12 targets 12.3/12.5, or
any future Category 13 "Performance optimization" work) has an actual
measured baseline to justify itself against, the same discipline
`decisions.md` ADR-016 already used for the Category 7 Toeplitz-matrix
cache — never optimize on a hunch.

Run with:

```bash
python profiling/baseline_profile.py
```

Output is printed only (no CSV/`outputs/` folder — this is diagnostic
tooling, not a `structures/` physics run, so it doesn't follow
ADR-009/010's raw-data-output convention).
