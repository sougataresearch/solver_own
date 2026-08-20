# Deployment — sougata_solver

Scope note (see `decisions.md` ADR-007): this is a **solo research tool
today**, run locally on Windows. This document is deliberately light —
no PyPI publishing, no Docker, no production servers — and says so
explicitly rather than padding with unneeded process. Revisit if this ever
becomes a shared/public project.

**Update (2026-08-19, ADR-038)**: ADR-007's own revisit condition — "used
by a second person" — was triggered. The resolution stays inside this
document's existing scope: distribution is by `git clone` from the
already-public `github.com/sougataresearch/solver_own` repo, not PyPI or
Docker. `GETTING_STARTED.md` (new) and `setup.ps1` (new, repo root) give a
new user a one-command path to a working, test-verified install. Full
PyPI publishing (version-pinning policy, `CHANGELOG.md`, release tagging)
remains a separate, not-yet-made decision — see ADR-038.

## Environment Setup

```powershell
cd C:\Users\sougata.bhunia\Desktop\Solver_own\sougata_solver
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

Requirements (from `pyproject.toml`): Python ≥ 3.10, `numpy>=1.24`,
`scipy>=1.10`; dev extra adds `pytest>=7.0`, `matplotlib>=3.7` (field/
result plotting, Phase 7 and Category 16), and `ruff>=0.16` (static
analysis, Category 17 target 17.5).

## Build Steps

`sougata_solver` builds as a standard `setuptools` src-layout package
(`[tool.setuptools.packages.find] where = ["src"]`). To build a wheel
locally (not currently published anywhere):

```bash
pip install build
python -m build
```

No build step is required for day-to-day development — `pip install -e .`
(editable install) is the normal workflow.

## Docker

**Not used.** No `Dockerfile` exists and none is planned at current scope
— a local virtualenv is sufficient for a single-developer numerical
library with no service component. If this changes (e.g. needing a
reproducible environment for sharing a specific validation run), prefer a
minimal `python:3.12-slim` image with `pip install -e ".[dev]"`, not a
speculative multi-stage build built ahead of actual need.

## CI/CD

**Set up as of Category 17 targets 17.2/17.3 (2026-08-07)** — two GitHub
Actions workflows under `.github/workflows/`:

- **`ci.yml`** ("CI (fast suite)"): on every push/PR to `main` and on
  manual dispatch, on `windows-latest` across a Python 3.10/3.11/3.12
  matrix (matching `pyproject.toml`'s `requires-python`) — `ruff check .`
  (target 17.5's static-analysis gate) then `pytest -m "not slow" -q`.
  `windows-latest`, not `ubuntu-latest` (the earlier placeholder above
  this section used to show), since this project has been developed and
  tested exclusively on Windows so far — matches this target's own
  explicit "Windows CI" wording.
- **`slow-tests.yml`**: `pytest -m slow -q` on a weekly cron schedule
  (Monday 06:00 UTC) plus manual dispatch, per target 17.3's own
  "schedule or manually trigger" wording — kept separate from `ci.yml` so
  convergence/benchmark studies (`tests/test_harmonic_convergence_matrix.py`'s
  4 `slow` cases alone take ~450s) don't gate every push.

No CD (continuous deployment) is applicable — there is no deployment
target (no server, no package registry) at current scope.

## Production Deployment

**Not applicable.** There is no "production" for this project today — it
is invoked directly via `python structures/*.py` / `python postprocessing/*.py` or imported into ad hoc
scripts. If/when public PyPI distribution is ever desired, that would be a
new, explicitly-scoped decision (see `decisions.md` ADR-007) requiring:
version pinning strategy, a `CHANGELOG.md`, and a release-tagging
convention — none of which exist or are needed yet.

## Rollback Strategy

Since there is no deployed service, "rollback" reduces to standard git
practice: never force-push a shared branch, tag known-good states once
Phase milestones complete (e.g. `git tag phase1-done` after this
documentation-creation commit, if useful), and rely on `pytest` passing as
the gate for trusting any given commit — a failing oracle-comparison test
is the actual "this isn't safe" signal in a physics codebase, more so than
any deployment mechanism.

## Monitoring

**Not applicable** — no running service to monitor. The closest analogue
is the test suite itself (`testing.md`) and the convergence studies
(Phase 5/8) acting as ongoing correctness monitors, run manually per
`rules.md`'s Code Review Checklist before trusting a commit.

## Logging

See `design.md`'s Logging Strategy section — the library itself does not
log (pure functions, raise-don't-log on error); `logging` module usage is
reserved for future numerically-concerning-but-not-fatal warnings (e.g.
ill-conditioned Toeplitz matrices in Phase 4), emitted at `WARNING` level,
never routine per-solve chatter.
