# Deployment — pyrcwa

Scope note (see `decisions.md` ADR-007): this is a **solo research tool
today**, run locally on Windows. This document is deliberately light —
no PyPI publishing, no Docker, no production servers — and says so
explicitly rather than padding with unneeded process. Revisit if this ever
becomes a shared/public project.

## Environment Setup

```powershell
cd c:\Users\d14k4\Desktop\Solver_own\pyrcwa
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

Requirements (from `pyproject.toml`): Python ≥ 3.10, `numpy>=1.24`,
`scipy>=1.10`; dev extra adds `pytest>=7.0`. `matplotlib` will need to be
added as a dev/example dependency once Phase 7 (field visualization) lands
— see `tasks.md`.

## Build Steps

`pyrcwa` builds as a standard `setuptools` src-layout package
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

**Not yet set up** — `pyrcwa` was just initialized as its own git
repository in this session (see `decisions.md` ADR-008) and has no remote
yet. Once a remote (e.g. GitHub) exists, the minimum useful CI is a single
GitHub Actions workflow:

```yaml
# .github/workflows/test.yml (add once a remote exists)
name: test
on: [push, pull_request]
jobs:
  pytest:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - run: pip install -e ".[dev]"
      - run: pytest -m "not slow"
```

Keep `slow`-marked convergence studies out of the default CI run (they're
for local investigation, not every-push validation) — run them manually or
in a separate, manually-triggered workflow if they ever become expensive
enough to matter.

No CD (continuous deployment) is applicable — there is no deployment
target (no server, no package registry) at current scope.

## Production Deployment

**Not applicable.** There is no "production" for this project today — it
is invoked directly via `python examples/NN_*.py` or imported into ad hoc
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
