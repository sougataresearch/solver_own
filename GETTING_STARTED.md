# Getting Started (new users)

This is the recipe for getting `sougata_solver` running on a machine that
has never touched it before. It assumes nothing about where you put the
code or what your username is.

Scope note (see `decisions.md` ADR-007/ADR-038): this project is
distributed by **git clone from its public GitHub repo**, not via PyPI or
Docker. `pip install -e .` still gives you the "it downloads everything it
needs" experience for `sougata_solver`'s own dependencies (`numpy`,
`scipy`, and — for the `dev` extra — `pytest`/`matplotlib`/`ruff`); PyPI
publishing of `sougata_solver` itself is a separate, not-yet-made decision.

## 1. Install Python

You need **Python 3.10, 3.11, or 3.12** (matches `pyproject.toml`'s
`requires-python` and the versions this project's CI actually tests
against).

- Download from [python.org/downloads](https://www.python.org/downloads/).
- During install, check **"Add python.exe to PATH"**.
- Verify it worked by opening a new PowerShell window and running:

```powershell
python --version
```

## 2. Get the code

Pick any folder you want (this does **not** need to match any particular
path — the setup script below works from wherever you clone it):

```powershell
git clone https://github.com/sougataresearch/solver_own.git
cd solver_own\sougata_solver
```

Don't have `git`? Install it from [git-scm.com](https://git-scm.com/downloads),
or download the repo as a ZIP from the GitHub page ("Code" → "Download ZIP")
and extract it instead.

## 3. Run the setup script

From inside the `sougata_solver` folder:

```powershell
powershell -ExecutionPolicy Bypass -File setup.ps1
```

This automatically:
1. Checks your Python version.
2. Creates an isolated virtual environment in `.venv` (so this doesn't
   touch any Python packages already on your system).
3. Runs `pip install -e ".[dev]"`, which downloads and installs
   `sougata_solver` itself plus `numpy`, `scipy`, `pytest`, `matplotlib`,
   and `ruff` — nothing to install by hand.
4. Runs the fast test suite (`pytest -m "not slow"`) to prove the install
   actually works before you rely on it.

If it ends with **"Setup complete: all tests passed."**, you're done. If
tests fail, the script tells you and stops — don't proceed on a red
install; see `troubleshooting.md` or share the output with the project
owner.

## 4. Activate the environment (every new terminal session)

The setup script only installs things — it doesn't leave the environment
"active" once it exits. Each time you open a new terminal to work on this
project:

```powershell
cd solver_own\sougata_solver
.venv\Scripts\Activate.ps1
```

You'll know it worked because your prompt gets a `(.venv)` prefix.

## 5. Run something

With the environment activated:

```powershell
python structures\thin_film\<some_example>.py
```

Browse `structures/` for worked examples by structure type (thin film,
trench/grating, via/pillar), each with its own `README.md`. For the
underlying library API, start with `README.md`'s Features section, which
links directly to the relevant module under `src/sougata_solver/`.

## 6. Run the tests yourself (optional, but recommended)

```powershell
pytest -m "not slow" -q     # fast suite -- same one setup.ps1 already ran
pytest -m slow -q           # slow convergence/benchmark studies (minutes)
ruff check .                # static analysis
```

## Troubleshooting

- **`python` not recognized** — Python wasn't added to PATH during install;
  re-run the Python installer and check that box, or add it manually.
- **`setup.ps1` won't run / "running scripts is disabled"** — that's
  Windows' default PowerShell execution policy blocking unsigned scripts.
  The `-ExecutionPolicy Bypass` flag in the command above already works
  around this for just this one invocation; it does not change any
  system-wide setting.
- **Anything else** — see `troubleshooting.md` for known failure modes, or
  `rules.md` for how this project expects bugs to be diagnosed (never
  silently loosen a tolerance or delete a failing oracle test to "fix" it).

## What you do *not* need

- The `REFERENCE/` folder (vendored oracle repos, ~1GB, one directory up
  from `sougata_solver/`) — it's cited in docstrings as the source of
  formulas for human review, but nothing in `sougata_solver` imports from
  it at runtime or in tests. You can develop and run the full test suite
  without it.
- Docker — none is used or required (`deployment.md`).
- A PyPI account — the package isn't published there; `git clone` +
  `pip install -e .` is the whole distribution mechanism today.

## Formula sources (checking a citation without `REFERENCE/`)

Docstrings and `decisions.md`/`references.md`/`theory.md` cite formulas by
exact repo/file/line (e.g. `S4/S4/fmm/fmm_closed.cpp:165-256`). You don't
need the vendored `REFERENCE/` copy to check one — clone the one original
public repo the citation names instead:

| Cited as | Public repo |
|---|---|
| `S4/...` | [github.com/victorliu/S4](https://github.com/victorliu/S4) |
| `EMpy/...` | [github.com/lbolla/EMpy](https://github.com/lbolla/EMpy) |
| `RigorousCoupledWaveAnalysis.jl/...` | [github.com/jonschlipf/RigorousCoupledWaveAnalysis.jl](https://github.com/jonschlipf/RigorousCoupledWaveAnalysis.jl) |
| `Rigorous-Coupled-Wave-Analysis/...` | [github.com/zhaonat/Rigorous-Coupled-Wave-Analysis](https://github.com/zhaonat/Rigorous-Coupled-Wave-Analysis) |
| `PyRCWA/...` | [github.com/vitamingcheng/PyRCWA](https://github.com/vitamingcheng/PyRCWA) |

(All five confirmed publicly reachable as of 2026-08-20.) `EMTutorial` is
the one exception — it's vendored JCMsuite project files, not an
independent public code repo, so there's nothing to link to; those
citations are illustrative geometry references, not formula sources, and
aren't needed to verify any equation. See `references.md` for the full
per-repo citation index and why each was chosen over alternatives.
