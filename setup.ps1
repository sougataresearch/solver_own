# Onboarding setup script for a fresh clone of sougata_solver.
#
# What this does (mirrors .github/workflows/ci.yml's own install/verify
# steps exactly, so "it worked here" means the same thing as "CI is green"):
#   1. Check the Python version meets pyproject.toml's requires-python.
#   2. Create a local .venv (skipped if one already exists).
#   3. pip install -e ".[dev]" -- pulls numpy/scipy/pytest/matplotlib/ruff
#      from PyPI automatically; no manual dependency list to maintain here.
#   4. Run the fast test suite (pytest -m "not slow") to confirm the
#      install actually works before the user starts editing anything.
#
# Usage (from anywhere -- this script cd's to its own folder first):
#   powershell -ExecutionPolicy Bypass -File setup.ps1
#
# Re-running is safe: an existing .venv is reused, not recreated.

$ErrorActionPreference = "Stop"

# Run relative to this script's own location, not the caller's cwd, so it
# works regardless of where the repo was cloned.
Set-Location -Path $PSScriptRoot

Write-Host "== sougata_solver setup ==" -ForegroundColor Cyan

# --- 1. Python version check -------------------------------------------
$pythonCmd = Get-Command python -ErrorAction SilentlyContinue
if (-not $pythonCmd) {
    Write-Host "ERROR: 'python' was not found on PATH." -ForegroundColor Red
    Write-Host "Install Python 3.10 or later from https://www.python.org/downloads/ (check 'Add python.exe to PATH' during install), then re-run this script."
    exit 1
}

$versionOutput = & python --version 2>&1
if ($versionOutput -match "Python (\d+)\.(\d+)") {
    $major = [int]$Matches[1]
    $minor = [int]$Matches[2]
    if ($major -lt 3 -or ($major -eq 3 -and $minor -lt 10)) {
        Write-Host "ERROR: found $versionOutput, but pyproject.toml requires Python >= 3.10." -ForegroundColor Red
        exit 1
    }
    Write-Host "Found $versionOutput (OK, >= 3.10 required)."
} else {
    Write-Host "WARNING: could not parse '$versionOutput' -- continuing anyway." -ForegroundColor Yellow
}

# --- 2. Create venv (skip if present) -----------------------------------
if (Test-Path ".venv") {
    Write-Host "Existing .venv found -- reusing it."
} else {
    Write-Host "Creating virtual environment in .venv ..."
    python -m venv .venv
}

$venvPython = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $venvPython)) {
    Write-Host "ERROR: venv creation appears to have failed -- $venvPython not found." -ForegroundColor Red
    exit 1
}

# --- 3. Install package + dev extras ------------------------------------
Write-Host "Installing sougata_solver + dependencies (numpy, scipy, pytest, matplotlib, ruff) ..."
& $venvPython -m pip install --upgrade pip
& $venvPython -m pip install -e ".[dev]"

# --- 4. Verify with the same fast test suite CI runs --------------------
Write-Host "Running the fast test suite to confirm the install works ..."
& $venvPython -m pytest -m "not slow" -q
$testExitCode = $LASTEXITCODE

Write-Host ""
if ($testExitCode -eq 0) {
    Write-Host "== Setup complete: all tests passed. ==" -ForegroundColor Green
    Write-Host "Activate the environment in future sessions with:"
    Write-Host "    .venv\Scripts\Activate.ps1" -ForegroundColor Cyan
} else {
    Write-Host "== Setup finished, but tests FAILED (exit code $testExitCode). ==" -ForegroundColor Red
    Write-Host "The environment was created and dependencies were installed, but something"
    Write-Host "is wrong -- do not assume the solver works. See troubleshooting.md, or share"
    Write-Host "the test output above with the project owner."
    exit $testExitCode
}
