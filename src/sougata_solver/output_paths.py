"""Date-stamped output folder shared by structures/ and postprocessing/ scripts.

Every script-generated file (CSV, PNG, ...) lands under
`outputs/YYYY-MM-DD/HH-MM-SS_<run_name>/`, so runs from the same day stay
grouped by date but each individual run -- whether it's a different script or
a rerun of the same one -- gets its own subfolder and never overwrites or
gets mixed up with another run's files.
"""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

OUTPUTS_ROOT = Path(__file__).resolve().parents[2] / "outputs"


def dated_output_dir(run_date: date | None = None) -> Path:
    """`outputs/YYYY-MM-DD/` for `run_date` (default: today), created if needed."""
    d = OUTPUTS_ROOT / (run_date or date.today()).isoformat()
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_output_dir(run_name: str) -> Path:
    """A fresh subfolder for one script invocation:
    `outputs/YYYY-MM-DD/HH-MM-SS_<run_name>/`, created if needed -- so two
    different scripts, or two runs of the same script, never collide.
    """
    timestamp = datetime.now().strftime("%H-%M-%S")
    d = dated_output_dir() / f"{timestamp}_{run_name}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_output_path(run_name: str, filename: str) -> Path:
    """Full path for `filename` inside a fresh subfolder for this run
    (see `run_output_dir`)."""
    return run_output_dir(run_name) / filename


def find_latest_output(filename: str) -> Path:
    """Most recently written `outputs/*/*/filename` (newest run subfolder
    first -- folder names sort chronologically since they're `HH-MM-SS_...`
    inside a `YYYY-MM-DD` date folder).

    For scripts (e.g. postprocessing) that read a file an earlier script
    wrote, possibly on a previous day -- so a same-day workflow still works
    with no path editing, and a next-day rerun of just the postprocessing
    step still finds the most recent matching run.
    """
    candidates = sorted(OUTPUTS_ROOT.glob(f"*/*/{filename}"), reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No {filename!r} found under {OUTPUTS_ROOT}")
    return candidates[0]
