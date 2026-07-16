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
    (see `run_output_dir`). Only safe to call once per run -- a second
    call creates a *different* timestamped subfolder. A script that writes
    more than one file (e.g. a CSV and a metadata file) should call
    `run_output_dir` once instead and build every path under that same
    directory (see `write_run_metadata`)."""
    return run_output_dir(run_name) / filename


def write_run_metadata(output_dir: Path, script_path: str, **params: object) -> Path:
    """Write a human-readable `run_metadata.txt` into `output_dir` recording
    which script produced this run and its key parameters (materials,
    thicknesses, angle, wavelength range, ...) -- so a run folder can be
    told apart from another run of the same script (or a different script)
    without re-reading code or guessing from the timestamp alone.

    Call once per run, alongside whatever data file(s) the script writes
    into the same `output_dir` (from `run_output_dir`).
    """
    lines = [f"script: {script_path}", f"run_at: {datetime.now().isoformat(timespec='seconds')}", ""]
    lines.extend(f"{key}: {value}" for key, value in params.items())
    path = output_dir / "run_metadata.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


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
