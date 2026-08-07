"""Category 15 targets 15.5/15.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
command-line interface for the validated `config.py` workflow.

## CLI design (target 15.5)

One subcommand, `run`, for one job: load a JSON configuration
(`config.simulation_from_dict`'s schema, see `config.py`'s module
docstring), validate it, solve it, and report reflectance/transmittance.

    sougata-solver run <config.json> [--output-dir DIR]

**Exit codes** (checked by `tests/test_cli.py`, not just documented here):

- `0` -- solved successfully; reflectance/transmittance printed to stdout
  and written to the output directory.
- `2` -- the configuration file could not be read or was invalid: missing
  file, malformed JSON, or a schema violation caught by `config.py`'s
  construction-time validation (target 15.3). This never reaches
  `Simulation.solve()`, matching `config.py`'s own "validated before any
  numerical calculation" guarantee.
- `1` -- any other, unexpected failure (e.g. a `LinAlgError` from the
  solver itself). Kept distinct from `2` so a caller can tell "your config
  was wrong" apart from "the solve itself failed."

**Output location**: reuses `output_paths.py`'s existing
`outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/` convention (already used by every
`structures/*.py` example script) rather than inventing a second one --
`--output-dir` overrides it with a caller-chosen directory instead, for
scripted/CI use where a fixed, predictable path is wanted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sougata_solver import config
from sougata_solver.output_paths import run_output_dir, write_run_metadata

EXIT_OK = 0
EXIT_RUNTIME_ERROR = 1
EXIT_CONFIG_ERROR = 2


def _run(config_path: str, output_dir: str | None) -> int:
    path = Path(config_path)
    try:
        text = path.read_text(encoding="utf-8")
        data = json.loads(text)
        simulation, excitation = config.simulation_from_dict(data)
    except OSError as exc:
        print(f"error: could not read {config_path!r}: {exc}", file=sys.stderr)
        return EXIT_CONFIG_ERROR
    except (json.JSONDecodeError, ValueError) as exc:
        print(f"error: invalid configuration {config_path!r}: {exc}", file=sys.stderr)
        return EXIT_CONFIG_ERROR

    try:
        result = simulation.solve(excitation)
        r, t = result.reflectance(), result.transmittance()
    except Exception as exc:  # noqa: BLE001 -- deliberately broad: any solver failure maps to EXIT_RUNTIME_ERROR
        print(f"error: solve failed: {exc}", file=sys.stderr)
        return EXIT_RUNTIME_ERROR

    print(f"reflectance:   {r:.6f}")
    print(f"transmittance: {t:.6f}")
    print(f"R + T:         {r + t:.6f}")

    out_dir = Path(output_dir) if output_dir else run_output_dir(path.stem)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_run_metadata(out_dir, str(path), reflectance=r, transmittance=t)
    (out_dir / "result.json").write_text(json.dumps({"reflectance": r, "transmittance": t}), encoding="utf-8")
    print(f"wrote results to {out_dir}")
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="sougata-solver")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="solve a validated JSON simulation configuration")
    run_parser.add_argument("config_path", help="path to a config.py-schema JSON file")
    run_parser.add_argument("--output-dir", default=None, help="write results here instead of the dated outputs/ folder")

    args = parser.parse_args(argv)
    if args.command == "run":
        return _run(args.config_path, args.output_dir)
    parser.error(f"unknown command {args.command!r}")  # pragma: no cover -- unreachable, argparse enforces choices
    return EXIT_RUNTIME_ERROR


if __name__ == "__main__":
    sys.exit(main())
