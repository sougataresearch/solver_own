"""Category 15 target 15.6: CLI `run` command tests -- exit codes and
output-file behavior per `cli.py`'s own documented design (target 15.5)."""

from __future__ import annotations

import json

import pytest

from sougata_solver import cli

_VALID_CONFIG = {
    "unit": "m",
    "lattice": {"type": "2d", "a": [1e-6, 0.0], "b": [0.0, 1e-6]},
    "materials": {"air": {"eps_re": 1.0}, "glass": {"eps_re": 2.25}, "MgF2": {"eps_re": 1.9044}},
    "layers": [{"name": "MgF2", "thickness": 0.55e-6 / (4 * 1.38), "material": "MgF2"}],
    "incidence": "air",
    "transmission": "glass",
    "num_orders": 1,
    "excitation": {"wavelength": 0.55e-6, "theta_deg": 30.0, "phi_deg": 0.0, "s_amplitude_re": 1.0},
}


def test_run_valid_config_exits_zero_and_writes_output(tmp_path, capsys):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(_VALID_CONFIG), encoding="utf-8")
    out_dir = tmp_path / "out"

    exit_code = cli.main(["run", str(config_path), "--output-dir", str(out_dir)])

    assert exit_code == cli.EXIT_OK
    captured = capsys.readouterr()
    assert "reflectance:" in captured.out
    assert "transmittance:" in captured.out
    assert (out_dir / "result.json").exists()
    assert (out_dir / "run_metadata.txt").exists()

    result = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert 0.0 <= result["reflectance"] <= 1.0
    assert result["reflectance"] + result["transmittance"] == pytest.approx(1.0, abs=1e-9)


def test_run_missing_file_exits_config_error(tmp_path, capsys):
    exit_code = cli.main(["run", str(tmp_path / "does_not_exist.json")])
    assert exit_code == cli.EXIT_CONFIG_ERROR
    assert "error:" in capsys.readouterr().err


def test_run_malformed_json_exits_config_error(tmp_path, capsys):
    config_path = tmp_path / "bad.json"
    config_path.write_text("{not valid json", encoding="utf-8")
    exit_code = cli.main(["run", str(config_path)])
    assert exit_code == cli.EXIT_CONFIG_ERROR


def test_run_schema_violation_exits_config_error(tmp_path, capsys):
    bad_config = json.loads(json.dumps(_VALID_CONFIG))
    bad_config["num_orders"] = -1
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(bad_config), encoding="utf-8")

    exit_code = cli.main(["run", str(config_path)])

    assert exit_code == cli.EXIT_CONFIG_ERROR
    assert "num_orders" in capsys.readouterr().err


def test_no_command_raises_systemexit():
    with pytest.raises(SystemExit):
        cli.main([])
