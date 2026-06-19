from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def read_key_value_summary(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["key", "value"]
        return {key: value for key, value in reader}


def test_solve_fd_heat_script_outputs(tmp_path: Path) -> None:
    outdir = tmp_path / "fd2"

    subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "solve_fd_heat.py"),
            "--n",
            "8",
            "--source",
            "square",
            "--outdir",
            str(outdir),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )

    summary = read_key_value_summary(outdir / "summary.csv")
    assert (outdir / "temperature.png").is_file()
    assert (outdir / "solution.npz").is_file()
    assert int(summary["dofs"]) == 64
    assert float(summary["relative_residual"]) < 1e-12
    assert float(summary["t_max"]) > 0.0


def test_solve_fd_heat_3d_script_outputs(tmp_path: Path) -> None:
    outdir = tmp_path / "fd3"

    subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "solve_fd_heat_3D.py"),
            "--n",
            "4",
            "--source",
            "gaussian",
            "--outdir",
            str(outdir),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )

    summary = read_key_value_summary(outdir / "summary.csv")
    assert (outdir / "temperature_midplane.png").is_file()
    assert (outdir / "temperature_3d.png").is_file()
    assert (outdir / "solution.npz").is_file()
    assert int(summary["dofs"]) == 64
    assert float(summary["relative_residual"]) < 1e-12
    assert float(summary["t_max"]) > 0.0
