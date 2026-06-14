#!/usr/bin/env python3
"""Run finite-difference mesh and parameter experiments."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from time import perf_counter

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

cache_dir = Path(__file__).resolve().parents[1] / ".cache"
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

import matplotlib.pyplot as plt
import numpy as np

from solve_fd_heat import build_laplacian_matrix, build_source, residual_metadata, solve_system, source_metadata


def run_case(n: int, source_kind: str, strength: float, solver: str, rtol: float) -> dict[str, float | int | str | None]:
    h = 1.0 / (n + 1)
    xx, yy, source = build_source(n, source_kind, strength)
    a = build_laplacian_matrix(n)
    rhs = source.reshape(-1)

    t0 = perf_counter()
    u, iterations = solve_system(a, rhs, solver, rtol)
    elapsed = perf_counter() - t0
    temp = u.reshape(n, n)

    return {
        "n": n,
        "dofs": n * n,
        "source_kind": source_kind,
        "strength": strength,
        **source_metadata(source, h),
        "solver": solver,
        "cg_iterations": iterations,
        "elapsed_seconds": elapsed,
        **residual_metadata(a, u, rhs),
        "matrix_nnz": a.nnz,
        "t_min": float(temp.min()),
        "t_max": float(temp.max()),
        "t_mean": float(temp.mean()),
    }


def write_csv(path: Path, rows: list[dict[str, float | int | str | None]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def plot_mesh_study(path: Path, rows: list[dict[str, float | int | str | None]]) -> None:
    ns = np.array([row["n"] for row in rows], dtype=float)
    dofs = np.array([row["dofs"] for row in rows], dtype=float)
    tmax = np.array([row["t_max"] for row in rows], dtype=float)
    elapsed = np.array([row["elapsed_seconds"] for row in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    axes[0].plot(dofs, tmax, marker="o")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("degrees of freedom")
    axes[0].set_ylabel("max temperature")
    axes[0].set_title("Mesh convergence indicator")

    axes[1].plot(ns, elapsed, marker="o", color="tab:red")
    axes[1].set_xlabel("n")
    axes[1].set_ylabel("seconds")
    axes[1].set_title("Runtime growth")

    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_strength_study(path: Path, rows: list[dict[str, float | int | str | None]]) -> None:
    strengths = np.array([row["strength"] for row in rows], dtype=float)
    tmax = np.array([row["t_max"] for row in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
    ax.plot(strengths, tmax, marker="o")
    ax.set_xlabel("source strength")
    ax.set_ylabel("max temperature")
    ax.set_title("Linearity check")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_residual_study(
    path: Path,
    mesh_rows: list[dict[str, float | int | str | None]],
    strength_rows: list[dict[str, float | int | str | None]],
) -> None:
    mesh_dofs = np.array([row["dofs"] for row in mesh_rows], dtype=float)
    mesh_residuals = np.array([row["relative_residual"] for row in mesh_rows], dtype=float)
    strengths = np.array([row["strength"] for row in strength_rows], dtype=float)
    strength_residuals = np.array([row["relative_residual"] for row in strength_rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    axes[0].plot(mesh_dofs, mesh_residuals, marker="o")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("degrees of freedom")
    axes[0].set_ylabel("relative residual")
    axes[0].set_title("Residual vs mesh size")

    axes[1].plot(strengths, strength_residuals, marker="o", color="tab:green")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("source strength")
    axes[1].set_ylabel("relative residual")
    axes[1].set_title("Residual vs source strength")

    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/experiments/fd_experiments"))
    parser.add_argument("--solver", choices=["direct", "cg"], default="direct")
    parser.add_argument("--source", choices=["square", "gaussian", "two_gaussian", "two_hotspots"], default="square")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--mesh-n", nargs="+", type=int, default=[16, 24, 32, 48, 64, 96, 128])
    parser.add_argument("--strengths", nargs="+", type=float, default=[25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    mesh_rows = [run_case(n, args.source, 100.0, args.solver, args.rtol) for n in args.mesh_n]
    write_csv(args.outdir / "mesh_study.csv", mesh_rows)
    plot_mesh_study(args.outdir / "mesh_study.png", mesh_rows)

    strength_rows = [run_case(64, args.source, s, args.solver, args.rtol) for s in args.strengths]
    write_csv(args.outdir / "strength_study.csv", strength_rows)
    plot_strength_study(args.outdir / "strength_study.png", strength_rows)
    plot_residual_study(args.outdir / "residual_study.png", mesh_rows, strength_rows)

    print("Experiment outputs:")
    print(args.outdir.resolve())
    print("mesh study:")
    for row in mesh_rows:
        print(f"  n={row['n']}, dofs={row['dofs']}, t_max={row['t_max']:.8g}, time={row['elapsed_seconds']:.4g}s")
    print("strength study:")
    for row in strength_rows:
        print(f"  strength={row['strength']}, t_max={row['t_max']:.8g}")


if __name__ == "__main__":
    main()
