#!/usr/bin/env python3
"""Solve a 2D steady heat diffusion problem with finite differences.

Model:
    -Delta T = f on [0, 1] x [0, 1]
    T = 0 on the boundary

This script is intentionally small and explicit. It is meant for numerical
experiments before moving to the finite element implementation.
"""

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
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import cg, spsolve


def build_source(n: int, kind: str, strength: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return grid coordinates and source values on interior points."""
    h = 1.0 / (n + 1)
    x = np.linspace(h, 1.0 - h, n)
    y = np.linspace(h, 1.0 - h, n)
    xx, yy = np.meshgrid(x, y, indexing="xy")

    if kind == "square":
        f = np.zeros_like(xx)
        mask = (xx > 0.4) & (xx < 0.6) & (yy > 0.4) & (yy < 0.6)
        f[mask] = strength
    elif kind == "gaussian":
        sigma = 0.06
        r2 = (xx - 0.5) ** 2 + (yy - 0.5) ** 2
        f = strength * np.exp(-r2 / (2.0 * sigma**2))
    elif kind in {"two_gaussian", "two_hotspots"}:
        sigma = 0.06
        r2_a = (xx - 0.35) ** 2 + (yy - 0.45) ** 2
        r2_b = (xx - 0.68) ** 2 + (yy - 0.58) ** 2
        f = strength * np.exp(-r2_a / (2.0 * sigma**2))
        f += strength * np.exp(-r2_b / (2.0 * sigma**2))
    else:
        raise ValueError(f"Unknown source kind: {kind}")

    return xx, yy, f


def build_laplacian_matrix(n: int):
    """Build the sparse matrix for -Delta with zero Dirichlet boundaries."""
    h = 1.0 / (n + 1)
    main = 2.0 * np.ones(n)
    off = -1.0 * np.ones(n - 1)
    t = diags([off, main, off], [-1, 0, 1], format="csr")
    i = eye(n, format="csr")
    a = (kron(i, t) + kron(t, i)) / h**2
    return a


def source_metadata(source: np.ndarray, h: float) -> dict[str, float | int]:
    """Return simple grid/source scale diagnostics."""
    active_points = int(np.count_nonzero(source))
    cell_area = h**2
    return {
        "grid_size_h": h,
        "grid_cell_area": cell_area,
        "heat_source_active_points": active_points,
        "heat_source_active_area": active_points * cell_area,
        "heat_source_integral": float(source.sum() * cell_area),
    }


def residual_metadata(a, u: np.ndarray, rhs: np.ndarray) -> dict[str, float]:
    """Return absolute and relative residual diagnostics for A u = rhs."""
    residual = rhs - a @ u
    abs_residual = float(np.linalg.norm(residual))
    rhs_norm = float(np.linalg.norm(rhs))
    rel_residual = abs_residual / rhs_norm if rhs_norm > 0.0 else abs_residual
    return {
        "absolute_residual": abs_residual,
        "relative_residual": rel_residual,
    }


def solve_system(a, rhs: np.ndarray, solver: str, rtol: float) -> tuple[np.ndarray, int | None]:
    """Solve A u = rhs and return solution plus CG iteration count if relevant."""
    if solver == "direct":
        return spsolve(a, rhs), None

    iterations = 0

    def count_iteration(_):
        nonlocal iterations
        iterations += 1

    try:
        u, info = cg(a, rhs, rtol=rtol, atol=0.0, callback=count_iteration)
    except TypeError:
        u, info = cg(a, rhs, tol=rtol, atol=0.0, callback=count_iteration)
    if info != 0:
        raise RuntimeError(f"CG did not converge, scipy returned info={info}")
    return u, iterations


def save_outputs(
    outdir: Path,
    xx: np.ndarray,
    yy: np.ndarray,
    source: np.ndarray,
    temp: np.ndarray,
    metadata: dict[str, str | int | float | None],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    np.savez(
        outdir / "solution.npz",
        x=xx,
        y=yy,
        heat_source=source,
        temperature=temp,
        **metadata,
    )

    with (outdir / "summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["key", "value"])
        for key, value in metadata.items():
            writer.writerow([key, value])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)

    im0 = axes[0].pcolormesh(xx, yy, source, shading="auto", cmap="inferno")
    axes[0].set_title("Heat source f(x, y)")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].pcolormesh(xx, yy, temp, shading="auto", cmap="viridis")
    axes[1].contour(xx, yy, temp, colors="white", linewidths=0.6, alpha=0.7)
    axes[1].set_title("Temperature T(x, y)")
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    fig.colorbar(im1, ax=axes[1])

    fig.savefig(outdir / "temperature.png", dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=64, help="Interior grid points per direction.")
    parser.add_argument(
        "--source",
        choices=["square", "gaussian", "two_gaussian", "two_hotspots"],
        default="square",
        help="Heat source type.",
    )
    parser.add_argument("--strength", type=float, default=100.0, help="Source strength.")
    parser.add_argument("--solver", choices=["direct", "cg"], default="direct")
    parser.add_argument("--rtol", type=float, default=1e-10, help="CG relative tolerance.")
    parser.add_argument("--outdir", type=Path, default=Path("results/baselines/fd_single"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n < 3:
        raise ValueError("--n must be at least 3")

    h = 1.0 / (args.n + 1)
    xx, yy, source = build_source(args.n, args.source, args.strength)
    a = build_laplacian_matrix(args.n)
    rhs = source.reshape(-1)

    t0 = perf_counter()
    u, iterations = solve_system(a, rhs, args.solver, args.rtol)
    elapsed = perf_counter() - t0
    temp = u.reshape(args.n, args.n)

    metadata = {
        "n": args.n,
        "dofs": args.n * args.n,
        "source_kind": args.source,
        "strength": args.strength,
        **source_metadata(source, h),
        "solver": args.solver,
        "cg_iterations": iterations,
        "elapsed_seconds": elapsed,
        **residual_metadata(a, u, rhs),
        "matrix_nnz": a.nnz,
        "t_min": float(temp.min()),
        "t_max": float(temp.max()),
    }
    save_outputs(args.outdir, xx, yy, source, temp, metadata)

    print("Solved 2D heat diffusion problem")
    for key, value in metadata.items():
        print(f"{key}: {value}")
    print(f"outputs: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
