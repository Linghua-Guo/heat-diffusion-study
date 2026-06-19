#!/usr/bin/env python3
"""FEniCSx version of the 2D steady heat diffusion solver.

This script requires dolfinx, mpi4py, and PETSc through petsc4py. It is kept
separate from the finite-difference scripts because the current local
environment may not have FEniCSx installed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

try:
    from mpi4py import MPI
    from dolfinx import fem, io, mesh
    from dolfinx.fem.petsc import LinearProblem
    import ufl
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing FEniCSx dependencies. Install dolfinx/mpi4py first, then rerun this script."
    ) from exc


SOURCE_CHOICES = ("square", "gaussian", "two_gaussian", "two_hotspots")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument(
        "--source",
        choices=SOURCE_CHOICES,
        default="square",
        help="Heat source type.",
    )
    parser.add_argument("--strength", type=float, default=100.0)
    parser.add_argument("--outdir", type=Path, default=Path("results/baselines/fenicsx_single"))
    parser.add_argument("--ksp-type", default="cg", help="PETSc KSP solver type.")
    parser.add_argument("--pc-type", default="hypre", help="PETSc preconditioner type.")
    parser.add_argument("--ksp-rtol", type=float, default=1e-10, help="PETSc KSP relative tolerance.")
    parser.add_argument(
        "--no-png",
        action="store_true",
        help="Skip the matplotlib quick-look PNG output.",
    )
    return parser.parse_args()


def evaluate_source(x: np.ndarray, kind: str, strength: float) -> np.ndarray:
    """Return source values at FEniCSx interpolation points."""
    values = np.zeros(x.shape[1], dtype=np.float64)

    if kind == "square":
        inside = (x[0] > 0.4) & (x[0] < 0.6) & (x[1] > 0.4) & (x[1] < 0.6)
        values[inside] = strength
    elif kind == "gaussian":
        sigma = 0.06
        r2 = (x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2
        values = strength * np.exp(-r2 / (2.0 * sigma**2))
    elif kind in {"two_gaussian", "two_hotspots"}:
        sigma = 0.06
        r2_a = (x[0] - 0.35) ** 2 + (x[1] - 0.45) ** 2
        r2_b = (x[0] - 0.68) ** 2 + (x[1] - 0.58) ** 2
        values = strength * np.exp(-r2_a / (2.0 * sigma**2))
        values += strength * np.exp(-r2_b / (2.0 * sigma**2))
    else:
        raise ValueError(f"Unknown source kind: {kind}")

    return values


def create_mesh_and_space(n: int):
    """Create the unit-square triangular mesh and P1 Lagrange space."""
    domain = mesh.create_unit_square(MPI.COMM_WORLD, n, n, mesh.CellType.triangle)
    v_space = fem.functionspace(domain, ("Lagrange", 1))
    return domain, v_space


def create_source(v_space, kind: str, strength: float):
    """Interpolate a configured source into the finite element space."""
    source = fem.Function(v_space)
    # Caveat for square: the intended source is discontinuous, but interpolation
    # into continuous P1 Lagrange gives a nodal approximation, not a cellwise
    # constant source. Smooth sources such as Gaussian are naturally better
    # matched to this representation.
    source.interpolate(lambda x: evaluate_source(x, kind, strength))
    return source


def create_dirichlet_bc(domain, v_space):
    """Create zero-temperature Dirichlet boundary conditions."""
    def boundary(x):
        return (
            np.isclose(x[0], 0.0)
            | np.isclose(x[0], 1.0)
            | np.isclose(x[1], 0.0)
            | np.isclose(x[1], 1.0)
        )

    boundary_dofs = fem.locate_dofs_geometrical(v_space, boundary)
    return fem.dirichletbc(fem.Constant(domain, 0.0), boundary_dofs, v_space)


def solve_heat_problem(v_space, source, bc, ksp_type: str, pc_type: str, ksp_rtol: float):
    """Assemble and solve the FEniCSx weak form."""
    u = ufl.TrialFunction(v_space)
    v = ufl.TestFunction(v_space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    l_form = source * v * ufl.dx

    problem = LinearProblem(
        a,
        l_form,
        bcs=[bc],
        petsc_options_prefix="heat_",
        petsc_options={"ksp_type": ksp_type, "pc_type": pc_type, "ksp_rtol": ksp_rtol},
    )
    uh = problem.solve()
    uh.name = "temperature"
    return uh


def global_min_max(domain, values: np.ndarray) -> tuple[float, float]:
    """Return global min/max across MPI ranks."""
    local_max = np.array(values.max(), dtype=np.float64)
    local_min = np.array(values.min(), dtype=np.float64)
    global_max = np.array(0.0, dtype=np.float64)
    global_min = np.array(0.0, dtype=np.float64)
    domain.comm.Allreduce(local_max, global_max, op=MPI.MAX)
    domain.comm.Allreduce(local_min, global_min, op=MPI.MIN)
    return float(global_min), float(global_max)


def global_source_integral(domain, source) -> float:
    """Return the global integral of the interpolated source."""
    local_integral = np.array(fem.assemble_scalar(fem.form(source * ufl.dx)), dtype=np.float64)
    global_integral = np.array(0.0, dtype=np.float64)
    domain.comm.Allreduce(local_integral, global_integral, op=MPI.SUM)
    return float(global_integral)


def write_xdmf(domain, uh, path: Path) -> None:
    """Write mesh and temperature field to XDMF/HDF5."""
    with io.XDMFFile(domain.comm, path, "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)


def build_summary(
    args: argparse.Namespace,
    domain,
    v_space,
    source_integral: float,
    t_min: float,
    t_max: float,
    xdmf_path: Path,
    png_path: Path,
) -> dict[str, str | int | float]:
    """Collect reproducibility metadata for the solve."""
    png_value = str(png_path.resolve()) if not args.no_png and domain.comm.size == 1 else ""
    return {
        "n": args.n,
        "dofs": v_space.dofmap.index_map.size_global,
        "source_kind": args.source,
        "strength": args.strength,
        "mesh_size_h": 1.0 / args.n,
        "source_integral": source_integral,
        "ksp_type": args.ksp_type,
        "pc_type": args.pc_type,
        "ksp_rtol": args.ksp_rtol,
        "t_min": t_min,
        "t_max": t_max,
        "xdmf_path": str(xdmf_path.resolve()),
        "png_path": png_value,
    }


def write_summary(path: Path, summary: dict[str, str | int | float]) -> None:
    """Write summary metadata as one CSV row."""
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary.keys())
        writer.writeheader()
        writer.writerow(summary)


def print_summary(summary: dict[str, str | int | float]) -> None:
    """Print a compact solve report."""
    print("Solved FEniCSx heat diffusion problem")
    for key, value in summary.items():
        print(f"{key}: {value}")


def save_temperature_png(domain, v_space, uh, path: Path) -> None:
    """Save a quick-look plot for a scalar P1 field on a triangular mesh."""
    if domain.comm.size != 1:
        if domain.comm.rank == 0:
            print("Skipping PNG output: parallel plotting is not implemented.")
        return

    try:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "PNG output requires matplotlib in the FEniCSx environment. "
            "Install it with: conda install -n fenicsx matplotlib"
        ) from exc

    tdim = domain.topology.dim
    num_cells = domain.topology.index_map(tdim).size_local
    cells = np.array(
        [v_space.dofmap.cell_dofs(cell) for cell in range(num_cells)],
        dtype=np.int32,
    )
    coordinates = v_space.tabulate_dof_coordinates()[:, :2]
    values = uh.x.array[: coordinates.shape[0]]

    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1], cells)

    fig, ax = plt.subplots(figsize=(6.0, 5.2), constrained_layout=True)
    plot = ax.tripcolor(triangulation, values, shading="gouraud", cmap="inferno")
    ax.triplot(triangulation, color="white", linewidth=0.2, alpha=0.25)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("FEniCSx steady temperature")
    fig.colorbar(plot, ax=ax, label="T")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    domain, v_space = create_mesh_and_space(args.n)
    source = create_source(v_space, args.source, args.strength)
    bc = create_dirichlet_bc(domain, v_space)
    uh = solve_heat_problem(v_space, source, bc, args.ksp_type, args.pc_type, args.ksp_rtol)

    xdmf_path = args.outdir / "temperature.xdmf"
    write_xdmf(domain, uh, xdmf_path)

    png_path = args.outdir / "temperature.png"
    if not args.no_png:
        save_temperature_png(domain, v_space, uh, png_path)

    t_min, t_max = global_min_max(domain, uh.x.array)
    source_integral = global_source_integral(domain, source)

    if domain.comm.rank == 0:
        summary = build_summary(args, domain, v_space, source_integral, t_min, t_max, xdmf_path, png_path)
        write_summary(args.outdir / "summary.csv", summary)
        print_summary(summary)


if __name__ == "__main__":
    main()
