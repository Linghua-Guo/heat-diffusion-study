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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--strength", type=float, default=100.0)
    parser.add_argument("--outdir", type=Path, default=Path("results/baselines/fenicsx_single"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    domain = mesh.create_unit_square(MPI.COMM_WORLD, args.n, args.n, mesh.CellType.triangle)
    v_space = fem.functionspace(domain, ("Lagrange", 1))

    source = fem.Function(v_space)

    def heat_source(x):
        values = np.zeros(x.shape[1], dtype=np.float64)
        inside = (x[0] > 0.4) & (x[0] < 0.6) & (x[1] > 0.4) & (x[1] < 0.6)
        values[inside] = args.strength
        return values

    source.interpolate(heat_source)

    def boundary(x):
        return (
            np.isclose(x[0], 0.0)
            | np.isclose(x[0], 1.0)
            | np.isclose(x[1], 0.0)
            | np.isclose(x[1], 1.0)
        )

    boundary_dofs = fem.locate_dofs_geometrical(v_space, boundary)
    bc = fem.dirichletbc(fem.Constant(domain, 0.0), boundary_dofs, v_space)

    u = ufl.TrialFunction(v_space)
    v = ufl.TestFunction(v_space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    l_form = source * v * ufl.dx

    problem = LinearProblem(
        a,
        l_form,
        bcs=[bc],
        petsc_options_prefix="heat_",
        petsc_options={"ksp_type": "cg", "pc_type": "hypre", "ksp_rtol": 1e-10},
    )
    uh = problem.solve()
    uh.name = "temperature"

    xdmf_path = args.outdir / "temperature.xdmf"
    with io.XDMFFile(domain.comm, xdmf_path, "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)

    local_values = uh.x.array
    local_max = np.array(local_values.max(), dtype=np.float64)
    local_min = np.array(local_values.min(), dtype=np.float64)
    global_max = np.array(0.0, dtype=np.float64)
    global_min = np.array(0.0, dtype=np.float64)
    domain.comm.Allreduce(local_max, global_max, op=MPI.MAX)
    domain.comm.Allreduce(local_min, global_min, op=MPI.MIN)

    local_source_integral = np.array(fem.assemble_scalar(fem.form(source * ufl.dx)), dtype=np.float64)
    global_source_integral = np.array(0.0, dtype=np.float64)
    domain.comm.Allreduce(local_source_integral, global_source_integral, op=MPI.SUM)

    if domain.comm.rank == 0:
        summary = {
            "n": args.n,
            "dofs": v_space.dofmap.index_map.size_global,
            "source_kind": "square",
            "strength": args.strength,
            "mesh_size_h": 1.0 / args.n,
            "source_integral": float(global_source_integral),
            "t_min": float(global_min),
            "t_max": float(global_max),
            "xdmf_path": str(xdmf_path.resolve()),
        }
        with (args.outdir / "summary.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary.keys())
            writer.writeheader()
            writer.writerow(summary)

        print("Solved FEniCSx heat diffusion problem")
        print(f"n: {args.n}")
        print(f"dofs: {v_space.dofmap.index_map.size_global}")
        print(f"source_integral: {float(global_source_integral)}")
        print(f"t_min: {float(global_min)}")
        print(f"t_max: {float(global_max)}")
        print(f"output: {xdmf_path.resolve()}")


if __name__ == "__main__":
    main()
