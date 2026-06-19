from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import solve_fd_heat as fd2
import solve_fd_heat_3D as fd3


def test_2d_laplacian_shape_and_stencil() -> None:
    n = 3
    h = 1.0 / (n + 1)
    a = fd2.build_laplacian_matrix(n).toarray()

    assert a.shape == (n * n, n * n)
    assert np.isclose(a[4, 4], 4.0 / h**2)

    neighbors_of_center = {1, 3, 5, 7}
    for index in neighbors_of_center:
        assert np.isclose(a[4, index], -1.0 / h**2)

    assert np.allclose(a, a.T)


def test_2d_square_source_metadata_and_solution_residual() -> None:
    n = 7
    strength = 100.0
    h = 1.0 / (n + 1)
    _, _, source = fd2.build_source(n, "square", strength)
    a = fd2.build_laplacian_matrix(n)
    rhs = source.reshape(-1)
    u, iterations = fd2.solve_system(a, rhs, "direct", 1e-10)
    residual = fd2.residual_metadata(a, u, rhs)
    metadata = fd2.source_metadata(source, h)

    assert iterations is None
    assert metadata["heat_source_active_points"] == 1
    assert np.isclose(metadata["heat_source_integral"], strength * h**2)
    assert residual["relative_residual"] < 1e-12
    assert float(u.max()) > 0.0


def test_2d_gaussian_source_is_positive_and_centered() -> None:
    n = 9
    _, _, source = fd2.build_source(n, "gaussian", 100.0)

    center = n // 2
    assert np.isclose(source[center, center], source.max())
    assert np.all(source > 0.0)


def test_3d_laplacian_shape_and_stencil() -> None:
    n = 3
    h = 1.0 / (n + 1)
    a = fd3.build_laplacian_matrix(n).toarray()
    center = 13

    assert a.shape == (n**3, n**3)
    assert np.isclose(a[center, center], 6.0 / h**2)

    for index in {4, 10, 12, 14, 16, 22}:
        assert np.isclose(a[center, index], -1.0 / h**2)

    assert np.allclose(a, a.T)


def test_3d_gaussian_small_solve_residual() -> None:
    n = 4
    _, _, _, source = fd3.build_source(n, "gaussian", 100.0)
    a = fd3.build_laplacian_matrix(n)
    rhs = source.reshape(-1)
    u, _ = fd3.solve_system(a, rhs, "direct", 1e-10)
    residual = fd3.residual_metadata(a, u, rhs)

    assert residual["relative_residual"] < 1e-12
    assert float(u.max()) > 0.0
