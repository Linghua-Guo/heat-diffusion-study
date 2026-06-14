# 有限元离散与 FEniCSx 基础

## 5. 有限元离散

### 5.1 网格

将区域 $\Omega$ 划分成若干小三角形单元：

$$
\Omega \approx \bigcup_e K_e.
$$

在每个三角形单元 $K_e$ 上，使用局部多项式近似温度场。最常见的入门选择是一次 Lagrange 元，也称为 $P_1$ 元。

### 5.2 基函数展开

用有限个基函数 $\phi_i(x,y)$ 近似温度：

$$
T_h(x,y)=\sum_{j=1}^{N}T_j\phi_j(x,y).
$$

其中：

- $T_h$：离散后的近似温度场；
- $\phi_j$：第 $j$ 个有限元基函数；
- $T_j$：第 $j$ 个节点上的未知温度；
- $N$：自由度数量。

对于 $P_1$ 三角形单元，基函数具有局部支撑性质：一个节点的基函数只在与该节点相邻的单元上非零。由此得到的全局矩阵具有稀疏结构。

### 5.3 离散线性系统

把 $T_h=\sum_j T_j\phi_j$ 代入弱形式，并取测试函数 $v=\phi_i$，得到：

$$
\sum_{j=1}^{N}T_j\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dx
=\int_\Omega f\phi_i\,dx.
$$

写成矩阵形式：

$$
AU=b.
$$

其中：

$$
A_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dx,
$$

$$
b_i=\int_\Omega f\phi_i\,dx,
$$

$$
U=(T_1,T_2,\dots,T_N)^T.
$$

矩阵 $A$ 通常称为刚度矩阵。对于该 Poisson 问题，在 Dirichlet 边界条件处理后，$A$ 是稀疏、对称、正定矩阵，因此适合使用共轭梯度法 CG 和预条件方法求解。

---

## 6. FEniCSx 实现思路

FEniCSx 是新版 FEniCS 的主要实现。其优势是代码形式非常接近数学表达。数学中的弱形式：

$$
\int_\Omega \nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx
$$

离散后，在 FEniCSx/UFL 中可写为：

```python
a = inner(grad(u), grad(v)) * dx
L = f * v * dx
```

这里代码变量 `u` 表示 trial function，对应数学中的离散未知函数 $T_h$；代码变量 `v` 表示 test function，对应离散测试函数 $v_h$。代码中的 `V` 是 FEniCSx 变量名，对应数学上的有限元空间 $V_h$。

对应关系如下：

| 数学对象 | FEniCSx 对象 |
|---|---|
| 计算区域 $\Omega$ | `mesh` / `domain` |
| 有限元空间 $V_h$ | `fem.functionspace(domain, ...)`，代码变量常写作 `V` |
| 离散未知函数 $T_h$ | `ufl.TrialFunction(V)`，代码变量常写作 `u` |
| 离散测试函数 $v_h$ | `ufl.TestFunction(V)`，代码变量常写作 `v` |
| 热源 $f$ | `fem.Function` / `fem.Constant` / expression |
| 双线性型 $a(T_h,v_h)$ | UFL form，代码变量常写作 `a` |
| 线性型 $L(v_h)$ | UFL form，代码变量常写作 `L` 或 `l_form` |
| 边界条件 $T_h=0$ | `fem.dirichletbc(...)` |
| 稀疏线性系统 $AU=b$ | PETSc matrix/vector，由 `LinearProblem` 装配和求解 |
| 解向量 $U=(T_1,\dots,T_N)^T$ | 求解后的 `Function`，例如 `uh` 及其系数数组 `uh.x.array` |

---

## 7. FEniCSx 示例代码

下面的代码使用 FEniCSx 求解二维芯片稳态热扩散问题。

```python
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

from dolfinx import mesh, fem, plot
from dolfinx.fem.petsc import LinearProblem
import ufl


# 1. 创建二维单位方形网格
domain = mesh.create_unit_square(
    MPI.COMM_WORLD,
    64,
    64,
    cell_type=mesh.CellType.triangle,
)


# 2. 定义有限元函数空间：一次 Lagrange 元
V = fem.functionspace(domain, ("Lagrange", 1))


# 3. 定义热源 f
f = fem.Function(V)

def heat_source(x):
    values = np.zeros(x.shape[1], dtype=np.float64)
    inside = (
        (x[0] > 0.4) & (x[0] < 0.6) &
        (x[1] > 0.4) & (x[1] < 0.6)
    )
    values[inside] = 100.0
    return values

f.interpolate(heat_source)


# 4. 定义 Dirichlet 边界条件 T = 0
def boundary(x):
    return (
        np.isclose(x[0], 0.0) |
        np.isclose(x[0], 1.0) |
        np.isclose(x[1], 0.0) |
        np.isclose(x[1], 1.0)
    )

boundary_dofs = fem.locate_dofs_geometrical(V, boundary)
bc = fem.dirichletbc(
    value=fem.Constant(domain, 0.0),
    dofs=boundary_dofs,
    V=V,
)


# 5. 写出弱形式
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx


# 6. 求解线性系统
problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options={
        "ksp_type": "cg",
        "pc_type": "hypre",
        "ksp_rtol": 1e-10,
    },
)

uh = problem.solve()


# 7. 输出基本数值信息
temperature = uh.x.array
print("Number of DOFs:", V.dofmap.index_map.size_global)
print("T_min:", temperature.min())
print("T_max:", temperature.max())


# 8. 简单可视化
try:
    import pyvista

    topology, cell_types, geometry = plot.vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid.point_data["Temperature"] = temperature
    grid.set_active_scalars("Temperature")

    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=False)
    plotter.view_xy()
    plotter.show()

except ImportError:
    print("PyVista is not installed. Use XDMF output or install pyvista for visualization.")
```

这段代码对应的核心流程是：

1. 创建单位方形网格；
2. 定义 $P_1$ 有限元空间；
3. 定义中心局部热源；
4. 对四周边界施加 $T=0$；
5. 写出弱形式；
6. 调用 PETSc 求解稀疏线性系统；
7. 查看自由度数量、温度最小值和最大值；
8. 可视化温度场。

---

## 8. 代码与数学的逐行对应

### 8.1 网格

```python
domain = mesh.create_unit_square(MPI.COMM_WORLD, 64, 64)
```

对应：

$$
\Omega=[0,1]\times[0,1].
$$

`64, 64` 表示在两个方向上划分网格。网格越细，自由度越多，数值解通常越精细，但计算代价也越高。

### 8.2 有限元空间

```python
V = fem.functionspace(domain, ("Lagrange", 1))
```

对应：

$$
V_h=\text{span}\{\phi_1,\phi_2,\dots,\phi_N\}.
$$

`("Lagrange", 1)` 表示一次 Lagrange 有限元，也就是每个三角形单元上的线性函数。

### 8.3 热源

```python
inside = (
    (x[0] > 0.4) & (x[0] < 0.6) &
    (x[1] > 0.4) & (x[1] < 0.6)
)
values[inside] = 100.0
```

对应：

$$
f(x,y)=
\begin{cases}
100, & 0.4<x<0.6,\ 0.4<y<0.6,\\
0, & \text{otherwise}.
\end{cases}
$$

### 8.4 边界条件

```python
bc = fem.dirichletbc(value=fem.Constant(domain, 0.0), dofs=boundary_dofs, V=V)
```

对应：

$$
T=0,\quad \text{on }\partial\Omega.
$$

### 8.5 弱形式

```python
a = inner(grad(u), grad(v)) * dx
L = f * v * dx
```

对应：

$$
\int_\Omega \nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx.
$$

---

---

[返回目录](README.md) | [上一章：从强形式到弱形式](02_weak_form.md) | [下一章：数值结果、网格研究与扩展](04_numerical_intuition.md)
