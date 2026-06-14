# 二维热扩散问题讲义：从 Poisson 方程到 FEniCSx 实现

## 0. 本讲义目标

本讲义基于项目文件 `docs/reference_project_plan.md` 和配套脚本。

我们围绕一个二维芯片稳态热扩散问题，建立从物理模型、数学方程、有限元弱形式到 FEniCSx（新版 FEniCS）代码实现的完整链条。

完成本讲义后，读者应能够理解：

- 稳态热扩散为什么可以写成 Poisson 方程；
- Dirichlet 边界条件的物理意义；
- 有限元方法为什么要引入弱形式；
- 如何把 PDE 离散成稀疏线性方程组；
- 如何用 FEniCSx 编写一个二维热扩散求解器；
- 如何分析网格、热源、边界条件和线性求解器对结果的影响。

### 0.1 讲义使用方式

这份讲义可按两条路径阅读。

第一种是理论优先：

```text
物理背景
  -> PDE 模型
  -> 强形式
  -> 弱形式
  -> FEM 离散
  -> FEniCSx 实现
  -> 数值实验
```

第二种是实操优先：

```text
先运行有限差分脚本
  -> 看温度云图
  -> 回到 Poisson 方程
  -> 理解边界条件
  -> 理解稀疏矩阵
  -> 再进入 FEM 和 FEniCSx
```

若目标是尽快建立数值直觉，建议先运行脚本并观察结果。若目标是打好 FEM 数学基础，建议从弱形式开始逐步推导。

### 0.2 细分学习模块

本讲义可拆成 8 个学习模块：

| 模块 | 主题 | 核心问题 | 对应章节 |
|---|---|---|---|
| M1 | 物理建模 | 为什么芯片热扩散能写成 PDE？ | 第 1-3 节 |
| M2 | PDE 数学模型 | Poisson 方程每一项是什么意思？ | 第 2-3 节 |
| M3 | 弱形式 | 为什么 FEM 不直接解强形式？ | 第 4 节 |
| M4 | FEM 离散 | 基函数、自由度、刚度矩阵如何产生？ | 第 5 节 |
| M5 | 有限差分实操 | 如何先完成一个可执行模拟？ | 第 14 节 |
| M6 | FEniCSx 实现 | 数学弱形式如何变成 FEM 代码？ | 第 6-8、14.7 节 |
| M7 | 数值分析 | 网格、误差、求解器如何影响结果？ | 第 9-11、13 节 |
| M8 | 工程扩展 | 如何进入多材料、瞬态、多物理场？ | 第 12 节、第 16-18 节 |

### 0.3 建议学习节奏

按 6 次课推进时，可采用以下安排：

| 课次 | 内容 | 产出 |
|---|---|---|
| 第 1 次 | 物理背景、几何、热源、边界条件 | 能画出问题示意图 |
| 第 2 次 | Poisson 方程、强形式、弱形式 | 能独立写出弱形式 |
| 第 3 次 | 网格、基函数、矩阵系统 | 能解释 $AU=b$ 的来源 |
| 第 4 次 | 有限差分脚本实操 | 能生成温度云图和实验表 |
| 第 5 次 | FEniCSx 代码结构 | 能把弱形式翻译成 UFL |
| 第 6 次 | 网格实验、求解器、扩展方向 | 能设计自己的数值实验 |

---

## 1. 物理背景：芯片中的稳态导热

考虑一个简化二维芯片区域：

$$
\Omega=[0,1]\times[0,1].
$$

芯片内部可以导热，中间某个小区域持续发热，四周边界连接理想散热器。目标是求解系统达到热平衡后的温度分布，即每个空间位置 $(x,y)$ 上的稳态温度：

$$
T=T(x,y).
$$

在真实工程中，这类问题出现在：

- 半导体芯片热分析；
- 高功耗模块的 hotspot 识别；
- Electronic Design Automation, EDA；
- 热-结构耦合模拟；
- 工业多物理场仿真。

本项目的几何和方程较为简单，但工作流与工业有限元软件高度一致：

$$
\text{物理建模}\rightarrow \text{PDE}\rightarrow \text{弱形式}\rightarrow \text{网格}\rightarrow \text{矩阵装配}\rightarrow \text{线性求解}\rightarrow \text{可视化与分析}.
$$

---

## 2. 数学模型：二维 Poisson 方程

稳态热传导的一般形式可以写为：

$$
-\nabla\cdot(k\nabla T)=q.
$$

其中：

- $T(x,y)$：温度场；
- $k$：热导率；
- $q(x,y)$：体热源强度；
- $-\nabla\cdot(k\nabla T)$：热通量的散度。

如果假设材料均匀且 $k=1$，则方程简化为：

$$
-\nabla^2 T=f(x,y),
$$

也就是二维 Poisson 方程。

其中 Laplace 算子为：

$$
\nabla^2 T=\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}.
$$

因此本项目的模型为：

$$
\begin{cases}
-\nabla^2 T=f(x,y), & (x,y)\in\Omega,\\
T=0, & (x,y)\in\partial\Omega.
\end{cases}
$$

这里 $\partial\Omega$ 表示芯片的边界。

---

## 3. 几何、热源与边界条件

### 3.1 计算区域

我们取单位正方形：

$$
\Omega=[0,1]\times[0,1].
$$

这相当于对实际芯片尺寸进行无量纲化。真实工程中的长度尺度可能是毫米或微米，但只要方程和参数保持一致缩放，数值流程并不会改变。

### 3.2 热源区域

中心发热区定义为：

$$
0.4<x<0.6,\quad 0.4<y<0.6.
$$

热源函数设为：

$$
f(x,y)=
\begin{cases}
100, & 0.4<x<0.6,\ 0.4<y<0.6,\\
0, & \text{其他区域}.
\end{cases}
$$

该区域可理解为芯片中的高功耗单元，或者一个局部 hotspot。

### 3.3 Dirichlet 边界条件

边界条件取：

$$
T=0,\quad (x,y)\in\partial\Omega.
$$

也就是：

$$
T(x,0)=T(x,1)=T(0,y)=T(1,y)=0.
$$

物理含义是芯片四周连接理想散热器，边缘温度被固定。该边界条件虽然理想化，但非常适合作为有限元方法的入门模型。

---

## 4. 从强形式到弱形式

### 4.1 强形式

原始 PDE：

$$
-\nabla^2 T=f
$$

称为强形式。强形式要求 $T$ 至少具有二阶导数；在复杂几何、非光滑热源或分片材料中，这一要求会带来困难。

有限元方法通常不直接离散强形式，而是先转化为弱形式。

### 4.2 引入测试函数

选取任意测试函数 $v$，并要求该函数在边界上为 0：

$$
v=0,\quad \text{on }\partial\Omega.
$$

将 PDE 两边乘以 $v$，并在区域 $\Omega$ 上积分：

$$
\int_\Omega (-\nabla^2 T)v\,dx=\int_\Omega fv\,dx.
$$

### 4.3 分部积分

利用 Green 公式：

$$
\int_\Omega (-\nabla^2 T)v\,dx
=\int_\Omega \nabla T\cdot\nabla v\,dx-\int_{\partial\Omega}\frac{\partial T}{\partial n}v\,ds.
$$

由于 Dirichlet 问题中测试函数在边界上满足 $v=0$，边界积分项消失：

$$
\int_{\partial\Omega}\frac{\partial T}{\partial n}v\,ds=0.
$$

因此得到弱形式：

$$
\int_\Omega \nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx.
$$

### 4.4 弱形式的意义

弱形式的核心好处是：

- 二阶导数变成了一阶导数；
- 可以处理分片线性函数；
- 自然适合在三角形或四边形网格上离散；
- 最终会形成对称正定的稀疏线性系统。

用有限元语言写，就是寻找：

$$
T\in V,\quad T|_{\partial\Omega}=0,
$$

使得对所有测试函数 $v\in V_0$，都有：

$$
a(T,v)=L(v),
$$

其中：

$$
a(T,v)=\int_\Omega \nabla T\cdot\nabla v\,dx,
$$

$$
L(v)=\int_\Omega fv\,dx.
$$

---

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

## 9. 数值结果解读

求解后得到的温度场通常具有以下特征：

- 最高温度出现在中心热源区域附近；
- 温度从中心向边界逐渐降低；
- 四条边界上的温度严格等于 0；
- 等温线大致围绕中心区域分布；
- 热源区域边缘附近温度梯度较大。

网格加密后，温度场通常更平滑，热点边界附近的解析度也会提高。由于热源是分片常数，热源区域边界处存在不连续性，因此局部网格加密通常比全局均匀加密更有效。

---

## 10. 网格研究：精度与计算量

可设置不同网格尺寸：

```python
for n in [16, 32, 64, 128]:
    domain = mesh.create_unit_square(MPI.COMM_WORLD, n, n)
```

建议记录：

- 自由度数量；
- 最大温度 $T_{\max}$；
- 求解时间；
- 可视化结果；
- 稀疏矩阵非零元数量。

一般规律是：

- 网格越细，自由度越多；
- 空间细节解析度越高；
- 线性系统规模快速增长；
- 求解时间和内存占用增加；
- hotspot 附近更需要细网格。

对于二维 $P_1$ 三角形网格，自由度数量大致随 $n^2$ 增长。

---

## 11. 稀疏矩阵与线性求解器

有限元离散后得到：

$$
AU=b.
$$

矩阵 $A$ 的稀疏性来自基函数的局部支撑。一个节点只与邻近节点发生耦合，因此矩阵中大部分元素为 0。

这一性质对大规模仿真尤其重要，因为工业级仿真中的自由度可能达到：

- $10^6$；
- $10^8$；
- 更高规模。

此时不可能使用普通稠密矩阵存储和求解。

对于本问题：

- $A$ 对称；
- $A$ 正定；
- 适合使用 CG；
- 预条件器可以选 Jacobi、ILU、AMG 或 Hypre。

FEniCSx 中可以通过 PETSc 选项控制求解器：

```python
petsc_options={
    "ksp_type": "cg",
    "pc_type": "hypre",
    "ksp_rtol": 1e-10,
}
```

网格较小时，也可先使用直接求解器：

```python
petsc_options={
    "ksp_type": "preonly",
    "pc_type": "lu",
}
```

但直接求解器通常不适合大规模三维工业问题。

---

## 12. 常见扩展

### 12.1 Gaussian 热源

把分片常数热源改成平滑高斯热源：

$$
f(x,y)=A\exp\left(-\frac{(x-x_0)^2+(y-y_0)^2}{2\sigma^2}\right).
$$

代码中可写为：

```python
def gaussian_source(x):
    A = 100.0
    x0, y0 = 0.5, 0.5
    sigma = 0.06
    r2 = (x[0] - x0) ** 2 + (x[1] - y0) ** 2
    return A * np.exp(-r2 / (2.0 * sigma ** 2))
```

Gaussian 热源比方形热源更平滑，所得温度场通常也更光滑。

为了比较单个 Gaussian hotspot 与两个 Gaussian hotspot 的影响，建议保持每个 hotspot 的 $\sigma$ 和峰值强度一致。这样差异主要来自热源数量和位置，而不是热源宽度或单个 hotspot 振幅。

### 12.2 多个 hotspot

可定义多个发热区域：

$$
f(x,y)=f_1(x,y)+f_2(x,y)+\cdots+f_m(x,y).
$$

这种设置更接近芯片上多个功能模块同时发热的情形。

### 12.3 非均匀热导率

如果材料热导率不是常数，模型变为：

$$
-\nabla\cdot(k(x,y)\nabla T)=f.
$$

弱形式变为：

$$
\int_\Omega k(x,y)\nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx.
$$

代码中只需要把：

```python
a = inner(grad(u), grad(v)) * dx
```

改成：

```python
a = k * inner(grad(u), grad(v)) * dx
```

这是多材料热扩散问题的入口。

### 12.4 时间演化热方程

如果考虑瞬态导热，则方程变为：

$$
\frac{\partial T}{\partial t}-\nabla^2T=f.
$$

使用后向 Euler 离散：

$$
\frac{T^{n+1}-T^n}{\Delta t}-\nabla^2T^{n+1}=f.
$$

弱形式为：

$$
\int_\Omega T^{n+1}v\,dx+\Delta t\int_\Omega\nabla T^{n+1}\cdot\nabla v\,dx
=\int_\Omega T^nv\,dx+\Delta t\int_\Omega fv\,dx.
$$

这将进一步涉及质量矩阵、时间步进、稳定性和瞬态仿真。

---

## 13. 建议实验任务

### 实验 1：基础求解

完成单位方形、中心热源、零温边界的稳态求解，并输出：

- 温度最大值；
- 温度最小值；
- 温度云图；
- 等温线图。

### 实验 2：网格加密

分别使用：

$$
n=16,\ 32,\ 64,\ 128
$$

观察：

- $T_{\max}$ 如何变化；
- hotspot 边缘是否更清晰；
- 计算时间如何增长。

### 实验 3：热源强度变化

把热源强度从 100 改成：

$$
50,\ 100,\ 200,\ 500.
$$

观察 $T_{\max}$ 是否近似线性增长。

由于方程是线性的，当边界条件不变且 $f$ 成比例放大时，解 $T$ 也会成比例放大。

### 实验 4：多个发热区

构造两个或三个 hotspot，观察温度场叠加后的形态。

### 实验 5：求解器比较

比较：

- LU；
- CG；
- CG + Jacobi；
- CG + Hypre/AMG。

记录不同网格下的求解时间和迭代次数。

---

## 14. 可执行实操脚本

当前目录下已经新增 `scripts/` 文件夹，包含四份脚本：

- `scripts/solve_fd_heat.py`：有限差分版本，当前环境可直接运行；
- `scripts/run_fd_experiments.py`：批量实验脚本，用于网格加密和热源强度实验；
- `scripts/solve_fenicsx_heat.py`：FEniCSx 版本，使用 `fenicsx` conda 环境运行；
- `scripts/run_fenicsx_heat.sh`：FEniCSx wrapper，用于设置干净的 JIT 编译环境。

### 14.1 为什么先提供有限差分脚本

FEniCSx 是本项目最终要掌握的 FEM 工具。当前 conda base 环境已经具备有限差分脚本所需依赖：

- `numpy`；
- `scipy`；
- `matplotlib`。

FEniCSx 依赖放在独立 conda 环境 `fenicsx` 中，包括：

- `dolfinx`；
- `mpi4py`；
- `petsc4py`。

因此，讲义配套了一套有限差分脚本。其目的不是替代 FEM，而是先完成 PDE、边界条件、稀疏矩阵、线性求解器和可视化的完整验证。进入 FEM 部分时，再使用 `fenicsx` 环境运行 FEniCSx 实现。

### 14.2 有限差分离散

在单位正方形内部取 $n\times n$ 个未知点，网格间距为：

$$
h=\frac{1}{n+1}.
$$

对内部点 $(i,j)$，二维 Laplace 算子使用五点差分：

$$
-\nabla^2 T(x_i,y_j)
\approx
\frac{4T_{i,j}-T_{i-1,j}-T_{i+1,j}-T_{i,j-1}-T_{i,j+1}}{h^2}.
$$

于是离散方程为：

$$
\frac{4T_{i,j}-T_{i-1,j}-T_{i+1,j}-T_{i,j-1}-T_{i,j+1}}{h^2}=f_{i,j}.
$$

边界条件 $T=0$ 通过边界外侧的已知边界值处理。最终仍然得到一个稀疏线性系统：

$$
AU=b.
$$

### 14.3 Kronecker 积形式的推导

脚本中的 `build_laplacian_matrix` 使用如下形式构造二维离散 Laplace 矩阵：

```python
a = (kron(i, t) + kron(t, i)) / h**2
```

数学上对应：

$$
A=\frac{1}{h^2}(I\otimes T+T\otimes I).
$$

这里的 $\otimes$ 是 Kronecker product，中文常称为克罗内克积或张量积。它不是普通内积，而是一种构造块矩阵的运算。

先定义一维二阶差分矩阵：

$$
T=
\begin{bmatrix}
2&-1&0&\cdots&0\\
-1&2&-1&\cdots&0\\
0&-1&2&\cdots&0\\
\vdots&\vdots&\vdots&\ddots&-1\\
0&0&0&-1&2
\end{bmatrix}.
$$

一维中有：

$$
-u''(x_i)\approx \frac{2u_i-u_{i-1}-u_{i+1}}{h^2}.
$$

因此 $T$ 表示未除以 $h^2$ 的一维 $-\partial_{xx}$ 差分模板。

二维内部未知量写成矩阵：

$$
U=
\begin{bmatrix}
u_{1,1}&u_{1,2}&\cdots&u_{1,n}\\
u_{2,1}&u_{2,2}&\cdots&u_{2,n}\\
\vdots&\vdots&\ddots&\vdots\\
u_{n,1}&u_{n,2}&\cdots&u_{n,n}
\end{bmatrix}.
$$

对 $U$ 左乘 $T$：

$$
(TU)_{i,j}=2u_{i,j}-u_{i-1,j}-u_{i+1,j}.
$$

这表示第一个方向上的二阶差分。

对 $U$ 右乘 $T$：

$$
(UT)_{i,j}=2u_{i,j}-u_{i,j-1}-u_{i,j+1}.
$$

这表示第二个方向上的二阶差分。

两个方向相加：

$$
(TU+UT)_{i,j}
=
4u_{i,j}
-u_{i-1,j}
-u_{i+1,j}
-u_{i,j-1}
-u_{i,j+1}.
$$

因此二维离散 Laplace 算子可写为：

$$
-\Delta_h U=\frac{TU+UT}{h^2}.
$$

接下来把矩阵 $U$ 拉平成向量：

$$
\operatorname{vec}(U)
=
(u_{1,1},u_{2,1},\dots,u_{n,1},u_{1,2},\dots,u_{n,n})^T.
$$

采用列优先展开时，有标准恒等式：

$$
\operatorname{vec}(AUB)=(B^T\otimes A)\operatorname{vec}(U).
$$

对第一项：

$$
\operatorname{vec}(TU)
=
\operatorname{vec}(TUI)
=
(I^T\otimes T)\operatorname{vec}(U)
=
(I\otimes T)\operatorname{vec}(U).
$$

对第二项：

$$
\operatorname{vec}(UT)
=
\operatorname{vec}(IUT)
=
(T^T\otimes I)\operatorname{vec}(U).
$$

由于 $T$ 是对称矩阵，$T^T=T$，所以：

$$
\operatorname{vec}(UT)
=
(T\otimes I)\operatorname{vec}(U).
$$

于是：

$$
\operatorname{vec}(TU+UT)
=
(I\otimes T+T\otimes I)\operatorname{vec}(U).
$$

最终得到：

$$
\boxed{
A=\frac{1}{h^2}(I\otimes T+T\otimes I)
}
$$

这正是脚本中：

```python
a = (kron(i, t) + kron(t, i)) / h**2
```

的数学来源。

如果采用行优先展开，$I\otimes T$ 与 $T\otimes I$ 对应的方向会互换，但两项相加后的二维 Laplace 矩阵等价。

三维情况可以从基向量角度直接推广。若一维离散空间的标准基为 $e_i$，则三维网格点 $(i,j,k)$ 可对应张量积基：

$$
e_i\otimes e_j\otimes e_k.
$$

一维差分矩阵 $T$ 作用在哪个方向，就把 $T$ 放在哪个张量因子上；其他方向不变，因此放单位矩阵 $I$。于是三维离散 Laplace 矩阵为：

$$
\boxed{
A_{3D}
=
\frac{1}{h^2}
\left(
T\otimes I\otimes I
+
I\otimes T\otimes I
+
I\otimes I\otimes T
\right)
}
$$

它对应三维七点差分：

$$
-\Delta_h u_{i,j,k}
\approx
\frac{
6u_{i,j,k}
-u_{i-1,j,k}
-u_{i+1,j,k}
-u_{i,j-1,k}
-u_{i,j+1,k}
-u_{i,j,k-1}
-u_{i,j,k+1}
}{h^2}.
$$

代码形式可写成：

```python
A3 = (
    kron(kron(t, i), i)
    + kron(kron(i, t), i)
    + kron(kron(i, i), t)
) / h**2
```

因此，二维和三维的统一规律是：每个空间方向贡献一个一维差分矩阵 $T$，其余方向放单位矩阵 $I$，最后将各方向相加。

有限差分与有限元的关系可概括为：

- 有限差分直接离散微分算子；
- 有限元先写出弱形式，再离散函数空间；
- 两者最终都会形成稀疏矩阵；
- 对规则方形区域，有限差分实现更简洁；
- 对复杂几何和复杂材料，有限元更自然。

### 14.4 单次求解脚本

运行：

```bash
python3 scripts/solve_fd_heat.py --n 64 --source square --solver direct --outdir results/baselines/fd_single
```

参数含义：

- `--n 64`：每个方向 64 个内部网格点；
- `--source square`：中心方形热源；
- `--solver direct`：使用 SciPy 直接稀疏求解器；
- `--outdir results/baselines/fd_single`：输出目录。

输出文件：

- `results/baselines/fd_single/temperature.png`：热源和温度场图；
- `results/baselines/fd_single/summary.csv`：自由度、最大温度、求解时间等摘要；
- `results/baselines/fd_single/solution.npz`：网格、热源和温度数组。

还可尝试 Gaussian 热源：

```bash
python3 scripts/solve_fd_heat.py --n 64 --source gaussian --outdir results/baselines/fd_gaussian
```

或者两个 hotspot：

```bash
python3 scripts/solve_fd_heat.py --n 64 --source two_gaussian --outdir results/baselines/fd_two_gaussian
```

若要测试迭代求解器：

```bash
python3 scripts/solve_fd_heat.py --n 128 --source square --solver cg --outdir results/experiments/fd_cg
```

### 14.5 批量实验脚本

运行：

```bash
python3 scripts/run_fd_experiments.py --mesh-n 16 32 64 128 --strengths 50 100 200 500 --outdir results/experiments/fd_experiments
```

脚本会自动完成两组实验。

第一组是网格加密实验：

- $n=16$；
- $n=32$；
- $n=64$；
- $n=128$。

输出：

- `results/experiments/fd_experiments/mesh_study.csv`；
- `results/experiments/fd_experiments/mesh_study.png`。

第二组是热源强度实验：

- $f=50$；
- $f=100$；
- $f=200$；
- $f=500$。

输出：

- `results/experiments/fd_experiments/strength_study.csv`；
- `results/experiments/fd_experiments/strength_study.png`。

该实验可以验证线性 PDE 的一个重要性质：

$$
f\rightarrow \alpha f
\quad\Rightarrow\quad
T\rightarrow \alpha T.
$$

### 14.6 已验证的运行结果

使用当前环境测试：

```bash
python3 scripts/solve_fd_heat.py --n 32 --source square --solver direct --outdir results/baselines/fd_single_test
```

得到摘要：

```text
n: 32
dofs: 1024
source_kind: square
strength: 100.0
solver: direct
matrix_nnz: 4992
t_min: 0.0033158550721638733
t_max: 1.1278362573560161
```

批量实验：

```bash
python3 scripts/run_fd_experiments.py --mesh-n 16 32 64 --strengths 50 100 200 --outdir results/experiments/fd_experiments_test
```

其中热源强度实验给出：

```text
strength=50.0,  t_max=0.65457619
strength=100.0, t_max=1.3091524
strength=200.0, t_max=2.6183048
```

结果表明，最大温度随热源强度近似线性放大，这与 Poisson 方程的线性性质一致。

### 14.7 FEniCSx 脚本

安装好 FEniCSx 后，可以运行：

```bash
bash scripts/run_fenicsx_heat.sh --n 64 --outdir results/baselines/fenicsx_single
```

这里使用 wrapper 脚本，是为了避免旧 Anaconda 环境变量污染 FEniCSx 的 JIT 编译过程。FEniCSx 会把 UFL 弱形式即时编译成本地代码；如果 `CC`、`CXX`、`CFLAGS` 或 `LDFLAGS` 指向旧环境，可能导致 linker 错误。wrapper 默认把 JIT cache 放到 `/private/tmp`，避免项目目录被 cache 文件打乱。

输出：

- `results/baselines/fenicsx_single/temperature.xdmf`。

该文件可用 ParaView 打开，也可继续用 PyVista 读取并可视化。

FEniCSx 版本的核心仍然是弱形式：

```python
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
l_form = source * v * ufl.dx
```

#### 14.7.1 `solve_fenicsx_heat.py` 代码解读

当前脚本位于：

```text
scripts/solve_fenicsx_heat.py
```

它求解的问题仍然是单位正方形上的稳态热扩散问题：

$$
-\Delta T=f,\qquad T|_{\partial\Omega}=0.
$$

脚本的整体流程是：

```text
解析参数
  -> 创建 mesh
  -> 定义有限元空间
  -> 定义热源函数
  -> 定义 Dirichlet 边界条件
  -> 写出 UFL 弱形式
  -> 交给 PETSc 求解
  -> 输出 XDMF 和 summary.csv
```

第一部分是依赖导入：

```python
from mpi4py import MPI
from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import LinearProblem
import ufl
```

含义如下：

| Python 对象 | 作用 |
|---|---|
| `MPI` | 并行通信接口；即使单进程运行，FEniCSx 也按 MPI 程序组织 |
| `mesh` | 创建几何网格 |
| `fem` | 定义 function space、function、边界条件、积分装配 |
| `io` | 输出 XDMF/HDF5 文件 |
| `LinearProblem` | 自动装配线性系统并调用 PETSc 求解 |
| `ufl` | 用接近数学表达式的方式写弱形式 |

参数部分：

```python
parser.add_argument("--n", type=int, default=64)
parser.add_argument("--strength", type=float, default=100.0)
parser.add_argument("--outdir", type=Path, default=Path("results/baselines/fenicsx_single"))
```

这里 `--n` 在 FEniCSx 中表示每个方向的 mesh subdivision 数。单位正方形被切成 `n` 份，因此 mesh scale 近似为：

$$
h=\frac{1}{n}.
$$

这与有限差分脚本不同。有限差分脚本的 `n` 是内部未知点数量，其网格间距是：

$$
h=\frac{1}{n+1}.
$$

网格创建：

```python
domain = mesh.create_unit_square(
    MPI.COMM_WORLD,
    args.n,
    args.n,
    mesh.CellType.triangle,
)
```

这表示在单位正方形 $\Omega=[0,1]\times[0,1]$ 上创建三角形 mesh。`mesh.CellType.triangle` 表示每个小方格会被切成 triangle cells。FEM 的积分和 basis function 都定义在这些三角形单元上。

有限元空间：

```python
v_space = fem.functionspace(domain, ("Lagrange", 1))
```

这定义了 P1 Lagrange 有限元空间：

$$
V_h=\mathrm{span}\{\phi_1,\phi_2,\dots,\phi_N\}.
$$

其中每个 basis function $\phi_i$ 是分片一次多项式。温度近似写成：

$$
T_h(x,y)=\sum_{j=1}^{N}U_j\phi_j(x,y).
$$

这里的 $U_j$ 是第 $j$ 个自由度，即对应节点上的温度系数。

热源函数：

```python
source = fem.Function(v_space)
```

这里把热源 $f$ 也放在同一个 P1 Lagrange 空间中。也就是说，热源不是作为任意解析表达式直接参与积分，而是先被插值到有限元空间。

当前脚本的热源是 square source：

```python
inside = (x[0] > 0.4) & (x[0] < 0.6) & (x[1] > 0.4) & (x[1] < 0.6)
values[inside] = args.strength
```

数学上对应：

$$
f(x,y)=
\begin{cases}
\mathrm{strength}, & 0.4<x<0.6,\ 0.4<y<0.6,\\
0, & \text{otherwise}.
\end{cases}
$$

然后：

```python
source.interpolate(heat_source)
```

表示把这个函数插值到有限元节点上。粗网格时，square source 的边界不一定与 mesh 对齐，因此离散热源积分 $\int_\Omega f_h dx$ 会随 mesh 变化。

边界条件：

```python
boundary_dofs = fem.locate_dofs_geometrical(v_space, boundary)
bc = fem.dirichletbc(fem.Constant(domain, 0.0), boundary_dofs, v_space)
```

`boundary` 函数识别满足 $x=0$、$x=1$、$y=0$、$y=1$ 的边界节点。`dirichletbc` 将这些自由度固定为 0，即：

$$
T_h=0\quad\text{on }\partial\Omega.
$$

弱形式：

```python
u = ufl.TrialFunction(v_space)
v = ufl.TestFunction(v_space)
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
l_form = source * v * ufl.dx
```

对应数学表达：

$$
a(T_h,v_h)=\int_\Omega \nabla T_h\cdot\nabla v_h\,dx,
$$

$$
L(v_h)=\int_\Omega f_h v_h\,dx.
$$

求解目标是：

$$
a(T_h,v_h)=L(v_h),\qquad \forall v_h\in V_h^0.
$$

这里 $V_h^0$ 表示满足零 Dirichlet 边界条件的测试函数空间。

求解器：

```python
problem = LinearProblem(
    a,
    l_form,
    bcs=[bc],
    petsc_options_prefix="heat_",
    petsc_options={"ksp_type": "cg", "pc_type": "hypre", "ksp_rtol": 1e-10},
)
uh = problem.solve()
```

`LinearProblem` 会自动完成：

```text
UFL form
  -> local element matrix
  -> global sparse matrix
  -> apply Dirichlet BC
  -> PETSc linear solve
```

其中：

| PETSc option | 含义 |
|---|---|
| `ksp_type: cg` | 使用 conjugate gradient method |
| `pc_type: hypre` | 使用 Hypre/AMG preconditioner |
| `ksp_rtol: 1e-10` | 相对收敛阈值 |
| `petsc_options_prefix: "heat_"` | FEniCSx 0.10 要求的 PETSc option namespace |

这里 CG 是合理选择，因为 Poisson 问题在 Dirichlet 条件下得到的 stiffness matrix 通常是 sparse、symmetric、positive definite。

输出：

```python
with io.XDMFFile(domain.comm, xdmf_path, "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(uh)
```

XDMF 文件记录 mesh 和 temperature field，可以用 ParaView 打开。

脚本还计算并输出：

```python
source_integral = fem.assemble_scalar(fem.form(source * ufl.dx))
t_min = uh.x.array.min()
t_max = uh.x.array.max()
```

其中：

- `source_integral` 是离散热源总量 $\int_\Omega f_h dx$；
- `t_min` 用于检查 Dirichlet 边界是否为 0；
- `t_max` 是最常用的温度尺度 sanity check。

最后脚本写出：

```text
temperature.xdmf
temperature.h5
summary.csv
```

其中 `summary.csv` 用来记录可复现实验的核心数值。

该脚本与有限差分脚本的差别在于：

- 有限差分脚本手动构造五点差分矩阵；
- FEniCSx 脚本由框架根据弱形式自动装配有限元矩阵；
- FEniCSx 更适合后续扩展到复杂几何、多材料和三维问题。

### 14.8 FEniCSx 与有限差分结果对比

删除旧 Anaconda 并切换到新的 conda workflow 后，FEniCSx smoke test 可以正常运行。例如：

```bash
bash scripts/run_fenicsx_heat.sh --n 16 --outdir results/tmp_fenicsx_check
```

得到：

```text
dofs = 289
source_integral = 3.515625
t_min = 0
t_max = 1.158638793628342
```

这里 `t_min=0` 来自 Dirichlet 边界条件。FEniCSx 的自由度包括边界节点，因此当网格划分数为 `n` 时，单位正方形上的 P1 Lagrange 元自由度大约是 $(n+1)^2$。例如 `n=16` 时，dofs 为 $17^2=289$。

为了和有限差分结果做公平比较，需要注意二者对 `n` 的定义不同：

- 有限差分脚本中的 `n` 是内部未知点个数，因此网格间距为 $h=1/(n+1)$；
- FEniCSx 脚本中的 `n` 是每个方向的 mesh subdivision 数，因此网格尺度为 $h=1/n$；
- 因此若 FD 使用 `n=64`，对应的 FEniCSx 网格应使用 `n=65`，这样二者的 $h$ 相同。

对齐网格尺度后，得到如下比较：

| FD `n` | FEniCSx `n` | $h$ | FD $\int_\Omega f dx$ | FEniCSx $\int_\Omega f dx$ | FD $T_{\max}$ | FEniCSx $T_{\max}$ | 相对差 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 17 | 0.05882353 | 5.53633218 | 5.53633218 | 1.65596296 | 1.60044130 | -3.35% |
| 32 | 33 | 0.03030303 | 3.30578512 | 3.30578512 | 1.12783626 | 1.11279177 | -1.33% |
| 64 | 65 | 0.01538462 | 4.00000000 | 4.00000000 | 1.30915239 | 1.30521504 | -0.30% |

这个结果说明：

1. 当 $h$ 和离散热源积分对齐后，FD 和 FEniCSx 的结果非常接近；
2. 随着网格加密，$T_{\max}$ 的相对差从约 $3.35\%$ 降到约 $0.30\%$；
3. 二者收敛到同一个连续问题的解，但粗网格上会因为离散方式不同而有可见差异。

从原理上看，有限差分直接离散强形式：

$$
-\Delta T=f,
$$

在每个内部格点上用五点 stencil 近似二阶导数。它要求点值方程近似成立：

$$
\frac{4T_{i,j}-T_{i-1,j}-T_{i+1,j}-T_{i,j-1}-T_{i,j+1}}{h^2}
=f_{i,j}.
$$

FEniCSx/FEM 则先写弱形式：

$$
\int_\Omega \nabla T_h\cdot\nabla v_h\,dx
=
\int_\Omega f v_h\,dx,
$$

然后在 P1 Lagrange 有限元空间中求解。它不是逐点要求 PDE 成立，而是要求 residual 对所有测试函数 $v_h$ 的投影为零：

$$
\int_\Omega \left(-\Delta T_h-f\right)v_h\,dx=0.
$$

因此，FD 和 FEM 的矩阵来源不同：

- FD 矩阵来自五点差分 stencil；
- FEM 刚度矩阵来自 $\int_\Omega \nabla\phi_i\cdot\nabla\phi_j\,dx$；
- FD 的未知量是内部格点温度；
- FEM 的未知量是有限元节点上的系数，包括边界节点，但 Dirichlet 边界节点被固定为零。

两者的物理问题相同，离散路径不同。网格加密后，如果热源、边界条件和几何区域一致，二者应趋于同一个连续解。

---

## 15. 常见错误与调试方法

### 15.1 边界条件没有正确施加

现象：

- 边界温度不是 0；
- 解整体偏移；
- 矩阵可能奇异。

检查：

- `boundary` 函数是否正确识别四条边；
- `locate_dofs_geometrical` 是否找到了边界自由度；
- `bcs=[bc]` 是否传入 `LinearProblem`。

### 15.2 热源没有插值成功

现象：

- 温度场全为 0；
- 最大温度为 0。

检查：

- `f.interpolate(heat_source)` 是否执行；
- `heat_source` 返回数组长度是否正确；
- hotspot 判断条件是否过窄或写错。

### 15.3 网格太粗

现象：

- hotspot 区域形状不清楚；
- 最大温度不稳定；
- 等温线粗糙。

解决：

- 增大 `n`；
- 在 hotspot 附近做局部加密；
- 使用更高阶有限元。

### 15.4 求解器不收敛

现象：

- 迭代次数较多；
- 残差下降缓慢；
- 程序运行时间过长。

解决：

- 先使用 LU 验证模型是否正确；
- 使用 CG + AMG；
- 检查边界条件；
- 检查矩阵是否仍然对称正定。

---

## 16. 本项目的知识地图

二维热扩散项目连接了多个重要方向：

```text
热传导物理
    ↓
Poisson / Heat Equation
    ↓
弱形式与变分方法
    ↓
有限元基函数与网格
    ↓
稀疏矩阵装配
    ↓
PETSc 线性求解器
    ↓
科学计算与工业多物理场仿真
```

本项目的价值不只是写出一个可运行程序，而是建立一条完整的认知链条：

$$
\text{Physics}\rightarrow \text{PDE}\rightarrow \text{FEM}\rightarrow \text{Sparse Solver}\rightarrow \text{Simulation Workflow}.
$$

这正是从数据分析进入工业数值物理模拟的重要第一步。

---

## 17. 进一步学习路线

建议按以下顺序继续学习：

1. 理解 Poisson 方程和 Dirichlet 边界条件；
2. 手动推导一次弱形式，尤其是分部积分；
3. 运行并验证 FEniCSx 基础代码；
4. 改变网格尺寸，观察数值结果；
5. 改变热源，理解线性系统响应；
6. 研究稀疏矩阵和求解器；
7. 扩展到非均匀材料；
8. 扩展到瞬态热方程；
9. 最后考虑多物理场或 AI surrogate model。

推荐资料：

- FEniCSx Tutorial: https://jsdokken.com/dolfinx-tutorial/
- MIT 18.336 Numerical Methods for PDEs: https://ocw.mit.edu/courses/18-336-numerical-methods-for-partial-differential-equations-spring-2009/
- PETSc: https://petsc.org
- Zienkiewicz, The Finite Element Method

---

## 18. 分模块讲授方案

本节将前面的内容进一步细分为可逐节学习、逐节验收的课程单元。每个模块都包含：

- 学习目标；
- 必备前置知识；
- 讲授重点；
- 实操任务；
- 检查问题；
- 本节交付物。

### 18.1 模块 M1：从物理问题到计算问题

学习目标：

- 理解芯片热扩散问题的工程背景；
- 能够说明温度场、热源、边界散热器分别代表什么；
- 能够把真实物理问题抽象成二维单位方形上的计算问题。

前置知识：

- 基本函数概念；
- 二维坐标系；
- 温度、热源、边界的物理直觉。

讲授重点：

- 芯片内部导热；
- 中心 hotspot；
- 四周理想散热；
- 稳态温度分布。

实操任务：

1. 手动画出单位方形区域。
2. 标出中心热源区域。
3. 标出四条边界上的 $T=0$。
4. 用一句话解释最高温度为什么应出现在热源附近。

检查问题：

- 为什么该问题可以先用二维模型而不是三维模型？
- 为什么边界温度固定为 0 是一种理想化？
- 如果热源移到左上角，温度云图会如何变化？

本节交付物：

- 一张问题示意图；
- 一段对物理模型的文字描述。

### 18.2 模块 M2：Poisson 方程建模

学习目标：

- 理解稳态热扩散方程；
- 能解释 $-\nabla^2 T=f$ 中每一项的含义；
- 能区分未知量、已知量、区域和边界条件。

前置知识：

- 一阶、二阶导数；
- 偏导数；
- 梯度和散度的基本直觉。

讲授重点：

- 温度场 $T(x,y)$ 是未知函数；
- 热源 $f(x,y)$ 是已知输入；
- Laplace 算子刻画扩散和平滑；
- Dirichlet 边界条件给出边界温度。

核心公式：

$$
\begin{cases}
-\nabla^2T=f, & (x,y)\in\Omega,\\
T=0, & (x,y)\in\partial\Omega.
\end{cases}
$$

实操任务：

1. 写出单位方形区域 $\Omega$。
2. 写出中心方形热源函数。
3. 写出四条边界上的温度条件。
4. 判断方程中哪些是输入，哪些是输出。

检查问题：

- 为什么热源越大，温度通常越高？
- 如果 $f=0$ 且边界 $T=0$，解应该是什么？
- 如果只固定一条边界，问题会发生什么变化？

本节交付物：

- 一页数学模型说明。

### 18.3 模块 M3：强形式到弱形式

学习目标：

- 理解强形式和弱形式的区别；
- 能完成乘测试函数、积分、分部积分这三步；
- 理解二阶导数如何转化为一阶导数。

前置知识：

- 积分；
- 分部积分；
- 向量点乘；
- 边界项的基本概念。

讲授重点：

- 强形式要求更高光滑性；
- 测试函数用于“探测”方程是否成立；
- 分部积分把导数从 $T$ 转移到 $v$；
- Dirichlet 问题中测试函数在边界为 0。

推导主线：

$$
-\nabla^2T=f
$$

乘以测试函数 $v$ 并积分：

$$
\int_\Omega (-\nabla^2T)v\,dx=\int_\Omega fv\,dx.
$$

分部积分后得到：

$$
\int_\Omega \nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx.
$$

实操任务：

1. 手写一遍弱形式推导。
2. 标出哪一步使用了分部积分。
3. 标出边界项为什么消失。
4. 把 $a(T,v)$ 和 $L(v)$ 分别写出来。

检查问题：

- 为什么弱形式只需要一阶导数？
- 测试函数为什么要在 Dirichlet 边界上为 0？
- 如果是 Neumann 边界，边界项还会消失吗？

本节交付物：

- 一份完整弱形式推导。

### 18.4 模块 M4：有限元离散

学习目标：

- 理解网格、单元、节点、自由度；
- 理解基函数展开 $T_h=\sum_jT_j\phi_j$；
- 能解释刚度矩阵和右端向量的来源。

前置知识：

- 线性代数基础；
- 矩阵乘法；
- 多项式插值的直觉。

讲授重点：

- 用三角形网格近似区域；
- 用局部基函数近似未知函数；
- 每个节点温度是一个自由度；
- 弱形式离散后得到 $AU=b$。

核心公式：

$$
A_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dx,
$$

$$
b_i=\int_\Omega f\phi_i\,dx.
$$

实操任务：

1. 画一个由三角形组成的粗网格。
2. 在一个节点上画出局部支撑区域。
3. 写出 $T_h$ 的基函数展开。
4. 解释为什么矩阵 $A$ 是稀疏的。

检查问题：

- 为什么一个节点只和附近节点耦合？
- 一次 Lagrange 元的自由度在哪里？
- 网格变细时，自由度数量如何变化？

本节交付物：

- 一张 FEM 离散示意图；
- 一段对稀疏矩阵来源的解释。

### 18.5 模块 M5：有限差分实操

学习目标：

- 运行并验证当前环境可执行的热扩散脚本；
- 理解五点差分矩阵；
- 能生成温度云图、摘要 CSV 和实验曲线。

前置知识：

- Python 基础；
- NumPy 数组；
- Matplotlib 基础；
- 稀疏矩阵的概念。

讲授重点：

- 为什么可先用有限差分建立数值直觉；
- 五点差分如何离散 $-\nabla^2T$；
- `scipy.sparse` 如何构造矩阵；
- 直接求解器和 CG 求解器的区别。

实操任务：

运行单次求解：

```bash
python3 scripts/solve_fd_heat.py --n 64 --source square --outdir results/baselines/fd_single
```

运行 Gaussian 热源：

```bash
python3 scripts/solve_fd_heat.py --n 64 --source gaussian --outdir results/baselines/fd_gaussian
```

运行批量实验：

```bash
python3 scripts/run_fd_experiments.py --mesh-n 16 32 64 128 --strengths 50 100 200 500 --outdir results/experiments/fd_experiments
```

检查问题：

- 为什么热源强度翻倍时，最大温度也近似翻倍？
- 为什么不同网格下 $T_{\max}$ 会变化？
- 为什么直接求解器在大规模问题中不一定合适？

本节交付物：

- 一张温度云图；
- 一份 `summary.csv`；
- 一张网格实验曲线；
- 一张热源强度实验曲线。

### 18.6 模块 M6：FEniCSx 实现

学习目标：

- 理解 FEniCSx 中 mesh、function space、trial function、test function 的含义；
- 能把弱形式翻译成 UFL 代码；
- 能看懂 `solve_fenicsx_heat.py`。

前置知识：

- 弱形式；
- FEM 基函数；
- Python 函数和对象。

讲授重点：

- `mesh.create_unit_square` 对应计算区域；
- `fem.functionspace` 对应有限元空间；
- `ufl.TrialFunction` 和 `ufl.TestFunction` 对应未知函数和测试函数；
- `inner(grad(u), grad(v))*dx` 对应左端双线性型；
- `source*v*dx` 对应右端线性型。

实操任务：

安装 FEniCSx 后运行：

```bash
bash scripts/run_fenicsx_heat.sh --n 64 --outdir results/baselines/fenicsx_single
```

检查问题：

- FEniCSx 为什么不需要手动写矩阵元素？
- UFL 代码和弱形式之间如何对应？
- FEniCSx 相比有限差分更适合哪些问题？

本节交付物：

- 一个 `temperature.xdmf` 文件；
- 一段 FEniCSx 弱形式代码解释。

### 18.7 模块 M7：数值实验与误差直觉

学习目标：

- 理解网格加密对结果的影响；
- 理解热源形状对温度场的影响；
- 初步建立误差、收敛和计算成本的直觉。

前置知识：

- 已经运行并验证有限差分脚本；
- 知道自由度和网格大小的关系。

讲授重点：

- 粗网格可能低估或高估局部热点；
- 分片常数热源边缘附近梯度较大；
- Gaussian 热源更平滑；
- 自由度增长会带来计算成本增长。

实操任务：

1. 固定热源，改变网格。
2. 固定网格，改变热源强度。
3. 比较 square、gaussian、two_gaussian 三种热源。
4. 记录每次实验的 $T_{\max}$ 和运行时间。

检查问题：

- 为什么结果不是简单地随着网格单调变化？
- 热源不连续时，哪里最需要细网格？
- 如何判断一个网格已经“够用”？

本节交付物：

- 一份实验记录表；
- 一段对结果变化的解释。

### 18.8 模块 M8：工程扩展路线

学习目标：

- 知道如何从基础 Poisson 问题扩展到真实工程问题；
- 理解多材料、瞬态、各向异性和 AI surrogate 的入口；
- 形成后续项目路线。

前置知识：

- 已理解稳态 Poisson 模型；
- 已运行并验证至少一个数值脚本。

讲授重点：

- 多材料：把常数热导率变成 $k(x,y)$；
- 瞬态：加入时间导数 $\partial T/\partial t$；
- 各向异性：把热导率变成张量；
- AI surrogate：学习从参数到温度场的映射。

实操任务：

1. 把热源改成两个 hotspot。
2. 把方形热源改成 Gaussian 热源。
3. 设计一个多材料芯片区域。
4. 写出瞬态热方程的后向 Euler 弱形式。

检查问题：

- 多材料界面上温度和热流应该满足什么条件？
- 瞬态问题为什么需要初始条件？
- surrogate model 的输入和输出应该如何设计？

本节交付物：

- 一个后续扩展方案；
- 一个新的数值实验设计。

---

## 19. 分层练习题

### 19.1 基础题

1. 写出二维稳态热扩散问题的 PDE、区域和边界条件。

    答：
    - 二维PDE写作$\nabla T(x,y)=f(x,y)$，即$\frac{\partial^2 T(x,y)}{\partial x^2}+\frac{\partial^2 T(x,y)}{\partial y^2}=f(x,y)$.
    - 区域其实取决于实际情况，如上述例子则是一个方形；区域常记做$\Omega$，其边界常记为$\partial\Omega$.
    - 边界条件以上述例子为参照，则是区域边界为恒温0度，比如边界和冷却系统联通。

    批改：
    - 区域和边界条件的理解是对的。
    - PDE 这里需要改正。热扩散稳态 Poisson 方程不是 $\nabla T=f$，而是
      $$
      -\nabla^2T=f.
      $$
      其中 $\nabla T$ 是梯度，是一个向量；$\nabla^2T$ 或 $\Delta T$ 才是 Laplace 算子，是标量。
    - 如果展开，应写为
      $$
      -\left(\frac{\partial^2T}{\partial x^2}+\frac{\partial^2T}{\partial y^2}\right)=f(x,y).
      $$
    - 本题较完整的答案是：
      $$
      \begin{cases}
      -\nabla^2T=f, & (x,y)\in\Omega,\\
      T=0, & (x,y)\in\partial\Omega.
      \end{cases}
      $$

2. 解释 $T(x,y)$、$f(x,y)$、$\Omega$、$\partial\Omega$ 的含义。

    答：
    - $T(x,y)$ 表示在位置点$(x,y)$的温度（也可叫做温度场）。
    - $f(x,y)$ 表示热源在空间上的影响函数。
    - $\Omega$ 表示研究的区域，$\partial\Omega$ 表示区域的边界。

    批改：
    - $T(x,y)$、$\Omega$、$\partial\Omega$ 的解释正确。
    - $f(x,y)$ 建议说得更物理一些：它是热源项或单位热导率下的体热源强度，不是温度本身。
    - 更准确表述：
      $$
      f(x,y)\text{ 表示单位面积上的热源强度，决定区域内部哪里在发热、发热多强。}
      $$
    - 在无量纲模型 $-\nabla^2T=f$ 中，$f$ 的量纲相当于温度除以长度平方。

3. 为什么中心热源会导致中心附近温度更高？

    答：中心热源产生的热量会向边缘低温区域扩散，这个符合热力学第二定律（熵增），越靠近边缘温度越低，越靠近中心温度越高。

    批改：
    - 直觉是正确的：热源处持续输入热量，边界又被固定为低温，因此温度峰值会出现在热源附近。
    - “符合热力学第二定律”可以作为物理直觉，但在本 PDE 语境中，更直接的解释是：Poisson 方程把热源 $f>0$ 转化为温度场的局部抬升，Dirichlet 边界 $T=0$ 则把边界固定为散热端。
    - 可以补一句最大值原理的直觉：对于内部正热源和零边界，解通常在内部形成正温度峰值，并向边界衰减。

4. 为什么边界条件 $T=0$ 会让温度向边界降低？

    答：答案同问题3.

    批改：
    - 这题和第 3 题相关，但最好单独说明边界条件的作用。
    - 更准确的回答是：$T=0$ 把边界强制固定为最低参考温度。内部热源使区域内温度升高，而边界始终保持 0，因此从内部高温区域到边界会形成温度梯度。
    - 物理上，零温边界相当于理想散热器；数学上，Dirichlet 边界给解施加了固定边界值，使温度场必须在边界处降到 0。

### 19.2 推导题

1. 从 $-\nabla^2T=f$ 出发，推导弱形式。

    答：
    从PDE出发，首先定义基函数线性空间 $span\{v_i\}$，且如果真实解可以近似为$T\approx T_h=\sum_i T_iv_i$，那么近似解的残差$R=-\nabla^2T_h-f$在基函数线性空间中为0，即残差$R$和所有基函数的内积为0，$(R,v_i)=\iint_\Omega R\,v_i\,dS=0,\,\forall v_i$. 这样我们就可以写出以下方程：
    $$
    -\sum_i T_i \iint_\Omega (\nabla^2v_i)v_j\,dS-\iint_\Omega fv_j\,dS=0
    $$
    因为$(\nabla^2v_i)v_j=\nabla((\nabla v_i)v_j)-\nabla v_i\nabla v_j$，且第一项表示向量的散度在2维面积分，由高斯公式得到
    $$
    \iint_\Omega \nabla((\nabla v_i)v_j)\,dS=\int_{\partial\Omega} (\nabla v_i)v_j\cdot d\vec{l}
    $$
    其中$d\vec{l}$为方向向外的线积分元。
    无需进一步推导，因为基函数构造时需满足边界条件，即$v_i$在边界$\partial\Omega$恒等于0，所以积分为0，则方程可以写成：
    $$
    \sum_i T_i \iint_\Omega \nabla v_i \nabla v_j\,dS-\iint_\Omega fv_j\,dS=0
    $$
    至此我们通过构造基函数以及应用边界条件，成功将二次方程转化为一次方程。

    批改：
    - 你的主线是对的：从近似解 $T_h$ 出发，考虑残差，并要求残差在测试空间方向上为零。这就是 Galerkin 方法的核心直觉。
    - 需要修正几个关键符号。
    - 第一，测试函数空间建议写作 $V_h=\operatorname{span}\{\phi_i\}$，不要写成 $span\{v_i\}$ 后又把 $v_i$ 同时当测试函数和基函数。通常令
      $$
      T_h=\sum_jT_j\phi_j,\quad v_h=\phi_i.
      $$
    - 第二，散度恒等式应写成
      $$
      \nabla\cdot(v_j\nabla v_i)
      =
      \nabla v_j\cdot\nabla v_i+v_j\nabla^2v_i.
      $$
      因此
      $$
      (\nabla^2v_i)v_j
      =
      \nabla\cdot(v_j\nabla v_i)-\nabla v_i\cdot\nabla v_j.
      $$
      你写的方向是对的，但要把散度符号写成 $\nabla\cdot(\cdot)$，并注意是 $v_j\nabla v_i$。
    - 第三，边界项中消失的是测试函数 $v_j$。边界项应类似
      $$
      \int_{\partial\Omega}v_j\nabla v_i\cdot n\,ds.
      $$
      由于 Dirichlet 问题的测试函数满足 $v_j=0$ on $\partial\Omega$，所以边界项为 0。
    - 第四，最后一句“二次方程转化为一次方程”建议改为“把对 $T$ 的二阶导数要求转化为一阶导数要求”。这里不是代数意义上的二次方程，而是二阶微分算子变成了一阶梯度内积。
    - 更标准的最终弱形式是：
      $$
      \int_\Omega \nabla T_h\cdot\nabla v_h\,dS
      =
      \int_\Omega fv_h\,dS,
      \quad \forall v_h\in V_h.
      $$
      取 $v_h=\phi_i$，得到
      $$
      \sum_jT_j\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dS
      =
      \int_\Omega f\phi_i\,dS.
      $$

2. 说明分部积分后边界项为什么消失。

    答：边界项消失的原因是因为运用高斯公式后，积分含有基函数$v_i$，其在边界上很等于0。

    批改：
    - 方向正确，但有两个小问题。
    - “很等于 0”应为“恒等于 0”。
    - 更精确地说，边界项含有的是测试函数，而不是所有基函数都自动在边界为 0。对于齐次 Dirichlet 边界条件，测试函数空间取为
      $$
      V_0=\{v: v=0\text{ on }\partial\Omega\}.
      $$
      因此边界项
      $$
      \int_{\partial\Omega}v\frac{\partial T}{\partial n}\,ds
      $$
      因为 $v=0$ 而消失。
    - 如果是非齐次 Dirichlet 条件，通常会先做 lifting，把未知函数拆成满足边界值的部分加上零边界的测试空间；测试函数仍然在边界为 0。

3. 写出有限元离散后的矩阵元素 $A_{ij}$ 和 $b_i$。

    答：这里我们可以讲1.中最后公式的$ij$对换，以保证index的对应，那么对应的矩阵元素是：
    $$
    A_{ij}=\sum_i T_j \iint_\Omega \nabla v_j \nabla v_i\,dS\\
    b_i=\iint_\Omega fv_i\,dS
    $$

    批改：
    - 这里 $b_i$ 基本正确。
    - $A_{ij}$ 需要修正：矩阵元素本身不包含未知系数 $T_j$，也不包含外面的求和。
    - 正确写法是
      $$
      A_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dS,
      $$
      $$
      b_i=\int_\Omega f\phi_i\,dS.
      $$
    - 线性系统整体才是
      $$
      \sum_jA_{ij}T_j=b_i.
      $$
    - 如果用你的 $v_i$ 记号，也应写为
      $$
      A_{ij}=\int_\Omega \nabla v_j\cdot\nabla v_i\,dS.
      $$
      其中 $i$ 对应测试函数，$j$ 对应 trial basis function。

4. 解释矩阵 $A$ 为什么是稀疏矩阵。
    
    答：因为基函数定义为当前格点局部的影响函数，且不同基函数在线性空间中正交，这意味着矩阵 $A$ 是稀疏。不过我觉得我的解释还不够严密。

    批改：
    - 你觉得“不够严密”是对的；关键问题在于“正交”这句话不准确。
    - Lagrange 有限元基函数一般不是正交基。矩阵稀疏不是因为基函数正交，而是因为基函数具有局部支撑。
    - 矩阵元素是
      $$
      A_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,dS.
      $$
      如果 $\phi_i$ 和 $\phi_j$ 的支撑区域没有重叠，则它们的梯度也没有共同非零区域，因此
      $$
      A_{ij}=0.
      $$
    - 对于 $P_1$ Lagrange 元，一个节点的基函数只在围绕该节点的一圈三角形上非零。因此节点 $i$ 只和相邻节点产生矩阵耦合，远处节点对应的矩阵元素为 0。
    - 更准确的答案是：$A$ 稀疏来自有限元基函数的局部支撑和网格的局部连接结构，而不是来自基函数正交。

### 19.3 代码题

1. 运行 `solve_fd_heat.py`，生成一张温度云图。

    results/baselines/fd_single/temperature.png

    批改：
    - 完成。建议在答案里顺手记录 `summary.csv` 中的核心数值，例如 `t_max`、`relative_residual` 和 `heat_source_integral`。这样不仅有图，也有可复现实验记录。

2. 把 `--source square` 改成 `--source gaussian`，比较温度场。

    results/baselines/fd_gaussian/temperature.png

    批改：
    - 完成。建议补一句观察结论：Gaussian 热源的热源分布和温度分布都更平滑；square 热源在热源边界处有不连续跳变，因此温度梯度变化更明显。
    - 也可以比较两个 summary：
      - square source 的 `heat_source_integral` 约为 4；
      - gaussian source 的 `heat_source_integral` 约为 2.262；
      因此二者的 `t_max` 不能只从峰值强度 100 来比较，还要考虑总热源强度。

3. 把 `--n` 从 32 改到 128，记录最大温度和运行时间。

    - n=32:
        - t_max: 0.7637532352767834
        - elapsed_seconds: 0.001334624997980427
    - n=128:
        - t_max: 0.7693165760544349
        - elapsed_seconds: 0.0456337910000002
    最大温度几乎相等，但是n=128用时约为n=32的34倍。

    批改：
    - 记录方式很好，已经包含了关键的 accuracy/runtime 对比。
    - 这里你记录的是 Gaussian source 的结果，而不是 square source；从数值看，$T_{\max}$ 已经从 0.7638 收敛到 0.7693，说明 Gaussian 情况下网格收敛很平滑。
    - “最大温度几乎相等”可以更精确地说：相对变化约为
      $$
      \frac{0.7693-0.7638}{0.7693}\approx 0.7\%.
      $$
    - 用时增长明显是合理的，因为自由度从
      $$
      32^2=1024
      $$
      增加到
      $$
      128^2=16384,
      $$
      自由度增加了 16 倍；直接求解器的时间通常会比线性增长更快。

4. 使用 `--solver cg`，观察迭代求解器是否正常收敛。

    - t_max: 0.7693165760544942
    - elapsed_seconds: 0.08487920900006429

    奇怪的是cg求解反而慢，可能是因为查分情况下矩阵并不复杂，解析求解更快。

    批改：
    - 观察正确：在这个二维、规模不大的问题中，CG 可能比直接求解慢。
    - 这里“解析求解”建议改成“直接求解”。`spsolve` 不是解析解，而是稀疏直接线性代数求解。
    - CG 是否正常收敛，最好同时记录：
      - `cg_iterations`
      - `relative_residual`
    - 对小型 2D 问题，直接求解器常常很快；CG 的优势通常出现在更大规模、尤其是 3D 问题中，并且通常需要合适的预条件器。
    - 因此结论可以写成：当前 2D 测试规模下，直接求解器更快；但随着自由度增加，直接求解器内存和 fill-in 成本会迅速上升，迭代法更有优势。

### 19.4 分析题

1. 为什么热源强度翻倍时，温度场也会翻倍？

    根据泊松方程$-\nabla^2 T=f$, 热源在空间上扩散程度与$f$成正比，所以热源翻倍，温度场也翻倍。或着可以直接理解为温度是热源函数在空间上的2重积分，积分的结果和热源强度常数有线性关系。

    批改：
    - 结论正确，原因是 Poisson 方程和边界条件都是线性的。
    - “温度是热源函数在空间上的 2 重积分”这个说法有一定直觉意义，但不够严谨。更准确地说，温度是热源经过 Poisson 算子的逆作用得到的结果：
      $$
      T=(-\Delta)^{-1}f.
      $$
    - 由于 $(-\Delta)^{-1}$ 是线性算子，所以
      $$
      f\to \alpha f
      \quad\Rightarrow\quad
      T\to \alpha T.
      $$
    - 离散形式也一样：
      $$
      AU=b.
      $$
      若 $b\to \alpha b$，则
      $$
      U\to \alpha U.
      $$

2. 为什么 Gaussian 热源的温度场通常比方形热源更平滑？

    因为高斯函数在全平面可导，而方形热源在边缘并不连续。由柏松方程的弱形式我们得知温度的梯度与热源函数相关，如果热源连续，那温度也会更加平滑。

    批改：
    - 主要判断正确：Gaussian 热源平滑，square 热源有跳变，所以 Gaussian 产生的温度场通常更平滑。
    - “柏松”应写为“泊松”。
    - “弱形式得知温度的梯度与热源函数相关”这个表述可以更精确：弱形式中
      $$
      \int_\Omega \nabla T\cdot\nabla v\,dx=\int_\Omega fv\,dx
      $$
      说明热源 $f$ 通过整体积分关系驱动温度梯度，而不是点对点地直接决定 $\nabla T$。
    - 更严谨的说法是：椭圆型方程有平滑效应，右端项越平滑，解通常越平滑；方形热源的边界不连续，会在温度场梯度中留下更明显的局部变化。

3. 为什么热点边缘附近需要更细网格？

    热点边缘由于热源只在内部，天然形成温度梯度，更多的格点有助于更好的模拟温度的非线性行为。

    批改：
    - 方向正确：热点边缘附近确实需要更细网格。
    - 建议把“非线性行为”改掉。当前 Poisson 方程是线性的；需要细网格不是因为方程非线性，而是因为热源在边缘处变化快或不连续，导致温度梯度变化较强。
    - 更准确说法：热点边缘处 $f$ 从高值突然变为 0，局部解的曲率和梯度变化更明显。更细的网格可以更准确地解析这些局部变化，降低几何采样误差和离散误差。

4. 直接求解器和迭代求解器分别适合什么规模的问题？

    直接求解器适用于低自由度矩阵运算，迭代求解器适合高自由度。

    批改：
    - 结论方向正确，但可以更完整。
    - 直接求解器适合小到中等规模问题，优点是稳健、精度高、调参少；缺点是内存消耗大，稀疏矩阵分解时会产生 fill-in。
    - 迭代求解器适合大规模稀疏系统，尤其是 3D FEM/FD 问题；优点是内存更可控，能利用稀疏矩阵乘法；缺点是需要收敛判据和预条件器。
    - 对当前 Poisson 问题，矩阵是稀疏、对称、正定的，因此 CG 是合适的迭代方法。
    - 更完整答案可以写成：二维小规模实验中直接求解器通常更快；当自由度增大到几十万、几百万，尤其是三维问题时，迭代求解器配合预条件器更适合。

### 19.5 扩展题

1. 设计一个有两个不同热导率区域的芯片模型。
2. 写出非均匀热导率问题的弱形式。
3. 写出瞬态热方程的后向 Euler 离散。
4. 设计一个 AI surrogate model 的输入参数和输出温度场。
