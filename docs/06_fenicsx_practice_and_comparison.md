# FEniCSx 实操与 FD 对比

### 14.7 FEniCSx 脚本

安装好 FEniCSx 后，可以运行：

```bash
bash scripts/run_fenicsx_heat.sh --n 64 --outdir results/baselines/fenicsx_single
```

这里使用 wrapper 脚本，是为了避免旧 Anaconda 环境变量污染 FEniCSx 的 JIT 编译过程。FEniCSx 会把 UFL 弱形式即时编译成本地代码；如果 `CC`、`CXX`、`CFLAGS` 或 `LDFLAGS` 指向旧环境，可能导致 linker 错误。wrapper 默认把 JIT cache 放到 `/private/tmp`，避免项目目录被 cache 文件打乱。

输出：

- `results/baselines/fenicsx_single/temperature.xdmf`：给 ParaView 使用的 mesh + field metadata；
- `results/baselines/fenicsx_single/temperature.h5`：XDMF 引用的实际数值数据；
- `results/baselines/fenicsx_single/temperature.png`：用 matplotlib 生成的 quick-look 图；
- `results/baselines/fenicsx_single/summary.csv`：核心数值 metadata。

其中 `.xdmf` 和 `.h5` 需要放在一起。ParaView 打开 `.xdmf` 时会自动读取 `.h5`；`temperature.png` 用于在 VS Code 或 GitHub 中快速检查结果。

FEniCSx 版本的核心仍然是弱形式：

```python
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
l_form = source * v * ufl.dx
```

#### 14.7.1 UFL notation 小附录

UFL 的 notation 要分成两层理解：

```text
pointwise expression
integration measure
```

例如：

```python
ufl.inner(ufl.grad(u), ufl.grad(v))
```

只表示每个空间点上的向量内积：

$$
\nabla u(x)\cdot\nabla v(x).
$$

在二维中就是：

$$
\frac{\partial u}{\partial x}
\frac{\partial v}{\partial x}
+
\frac{\partial u}{\partial y}
\frac{\partial v}{\partial y}.
$$

它只是 integrand，还不是区域积分。乘上：

```python
* ufl.dx
```

才表示在 cell domain $\Omega$ 上积分：

$$
\int_\Omega \nabla u\cdot\nabla v\,dx.
$$

因此：

| UFL | 数学含义 |
|---|---|
| `ufl.grad(u)` | $\nabla u$ |
| `ufl.inner(a, b)` | pointwise inner product $a(x)\cdot b(x)$ |
| `* ufl.dx` | 对 cell domain 积分 |
| `* ufl.ds` | 对外边界 facet 积分 |
| `* ufl.dS` | 对内部 facet 积分 |

所以：

```python
ufl.inner(u, v) * ufl.dx
```

对应函数空间里的 $L^2$ inner product：

$$
\int_\Omega uv\,dx.
$$

而：

```python
ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
```

对应 Poisson/FEM 弱形式中的 gradient inner product：

$$
\int_\Omega \nabla u\cdot\nabla v\,dx.
$$

#### 14.7.2 `solve_fenicsx_heat.py` 代码解读

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
  -> 输出 XDMF/HDF5、PNG 和 summary.csv
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
parser.add_argument("--source", choices=SOURCE_CHOICES, default="square")
parser.add_argument("--strength", type=float, default=100.0)
parser.add_argument("--outdir", type=Path, default=Path("results/baselines/fenicsx_single"))
parser.add_argument("--ksp-type", default="cg")
parser.add_argument("--pc-type", default="hypre")
parser.add_argument("--ksp-rtol", type=float, default=1e-10)
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

热源函数被单独拆成 `evaluate_source`：

```python
def evaluate_source(x, kind, strength):
    ...
```

它的输入 `x` 是 FEniCSx 传入的 interpolation points 坐标。在二维中，shape 近似为：

```python
x.shape == (2, num_points)
```

其中 `x[0]` 是所有点的横坐标，`x[1]` 是所有点的纵坐标。函数返回 shape 为 `(num_points,)` 的一维数组，表示每个 interpolation point 处的热源值。

当前脚本支持：

```text
square
gaussian
two_gaussian
two_hotspots
```

square source 数学上对应：

$$
f(x,y)=
\begin{cases}
\mathrm{strength}, & 0.4<x<0.6,\ 0.4<y<0.6,\\
0, & \text{otherwise}.
\end{cases}
$$

然后：

```python
source = fem.Function(v_space)
source.interpolate(lambda x: evaluate_source(x, kind, strength))
```

表示把这个函数 interpolate 到有限元节点上。

这里有一个重要 caveat：真实的 square source 是不连续函数，但当前代码把它 interpolate 到连续的 P1 Lagrange 空间。也就是说，FEniCSx 实际参与积分的是 $f_h$，不是精确的分片常数 $f$。对 P1 Lagrange 元来说，`evaluate_source(x, kind, strength)` 在节点上取值，然后由 basis functions 在单元内线性延拓。

因此：

- 如果某个 triangle 的顶点一部分在 source 区域内、一部分在区域外，该单元内部的 $f_h$ 会线性过渡，而不是突然跳变；
- 粗网格时，square source 的边界不一定与 mesh 对齐，因此离散热源积分 $\int_\Omega f_h dx$ 会随 mesh 变化；
- 网格加密后，这个 nodal interpolation 会逐渐逼近原始 discontinuous source，但 sharp edge 附近仍是主要误差来源。

如果后续想更忠实地表示分片常数热源，可以改用 cellwise 的 `DG0` source，或在 UFL 积分表达式中使用 conditional expression。

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
temperature.png
summary.csv
```

其中 `temperature.xdmf`/`temperature.h5` 用于 ParaView，`temperature.png` 是 matplotlib 生成的 quick-look 图，`summary.csv` 用来记录可复现实验的核心数值。

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

---

[返回目录](README.md) | [上一章：有限差分实操](05_finite_difference_practice.md) | [下一章：调试、知识地图与学习路线](07_debugging_and_roadmap.md)
