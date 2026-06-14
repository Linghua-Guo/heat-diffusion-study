# 有限差分实操

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

---

[返回目录](README.md) | [上一章：数值结果、网格研究与扩展](04_numerical_intuition.md) | [下一章：FEniCSx 实操与 FD 对比](06_fenicsx_practice_and_comparison.md)
