# 调试、知识地图与学习路线

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

---

[返回目录](README.md) | [上一章：FEniCSx 实操与 FD 对比](06_fenicsx_practice_and_comparison.md) | [下一章：分层练习题与批改](08_exercises.md)
