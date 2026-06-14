# Project 1 —— 使用 FEniCS 实现二维芯片热扩散模拟

# 项目目标

本项目旨在作为进入以下领域的第一步：

- Scientific Computing（科学计算）
- Partial Differential Equations（偏微分方程）
- Finite Element Method（有限元方法，FEM）
- Sparse Linear Algebra（稀疏线性代数）
- Industrial Multiphysics Simulation（工业多物理场模拟）

项目核心目标：

> 使用 FEniCS 构建一个完整的二维稳态热扩散模拟器。

这个项目数学上并不复杂，但它的结构与真实工业仿真流程高度一致，因此非常适合作为进入工业 scientific computing 的第一个项目。

其背景广泛存在于：

- 半导体热分析
- Electronic Design Automation（EDA）
- 芯片散热设计
- 热-结构耦合分析
- 多物理场模拟

---

# 1. 物理背景

考虑一个简化二维芯片：

- 芯片内部能够导热
- 中间某区域持续发热
- 芯片边缘连接理想散热器

目标是计算：

## 稳态温度分布

即：

系统达到热平衡后，芯片内部每一点的温度是多少。

---

# 2. 数学模型

系统满足二维 Poisson 方程：

$$-\nabla^2 T = f(x,y)$$


其中：

- $T(x,y)$：温度场
- $f(x,y)$：热源分布

这个方程描述：

> 热量如何在材料内部扩散。

---

# 3. 几何模型

建立一个二维方形区域：

- 区域范围：$[0,1] \times [0,1]$

在中心区域定义一个发热区：

例如：

- $0.4 < x < 0.6$
- $0.4 < y < 0.6$

这个区域可以理解为：

- 芯片中的高功耗单元
- hotspot
- 发热晶体管区域

---

# 4. 热源定义

定义热源：

$$
f(x,y)=
\begin{cases}
100 & \text{热源区域内部} \\
0 & \text{其他区域}
\end{cases}
$$

后续你可以继续扩展：

- Gaussian 热源
- 多个 hotspot
- 非均匀功耗分布

---

# 5. 边界条件

采用 Dirichlet 边界条件：

$$
T=0
$$

即：

芯片四周边界温度固定为 0。

物理意义：

- 芯片边缘连接理想散热器
- 热量可以被无限快带走

这是最经典的 FEM 边界条件之一。

---

# 6. 数值方法

本项目采用：

# Finite Element Method（FEM）

你需要理解 FEM 的完整流程：

PDE

→ Weak Form（弱形式）

→ Basis Expansion（基函数展开）

→ Element Integral（单元积分）

→ Global Sparse Matrix（全局稀疏矩阵）

→ Linear Solver（线性求解器）

这是整个工业 FEM 软件的核心逻辑。

---

# 7. 本项目最重要的学习内容

# 7.1 Weak Formulation（弱形式）

这是 FEM 最核心的思想。

你需要研究：

- 为什么 FEM 要引入 weak form
- 为什么需要 test function
- 为什么二阶导数会消失
- integration by parts 的意义
- variational formulation

核心目标：

真正理解：

$$-\nabla^2 T = f$$

如何变成积分形式。

这是进入 FEM 世界最关键的一步。

---

# 7.2 Basis Function（基函数）

你需要理解：

$$
T(x,y)=\sum_i T_i \phi_i(x,y)
$$

其中：

- $\phi_i$：basis / shape function
- $T_i$：节点自由度

重点学习：

- local support
- linear basis function
- degree of freedom (DOF)
- interpolation

---

# 7.3 Mesh（网格）

研究：

- triangular mesh
- mesh refinement
- local resolution
- 高梯度区域为什么需要更细 mesh

建议尝试：

## 第一版

uniform mesh

## 第二版

在 hotspot 周围加密 mesh

这是真实工业仿真中的重要思想。

---

# 7.4 Sparse Matrix（稀疏矩阵）

理解为什么 FEM 天然产生：

# 稀疏矩阵

重点研究：

- sparsity pattern
- 邻接耦合
- memory efficiency
- sparse storage

工业 scientific computing 中：

solver 的核心几乎都围绕 sparse matrix 展开。

---

# 7.5 Linear Solver（线性求解器）

学习：

- Conjugate Gradient（CG）
- GMRES
- preconditioning

理解：

- 为什么 direct solver 不适合超大系统
- 为什么 iterative solver 是 HPC 核心

---

# 8. 软件栈

# 主框架

## FEniCS

官方网站：

https://fenicsproject.org

推荐原因：

- 与数学表达极其接近
- 对物理背景非常友好
- FEM 学习成本低
- 可以快速建立 PDE → FEM 的直觉

---

# 9. 推荐学习资料

# PDE 与数值方法

## MIT 18.336 Numerical Methods for PDEs

https://ocw.mit.edu/courses/18-336-numerical-methods-for-partial-differential-equations-spring-2009/

推荐原因：

- 数学严谨
- 非常适合物理背景
- 注重 numerical intuition

---

# FEM 入门

## FEniCS Tutorial

https://jsdokken.com/dolfinx-tutorial/

建议：

边学边实现。

---

# FEM 理论书籍

## The Finite Element Method

作者：

Zienkiewicz

这是 FEM 经典教材之一。

---

# Sparse Solver

## PETSc

https://petsc.org

工业和 HPC 中非常重要。

---

# 10. 项目 Milestone

# Milestone 1

安装：

- Python
- FEniCS
- numpy
- matplotlib

验证：

- 官方 Poisson example 可以运行

---

# Milestone 2

建立：

- geometry
- mesh

输出：

- mesh visualization

---

# Milestone 3

实现：

- heat source
- boundary condition
- Poisson equation

求解：

- 温度场

---

# Milestone 4

实现可视化：

- temperature heatmap
- contour plot
- mesh plot

---

# Milestone 5

研究 mesh：

- coarse mesh
- fine mesh
- local refinement

观察：

- 精度变化
- runtime 变化

---

# Milestone 6

研究 solver：

- direct solver
- iterative solver
- CG convergence
- sparse matrix size

---

# 11. 后续扩展方向

完成基础版本后，可以继续扩展。

---

# Extension A —— 时间演化

求解：

$$
\frac{\partial T}{\partial t}-\nabla^2 T=f
$$

研究：

- time stepping
- stability
- transient simulation

---

# Extension B —— 多材料系统

不同区域使用不同 thermal conductivity。

研究：

- material interface
- discontinuous coefficient

---

# Extension C —— 各向异性导热

研究：

- tensor conductivity
- direction-dependent diffusion

---

# Extension D —— AI Surrogate Model

训练神经网络：

输入：

- 热源参数
- boundary condition

输出：

- temperature field

这是未来工业 scientific AI 的重要方向：

- Physics-informed AI
- surrogate simulation
- neural operator

---

# 12. 项目最终目标

完成本项目后，你应该能够理解：

# Physics

- 热扩散 PDE
- boundary condition

# Mathematics

- weak form
- FEM discretization

# Numerical Methods

- sparse matrix
- iterative solver

# Scientific Computing

- 完整 PDE simulation workflow
- 工业 FEM 思维

# Software Engineering

- 仿真代码组织
- visualization
- numerical analysis

---

# 13. 本项目真正的意义

这个项目不仅仅是：

> “写一个热扩散程序”。

它真正的目标是：

建立：

- PDE
- FEM
- sparse solver
- HPC
- multiphysics
- scientific computing

之间的完整认知链条。

这是从：

> 数据分析

真正进入：

> 工业数值物理模拟

最关键的第一步。

