# 分层 Q&A

本章把前面练习题整理成 Q&A。所有问题按主题编码，答案只保留标准表述，便于复习和检索。

## 问题索引

| 编码 | 主题 |
|---|---|
| Q01-Q09 | 物理模型、PDE、边界条件与温度尺度 |
| Q10-Q17 | 弱形式、分部积分、Galerkin 与残差 |
| Q18-Q31 | 有限元空间、basis functions、stiffness matrix |
| Q32-Q40 | 有限差分、Kronecker product、residual |
| Q41-Q49 | CG、solver、tolerance 与数值尺度 |
| Q50-Q67 | FEniCSx 实现、UFL、source interpolation、ParaView |
| Q68-Q78 | 代码细节、结果分析、FD/FEM 对比 |
| Q79-Q86 | 非均匀热导率、瞬态问题、surrogate model |

## A. 物理模型与 PDE

### Q01. 二维稳态热扩散问题的标准 PDE 是什么？

标准模型是 Poisson equation：

$$
\begin{cases}
-\Delta T=f, & (x,y)\in\Omega,\\
T=0, & (x,y)\in\partial\Omega.
\end{cases}
$$

其中 $T$ 是温度场，$f$ 是热源项，$\Omega$ 是计算区域，$\partial\Omega$ 是区域边界。注意 $\nabla T$ 是梯度向量，$\Delta T=\nabla^2T$ 才是 Laplacian。

### Q02. 展开写出二维 Laplacian 后，方程是什么？

二维 Laplacian 为：

$$
\Delta T=
\frac{\partial^2T}{\partial x^2}
+
\frac{\partial^2T}{\partial y^2}.
$$

因此强形式为：

$$
-\left(
\frac{\partial^2T}{\partial x^2}
+
\frac{\partial^2T}{\partial y^2}
\right)=f(x,y).
$$

### Q03. $T(x,y)$ 的物理含义是什么？

$T(x,y)$ 是位置 $(x,y)$ 处的温度，也称为温度场。它是未知函数，需要通过 PDE 和边界条件求解。

### Q04. $f(x,y)$ 的物理含义是什么？

$f(x,y)$ 是热源项，表示单位面积上的发热强度。在无量纲模型 $-\Delta T=f$ 中，$f$ 的量纲相当于温度除以长度平方。

### Q05. $\Omega$ 和 $\partial\Omega$ 分别是什么？

$\Omega$ 是研究的空间区域。例如单位正方形：

$$
\Omega=[0,1]\times[0,1].
$$

$\partial\Omega$ 是区域边界，即正方形的四条边。

### Q06. 为什么中心热源会让中心附近温度更高？

内部热源 $f>0$ 持续向区域输入热量，而边界被固定为 $T=0$，相当于理想散热器。因此温度通常在热源附近抬升，并向边界衰减。数学上，Poisson 方程把正热源转化为内部的正温度响应。

### Q07. 为什么 $T=0$ 边界会让温度向边界降低？

Dirichlet 条件 $T=0$ 强制边界温度固定为参考低温。内部热源使区域内部升温，但边界不能升温，因此从内部到边界形成温度梯度。

### Q08. 这里的 $T=0$ 是否表示真实绝对零度？

不一定。通常 $T=0$ 是参考温度或相对温度。例如边界与恒温冷却系统相连，可把冷却温度设为零点。

### Q09. 如果热源强度翻倍，温度场为什么也翻倍？

Poisson 方程和齐次 Dirichlet 边界条件都是线性的。形式上：

$$
T=(-\Delta)^{-1}f.
$$

所以：

$$
f\to\alpha f
\quad\Rightarrow\quad
T\to\alpha T.
$$

离散后也是：

$$
AU=b,
$$

若 $b\to\alpha b$，则 $U\to\alpha U$。

## B. 弱形式与 Galerkin 思想

### Q10. 为什么需要弱形式？

强形式 $-\Delta T=f$ 要求 $T$ 至少有二阶导数。复杂几何、不连续热源、分片材料会让这个要求过强。弱形式通过分部积分把二阶导数转移为一阶梯度内积，使解空间要求降低，也更适合有限元离散。

### Q11. 从 $-\Delta T=f$ 推导弱形式的第一步是什么？

取测试函数 $v$，两边乘以 $v$ 并在区域上积分：

$$
\int_\Omega (-\Delta T)v\,dx
=
\int_\Omega fv\,dx.
$$

这里 $dx$ 表示面积积分元，在二维中也可写作 $dS$。

### Q12. 分部积分公式在这里如何使用？

使用恒等式：

$$
\nabla\cdot(v\nabla T)
=
\nabla v\cdot\nabla T
+
v\Delta T.
$$

整理得：

$$
-(\Delta T)v
=
\nabla T\cdot\nabla v
-
\nabla\cdot(v\nabla T).
$$

积分后：

$$
\int_\Omega (-\Delta T)v\,dx
=
\int_\Omega \nabla T\cdot\nabla v\,dx
-
\int_{\partial\Omega}v\nabla T\cdot n\,ds.
$$

### Q13. 为什么分部积分后的边界项消失？

对齐次 Dirichlet 问题，测试函数取自零边界空间：

$$
V_0=\{v: v=0\text{ on }\partial\Omega\}.
$$

因此边界项

$$
\int_{\partial\Omega}v\nabla T\cdot n\,ds
$$

因为 $v=0$ 而消失。

### Q14. 最终弱形式是什么？

寻找 $T\in V$，使得对任意测试函数 $v\in V_0$：

$$
\int_\Omega\nabla T\cdot\nabla v\,dx
=
\int_\Omega fv\,dx.
$$

齐次 Dirichlet 条件下，通常 $V=V_0$。

### Q15. 残差角度如何理解弱形式？

令近似解为 $T_h$，残差为：

$$
R=-\Delta T_h-f.
$$

Galerkin 方法不是要求 $R$ 在每个点都为零，而是要求残差对所有测试函数方向正交：

$$
\int_\Omega Rv_h\,dx=0,
\quad \forall v_h\in V_h.
$$

经过分部积分后，就得到有限元弱形式。

### Q16. 测试函数和基函数是否相同？

在标准 Galerkin FEM 中，测试空间和 trial space 通常取同一个有限元空间。因此测试函数可以取为基函数：

$$
v_h=\phi_i.
$$

但概念上二者角色不同：trial function 用来表示未知解，test function 用来投影残差。

### Q17. 非齐次 Dirichlet 条件下，测试函数还为零吗？

是。若边界条件为 $T=g$，通常把解写成：

$$
T=\tilde T+T_g,
$$

其中 $T_g$ 满足边界值，$\tilde T$ 在边界为零。测试函数仍取零边界空间，所以分部积分边界项仍因测试函数为零而消失。

## C. 有限元离散

### Q18. 有限元空间 $V_h$ 是什么？

$V_h$ 是有限维函数空间，由有限元 basis functions 张成：

$$
V_h=\mathrm{span}\{\phi_1,\phi_2,\dots,\phi_N\}.
$$

这里 $N$ 是自由度数量。对 P1 Lagrange 元，通常每个 mesh vertex 对应一个自由度。

### Q19. span 是什么意思？

$\mathrm{span}$ 表示所有线性组合的集合：

$$
\mathrm{span}\{\phi_1,\dots,\phi_N\}
=
\left\{
\sum_{j=1}^{N}U_j\phi_j
\right\}.
$$

因此 $V_h$ 中的任意函数都可以写成 basis functions 的线性组合。

### Q20. 近似温度场 $T_h$ 如何表示？

有限元近似写作：

$$
T_h(x,y)=\sum_{j=1}^{N}U_j\phi_j(x,y).
$$

$U_j$ 是第 $j$ 个自由度的系数。对 P1 Lagrange 元，$U_j$ 就是第 $j$ 个节点处的温度值。

### Q21. 节点是什么意思？

节点是 mesh 划分后用于定义自由度的位置。对三角形 P1 Lagrange 元，节点就是 triangle vertices。更高阶 Lagrange 元还会有边中点、单元内部点等额外节点。

### Q22. 自由度是什么意思？

自由度是离散未知量的数量和位置。对 P1 Lagrange 温度场，每个非边界节点有一个未知温度系数。Dirichlet 边界节点被固定，不再作为自由未知量参与求解。

### Q23. Lagrange basis function 有什么特征？

第 $i$ 个 Lagrange basis function 满足 nodal property：

$$
\phi_i(x_j)=\delta_{ij}.
$$

也就是在自己的节点为 1，在其他节点为 0。

### Q24. 为什么 P1 Lagrange basis 是一次 polynomial？

P1 表示单元上使用一次多项式。对 triangle cell，basis function 在每个三角形内是线性函数。许多 P1 basis functions 拼接起来，形成全局连续、分片线性的近似温度场。

### Q25. 基函数定义在整个区域还是局部？

全局 basis function $\phi_i$ 是定义在整个区域 $\Omega$ 上的函数，但它只在与节点 $i$ 相邻的一圈单元上非零。这叫局部支撑。

### Q26. 有限元矩阵元素 $A_{ij}$ 和 $b_i$ 是什么？

将

$$
T_h=\sum_jU_j\phi_j
$$

代入弱形式，并取 $v_h=\phi_i$，得到：

$$
\sum_jU_j\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,dx
=
\int_\Omega f\phi_i\,dx.
$$

因此：

$$
A_{ij}=\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,dx,
$$

$$
b_i=\int_\Omega f\phi_i\,dx.
$$

### Q27. 为什么 stiffness matrix 是 sparse？

稀疏性来自 basis functions 的局部支撑。若 $\phi_i$ 和 $\phi_j$ 的支撑不重叠，则：

$$
A_{ij}=\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,dx=0.
$$

所以一个节点只和相邻节点耦合，远处节点对应的矩阵元素为零。

### Q28. 稀疏性是否来自 basis functions 正交？

不是。Lagrange basis functions 一般不正交。矩阵稀疏来自局部支撑和 mesh 的局部连接结构。

### Q29. stiffness matrix 的英文是什么？

刚度矩阵是 stiffness matrix。对 Poisson/FEM 问题，它通常指：

$$
A_{ij}=\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,dx.
$$

### Q30. stiffness matrix 为什么通常是 symmetric？

因为：

$$
A_{ij}
=
\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,dx
=
\int_\Omega\nabla\phi_i\cdot\nabla\phi_j\,dx
=
A_{ji}.
$$

### Q31. stiffness matrix 为什么 positive definite？

对非零函数

$$
v_h=\sum_iU_i\phi_i,
$$

有：

$$
U^TAU=\int_\Omega |\nabla v_h|^2\,dx.
$$

在零 Dirichlet 边界条件下，若 $v_h$ 非零，则该积分大于零。因此矩阵是 positive definite。

## D. 有限差分与 Kronecker 结构

### Q32. FD 中 $n$ 和 $h$ 的关系是什么？

有限差分脚本中 $n$ 表示每个方向的内部未知点数。单位区间有两个边界点，因此：

$$
h=\frac{1}{n+1}.
$$

### Q33. 为什么不是 $h=1/n$？

因为 $n$ 是内部未知点数量，不是切分的区间数量。若有 $n$ 个内部点，加上两个边界点，总共有 $n+2$ 个点，因此区间数是 $n+1$。

### Q34. 二维五点差分格式是什么？

对内部点 $(i,j)$：

$$
-\Delta T_{i,j}
\approx
\frac{
4T_{i,j}
-T_{i-1,j}
-T_{i+1,j}
-T_{i,j-1}
-T_{i,j+1}
}{h^2}.
$$

### Q35. 二维 Laplacian 矩阵为什么可以写成 Kronecker product？

二维算子可拆成两个方向：

$$
-\Delta=-\partial_{xx}-\partial_{yy}.
$$

离散后，$x$ 方向作用在每一行上，$y$ 方向作用在每一列上。若一维矩阵为 $T$，单位矩阵为 $I$，则二维矩阵为：

$$
A=\frac{1}{h^2}(I\otimes T+T\otimes I).
$$

### Q36. `kron(I, T)` 和 `kron(T, I)` 分别代表什么？

`kron(I, T)` 表示在一个方向上对每条线应用一维 Laplacian，另一个方向保持不变。`kron(T, I)` 表示反过来，在另一个方向应用一维 Laplacian。

### Q37. 三维 Laplacian 的 Kronecker 形式是什么？

三维中：

$$
-\Delta=-\partial_{xx}-\partial_{yy}-\partial_{zz}.
$$

离散矩阵为：

$$
A=\frac{1}{h^2}
\left(
I\otimes I\otimes T
+
I\otimes T\otimes I
+
T\otimes I\otimes I
\right).
$$

这正对应 `solve_fd_heat_3D.py` 中的矩阵构造。

### Q38. `source.reshape(-1)` 做什么？

它把二维或三维数组按内存顺序 flatten 成一维向量。因为线性系统写成：

$$
AU=b,
$$

其中 $U$ 和 $b$ 都是一维向量，所以需要把网格上的 source array 展平。

### Q39. `matrix_nnz` 是什么？

`matrix_nnz` 是 sparse matrix 中非零元素的数量。它用于衡量矩阵稀疏程度和存储规模。

### Q40. residual 如何定义？

对线性系统：

$$
AU=b,
$$

残差为：

$$
r=b-AU.
$$

相对残差通常为：

$$
\frac{\|r\|}{\|b\|}.
$$

它用于检查数值解是否真的满足离散方程。

## E. 求解器与数值尺度

### Q41. CG 是什么？

CG 是 Conjugate Gradient method，共轭梯度法。它是求解 sparse、symmetric、positive definite 线性系统的经典迭代方法。

### Q42. 为什么 Poisson 问题适合 CG？

零 Dirichlet 条件下，Poisson 离散矩阵通常是 sparse、symmetric、positive definite。因此满足 CG 的基本适用条件。

### Q43. `rtol` 和 `atol` 分别是什么？

`rtol` 是 relative tolerance，控制相对残差停止条件。`atol` 是 absolute tolerance，控制绝对残差停止条件。SciPy CG 通常检查类似：

$$
\|b-Ax\|\leq \max(\mathrm{rtol}\|b\|,\mathrm{atol}).
$$

脚本中 `atol=0.0` 表示主要依赖相对残差。

### Q44. 直接求解器和迭代求解器如何选择？

直接求解器适合小到中等规模问题，优点是稳健、调参少；缺点是内存消耗大，会产生 fill-in。迭代求解器适合大规模 sparse systems，尤其是 3D 问题，但通常需要合适的预条件器。

### Q45. 为什么小规模 2D 问题中 CG 可能比 direct solver 慢？

小规模二维问题的矩阵不大，稀疏直接求解器开销很低。CG 虽然每步便宜，但需要多次迭代，并且 Python 调用、收敛判断也有开销。因此小规模测试中 direct solver 更快并不奇怪。

### Q46. 热源尺度如何估计温度尺度？

粗略地，若热源总量为：

$$
Q=\int_\Omega f\,dx,
$$

温度尺度会随 $Q$ 线性变化。更具体地：

$$
T=(-\Delta)^{-1}f.
$$

所以改变 source strength 或 source area 都会影响 $T_{\max}$。

### Q47. 为什么 Gaussian 和 square source 的 $T_{\max}$ 不能只看峰值强度？

因为峰值强度相同不代表总热源量相同。需要比较：

$$
\int_\Omega f\,dx.
$$

如果 Gaussian 的积分小于 square source，即使峰值都是 100，温度峰值也可能更低。

### Q48. 为什么 residual 的绝对值可能随自由度增加而变大？

绝对残差是向量范数，向量长度随自由度增加而增加。更多自由度意味着 residual vector 有更多分量，因此绝对残差不一定可直接跨网格比较。

### Q49. 为什么更应该看 relative residual？

relative residual 用 $\|b\|$ 做归一化：

$$
\frac{\|b-AU\|}{\|b\|}.
$$

它更适合比较不同网格、不同热源强度下的线性系统求解质量。

## F. FEniCSx 实现

### Q50. FEniCSx 中 mesh 对应什么数学对象？

`mesh` 对应计算区域 $\Omega$ 的离散划分。代码中：

```python
domain = mesh.create_unit_square(MPI.COMM_WORLD, n, n, mesh.CellType.triangle)
```

表示把单位正方形划分为 triangle cells。

### Q51. `cell_type=mesh.CellType.triangle` 是什么意思？

它表示二维 mesh 的单元是三角形。单位正方形中的小方格会被切分成 triangles，FEM basis functions 和积分都定义在这些 triangles 上。

### Q52. `v_space = fem.functionspace(domain, ("Lagrange", 1))` 是什么意思？

它定义 P1 Lagrange 有限元空间。该空间中的函数在每个三角形上是一次多项式，并且全局连续。

### Q53. `source = fem.Function(v_space)` 是什么意思？

它定义一个属于 $V_h$ 的有限元函数，用来表示离散热源 $f_h$。当前脚本把 source 也放入 P1 Lagrange 空间。

### Q54. `heat_source(x)` 中的 `x` shape 是什么？

在二维中，`x` 的形状通常是：

```python
(2, num_points)
```

其中 `x[0]` 是所有采样点的横坐标，`x[1]` 是纵坐标。返回值应为 shape 为 `(num_points,)` 的一维数组。

### Q55. `source.interpolate(heat_source)` 做什么？

它把 Python 函数 `heat_source` 插值到有限元空间中，得到 FEM 函数 $f_h$。对 P1 Lagrange 元来说，就是在节点上取 source 值，再在每个 triangle 内线性延拓。

### Q56. 为什么 square source 的 interpolation 需要 caveat？

真实 square source 是 discontinuous，但 P1 Lagrange 空间是连续、分片线性的。插值后参与积分的是 $f_h$，不是精确的分片常数 $f$。在 source 边界附近，$f_h$ 会线性过渡，粗网格下热源积分也会受 mesh 对齐影响。

### Q57. 如果想更忠实地表示分片常数热源，应怎么办？

可以使用 cellwise constant 的 `DG0` source，或在 UFL form 中使用 conditional expression。这样能更接近原始 discontinuous square source。

### Q58. FEniCSx 中 TrialFunction 和 TestFunction 分别是什么？

`TrialFunction` 表示未知函数所在方向，用来形成矩阵；`TestFunction` 表示测试函数，用来形成弱形式中的投影方程。在线性问题中：

```python
u = ufl.TrialFunction(v_space)
v = ufl.TestFunction(v_space)
```

### Q59. UFL 是什么？

UFL 是 Unified Form Language。它允许用接近数学表达式的方式写弱形式，例如：

```python
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
l_form = source * v * ufl.dx
```

### Q60. `a` 和 `l_form` 分别对应什么？

`a` 是 bilinear form：

$$
a(u,v)=\int_\Omega\nabla u\cdot\nabla v\,dx.
$$

`l_form` 是 linear form：

$$
L(v)=\int_\Omega fv\,dx.
$$

### Q61. Dirichlet boundary condition 在代码中如何施加？

代码先定位边界自由度：

```python
boundary_dofs = fem.locate_dofs_geometrical(v_space, boundary)
```

再施加：

```python
bc = fem.dirichletbc(fem.Constant(domain, 0.0), boundary_dofs, v_space)
```

这表示把边界节点温度固定为 0。

### Q62. PETSc options 中 `ksp_type="cg"` 是什么？

`ksp_type` 指 Krylov subspace solver 类型。`cg` 表示使用 Conjugate Gradient method。

### Q63. PETSc options 中 `pc_type="hypre"` 是什么？

`pc_type` 指 preconditioner 类型。`hypre` 通常用于代数多重网格预条件，适合 Poisson 类问题。

### Q64. `pc_type="lu"` 是什么？

`lu` 表示使用 LU factorization 作为求解方式或预条件方式。它更接近直接求解，通常稳健但内存消耗大，不适合很大规模 3D 问题。

### Q65. FEniCSx 输出的 `.xdmf` 和 `.h5` 分别是什么？

`.xdmf` 是 XML metadata，描述 mesh 和 field。`.h5` 存实际数值数据。ParaView 打开 `.xdmf` 时会自动读取配套 `.h5`。

### Q66. 为什么还要输出 `temperature.png`？

`temperature.png` 是 quick-look 图，用于在 VS Code、GitHub 或文档中快速检查结果。ParaView 仍然更适合交互式检查 mesh、field、contour 和切片。

### Q67. ParaView 中如何看 mesh？

把 representation 从 `Surface` 改成：

```text
Surface With Edges
```

即可同时显示温度场和 mesh edges。若只想看网格线，可用 `Wireframe`。

## G. 代码与结果分析

### Q68. `pcolormesh` 做什么？

`pcolormesh(x, y, value)` 根据网格坐标和对应数值画二维色彩图。它要求坐标数组和值数组在索引上对应，即同一个 index 表示同一个空间位置。

### Q69. `xx`、`yy`、`source` 的 shape 为什么要一致？

因为 `xx[i,j]`、`yy[i,j]`、`source[i,j]` 共同表示同一个网格点或 cell 附近的值。绘图时不是只看数值集合，而是看坐标和值之间的对应关系。

### Q70. `outdir.resolve()` 做什么？

`Path.resolve()` 返回路径的绝对形式。它常用于 summary 或 print 输出，方便以后知道结果文件具体写在哪里。

### Q71. `outdir.mkdir(parents=True, exist_ok=True)` 做什么？

它创建输出目录。`parents=True` 表示如果上级目录不存在，也一起创建；`exist_ok=True` 表示目录已存在时不报错。

### Q72. `from __future__ import annotations` 是什么？

它让类型注解延迟求值。这样可以减少运行时类型解析问题，也方便在类型注解中引用尚未定义的类或较新的 typing 语法。

### Q73. Gaussian source 为什么通常比 square source 产生更平滑的温度场？

Gaussian source 本身是光滑函数；square source 在边界处有跳变。椭圆型方程有平滑效应，但右端项的不连续仍会在温度梯度中留下局部变化。因此 Gaussian 的温度场通常更平滑。

### Q74. 为什么热点边缘附近需要更细 mesh？

因为热点边缘处 $f$ 变化快，甚至不连续。更细 mesh 可以更好解析局部梯度和曲率变化，降低几何采样误差和离散误差。

### Q75. 如何判断一次数值计算是否合理？

至少检查：

- 边界温度是否满足 Dirichlet 条件；
- $T_{\max}$ 是否出现在热源附近；
- `heat_source_integral` 是否符合预期；
- relative residual 是否足够小；
- 网格加密后关键 FoM 是否收敛；
- 结果图是否和物理直觉一致。

### Q76. 可以用哪些 FoM 比较结果？

常用 figures of merit 包括：

- $T_{\max}$；
- $T_{\min}$；
- mean temperature；
- heat source integral；
- relative residual；
- energy norm；
- 不同网格或方法之间的 field error。

### Q77. FD 和 FEM 结果为什么不会完全一样？

两者离散路径不同。FD 使用点值和 stencil；FEM 使用弱形式、basis functions 和数值积分。即使连续问题相同，有限网格下也会有不同离散误差。网格加密后，若热源和边界处理一致，两者应趋向同一个连续解。

### Q78. 为什么 FEniCSx 的 $h$ 和 FD 的 $h$ 定义不同？

FEniCSx 中 `n` 表示每个方向的 mesh subdivision 数，所以：

$$
h\approx\frac{1}{n}.
$$

FD 脚本中 `n` 表示内部未知点数，所以：

$$
h=\frac{1}{n+1}.
$$

比较 FD 和 FEniCSx 时要注意这个差异。

## H. 扩展问题

### Q79. 非均匀热导率问题的 PDE 是什么？

若热导率为 $k(x,y)$，稳态方程为：

$$
-\nabla\cdot(k\nabla T)=f.
$$

这比 $-\Delta T=f$ 更一般，适合多材料芯片模型。

### Q80. 非均匀热导率的弱形式是什么？

乘以测试函数并分部积分，得到：

$$
\int_\Omega k\nabla T\cdot\nabla v\,dx
=
\int_\Omega fv\,dx.
$$

因此 stiffness matrix 变为：

$$
A_{ij}=\int_\Omega k\nabla\phi_j\cdot\nabla\phi_i\,dx.
$$

### Q81. 瞬态热方程是什么？

常见形式为：

$$
\frac{\partial T}{\partial t}
-
\alpha\Delta T
=
f.
$$

其中 $\alpha$ 是热扩散率。

### Q82. 后向 Euler 如何离散瞬态热方程？

设时间步长为 $\Delta t$，则：

$$
\frac{T^{n+1}-T^n}{\Delta t}
-
\alpha\Delta T^{n+1}
=
f^{n+1}.
$$

有限元离散后通常得到：

$$
(M+\Delta t\,\alpha A)U^{n+1}
=
MU^n+\Delta t\,b^{n+1}.
$$

其中 $M$ 是 mass matrix，$A$ 是 stiffness matrix。

### Q83. 如何设计一个两种热导率区域的芯片模型？

可以令：

$$
k(x,y)=
\begin{cases}
k_1, & (x,y)\in\Omega_1,\\
k_2, & (x,y)\in\Omega_2.
\end{cases}
$$

其中 $\Omega_1$ 表示高导热材料区，$\Omega_2$ 表示低导热材料区。然后求解：

$$
-\nabla\cdot(k\nabla T)=f.
$$

### Q84. AI surrogate model 的输入和输出可以如何设计？

输入可以包括：

- 热源位置；
- 热源强度；
- 热源宽度；
- 材料热导率；
- 几何参数；
- 边界条件参数。

输出可以是：

- 全场温度 $T(x,y)$；
- $T_{\max}$；
- 热点位置；
- 若干 sensor points 的温度。

### Q85. surrogate model 为什么仍需要 PDE solver？

PDE solver 用来生成高质量训练数据，并作为验证基准。surrogate model 只能近似已学习参数范围内的输入输出映射，不能替代物理模型本身。

### Q86. 学完本项目后，下一步最自然的问题是什么？

最自然的扩展是：

- 用 FEniCSx 支持 Gaussian 和 two Gaussian source；
- 使用 `DG0` 表示 discontinuous square source；
- 做 FD/FEM field-level error；
- 加入非均匀热导率；
- 从 steady problem 走向 transient heat equation；
- 把 2D 模型扩展到 3D FEM。

---

[返回目录](README.md) | [上一章：调试、知识地图与学习路线](07_debugging_and_roadmap.md)
