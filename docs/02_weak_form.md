# 从强形式到弱形式

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

---

[返回目录](README.md) | [上一章：物理背景与数学模型](01_physics_and_model.md) | [下一章：有限元离散与 FEniCSx 基础](03_fem_discretization.md)
