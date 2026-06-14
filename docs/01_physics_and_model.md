# 物理背景与数学模型

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

---

[返回目录](README.md) | [上一章：总览与学习路线](00_overview.md) | [下一章：从强形式到弱形式](02_weak_form.md)
