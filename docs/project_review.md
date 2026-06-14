# 二维热扩散项目检阅

## 1. 当前目标

本项目围绕稳态热扩散/Poisson 问题：

$$
-\Delta T=f,\qquad T|_{\partial\Omega}=0
$$

建立从数学模型、有限差分、有限元弱形式到 FEniCSx 实现的完整学习链条。

当前项目已经进入 FEM/FEniCSx 阶段。有限差分部分主要作为可控参照系，用于验证尺度、边界条件、残差和收敛行为。

## 2. 环境状态

当前采用 conda workflow，不再维护 `.venv`。

- Anaconda: `/opt/anaconda3`
- conda: `26.1.1`
- FD 脚本运行环境: conda `base`
- FEniCSx 环境: `/opt/anaconda3/envs/fenicsx`
- FEniCSx wrapper: `scripts/run_fenicsx_heat.sh`

旧 Anaconda `/Users/glh/anaconda3` 已删除。项目中的 `.venv` 也已删除。

FEniCSx smoke test 已通过：

```text
n = 16
dofs = 289
source_integral = 3.515625
t_min = 0
t_max = 1.158638793628342
```

## 3. 代码结构

### 有限差分

- `scripts/solve_fd_heat.py`
  - 2D steady heat diffusion finite-difference solver。
  - 支持 `square`、`gaussian`、`two_gaussian`。
  - 输出 `solution.npz`、`summary.csv`、`temperature.png`。
  - 记录网格尺度、热源面积/积分、残差、矩阵非零元数、最大温度。

- `scripts/run_fd_experiments.py`
  - 批量网格和热源强度实验。
  - 输出 mesh/strength/residual studies。

- `scripts/solve_fd_heat_3D.py`
  - 3D finite-difference solver。
  - 支持 3D Gaussian 测试和 3D 可视化。

### FEniCSx

- `scripts/solve_fenicsx_heat.py`
  - 2D FEniCSx/FEM solver。
  - 当前实现 square source。
  - 使用 P1 Lagrange 元，弱形式为：

$$
\int_\Omega \nabla T_h\cdot\nabla v_h\,dx
=
\int_\Omega f v_h\,dx
$$

- `scripts/run_fenicsx_heat.sh`
  - FEniCSx wrapper。
  - 清理 `CC/CXX/CFLAGS/LDFLAGS`，默认设置 JIT cache 到 `/private/tmp`。

## 4. 当前核心结果

### FD baseline

| case | n | dofs | source integral | Tmax | relative residual |
|---|---:|---:|---:|---:|---:|
| square | 64 | 4096 | 4.000000 | 1.309152 | 4.99e-14 |
| gaussian | 64 | 4096 | 2.261947 | 0.768182 | 5.45e-14 |
| two_gaussian | 64 | 4096 | 4.523893 | 0.936249 | 6.20e-14 |
| 3D gaussian | 64 | 262144 | 0.340192 | 0.310906 | 9.35e-11 |

### FD 与 FEniCSx 对比

对齐网格尺度后，即 FD 的 `n` 对应 FEniCSx 的 `n+1`：

| FD `n` | FEniCSx `n` | h | FD Tmax | FEniCSx Tmax | relative difference |
|---:|---:|---:|---:|---:|---:|
| 16 | 17 | 0.05882353 | 1.65596296 | 1.60044130 | -3.35% |
| 32 | 33 | 0.03030303 | 1.12783626 | 1.11279177 | -1.33% |
| 64 | 65 | 0.01538462 | 1.30915239 | 1.30521504 | -0.30% |

结论：当 $h$ 和热源积分对齐后，两种方法结果接近；网格加密后差异明显下降。

## 5. 原理解读

有限差分直接离散强形式：

$$
-\Delta T=f
$$

在规则网格点上构造五点 stencil。

FEM/FEniCSx 先做弱形式，再在有限元空间中求解 residual 的正交条件：

$$
\int_\Omega \left(-\Delta T_h-f\right)v_h\,dx=0
$$

其中 $v_h$ 取自测试函数空间。

因此二者矩阵来源不同：

- FD 矩阵来自差分 stencil；
- FEM 刚度矩阵来自 $\int_\Omega\nabla\phi_i\cdot\nabla\phi_j\,dx$；
- FD 的未知量是内部格点温度；
- FEM 的自由度是有限元节点系数，Dirichlet 边界节点被固定为零。

两者离散路径不同，但在同一连续问题、同一边界条件、同一热源极限下应收敛到同一个解。

## 6. 已知限制

1. `solve_fenicsx_heat.py` 当前只支持 square source；还没有支持 gaussian/two_gaussian。
2. FEniCSx 当前只输出 XDMF/HDF5 和 summary，没有生成 PNG 图。
3. FD 与 FEM 对比目前主要比较 `Tmax` 和 source integral；后续应加入场级误差比较。
4. FEniCSx 中的 square source 是通过 nodal interpolation 表示的，粗网格下热源积分会随 mesh 对齐方式变化。
5. 目前 FEniCSx 与 FD 的结果组织已经统一，但还没有自动化 comparison runner。

## 7. 下一阶段建议

优先级从高到低：

1. 扩展 FEniCSx source：支持 `square`、`gaussian`、`two_gaussian`，与 FD 脚本一致。
2. 增加 FEniCSx 可视化：输出温度图和热源图。
3. 写一个统一 comparison 脚本，自动运行 FD/FEniCSx 对齐网格并生成表格。
4. 对 FEM 结果做残差/能量范数分析：

$$
\|T_h\|_E^2=\int_\Omega |\nabla T_h|^2 dx
$$

5. 用同一热源积分归一化后比较不同 source 的温度尺度。
6. 保持 `results/` 只存放 canonical outputs；临时 smoke test 输出放到 `/tmp` 或 `results/tmp_*` 并及时删除。

## 8. 建议保留的核心结果

建议作为后续分析基准保留：

- `results/baselines/fd_single`
- `results/baselines/fd_gaussian`
- `results/baselines/fd_two_gaussian`
- `results/baselines/fd_3d_gaussian`
- `results/experiments/fd_experiment_single`
- `results/experiments/fd_experiment_gaussian`
- `results/comparisons/fd_fenicsx_comparison`
- `results/comparisons/fenicsx_compare_n17`
- `results/comparisons/fenicsx_compare_n33`
- `results/comparisons/fenicsx_compare_n65`

环境验证类结果已经从 canonical results tree 中删除。后续如果需要重新做 smoke test，可临时输出到 `results/tmp_*` 或 `/tmp`，确认通过后删除。
