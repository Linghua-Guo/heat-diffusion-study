# 二维热扩散问题学习项目

本项目围绕稳态热扩散/Poisson 方程：

$$
-\Delta T=f,\qquad T|_{\partial\Omega}=0
$$

目标是从物理建模、有限差分、有限元弱形式一路走到 FEniCSx 实现与结果分析。

## 学习文档

按顺序阅读：

- [学习文档总览](docs/README.md)
- [项目检阅](docs/project_review.md)
- [结果索引](results/INDEX.md)

## 环境

本项目使用 conda。

- FD 脚本：conda `base`
- FEniCSx 脚本：conda `fenicsx`

检查：

```bash
conda --version
conda env list
```

## 运行 FD

```bash
conda activate base
python scripts/solve_fd_heat.py --n 64 --source square --outdir results/baselines/fd_single
```

批量实验：

```bash
python scripts/run_fd_experiments.py --mesh-n 16 32 64 128 --strengths 50 100 200 500 --outdir results/experiments/fd_experiments
```

## 运行 FEniCSx

```bash
bash scripts/run_fenicsx_heat.sh --n 64 --outdir results/baselines/fenicsx_single
```

若 FEniCSx 环境路径不同：

```bash
FENICSX_PYTHON=/path/to/fenicsx/bin/python bash scripts/run_fenicsx_heat.sh --n 64
```

## 目录

```text
docs/       学习讲义与项目文档
scripts/    可执行脚本
results/    已整理的结果
```

