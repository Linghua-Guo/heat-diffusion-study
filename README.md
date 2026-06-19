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

`requirements.txt` 只覆盖有限差分脚本所需的轻量 Python 依赖：

```bash
conda activate base
pip install -r requirements.txt
```

FEniCSx 依赖 PETSc、MPI 和 C++/JIT 编译工具链，不建议放进普通 `pip`/`venv`。本项目默认使用独立 conda 环境：

```bash
conda create -n fenicsx -c conda-forge fenics-dolfinx matplotlib
```

创建后可以检查：

```bash
conda activate fenicsx
python -c "import dolfinx; print(dolfinx.__version__)"
python -c "import mpi4py, petsc4py, matplotlib"
```

当前 wrapper 默认使用：

```text
/opt/anaconda3/envs/fenicsx/bin/python
```

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
conda activate fenicsx
bash scripts/run_fenicsx_heat.sh --n 64 --outdir results/baselines/fenicsx_single
```

该脚本会输出 ParaView 可读的 `temperature.xdmf`/`temperature.h5`，并额外保存 matplotlib quick-look 图 `temperature.png`。

若 FEniCSx 环境路径不同：

```bash
FENICSX_PYTHON=/path/to/fenicsx/bin/python bash scripts/run_fenicsx_heat.sh --n 64
```

## 测试

本地运行轻量测试：

```bash
conda activate base
pytest -q
```

GitHub Actions 会在 push 和 pull request 时运行这些 FD smoke tests。默认 CI 不安装 FEniCSx；FEniCSx 仍建议在本地 `fenicsx` conda 环境中验证。

## 目录

```text
docs/       学习讲义与项目文档
scripts/    可执行脚本
results/    已整理的结果
```
