# 肠菌靶点多肽结合剂数据库
## Gut Microbiota-Targeted Peptide Binder Database

本目录包含从 binder 初步设计、候选序列筛选、Rosetta 界面评估、GROMACS 分子动力学模拟、轨迹后处理、指标分析到结果绘图的脚本。脚本中的 `#SBATCH` 调度声明、个人账户和邮箱信息已移除。

## 1. 流程概览

```text
目标蛋白 input.pdb
        |
        v
[01] RFdiffusion 生成 binder 骨架
        |
        v
[01] ProteinMPNN 序列设计
        |
        v
[01] AlphaFold2 initial guess 预测
        |
        v
[02] 提取 AF2 候选 PDB，Rosetta 打分与 InterfaceAnalyzer 界面分析
        |
        v
选择候选 binder，并放入 GROMACS 输入目录
        |
        v
[03] GROMACS 预处理、溶剂化、加离子、能量最小化、NVT/NPT/正式 MD
        |
        v
[04] 轨迹去跳跃、居中、旋转平移拟合、导出结构帧
        |
        v
[05] RMSD、Rg、RMSF、最小距离、氢键、温度、势能、压力分析
        |
        v
[06] 将 .xvg 结果批量绘制成 PNG 图片
```

## 2. 文件说明

| 文件 | 作用 |
| --- | --- |
| `01_rfdiffusion_binder_design.sh` | 调用 RFdiffusion 生成 2000 个 binder 设计，随后进行 ProteinMPNN 序列设计和 AlphaFold2 initial guess 预测。 |
| `02_rosetta_binder_screening.sh` | 从 AF2 输出的 `out.silent` 中提取候选 PDB，使用 Rosetta `score_jd2` 和 `InterfaceAnalyzer` 评估结构与界面。 |
| `03_gromacs_md_simulation_gpu.sh` | 使用 GROMACS GPU 版本完成体系预处理、溶剂化、加离子、能量最小化、NVT/NPT 平衡和正式分子动力学模拟。 |
| `04_gromacs_trajectory_postprocess.sh` | 对 MD 轨迹执行 `nojump`、居中、旋转平移拟合，并导出 `md_fit.gro`。 |
| `05_gromacs_md_analysis.sh` | 自动解析 PDB 链范围，计算热力学、结构稳定性和 binder-蛋白结合指标。 |
| `06_plot_xvg_results.py` | 可选可视化脚本，批量读取 8 个 `.xvg` 文件并输出 PNG；依赖 Python 3、NumPy、Matplotlib。 |

## 3. 软件环境

| 软件环境 | 版本 | 官方地址 |
| --- | --- | --- |
| RFdiffusion | 源码快照，无 `.git` 元数据；功能对应 `v1.1.0` | https://github.com/RosettaCommons/RFdiffusion |
| dl_binder_design | `cafa3853ac94dceb1b908c8d9e6954d71749871a` | https://github.com/nrbennet/dl_binder_design |
| ProteinMPNN | `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57` | https://github.com/dauparas/ProteinMPNN |
| silent_tools | `1.0.0 / ced9d7051b4ecb0bc8e55f94ae2c76fd29a3372c` | https://github.com/bcov77/silent_tools |
| Rosetta | `2024.09+release.06b3cf8 / release-371` | https://rosettacommons.org/software/download/ |
| GROMACS | `2025.0` | https://ftp.gromacs.org/gromacs/gromacs-2025.0.tar.gz |
| AlphaFold2 | `alphafold_params_2022-12-06` | https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar |

## 4. 主要输出

### 01. RFdiffusion + ProteinMPNN + AF2

- `project_output/001_*.pdb`：RFdiffusion 设计结构；  # 001 为示例编号，具体以实际项目为准
- `project_output/test/design.silent`：ProteinMPNN 输入；
- `project_output/test/design2.silent`：ProteinMPNN 序列设计结果；
- AF2 initial guess 预测结果及 `out.silent`。

### 02. Rosetta 候选结构筛选

- `select.txt`：候选结构名称；
- `select/`：从 `out.silent` 提取出的 PDB；
- `last/`：Rosetta 打分和 InterfaceAnalyzer 结果。

### 03. GROMACS GPU 分子动力学模拟

主要输出位于 `MD/001/output/`：

- `md.tpr`、`md.xtc`、`md.edr`、`md.gro`、`md.cpt`、`md.log`；
- `em.gro`、`nvt.gro`、`npt.gro`；
- topology 文件；
- `atom_count.out`。

### 04. GROMACS 轨迹后处理

- `md_nojump.xtc`；
- `md_center.xtc`；
- `md_fit.xtc`；
- `md_fit.gro`；
- `index.ndx`；
- `postprocess.log`。

### 05. GROMACS 轨迹分析

分析结果输出到 `MD/001/output/xvg/`：

| 文件 | 内容 |
| --- | --- |
| `temp.xvg` | 温度 |
| `potential.xvg` | 势能 |
| `pressure.xvg` | 压力 |
| `rmsd.xvg` | RMSD |
| `rg.xvg` | 回转半径 Rg |
| `rmsf.xvg` | 残基 RMSF |
| `mindist.xvg` | binder-蛋白最小距离 |
| `hbond.xvg` | binder-蛋白氢键数目 |

### 06. xvg 结果可视化

图片默认输出到 `MD/001/output/xvg/plots/`：

- `temperature.png`、`potential.png`、`pressure.png`；
- `rmsd.png`、`rg.png`、`rmsf.png`；
- `mindist.png`、`hbond.png`。

## 5. 已保留的关键复现参数

- RFdiffusion 热点残基：`A182,A155,A169,A175,A189,A204`；  # 备注：此为实例热点，具体以实际为准
- RFdiffusion 设计数量：`inference.num_designs=2000`；
- RFdiffusion 噪声参数：`noise_scale_ca=0`、`noise_scale_frame=0`；
- ProteinMPNN 权重：`v_48_020.pt`；
- GROMACS 力场：`charmm27`；
- 水模型：`tip3p`；
- 模拟盒子：十二面体，`-d 3.0`；
- 离子：`-neutral`，阳离子 `NA`，阴离子 `CL`；
- 最小距离接触截断：`0.4 nm`；
- 轨迹结构导出时间：`50000 ps`。

## 6. 脱敏与修改说明

1. 删除所有脚本开头的 `#SBATCH` 调度声明；
2. 删除个人账户名、邮箱等调度平台信息；
3. 将本机绝对路径替换为集中配置的占位路径；
4. 将步骤 02 对 `SLURM_SUBMIT_DIR` 的依赖改为自动定位脚本所在目录；
5. 在步骤 03 的结果归档部分补充移动 `md.tpr`，以保证步骤 04 和步骤 05 能在 `output/` 目录中找到该文件；
6. 未修改科学计算命令、关键参数、热点残基、设计数量和分析指标。

