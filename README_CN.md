# meloculardftcal：分子 DFT / TD-DFT / SOC 高通量工作流

[English README](README_EN.md) | [中文智能体接口](AGENT_INTERFACE_CN.md) | [English agent interface](AGENT_INTERFACE_EN.md) | [第三方说明](THIRD_PARTY_NOTICES.md)

仓库地址：<https://github.com/letsgo-2/meloculardftcal>

本项目是一套以 ORCA 为核心的批量分子量化计算与后处理工作流，覆盖三维结构准备、基态几何优化与频率计算、TD-DFT 激发态与自旋轨道耦合（SOC）计算、Slurm 调度、结果解析、性质合并、分子描述符计算，以及统计和结构可视化。

仓库保存代码、参数模板、说明文档、小型示例，以及唯一一套
`pic16/S1` 核心表格数据。ORCA 安装包、OpenMPI 源码、原始量化日志、GBW
波函数、其他完整批次和批量图片不应提交到 Git。

## 1. 工作流与科学目标

```text
分子生成 CSV / SMILES
        │
        ▼
Open Babel 生成三维结构
        │
        ▼
ORCA 基态 Opt + Freq
        │
        ├── 优化结构 .xyz
        ├── 分子轨道 .gbw
        └── 基态/热力学/频率性质
        │
        ▼
TD-DFT + Triplet + SOC
        │
        ├── S1/T1 能量
        ├── S1 振子强度与波长
        ├── S1−Ti 能隙
        └── SOCME 与 SOCME/|ΔE|
        │
        ▼
原始生成指标 + DFT + RDKit 描述符
        │
        ▼
分布比较、Pareto 筛选、化学空间和结构可视化
```

当前主要优化方向是寻找同时具有以下性质的候选分子：

- 较小的正 `S1−Ti` 能隙；
- 较大的 `SOCME(S1,Ti)`；
- 较大的 `SOCME/|ΔE|`；
- 适合目标应用的 S1 能量、吸收波长和振子强度；
- 无虚频的稳定优化结构；
- 可接受的分子量、LogP、TPSA、可旋转键和 SA Score；
- 在候选库中保持足够的结构多样性。

这些是多目标优化方向，不是固定的普适阈值。项目中的 LogP/TPSA 窗口、电荷阈值和 Pareto 轴范围均可调整。

## 2. 仓库结构

```text
.
├── README.md                         # 语言入口
├── README_CN.md                      # 完整中文说明
├── README_EN.md                      # Complete English documentation
├── AGENT_INTERFACE.md                # 智能体接口语言入口
├── AGENT_INTERFACE_CN.md             # 完整中文智能体接口
├── AGENT_INTERFACE_EN.md             # Complete English agent interface
├── workflow.py                       # 机器可读接口查询 CLI
├── THIRD_PARTY_NOTICES.md            # SA Score 等第三方组件说明
├── environment.yml                   # 推荐 Conda 环境
├── requirements.txt                  # pip 依赖
├── config/
│   ├── agent_interface.json          # 阶段、路径、版本和默认参数注册表
│   ├── cluster.env.example           # 集群资源和 ORCA 路径模板
│   └── pictos.txt                    # pic 批次到展示标签映射
├── 01_structure_preparation/
│   ├── convertInputs.py              # SMILES → 3D GJF
│   ├── gjf_to_inp.sh                 # 普通体系 Opt/Freq 输入
│   ├── gjf_to_inpp13.sh              # 重元素/ZORA Opt/Freq 输入
│   ├── gjf_to_inptd.sh               # 普通体系 TD-DFT/SOC 输入
│   └── gjf_to_inptdp13.sh            # 重元素/ZORA TD-DFT/SOC 输入
├── 02_cluster_execution/
│   ├── suboptsbatch.sh                # 批量提交 Opt/Freq
│   ├── subtdsbatch.sh                 # 批量提交 TD/SOC
│   ├── find_missing_opt_ids.sh        # 基态缺失检测
│   ├── find_missing_td_ids.sh         # TD 缺失检测
│   ├── list_opt_jobs.sh               # 查询基态队列任务
│   ├── list_td_jobs.sh                # 查询 TD 队列任务
│   ├── auto_submit_opt.sh             # 基态自动补交
│   └── auto_submit_td.sh              # TD 自动补交
├── 03_result_parsing/
│   ├── parse_orca_opt.py              # 基态日志解析
│   ├── parse_orca_td_socv2.py         # TD/SOC 主解析器
│   └── merge_dft_csv.py               # TD 与基态结果合并
├── 04_data_processing/
│   ├── merge_summary.py               # 与原始生成表合并
│   ├── addrdkit.py                    # RDKit 描述符和 SA Score
│   ├── addcenter.py                   # attachment-point 对称性
│   └── calculate_molecular_weight.py  # 分子量和分布
├── 05_analysis/
│   ├── plotters/                      # Python 分析和绘图脚本
│   └── data/                          # 标签和小型参考点
├── 06_visualization/
│   ├── generate_molecule_images.py
│   ├── generate_molecule_images_with_props.py
│   ├── generate_top_selected_images.py
│   ├── execlwithimage.py
│   └── scaffold/                      # 骨架、装饰基团和对称性图
├── 07_validation/
│   ├── check_environment.py           # 软件/硬件能力报告
│   ├── collect_metrics.py             # 数据质量和科学性能指标
│   └── UnitTestSAScore.py
├── data/
│   ├── core/pic16_s1/                 # 随仓库提供的 222 条 S1 核心数据
│   ├── raw/                           # 原始 CSV；默认不提交
│   ├── processed/                     # v20/v21/v22；默认不提交
│   └── example_molecules.csv
├── examples/example_v22.csv
├── reports/                           # 环境、指标和运行清单；可提交小型 JSON
├── work/                              # ORCA 输入和运行产物；默认不提交
└── results/figures/                   # 图片输出；默认不提交
```

所有命令都应从仓库根目录运行。文档中的路径均为相对仓库根目录的相对路径。

### 2.1 机器可读接口

中文智能体应先读取 [`AGENT_INTERFACE_CN.md`](AGENT_INTERFACE_CN.md)，
英文智能体应读取 [`AGENT_INTERFACE_EN.md`](AGENT_INTERFACE_EN.md)，再通过
`config/agent_interface.json` 获取稳定的阶段、入口、版本、路径模板和默认
参数。只读查询命令：

```bash
python3 workflow.py stages
python3 workflow.py stage S41_PARSE_TD
python3 workflow.py software
python3 workflow.py parameters
python3 workflow.py paths --check
```

输出均为 JSON。`config/agent_interface.json` 是接口注册表；真正写入 ORCA
输入的参数仍以四个 `01_structure_preparation/gjf_to_*.sh` 为执行权威，
修改后必须重新生成输入。

## 3. 参考运行环境

### 3.1 历史集群环境

当前脚本最初在以下环境中使用。该表描述已记录的参考环境，不表示新集群必须完全相同。

| 项目 | 参考配置 |
|---|---|
| 操作系统 | CentOS Linux 7，Linux 3.10 |
| 调度器 | Slurm |
| ORCA | 6.0.1，shared OpenMPI 4.1.6 AVX2 发行包 |
| MPI | OpenMPI 4.1.6 |
| 通信库 | UCX 1.15.0 |
| CPU | 2 × AMD EPYC 7H12 64-Core |
| 物理核心 | 128 核/节点 |
| 逻辑 CPU | 256 线程/节点，SMT2 |
| 内存 | 约 503 GiB/节点 |
| Swap | 无 |
| GPU | 工作流不依赖 GPU |
| 文件系统 | 集群共享并行文件系统 |

集群系统默认 Python 曾为 2.7，但本仓库的数据处理和绘图统一要求 Python 3.9 或更高版本。建议在登录节点或分析节点创建独立 Conda 环境，不要依赖系统 Python。

### 3.2 历史计算规模

整理前工作区包含的输入规模约为：

| 批次 | 分子数 |
|---|---:|
| pic3 | 154 |
| pic4 | 300 |
| pic8 | 601 |
| pic9 | 77 |
| pic10 | 150 |
| pic13 | 537 |
| pic16 | 222 |
| pic24 | 148 |
| pic25 | 400 |
| pic26 | 400 |
| publish 参考集 | 9 |
| 合计 | 2,998 |

这是输入文件规模，不等同于成功完成的计算数。每个候选通常对应基态 `.inp/.sh` 和 TD `.inp/.sh`，运行后还会产生日志、XYZ、GBW 和 Slurm 输出。

### 3.3 单任务资源参数

资源由 [`config/cluster.env.example`](config/cluster.env.example) 提供：

| 参数 | Opt/Freq 默认 | TD/SOC 默认 | 含义 |
|---|---:|---:|---|
| `SBATCH_PARTITION` | `amd_512q` | `amd_512q` | Slurm 分区 |
| `SBATCH_NODES` | 1 | 1 | 单任务节点数 |
| `SBATCH_NTASKS` | 64 | 32 | MPI rank 数，同时写入 `%pal nprocs` |
| `OMP_NUM_THREADS` | 1 | 1 | 每个 MPI rank 的 OpenMP 线程 |
| `%maxcore` | 3000 MB | 3000 MB | ORCA 每个处理进程的近似内存上限 |

按 `maxcore × nprocs` 粗略估算：

- 64 ranks：约 192,000 MB，即约 187.5 GiB；
- 32 ranks：约 96,000 MB，即约 93.75 GiB。

实际峰值可能高于简单乘积。调整 `SBATCH_NTASKS` 或 `%maxcore` 时必须同时考虑节点物理内存、调度器内存策略、ORCA 模块的并行效率和分子大小。

自动补交脚本默认：

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `MAX_JOBS` | 20 | 当前用户整个 Slurm 队列的总任务上限，不仅是当前批次 |
| `SLEEP_SEC` | 1800 s | 每 30 分钟检查一次 |

为避免智能体误取消任务或删除提交脚本，历史超时清理脚本没有收入公开整理版。
需要取消任务时应显式检查 `squeue` 后调用 `scancel <job_id>`，并单独获得外部
调度器修改权限。

## 4. 软件安装

### 4.1 Conda 环境

```bash
conda env create -f environment.yml
conda activate dft-workflow
```

环境包含：

- NumPy、pandas、SciPy；
- Matplotlib、Seaborn、Pillow；
- RDKit、Open Babel；
- scikit-learn、UMAP；
- openpyxl、tqdm。

### 4.2 pip 环境

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

结构生成和部分绘图使用 `obabel` 命令行。pip 依赖不会安装该可执行文件，
因此推荐从 conda-forge 安装 Open Babel。

### 4.3 集群 ORCA 配置

```bash
cp config/cluster.env.example config/cluster.env
```

编辑 `config/cluster.env`：

```bash
export SBATCH_PARTITION="amd_512q"
export SBATCH_NODES="1"
export SBATCH_NTASKS="32"
export WORK_DIR="<集群上的仓库绝对路径>"
export ORCA_BIN="<集群上的ORCA可执行文件绝对路径>"
export MODULE_INIT="<集群上的modules初始化脚本绝对路径>"
export ORCA_MODULE="orca/6_0_1_linux_x86-64_shared_openmpi416_avx2"
```

然后加载：

```bash
source config/cluster.env
```

`config/cluster.env` 含本机路径，已被 `.gitignore` 排除。
模板中的 `SBATCH_NTASKS=32` 是共享初值；按照本仓库的参考资源配置，生成
基态任务前显式设为 64，生成 TD 任务前再设为 32。修改该值后必须重新生成
`.inp/.sh`，确保 Slurm rank 数与 ORCA `nprocs` 一致。

## 5. 量化计算协议

### 5.1 普通元素体系：基态 Opt/Freq

入口：

```text
01_structure_preparation/gjf_to_inp.sh
```

默认 ORCA 关键字：

```text
! B3LYP def2-SVP Opt Freq RIJCOSX CPCM(Water)
%maxcore 3000
%pal nprocs <SBATCH_NTASKS> end
```

含义：

- 泛函：B3LYP；
- 基组：def2-SVP；
- 任务：几何优化和谐振频率；
- RIJCOSX：加速混合泛函计算；
- 溶剂：CPCM 水；
- 输出用于提取稳定结构、HOMO/LUMO、偶极矩、Gibbs 自由能和虚频。

### 5.2 重元素/ZORA 体系：基态 Opt/Freq

入口：

```text
01_structure_preparation/gjf_to_inpp13.sh
```

默认关键字：

```text
! PBE0 D3BJ ZORA ZORA-def2-TZVP SARC/J Opt Freq RIJCOSX TightSCF CPCM(Water)
```

与普通体系相比：

- 泛函改为 PBE0；
- 使用 D3BJ 色散修正；
- 使用 ZORA 标量相对论处理；
- 使用 ZORA-def2-TZVP 和 SARC/J；
- 使用 TightSCF。

该模板最初为特定重元素批次准备。是否适用于新的元素和电荷状态必须由计算化学方案决定。

### 5.3 普通体系：TD-DFT/SOC

入口：

```text
01_structure_preparation/gjf_to_inptd.sh
```

默认关键字：

```text
! wB97X-D3 def2-TZVP RIJCOSX TightSCF CPCM(Water)
%tddft
    NRoots   10
    Triplets true
    dosoc true
end
```

此外：

```text
%scf
  Guess MORead
end
%moinp "<base>.gbw"
* xyzfile <charge> <multiplicity> <base>.xyz
```

即 TD 阶段读取基态优化后的几何和 GBW 波函数。必须先确认基态任务正常结束，并且 `.xyz` 与 `.gbw` 存在。

### 5.4 重元素/ZORA：TD-DFT/SOC

入口：

```text
01_structure_preparation/gjf_to_inptdp13.sh
```

默认关键字：

```text
! wB97X-D3 ZORA ZORA-def2-TZVP SARC/J RIJCOSX RI-SOMF(1X) TightSCF CPCM(Water)
```

TD 块仍为 10 个根、包含三重态并计算 SOC；另外使用 ZORA 和 `RI-SOMF(1X)`。

### 5.5 可调计算参数

| 参数 | 修改位置 | 修改后动作 |
|---|---|---|
| 泛函、基组、色散、溶剂、SCF 关键字 | 四个 `01_structure_preparation/gjf_to_*.sh` 中的 `echo "! ..."` | 删除或归档旧生成输入，重新生成 `.inp/.sh` |
| TD 根数 | `01_structure_preparation/gjf_to_inptd.sh`、`01_structure_preparation/gjf_to_inptdp13.sh` 中 `NRoots` | 重新生成 TD 输入并重新计算 |
| 是否计算 Triplet/SOC | 同上，`Triplets`、`dosoc` | 重新生成并重新计算 |
| 每进程内存 | 四个生成脚本中的 `%maxcore` | 与节点内存和 rank 数共同复核 |
| MPI ranks | `config/cluster.env` 的 `SBATCH_NTASKS` | 重新 `source` 配置并重新生成提交脚本 |
| 分区和节点数 | `SBATCH_PARTITION`、`SBATCH_NODES` | 重新生成提交脚本 |
| 总电荷与多重度 | 原始 GJF 的 charge/multiplicity 行 | 重新生成全部 ORCA 输入 |
| 水溶剂模型 | `CPCM(Water)` | 改动后必须作为新的计算协议保存，不能与旧结果混合 |

建议每次协议变化都保存一个运行清单，至少记录 ORCA 版本、MPI 版本、完整关键字、资源参数、日期和 Git commit。

## 6. 数据和命名契约

原始 CSV 至少需要：

```csv
MOLNAME,SMILES
pic16_0001,CCO
pic16_0002,c1ccccc1
```

推荐任务命名：

```text
pic16_0001.inp
pic16_0001.sh
pic16_0001.log
pic16_0001.xyz
pic16_0001.gbw
pic16_0001td.inp
pic16_0001td.sh
pic16_0001td.log
```

关键约束：

1. `<batch>_<四位编号>` 必须在所有阶段保持一致；
2. TD 解析器通过日志编号寻找原始 CSV 的对应数据行；
3. `0001` 对应 CSV 第一条数据行，`0012` 对应第十二条数据行；
4. 如果 CSV 被重新排序，SMILES 与量化结果会错配；
5. `03_result_parsing/merge_dft_csv.py` 按 `filename` 合并；
6. `04_data_processing/merge_summary.py` 将 `filename` 转成 `ID` 后与原始生成表合并。

原始生成表还可以包含：

- `Scaffold`
- `st_gap`
- `charge_score`
- `symmetry_score`
- `decoration_ring_penalty`
- `total_score`
- `Step`
- `Molecular_Weight`

### 6.1 随仓库提供的 pic16 / S1 核心数据

`config/pictos.txt` 将 `pic16` 映射为展示标签 `S1`。公开整理版只保留这一
套完整的 222 条核心数据：

| 路径 | 内容 |
|---|---|
| `data/core/pic16_s1/input/pic16.csv` | 生成指标、Scaffold 和 SMILES |
| `data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv` | TD-DFT/SOC 汇总 |
| `data/core/pic16_s1/intermediate/pic16inpdftv20.csv` | 基态与 TD/SOC 合并 |
| `data/core/pic16_s1/intermediate/pic16inpdftv21.csv` | 加回原始生成指标 |
| `data/core/pic16_s1/final/pic16inpdftv22.csv` | 最终 RDKit/SA Score 数据 |
| `05_analysis/data/publishinpdftv22.csv` | 用于相对性能比较的小型参考分子表 |
| `reports/pic16_s1/metrics.json` | 数据质量和科学性能快照 |

文件来源、行数和校验值见
[`data/core/pic16_s1/DATASET.md`](data/core/pic16_s1/DATASET.md)。原始 ORCA
日志、GBW、XYZ、Excel 和批量图片不放入仓库。

## 7. 完整运行流程

以下命令均从仓库根目录执行。

### 7.1 SMILES → 三维 GJF

```bash
python3 01_structure_preparation/convertInputs.py \
  1 data/raw/pic16.csv work/pic16
```

三个位置参数依次为模式、源 CSV 和可选输出目录。省略输出目录时，为兼容历史
用法，脚本写入仓库根目录下与 CSV 同名的目录；推荐始终显式使用
`work/<batch>`。

`01_structure_preparation/convertInputs.py` 的模式：

| 模式 | 功能 |
|---:|---|
| 1 | SMILES → 3D GJF，结构优化准备 |
| 2 | 历史 Gaussian LOG → TD GJF |
| 3 | 历史 Gaussian Guess=Read TD 输入 |

当前 ORCA 主链路主要使用模式 1。也可以直接调用：

```bash
mkdir -p work/pic16
obabel -:"CCO" -O work/pic16/pic16_0001.gjf --gen3d
```

### 7.2 生成和提交基态任务

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh work/pic16 work/pic16inp
bash 02_cluster_execution/suboptsbatch.sh work/pic16inp 1-20
```

支持：

```bash
bash 02_cluster_execution/suboptsbatch.sh work/pic16inp 1-10 20 30-35
bash 02_cluster_execution/suboptsbatch.sh work/pic16inp "[1-10],[20-30],78"
```

检查缺失：

```bash
bash 02_cluster_execution/find_missing_opt_ids.sh work/pic16inp
```

自动补交：

```bash
bash 02_cluster_execution/auto_submit_opt.sh work/pic16inp 20 1800
```

脚本名保留历史 `pic16`，实际会从目录名推断 `pic3`、`pic13` 等前缀。

### 7.3 基态质量门控

在进入 TD 前至少检查：

- 日志包含 `ORCA TERMINATED NORMALLY`；
- 优化完成；
- `.xyz` 和 `.gbw` 存在且非空；
- `Imaginary == No`，或对预期过渡态另行处理；
- 最终结构没有解离、异常键长或错误质子化；
- 电荷和多重度正确。

### 7.4 生成和提交 TD/SOC

```bash
export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh work/pic16 work/pic16inp
bash 02_cluster_execution/subtdsbatch.sh work/pic16inp 1-20
```

TD 缺失检测：

```bash
bash 02_cluster_execution/find_missing_td_ids.sh work/pic16inp
```

自动补交：

```bash
bash 02_cluster_execution/auto_submit_td.sh work/pic16inp 20 1800
```

TD 待运行集合定义为：

```text
基态正常结束编号 − TD 正常结束编号
```

### 7.5 解析基态

```bash
python3 03_result_parsing/parse_orca_opt.py \
  work/pic16inp \
  -o work/pic16inpout
```

提取：

- 最终单点能、核排斥能；
- HOMO、LUMO 和 gap；
- 偶极矩；
- Gibbs 自由能；
- 虚频标记。

注意：当前解析器以能否提取 HOMO 作为是否写入汇总的重要条件，并不会独立强制验证 `ORCA TERMINATED NORMALLY`。自动工作流应在解析前先做正常结束检查。

### 7.6 解析 TD/SOC

```bash
python3 03_result_parsing/parse_orca_td_socv2.py \
  work/pic16inp \
  -o work/pic16inpout \
  --csv data/raw/pic16.csv
```

批量解析会直接生成
`work/pic16inpout/pic16inpout_tdv2.csv`，名称与合并阶段契约一致。

解析定义：

- `S1-T1 = E(S1) − E(T1)`，单位 eV；
- `min_positive_S1-Ti` 是所有正 `E(S1) − E(Ti)` 中的最小值；
- SOCME 将 ORCA 输出的 X/Y/Z 三个复分量做模长合成，单位 cm⁻¹；
- `SOCME/|ΔE|` 使用 SOCME 除以绝对能隙，单位 cm⁻¹/eV；
- 汇总最大 SOCME 和最大 `SOCME/|ΔE|` 时，只考虑不高于最接近 S1 的下方三重态范围。

### 7.7 合并为 v20/v21/v22

```bash
python3 03_result_parsing/merge_dft_csv.py \
  work/pic16inp \
  -o work/pic16inpout
```

得到：

```text
work/pic16inpout/pic16inpdftv20.csv
```

合并原始生成指标：

```bash
python3 04_data_processing/merge_summary.py \
  --piccsvfile data/raw/pic16.csv \
  --summary work/pic16inpout/pic16inpdftv20.csv \
  --output data/processed/pic16inpdftv21.csv
```

增加 RDKit 指标：

```bash
python3 04_data_processing/addrdkit.py \
  data/processed/pic16inpdftv21.csv \
  --output data/processed/pic16inpdftv22.csv
```

版本含义：

| 版本 | 内容 |
|---|---|
| v20 | 基态 + TD/SOC |
| v21 | v20 + 原始生成/骨架指标 |
| v22 | v21 + LogP、TPSA、可旋转键、形式电荷、SA Score |
| v23 | v22 的带结构图片 Excel 展示版 |

## 8. 结果字段与优化方向

| 字段 | 单位 | 方向/用途 |
|---|---|---|
| `S1_energy_eV` | eV | 目标激发能；根据应用设窗口 |
| `T1_energy_eV` | eV | 三重态能量 |
| `S1_wavelength_nm` | nm | 吸收波长；根据光源/治疗窗口设目标 |
| `S1_fosc_D2` | 无量纲 | S1 振子强度，通常越大表示允许跃迁更强 |
| `S1-T1` | eV | S1 与 T1 能隙 |
| `min_positive_S1-Ti` | eV | 越小通常越有利于接近共振 |
| `max_SOCME_S1_Ti_value` | cm⁻¹ | 越大表示自旋轨道耦合更强 |
| `max_SOCME_over_deltaE_S1_Ti_value` | cm⁻¹/eV | 耦合与能隙联合代理指标，通常越大越优 |
| `HOMO(eV)` / `LUMO(eV)` | eV | 前线轨道位置 |
| `Gap(eV)` | eV | HOMO−LUMO gap |
| `Dipole(D)` | Debye | 极性代理指标 |
| `Gibbs(Eh)` | Hartree | 热力学结果；只宜在同协议下比较 |
| `Imaginary` | Yes/No | 最低点质量检查 |
| `charge` | 整数 | 形式电荷 |
| `LogP` | 无量纲 | 脂溶性 |
| `TPSA` | Å² | 极性表面积 |
| `Rotatable_Bonds` | 个 | 柔性 |
| `SA_Score` | 约 1–10 | 合成可及性；越低通常越容易合成 |
| `Molecular_Weight` | g/mol | 分子规模 |

不同泛函、基组、溶剂模型或电荷状态的绝对能量不可直接混合比较。

## 9. 分析与绘图脚本

所有批次都直接调用 `05_analysis/plotters/` 中的 Python 程序，并显式传入
数据和输出路径。克隆仓库后可直接使用
`data/core/pic16_s1/final/pic16inpdftv22.csv`；新计算使用
`data/processed/<batch>inpdftv22.csv`。历史固定批次白名单的包装脚本已
移除。

### 9.1 绘图指标和可调参数总表

| 脚本 | 展示内容 | 主要 CLI 参数 |
|---|---|---|
| `05_analysis/plotters/plot_boxplot_HL.py` | HOMO、LUMO、HOMO-LUMO gap、Dipole 的小提琴图/箱线图；HOMO-LUMO、gap-dipole、gap-S1T1、charge-dipole 关系 | `--columns`、`--output`、`--width`、`--height` |
| `05_analysis/plotters/plot_boxplot_s1t1.py` | S1/T1 能量分布；可叠加 S1 振子强度 | `--columns`、`--output`、`--width`、`--height`、`--enable-osc`、`--osc-column` |
| `05_analysis/plotters/plot_boxplot_s1_t1ti.py` | `S1-T1` 与最小正 `S1-Ti` 分布 | `--columns`、`--output`、`--width`、`--height`、`--filter-imaginary` |
| `05_analysis/plotters/plot_boxplot_socme.py` | 最大 SOCME 与最大 `SOCME/|ΔE|` 分布 | `--columns`、`--output`、`--width`、`--height` |
| `05_analysis/plotters/plot_pareto_s1ti_socme.py` | X=`min_positive_S1-Ti`，Y=SOCME；Pareto 前沿 | `--output`、`--width`、`--height`、`--xmin/xmax`、`--ymin/ymax`、`--show-pareto` |
| `05_analysis/plotters/plot_pareto_s1ti_socmegap.py` | X=`min_positive_S1-Ti`，Y=`SOCME/|ΔE|`；总图、各骨架图和边际分布 | 与上项相同 |
| `05_analysis/plotters/plot_ridge_socmev2.py` | SOCME 和 `SOCME/|ΔE|` 双栏山脊图，标注全量 median/max | `--columns`、`--output`、`--width`、`--height`、`--overlap`、`--scale`、`--xlabels`、`--xmax` |
| `05_analysis/plotters/plot_summary.py` | SA Score、MW、LogP、S1/波长的 KDE 与联合散点 | `--col`、`--output`；更多范围在代码中设置 |
| `05_analysis/plotters/plot_summary2.py` | X=`S1-T1`，Y=`T1/S1`，颜色=S1；参考集与生成集叠加 | 仅输入文件；坐标范围和色标在代码中设置 |
| `05_analysis/plotters/plot_targeting.py` | 电荷组成、正电比例、LogP-TPSA、电荷-TPSA、分面和组合图 | `--output-dir`、`--pictos`、`--filter-imaginary`、`--charge-max`、`--positive-threshold`、LogP/TPSA 窗口、`--save-pdf`、`--no-facets` |
| `05_analysis/plotters/plot_decorations_statistics.py` | Scaffold 中装饰片段频率和片段结构 | `--scaffold-col`、`--output`、`--top-n`、`--width`、`--height` |
| `05_analysis/plotters/plotsample.py` | Morgan 指纹化学空间；MaxMin 与随机采样对比；UMAP/t-SNE/MDS；NND 和来源采样率 | 样本数、SMILES 列、样本 CSV、绘图目录、映射表、投影方法、预选池、指纹、随机种子和降维参数 |

### 9.2 前线轨道和偶极矩

```bash
python3 05_analysis/plotters/plot_boxplot_HL.py \
  data/processed/pic3inpdftv22.csv \
  data/processed/pic16inpdftv22.csv \
  --output results/figures/frontier \
  --width 14 \
  --height 8
```

会以 `frontier` 为前缀生成：

- `frontier_violin.png`
- `frontier_box.png`
- `frontier_scatter.png`
- `frontier_gap_dipole.png`
- `frontier_gap_s1t1.png`
- `frontier_charge_dipole.png`

### 9.3 S1/T1 和振子强度

```bash
python3 05_analysis/plotters/plot_boxplot_s1t1.py \
  data/processed/pic16inpdftv22.csv \
  --enable-osc \
  --osc-column S1_fosc_D2 \
  --output results/figures/s1_t1.png
```

### 9.4 能隙和虚频过滤

```bash
python3 05_analysis/plotters/plot_boxplot_s1_t1ti.py \
  data/processed/pic16inpdftv22.csv \
  --filter-imaginary \
  --output results/figures/st_gap.png
```

`--filter-imaginary` 只保留 `Imaginary == No`。用于正式统计时建议启用。

### 9.5 Pareto 分析

```bash
python3 05_analysis/plotters/plot_pareto_s1ti_socmegap.py \
  data/processed/pic3inpdftv22.csv \
  data/processed/pic16inpdftv22.csv \
  --xmin 0 \
  --xmax 0.6 \
  --ymin 0 \
  --ymax 20 \
  --show-pareto pic3 pic16 \
  --output results/figures/pareto_soc_gap.png
```

Pareto 方向：

```text
X 越小越好：min_positive_S1-Ti
Y 越大越好：SOCME 或 SOCME/|ΔE|
```

轴范围只控制展示，不应在无记录的情况下用于删除原始数据。

### 9.6 Targeting 代理指标

```bash
python3 05_analysis/plotters/plot_targeting.py \
  data/processed/pic3inpdftv22.csv \
  data/processed/pic16inpdftv22.csv \
  --output-dir results/figures/targeting \
  --pictos config/pictos.txt \
  --filter-imaginary \
  --charge-max 6 \
  --positive-threshold 0 \
  --logp-min 1 \
  --logp-max 5 \
  --tpsa-min 20 \
  --tpsa-max 120 \
  --save-pdf
```

输出包括：

- 电荷组成热图；
- 正电荷比例/电荷组成；
- LogP-TPSA 平面和各电荷层；
- charge-TPSA 抖动散点；
- 各批次分面图；
- a–d 组合图；
- 部分额外的单批次和三维理化空间图。

`LogP 1–5`、`TPSA 20–120 Å²` 和 `charge > 0` 是当前脚本的可调代理窗口，不应写成普适生物学结论。

### 9.7 化学空间和采样

```bash
python3 05_analysis/plotters/plotsample.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --n_pick 20 \
  --method tsne \
  --pre_sample 50000 \
  --pictos config/pictos.txt \
  --output data/processed/pic16_s1_diverse20.csv \
  --plot-dir results/figures/pic16_s1/chemical_space
```

算法默认值和对应 CLI：

- Morgan 指纹半径 2、长度 2048：`--radius`、`--fp-size`；
- 随机种子 42：`--seed`；
- UMAP：Jaccard、默认邻居数 30、最小距离 0.1：
  `--umap-neighbors`、`--umap-min-dist`；
- t-SNE：先 PCA 到最多 50 维、默认 perplexity 30：`--tsne-perplexity`；
- MDS：Jaccard 距离，数据量大时内存和时间开销显著。

脚本把 MaxMin 选择结果写到 `--output`，并在 `--plot-dir` 下写入评估图和
投影坐标 CSV，供智能体计算覆盖度或继续聚类。仅在选择 `--method umap` 时
需要 `umap-learn`；查看帮助或使用 t-SNE/MDS 不依赖 UMAP。

## 10. 分子、骨架和结果展示

### 10.1 普通二维分子图

```bash
python3 06_visualization/generate_molecule_images.py \
  data/processed/pic16inpdftv22.csv \
  --tool rdkit \
  --output results/figures/molecules
```

可调：

- `--tool rdkit|openbabel`
- `--output`
- `--test INDEX`：只画某一行。

### 10.2 带性质标注的分子图

```bash
python3 06_visualization/generate_molecule_images_with_props.py \
  data/processed/pic16inpdftv22.csv \
  --tool rdkit \
  --output results/figures/properties \
  --verbose
```

图中展示：

- HOMO-LUMO gap；
- S1、T1；
- `S1-T1`；
- `min(S1-Ti)`；
- SOCME；
- `SOCME/δE`；
- 分子 ID、charge、S1 波长、S1 振子强度等可用字段。

可调参数：`--tool`、`--output`、`--test`、`--verbose`。

### 10.3 相对参考分子筛选 Top N

```bash
python3 06_visualization/generate_top_selected_images.py \
  data/processed/pic16inpdftv22.csv \
  --reference-csv 05_analysis/data/publishinpdftv22.csv \
  --top-n 20 \
  --tool rdkit \
  --output-dir results/figures/top_selected
```

筛选条件：

- 候选 `min(S1-Ti)` 小于对应参考分子；
- 候选 `SOCME/δE` 大于对应参考分子；
- 再按脚本排序规则选前 N。

### 10.4 分子、骨架和装饰片段

```bash
python3 06_visualization/scaffold/generate_mol_and_scaffold_figs.py \
  --csv data/processed/pic16inpdftv22.csv \
  --row 1 \
  --output-dir results/figures/scaffold \
  --w 900 \
  --h 700 \
  --show-properties
```

输出：

- 整分子图；
- scaffold 图；
- 各 attachment point 装饰片段图；
- 可选的量化性质文字。

可调：`--smiles` 或 `--csv`、`--scaffold`、列名、行号、输出名、宽高和 `--show-properties`。

`06_visualization/scaffold/generate_mol_and_scaffold_figsv2.py` 使用更宽松的骨架匹配，最初用于通用匹配失败的 pic4；不应在未检查匹配结果时自动替代主版本。

### 10.5 Excel 图片表

```bash
python3 06_visualization/execlwithimage.py \
  --csv data/processed/pic16inpdftv22.csv \
  --output results/pic16inpdftv23.xlsx \
  --smiles-col SMILES \
  --img-col-name Structure \
  --img-width 200 \
  --img-height 200
```

该脚本使用 Open Babel 生成结构图并插入 Excel。`--img-height` 主要控制 Excel 行高，实际图片高度由 Open Babel 比例决定。

## 11. 质量控制和量化运行指标

建议每个批次至少记录：

| 指标 | 定义 |
|---|---|
| 输入数 | 原始 CSV 数据行数 |
| 基态正常结束率 | 含 `ORCA TERMINATED NORMALLY` 的基态日志数 / 输入数 |
| 基态解析覆盖率 | `opt.csv` 行数 / 基态正常日志数 |
| 无虚频率 | `Imaginary == No` 数 / 基态解析数 |
| TD 正常结束率 | 正常 TD 日志数 / 基态正常数 |
| TD 解析覆盖率 | TD 汇总行数 / 正常 TD 日志数 |
| 完整 v22 率 | 核心字段均非空的 v22 行数 / 输入数 |
| Pareto 数量和比例 | 非支配候选数 / 完整候选数 |
| 目标窗口命中率 | 满足用户定义能隙、SOC、波长和药化约束的比例 |
| 多样性 | MaxMin 样本最近邻距离、指纹覆盖或聚类覆盖 |

生成机器可读报告：

```bash
python3 07_validation/collect_metrics.py
```

新批次应显式提供每个阶段的数据、身份和同协议参考分子：

```bash
python3 07_validation/collect_metrics.py \
  --source data/raw/<batch>.csv \
  --dataset data/processed/<batch>inpdftv22.csv \
  --td-summary work/<batch>inpout/<batch>inpout_tdv2.csv \
  --v20 work/<batch>inpout/<batch>inpdftv20.csv \
  --v21 data/processed/<batch>inpdftv21.csv \
  --reference 05_analysis/data/publishinpdftv22.csv \
  --reference-index <reference-id> \
  --dataset-id <batch> \
  --display-label <label> \
  --output reports/<batch>/metrics.json
```

没有同一计算协议下的合适参考分子时，使用 `--no-reference`，不要将不同泛函、
基组、溶剂或相对论方案的数值直接比较。

当前 `reports/pic16_s1/metrics.json` 快照显示：

- 输入、TD 汇总、v20、v21、v22 均为 222 行，v22 核心字段完整率 100%；
- 无虚频 191 条，无虚频率约 86.0%；
- gap–`SOCME/|ΔE|` 二目标 Pareto 为 2 条；
- gap–SOCME–ratio 三目标 Pareto 为 14 条；
- 相对 `05_analysis/data/publishinpdftv22.csv` 中的 `pic16` 参考行，
  gap 与 ratio 同时达到或优于参考的候选为 175/222；排除虚频后为 156/191；
- 再要求 SOCME 本身达到或优于参考时为 155/222；排除虚频后为 138/191；
- 默认 LogP 1–5 与 TPSA 20–120 Å² 联合窗口命中 0 条。

该 `pic16` 参考行自身为 `Imaginary=Yes`，因此上述计数只表示电子性质阈值
比较，不能解释为相对一个已确认几何极小值的性能提升。

最后一项说明默认 targeting 窗口不适合 S1 数据分布，不能把该绘图默认值当成
通用筛选标准。S1 的 LogP 中位数约 14.08，应先根据具体应用、分子尺寸和参考
体系重新定义窗口并做敏感性分析。

优化顺序建议：

1. 先修复正常结束率和解析覆盖率；
2. 再排除虚频、错电荷、异常结构和不完整记录；
3. 再做电子性质 Pareto；
4. 加入 LogP/TPSA/SA/MW 等可行性约束；
5. 最后通过化学空间采样控制结构冗余；
6. 对 Pareto 或目标窗口中的代表分子使用更高层级方法复算。

## 12. 已知限制

1. 当前没有统一的 Python 包或单一工作流引擎，阶段之间通过文件命名和 CSV 连接。
2. `03_result_parsing/parse_orca_opt.py` 不会单独强制验证正常结束，自动化代理必须先检查日志。
3. TD 解析器依赖 CSV 行序与任务编号一致。
4. `03_result_parsing/merge_dft_csv.py` 使用两个表的 ID 并集，可能产生只有基态或只有 TD 数据的部分行。
5. 部分绘图程序存在代码内固定坐标范围、DPI、调色板和随机种子；CLI 未暴露时需要修改源文件。
6. `config/agent_interface.json` 是声明式接口，不会自动改写四个 ORCA 生成脚本；生成脚本是执行权威。
7. 当前参考环境使用 ORCA 6.0.1。不同 ORCA 大版本的输出格式可能导致正则解析失效。
8. `SA_Score` 使用随仓库提供的 RDKit SA Score 数据，第三方说明见 [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)。
9. 环境和指标报告可以自动生成，但每次计算的完整 provenance 仍需写入 `reports/runs/`。

## 13. 验证

Python 语法：

```bash
find . -type f -name '*.py' \
  -exec python3 -c 'import pathlib,sys; p=pathlib.Path(sys.argv[1]); compile(p.read_text(encoding="utf-8"), str(p), "exec")' {} \;
```

Shell 语法：

```bash
find . -type f -name '*.sh' -exec bash -n {} \;
```

SA Score 测试：

```bash
python3 07_validation/UnitTestSAScore.py
```

接口、环境和核心数据：

```bash
python3 workflow.py paths --check
python3 07_validation/check_environment.py --output reports/environment.json
sha256sum -c data/core/pic16_s1/SHA256SUMS
python3 07_validation/collect_metrics.py
```

绘图烟雾测试：

```bash
python3 05_analysis/plotters/plot_pareto_s1ti_socmegap.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --output results/figures/pic16_s1_pareto.png
```

## 14. GitHub 发布

```bash
git init
git remote add origin https://github.com/letsgo-2/meloculardftcal.git
git add .
git status
git commit -m "Initial public workflow release"
git branch -M main
git push -u origin main
```

公开前必须：

- 选择并添加合适的 `LICENSE`；
- 确认参考分子和数据允许公开；
- 检查用户名、账号、令牌和集群绝对路径；
- 不上传 ORCA 安装包和受限软件；
- 不上传完整 GBW、日志、批量图片和敏感研究数据；
- 用小批次重新跑通完整链路；
- 保存并检查计算协议和运行清单。
