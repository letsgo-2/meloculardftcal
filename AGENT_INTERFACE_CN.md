# meloculardftcal 智能体接口（中文）

本文件是自动化智能体准备、运行、检查、分析和可视化本仓库时应遵守的
操作契约。

相关入口：

- 中文说明：[`README_CN.md`](README_CN.md)
- English README：[`README_EN.md`](README_EN.md)
- English agent interface：[`AGENT_INTERFACE_EN.md`](AGENT_INTERFACE_EN.md)
- GitHub 仓库：<https://github.com/letsgo-2/meloculardftcal>

机器可读接口：

- 阶段、版本、路径和参数注册表：`config/agent_interface.json`
- 只读查询命令：`workflow.py`
- 当前运行环境报告：`reports/environment.json`
- pic16/S1 数据和性能报告：`reports/pic16_s1/metrics.json`

智能体在构造命令前应先查询接口：

```bash
python3 workflow.py stages
python3 workflow.py stage S41_PARSE_TD
python3 workflow.py software
python3 workflow.py parameters
python3 workflow.py paths --check
```

所有命令输出 JSON。`config/agent_interface.json` 是声明式注册表；真正执行
计算、解析或绘图时，以注册表指向的脚本为准。

## 1. 范围和不变量

1. 所有命令从仓库根目录运行。
2. 项目内路径必须按仓库相对路径解析，不能假定仓库位于固定绝对路径。
3. `ORCA_BIN`、`MODULE_INIT` 和 `WORK_DIR` 是集群外部路径，只能写入被
   Git 忽略的 `config/cluster.env`。
4. 不得提交 ORCA 程序、MPI 安装包、GBW、大型日志、Slurm 输出、私有数据、
   密钥、令牌和批量图片。
5. `data/raw/`、`data/processed/`、`work/` 和 `results/figures/` 是运行时
   或生成数据目录。
6. `data/core/pic16_s1/`、`05_analysis/data/publishinpdftv22.csv` 以及
   `reports/` 下的小型 JSON 是有意纳入版本控制的参考数据。
7. 生成任务后不得改变源 CSV 行序。TD 解析器将 `*_<id>td.log` 中的编号
   映射到源 CSV 的第 `id - 1` 行。
8. 不得在没有记录协议差异的情况下，比较不同泛函、基组、溶剂、相对论方案
   或 ORCA 大版本产生的数值。
9. 运行成功不能只根据进程返回码判断，必须同时应用阶段质量门控。

## 2. 软件、环境和机器接口

### 2.1 参考软件版本

| 组件 | 要求或已记录值 |
|---|---|
| Python | `>=3.9` |
| ORCA 解析目标 | ORCA `6.0.1` 输出 |
| ORCA 参考构建 | `linux_x86-64_shared_openmpi416_avx2` |
| MPI | OpenMPI `4.1.6` |
| 通信层 | UCX `1.15.0` |
| 调度器 | Slurm |
| 三维结构/图片后端 | Open Babel `obabel` |
| 化学信息学 | RDKit |
| 数据和绘图 | NumPy、pandas、SciPy、Matplotlib、Seaborn |
| 化学空间 | scikit-learn；UMAP 方法额外需要 `umap-learn` |

推荐环境：

```bash
conda env create -f environment.yml
conda activate dft-workflow
```

运行前检查：

```bash
python3 07_validation/check_environment.py \
  --output reports/environment.json
python3 workflow.py paths --check
```

集群上还应检查：

```bash
source config/cluster.env
"${ORCA_BIN}" --help
sbatch --version
squeue --version
```

智能体遇到能力缺失时必须报告，不得静默更换 ORCA 版本、计算方法、基组或
溶剂模型。

### 2.2 参考硬件和规模

历史节点配置：

- 2 × AMD EPYC 7H12；
- 128 个物理核心、256 个逻辑 CPU；
- 约 503 GiB 内存；
- 不需要 GPU。

历史多批次输入约 2,998 个分子。仓库只附带 pic16/S1 的 222 条核心数据。

默认任务资源：

| 配置 | 节点 | MPI ranks | 每 rank OMP 线程 | ORCA `%maxcore` | 估算声明内存 |
|---|---:|---:|---:|---:|---:|
| `opt_standard` | 1 | 64 | 1 | 3000 MB/rank | 187.5 GiB |
| `opt_heavy_zora` | 1 | 64 | 1 | 3000 MB/rank | 187.5 GiB |
| `td_standard` | 1 | 32 | 1 | 3000 MB/rank | 93.75 GiB |
| `td_heavy_zora` | 1 | 32 | 1 | 3000 MB/rank | 93.75 GiB |

估算值为 `nprocs × maxcore`，实际峰值可能更高。增加 rank 或 `%maxcore`
前必须复核节点内存、调度策略、临时空间、分子规模和 ORCA 并行效率。

## 3. 集群配置接口

创建私有配置：

```bash
cp config/cluster.env.example config/cluster.env
```

`config/cluster.env` 已被 `.gitignore` 排除。

| 变量 | 模板值 | 使用者 | 含义 |
|---|---|---|---|
| `SBATCH_PARTITION` | `amd_512q` | 四个 ORCA 生成脚本 | Slurm 分区 |
| `SBATCH_NODES` | `1` | 四个生成脚本 | 每任务节点数 |
| `SBATCH_NTASKS` | `32` | 四个生成脚本 | MPI ranks 和 ORCA `nprocs` |
| `WORK_DIR` | 站点相关 | 生成的 `.sh` | 集群仓库/工作目录 |
| `ORCA_BIN` | 站点相关 | 生成的 `.sh` | ORCA 可执行文件 |
| `MODULE_INIT` | 站点相关 | 生成的 `.sh` | Modules 初始化脚本 |
| `ORCA_MODULE` | ORCA 6.0.1 模块名 | 生成的 `.sh` | 要加载的软件模块 |

基态生成器在未设置 `SBATCH_NTASKS` 时默认 64 ranks，TD 生成器默认
32 ranks。若 `config/cluster.env` 设置了共享值，应在每个阶段生成前显式覆盖：

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh \
  work/pic16 work/pic16inp

export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh \
  work/pic16 work/pic16inp
```

生成的提交脚本固定写入 `OMP_NUM_THREADS=1`，当前未通过
`config/cluster.env` 开放。

## 4. 稳定阶段注册表

| 阶段 ID | 入口 | 主要输入 | 主要输出 | 外部副作用 |
|---|---|---|---|---|
| `S00_PREFLIGHT` | `07_validation/check_environment.py` | 环境与配置 | 环境 JSON | 无 |
| `S10_STRUCTURE` | `01_structure_preparation/convertInputs.py` | 源 CSV | 编号 GJF | 写工作目录 |
| `S20_OPT_GENERATE` | `gjf_to_inp.sh` / `gjf_to_inpp13.sh` | GJF、集群配置 | 基态 `.inp/.sh` | 写工作目录 |
| `S21_OPT_SUBMIT` | `02_cluster_execution/suboptsbatch.sh` | 基态任务、ID | Slurm 任务 | 提交外部任务 |
| `S22_OPT_QC` | `find_missing_opt_ids.sh` | 基态日志 | 缺失 ID/QC | 无 |
| `S30_TD_GENERATE` | `gjf_to_inptd.sh` / `gjf_to_inptdp13.sh` | GJF、XYZ、GBW | TD `.inp/.sh` | 写工作目录 |
| `S31_TD_SUBMIT` | `02_cluster_execution/subtdsbatch.sh` | TD 任务、ID | Slurm 任务 | 提交外部任务 |
| `S32_TD_QC` | `find_missing_td_ids.sh` | TD 日志 | 缺失 ID/QC | 无 |
| `S40_PARSE_OPT` | `03_result_parsing/parse_orca_opt.py` | 基态日志 | 基态 CSV | 写解析目录 |
| `S41_PARSE_TD` | `03_result_parsing/parse_orca_td_socv2.py` | TD 日志、源 CSV | TD/SOC CSV | 写解析目录 |
| `S42_MERGE_DFT` | `03_result_parsing/merge_dft_csv.py` | 基态/TD 汇总 | v20 | 写解析目录 |
| `S50_MERGE_SOURCE` | `04_data_processing/merge_summary.py` | 源 CSV、v20 | v21 | 写处理目录 |
| `S51_ENRICH` | `04_data_processing/addrdkit.py` | v21 | v22 | 写处理目录 |
| `S52_SYMMETRY` | `04_data_processing/addcenter.py` | 含骨架的 CSV | 对称性字段 | 写指定 CSV |
| `S60_METRICS` | `07_validation/collect_metrics.py` | 各阶段 CSV | 指标 JSON | 写 `reports/` |
| `S70_ANALYZE` | `05_analysis/plotters/` | v22 | 图、筛选 CSV | 写指定输出 |
| `S80_VISUALIZE` | `06_visualization/` | CSV/SMILES | 分子和骨架图片 | 写指定输出 |

调用阶段前使用以下命令读取完整参数模板：

```bash
python3 workflow.py stage <STAGE_ID>
```

`S21_OPT_SUBMIT` 和 `S31_TD_SUBMIT` 会改变外部调度器状态，必须获得明确的
任务提交权限。其他标记为无外部副作用的阶段仍可能写入仓库运行目录。

## 5. 标准工作流命令

以下 `<batch>` 使用稳定批次名，例如 `pic16`。

### 5.1 结构生成

```bash
python3 01_structure_preparation/convertInputs.py \
  1 data/raw/<batch>.csv work/<batch>
```

源 CSV 至少包含 `SMILES`；骨架分析还需要 `Scaffold`。必须显式传入输出
目录，以避免兼容模式把 GJF 写到仓库根目录。

### 5.2 基态生成、提交和检查

普通元素：

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh \
  work/<batch> work/<batch>inp
bash 02_cluster_execution/suboptsbatch.sh \
  work/<batch>inp 1-20
```

重元素/ZORA：

```bash
bash 01_structure_preparation/gjf_to_inpp13.sh \
  work/<batch> work/<batch>inp
```

ID 可以是单值或闭区间：

```bash
bash 02_cluster_execution/suboptsbatch.sh \
  work/<batch>inp 1-10 20 30-35
```

只读检查：

```bash
bash 02_cluster_execution/find_missing_opt_ids.sh \
  work/<batch>inp
```

### 5.3 TD/SOC 生成、提交和检查

仅当对应基态任务正常结束且 `.xyz/.gbw` 完整时继续：

```bash
source config/cluster.env
export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh \
  work/<batch> work/<batch>inp
bash 02_cluster_execution/subtdsbatch.sh \
  work/<batch>inp 1-20
bash 02_cluster_execution/find_missing_td_ids.sh \
  work/<batch>inp
```

重元素体系使用 `01_structure_preparation/gjf_to_inptdp13.sh`。

### 5.4 解析、合并和描述符

```bash
python3 03_result_parsing/parse_orca_opt.py \
  work/<batch>inp \
  --output work/<batch>inpout

python3 03_result_parsing/parse_orca_td_socv2.py \
  work/<batch>inp \
  --csv data/raw/<batch>.csv \
  --output work/<batch>inpout

python3 03_result_parsing/merge_dft_csv.py \
  work/<batch>inp \
  --output work/<batch>inpout

python3 04_data_processing/merge_summary.py \
  --piccsvfile data/raw/<batch>.csv \
  --summary work/<batch>inpout/<batch>inpdftv20.csv \
  --output data/processed/<batch>inpdftv21.csv

python3 04_data_processing/addrdkit.py \
  data/processed/<batch>inpdftv21.csv \
  --output data/processed/<batch>inpdftv22.csv
```

对 `work/<batch>inp`，合并器从 `work/<batch>inpout/` 读取
`<batch>inpout_opt.csv` 和 `<batch>inpout_tdv2.csv`。

## 6. 量化计算协议和参数修改

### 6.1 当前协议

| 配置 | ORCA 关键词/控制 |
|---|---|
| 普通 Opt/Freq | `B3LYP def2-SVP Opt Freq RIJCOSX CPCM(Water)` |
| 重元素 Opt/Freq | `PBE0 D3BJ ZORA ZORA-def2-TZVP SARC/J Opt Freq RIJCOSX TightSCF CPCM(Water)` |
| 普通 TD/SOC | `wB97X-D3 def2-TZVP RIJCOSX TightSCF CPCM(Water)`；`NRoots 10`、`Triplets true`、`dosoc true`、`MORead` |
| 重元素 TD/SOC | `wB97X-D3 ZORA ZORA-def2-TZVP SARC/J RIJCOSX RI-SOMF(1X) TightSCF CPCM(Water)`，TD 控制同上 |

TD 任务读取收敛的基态 GBW 和优化 XYZ。电荷与多重度来自 GJF 的
charge/multiplicity 行。

### 6.2 参数修改接口

`config/agent_interface.json` 公开当前结构化协议，但不会自动改写 ORCA 输入。
修改下列执行脚本后必须重新生成 `.inp/.sh`。

| 参数 | 普通体系文件 | 重元素文件 | 修改位置 |
|---|---|---|---|
| 基态方法、基组、溶剂 | `gjf_to_inp.sh` | `gjf_to_inpp13.sh` | 输出的 `!` 关键字行 |
| TD 方法、基组、溶剂 | `gjf_to_inptd.sh` | `gjf_to_inptdp13.sh` | 输出的 `!` 关键字行 |
| TD 根数 | 两个 TD 生成器 | 两个 TD 生成器 | `NRoots 10` |
| 三重态 | 两个 TD 生成器 | 两个 TD 生成器 | `Triplets true` |
| SOC | 两个 TD 生成器 | 两个 TD 生成器 | `dosoc true` |
| 每 rank 内存 | 四个生成器 | 四个生成器 | `%maxcore 3000` |
| MPI ranks | `config/cluster.env` | `config/cluster.env` | `SBATCH_NTASKS` |
| 分区/节点 | `config/cluster.env` | `config/cluster.env` | `SBATCH_PARTITION`、`SBATCH_NODES` |
| 电荷/多重度 | 源 GJF | 源 GJF | charge/multiplicity 行 |

变更流程：

1. 创建 Git 分支或复制当前协议生成器。
2. 基准测试时一次只改变一个协议维度。
3. 在新目录重新生成输入，不得混合旧、新 `.inp`。
4. 人工检查至少一个新 `.inp` 和 `.sh`。
5. 先计算小型参考子集，再提交完整批次。
6. 在运行清单记录 ORCA/MPI 版本、完整关键字、资源、Git commit 和输入哈希。
7. 分析时隔离不兼容协议的数据。

## 7. 调度控制和安全

只读状态查询：

```bash
bash 02_cluster_execution/find_missing_opt_ids.sh work/<batch>inp
bash 02_cluster_execution/find_missing_td_ids.sh work/<batch>inp
bash 02_cluster_execution/list_opt_jobs.sh <batch>
bash 02_cluster_execution/list_td_jobs.sh <batch>
```

持续提交：

```bash
bash 02_cluster_execution/auto_submit_opt.sh \
  work/<batch>inp 20 1800
bash 02_cluster_execution/auto_submit_td.sh \
  work/<batch>inp 20 1800
```

两个数字参数为 `MAX_JOBS` 和 `SLEEP_SEC`。`MAX_JOBS` 统计当前 Slurm 用户
的全部任务，而不仅是选定批次。持续提交会长期运行并改变外部调度器状态，
智能体必须获得明确权限。

历史批量超时终止脚本已从仓库删除。取消任务必须使用经复核的明确 job ID，
并单独取得执行 `scancel` 的权限。

## 8. 文件、路径和字段契约

### 8.1 命名

TD 日志命名必须匹配 `*_<id>td.log`，其中 `<id>` 为正整数。基态文件必须
保留相同批次前缀和编号。批次名不要包含空格。

| 名称 | 含义 |
|---|---|
| `<batch>inpdftv20.csv` | 合并基态和 TD 量化结果 |
| `<batch>inpdftv21.csv` | v20 加回源生成字段和稳定索引 |
| `<batch>inpdftv22.csv` | v21 加 RDKit 描述符和 SA Score |

### 8.2 可直接查看的数据和结果路径

| 内容 | 仓库相对路径 |
|---|---|
| pic16/S1 源分子 | `data/core/pic16_s1/input/pic16.csv` |
| TD/SOC 汇总 | `data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv` |
| v20 | `data/core/pic16_s1/intermediate/pic16inpdftv20.csv` |
| v21 | `data/core/pic16_s1/intermediate/pic16inpdftv21.csv` |
| v22 | `data/core/pic16_s1/final/pic16inpdftv22.csv` |
| 性能参考表 | `05_analysis/data/publishinpdftv22.csv` |
| 核心数据校验和 | `data/core/pic16_s1/SHA256SUMS` |
| 数据来源说明 | `data/core/pic16_s1/DATASET.md` |
| 量化性能快照 | `reports/pic16_s1/metrics.json` |
| 当前环境检查 | `reports/environment.json` |
| 生成图片 | `results/figures/<batch>/` |
| 单次运行清单 | `reports/runs/<run_id>.json` |

使用以下命令机器化发现路径：

```bash
python3 workflow.py paths --check
```

仓库未附带原始 ORCA 日志和波函数，因此核心指标报告中的正常结束率为
`null`，不能从提交的数据反推。

### 8.3 主要量化字段

| 字段 | 单位/类型 | 目标或作用 |
|---|---|---|
| `S1_energy_eV` | eV | 目标/窗口 |
| `T1_energy_eV` | eV | 上下文 |
| `S1_wavelength_nm` | nm | 目标/窗口 |
| `S1_fosc_D2` | 无量纲 | 通常最大化 |
| `S1-T1` | eV | 通常最小化绝对或正能隙 |
| `min_positive_S1-Ti` | eV | 最小化 |
| `max_SOCME_S1_Ti_value` | cm^-1 | 最大化 |
| `max_SOCME_over_deltaE_S1_Ti_value` | cm^-1/eV | 最大化 |
| `HOMO(eV)` / `LUMO(eV)` | eV | 目标/上下文 |
| `Gap(eV)` | eV | 目标/上下文 |
| `Dipole(D)` | Debye | 目标/上下文 |
| `Gibbs(Eh)` | Hartree | 只在同协议内比较 |
| `Imaginary` | `Yes`/`No` | 极小值要求 `No` |
| `LogP` | 无量纲 | 可行性窗口 |
| `TPSA` | Å² | 可行性窗口 |
| `Rotatable_Bonds` | 数量 | 柔性约束 |
| `SA_Score` | 约 1–10 | 最小化 |
| `Molecular_Weight` | g/mol | 尺寸约束 |

SOCME 是 X/Y/Z 复分量模长，ratio 定义为
`SOCME / abs(delta_E_eV)`。状态筛选规则以
`03_result_parsing/parse_orca_td_socv2.py` 为准。

## 9. 质量门控和性能指标

每个批次至少记录：

```text
N_input            = 源 CSV 行数
N_opt_normal       = 正常结束的基态日志数
N_opt_parsed       = 基态汇总行数
N_no_imaginary     = Imaginary == "No" 行数
N_td_normal        = 正常结束的 TD 日志数
N_td_parsed        = TD 汇总行数
N_complete_v22     = 所有核心字段完整的 v22 行数
```

建议派生率：

```text
opt_success_rate       = N_opt_normal / N_input
opt_parse_rate         = N_opt_parsed / N_opt_normal
minimum_rate           = N_no_imaginary / N_opt_parsed
td_success_rate        = N_td_normal / N_no_imaginary
td_parse_rate          = N_td_parsed / N_td_normal
pipeline_complete_rate = N_complete_v22 / N_input
```

分母为零或数据不可获得时使用 `null` 并给出原因，不得用零替代缺失值。

质量链：

```text
源行和结构有效
  -> 基态日志正常结束
  -> 优化 XYZ 和 GBW 存在
  -> 频率已解析且 Imaginary == "No"
  -> TD 日志正常结束
  -> S1、三重态、振子强度和 SOC 已解析
  -> 源数据与 DFT 一对一合并
  -> RDKit 描述符有效
  -> 候选可进入排序
```

生成默认 pic16/S1 报告：

```bash
python3 07_validation/collect_metrics.py
```

新批次：

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

没有同协议参考分子时使用 `--no-reference`。不同量化协议的参考值不得直接
比较。

当前 pic16/S1 快照：

- 输入、TD、v20、v21、v22 均为 222 行；
- 222 行核心字段完整；
- 无虚频 191 行，约 86.0%；
- gap–ratio 二目标 Pareto 为 2 行；
- gap–SOCME–ratio 三目标 Pareto 为 14 行；
- 相对 `pic16` 参考，gap 和 ratio 达标 175/222，去虚频后 156/191；
- 再要求 SOCME 达标为 155/222，去虚频后 138/191；
- 默认 LogP 1–5、TPSA 20–120 Å² 联合窗口命中 0 行。

参考表中的 `pic16` 行自身为 `Imaginary=Yes`，因此比较结果只是电子性质
阈值，不能视为相对已验证几何极小值的性能提升。

## 10. 优化目标接口

建议智能体使用显式目标结构：

```yaml
objective_version: 1
quality_constraints:
  orca_opt_normal: true
  imaginary: "No"
  orca_td_normal: true
objectives:
  - column: min_positive_S1-Ti
    direction: minimize
  - column: max_SOCME_S1_Ti_value
    direction: maximize
  - column: max_SOCME_over_deltaE_S1_Ti_value
    direction: maximize
optional_objectives:
  - column: S1_fosc_D2
    direction: maximize
feasibility_windows:
  LogP: [1.0, 5.0]
  TPSA: [20.0, 120.0]
targets:
  S1_wavelength_nm: null
constraints:
  Molecular_Weight: null
  SA_Score: null
diversity:
  fingerprint: Morgan
  selection: MaxMin
```

LogP/TPSA 是当前绘图默认值，不是普适科学约束。智能体必须公开所有阈值并
执行敏感性分析。

建议下一步策略：

1. 优先修复正常结束率和解析覆盖率。
2. 从科学排序中排除虚频和不完整记录。
3. 使用 Pareto 集，不得隐藏为未记录的单一分数。
4. 添加应用相关的波长、电荷、LogP、TPSA、MW 和 SA 约束。
5. 使用 `05_analysis/plotters/plotsample.py` 选择结构多样子集。
6. 对代表候选使用更高层级或替代协议复算。
7. 比较排序稳定性和不确定性后再提出新分子。

## 11. 分析和可视化接口

### 11.1 绘图入口

| 脚本 | 展示内容 | 可调参数 |
|---|---|---|
| `plot_boxplot_HL.py` | HOMO、LUMO、Gap、Dipole | 列、输出、宽高 |
| `plot_boxplot_s1t1.py` | S1/T1 和可选振子强度 | 列、输出、尺寸、振子选项 |
| `plot_boxplot_s1_t1ti.py` | `S1-T1`、最小正 `S1-Ti` | 列、输出、尺寸、虚频过滤 |
| `plot_boxplot_socme.py` | SOCME 和 `SOCME/abs(deltaE)` | 列、输出、尺寸 |
| `plot_pareto_s1ti_socme.py` | gap–SOCME Pareto | 输出、尺寸、轴范围、曲线选择 |
| `plot_pareto_s1ti_socmegap.py` | gap–SOCME/gap Pareto | 输出、尺寸、轴范围、曲线选择 |
| `plot_ridge_socmev2.py` | SOCME 与 ratio 山脊分布 | 列、输出、尺寸、重叠、缩放 |
| `plot_summary.py` | SA、MW、LogP、S1、波长 | 输入、SMILES 列、输出；部分源码常量 |
| `plot_summary2.py` | `S1-T1` 对 `T1/S1`，按 S1 着色 | 输入；部分源码轴和色标 |
| `plot_targeting.py` | 电荷、LogP–TPSA、charge–TPSA | 输出目录、过滤、窗口、PDF、分面 |
| `plot_decorations_statistics.py` | 装饰片段数量和结构 | Scaffold 列、Top N、输出尺寸 |
| `plotsample.py` | 指纹空间、MaxMin/随机、NND、来源采样率 | 样本数、列、输出、指纹、随机种子、降维 |

查看真实 CLI：

```bash
python3 05_analysis/plotters/<script>.py --help
```

只在源码中开放的坐标范围或图形常量，修改时必须记录具体变量。裁剪坐标轴
不能解释为过滤了数据。

### 11.2 标准分析命令

```bash
python3 05_analysis/plotters/plot_pareto_s1ti_socmegap.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --output results/figures/pic16_s1/pareto_ratio.png

python3 05_analysis/plotters/plot_targeting.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --output-dir results/figures/pic16_s1/targeting \
  --pictos config/pictos.txt \
  --filter-imaginary \
  --logp-min 1 --logp-max 5 \
  --tpsa-min 20 --tpsa-max 120

python3 05_analysis/plotters/plotsample.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --n_pick 20 \
  --method tsne \
  --pictos config/pictos.txt \
  --output data/processed/pic16_s1_diverse20.csv \
  --plot-dir results/figures/pic16_s1/chemical_space
```

`plotsample.py` 将 MaxMin 结果写到 `--output`，并在 `--plot-dir` 写入投影
坐标 CSV 和评估图。UMAP 仅在 `--method umap` 时加载。

### 11.3 分子和骨架图

| 脚本 | 用途 | 主要参数 |
|---|---|---|
| `06_visualization/generate_molecule_images.py` | 普通二维分子图 | CSV、SMILES 列、后端、输出 |
| `generate_molecule_images_with_props.py` | 带性质标注的结构图 | CSV、RDKit/Open Babel、输出、单行测试 |
| `generate_top_selected_images.py` | 相对参考筛选 Top N 并绘图 | 候选、参考、Top N、后端、输出 |
| `execlwithimage.py` | 生成嵌入图片的 Excel | CSV、图片目录、输出 |
| `scaffold/generate_mol_and_scaffold_figsv2.py` | 分子、核心和装饰片段 | 行/列、输出、宽高 |
| `scaffold/symmetry_and_drawv2.py` | 连接点对称性和绘图 | 手性、匹配/移动限制、输出前缀 |

图片输出默认放到 `results/figures/<batch>/`，不得把大规模生成图提交到 Git。

## 12. 运行来源清单

建议每次计算写入：

```text
reports/runs/<run_id>.json
```

至少包括：

```yaml
schema_version: 1
run_id: "<UTC timestamp or stable ID>"
git_commit: "<commit>"
batch: "<batch>"
stage_id: "<stage>"
software:
  orca: "6.0.1"
  mpi: "OpenMPI 4.1.6"
hardware:
  cpu_model: "<model>"
  nodes: 1
resources:
  mpi_ranks: 64
  omp_threads_per_rank: 1
  maxcore_mb_per_rank: 3000
protocol:
  method: "B3LYP"
  basis: "def2-SVP"
  solvent: "CPCM(Water)"
inputs: []
outputs: []
input_sha256: {}
quality_metrics: {}
warnings: []
```

不得在清单中保存用户名、凭据、私有主机名或授权软件绝对路径。

## 13. 已知自动化缺口

1. 量化关键词仍嵌入四个 Shell 生成器，没有统一写入结构化协议配置。
2. 部分绘图轴范围、标签和图片尺寸仍是源码常量。
3. `convertInputs.py` 使用位置参数，没有完整 argparse 接口。
4. TD–SMILES 合并仍依赖源行序和文件名 ID。
5. ORCA 解析器针对 6.0.1 文本格式，支持新大版本前需要回归测试。
6. 基态解析器本身不强制正常结束标记。
7. 合并步骤可能保留部分记录。
8. `config/agent_interface.json` 是只读注册表，不是工作流引擎或参数写入器。
9. 多数科研绘图仍通过 stdout 和图片报告；稳定 JSON 边界由
   `workflow.py` 和 `collect_metrics.py` 提供。

优先改进方向是统一协议配置、非交互结构生成、显式 join key、解析器测试
夹具，以及继续把绘图常量开放为 CLI。

## 14. 智能体完成报告

每次自主运行应返回：

```yaml
status: success | partial | failed
stage_id: "<stage>"
commands_executed: []
inputs: []
outputs: []
parameters_changed: {}
quality_metrics: {}
warnings: []
recommended_next_action: null
```

只有命令返回码为零不能标记为 `success`。必须验证该阶段预期文件并应用对应
质量门控。
