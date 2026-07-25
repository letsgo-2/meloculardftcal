# meloculardftcal

[中文 README](README_CN.md) | [English README](README_EN.md) | [中文智能体接口](AGENT_INTERFACE_CN.md) | [English agent interface](AGENT_INTERFACE_EN.md)

Repository: <https://github.com/letsgo-2/meloculardftcal>

This repository organizes a high-throughput molecular workflow for ORCA
ground-state DFT, TD-DFT/SOC calculations, Slurm execution, result parsing,
RDKit enrichment, quantitative analysis, and molecular/scaffold visualization.

本仓库整理了以 ORCA 为核心的分子高通量工作流，包括基态 DFT、
TD-DFT/SOC、Slurm 集群运行、结果解析、RDKit 描述符、量化分析以及分子和
骨架可视化。

```text
SMILES/CSV
  -> 3D structure
  -> Opt + Freq
  -> TD-DFT + SOC
  -> parsed v20/v21/v22 datasets
  -> QC, Pareto analysis, chemical-space sampling, figures
```

Start with [`README_CN.md`](README_CN.md) or
[`README_EN.md`](README_EN.md). Autonomous tools should read
[`AGENT_INTERFACE_CN.md`](AGENT_INTERFACE_CN.md) or
[`AGENT_INTERFACE_EN.md`](AGENT_INTERFACE_EN.md) before changing parameters,
submitting jobs, or ranking candidates. [`AGENT_INTERFACE.md`](AGENT_INTERFACE.md)
is the language landing page.

Machine-readable discovery:

```bash
python3 workflow.py stages
python3 workflow.py software
python3 workflow.py parameters
python3 workflow.py paths --check
```

The repository includes the 222-row pic16/S1 core dataset under
`data/core/pic16_s1/` and its quantitative snapshot at
`reports/pic16_s1/metrics.json`. Large ORCA logs, GBW files, Excel workbooks,
and bulk figures are intentionally excluded.

All repository file paths in the documentation are relative to the repository
root. Site-specific ORCA, module, and checkout paths belong only in the
Git-ignored `config/cluster.env`.
