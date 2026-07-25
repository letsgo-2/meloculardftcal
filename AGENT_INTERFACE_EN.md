# meloculardftcal Agent Interface (English)

This document is the operational contract for autonomous agents that prepare,
run, inspect, analyze, or visualize this repository.

Human documentation:

- Chinese: [`README_CN.md`](README_CN.md)
- English: [`README_EN.md`](README_EN.md)
- Chinese agent interface: [`AGENT_INTERFACE_CN.md`](AGENT_INTERFACE_CN.md)
- Repository: <https://github.com/letsgo-2/meloculardftcal>

Machine-readable companions:

- Interface registry: `config/agent_interface.json`
- Read-only query CLI: `workflow.py`
- Environment report: `reports/environment.json`
- Core data/performance report: `reports/pic16_s1/metrics.json`

Agents should query the interface before constructing commands:

```bash
python3 workflow.py stages
python3 workflow.py stage S41_PARSE_TD
python3 workflow.py software
python3 workflow.py parameters
python3 workflow.py paths --check
```

All commands return JSON. `config/agent_interface.json` is declarative;
the referenced scripts remain authoritative for actual execution.
`python3 workflow.py parameters` includes calculation profiles, scheduler
controls, plotting defaults, parameter-edit locations, and output templates.

## 1. Scope and invariants

1. Use the repository root as the working directory.
2. Resolve every repository file from the relative paths in this document.
3. Do not assume that the repository directory has a fixed absolute path.
4. `ORCA_BIN`, `MODULE_INIT`, and `WORK_DIR` are external cluster locations and
   are the only paths that normally need to be absolute in
   `config/cluster.env`.
5. Do not commit ORCA binaries, MPI distributions, `.gbw`, large `.log`,
   Slurm output, private datasets, or generated figure collections.
6. Treat `data/raw/`, `data/processed/`, `work/`, and `results/figures/` as
   generated-data locations.
   `data/core/pic16_s1/` and small JSON files under `reports/` are intentionally
   versioned reference artifacts.
7. Preserve the source CSV row order after jobs have been generated. The TD
   parser maps the integer in `*_<id>td.log` to source CSV row `id - 1`.
8. Do not compare energies or rank candidates across different quantum
   protocols without recording the protocol difference.

## 2. Capability and environment contract

### 2.1 Reference versions

| Component | Required or recorded value |
|---|---|
| Python | `>=3.9` |
| ORCA parser target | ORCA `6.0.1` output |
| Recorded ORCA package | shared OpenMPI 4.1.6 AVX2 build |
| Recorded MPI | OpenMPI `4.1.6` |
| Recorded communication layer | UCX `1.15.0` |
| Scheduler | Slurm |
| Structure generation | Open Babel `obabel` CLI |
| Cheminformatics | RDKit |
| Data/plots | NumPy, pandas, SciPy, Matplotlib, Seaborn |
| Chemical-space plots | scikit-learn and UMAP |

Create the recommended local environment with:

```bash
conda env create -f environment.yml
conda activate dft-workflow
```

Minimum preflight:

```bash
python3 07_validation/check_environment.py \
  --output reports/environment.json
python3 workflow.py paths --check
```

On a compute cluster, also verify:

```bash
source config/cluster.env
"${ORCA_BIN}" --help
sbatch --version
squeue --version
```

An agent must report a failed capability check instead of silently switching
software, method, basis set, or solvent model.

### 2.2 Reference hardware and workload

The recorded node has 2 × AMD EPYC 7H12 CPUs, 128 physical cores, 256 logical
CPUs, and approximately 503 GiB RAM. GPU resources are not used. The historical
input collection contained about 2,998 molecules across multiple batches.

The default task profiles are:

| Profile | Nodes | MPI ranks | OMP threads/rank | ORCA `%maxcore` | Approximate declared memory |
|---|---:|---:|---:|---:|---:|
| `opt_standard` | 1 | 64 | 1 | 3000 MB/rank | 187.5 GiB |
| `opt_heavy_zora` | 1 | 64 | 1 | 3000 MB/rank | 187.5 GiB |
| `td_standard` | 1 | 32 | 1 | 3000 MB/rank | 93.75 GiB |
| `td_heavy_zora` | 1 | 32 | 1 | 3000 MB/rank | 93.75 GiB |

The memory estimate is `nprocs × maxcore`; actual peak usage may be higher.
Before increasing either value, verify physical memory, scheduler policy,
scratch capacity, molecule size, and ORCA parallel scaling.

## 3. Configuration interface

Create the private cluster configuration:

```bash
cp config/cluster.env.example config/cluster.env
```

`config/cluster.env` is ignored by Git. Its public interface is:

| Variable | Default/template | Consumer | Meaning |
|---|---|---|---|
| `SBATCH_PARTITION` | `amd_512q` | all four ORCA input generators | Slurm partition |
| `SBATCH_NODES` | `1` | all four ORCA input generators | nodes per job |
| `SBATCH_NTASKS` | `32` in template | all four generators | MPI ranks and generated ORCA `nprocs` |
| `WORK_DIR` | site-specific | generated `.sh` jobs | cluster checkout/work directory |
| `ORCA_BIN` | site-specific | generated `.sh` jobs | ORCA executable |
| `MODULE_INIT` | site-specific | generated `.sh` jobs | environment-modules initialization |
| `ORCA_MODULE` | ORCA 6.0.1 module name | generated `.sh` jobs | module to load |

Important: the ground-state generators default to 64 ranks and the TD
generators default to 32 ranks when `SBATCH_NTASKS` is unset. Sourcing the
example configuration sets one shared value. For stage-specific ranks, export
the intended value before each generation step:

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh work/pic16 work/pic16inp

export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh work/pic16 work/pic16inp
```

`OMP_NUM_THREADS=1` is written into generated job scripts and is not currently
exposed by `config/cluster.env`.

## 4. Stage registry

The following IDs are stable names for orchestration:

| Stage ID | Entry point | Primary inputs | Primary outputs | Side effects |
|---|---|---|---|---|
| `S00_PREFLIGHT` | commands in §2.1 | environment | capability report | read-only |
| `S10_STRUCTURE` | `01_structure_preparation/convertInputs.py` | source CSV | numbered `.gjf` files | writes selected output directory |
| `S20_OPT_GENERATE` | `01_structure_preparation/gjf_to_inp.sh` or `01_structure_preparation/gjf_to_inpp13.sh` | `.gjf` files | `.inp`, `.sh` | writes `work/<batch>inp/` |
| `S21_OPT_SUBMIT` | `02_cluster_execution/suboptsbatch.sh` | `.sh`, ID selection | Slurm jobs | external scheduler mutation |
| `S22_OPT_QC` | `02_cluster_execution/find_missing_opt_ids.sh`, logs, parser output | Opt logs | missing IDs and QC counts | read-only |
| `S30_TD_GENERATE` | `01_structure_preparation/gjf_to_inptd.sh` or `01_structure_preparation/gjf_to_inptdp13.sh` | `.gjf`, optimized `.xyz/.gbw` | `*td.inp`, `*td.sh` | writes job files |
| `S31_TD_SUBMIT` | `02_cluster_execution/subtdsbatch.sh` | TD job files, IDs | Slurm jobs | external scheduler mutation |
| `S32_TD_QC` | `02_cluster_execution/find_missing_td_ids.sh`, TD logs | TD logs | missing IDs and QC counts | read-only |
| `S40_PARSE_OPT` | `03_result_parsing/parse_orca_opt.py` | Opt logs | per-log CSV and summary | writes output directory |
| `S41_PARSE_TD` | `03_result_parsing/parse_orca_td_socv2.py` | TD logs, source CSV | per-log CSV and summary | writes output directory |
| `S42_MERGE_DFT` | `03_result_parsing/merge_dft_csv.py` | Opt/TD summaries | DFT v20-style CSV | writes output directory |
| `S50_MERGE_SOURCE` | `04_data_processing/merge_summary.py` | source CSV, DFT summary | v21-style indexed CSV | writes one CSV |
| `S51_ENRICH` | `04_data_processing/addrdkit.py` | v21 CSV | v22 CSV | writes one CSV |
| `S52_SYMMETRY` | `04_data_processing/addcenter.py` | CSV with `SMILES`, `Scaffold` | symmetry-enriched CSV | writes one CSV |
| `S60_METRICS` | `07_validation/collect_metrics.py` | source, TD summary, v20/v21/v22 | metrics JSON | writes `reports/` |
| `S70_ANALYZE` | scripts under `05_analysis/plotters/` | v22 CSV files | plots, selected CSVs | writes requested outputs |
| `S80_VISUALIZE` | scripts under `06_visualization/` | CSV/SMILES | molecular/scaffold images | writes requested outputs |

## 5. Canonical workflow calls

Replace `<batch>` with a stable batch token such as `pic16`. The source CSV
should be `data/raw/<batch>.csv`.

### 5.1 Structure generation

`01_structure_preparation/convertInputs.py` takes a mode, source CSV, and
optional output directory:

```bash
python3 01_structure_preparation/convertInputs.py \
  1 data/raw/<batch>.csv work/<batch>
```

Expected source columns include `SMILES`; `Scaffold` is needed by later
scaffold analysis. Output files are numbered in source-row order. Always pass
the output directory: omitting it preserves a historical repository-root
output mode.

### 5.2 Ground-state generation and submission

Standard elements:

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh \
  work/<batch> work/<batch>inp
bash 02_cluster_execution/suboptsbatch.sh \
  work/<batch>inp 1-20
```

Heavy-element/ZORA route:

```bash
bash 01_structure_preparation/gjf_to_inpp13.sh \
  work/<batch> work/<batch>inp
```

ID selectors may be individual values or inclusive ranges:

```bash
bash 02_cluster_execution/suboptsbatch.sh \
  work/<batch>inp 1-10 20 30-35
```

### 5.3 TD generation and submission

Proceed only for ground-state jobs that terminated normally and produced the
optimized geometry and GBW expected by the generator.

```bash
source config/cluster.env
export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh \
  work/<batch> work/<batch>inp
bash 02_cluster_execution/subtdsbatch.sh \
  work/<batch>inp 1-20
```

Use `01_structure_preparation/gjf_to_inptdp13.sh` for the
heavy-element/ZORA route.

### 5.4 Parsing and enrichment

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

`03_result_parsing/merge_dft_csv.py` derives expected filenames from the input
directory basename
and reads both summaries from its output directory. For
`work/<batch>inp`, these are `<batch>inpout_opt.csv` and
`<batch>inpout_tdv2.csv` under `work/<batch>inpout/`. The TD parser now writes
the `_tdv2.csv` summary directly.

## 6. Quantum-protocol interface

### 6.1 Current protocols

| Profile | ORCA keywords and blocks |
|---|---|
| standard Opt/Freq | `B3LYP def2-SVP Opt Freq RIJCOSX CPCM(Water)` |
| heavy Opt/Freq | `PBE0 D3BJ ZORA ZORA-def2-TZVP SARC/J Opt Freq RIJCOSX TightSCF CPCM(Water)` |
| standard TD/SOC | `wB97X-D3 def2-TZVP RIJCOSX TightSCF CPCM(Water)`; `%tddft NRoots 10, Triplets true, dosoc true`; `MORead`; optimized `xyzfile` |
| heavy TD/SOC | `wB97X-D3 ZORA ZORA-def2-TZVP SARC/J RIJCOSX RI-SOMF(1X) TightSCF CPCM(Water)` with the same TD controls |

TD jobs read the converged ground-state GBW and optimized XYZ. Charge and
multiplicity originate from the GJF charge/multiplicity line.

### 6.2 Where an agent may change calculation parameters

`config/agent_interface.json` exposes the current protocol as structured data,
but it does not rewrite ORCA inputs. Parameters must be edited in the generator
named below, followed by regeneration of `.inp` and `.sh` files.

| Parameter | Standard file | Heavy-element file | Search/edit target |
|---|---|---|---|
| Opt method/basis/solvent | `01_structure_preparation/gjf_to_inp.sh` | `01_structure_preparation/gjf_to_inpp13.sh` | emitted line beginning `!` |
| TD method/basis/solvent | `01_structure_preparation/gjf_to_inptd.sh` | `01_structure_preparation/gjf_to_inptdp13.sh` | emitted line beginning `!` |
| TD root count | `01_structure_preparation/gjf_to_inptd.sh` | `01_structure_preparation/gjf_to_inptdp13.sh` | `NRoots 10` |
| triplet states | `01_structure_preparation/gjf_to_inptd.sh` | `01_structure_preparation/gjf_to_inptdp13.sh` | `Triplets true` |
| SOC calculation | `01_structure_preparation/gjf_to_inptd.sh` | `01_structure_preparation/gjf_to_inptdp13.sh` | `dosoc true` |
| memory/rank | all four generators | all four generators | `%maxcore 3000`; generated `%pal`/`nprocs`; `SBATCH_NTASKS` |
| Slurm partition/node count | all four generators | all four generators | `SBATCH_PARTITION`, `SBATCH_NODES` |
| charge/multiplicity | source `.gjf` | source `.gjf` | charge/multiplicity line |

Mutation protocol:

1. Copy the current generator or create a Git branch.
2. Change one documented protocol dimension at a time when benchmarking.
3. Regenerate inputs in a new directory; never mix old and new `.inp` files.
4. Inspect at least one generated `.inp` and `.sh`.
5. Run a small reference subset before a full batch.
6. Record the exact ORCA/MPI versions, full keywords, resources, source commit,
   and input hashes in a run manifest.
7. Keep datasets from incompatible protocols separate during analysis.

## 7. Scheduler-control interface and safety

Read-only inspection:

```bash
bash 02_cluster_execution/find_missing_opt_ids.sh work/<batch>inp
bash 02_cluster_execution/find_missing_td_ids.sh work/<batch>inp
bash 02_cluster_execution/list_opt_jobs.sh <batch>
bash 02_cluster_execution/list_td_jobs.sh <batch>
```

Continuous submission:

```bash
bash 02_cluster_execution/auto_submit_opt.sh \
  work/<batch>inp 20 1800
bash 02_cluster_execution/auto_submit_td.sh \
  work/<batch>inp 20 1800
```

The two numeric arguments are `MAX_JOBS` and `SLEEP_SEC`. `MAX_JOBS` counts all
jobs belonging to the current Slurm user, not only the selected batch. These
commands are long-running and submit external jobs. An agent needs explicit
authority to invoke them.

The historical timeout cleanup script is excluded from this repository.
Cancellation must use an explicit reviewed job ID and separate authority for
`scancel`.

## 8. File and schema contract

### 8.1 Naming

The supported TD log pattern is `*_<id>td.log` with a positive integer ID.
Ground-state names are expected to preserve the same batch prefix and ID.
Avoid spaces in batch names.

Recommended version meanings:

| Name | Meaning |
|---|---|
| `<batch>inpdftv20.csv` | merged Opt and TD quantum results |
| `<batch>inpdftv21.csv` | source generation table plus DFT properties and stable index |
| `<batch>inpdftv22.csv` | v21 plus RDKit descriptors and SA Score |

### 8.2 Bundled data and inspectable result paths

The only complete research batch committed to Git is pic16, displayed as S1:

| Artifact | Repository-relative path |
|---|---|
| source molecules | `data/core/pic16_s1/input/pic16.csv` |
| TD/SOC summary | `data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv` |
| v20 | `data/core/pic16_s1/intermediate/pic16inpdftv20.csv` |
| v21 | `data/core/pic16_s1/intermediate/pic16inpdftv21.csv` |
| v22 | `data/core/pic16_s1/final/pic16inpdftv22.csv` |
| electronic-property reference table | `05_analysis/data/publishinpdftv22.csv` |
| checksums | `data/core/pic16_s1/SHA256SUMS` |
| dataset provenance | `data/core/pic16_s1/DATASET.md` |
| quantitative snapshot | `reports/pic16_s1/metrics.json` |
| current environment check | `reports/environment.json` |
| generated figures | `results/figures/pic16_s1/` |
| per-run provenance | `reports/runs/<run_id>.json` |

Use `python3 workflow.py paths --check` to discover these paths without parsing
Markdown. Raw ORCA logs and wavefunctions are not bundled, so normal
termination rates in the committed metric report are `null`.

### 8.3 Quantitative columns

| Column | Unit/type | Direction or role |
|---|---|---|
| `S1_energy_eV` | eV | target/window |
| `T1_energy_eV` | eV | context |
| `S1_wavelength_nm` | nm | target/window |
| `S1_fosc_D2` | dimensionless | usually maximize |
| `S1-T1` | eV | usually minimize absolute/positive gap |
| `min_positive_S1-Ti` | eV | minimize |
| `max_SOCME_S1_Ti_value` | cm^-1 | maximize |
| `max_SOCME_S1_Ti_index` | state label/index | provenance |
| `max_SOCME_over_deltaE_S1_Ti_value` | cm^-1/eV | maximize |
| `max_SOCME_over_deltaE_S1_Ti_index` | state label/index | provenance |
| `HOMO(eV)` / `LUMO(eV)` | eV | target/context |
| `Gap(eV)` | eV | target/context |
| `Dipole(D)` | Debye | target/context |
| `Gibbs(Eh)` | Hartree | compare only within one protocol |
| `Imaginary` | `Yes`/`No` | require `No` for a minimum |
| `charge` | integer | constraint/targeting proxy |
| `LogP` | dimensionless | feasibility window |
| `TPSA` | Å^2 | feasibility window |
| `Rotatable_Bonds` | count | flexibility constraint |
| `SA_Score` | approximately 1–10 | minimize |
| `Molecular_Weight` | g/mol | size constraint |

SOCME is parsed as the magnitude of the complex X/Y/Z components. The ratio is
`SOCME / abs(delta_E_eV)`. The summary maxima are restricted by the parser's
triplet-state eligibility rule relative to S1; inspect
`03_result_parsing/parse_orca_td_socv2.py` before changing that rule.

## 9. Quality gates and run metrics

An agent must calculate and report these counts for every batch:

```text
N_input            = source CSV row count
N_opt_normal       = ground-state logs containing "ORCA TERMINATED NORMALLY"
N_opt_parsed       = rows in the Opt summary
N_no_imaginary     = parsed rows where Imaginary == "No"
N_td_normal        = TD logs containing "ORCA TERMINATED NORMALLY"
N_td_parsed        = rows in the TD summary
N_complete_v22     = v22 rows with all required optimization metrics
```

Derived rates:

```text
opt_success_rate       = N_opt_normal / N_input
opt_parse_rate         = N_opt_parsed / N_opt_normal
minimum_rate           = N_no_imaginary / N_opt_parsed
td_success_rate        = N_td_normal / N_no_imaginary
td_parse_rate          = N_td_parsed / N_td_normal
pipeline_complete_rate = N_complete_v22 / N_input
```

Use `null` and an explanatory reason when a denominator is zero. Do not replace
missing values with zero.

Quality gate:

```text
source row and structure valid
  -> Opt log terminated normally
  -> optimized XYZ and GBW exist
  -> frequency result parsed and Imaginary == "No"
  -> TD log terminated normally
  -> S1, triplets, oscillator strength, and SOC sections parsed
  -> source/DFT join has one unambiguous row
  -> RDKit descriptors valid
  -> candidate is eligible for ranking
```

The Opt parser does not itself enforce normal ORCA termination. Check logs
before accepting its CSV. The merge is a union and can retain partial records;
row count alone is not a completeness test.

Generate or refresh the metric report:

```bash
python3 07_validation/collect_metrics.py
```

For a new batch:

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

Use `--no-reference` instead of `--reference` and `--reference-index` when no
scientifically compatible reference exists. Reference comparisons are valid
only when the candidate and reference quantum protocols are compatible.

The committed pic16/S1 snapshot contains 222 complete v22 rows, 191 minima
without imaginary frequencies (86.0%), 2 gap-versus-ratio Pareto rows, and 14
three-objective Pareto rows. Its default joint LogP/TPSA window has zero hits,
so an agent must retune that visualization window before using it as a
feasibility constraint.

Against row `pic16` in `05_analysis/data/publishinpdftv22.csv`, 175/222
candidates meet or improve both the S1–Ti gap and SOCME/gap ratio thresholds;
156/191 do so after the no-imaginary-frequency filter. Requiring SOCME itself
to also meet the reference gives 155/222 and 138/191, respectively. The
reference row has `Imaginary=Yes`, so these are electronic-property benchmark
counts—not comparisons against a validated minimum.

## 10. Optimization interface

A suggested agent-side objective specification is:

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

The displayed LogP/TPSA values are current plotting defaults, not scientifically
universal constraints. An agent must disclose all selected thresholds and
perform sensitivity checks before presenting a ranking as robust.

Recommended next-step logic:

1. Fix pipeline failures if success or parse rates are low.
2. Exclude non-minima and incomplete records from scientific ranking.
3. Compute a Pareto set instead of collapsing all objectives into an
   undocumented scalar.
4. Apply application-specific wavelength, charge, LogP, TPSA, MW, and SA
   constraints.
5. Select a structurally diverse subset with
   `05_analysis/plotters/plotsample.py`.
6. Recalculate top candidates with a higher-level or alternative protocol.
7. Compare rank stability and uncertainty before proposing new chemistry.

## 11. Analysis and visualization interface

### 11.1 Plot entry points

| Script | Quantities shown | Runtime controls |
|---|---|---|
| `05_analysis/plotters/plot_boxplot_HL.py` | HOMO, LUMO, gap, dipole and their relations | columns, output, width, height |
| `05_analysis/plotters/plot_boxplot_s1t1.py` | S1/T1 distributions and optional oscillator strength | columns, output, size, oscillator options |
| `05_analysis/plotters/plot_boxplot_s1_t1ti.py` | `S1-T1`, minimum positive `S1-Ti` | columns, output, size, imaginary filter |
| `05_analysis/plotters/plot_boxplot_socme.py` | SOCME and `SOCME/abs(deltaE)` distributions | columns, output, size |
| `05_analysis/plotters/plot_pareto_s1ti_socme.py` | gap versus SOCME Pareto front | output, size, axis limits, Pareto toggle |
| `05_analysis/plotters/plot_pareto_s1ti_socmegap.py` | gap versus SOCME/gap Pareto front | output, size, axis limits, Pareto toggle |
| `05_analysis/plotters/plot_ridge_socmev2.py` | paired SOCME and ratio ridges | columns, output, size, overlap, scale, labels, x maxima |
| `05_analysis/plotters/plot_summary.py` | SA, MW, LogP, S1, wavelength KDE/joint views | input SMILES column and output; code-level limits |
| `05_analysis/plotters/plot_summary2.py` | `S1-T1` vs `T1/S1`, colored by S1 | inputs; code-level axes/colors |
| `05_analysis/plotters/plot_targeting.py` | charge, positive fraction, LogP-TPSA and charge-TPSA | output directory, filters, windows, PDF/facets |
| `05_analysis/plotters/plot_decorations_statistics.py` | decoration-fragment counts and structures | scaffold column, top N, output size |
| `05_analysis/plotters/plotsample.py` | fingerprint space, MaxMin/random samples, NND and source sampling rates | pick count, column, sample CSV, plot directory, labels, method, pool, fingerprint, seed, projection settings |

Discover exact CLI controls without editing source:

```bash
python3 05_analysis/plotters/<script>.py --help
```

Scripts with code-level-only plot limits must be edited in
`05_analysis/plotters/<script>.py`; record the changed constant and do not
interpret an axis crop as a data filter. Prefer a new CLI option when making
reusable changes.

### 11.2 Standard commands

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

### 11.3 Molecular images

| Script | Purpose | Main adjustable settings |
|---|---|---|
| `06_visualization/generate_molecule_images.py` | plain 2D molecule grids/files | CSV, SMILES column, output location, image size in source |
| `06_visualization/generate_molecule_images_with_props.py` | structures annotated with selected properties | property columns, labels, output/image settings |
| `06_visualization/generate_top_selected_images.py` | top candidates relative to references | ranking columns, direction, top N, output |
| `06_visualization/execlwithimage.py` | Excel workbook with embedded images | source CSV, image folder, output workbook |
| `06_visualization/scaffold/generate_mol_and_scaffold_figsv2.py` | molecule, core, and decorations | CSV row/columns, output, width, height |
| `06_visualization/scaffold/symmetry_and_drawv2.py` | attachment-point symmetry and drawings | chirality, match/movement limits, output prefix |

For scripts without an argparse interface, edit their `if __name__ ==
"__main__"` configuration block or add a CLI before automation. Do not assume
that a hard-coded historical path is a supported parameter.

## 12. Provenance manifest

For each calculation campaign, an agent should create a non-sensitive manifest,
for example `reports/runs/<run_id>.json`. The following YAML-shaped example
shows the required fields; JSON is preferred for the saved artifact:

```yaml
run_id: "<batch>-<date>-<protocol-tag>"
repository_commit: "<git-commit>"
source_csv: "data/raw/<batch>.csv"
source_sha256: "<sha256>"
input_count: null
software:
  orca: "6.0.1"
  mpi: "4.1.6"
  python: null
machine:
  cpu_model: null
  physical_cores: null
  memory_gib: null
protocol:
  generator: "01_structure_preparation/gjf_to_inptd.sh"
  keywords: null
  nroots: 10
  triplets: true
  dosoc: true
resources:
  partition: null
  nodes: 1
  mpi_ranks: 32
  omp_threads: 1
  maxcore_mb_per_rank: 3000
quality_counts:
  N_opt_normal: null
  N_no_imaginary: null
  N_td_normal: null
  N_complete_v22: null
outputs:
  dataset: "data/processed/<batch>inpdftv22.csv"
  figures: "results/figures/"
plot_parameters: {}
```

Do not store usernames, credentials, private hostnames, or licensed-software
locations in a committed manifest.

## 13. Known automation gaps

1. Quantum keywords are embedded in four shell generators rather than one
   structured config.
2. Some plotting limits, labels, and image sizes are source constants rather
   than CLI parameters.
3. `01_structure_preparation/convertInputs.py` uses positional arguments and
   has no argparse-style `--help` interface.
4. The TD-to-SMILES join depends on source row order and filename IDs.
5. ORCA parsers target the recorded ORCA 6.0.1 text layout and need regression
   tests before supporting another major version.
6. Opt parsing does not itself require a normal-termination marker.
7. Merge steps can preserve partial records.
8. `config/agent_interface.json` is a read-only registry, not a workflow engine
   and not an automatic parameter writer.
9. Most scientific plotting scripts report through stdout and files rather
   than stable JSON; `workflow.py` and `collect_metrics.py` provide JSON at the
   orchestration boundary.

When improving automation, prioritize a structured protocol configuration,
non-interactive structure generation, explicit join keys, parser fixtures,
and CLI exposure of remaining plot constants.

## 14. Agent completion report

Every autonomous run should return:

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

Never label a run `success` only because a process exited with code zero. Apply
the stage-specific quality gates and verify the expected output files.
