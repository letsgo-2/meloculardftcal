# meloculardftcal: High-throughput Molecular DFT / TD-DFT / SOC Workflow

[中文 README](README_CN.md) | [English agent interface](AGENT_INTERFACE_EN.md) | [中文智能体接口](AGENT_INTERFACE_CN.md) | [Third-party notices](THIRD_PARTY_NOTICES.md)

Repository: <https://github.com/letsgo-2/meloculardftcal>

This repository provides a file-based, high-throughput molecular quantum-chemistry workflow centered on ORCA. It covers 3D structure preparation, ground-state geometry optimization and frequency analysis, TD-DFT excited states and spin-orbit coupling (SOC), Slurm scheduling, result parsing, property enrichment, statistical analysis, and molecular/scaffold visualization.

The repository contains source code, configuration templates, documentation,
small examples, and one pic16/S1 core tabular dataset. ORCA distributions,
OpenMPI source trees, raw quantum-chemistry logs, GBW wavefunctions, other
complete batches, and bulk figures should not be committed.

## 1. Workflow and scientific objectives

```text
Generated-molecule CSV / SMILES
        │
        ▼
Open Babel 3D structure generation
        │
        ▼
ORCA ground-state Opt + Freq
        │
        ├── optimized .xyz geometry
        ├── .gbw molecular orbitals
        └── ground-state, thermal, and frequency properties
        │
        ▼
TD-DFT + triplets + SOC
        │
        ├── S1/T1 energies
        ├── S1 oscillator strength and wavelength
        ├── S1−Ti energy gaps
        └── SOCME and SOCME/|ΔE|
        │
        ▼
generation metrics + DFT + RDKit descriptors
        │
        ▼
distribution comparison, Pareto selection,
chemical-space analysis, and structure visualization
```

The current multi-objective search favors candidates with:

- a small positive `S1−Ti` energy gap;
- a large `SOCME(S1,Ti)`;
- a large `SOCME/|ΔE|`;
- application-appropriate S1 energy, absorption wavelength, and oscillator strength;
- an optimized structure without imaginary frequencies;
- acceptable molecular weight, LogP, TPSA, rotatable-bond count, and SA Score;
- sufficient structural diversity across the candidate set.

These are optimization directions, not universal fixed thresholds. Charge rules, LogP/TPSA windows, wavelength targets, and Pareto plotting ranges must be chosen for the intended application.

## 2. Repository layout

```text
.
├── README.md                         # language landing page
├── README_CN.md                      # complete Chinese documentation
├── README_EN.md                      # this document
├── AGENT_INTERFACE.md                # agent-interface language landing page
├── AGENT_INTERFACE_CN.md             # complete Chinese agent interface
├── AGENT_INTERFACE_EN.md             # complete English agent interface
├── workflow.py                       # read-only machine interface CLI
├── THIRD_PARTY_NOTICES.md
├── environment.yml
├── requirements.txt
├── config/
│   ├── agent_interface.json          # stage/path/version/parameter registry
│   ├── cluster.env.example           # cluster resource and ORCA template
│   └── pictos.txt                    # batch-to-display-label mapping
├── 01_structure_preparation/
│   ├── convertInputs.py              # SMILES to 3D GJF
│   ├── gjf_to_inp.sh                 # standard Opt/Freq
│   ├── gjf_to_inpp13.sh              # heavy-element/ZORA Opt/Freq
│   ├── gjf_to_inptd.sh               # standard TD-DFT/SOC
│   └── gjf_to_inptdp13.sh            # heavy-element/ZORA TD-DFT/SOC
├── 02_cluster_execution/
│   ├── suboptsbatch.sh
│   ├── subtdsbatch.sh
│   ├── find_missing_opt_ids.sh
│   ├── find_missing_td_ids.sh
│   ├── list_opt_jobs.sh
│   ├── list_td_jobs.sh
│   ├── auto_submit_opt.sh
│   └── auto_submit_td.sh
├── 03_result_parsing/
│   ├── parse_orca_opt.py
│   ├── parse_orca_td_socv2.py
│   └── merge_dft_csv.py
├── 04_data_processing/
│   ├── merge_summary.py
│   ├── addrdkit.py
│   ├── addcenter.py
│   └── calculate_molecular_weight.py
├── 05_analysis/
│   ├── plotters/
│   └── data/
├── 06_visualization/
│   └── scaffold/
├── 07_validation/
│   ├── check_environment.py
│   ├── collect_metrics.py
│   └── UnitTestSAScore.py
├── data/
│   ├── core/pic16_s1/                # bundled 222-row S1 core dataset
│   ├── raw/
│   ├── processed/
│   └── example_molecules.csv
├── examples/example_v22.csv
├── reports/
├── work/
└── results/figures/
```

Run commands from the repository root. All file paths in this documentation are relative to that root.

### 2.1 Machine-readable interface

English-speaking agents should read
[`AGENT_INTERFACE_EN.md`](AGENT_INTERFACE_EN.md); Chinese-speaking agents
should read [`AGENT_INTERFACE_CN.md`](AGENT_INTERFACE_CN.md), then use
`config/agent_interface.json` for stable stage IDs, entry points, versions,
path templates, and defaults. Read-only queries are:

```bash
python3 workflow.py stages
python3 workflow.py stage S41_PARSE_TD
python3 workflow.py software
python3 workflow.py parameters
python3 workflow.py paths --check
```

Every command emits JSON. `config/agent_interface.json` is the interface
registry; the four `01_structure_preparation/gjf_to_*.sh` generators remain
the execution authority for ORCA keywords. Regenerate inputs after changing
them.

## 3. Reference execution environment

### 3.1 Recorded cluster environment

The scripts were originally used in the following environment. This is a reproducibility reference, not a strict requirement for another cluster.

| Component | Recorded configuration |
|---|---|
| Operating system | CentOS Linux 7, Linux 3.10 |
| Scheduler | Slurm |
| ORCA | 6.0.1 shared OpenMPI 4.1.6 AVX2 distribution |
| MPI | OpenMPI 4.1.6 |
| Communication layer | UCX 1.15.0 |
| CPU | 2 × AMD EPYC 7H12 64-Core |
| Physical cores | 128 per node |
| Logical CPUs | 256 per node, SMT2 |
| Memory | approximately 503 GiB per node |
| Swap | none |
| GPU | not required |
| Storage | shared parallel cluster filesystem |

The cluster's system Python was historically Python 2.7. This repository requires Python 3.9 or newer for data processing and plotting. Use the supplied Conda environment instead of the cluster system Python.

### 3.2 Historical workload scale

The pre-cleanup workspace contained approximately:

| Batch | Molecules |
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
| publish reference set | 9 |
| Total | 2,998 |

These are input counts, not successful-completion counts. Each molecule normally has ground-state and TD input/job scripts and may produce logs, XYZ files, GBW files, and scheduler output.

### 3.3 Per-job resource settings

Resource variables are defined through [`config/cluster.env.example`](config/cluster.env.example).

| Parameter | Opt/Freq default | TD/SOC default | Meaning |
|---|---:|---:|---|
| `SBATCH_PARTITION` | `amd_512q` | `amd_512q` | Slurm partition |
| `SBATCH_NODES` | 1 | 1 | nodes per job |
| `SBATCH_NTASKS` | 64 | 32 | MPI ranks and ORCA `%pal nprocs` |
| `OMP_NUM_THREADS` | 1 | 1 | OpenMP threads per MPI rank |
| `%maxcore` | 3000 MB | 3000 MB | approximate ORCA memory limit per process |

A simple `maxcore × nprocs` estimate gives:

- 64 ranks: about 192,000 MB, or 187.5 GiB;
- 32 ranks: about 96,000 MB, or 93.75 GiB.

Actual peak memory can exceed this estimate. Any change to rank count or `maxcore` must be checked against node memory, scheduler policy, molecule size, and ORCA parallel scaling.

Automatic submission defaults:

| Parameter | Default | Meaning |
|---|---:|---|
| `MAX_JOBS` | 20 | maximum total Slurm jobs for the current user, not just this batch |
| `SLEEP_SEC` | 1800 s | polling interval |

The historical timeout cleanup script is intentionally excluded from the
public, streamlined repository because it could cancel jobs and delete job
scripts. Inspect `squeue` and use an explicit `scancel <job_id>` only with
separate scheduler-mutation authorization.

## 4. Software setup

### 4.1 Conda

```bash
conda env create -f environment.yml
conda activate dft-workflow
```

The environment includes NumPy, pandas, SciPy, Matplotlib, Seaborn, Pillow, RDKit, Open Babel, scikit-learn, UMAP, openpyxl, and tqdm.

### 4.2 pip

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Structure generation and some image tools use the `obabel` executable. The pip
dependencies do not install this CLI, so installing Open Babel from conda-forge
is recommended.

### 4.3 ORCA and cluster configuration

```bash
cp config/cluster.env.example config/cluster.env
```

Edit the local file:

```bash
export SBATCH_PARTITION="amd_512q"
export SBATCH_NODES="1"
export SBATCH_NTASKS="32"
export WORK_DIR="<ABSOLUTE_REPOSITORY_PATH_ON_CLUSTER>"
export ORCA_BIN="<ABSOLUTE_ORCA_EXECUTABLE_PATH_ON_CLUSTER>"
export MODULE_INIT="<ABSOLUTE_MODULE_INITIALIZATION_PATH_ON_CLUSTER>"
export ORCA_MODULE="orca/6_0_1_linux_x86-64_shared_openmpi416_avx2"
```

Load it before generating jobs:

```bash
source config/cluster.env
```

`config/cluster.env` is ignored by Git because it contains installation-specific paths.
The template's `SBATCH_NTASKS=32` is a shared initial value. Under the recorded
resource profile, explicitly set it to 64 before ground-state generation and
back to 32 before TD generation. Regenerate `.inp/.sh` after changing it so
Slurm ranks and ORCA `nprocs` remain consistent.

## 5. Quantum-chemistry protocols

### 5.1 Standard ground-state Opt/Freq

Generator:

```text
01_structure_preparation/gjf_to_inp.sh
```

Default input:

```text
! B3LYP def2-SVP Opt Freq RIJCOSX CPCM(Water)
%maxcore 3000
%pal nprocs <SBATCH_NTASKS> end
```

This uses B3LYP/def2-SVP for geometry optimization and harmonic frequencies, RIJCOSX acceleration, and CPCM water solvation.

### 5.2 Heavy-element/ZORA ground-state Opt/Freq

Generator:

```text
01_structure_preparation/gjf_to_inpp13.sh
```

Default input:

```text
! PBE0 D3BJ ZORA ZORA-def2-TZVP SARC/J Opt Freq RIJCOSX TightSCF CPCM(Water)
```

This protocol adds PBE0, D3BJ dispersion, scalar-relativistic ZORA, ZORA-def2-TZVP, SARC/J, and TightSCF. It originated from a specific heavy-element batch and must be reassessed for new elements or charge states.

### 5.3 Standard TD-DFT/SOC

Generator:

```text
01_structure_preparation/gjf_to_inptd.sh
```

Default input:

```text
! wB97X-D3 def2-TZVP RIJCOSX TightSCF CPCM(Water)
%tddft
    NRoots   10
    Triplets true
    dosoc true
end
```

The TD job reads the optimized geometry and wavefunction:

```text
%scf
  Guess MORead
end
%moinp "<base>.gbw"
* xyzfile <charge> <multiplicity> <base>.xyz
```

Do not submit TD calculations until the corresponding ground-state `.xyz` and `.gbw` files are complete and validated.

### 5.4 Heavy-element/ZORA TD-DFT/SOC

Generator:

```text
01_structure_preparation/gjf_to_inptdp13.sh
```

Default input:

```text
! wB97X-D3 ZORA ZORA-def2-TZVP SARC/J RIJCOSX RI-SOMF(1X) TightSCF CPCM(Water)
```

The calculation still requests ten roots, triplets, and SOC, while adding ZORA and `RI-SOMF(1X)`.

### 5.5 Changing calculation parameters

| Parameter | Where to change it | Required follow-up |
|---|---|---|
| functional, basis set, dispersion, solvent, SCF keywords | `echo "! ..."` in `01_structure_preparation/gjf_to_*.sh` | regenerate `.inp/.sh`; keep old results separate |
| TD root count | `NRoots` in both TD generator scripts | regenerate TD inputs and recompute |
| triplet/SOC calculation | `Triplets`, `dosoc` | regenerate and recompute |
| memory per process | `%maxcore` in the four generators | recheck node memory and rank count |
| MPI ranks | `SBATCH_NTASKS` in `config/cluster.env` | reload the environment and regenerate jobs |
| partition/nodes | `SBATCH_PARTITION`, `SBATCH_NODES` | regenerate job scripts |
| charge/multiplicity | charge/multiplicity line in the source GJF | regenerate all ORCA inputs |
| solvent | `CPCM(Water)` | treat as a new protocol; do not merge with old results silently |

For every protocol change, record ORCA version, MPI version, full keywords, resource settings, date, and Git commit.

## 6. Data and naming contract

The minimum source CSV is:

```csv
MOLNAME,SMILES
pic16_0001,CCO
pic16_0002,c1ccccc1
```

Recommended task names:

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

Important invariants:

1. `<batch>_<four-digit ID>` must be stable across all stages.
2. The TD parser maps a log ID back to the source CSV row.
3. `0001` maps to the first data row and `0012` maps to the twelfth data row.
4. Reordering the source CSV after job generation can mismatch SMILES and quantum results.
5. `03_result_parsing/merge_dft_csv.py` joins on `filename`.
6. `04_data_processing/merge_summary.py` converts `filename` to `ID` before merging with generation data.

Optional generation columns include `Scaffold`, `st_gap`, `charge_score`, `symmetry_score`, `decoration_ring_penalty`, `total_score`, `Step`, and `Molecular_Weight`.

### 6.1 Bundled pic16 / S1 core data

`config/pictos.txt` maps batch `pic16` to display label `S1`. The public
repository includes only this complete 222-row research batch:

| Path | Contents |
|---|---|
| `data/core/pic16_s1/input/pic16.csv` | generation metrics, scaffolds, and SMILES |
| `data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv` | TD-DFT/SOC summary |
| `data/core/pic16_s1/intermediate/pic16inpdftv20.csv` | merged ground-state and TD/SOC data |
| `data/core/pic16_s1/intermediate/pic16inpdftv21.csv` | source-generation fields restored |
| `data/core/pic16_s1/final/pic16inpdftv22.csv` | final RDKit/SA Score dataset |
| `05_analysis/data/publishinpdftv22.csv` | small reference-molecule table for relative performance |
| `reports/pic16_s1/metrics.json` | quality and scientific-performance snapshot |

See [`data/core/pic16_s1/DATASET.md`](data/core/pic16_s1/DATASET.md) for
lineage, row counts, and checksums. Raw ORCA logs, GBW/XYZ files, Excel
workbooks, and bulk images are excluded.

## 7. End-to-end execution

Run every command below from the repository root.

### 7.1 SMILES to 3D GJF

```bash
python3 01_structure_preparation/convertInputs.py \
  1 data/raw/pic16.csv work/pic16
```

The positional arguments are mode, source CSV, and an optional output
directory. Omitting the output directory retains the historical behavior of
writing a same-named directory at the repository root. Explicit
`work/<batch>` output is recommended.

`01_structure_preparation/convertInputs.py` modes:

| Mode | Function |
|---:|---|
| 1 | SMILES to 3D GJF for optimization preparation |
| 2 | legacy Gaussian LOG to TD GJF |
| 3 | legacy Gaussian Guess=Read TD input |

The current ORCA path mainly uses mode 1. Direct Open Babel usage is also possible:

```bash
mkdir -p work/pic16
obabel -:"CCO" -O work/pic16/pic16_0001.gjf --gen3d
```

### 7.2 Generate and submit ground-state jobs

```bash
source config/cluster.env
export SBATCH_NTASKS=64
bash 01_structure_preparation/gjf_to_inp.sh work/pic16 work/pic16inp
bash 02_cluster_execution/suboptsbatch.sh work/pic16inp 1-20
```

Multiple ranges are accepted:

```bash
bash 02_cluster_execution/suboptsbatch.sh work/pic16inp 1-10 20 30-35
```

Find unfinished IDs:

```bash
bash 02_cluster_execution/find_missing_opt_ids.sh work/pic16inp
```

Run the submission controller:

```bash
bash 02_cluster_execution/auto_submit_opt.sh work/pic16inp 20 1800
```

The `pic16` string in utility filenames is historical; the scripts infer the actual prefix from the supplied directory.

### 7.3 Ground-state quality gate

Before TD submission, verify:

- `ORCA TERMINATED NORMALLY` is present;
- geometry optimization completed;
- `.xyz` and `.gbw` exist and are non-empty;
- `Imaginary == No`, unless a transition state was intentionally requested;
- no dissociation, implausible bond lengths, or incorrect protonation;
- charge and multiplicity are correct.

### 7.4 Generate and submit TD/SOC jobs

```bash
export SBATCH_NTASKS=32
bash 01_structure_preparation/gjf_to_inptd.sh work/pic16 work/pic16inp
bash 02_cluster_execution/subtdsbatch.sh work/pic16inp 1-20
```

Find missing TD jobs:

```bash
bash 02_cluster_execution/find_missing_td_ids.sh work/pic16inp
```

Run the TD submission controller:

```bash
bash 02_cluster_execution/auto_submit_td.sh work/pic16inp 20 1800
```

The TD pending set is:

```text
normally completed ground-state IDs − normally completed TD IDs
```

### 7.5 Parse ground-state results

```bash
python3 03_result_parsing/parse_orca_opt.py \
  work/pic16inp \
  -o work/pic16inpout
```

The parser extracts final energy, nuclear repulsion, HOMO, LUMO, orbital gap, dipole, Gibbs free energy, and an imaginary-frequency flag.

The current parser does not independently require `ORCA TERMINATED NORMALLY`; it primarily requires a parseable HOMO block. An automated controller must run the normal-termination check before parsing.

### 7.6 Parse TD/SOC results

```bash
python3 03_result_parsing/parse_orca_td_socv2.py \
  work/pic16inp \
  -o work/pic16inpout \
  --csv data/raw/pic16.csv
```

Batch parsing directly creates
`work/pic16inpout/pic16inpout_tdv2.csv`, matching the merge-stage contract.

Definitions:

- `S1-T1 = E(S1) − E(T1)` in eV;
- `min_positive_S1-Ti` is the smallest positive `E(S1) − E(Ti)`;
- SOCME is the combined magnitude of ORCA's complex X/Y/Z components, in cm⁻¹;
- `SOCME/|ΔE|` is reported in cm⁻¹/eV;
- the maximum SOC metrics are summarized over triplet states up to the closest state lying below S1.

### 7.7 Build v20, v21, and v22

Merge ground-state and TD/SOC properties:

```bash
python3 03_result_parsing/merge_dft_csv.py \
  work/pic16inp \
  -o work/pic16inpout
```

Output:

```text
work/pic16inpout/pic16inpdftv20.csv
```

Merge generation data:

```bash
python3 04_data_processing/merge_summary.py \
  --piccsvfile data/raw/pic16.csv \
  --summary work/pic16inpout/pic16inpdftv20.csv \
  --output data/processed/pic16inpdftv21.csv
```

Add RDKit descriptors:

```bash
python3 04_data_processing/addrdkit.py \
  data/processed/pic16inpdftv21.csv \
  --output data/processed/pic16inpdftv22.csv
```

| Version | Contents |
|---|---|
| v20 | ground-state + TD/SOC properties |
| v21 | v20 + generation/scaffold metrics |
| v22 | v21 + LogP, TPSA, rotatable bonds, formal charge, SA Score |
| v23 | image-embedded Excel representation of v22 |

## 8. Result fields and optimization meaning

| Field | Unit | Typical use/direction |
|---|---|---|
| `S1_energy_eV` | eV | target excited-state energy |
| `T1_energy_eV` | eV | lowest triplet energy |
| `S1_wavelength_nm` | nm | absorption wavelength target |
| `S1_fosc_D2` | dimensionless | larger often means a stronger allowed transition |
| `S1-T1` | eV | singlet-triplet gap |
| `min_positive_S1-Ti` | eV | smaller favors near-resonance |
| `max_SOCME_S1_Ti_value` | cm⁻¹ | larger indicates stronger SOC |
| `max_SOCME_over_deltaE_S1_Ti_value` | cm⁻¹/eV | combined coupling/gap proxy; larger is usually favored |
| `HOMO(eV)`, `LUMO(eV)` | eV | frontier orbital positions |
| `Gap(eV)` | eV | HOMO-LUMO gap |
| `Dipole(D)` | Debye | polarity proxy |
| `Gibbs(Eh)` | Hartree | compare only under an identical protocol |
| `Imaginary` | Yes/No | stationary-point quality check |
| `charge` | integer | formal charge |
| `LogP` | dimensionless | lipophilicity |
| `TPSA` | Å² | polar surface area |
| `Rotatable_Bonds` | count | flexibility |
| `SA_Score` | approximately 1–10 | lower generally indicates easier synthesis |
| `Molecular_Weight` | g/mol | molecular size |

Absolute energies from different functionals, basis sets, solvent models, or charge states must not be pooled without a deliberate calibration.

## 9. Analysis and plotting

Invoke the Python programs under `05_analysis/plotters/` directly and pass
explicit input/output paths. A fresh clone can use
`data/core/pic16_s1/final/pic16inpdftv22.csv`; new runs normally use
`data/processed/<batch>inpdftv22.csv`. Historical fixed-whitelist wrappers
have been removed.

### 9.1 Plot catalog and adjustable settings

| Script | What it displays | Main CLI controls |
|---|---|---|
| `05_analysis/plotters/plot_boxplot_HL.py` | HOMO, LUMO, gap, and dipole distributions; HOMO-LUMO, gap-dipole, gap-S1T1, and charge-dipole relations | `--columns`, `--output`, `--width`, `--height` |
| `05_analysis/plotters/plot_boxplot_s1t1.py` | S1/T1 energy distributions and optional S1 oscillator strength | `--columns`, `--output`, size, `--enable-osc`, `--osc-column` |
| `05_analysis/plotters/plot_boxplot_s1_t1ti.py` | `S1-T1` and minimum positive `S1-Ti` distributions | columns, output, size, `--filter-imaginary` |
| `05_analysis/plotters/plot_boxplot_socme.py` | maximum SOCME and maximum `SOCME/|ΔE|` distributions | columns, output, size |
| `05_analysis/plotters/plot_pareto_s1ti_socme.py` | X=min positive S1−Ti, Y=SOCME, Pareto front | output, size, x/y limits, `--show-pareto` |
| `05_analysis/plotters/plot_pareto_s1ti_socmegap.py` | X=min positive S1−Ti, Y=`SOCME/|ΔE|`; combined and per-batch plots | same controls as above |
| `05_analysis/plotters/plot_ridge_socmev2.py` | two-panel SOCME and `SOCME/|ΔE|` ridges with full-data median/max | columns, output, size, overlap, scale, x labels/maxima |
| `05_analysis/plotters/plot_summary.py` | KDE and joint plots for SA, MW, LogP, S1, and wavelength | SMILES column and first output name; some limits are code constants |
| `05_analysis/plotters/plot_summary2.py` | X=`S1-T1`, Y=`T1/S1`, color=S1, with reference/generation overlays | input files only; axes and colors are code constants |
| `05_analysis/plotters/plot_targeting.py` | charge composition, positive fraction, LogP-TPSA, charge-TPSA, facets, and a multi-panel figure | output dir, labels, imaginary filter, charge settings, LogP/TPSA window, PDF/facet toggles |
| `05_analysis/plotters/plot_decorations_statistics.py` | decoration-fragment frequencies and fragment depictions | scaffold column, output, top N, width, height |
| `05_analysis/plotters/plotsample.py` | Morgan-fingerprint space; MaxMin vs random sampling; UMAP/t-SNE/MDS; NND and source sampling rates | sample count, SMILES column, sample CSV, plot directory, labels, projection, pool, fingerprint, seed, and projection settings |

### 9.2 Frontier orbitals and dipole

```bash
python3 05_analysis/plotters/plot_boxplot_HL.py \
  data/processed/pic3inpdftv22.csv \
  data/processed/pic16inpdftv22.csv \
  --output results/figures/frontier \
  --width 14 \
  --height 8
```

This produces `_violin`, `_box`, `_scatter`, `_gap_dipole`, `_gap_s1t1`, and `_charge_dipole` PNG files.

### 9.3 S1/T1 and oscillator strength

```bash
python3 05_analysis/plotters/plot_boxplot_s1t1.py \
  data/processed/pic16inpdftv22.csv \
  --enable-osc \
  --osc-column S1_fosc_D2 \
  --output results/figures/s1_t1.png
```

### 9.4 Energy gaps with frequency filtering

```bash
python3 05_analysis/plotters/plot_boxplot_s1_t1ti.py \
  data/processed/pic16inpdftv22.csv \
  --filter-imaginary \
  --output results/figures/st_gap.png
```

For formal statistics, enabling `--filter-imaginary` is recommended.

### 9.5 Pareto analysis

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

Pareto directions:

```text
minimize X: min_positive_S1-Ti
maximize Y: SOCME or SOCME/|ΔE|
```

Plot ranges are display settings and should not silently filter the source dataset.

### 9.6 Targeting proxy analysis

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

The current `LogP=1–5`, `TPSA=20–120 Å²`, and positive-charge settings are adjustable proxy windows, not universal biological rules.

### 9.7 Chemical space and sampling

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

Defaults and their CLI controls:

- Morgan radius 2 and 2048 bits: `--radius`, `--fp-size`;
- random seed 42: `--seed`;
- UMAP uses Jaccard with 30 neighbors and 0.1 minimum distance:
  `--umap-neighbors`, `--umap-min-dist`;
- t-SNE uses PCA up to 50 dimensions and perplexity 30:
  `--tsne-perplexity`;
- MDS: pairwise Jaccard distances and high memory/time cost for large sets.

The selected MaxMin rows are written to `--output`; the assessment figure and
projected coordinates are written under `--plot-dir`. `umap-learn` is imported
only for `--method umap`, so help, t-SNE, and MDS remain available without it.

## 10. Molecular and scaffold visualization

### 10.1 Plain molecular images

```bash
python3 06_visualization/generate_molecule_images.py \
  data/processed/pic16inpdftv22.csv \
  --tool rdkit \
  --output results/figures/molecules
```

Controls: `--tool rdkit|openbabel`, `--output`, and `--test INDEX`.

### 10.2 Molecular images with properties

```bash
python3 06_visualization/generate_molecule_images_with_props.py \
  data/processed/pic16inpdftv22.csv \
  --tool rdkit \
  --output results/figures/properties \
  --verbose
```

Available labels include the HOMO-LUMO gap, S1, T1, `S1-T1`, minimum `S1-Ti`, SOCME, `SOCME/δE`, ID, charge, S1 wavelength, and S1 oscillator strength.

Controls: `--tool`, `--output`, `--test`, and `--verbose`.

### 10.3 Top-N candidates relative to references

```bash
python3 06_visualization/generate_top_selected_images.py \
  data/processed/pic16inpdftv22.csv \
  --reference-csv 05_analysis/data/publishinpdftv22.csv \
  --top-n 20 \
  --tool rdkit \
  --output-dir results/figures/top_selected
```

The script requires a smaller candidate `min(S1-Ti)` and a larger candidate `SOCME/δE` than the corresponding reference molecule before ranking the top N.

### 10.4 Molecule, scaffold, and decorations

```bash
python3 06_visualization/scaffold/generate_mol_and_scaffold_figs.py \
  --csv data/processed/pic16inpdftv22.csv \
  --row 1 \
  --output-dir results/figures/scaffold \
  --w 900 \
  --h 700 \
  --show-properties
```

Outputs include the complete molecule, scaffold, and individual attachment-point decorations. Controls include direct SMILES or CSV mode, scaffold text, column names, row number, output name, dimensions, and property labels.

`06_visualization/scaffold/generate_mol_and_scaffold_figsv2.py` uses more permissive scaffold matching and was originally intended for pic4 cases that fail the standard matcher. It should not be selected automatically without validating the match.

### 10.5 Image-embedded Excel output

```bash
python3 06_visualization/execlwithimage.py \
  --csv data/processed/pic16inpdftv22.csv \
  --output results/pic16inpdftv23.xlsx \
  --smiles-col SMILES \
  --img-col-name Structure \
  --img-width 200 \
  --img-height 200
```

The script uses Open Babel images. `--img-height` primarily controls the Excel row height; Open Babel keeps the actual molecular image aspect ratio.

## 11. Quality and quantitative run metrics

Track at least:

| Metric | Definition |
|---|---|
| input count | source CSV rows |
| ground-state normal-termination rate | normally terminated ground logs / inputs |
| ground-state parse coverage | rows in opt summary / normally terminated ground logs |
| no-imaginary-frequency rate | `Imaginary == No` / parsed ground states |
| TD normal-termination rate | normally terminated TD logs / successful ground states |
| TD parse coverage | TD summary rows / normally terminated TD logs |
| complete-v22 rate | rows with all core fields / inputs |
| Pareto count and fraction | nondominated candidates / complete candidates |
| target-window hit rate | candidates meeting the selected electronic and feasibility constraints |
| diversity | nearest-neighbor distance, fingerprint coverage, or cluster coverage |

Generate the machine-readable report with:

```bash
python3 07_validation/collect_metrics.py
```

For a new batch, explicitly identify every stage artifact, the dataset, and a
same-protocol reference:

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

Use `--no-reference` when no scientifically compatible reference exists. Do
not directly compare values from different functionals, basis sets, solvents,
or relativistic protocols.

The committed `reports/pic16_s1/metrics.json` snapshot reports:

- 222 input, TD-summary, v20, v21, and v22 rows, with 100% required-field completeness;
- 191 no-imaginary-frequency rows, approximately 86.0%;
- 2 nondominated rows for the gap–`SOCME/|ΔE|` two-objective front;
- 14 nondominated rows for the gap–SOCME–ratio three-objective front;
- 175/222 candidates meeting or improving both the gap and ratio thresholds
  of reference row `pic16`, or 156/191 after the no-imaginary filter;
- 155/222 candidates also meeting the SOCME threshold, or 138/191 after the
  no-imaginary filter;
- zero rows in the joint default LogP 1–5 and TPSA 20–120 Å² window.

The selected reference row in `05_analysis/data/publishinpdftv22.csv` has
`Imaginary=Yes`. These counts are therefore electronic-property threshold
comparisons, not improvements over a validated geometric minimum.

The last result shows that the targeting plot defaults are unsuitable for the
S1 distribution and must not be treated as universal filters. The S1 LogP
median is approximately 14.08; select application-specific windows and run a
sensitivity analysis.

Recommended optimization order:

1. improve normal termination and parse coverage;
2. remove imaginary-frequency, wrong-charge, structurally invalid, and incomplete records;
3. optimize the electronic-property Pareto front;
4. apply LogP/TPSA/SA/MW feasibility constraints;
5. control redundancy through chemical-space sampling;
6. recompute representative Pareto candidates at a higher level of theory.

## 12. Known limitations

1. The project is a file-connected script workflow, not a unified Python package or workflow engine.
2. `03_result_parsing/parse_orca_opt.py` does not independently enforce normal termination.
3. TD-to-SMILES mapping depends on source CSV row order and task IDs.
4. `03_result_parsing/merge_dft_csv.py` uses the union of IDs and can produce partial rows.
5. Some plots have code-level axis ranges, DPI, palettes, and random seeds that are not exposed through CLI options.
6. `config/agent_interface.json` is declarative and does not rewrite the four ORCA generators; the generators are the execution authority.
7. Parsers target ORCA 6.0.1 output conventions; another major ORCA version may require regex updates.
8. SA Score components have separate third-party provenance; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
9. Environment and metric reports can be generated automatically, but complete run provenance must still be stored under `reports/runs/`.

## 13. Validation

Python syntax:

```bash
find . -type f -name '*.py' \
  -exec python3 -c 'import pathlib,sys; p=pathlib.Path(sys.argv[1]); compile(p.read_text(encoding="utf-8"), str(p), "exec")' {} \;
```

Shell syntax:

```bash
find . -type f -name '*.sh' -exec bash -n {} \;
```

SA Score test:

```bash
python3 07_validation/UnitTestSAScore.py
```

Interface, environment, and core data:

```bash
python3 workflow.py paths --check
python3 07_validation/check_environment.py --output reports/environment.json
sha256sum -c data/core/pic16_s1/SHA256SUMS
python3 07_validation/collect_metrics.py
```

Plot smoke test:

```bash
python3 05_analysis/plotters/plot_pareto_s1ti_socmegap.py \
  data/core/pic16_s1/final/pic16inpdftv22.csv \
  --output results/figures/pic16_s1_pareto.png
```

## 14. GitHub publication

```bash
git init
git remote add origin https://github.com/letsgo-2/meloculardftcal.git
git add .
git status
git commit -m "Initial public workflow release"
git branch -M main
git push -u origin main
```

Before publishing:

- choose and add an appropriate `LICENSE`;
- confirm that reference molecules and research data can be released;
- scan for usernames, credentials, tokens, and cluster-specific absolute paths;
- do not upload ORCA binaries or other restricted software;
- do not upload full GBW/log collections, bulk figures, or sensitive datasets;
- rerun a small end-to-end batch;
- preserve the exact computational protocol and run manifest.
