# pic16 / S1 core dataset

This directory contains the only complete research batch bundled with the
public workflow. In `config/pictos.txt`, batch `pic16` maps to display label
`S1`.

All paths below are relative to the repository root.

| Path | Rows | Role |
|---|---:|---|
| `data/core/pic16_s1/input/pic16.csv` | 222 | generated molecules, scaffolds, SMILES, and generation scores |
| `data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv` | 222 | parsed TD-DFT/SOC summary |
| `data/core/pic16_s1/intermediate/pic16inpdftv20.csv` | 222 | merged ground-state and TD/SOC properties |
| `data/core/pic16_s1/intermediate/pic16inpdftv21.csv` | 222 | v20 merged back into the generation table |
| `data/core/pic16_s1/final/pic16inpdftv22.csv` | 222 | v21 plus RDKit descriptors and SA Score |

A historical S1-labelled export in the source workspace has the same
scientific values as the bundled v22 table; only its ID prefix was renamed
from `pic16_` to `s1_`. The `pic16_` IDs are retained here because they preserve
direct source-row and calculation-filename provenance.

Raw ORCA logs, GBW wavefunctions, optimized geometries, Slurm outputs, Excel
workbooks, and bulk images are intentionally excluded. Consequently, normal
ORCA termination counts cannot be independently reconstructed from this
dataset alone.

The optional electronic-property benchmark is row `pic16` in
`05_analysis/data/publishinpdftv22.csv`. That row has `Imaginary=Yes`; use it
only as an electronic threshold reference, not as a validated minimum.

Verify copied files:

```bash
sha256sum -c data/core/pic16_s1/SHA256SUMS
```

Regenerate the quantitative report:

```bash
python3 07_validation/collect_metrics.py
```

The report is written to `reports/pic16_s1/metrics.json`.
