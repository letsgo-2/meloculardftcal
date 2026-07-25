# pic16 / S1 reports

- `reports/pic16_s1/metrics.json`: committed data-quality and scientific
  performance snapshot, including Pareto sizes, property distributions, and
  comparison with reference row `pic16`.
- `05_analysis/data/publishinpdftv22.csv`: electronic-property reference table.
  Its `pic16` row has `Imaginary=Yes` and is not a validated minimum.
- `results/figures/pic16_s1/`: generated figures; ignored by Git.
- `reports/runs/`: optional per-run provenance manifests.

Regenerate metrics from the repository root:

```bash
python3 07_validation/collect_metrics.py
```

Inspect configured paths:

```bash
python3 workflow.py paths --check
```
