# saturation

Implementation notes for [`saturation`](../contracts.md#saturation). Not normative for I/O.

## Behavior

- Single sample-level job (no per-chunk parallelism); gathers pre-dedup per-cell BAMs across all analysis chunks.
- Build per-cell per-base depth histogram from aligned primary reads; subsample fractions 1%–100% with median + IQR aggregation across sampled HQ cells.
- Plot: observed median genome fraction (IQR error bars) + fitted line/curve + 2× prediction; x-axis sequencing depth (Gbp) from fastp `summary.after_filtering.total_bases` (fallback `before_filtering`).
- HQ cell list sorted before random subsampling when count exceeds `saturation_max_cells`.

Formulas, extrapolation rules, and output column definitions: [`qc_metrics.md`](../qc_metrics.md#saturation-model).

## Defaults and workflow keys

See [`qc_metrics.md`](../qc_metrics.md#saturation-model) for `saturation_reads_threshold`, `saturation_max_cells`, `saturation_sample_seed`, `saturation_linear_r2_threshold`.

## Skip / incremental re-run

Supports `--dry-run`. Re-run overwrites `qc/saturation/` outputs.

## Toolchain and environment

`scripts/saturation.py`; uses `chrom_size_path` for genome size `G`.

## SeekSoulMethyl / dbit-matrix alignment

dbit-matrix CpG-coverage saturation style, adapted for DD-MET5 single-cell genome-fraction breadth.

## Out of scope

Does not depend on `merge_sc_metrics` or ALLCools count sidecars.
