# estimated_cells

Implementation notes for [`estimated_cells`](../contracts.md#estimated_cells). Not normative for I/O.

## Behavior

- Merges per-chunk/strand `*_cb_aligned_reads_counts.csv` into sample-level barcode tables.
- Default filtering: `expected_cell_num=3000`; 99th-percentile index × 0.1 read threshold (aligned with SeekSoulMethyl `step3_estimated_cells.py`).
- Optional `force_cell_num`: when set, take top N barcodes by `aligned_reads` among barcodes with reads > 0 (aligned with SeekSoulMethyl `cell_identify.R` force-cell mode). `force_cell_num` overrides `expected_cell_num` filtering; both keys may appear in workflow JSON.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `expected_cell_num` | `3000` | target cell count for percentile threshold |
| `force_cell_num` | (optional) | force top-N cells by read count |

## Skip / incremental re-run

Sample-level single job; re-run overwrites `cells/` outputs.

## Toolchain and environment

`scripts/estimated_cells.py`.

## SeekSoulMethyl / dbit-matrix alignment

Threshold logic aligned with SeekSoulMethyl `step3_estimated_cells.py` and `cell_identify.R`.

## Out of scope

Skipped when workflow uses `gexcb` mode (see [`chunk_model.md`](../chunk_model.md#barcode-selection)).
