# qc_summary

Implementation notes for [`qc_summary`](../contracts.md#qc_summary). Not normative for I/O.

## Behavior

- Glob `allcools/*_merged_fr_bam_allcools/*_allc.gz.count.csv` across analysis chunks; gather-only (barcodes are disjoint by prefix; see [`chunk_model.md`](../chunk_model.md)).
- Median-cell stats: among called cells with ALLCools metrics, sort by `genome_cov` descending and take the median index.
- `fraction_reads_in_cells` = `reads_in_cells` / pooled Bismark unique aligned reads.
- Optional workflow key `cbcsv`: path to `m_cb,gex_cb` map for `gex_cb` column (not the `gexcb` barcode list).

Column definitions for output tables: [`qc_metrics.md`](../qc_metrics.md#qc_summary-output-tables).

## Defaults and workflow keys

| Key | Role |
|-----|------|
| `cbcsv` | optional methylation ↔ GEX barcode map for `gex_cb` column |

Cell read-count input path depends on barcode-selection mode — see [`chunk_model.md`](../chunk_model.md#barcode-selection).

## Skip / incremental re-run

Single sample-level job; supports `--dry-run`.

## Toolchain and environment

`scripts/qc_summary.py`.

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `step3_merge_sc_metrics` + `step4_wgs_summary`.

## Out of scope

Per-cell JSON, HTML reports, or bulk merged ALLC metrics.
