# qc_summary

Implementation notes for [`qc_summary`](../contracts.md#qc_summary). Not normative for I/O.

## Behavior

- Glob `allcools/*_merged_fr_bam_allcools/*_allc.gz.count.csv` across analysis chunks; gather-only (barcodes are disjoint by prefix; see [`chunk_model.md`](../chunk_model.md)).
- **`mito_CG_mc_rate`:** parse sibling `*_allc.gz` per cell; keep rows with `chrom ∈ mito_chromosomes` and CG context (`context[1] == G`). Sample row pools mc/cov across called cells only.
- Median-cell stats: among called cells with ALLCools metrics, sort by `genome_cov` descending and take the median index.
- `fraction_reads_in_cells` = `reads_in_cells` / pooled Bismark unique aligned reads.
- Optional workflow key `cbcsv`: path to `m_cb,gex_cb` map for `gex_cb` column (not the `gexcb` barcode list).

Column definitions for output tables: [`qc_metrics.md`](../qc_metrics.md#qc_summary-output-tables).

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `cbcsv` | — | optional methylation ↔ GEX barcode map for `gex_cb` column |
| `mito_chromosomes` | `chrM` | comma-separated mitochondrial contigs for `mito_CG_mc_rate` |

Cell read-count input path depends on barcode-selection mode — see [`chunk_model.md`](../chunk_model.md#barcode-selection).

## Skip / incremental re-run

Single sample-level job; supports `--dry-run`.

## Toolchain and environment

`scripts/qc_summary.py`.

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `step3_merge_sc_metrics` + `step4_wgs_summary`.

## Validation

- Local: `work/dd-met5-example` (50 cells, `force_cell_num=50`); sixteen-stage e2e includes `qc_summary` (2026-06-25).
- Biological sanity **passed**: CtoT≈0.997; sample CpG mc≈77%; demux funnel and align rates consistent with DD-MET5 test chemistry; `cells_summary.tsv` metrics plausible for mouse bisulfite data.
- Slurm: aggregate `qc_summary` job in `dd_met5_test.json` DAG cluster-tested (HPC, 2026-06-25).

## Out of scope

Per-cell JSON, HTML reports, or bulk merged ALLC metrics.
