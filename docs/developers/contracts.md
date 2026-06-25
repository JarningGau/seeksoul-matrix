# Stage contracts

Normative input/output contracts for seeksoul-matrix workflow stages. Implementation details, QC field definitions, and chunk semantics live in [`doc-system.md`](doc-system.md) and linked pages.

## Main stage order

`fastp_split -> demux_extract_bc -> regroup_shards -> bismark_align -> bam_sort -> (count_mapped_reads -> estimated_cells)? -> split_bams -> merge_fr_bams -> bam_to_allc -> saturation -> qc_summary`

All twelve stages are implemented in this repository revision. Read-order vs analysis chunks and barcode-selection routing: [`chunk_model.md`](chunk_model.md).

### `fastp_split` {#fastp_split}

Purpose: quality control and chunking of paired FASTQ input.

Inputs:

- raw `R1 FASTQ`
- raw `R2 FASTQ`

Outputs:

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

Contract:

- Downstream `demux_extract_bc` reads paired shards from `shard_fastq/`.

See also: [stage notes](stage_notes/fastp_split.md) · [chunk model](chunk_model.md)

### `demux_extract_bc` {#demux_extract_bc}

Purpose: DD-MET5 barcode extraction, C→T QC, adapter trimming, and forward/reverse paired FASTQ output with barcode-prefix sub-sharding.

Inputs:

- `work/<sample>/shard_fastq/<readchunk>.R1.fq.gz` and paired R2 from `fastp_split`
- cell barcode whitelist: `whitelist/DD-MET5/U3CB_methylation.txt.gz` (default)
- workflow key `split_fastq_prefix_bases` (default `1`)

Per-read-chunk outputs under `work/<sample>/demux/`:

| File | Columns |
|------|---------|
| `<readchunk>.linker.tsv` | `CR`, `UB`, `C`, `T` |
| `<readchunk>.stats.json` | demux and C→T QC metrics |

When `split_fastq_prefix_bases > 0`, sub-shards under `work/<sample>/demux/shards/`:

| File |
|------|
| `<readchunk>__<prefix>.forward_1.fq.gz` / `.forward_2.fq.gz` |
| `<readchunk>__<prefix>.reverse_1.fq.gz` / `.reverse_2.fq.gz` |

Sample-level output (after all read-chunks):

| File | Columns |
|------|---------|
| `qc.CtoT.tsv` | `CR`, `C`, `T`, `CtoT` |

Contract:

- When prefix sub-sharding is enabled, every read-order chunk writes sub-shards for each observed barcode prefix.

See also: [stage notes](stage_notes/demux_extract_bc.md) · [chunk model](chunk_model.md) · [QC metrics](qc_metrics.md)

### `regroup_shards` {#regroup_shards}

Purpose: concatenate demux prefix sub-shards into one analysis shard per barcode prefix.

Inputs:

- `work/<sample>/demux/shards/<readchunk>__<prefix>.{forward,reverse}_{1,2}.fq.gz`

Outputs under `work/<sample>/demux/`:

| File |
|------|
| `<prefix>.forward_1.fq.gz` / `.forward_2.fq.gz` |
| `<prefix>.reverse_1.fq.gz` / `.reverse_2.fq.gz` |
| `chunks.tsv` |

Contract:

- Analysis `chunk_id` equals barcode prefix (e.g. `A`, `G`, `T` when `split_fastq_prefix_bases=1`).
- Downstream stages use top-level `demux/<prefix>.*` only (not `demux/shards/`).

See also: [stage notes](stage_notes/regroup_shards.md) · [chunk model](chunk_model.md)

### `bismark_align` {#bismark_align}

Purpose: Bismark bisulfite alignment of demux forward/reverse paired FASTQ.

Inputs:

- `work/<sample>/demux/<chunk>.forward_1.fq.gz` / `.forward_2.fq.gz`
- `work/<sample>/demux/<chunk>.reverse_1.fq.gz` / `.reverse_2.fq.gz`
- `bismark_ref`: parent directory of `Bisulfite_Genome/`

Per-chunk outputs under `work/<sample>/align/`:

| File |
|------|
| `<chunk>.forward_1_bismark_bt2_pe.bam` |
| `<chunk>.forward_1_bismark_bt2_PE_report.txt` |
| `<chunk>.reverse_1_bismark_bt2_pe.bam` |
| `<chunk>.reverse_1_bismark_bt2_PE_report.txt` |

Contract:

- One forward and one reverse BAM (plus reports) per analysis chunk.

See also: [stage notes](stage_notes/bismark_align.md)

### `bam_sort` {#bam_sort}

Purpose: name-sort Bismark paired-end BAMs for downstream per-cell splitting.

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe.bam`

Per-chunk outputs under `work/<sample>/align/`:

| File |
|------|
| `<chunk>.forward_1_bismark_bt2_pe_sortbyname.bam` |
| `<chunk>.reverse_1_bismark_bt2_pe_sortbyname.bam` |

Contract:

- Input BAMs must be name-sorted for `split_bams`; unsorted source BAMs are retained alongside sortbyname outputs.

See also: [stage notes](stage_notes/bam_sort.md)

### `count_mapped_reads` {#count_mapped_reads}

Purpose: count aligned reads per cell barcode from unsorted Bismark BAMs (methylation-only path).

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe.bam`

Per-strand outputs under `work/<sample>/align/`:

| File | Columns |
|------|---------|
| `<chunk>.{forward,reverse}_1_bismark_bt2_pe_cb_aligned_reads_counts.csv` | `barcode`, `aligned_reads` |

Contract:

- Uses unsorted BAMs (not sortbyname outputs).
- Skipped when workflow uses `gexcb` mode ([chunk model](chunk_model.md#barcode-selection)).

See also: [stage notes](stage_notes/count_mapped_reads.md) · [chunk model](chunk_model.md)

### `estimated_cells` {#estimated_cells}

Purpose: merge per-chunk/strand barcode counts and filter to called cells (methylation-only path).

Inputs:

- `work/<sample>/align/*_cb_aligned_reads_counts.csv`
- workflow key `expected_cell_num` (default `3000`; threshold filtering)
- optional workflow key `force_cell_num`: top N barcodes by `aligned_reads` among barcodes with reads > 0; when set, overrides `expected_cell_num` threshold filtering

Sample-level outputs under `work/<sample>/cells/`:

| File | Columns |
|------|---------|
| `merged_barcode_counts.csv` | `barcode`, `aligned_reads` (summed) |
| `filtered_barcode_read_counts.csv` | `aligned_reads`, `barcode` |
| `filtered_barcode` | one barcode per line |

Contract:

- Skipped when workflow uses `gexcb` mode ([chunk model](chunk_model.md#barcode-selection)).

See also: [stage notes](stage_notes/estimated_cells.md) · [chunk model](chunk_model.md)

### `split_bams` {#split_bams}

Purpose: split name-sorted Bismark BAMs into per-cell BAM files.

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe_sortbyname.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe_sortbyname.bam`
- barcode list (**one of**):
  - `work/<sample>/cells/filtered_barcode` (methylation-only path)
  - `gexcb` workflow path (RNA filtered barcodes)

Per-strand outputs under `work/<sample>/split_bams/`:

| Path |
|------|
| `<chunk>.forward_1/<barcode>.bam` |
| `<chunk>.reverse_1/<barcode>.bam` |
| `<chunk>.{forward,reverse}_1/<chunk>.{forward,reverse}_1_filtered_barcode` |
| `<chunk>.{forward,reverse}_1/<chunk>.{forward,reverse}_1_filtered_barcode_reads_counts.csv` |

Contract:

- Barcode list source is mutually exclusive; see [chunk model](chunk_model.md#barcode-selection).
- Input BAM must be sorted by read name.

See also: [stage notes](stage_notes/split_bams.md) · [chunk model](chunk_model.md)

### `merge_fr_bams` {#merge_fr_bams}

Purpose: merge forward and reverse per-cell split BAMs into one BAM per barcode.

Inputs (per chunk):

- `work/<sample>/split_bams/<chunk>.forward_1/*.bam`
- `work/<sample>/split_bams/<chunk>.reverse_1/*.bam`
- `*_filtered_barcode` and `*_filtered_barcode_reads_counts.csv` under both strand directories

Outputs under `work/<sample>/split_bams/merged/`:

| Path | Columns (where applicable) |
|------|----------------------------|
| `<chunk>_merged_fr_bam/<barcode>.bam` | |
| `<chunk>_merge_filtered_barcode` | one barcode per line |
| `<chunk>_merge_filtered_barcode_reads_counts.csv` | `reads_counts`, `barcode` |

Contract:

- One merged BAM per barcode per analysis chunk when both strands have reads.

See also: [stage notes](stage_notes/merge_fr_bams.md) · [chunk model](chunk_model.md)

### `bam_to_allc` {#bam_to_allc}

Purpose: convert merged per-cell BAMs to ALLC format.

Inputs (per chunk):

- `work/<sample>/split_bams/merged/<chunk>_merged_fr_bam/*.bam`
- `work/<sample>/split_bams/merged/<chunk>_merge_filtered_barcode`
- `genome_fa` (workflow JSON)
- `chrom_size_path` (workflow JSON)

Outputs under `work/<sample>/allcools/`:

| Path |
|------|
| `<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz` |
| `<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz.count.csv` |

Contract:

- Only barcodes listed in `<chunk>_merge_filtered_barcode` are processed.

See also: [stage notes](stage_notes/bam_to_allc.md)

### `saturation` {#saturation}

Purpose: estimate sample-level genome-fraction saturation curve from pre-dedup per-cell BAMs.

Inputs:

- `work/<sample>/shard_fastq/fastp.json`
- `work/<sample>/split_bams/merged/<chunk>_merged_fr_bam/<barcode>.bam` (all chunks)
- cell read-count table:
  - methylation-only: `work/<sample>/cells/filtered_barcode_read_counts.csv`
  - gexcb: `work/<sample>/split_bams/merged/*_merge_filtered_barcode_reads_counts.csv`
- `chrom_size_path`

Outputs under `work/<sample>/qc/saturation/`:

| File |
|------|
| `saturation_curve.png` |
| `saturation_summary.tsv` |

Contract:

- Single sample-level job; gathers per-cell BAMs across all analysis chunks.

See also: [stage notes](stage_notes/saturation.md) · [chunk model](chunk_model.md) · [QC metrics](qc_metrics.md)

### `qc_summary` {#qc_summary}

Purpose: gather per-cell and sample-level QC metrics into summary tables.

Inputs:

- `work/<sample>/shard_fastq/fastp.json`
- `work/<sample>/demux/*.stats.json`, `work/<sample>/demux/qc.CtoT.tsv`
- `work/<sample>/align/*_{forward,reverse}_1_bismark_bt2_PE_report.txt`
- `work/<sample>/qc/saturation/saturation_summary.tsv`
- `work/<sample>/allcools/*_merged_fr_bam_allcools/*_allc.gz.count.csv`
- `work/<sample>/allcools/*_merged_fr_bam_allcools/*_allc.gz` (sibling of each count sidecar; used for `mito_CG_mc_rate`)
- cell read-count table (methylation-only or gexcb; see [chunk model](chunk_model.md#barcode-selection))
- optional `cbcsv` (workflow JSON or `make_cmd.py --cbcsv`): methylation ↔ GEX barcode map for `gex_cb` column in `cells_summary.tsv`
- optional `mito_chromosomes` (workflow JSON or `make_cmd.py --mito-chromosomes`; default `chrM`): mitochondrial contig names for `mito_CG_mc_rate`

Outputs under `work/<sample>/summary/`:

| File |
|------|
| `cells_summary.tsv` |
| `sample_summary.tsv` |
| `wgs_summary.csv` |

Contract:

- Single sample-level job; globs ALLCools count sidecars across analysis chunks (gather-only).
- `cells_summary.tsv`: one row per called cell; `sample_summary.tsv` and `wgs_summary.csv`: one row per sample.

See also: [stage notes](stage_notes/qc_summary.md) · [chunk model](chunk_model.md) · [QC metrics](qc_metrics.md)

### `allc_to_matrix` {#allc_to_matrix}

Purpose: convert per-cell ALLC into a MethSCAn-compatible CSR sparse matrix store for downstream methylation analysis.

Inputs:

- `work/<sample>/allcools/<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz` (gather across analysis chunks)
- optional barcode list (**one of**):
  - `--cell-names` / explicit file
  - `work/<sample>/cells/filtered_barcode` (methylation-only path)
  - `work/<sample>/split_bams/merged/*_merge_filtered_barcode` (gexcb path)
  - if none present: all discovered ALLC barcodes (warning logged)

Workflow keys (when enabled via `run_meth_analysis`):

| Key | Default |
|-----|---------|
| `meth_context` | `CG` |
| `meth_chunksize` | `10000000` |
| `meth_round_sites` | `false` |
| `meth_main_chroms_only` | `false` |
| `meth_exclude_contigs` | `""` |

Outputs under `work/<sample>/meth/matrix/`:

| Path | Description |
|------|-------------|
| `{chrom}.npz` | CSR sparse matrix (`int8`: +1 methylated, -1 unmethylated) |
| `column_header.txt` | one cell barcode per line |
| `cell_stats.csv` | `cell_name`, `n_obs`, `n_meth`, `global_meth_frac` |
| `run_info.json` | stage parameters and runtime metadata |

Contract:

- Single sample-level job; globs ALLC across analysis chunks; rejects duplicate barcodes across chunks.
- Optional post-`qc_summary` stage; gated by workflow key `run_meth_analysis` (default `false`).
- MethSCAn `prepare` parity validated on `work/dd-met5-example` (50 cells, all-context comparison passed).

See also: [stage notes](stage_notes/allc_to_matrix.md) · [chunk model](chunk_model.md)

### `meth_smooth` {#meth_smooth}

Purpose: tricube-weighted pseudobulk smoothing over per-position methylation fractions (MethSCAn `smooth` equivalent).

Inputs:

- `work/<sample>/meth/matrix/{chrom}.npz` (CSR store from `allc_to_matrix`)

Workflow keys (when enabled via `run_meth_analysis`):

| Key | Default |
|-----|---------|
| `meth_smooth_bandwidth` | `1000` |
| `meth_smooth_use_weights` | `false` |

Outputs under `work/<sample>/meth/matrix/smoothed/`:

| Path | Description |
|------|-------------|
| `{chrom}.csv.gz` | two columns: genomic position, smoothed methylation fraction |
| `run_info.json` | stage parameters and runtime metadata |

Contract:

- Single sample-level job; reads the active CSR matrix store under `meth/matrix/`.
- Requires prior `allc_to_matrix` output.

See also: [stage notes](stage_notes/meth_smooth.md) · [chunk model](chunk_model.md)

### `meth_scan` {#meth_scan}

Purpose: sliding-window scan for variably methylated regions (VMRs; MethSCAn `scan` equivalent).

Inputs:

- `work/<sample>/meth/matrix/{chrom}.npz`
- `work/<sample>/meth/matrix/smoothed/{chrom}.csv.gz`

Workflow keys (when enabled via `run_meth_analysis`):

| Key | Default |
|-----|---------|
| `meth_scan_bandwidth` | `2000` |
| `meth_scan_stepsize` | `100` |
| `meth_scan_var_threshold` | `0.02` |
| `meth_scan_min_cells` | `6` |
| `meth_scan_bridge_gaps` | `0` |
| `meth_matrix_cores` | `8` |

Outputs under `work/<sample>/meth/vmr/`:

| Path | Description |
|------|-------------|
| `vmrs.bed` | tab-separated: `chromosome`, `VMR_start`, `VMR_end`, `variance`, `n_sites`, `n_cells` (no header by default) |
| `run_info.json` | stage parameters, variance threshold, VMR counts |

Contract:

- Single sample-level job; requires prior `meth_smooth` output.
- Global variance threshold is determined on the largest chromosome (by CSR file size), then applied genome-wide.
- Fails if no VMR passes `--min-cells`.

See also: [stage notes](stage_notes/meth_scan.md) · [chunk model](chunk_model.md)
