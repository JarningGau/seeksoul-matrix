# QC metrics

Field definitions, demux statistics schema, saturation model, and summary-table columns. For input/output paths, see [`contracts.md`](contracts.md).

## Demux per-read-chunk outputs

### `linker.tsv`

One row per UMI-deduped QC read under `work/<sample>/demux/<readchunk>.linker.tsv`:

| Column | Description |
|--------|-------------|
| `CR` | corrected cell barcode |
| `UB` | 12 bp UMI |
| `C` | convertible-base count |
| `T` | converted-base count |

### `qc.CtoT.tsv`

Sample-level file after all read-chunks complete:

| Column | Description |
|--------|-------------|
| `CR` | corrected cell barcode |
| `C` | aggregated convertible bases |
| `T` | aggregated converted bases |
| `CtoT` | `T / (C+T)` rounded to 3 decimals |

Merged from all `<readchunk>.linker.tsv`, aggregated by `CR`.

### Read name format

`{CB}_{UB}_{forward|reverse}_{alt}_{orig_name}:{CB}` where `alt` is `M` for exact whitelist match or HD=1 correction tag (e.g. `3A`).

## Demux `stats.json` schema

Per read-chunk file `<readchunk>.stats.json`. Nested **funnel** (read retention) plus parallel **ct** (C→T QC on TTT-insert reads):

```
funnel.total
├─ funnel.barcode_rejected
├─ funnel.barcode_ambiguous
└─ funnel.barcode_passed.total  (= exact + corrected)
   ├─ funnel.barcode_passed.exact
   ├─ funnel.barcode_passed.corrected
   ├─ funnel.barcode_passed.unknown_chain
   ├─ funnel.barcode_passed.too_short
   ├─ funnel.barcode_passed.chimeric_filtered.total
   │  ├─ funnel.barcode_passed.chimeric_filtered.forward
   │  └─ funnel.barcode_passed.chimeric_filtered.reverse
   └─ funnel.barcode_passed.valid.total
      ├─ funnel.barcode_passed.valid.forward
      └─ funnel.barcode_passed.valid.reverse
```

| Section | Fields | Description |
|---------|--------|-------------|
| (root) | `chunk_id`, `chemistry`, `filter_ch`, `input_r1`, `input_r2` | metadata |
| `funnel` | nested object above | demux read retention funnel (all counts are reads) |
| `ct` | `num_17lme`, `rate_17lme`, `ct_reads`, `ct_umi_dedup`, `ct_convertible_bases`, `ct_converted_bases`, `CtoT` | C→T QC funnel (parallel to FASTQ output; TTT-insert only; fields in processing order) |

**Accounting:** `funnel.total = barcode_rejected + barcode_ambiguous + barcode_passed.total`; `barcode_passed.total = unknown_chain + too_short + chimeric_filtered.total + valid.total` (when `unknown_chain` reads leave before length filter).

`CtoT` is rounded to 3 decimal places; fractions use 6 decimal places.

## Saturation model

Workflow parameters (defaults in parentheses):

| Key | Default | Role |
|-----|---------|------|
| `saturation_reads_threshold` | `100` | minimum aligned reads for HQ cell inclusion |
| `saturation_max_cells` | `100` | maximum HQ cells used for curve estimation |
| `saturation_sample_seed` | `42` | random seed when sampling HQ cells above `saturation_max_cells` |
| `saturation_linear_r2_threshold` | `0.99` | linear-fit R² at/above which linear extrapolation is used |

### Method

- Build per-cell per-base depth histogram from aligned primary reads across **all chunks** for each barcode; do not use deduplicated ALLC (`cov` is post-UMI-dedup and unsuitable for subsampling simulation).
- Depth is counted per **sequencing fragment** (read pair): PE mates share a query name and contribute at most once to each reference position's depth (overlapping mates are de-duplicated).
- Genome fraction at subsample depth `f`: `sum_d count(d) × (1 − (1−f)^d) / G`, where `count(d)` is the number of reference positions with fragment depth `d`, and `G` is genome size from `chrom_size_path`.
- HQ cells: `aligned_reads > saturation_reads_threshold` and per-cell BAM present. If HQ count exceeds `saturation_max_cells`, randomly sample `saturation_max_cells` cells (fixed `saturation_sample_seed`; HQ list sorted before sampling).
- For each subsample fraction (1%–100%), compute per-cell genome fraction; aggregate with **median** and **IQR** (asymmetric `Q1`/`Q3` error bars) across sampled cells.
- Extrapolation (`f` beyond 1×): fit both a through-origin linear model `y = m × f` and a saturation curve `y = a × (1 − exp(−b × f))` to the median curve.
  - If linear fit `R² ≥ saturation_linear_r2_threshold`: `extrapolation_model = linear`, `predicted@2× = m × 2`, `theoretical_max = saturation_rate = NA`.
  - Otherwise: `extrapolation_model = saturation`, `theoretical_max = a`, `predicted@2× = a × (1 − exp(−2b))`, `saturation_rate = median@100% / a × 100`.
- Plot x-axis: subsample fraction × fastp `total_bases` / 1e9 (Gbp); 2× prediction at twice the observed sequencing depth.
- Y-axis: median **genome fraction** (0–1 in TSV; plot as %).

### `saturation_summary.tsv` columns

`sample_id`, `observed_median_genome_fraction`, `theoretical_max_median_genome_fraction`, `predicted_median_genome_fraction_at_2x`, `saturation_rate`, `extrapolation_model`, `hq_cell_count`, `sampled_cell_count`, `sample_seed`.

## `qc_summary` output tables

### `cells_summary.tsv`

One row per cell. Core columns:

`cell_barcode`, `aligned_reads`, `CtoT`, `total_cpg_number`, `genome_cov`, `genome_cov_raw_umi`, `genome_cov_new_umi`, `cell_saturation`, `CG_mc_rate`, `mito_CG_mc_rate`, `CH_mc_rate`, `CHG_mc_rate`, `CHH_mc_rate`, `CA_mc_rate`, `CT_mc_rate`, `CC_mc_rate`, optional `gex_cb`.

Aggregated `*_mc_rate` values pool `mc` and `cov` across trinucleotide contexts:

- **CG:** `C?G` with second base `G`
- **CH:** non-CG
- **CHG / CHH:** third-base split within CH
- **CA / CT / CC:** second-base split within CH

**`mito_CG_mc_rate`:** per-cell mitochondrial CG methylation from the sibling `*_allc.gz` (not `count.csv`). Filter `chrom ∈ mito_chromosomes` (default `chrM`) and CG contexts (`CGA`, `CGC`, `CGG`, `CGT`); `sum(mc) / sum(cov)`. `NA` when the cell has no mitochondrial CG coverage.

Omits per-context `{context}_mc_rate`, `weighted_mc_rate`, `total_mc`, `total_cov`, `total_number`, and `{context}_mc/cov/number`.

### `sample_summary.tsv`

One-row internal QC summary (selected metric subset):

- **fastp:** `raw_reads`, `total_bases`, `duplication_rate`
- **demux:** `valid_barcode_rate` (`barcode_passed.total` / `funnel.total`), `valid_demux_rate` (`barcode_passed.valid.total` / `funnel.total`), `barcode_corrected_fraction`, `dropped_too_short`, `dropped_chimeric`, `forward_reads`, `reverse_reads`, `rate_17lme`, `CtoT`, plus `rate_7f`, `rate_7f17lme`, `cc_mean` as `NA` until demux stats extended
- **bismark:** `mapped_to_genome`, `confidently_mapped`, `cpg_methylation_rate`, `chg_methylation_rate`, `chh_methylation_rate` (no `unknown_methylation_rate`)
- **mito:** `mito_CG_mc_rate` — pooled `sum(mc) / sum(cov)` across mitochondrial CG sites in **called cells** only (same contig/context filter as per-cell column); `NA` when no coverage
- **cells:** `estimated_cells`, `reads_in_cells`, `fraction_reads_in_cells`
- **saturation:** `observed_median_genome_fraction`, `theoretical_max_median_genome_fraction`, `saturation_rate`, `extrapolation_model`, `sampled_cell_count`, `sample_seed` (no `hq_cell_count`, no `predicted_median_genome_fraction_at_2x`)
- **median called-cell QC:** `median_genome_cov`, `median_total_cpg_number`, `median_aligned_reads`, `median_cell_saturation`

### `wgs_summary.csv`

SeekSoul-compatible one-row CSV; max-cell and bulk-CpG columns are `NA`.
