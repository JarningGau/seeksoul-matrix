# Stage contracts

Normative input/output contracts for seeksoul-matrix workflow stages.

## Main stage order (planned)

`fastp_split -> demux_extract_bc -> regroup_shards -> bismark_align -> bam_sort -> (count_mapped_reads -> estimated_cells)? -> split_bams -> merge_fr_bams -> bam_to_allc -> saturation -> qc_summary`

`fastp_split`, `demux_extract_bc`, `regroup_shards`, `bismark_align`, `bam_sort`, `count_mapped_reads`, `estimated_cells`, `split_bams`, `merge_fr_bams`, `bam_to_allc`, `saturation`, and `qc_summary` are implemented in this repository revision.

### Chunk model (barcode-prefix shards)

| Concept | Key / ID | Role |
|---------|----------|------|
| Read-order chunk | numeric `0001..N` from `number_of_split_parts` | fastp split + parallel demux input only |
| Analysis chunk | barcode prefix from `split_fastq_prefix_bases` (default `1`) | all stages from `regroup_shards` through `qc_summary` |

With DD-MET5 barcodes (no `C`), `split_fastq_prefix_bases=1` yields up to 3 analysis chunks (`A`, `G`, `T`). Each cell barcode maps to exactly one prefix; regrouped shards are disjoint by barcode.

Stages that consume all analysis chunks (e.g. `saturation`, `qc_summary`) **gather** per-chunk outputs by globbing chunk directories. They do **not** deduplicate barcodes across chunks — a barcode never appears in more than one prefix shard.

Manifest: `work/<sample>/demux/chunks.tsv` (`chunk_id`, `prefix`, `stream`, `subshard_count`, `subshard_paths`).

### Barcode selection mode (mutually exclusive)

Downstream step3 stages require a cell barcode list for `split_bams`. Configure **exactly one** of:

| Mode | Workflow key | Stage path |
|------|--------------|------------|
| Methylation-only (default) | `expected_cell_num` (default `3000`) | `count_mapped_reads` → `estimated_cells` → `split_bams` |
| RNA + methylation | `gexcb` (path to RNA filtered barcodes) | `split_bams` only |

If neither key is set, `expected_cell_num=3000` applies. Setting both `gexcb` and `expected_cell_num` is an error.

### `fastp_split`

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

- Adapter trimming is disabled (`--disable_adapter_trimming`).
- Chunk count is controlled by `number_of_split_parts` (passed to `fastp --split`); this controls **read-order parallelism** for demux only, not analysis chunk boundaries.
- Downstream `demux_extract_bc` will read paired shards from `shard_fastq/`.

### `demux_extract_bc`

Purpose: DD-MET5 barcode extraction, C→T QC, adapter trimming, and forward/reverse paired FASTQ output with optional barcode-prefix sub-sharding.

Inputs:

- `work/<sample>/shard_fastq/R1*.fq.gz` or `<readchunk>.R1.fq.gz` and paired R2 from `fastp_split`
- cell barcode whitelist: `whitelist/DD-MET5/U3CB_methylation.txt.gz` (default)
- `split_fastq_prefix_bases` (workflow key, default `1`; SeekSoul `split_fastq`): first N bases of corrected barcode for prefix sub-shards

Per-read-chunk outputs under `work/<sample>/demux/`:

| File | Columns / fields | Description |
|------|------------------|-------------|
| `<readchunk>.linker.tsv` | `CR`, `UB`, `C`, `T` | one row per UMI-deduped QC read; `CR` = corrected CB, `UB` = 12 bp UMI, `C`/`T` = convertible-base counts |
| `<readchunk>.stats.json` | grouped JSON | read-chunk demux and C→T QC metrics (see schema below) |

When `split_fastq_prefix_bases > 0`, demux FASTQ sub-shards under `work/<sample>/demux/shards/`:

| File | Description |
|------|-------------|
| `<readchunk>__<prefix>.forward_1.fq.gz` / `.forward_2.fq.gz` | forward-strand sub-shard for barcode prefix |
| `<readchunk>__<prefix>.reverse_1.fq.gz` / `.reverse_2.fq.gz` | reverse-strand sub-shard for barcode prefix |

Sample-level output (after all read-chunks complete):

| File | Columns | Description |
|------|---------|-------------|
| `qc.CtoT.tsv` | `CR`, `C`, `T`, `CtoT` | merge all `<chunk>.linker.tsv`, aggregate by `CR`; `CtoT = T / (C+T)` rounded to 3 decimals |

Read name format: `{CB}_{UB}_{forward|reverse}_{alt}_{orig_name}:{CB}` where `alt` is `M` for exact whitelist match or HD=1 correction tag (e.g. `3A`).

Contract:

- Structure `B17U12` on R1; CB whitelist with Hamming distance 1; ambiguous multi-match CBs are discarded (no multi rescue).
- C→T QC uses ME5 positions and 17L+ME validation aligned with SeekSoulMethyl `calculate_ct_conversion_rate`.
- Only reads with **TTT** insert signature (forward chain) contribute to `linker.tsv` and CtoT statistics; reverse (`CCC`) reads are excluded from conversion QC.
- Each `(CR, UB)` pair contributes at most one row to `linker.tsv`.
- **CH chimeric filtering** (`--filter-ch`, default `2`, workflow key `filter_ch`): after adapter trim and length check, drop read pairs where trimmed R1 or R2 exceeds the threshold for strand-specific CH patterns (SeekSoulMethyl `should_filter_read_ch_pattern`). `0` disables filtering. Does not affect `linker.tsv` / C→T QC (computed on raw R1 before trim).
- Spike-in output and multi-barcode rescue are out of scope.

**`<chunk>.stats.json` schema:**

Nested **funnel** (read retention, top to bottom) plus parallel **ct** (C→T QC on TTT-insert reads):

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

Accounting: `funnel.total = barcode_rejected + barcode_ambiguous + barcode_passed.total`; `barcode_passed.total = unknown_chain + too_short + chimeric_filtered.total + valid.total` (when `unknown_chain` reads leave before length filter).

`CtoT` is rounded to 3 decimal places; fractions use 6 decimal places.

### `regroup_shards`

Purpose: concatenate demux prefix sub-shards into one analysis shard per barcode prefix (SeekSoul `split_fastq` alignment).

Inputs:

- `work/<sample>/demux/shards/<readchunk>__<prefix>.{forward,reverse}_{1,2}.fq.gz`

Outputs under `work/<sample>/demux/`:

| File | Description |
|------|-------------|
| `<prefix>.forward_1.fq.gz` / `.forward_2.fq.gz` | regrouped forward-strand paired FASTQ |
| `<prefix>.reverse_1.fq.gz` / `.reverse_2.fq.gz` | regrouped reverse-strand paired FASTQ |
| `chunks.tsv` | manifest of analysis chunks and source sub-shards |

Contract:

- Gzip member concatenation (binary append) merges sub-shards per stream.
- Analysis `chunk_id` equals barcode prefix (e.g. `A`, `G`, `T` when `split_fastq_prefix_bases=1`).
- Downstream stages discover chunks from top-level `demux/<prefix>.*` files only (not `demux/shards/`).

### `bismark_align`

Purpose: Bismark bisulfite alignment of demux forward/reverse paired FASTQ using the seekgene [Bismark](https://github.com/seekgene/Bismark) fork (`--add_barcode`, `--add_umi`).

Inputs:

- `work/<sample>/demux/<chunk>.forward_1.fq.gz` / `.forward_2.fq.gz`
- `work/<sample>/demux/<chunk>.reverse_1.fq.gz` / `.reverse_2.fq.gz`
- `bismark_ref`: parent directory of `Bisulfite_Genome/` (passed to `bismark --genome`)

Per-chunk outputs under `work/<sample>/align/`:

| File | Description |
|------|-------------|
| `<chunk>.forward_1_bismark_bt2_pe.bam` | forward-strand Bismark paired-end BAM |
| `<chunk>.forward_1_bismark_bt2_PE_report.txt` | forward alignment report |
| `<chunk>.reverse_1_bismark_bt2_pe.bam` | reverse-strand BAM (`--pbat`) |
| `<chunk>.reverse_1_bismark_bt2_PE_report.txt` | reverse alignment report |

Contract:

- Uses seekgene Bismark (stock bioconda Bismark lacks `--add_barcode` / `--add_umi`); requires `bowtie2` on PATH (`check-bismark-env` verifies).
- `bismark_align.py` prepends the Bismark executable directory to subprocess `PATH` so Bismark can invoke `bowtie2` on Slurm compute nodes without pixi activation.
- Forward: `bismark --genome <bismark_ref> --parallel <N> -1 <fwd_r1> -2 <fwd_r2> -o <align_dir> -X <max_insert> --add_barcode --add_umi`.
- Reverse: same with `--pbat`.
- Default `bismark_parallel=8`, `bismark_max_insert=1000` (aligned with SeekSoulMethyl `step2.nf`).
- BAM read names retain demux format; Bismark writes `CB` and `UR` tags from read names.
- `samtools sort`, UMI dedup, and per-cell BAM splitting are out of scope for this stage.

### `bam_sort`

Purpose: name-sort Bismark paired-end BAMs so read pairs are adjacent for downstream per-cell splitting.

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe.bam`

Per-chunk outputs under `work/<sample>/align/`:

| File | Description |
|------|-------------|
| `<chunk>.forward_1_bismark_bt2_pe_sortbyname.bam` | forward-strand name-sorted BAM |
| `<chunk>.reverse_1_bismark_bt2_pe_sortbyname.bam` | reverse-strand name-sorted BAM |

Contract:

- Uses `samtools sort -n -@ <sort_threads> -o <out> <in>` (aligned with SeekSoulMethyl `SORT_BAM_BY_NAME` in `step2.nf`).
- Default `sort_threads=6`; unsorted source BAMs are retained.
- If a sortbyname output exists and is newer than its input, that BAM is skipped on re-run.
- BAM indexing, UMI dedup, and per-cell splitting are out of scope for this stage.

### `count_mapped_reads` (methylation-only path only)

Purpose: count aligned reads per cell barcode from unsorted Bismark BAMs (SeekSoulMethyl `COUNTS_MAPPED_READS`).

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe.bam`

Per-strand outputs under `work/<sample>/align/`:

| File | Columns | Description |
|------|---------|-------------|
| `<chunk>.{forward,reverse}_1_bismark_bt2_pe_cb_aligned_reads_counts.csv` | `barcode`, `aligned_reads` | counts from BAM `CB:Z:` tag |

Contract:

- Uses unsorted BAMs (not sortbyname outputs).
- If a counts CSV exists and is newer than its input BAM, that BAM is skipped on re-run.
- Skipped when workflow uses `gexcb` mode.

### `estimated_cells` (methylation-only path only)

Purpose: merge per-chunk/strand barcode counts and filter to called cells (SeekSoulMethyl `ESTIMATED_CELLS`).

Inputs:

- `work/<sample>/align/*_cb_aligned_reads_counts.csv`

Sample-level outputs under `work/<sample>/cells/`:

| File | Description |
|------|-------------|
| `merged_barcode_counts.csv` | all barcodes with summed `aligned_reads` |
| `filtered_barcode_read_counts.csv` | barcodes passing threshold (`aligned_reads`, `barcode`) |
| `filtered_barcode` | one barcode per line |

Contract:

- Default `expected_cell_num=3000`; 99th-percentile index × 0.1 read threshold (aligned with SeekSoulMethyl `step3_estimated_cells.py`).
- Optional `force_cell_num`: when set, take top N barcodes by `aligned_reads` among barcodes with reads > 0 (aligned with SeekSoulMethyl `cell_identify.R` force-cell mode). `force_cell_num` overrides `expected_cell_num` filtering; both keys may appear in workflow JSON.
- Skipped when workflow uses `gexcb` mode.

### `split_bams`

Purpose: split name-sorted Bismark BAMs into per-cell BAM files (SeekSoulMethyl `SPLIT_BAM_FILES`).

Inputs:

- `work/<sample>/align/<chunk>.forward_1_bismark_bt2_pe_sortbyname.bam`
- `work/<sample>/align/<chunk>.reverse_1_bismark_bt2_pe_sortbyname.bam`
- barcode list (**one of**):
  - `work/<sample>/cells/filtered_barcode` (methylation-only path)
  - `gexcb` workflow path (RNA filtered barcodes file)

Per-strand outputs under `work/<sample>/split_bams/`:

| File / dir | Description |
|------------|-------------|
| `<chunk>.forward_1/<barcode>.bam` | forward-strand single-cell BAMs |
| `<chunk>.reverse_1/<barcode>.bam` | reverse-strand single-cell BAMs |
| `<chunk>.{forward,reverse}_1/<chunk>.{forward,reverse}_1_filtered_barcode` | barcodes with reads > 0 |
| `<chunk>.{forward,reverse}_1/<chunk>.{forward,reverse}_1_filtered_barcode_reads_counts.csv` | per-barcode read counts |

Contract:

- Input BAM must be sorted by read name; groups reads by `{CB}` prefix of QNAME (`qname.split("_")[0]`).
- Default `split_bams_cores=8` for parallel batch splitting.
- allcools and UMI dedup are out of scope for this stage.

### `merge_fr_bams`

Purpose: merge forward and reverse per-cell split BAMs into one BAM per barcode (SeekSoulMethyl `MERGE_BISMARK_BAM`).

Inputs (per chunk):

- `work/<sample>/split_bams/<chunk>.forward_1/*.bam`
- `work/<sample>/split_bams/<chunk>.reverse_1/*.bam`
- `*_filtered_barcode` and `*_filtered_barcode_reads_counts.csv` under both strand directories

Outputs under `work/<sample>/split_bams/merged/`:

| File / dir | Description |
|------------|-------------|
| `<chunk>_merged_fr_bam/<barcode>.bam` | merged single-cell BAM (e.g. `0001_merged_fr_bam`) |
| `<chunk>_merge_filtered_barcode` | union of F/R filtered barcodes |
| `<chunk>_merge_filtered_barcode_reads_counts.csv` | `reads_counts,barcode`; summed across strands |

Contract:

- Both strands present: `samtools merge -n -@ <threads> -o <out> <f> <r>`.
- Single strand only: copy input BAM to output path.
- Default `merge_fr_bams_cores=8`; one samtools thread per barcode merge, parallelized by process pool.
- Skip output BAM when existing file passes `samtools quickcheck`.
- Sample-wide union of filtered barcodes or merged BAMs across analysis chunks (barcodes are disjoint by prefix; see Chunk model), allcools conversion, and UMI dedup are out of scope for this stage.

### `bam_to_allc`

Purpose: convert merged per-cell BAMs to ALLC format via seekgene ALLCools `bam-to-allc` (SeekSoulMethyl `ALLCOOLS_BAM_TO_ALLC` / `step3_bam_to_allc.py`).

Inputs (per chunk):

- `work/<sample>/split_bams/merged/<chunk>_merged_fr_bam/*.bam`
- `work/<sample>/split_bams/merged/<chunk>_merge_filtered_barcode`
- `genome_fa`: reference genome FASTA (workflow JSON)
- `chrom_size_path`: chromosome sizes BED (workflow JSON; required for workflow parity; not passed to bam-to-allc CLI)

Outputs under `work/<sample>/allcools/`:

| File / dir | Description |
|------------|-------------|
| `<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz` | per-cell ALLC (gzip) |
| `<chunk>_merged_fr_bam_allcools/<barcode>_allc.count.csv` | optional count sidecar from `--save_count_df` |

Contract:

- Per barcode: `samtools sort` → `samtools index` → `allcools bam-to-allc --reference_fasta <genome_fa> --convert_bam_strandness --tag UR --save_count_df`; intermediate sorted BAM removed on success.
- Only barcodes listed in `<chunk>_merge_filtered_barcode` are processed.
- Default `bam_to_allc_cores=8`, `allcools_tag=UR`, `align_method=bismark`.
- Skip barcode when `<barcode>_allc.gz` exists and is non-empty.
- `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` during parallel workers.
- Requires seekgene ALLCools fork (`pixi run setup-allcools`) and `tabix` on PATH (`htslib` conda dep; `check-allcools-env` verifies). `bam_to_allc.py` prepends the pixi `bin` directory to subprocess `PATH` so ALLCools can invoke `tabix` on Slurm compute nodes without pixi activation (same pattern as `bismark_align` / bowtie2).
- `generate-dataset`, merge/extract allc, and sample-wide per-cell QC tables that glob `*_allc.gz.count.csv` from all analysis-chunk directories (gather-only; no barcode dedup across chunks) are out of scope for this stage.

### `saturation`

Purpose: estimate sample-level genome-fraction saturation curve via subsampling simulation on pre-dedup per-cell BAM per-base depth histograms (dbit-matrix CpG-coverage style, adapted for DD-MET5 single-cell breadth fraction).

Inputs:

- `work/<sample>/shard_fastq/fastp.json` — fastp report; `summary.after_filtering.total_bases` (fallback `before_filtering`) supplies plot x-axis sequencing depth in Gbp
- `work/<sample>/split_bams/merged/<chunk>_merged_fr_bam/<barcode>.bam` — pre-dedup per-cell BAMs
- Cell read-count table (HQ filter):
  - methylation-only: `work/<sample>/cells/filtered_barcode_read_counts.csv` (`aligned_reads`, `barcode`)
  - gexcb: `work/<sample>/split_bams/merged/*_merge_filtered_barcode_reads_counts.csv` (`reads_counts`, `barcode`; summed across chunks)
- `chrom_size_path`: chromosome sizes BED; summed to genome size `G` for normalization

Workflow parameters:

- `saturation_reads_threshold` (default `100`): minimum aligned reads for HQ cell inclusion
- `saturation_max_cells` (default `100`): maximum HQ cells used for curve estimation
- `saturation_sample_seed` (default `42`): random seed when sampling HQ cells above `saturation_max_cells`
- `saturation_linear_r2_threshold` (default `0.99`): linear-fit R² at/above which the linear extrapolation is used instead of the saturation curve

Outputs under `work/<sample>/qc/saturation/`:

| File | Description |
|------|-------------|
| `saturation_curve.png` | Observed median genome fraction (IQR error bars) + fitted line/curve + 2× prediction; x-axis = sequencing depth (Gbp) from fastp `total_bases` |
| `saturation_summary.tsv` | One-row TSV: `sample_id`, `observed_median_genome_fraction`, `theoretical_max_median_genome_fraction`, `predicted_median_genome_fraction_at_2x`, `saturation_rate`, `extrapolation_model`, `hq_cell_count`, `sampled_cell_count`, `sample_seed` |

Contract:

- Build per-cell per-base depth histogram from aligned primary reads across **all chunks** for each barcode; do not use deduplicated ALLC (`cov` is post-UMI-dedup and unsuitable for subsampling simulation).
- Depth is counted per **sequencing fragment** (read pair): PE mates share a query name and contribute at most once to each reference position's depth (overlapping mates are de-duplicated), so the subsampling unit matches a real fragment.
- Genome fraction at subsample depth `f`: `sum_d count(d) × (1 − (1−f)^d) / G`, where `count(d)` is the number of reference positions with fragment depth `d`.
- HQ cells: `aligned_reads > saturation_reads_threshold` and per-cell BAM present. If HQ count exceeds `saturation_max_cells`, randomly sample `saturation_max_cells` cells (fixed `saturation_sample_seed`; HQ list sorted before sampling).
- For each subsample fraction (1%–100%), compute per-cell genome fraction; aggregate with **median** and **IQR** (asymmetric `Q1`/`Q3` error bars) across sampled cells.
- Extrapolation (`f` beyond 1×): fit both a through-origin linear model `y = m × f` and a saturation curve `y = a × (1 − exp(−b × f))` to the median curve.
  - If the linear fit `R² ≥ saturation_linear_r2_threshold` (data essentially unsaturated), report `extrapolation_model = linear`, `predicted@2× = m × 2`, and `theoretical_max = saturation_rate = NA` (the asymptote is unidentifiable in the linear regime).
  - Otherwise report `extrapolation_model = saturation`, `theoretical_max = a`, `predicted@2× = a × (1 − exp(−2b))`, and `saturation_rate = median@100% / a × 100`.
- Plot x-axis: subsample fraction × fastp `total_bases` / 1e9 (Gbp); 2× prediction at twice the observed sequencing depth.
- Y-axis semantics: median **genome fraction** (0–1 in TSV; plot as %).
- Single sample-level job (no chunking); supports `--dry-run`.
- Does not depend on `merge_sc_metrics` or ALLCools count sidecars.

### `qc_summary`

Purpose: gather per-cell and sample-level QC metrics into summary tables (SeekSoulMethyl `step3_merge_sc_metrics` + `step4_wgs_summary` alignment).

Inputs:

- `work/<sample>/shard_fastq/fastp.json`
- `work/<sample>/demux/*.stats.json`, `work/<sample>/demux/qc.CtoT.tsv`
- `work/<sample>/align/*_{forward,reverse}_1_bismark_bt2_PE_report.txt`
- `work/<sample>/qc/saturation/saturation_summary.tsv`
- `work/<sample>/allcools/*_merged_fr_bam_allcools/*_allc.gz.count.csv`
- Cell read-count table:
  - methylation-only: `work/<sample>/cells/filtered_barcode_read_counts.csv`
  - gexcb: `work/<sample>/split_bams/merged/*_merge_filtered_barcode_reads_counts.csv` (summed across chunks)
- Optional workflow key `cbcsv`: path to `m_cb,gex_cb` map for `gex_cb` column (not the `gexcb` barcode list)

Outputs under `work/<sample>/summary/`:

| File | Description |
|------|-------------|
| `cells_summary.tsv` | One row per cell with ALLCools metrics, aligned reads, per-barcode `CtoT`, and `{context}_mc_rate` columns |
| `sample_summary.tsv` | One-row internal QC summary (selected metric subset) |
| `wgs_summary.csv` | SeekSoul-compatible one-row CSV; max-cell and bulk-CpG columns are `NA` |

`sample_summary.tsv` columns (selected subset):

- fastp: `raw_reads`, `total_bases`, `duplication_rate`
- demux: `valid_barcode_rate`, `barcode_corrected_fraction`, `dropped_too_short`, `dropped_chimeric`, `forward_reads`, `reverse_reads`, `rate_17lme`, `CtoT`, plus `rate_7f`, `rate_7f17lme`, `cc_mean` as `NA` until demux stats extended
- bismark: `mapped_to_genome`, `confidently_mapped`, `cpg_methylation_rate`, `chg_methylation_rate`, `chh_methylation_rate` (no `unknown_methylation_rate`)
- cells: `estimated_cells`, `reads_in_cells`, `fraction_reads_in_cells`
- saturation: `observed_median_genome_fraction`, `theoretical_max_median_genome_fraction`, `saturation_rate`, `extrapolation_model`, `sampled_cell_count`, `sample_seed` (no `hq_cell_count`, no `predicted_median_genome_fraction_at_2x`)
- median called-cell QC: `median_genome_cov`, `median_total_cpg_number`, `median_aligned_reads`, `median_cell_saturation`

`cells_summary.tsv` core columns: `cell_barcode`, `aligned_reads`, `CtoT`, `total_cpg_number`, `genome_cov`, `genome_cov_raw_umi`, `genome_cov_new_umi`, `cell_saturation`, `weighted_mc_rate`, optional `gex_cb`; plus `{context}_mc_rate` per ALLCools context. Omits `total_mc`, `total_cov`, `total_number`, and `{context}_mc/cov/number`.

Contract:

- Glob `allcools/*_merged_fr_bam_allcools/*_allc.gz.count.csv` across analysis chunks; gather-only (barcodes are disjoint by prefix; see Chunk model).
- Median-cell stats: among called cells with ALLCools metrics, sort by `genome_cov` descending and take the median index.
- `fraction_reads_in_cells` = `reads_in_cells` / pooled Bismark unique aligned reads.
- Single sample-level job; supports `--dry-run`.
- Does not produce per-cell JSON, HTML reports, or bulk merged ALLC metrics.
