# Stage contracts

Normative input/output contracts for seeksoul-matrix workflow stages.

## Main stage order (planned)

`fastp_split -> demux_extract_bc -> bismark_align -> bam_sort -> (count_mapped_reads -> estimated_cells)? -> split_bams -> merge_fr_bams -> bam_to_allc -> ...`

`fastp_split`, `demux_extract_bc`, `bismark_align`, `bam_sort`, `count_mapped_reads`, `estimated_cells`, `split_bams`, `merge_fr_bams`, and `bam_to_allc` are implemented in this repository revision.

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
- Chunk count is controlled by `number_of_split_parts` (passed to `fastp --split`).
- Downstream `demux_extract_bc` will read paired shards from `shard_fastq/`.

### `demux_extract_bc`

Purpose: DD-MET5 barcode extraction, C→T QC, adapter trimming, and forward/reverse paired FASTQ output.

Inputs:

- `work/<sample>/shard_fastq/R1*.fq.gz` or `<chunk>.R1.fq.gz` and paired R2 from `fastp_split`
- cell barcode whitelist: `whitelist/DD-MET5/U3CB_methylation.txt.gz` (default)

Per-chunk outputs under `work/<sample>/demux/`:

| File | Columns / fields | Description |
|------|------------------|-------------|
| `<chunk>.forward_1.fq.gz` / `.forward_2.fq.gz` | — | forward-strand paired FASTQ |
| `<chunk>.reverse_1.fq.gz` / `.reverse_2.fq.gz` | — | reverse-strand paired FASTQ |
| `<chunk>.linker.tsv` | `CR`, `UB`, `C`, `T` | one row per UMI-deduped QC read; `CR` = corrected CB, `UB` = 12 bp UMI, `C`/`T` = convertible-base counts |
| `<chunk>.stats.json` | grouped JSON | chunk-level demux and C→T QC metrics (see schema below) |

Sample-level output (after all chunks complete):

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

- Uses seekgene Bismark (stock bioconda Bismark lacks `--add_barcode` / `--add_umi`).
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
- Cross-chunk barcode consolidation, allcools, and UMI dedup are out of scope for this stage.

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
- Requires seekgene ALLCools fork (`pixi run setup-allcools`).
- `generate-dataset`, merge/extract allc, and cross-chunk metric consolidation are out of scope for this stage.
