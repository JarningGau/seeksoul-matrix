# Stage contracts

Normative input/output contracts for seeksoul-matrix workflow stages.

## Main stage order (planned)

`fastp_split -> demux_extract_bc -> bismark_align -> ...`

`fastp_split`, `demux_extract_bc`, and `bismark_align` are implemented in this repository revision.

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
- Spike-in output, multi-barcode rescue, and CH chimeric filtering are out of scope.

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
   └─ funnel.barcode_passed.valid.total
      ├─ funnel.barcode_passed.valid.forward
      └─ funnel.barcode_passed.valid.reverse
```

| Section | Fields | Description |
|---------|--------|-------------|
| (root) | `chunk_id`, `chemistry`, `input_r1`, `input_r2` | metadata |
| `funnel` | nested object above | demux read retention funnel (all counts are reads) |
| `ct` | `num_17lme`, `rate_17lme`, `ct_reads`, `ct_umi_dedup`, `ct_convertible_bases`, `ct_converted_bases`, `CtoT` | C→T QC funnel (parallel to FASTQ output; TTT-insert only; fields in processing order) |

Accounting: `funnel.total = barcode_rejected + barcode_ambiguous + barcode_passed.total`; `barcode_passed.total = unknown_chain + too_short + valid.total` (when `unknown_chain` reads leave before length filter).

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
