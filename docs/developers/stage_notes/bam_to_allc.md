# bam_to_allc

Implementation notes for [`bam_to_allc`](../contracts.md#bam_to_allc). Not normative for I/O.

## Behavior

Per barcode pipeline:

1. `samtools sort`
2. `samtools index`
3. `allcools bam-to-allc --reference_fasta <genome_fa> --convert_bam_strandness --tag UR --save_count_df`

Intermediate sorted BAM removed on success. Only barcodes listed in `<chunk>_merge_filtered_barcode` are processed.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `bam_to_allc_cores` | `8` | parallel workers |
| `allcools_tag` | `UR` | UMI tag for dedup |
| `align_method` | `bismark` | ALLCools align method |
| `genome_fa` | (required) | reference FASTA |
| `chrom_size_path` | (required in workflow) | chromosome sizes BED; workflow parity only, not passed to bam-to-allc CLI |

## Skip / incremental re-run

Skip barcode when `<barcode>_allc.gz` exists and is non-empty.

`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` during parallel workers.

## Toolchain and environment

- Requires seekgene ALLCools fork (`pixi run setup-allcools`) and `tabix` on PATH (`pixi run check-allcools-env`).
- `bam_to_allc.py` prepends the pixi `bin` directory to subprocess `PATH` so ALLCools can invoke `tabix` on Slurm compute nodes without pixi activation (same pattern as `bismark_align` / bowtie2).

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `ALLCOOLS_BAM_TO_ALLC` / `step3_bam_to_allc.py`.

## Out of scope

`generate-dataset`, merge/extract allc, and sample-wide per-cell QC tables that glob `*_allc.gz.count.csv` from all analysis-chunk directories (gather-only; no barcode dedup across chunks).
