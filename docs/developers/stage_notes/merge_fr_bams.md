# merge_fr_bams

Implementation notes for [`merge_fr_bams`](../contracts.md#merge_fr_bams). Not normative for I/O.

## Behavior

- Both strands present: `samtools merge -n -@ <threads> -o <out> <f> <r>`.
- Single strand only: copy input BAM to output path.
- One samtools thread per barcode merge, parallelized by process pool.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `merge_fr_bams_cores` | `8` | parallel merge workers |

## Skip / incremental re-run

Skip output BAM when existing file passes `samtools quickcheck`.

## Toolchain and environment

`samtools merge` via `scripts/merge_fr_bams.py`.

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `MERGE_BISMARK_BAM`.

## Out of scope

Sample-wide union of filtered barcodes or merged BAMs across analysis chunks (barcodes are disjoint by prefix; see [`chunk_model.md`](../chunk_model.md)), allcools conversion, and UMI dedup.
