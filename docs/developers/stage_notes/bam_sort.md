# bam_sort

Implementation notes for [`bam_sort`](../contracts.md#bam_sort). Not normative for I/O.

## Behavior

Uses `samtools sort -n -@ <sort_threads> -o <out> <in>` (aligned with SeekSoulMethyl `SORT_BAM_BY_NAME` in `step2.nf`). Unsorted source BAMs are retained.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `sort_threads` | `6` | `samtools sort -@` threads |

## Skip / incremental re-run

If a sortbyname output exists and is newer than its input, that BAM is skipped on re-run.

## Toolchain and environment

`samtools` from root pixi environment.

## SeekSoulMethyl / dbit-matrix alignment

Name-sort step aligned with SeekSoulMethyl `step2.nf`.

## Out of scope

BAM indexing, UMI dedup, and per-cell splitting.
