# count_mapped_reads

Implementation notes for [`count_mapped_reads`](../contracts.md#count_mapped_reads). Not normative for I/O.

## Behavior

- Counts aligned reads per cell from BAM `CB:Z:` tag.
- Uses **unsorted** Bismark BAMs (not sortbyname outputs).

## Defaults and workflow keys

None beyond per-chunk strand inputs.

## Skip / incremental re-run

If a counts CSV exists and is newer than its input BAM, that BAM is skipped on re-run.

## Toolchain and environment

`scripts/count_mapped_reads.py`; reads `CB:Z:` from unsorted BAMs.

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `COUNTS_MAPPED_READS`.

## Out of scope

Skipped when workflow uses `gexcb` mode (see [`chunk_model.md`](../chunk_model.md#barcode-selection)).
