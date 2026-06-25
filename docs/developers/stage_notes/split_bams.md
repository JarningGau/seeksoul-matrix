# split_bams

Implementation notes for [`split_bams`](../contracts.md#split_bams). Not normative for I/O.

## Behavior

- Input BAM must be sorted by read name.
- Groups reads by cell barcode: `{CB}` prefix of QNAME (`qname.split("_")[0]`).
- Barcode list source depends on workflow mode — see [`chunk_model.md`](../chunk_model.md#barcode-selection).

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `split_bams_cores` | `8` | parallel batch splitting |

## Skip / incremental re-run

Per analysis chunk and strand.

## Toolchain and environment

`scripts/split_bams.py`; process pool for parallel batches.

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoulMethyl `SPLIT_BAM_FILES`.

## Out of scope

allcools conversion and UMI dedup.
