# fastp_split

Implementation notes for [`fastp_split`](../contracts.md#fastp_split). Not normative for I/O.

## Behavior

- Adapter trimming is disabled (`--disable_adapter_trimming`).
- Chunk count is controlled by `number_of_split_parts` (passed to `fastp --split`); this controls **read-order parallelism** for demux only, not analysis chunk boundaries (see [`chunk_model.md`](../chunk_model.md)).

## Defaults and workflow keys

| Key | Role |
|-----|------|
| `number_of_split_parts` | read-order chunk count for `fastp --split` |

## Skip / incremental re-run

Not documented; re-run overwrites or regenerates shard FASTQs per script behavior.

## Toolchain and environment

Uses `fastp` from the root pixi environment.

## SeekSoulMethyl / dbit-matrix alignment

Read-order chunking pattern aligned with dbit-matrix `fastp_split`.

## Out of scope

Analysis-chunk (barcode-prefix) boundaries are determined in `demux_extract_bc` / `regroup_shards`.
