# regroup_shards

Implementation notes for [`regroup_shards`](../contracts.md#regroup_shards). Not normative for I/O.

## Behavior

- Gzip member concatenation (binary append) merges read-order sub-shards per stream into one analysis shard per barcode prefix.
- Analysis `chunk_id` equals barcode prefix (e.g. `A`, `G`, `T` when `split_fastq_prefix_bases=1`).
- Writes `chunks.tsv` manifest (see [`chunk_model.md`](../chunk_model.md)).

## Defaults and workflow keys

Inherits `split_fastq_prefix_bases` from demux workflow config.

## Skip / incremental re-run

Re-run regenerates regrouped FASTQs and manifest when inputs change.

## Toolchain and environment

`scripts/regroup_shards.py`; binary gzip append (no recompression).

## SeekSoulMethyl / dbit-matrix alignment

Aligned with SeekSoul `split_fastq` regroup semantics.

## Out of scope

Downstream discovery uses top-level `demux/<prefix>.*` only, not `demux/shards/`.
