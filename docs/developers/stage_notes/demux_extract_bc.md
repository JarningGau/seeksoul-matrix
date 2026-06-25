# demux_extract_bc

Implementation notes for [`demux_extract_bc`](../contracts.md#demux_extract_bc). Not normative for I/O.

## Behavior

- Structure `B17U12` on R1; CB whitelist with Hamming distance 1; ambiguous multi-match CBs are discarded (no multi rescue).
- C→T QC uses ME5 positions and 17L+ME validation aligned with SeekSoulMethyl `calculate_ct_conversion_rate`.
- Only reads with **TTT** insert signature (forward chain) contribute to `linker.tsv` and CtoT statistics; reverse (`CCC`) reads are excluded from conversion QC.
- Each `(CR, UB)` pair contributes at most one row to `linker.tsv`.
- **CH chimeric filtering** (`--filter-ch`, default `2`, workflow key `filter_ch`): after adapter trim and length check, drop read pairs where trimmed R1 or R2 exceeds the threshold for strand-specific CH patterns (SeekSoulMethyl `should_filter_read_ch_pattern`). `0` disables filtering. Does not affect `linker.tsv` / C→T QC (computed on raw R1 before trim).

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `split_fastq_prefix_bases` | `1` | first N bases of corrected barcode for prefix sub-shards (SeekSoul `split_fastq`) |
| `filter_ch` | `2` | CH chimeric filter threshold; `0` disables |

Default whitelist: `whitelist/DD-MET5/U3CB_methylation.txt.gz`.

## Skip / incremental re-run

Per read-order chunk; aggregate `qc.CtoT.tsv` runs after all chunks complete.

## Toolchain and environment

`cutadapt`, Python demux logic in `scripts/demux_extract_bc.py`; sample-level aggregation in `scripts/aggregate_ct_qc.py`.

## SeekSoulMethyl / dbit-matrix alignment

C→T QC and CH filtering aligned with SeekSoulMethyl demux steps; prefix sub-sharding aligned with SeekSoul `split_fastq`.

## Out of scope

Spike-in output and multi-barcode rescue.

Field definitions for `linker.tsv`, `qc.CtoT.tsv`, and `stats.json`: [`qc_metrics.md`](../qc_metrics.md).
