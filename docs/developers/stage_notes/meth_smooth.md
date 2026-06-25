# meth_smooth

Implementation notes for [`meth_smooth`](../contracts.md#meth_smooth). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `smooth`: tricube-weighted smoothing of per-position methylation fractions across a CSR matrix store.

1. Load each `{chrom}.npz` CSR from `meth/matrix/`.
2. Per genomic row: `mfrac = n_methylated / n_observed` (rows with zero coverage are skipped).
3. For each CpG position, smooth `mfrac` over a tricube kernel of width `bandwidth` (default 1000 bp; half-bandwidth = `bandwidth // 2`).
4. Write `matrix/smoothed/{chrom}.csv.gz` with columns `pos,smoothed_mfrac`.
5. Emit `matrix/smoothed/run_info.json`.

Optional `--use-weights` applies `log1p(coverage)` weights per site (MethSCAn `--use-weights`; default off).

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `run_meth_analysis` | `false` | append `meth_smooth` after `allc_to_matrix` when true |
| `meth_smooth_bandwidth` | `1000` | tricube kernel width (bp) |
| `meth_smooth_use_weights` | `false` | weight sites by log1p(coverage) |

## CLI flags (`scripts/meth_smooth.py`)

| Flag | Role |
|------|------|
| `--work-path` | read `<work_path>/meth/matrix/` |
| `--data-dir` | explicit matrix directory (test override) |
| `--bandwidth` | smoothing bandwidth (bp) |
| `--use-weights` | enable log1p coverage weights |
| `--dry-run` | print resolved paths and exit |

## Toolchain

Root pixi environment: `numpy`, `scipy`, `numba` (no `methscan` package).

Shared library: `scripts/lib/meth_matrix/smooth.py`.

## Format note

MethSCAn writes uncompressed `smoothed/{chrom}.csv`; seeksoul-matrix uses `{chrom}.csv.gz`. Numeric parity compares decoded rows only.

## Validation

MethSCAn `smooth` parity **passed** on `work/dd-met5-example` (50 cells, CG): max abs diff `0.0` vs v1.1.0 on all main chromosomes.

## License note

MethSCAn algorithms reimplemented clean-room with citation. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).
