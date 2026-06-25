# meth_scan

Implementation notes for [`meth_scan`](../contracts.md#meth_scan). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `scan`: sliding-window VMR detection on shrunken residuals.

1. Load CSR `{chrom}.npz` and smoothed `{chrom}.csv.gz` per chromosome.
2. Sort chromosomes by CSR file size (largest first).
3. Slide a window of width `bandwidth` (default 2000 bp) with step `stepsize` (default 100 bp).
4. Per window: compute per-cell mean shrunken residuals, then `nanvar` across cells.
5. On the **largest chromosome**, set global variance threshold = `nanquantile(window_variances, 1 - var_threshold)` (default top 2%).
6. Merge windows above threshold into VMR intervals (`bridge_gaps` optional).
7. Re-filter VMRs by `min_cells` (default 6); write `meth/vmr/vmrs.bed` + `run_info.json`.

### BED columns (no header by default)

| Col | Field |
|-----|-------|
| 1 | chromosome |
| 2 | VMR_start |
| 3 | VMR_end |
| 4 | variance |
| 5 | n_sites (CpGs in region) |
| 6 | n_cells (cells with coverage) |

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `meth_scan_bandwidth` | `2000` | sliding-window width (bp) |
| `meth_scan_stepsize` | `100` | window step (bp) |
| `meth_scan_var_threshold` | `0.02` | top fraction merged as VMRs |
| `meth_scan_min_cells` | `6` | minimum cells with coverage per VMR |
| `meth_scan_bridge_gaps` | `0` | merge neighboring VMRs within gap (0 = off) |
| `meth_matrix_cores` | `8` | numba parallel threads |

## CLI flags (`scripts/meth_scan.py`)

| Flag | Role |
|------|------|
| `--work-path` | matrix at `<work_path>/meth/matrix/`; output `<work_path>/meth/vmr/vmrs.bed` |
| `--data-dir` | matrix override (requires `--output`) |
| `--output` | BED path override |
| `--bandwidth`, `--stepsize`, `--var-threshold`, `--min-cells`, `--bridge-gaps`, `--threads` | scan parameters |
| `--write-header` | optional BED header row |
| `--dry-run` | print resolved paths and exit |

## Toolchain

Root pixi: `numpy`, `scipy`, `numba`, `pandas`.

Shared library: `scripts/lib/meth_matrix/scan.py`, `numerics.py`.

## Validation

MethSCAn `scan` parity **passed** on `work/dd-met5-example` (50 cells, default params): 2 VMRs; exact BED coordinate and statistic match vs v1.1.0.

## License note

MethSCAn algorithms reimplemented clean-room with citation. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).
