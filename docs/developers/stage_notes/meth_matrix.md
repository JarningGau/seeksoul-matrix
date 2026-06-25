# meth_matrix

Implementation notes for [`meth_matrix`](../contracts.md#meth_matrix). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `matrix` / `matrix_sparse`: aggregate per-cell methylation statistics over BED regions.

1. Resolve BED (`--regions-bed` or fallback `meth/vmr/vmrs.bed`).
2. Per chromosome with CSR data: load `{chrom}.npz` + `smoothed/{chrom}.csv.gz`.
3. Numba kernel `_calc_mean_mfracs`: per region × cell → methylated count, total count, shrunken-residual sum.
4. **Sparse (default):** emit non-zero coverage pairs to `matrix.mtx.gz`; finalize `features.tsv.gz` + `barcodes.tsv.gz`.
5. **Dense (`--dense`):** write four cell × region `.csv.gz` tables.

Region names: `{chrom}:{start}-{end}` using BED columns 1–3 (additional columns ignored).

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `run_meth_analysis` | `false` | prerequisite for meth stages |
| `run_meth_matrix` | `false` | append `meth_matrix` after `meth_scan` |
| `meth_regions_bed` | `""` | explicit BED; else `meth/vmr/vmrs.bed` |
| `meth_regions_label` | `""` | output dir under `meth/regions/` |
| `meth_matrix_dense` | `false` | dense four-table output |
| `meth_matrix_cores` | `8` | numba parallel threads |

## CLI flags (`scripts/meth_matrix.py`)

| Flag | Role |
|------|------|
| `--work-path` | matrix at `<work_path>/meth/matrix/` |
| `--data-dir` | matrix override (requires `--regions-bed` + `--output-dir`) |
| `--regions-bed` | BED path (optional with `--work-path`) |
| `--output-dir` | default `meth/regions/<label>/` |
| `--dense` | opt-in dense output |
| `--threads` | parallel threads |
| `--dry-run` | print resolved paths and exit |

## Skipped stage: `meth_matrix_filter`

MethSCAn `filter` (column subset / QC thresholds) is **not** implemented. `allc_to_matrix` already restricts cells via `cells/filtered_barcode` (methylation-only) or gexcb merge lists.

## Toolchain

Root pixi: `numpy`, `scipy`, `pandas`, `numba`.

Shared library: `scripts/lib/meth_matrix/region_matrix.py`.

## Validation

MethSCAn `matrix --sparse` parity **passed** on `work/dd-met5-example` (50 cells, 2 VMRs from `vmrs.bed`): 19 non-zero entries; exact match on `(row_i, col_i, mfrac)` vs v1.1.0 reference (MethSCAn requires uncompressed `smoothed/*.csv` for its CLI check; parity run used decompressed sidecar).

## License note

MethSCAn algorithms reimplemented clean-room with citation. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).
