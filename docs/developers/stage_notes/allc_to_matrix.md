# allc_to_matrix

Implementation notes for [`allc_to_matrix`](../contracts.md#allc_to_matrix). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `prepare`: per-cell ALLC → per-chromosome CSR sparse matrices.

1. Gather `allcools/*_merged_fr_bam_allcools/*_allc.gz` across analysis chunks.
2. Filter rows by context prefix (`meth_context`, default `CG`).
3. Encode sites: `mc > 0` → `+1`, `mc == 0` → `-1`; ambiguous sites (`0 < mc < cov`) discarded unless `meth_round_sites` (ties always discarded).
4. Write temporary COO chunks per `(chrom, pos // chunksize)`, convert to CSR `{chrom}.npz`.
5. Emit `column_header.txt`, `cell_stats.csv`, `run_info.json`.

ALLC `context` column (col 4) is the full trinucleotide from ALLCools (`CGA`, `CTT`, …). CG filtering uses **prefix match** (`context.startswith("CG")`), not equality.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `run_meth_analysis` | `false` | append `allc_to_matrix` after `qc_summary` in workflow driver |
| `meth_context` | `CG` | context prefix filter (`CG`, `CHG`, `CHH`, `CH`, `all`) |
| `meth_chunksize` | `10000000` | COO temp-file chunk size (bp) |
| `meth_round_sites` | `false` | round ambiguous sites to majority vote |
| `meth_main_chroms_only` | `false` | keep only chr1–19, chrX, chrY, chrM |
| `meth_exclude_contigs` | `""` | comma-separated contigs to skip |

## CLI flags (`scripts/allc_to_matrix.py`)

| Flag | Role |
|------|------|
| `--work-path` | production gather from `allcools/` |
| `--allc-dir` | flat `*_allc.gz` test override (requires `--output-dir`) |
| `--cell-names` | explicit barcode list file |
| `--barcode-mode` | `methylation_only` or `gexcb` barcode source |
| `--dry-run` | print resolved paths and exit |

## Barcode selection

Priority: `--cell-names` → contract barcode file (`cells/filtered_barcode` or gexcb merge lists) → all discovered ALLC barcodes with `warning=no_barcode_list_using_all`.

## Toolchain and environment

Root pixi environment: `numpy`, `scipy`, `pandas`, `numba` (no `methscan` PyPI package).

Shared library: `scripts/lib/meth_matrix/` (`allc.py`, `store.py`).

## License note

MethSCAn algorithms are reimplemented clean-room with citation; no MethSCAn GPL source is copied into seeksoul-matrix. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).

## Out of scope

- Foreign input formats (Bismark `.cov`, methylpy, biscuit).
- `meth_matrix_filter` (cell filtering done in `allc_to_matrix`).
- `meth_diff`, `meth_profile`, and other Phase 3 meth stages.

## Validation

MethSCAn `prepare` parity **passed** on `work/dd-met5-example` (50 cells, `meth_context=all`): CSR output matches reference. Default `meth_context=CG` is intentional pipeline scope (CpG-focused analysis).
