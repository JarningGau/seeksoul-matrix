# allc_to_matrix

Implementation notes for [`allc_to_matrix`](../contracts.md#allc_to_matrix). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `prepare`: per-cell ALLC → per-chromosome CSR sparse matrices.

1. Gather `allcools/*_merged_fr_bam_allcools/*_allc.gz` across analysis chunks.
2. Read cells in a process pool (`--threads`, workflow key `meth_matrix_cores`). Each worker filters rows by context prefix (`meth_context`, default `CG`) and returns in-memory per-chromosome arrays.
3. Encode sites: `mc > 0` → `+1`, `mc == 0` → `-1`; ambiguous sites (`0 < mc < cov`) discarded unless `meth_round_sites` (ties always discarded).
4. Assemble each chromosome as CSC (one column per cell) and convert with `tocsr()` to `{chrom}.npz`. No temporary COO files.
5. Emit `column_header.txt`, `cell_stats.csv`, `run_info.json`.

ALLC `context` column (col 4) is the full trinucleotide from ALLCools (`CGA`, `CTT`, …). CG filtering uses **prefix match** (`context.startswith("CG")`), not equality.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `run_meth_analysis` | `false` | append `allc_to_matrix` after `qc_summary` in workflow driver |
| `meth_context` | `CG` | context prefix filter (`CG`, `CHG`, `CHH`, `CH`, `all`) |
| `meth_chunksize` | `10000000` | recorded in `run_info.json`; does not change the matrix or runtime |
| `meth_round_sites` | `false` | round ambiguous sites to majority vote |
| `meth_main_chroms_only` | `false` | keep only chr1–19, chrX, chrY, chrM |
| `meth_exclude_contigs` | `""` | comma-separated contigs to skip |
| `meth_matrix_cores` | `8` | worker processes for per-cell ALLC reads (`--threads`); also used by `meth_scan`, `meth_matrix`, and `meth_profile` |

## CLI flags (`scripts/allc_to_matrix.py`)

| Flag | Role |
|------|------|
| `--work-path` | production gather from `allcools/` |
| `--allc-dir` | flat `*_allc.gz` test override (requires `--output-dir`) |
| `--cell-names` | explicit barcode list file |
| `--barcode-mode` | `methylation_only` or `gexcb` barcode source |
| `--threads` | worker processes for per-cell ALLC reads (default 1) |
| `--dry-run` | print resolved paths and exit |

## Barcode selection

Priority: `--cell-names` → contract barcode file (`cells/filtered_barcode` or gexcb merge lists) → all discovered ALLC barcodes with `warning=no_barcode_list_using_all`.

## Toolchain and environment

Root pixi environment: `numpy`, `scipy`, `pandas` (no `methscan` PyPI package). This stage does not use numba.

Shared library: `scripts/lib/meth_matrix/` (`allc.py`, `store.py`).

Context prefilter: for a prefix filter other than `all` / `*`, lines are skipped in the gzip byte stream when they do not contain a tab plus that prefix (both cases). `context_matches` still decides the keep. `CG` on this chemistry keeps about 4% of ALLC rows.

Workers return arrays to the parent. There is no spill-to-disk path. A `CG` run of about 300 cells is a few hundred MB of site arrays. `meth_context=all` keeps most ALLC rows and is roughly 20× that (on the order of 8 GB for 300 cells) before the CSR `indptr` arrays, which are one `int64` per genomic position on each chromosome and dominate memory either way.

Slurm `cpus_per_task` for this stage should be at least `meth_matrix_cores`. `workflow/dd_met5_slurm_large.json` sets both to 16.

## License note

MethSCAn algorithms are reimplemented clean-room with citation; no MethSCAn GPL source is copied into seeksoul-matrix. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).

## Out of scope

- Foreign input formats (Bismark `.cov`, methylpy, biscuit).
- `meth_matrix_filter` (cell filtering done in `allc_to_matrix`).
- `meth_diff` (MethSCAn `diff`) — not planned; see [`methscan_builtin_spec.md`](../../methscan_builtin_spec.md).
- `meth_profile` is a separate optional stage; see [`meth_profile.md`](meth_profile.md).

## Validation

MethSCAn `prepare` parity **passed** on `work/dd-met5-example` (50 cells, `meth_context=all`): CSR output matches reference. That golden directory is not in the repo. The parallel rewrite was checked against the previous COO implementation, which is the code that passed that comparison: 60-cell `CG` matrices, headers, and `cell_stats.csv` match exactly (`--threads 1` and `--threads 8`); 2-cell `meth_context=all` and `round_sites` + `main_chroms_only` + `exclude chrM` also match. Default `meth_context=CG` is intentional pipeline scope (CpG-focused analysis).

On `work/C283_Brain_DNAme_S1` (300 cells, `CG`), wall time went from 613 s single-process to 134 s at 8 workers and 112 s at 16 workers. The chromosome CSR build stays single-process, so extra workers past that mostly wait on it.
