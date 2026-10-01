# meth_profile

Implementation notes for [`meth_profile`](../contracts.md#meth_profile). Not normative for I/O.

## Behavior

Clean-room port of MethSCAn `profile`: average per-cell methylation around BED features. The stage writes a long CSV. It does not plot.

1. Read `column_header.txt` and the user BED (`--regions-bed` / `meth_profile_bed`). There is no `vmrs.bed` fallback.
2. For each region, center `c = (start + end) // 2` and take `width` bp starting at `c - width // 2`. Shorter features are extended; longer features are cut.
3. Walk CSR rows in that window. `+1` counts as methylated, `-1` as unmethylated. Both increment `n_total`.
4. Relative position is `pos - c`. A `-` strand (optional 1-indexed `--strand-column`) negates it so upstream stays negative. `.` is treated as `+`.
5. Sum counts across regions and chromosomes. `meth_frac = n_meth / n_total`. Write rows with `n_total > 0`.

The CSR row index is the BED coordinate, same as `meth_matrix`. The BED does not need to be sorted; regions are grouped by chromosome.

Stored relative positions are `[-width//2, width - width//2)`. On an even width, a minus-strand site at the genomic start of the window lands on `+width//2`, which is outside that axis, and that one base is dropped.

## Memory

Thread-local `int32` buffers are shape `(threads, n_cells, width)` for methylated and total counts. Peak size is about `threads * n_cells * width * 8` bytes, then summed into `int64` arrays of shape `(n_cells, width)`. Chunking over cells is not implemented.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `run_meth_analysis` | `false` | prerequisite for meth stages |
| `run_meth_profile` | `false` | append `meth_profile` after `meth_matrix` when that stage is also on |
| `meth_profile_bed` | `""` | required feature BED |
| `meth_profile_label` | `""` | output dir under `meth/profile/` |
| `meth_profile_width` | `4000` | profile width in bp |
| `meth_profile_strand_column` | unset | 1-indexed strand column |
| `meth_matrix_cores` | `8` | numba threads (`--threads`) |

The constant CSV `label` column is only the stage-script `--label` flag. It is not a workflow key.

## CLI flags (`scripts/meth_profile.py`)

| Flag | Role |
|------|------|
| `--work-path` | matrix at `<work_path>/meth/matrix/` |
| `--data-dir` | matrix override (requires `--output-dir`) |
| `--regions-bed` | feature BED (required) |
| `--output-dir` | default `meth/profile/<label>/` |
| `--regions-label` | subdirectory name |
| `--width` | profile width |
| `--strand-column` | 1-indexed strand column |
| `--label` | constant output column |
| `--threads` | parallel threads |
| `--dry-run` | print resolved paths and exit |

## Toolchain

Root pixi: `numpy`, `numba`. No `methscan` package.

Shared library: `scripts/lib/meth_matrix/profile.py`.

## License note

MethSCAn algorithms reimplemented clean-room with citation. Reference: Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)).
