# MethSCAn-native analysis — implementation spec

**Status:** Phase 1 in progress — `allc_to_matrix` implemented; MethSCAn `prepare` parity passed; `meth_smooth` / `meth_scan` pending.  
**Last updated:** 2026-06-25

## Summary

Integrate [MethSCAn](https://anders-biostat.github.io/MethSCAn/) single-cell methylation **analysis** capabilities into seeksoul-matrix as first-party `scripts/` stages, without a runtime dependency on the `methscan` PyPI package and without supporting MethSCAn’s generic per-tool input adapters (Bismark `.cov`, methylpy, biscuit, custom `--input-format` strings).

The pipeline already produces per-cell ALLC via `bam_to_allc`. New stages read **only** seeksoul-matrix contract paths and ALLC column semantics. Internal sparse-matrix layout and downstream algorithms follow MethSCAn’s published methods; implementation is a clean-room port aligned with dbit-matrix stage conventions.

When this work ships, update [`docs/developers/contracts.md`](developers/contracts.md), add matching [`docs/developers/stage_notes/`](developers/stage_notes/) pages, and bump [`docs/developers/status.md`](developers/status.md) per [`doc-system.md`](developers/doc-system.md).

---

## Goals

| Goal | Detail |
|------|--------|
| **No `methscan` package** | Do not add MethSCAn to `pixi.toml` / workflow env; no subprocess calls to `methscan`. |
| **No foreign input formats** | Single reader: gzipped ALLC from `bam_to_allc` (`<barcode>_allc.gz`). No Bismark `.cov` path. |
| **Pipeline-native I/O** | Inputs from [`contracts.md`](developers/contracts.md) (`allcools/`, `cells/`, optional `summary/`); outputs under a new `work/<sample>/meth/` tree. |
| **Thin stage scripts** | One script per logical step; JSON workflow keys; `make_cmd.py` driver; `--dry-run` on every stage. |
| **Method fidelity** | VMR scan, region matrix, smoothing, and (later) DMR/profile behavior should match MethSCAn defaults unless explicitly documented. |
| **Cite MethSCAn** | Document Kremer et al., *Nature Methods* 2024 ([doi:10.1038/s41592-024-02347-x](https://doi.org/10.1038/s41592-024-02347-x)) in user-facing docs and optional `--cite` on analysis scripts. |

## Non-goals (initial phases)

- Replacing `bam_to_allc` or ALLCools.
- Merged sample-wide ALLC matrix export ([`status.md`](developers/status.md) “generate-dataset” gap) — may remain separate.
- Supporting non–seeksoul-matrix ALLC variants or external barcode lists outside existing methylation-only / gexcb contracts.
- Interactive plotting inside the pipeline (MethSCAn `profile` produces CSV for external R/Python plotting).
- Default-on inclusion in the twelve-stage `fastp_split → qc_summary` path before validation; first delivery is **optional post-QC** stages.

---

## Background: what MethSCAn provides

Reference checkout: local `MethSCAn/` (gitignored; algorithm reference only). Upstream v1.1.0 commands:

| MethSCAn command | Role |
|------------------|------|
| `prepare` | Per-cell coverage files → per-chromosome CSR sparse matrices (`{chrom}.npz`), `column_header.txt`, `cell_stats.csv` |
| `filter` | Subset matrices by site count / global methylation or explicit cell list |
| `smooth` | Tricube-weighted pseudobulk smoothing → `smoothed/{chrom}.csv.gz` |
| `scan` | Sliding-window variance on shrunken residuals → VMR BED |
| `diff` | Two-group DMRs with permutation FDR |
| `matrix` | Cell × region methylation tables (counts, fractions, shrunken residuals) |
| `profile` | Mean methylation profile around BED features |

MethSCAn’s `prepare` accepts many input formats; seeksoul-matrix will **skip that layer** and ingest ALLC directly.

### seeksoul-matrix ALLC (canonical input)

From `bam_to_allc` ([`stage_notes/bam_to_allc.md`](developers/stage_notes/bam_to_allc.md)):

```
work/<sample>/allcools/<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz
```

Tab-separated rows (ALLCools `bam-to-allc`):

| Col | Field | Example |
|-----|-------|---------|
| 1 | `chrom` | `chr1` |
| 2 | `pos` | `34170091` |
| 3 | `strand` | `+` / `-` |
| 4 | `context` | `CG`, `CHG`, `CHH`, … |
| 5 | `mc` | methylated read count |
| 6 | `cov` | total coverage |
| 7+ | (optional ALLCools fields) | ignored unless needed later |

**Site encoding for the sparse matrix** (MethSCAn-compatible semantics):

- For each genomic position and cell, after context filter: `+1` methylated, `-1` unmethylated, `0` missing / ambiguous.
- **Context filter (default):** `context == "CG"` (workflow key `meth_context`, default `CG`). Document whether CHG/CHH modes are in scope for v1.
- **Ambiguous sites:** `mc > 0` and `mc < cov` → discard unless `--round-sites` (same rule as MethSCAn `prepare`).
- **Ties:** `mc == cov - mc` → always discard.

Cell names: barcode string from filename (`<barcode>_allc.gz` → `<barcode>`). Barcodes are unique across analysis chunks ([`logs.md`](developers/logs.md) — zero overlap in validated runs).

### Barcode / cell selection

| Mode | Barcode list source |
|------|---------------------|
| methylation-only | `work/<sample>/cells/filtered_barcode` (+ optional read counts from `filtered_barcode_read_counts.csv`) |
| gexcb | Union of `work/<sample>/split_bams/merged/*_merge_filtered_barcode` |

Gather pattern: same as `saturation` / `qc_summary` — sample-level job globs `allcools/*_merged_fr_bam_allcools/*_allc.gz` and intersects with the active barcode list. Reject or warn on duplicate barcodes across chunks.

---

## Proposed architecture

### Stage decomposition

Map MethSCAn commands to seeksoul-matrix stages (names tentative):

```
bam_to_allc  →  [existing]
     ↓
allc_to_matrix     # prepare equivalent: ALLC → CSR store
     ↓
meth_matrix_filter # optional; may default to passthrough if cells/ already filtered
     ↓
meth_smooth        # pseudobulk smoothing
     ↓
meth_scan          # VMRs (requires smooth)
     ↓
meth_matrix        # region matrix (requires smooth; needs external BED)
```

**Phase 2 (optional):** `meth_diff`, `meth_profile` — same internal store, extra CLI inputs (cell groups CSV, regions BED).

Default workflow extension (not enabled until validated):

```
… → qc_summary → allc_to_matrix → meth_smooth → meth_scan
```

`meth_matrix_filter` and `meth_matrix` are on-demand or workflow-flagged (`run_meth_matrix: true`, `vmr_regions_bed: …`).

### Output layout (draft contract)

Under `work/<sample>/meth/`:

| Path | Producer | Description |
|------|----------|-------------|
| `matrix/` | `allc_to_matrix` | `{chrom}.npz` (CSR), `column_header.txt`, `cell_stats.csv`, `run_info.json` |
| `matrix_filtered/` | `meth_matrix_filter` | Filtered CSR store (same layout) |
| `matrix/smoothed/` or `matrix_filtered/smoothed/` | `meth_smooth` | `{chrom}.csv.gz` smoothed pseudobulk |
| `vmr/vmrs.bed` | `meth_scan` | VMR intervals |
| `regions/<label>/` | `meth_matrix` | `methylated_sites.csv.gz`, `total_sites.csv.gz`, `methylation_fractions.csv.gz`, `mean_shrunken_residuals.csv.gz` |
| `dmr/` | `meth_diff` | DMR BED (phase 2) |
| `profile/` | `meth_profile` | Profile CSV (phase 2) |

Use `matrix_filtered/` when filtering runs; downstream stages take `--data-dir` pointing at the active store. Exact paths become normative in `contracts.md` at implementation time.

### Shared library vs monolithic scripts

Follow dbit-matrix “thin stage script” pattern:

```
scripts/
  allc_to_matrix.py
  meth_matrix_filter.py
  meth_smooth.py
  meth_scan.py
  meth_matrix.py
  meth_diff.py          # phase 2
  meth_profile.py       # phase 2
  lib/
    meth_matrix/        # shared: allc_reader, csr_build, smooth, scan, numerics
      __init__.py
      allc.py
      store.py
      smooth.py
      scan.py
      numerics.py
      ...
```

- Stage scripts: argparse, path validation, `workflow_input_checks`, `--dry-run`.
- `lib/meth_matrix/`: numba-heavy kernels; unit-testable without full workflow.

### Dependencies (pixi)

MethSCAn uses: `numpy`, `scipy`, `pandas`, `numba`, `statsmodels` (DMR only). Add to root `pixi.toml` when implementing — prefer conda-forge pins consistent with Python 3.11. No `click` / `methscan` CLI stack required (argparse only, matching existing stages).

### License note

MethSCAn is **GPL-3.0-or-later**. This spec assumes a **clean-room reimplementation** of algorithms with citation, not copying MethSCAn source into the repo. If any GPL code is adapted verbatim, legal review and license compatibility with the project default must be resolved before merge.

---

## Algorithm checklist (port from MethSCAn reference)

Use `MethSCAn/methscan/*.py` as a behavioral spec; reimplement in `scripts/lib/meth_matrix/`.

- [x] **ALLC → COO chunks → CSR** (`prepare.py`): chromosome chunking (`chunksize`, default 10 Mbp), COO temp files, CSR `indptr` construction, `int8` data values.
- [x] **Cell stats** (`cell_stats.csv`): `n_obs`, `n_meth`, `global_meth_frac` per cell.
- [ ] **Filter** (`filter.py`): column subset on CSR; rewrite `column_header.txt` + stats.
- [ ] **Smooth** (`smooth.py`): tricube kernel, bandwidth default 1000 bp, optional `log1p(coverage)` weights.
- [ ] **Shrunken residuals** (`numerics.py`): `_calc_mean_shrunken_residuals` for windowed scan/matrix.
- [ ] **Scan** (`scan.py`): sliding window (default bw 2000, step 100), variance threshold (default 0.02), `min_cells` (default 6), optional `bridge_gaps`, parallel over chromosomes.
- [ ] **Matrix** (`matrix.py`): per-region counts / fractions / mean shrunken residuals; optional sparse `.mtx.gz` output.
- [ ] **Diff** (`diff.py`, phase 2): t-test windows, permutation FDR, two groups only.
- [ ] **Profile** (`profile.py`, phase 2): fixed-width profiles around sorted BED regions.

Document any intentional numeric deviations from MethSCAn in stage notes.

**Validated deviations (`allc_to_matrix`):**

- Default `meth_context=CG` filters ALLC trinucleotide context by prefix; MethSCAn `prepare` ingests all contexts.
- Cell barcodes strip `_allc` filename suffix; MethSCAn uses full basename (`<barcode>_allc`).
- ALLCools files have no header; MethSCAn `--input-format allc` skips the first row per file — use builtin reader or a no-header custom format.

---

## Workflow / `make_cmd.py` integration

- [x] `allc_to_matrix` in stage list; gated by `run_meth_analysis` (default `false`).
- [ ] `meth_smooth`, `meth_scan`, … stage names and `STAGE_REQUIRED_FIELDS`.
- [ ] Workflow keys (draft):

| Key | Default | Used by |
|-----|---------|---------|
| `run_meth_analysis` | `false` | gate optional stages |
| `meth_context` | `CG` | `allc_to_matrix` |
| `meth_chunksize` | `10000000` | `allc_to_matrix` |
| `meth_round_sites` | `false` | `allc_to_matrix` |
| `meth_smooth_bandwidth` | `1000` | `meth_smooth` |
| `meth_scan_bandwidth` | `2000` | `meth_scan` |
| `meth_scan_stepsize` | `100` | `meth_scan` |
| `meth_scan_var_threshold` | `0.02` | `meth_scan` |
| `meth_scan_min_cells` | `6` | `meth_scan` |
| `meth_matrix_cores` | `8` | parallel stages |
| `vmr_regions_bed` | — | `meth_matrix` (optional VMR-free regions BED) |

- [x] Pixi dry-run: `meth-allc-to-matrix-dry-run`
- [ ] Pixi dry-run: `meth-smooth-dry-run`, `meth-scan-dry-run`, `meth-e2e-dry-run`
- [ ] Slurm: single aggregate jobs for sample-wide stages (like `saturation`), not per analysis chunk.

---

## Implementation phases

### Phase 0 — Design lock (this document)

- [x] Capture goals, I/O, stage split, and open questions.
- [x] Phase-1 command set confirmed: **scan-only MVP** first (`allc_to_matrix` → `meth_smooth` → `meth_scan`); `meth_matrix` deferred to Phase 2.
- [ ] Review with stakeholders; resolve remaining open questions below.

### Phase 1 — Core store + VMR path (MVP)

1. [x] `scripts/lib/meth_matrix/` — ALLC reader + CSR build.
2. [x] `scripts/allc_to_matrix.py` — gather ALLC across chunks; write `work/<sample>/meth/matrix/`.
3. [ ] `scripts/meth_smooth.py`
4. [ ] `scripts/meth_scan.py`
5. [x] `make_cmd.py` wiring + `workflow/dd_met5_test.json` optional block (`allc_to_matrix` only).
6. [x] Docs: `contracts.md`, `stage_notes/allc_to_matrix.md`, `status.md`, `logs.md`; [ ] `stage_notes/meth_scan.md` (pending implementation).

**MVP validation (`prepare` / `allc_to_matrix`):** [x] MethSCAn `prepare` parity passed on `work/dd-met5-example` (50 cells, all-context). **Pending:** VMR comparison after `meth_smooth` + `meth_scan` ship.

### Phase 2 — Filter + region matrix

1. `scripts/meth_matrix_filter.py` — integrate with `cells/filtered_barcode` as default `--cell-names` source.
2. `scripts/meth_matrix.py` — user-supplied BED (promoters, VMRs from phase 1, etc.).
3. Workflow keys for filter thresholds vs passthrough.

### Phase 3 — DMR + profile

1. `scripts/meth_diff.py`, `scripts/meth_profile.py`
2. `statsmodels` dependency gated to these stages if desired.

### Phase 4 — HPC + production config

1. Slurm memory/CPU tiers in `dd_met5_slurm.json`
2. Cluster validation on production-scale cell count
3. README / examples (`examples/run_meth_analysis_example.sh`)

---

## Testing strategy

| Level | Approach |
|-------|----------|
| Unit | Synthetic ALLC fixtures (few cells, one chrom); CSR round-trip; ambiguous-site rules — **not yet added** |
| Golden | MethSCAn `prepare` parity on `dd-met5-example` — **passed** (one-off local check; not in repo) |
| Integration | `pixi run meth-allc-to-matrix-dry-run`; local run on `dd-met5-example` after `bam_to_allc` |
| Regression | Record validation outcomes in `logs.md` when CSR logic changes |

No automated test suite exists today ([`status.md`](developers/status.md)); first meth stages should add pytest under `tests/` following dbit-matrix regression style when implemented.

---

## Open questions

1. **Default context:** CG-only for v1, or configurable CHG/CHH for non-CpG assays?
2. **Filter duplication:** Skip `meth_matrix_filter` when `estimated_cells` already ran, or always expose MethSCAn-style site/meth thresholds for analysis QC?
3. **Stage naming:** Prefer `allc_to_matrix` vs `meth_prepare` for consistency with MethSCAn vocabulary?
4. **VMR regions for clustering:** Auto-run `meth_scan` → feed `vmr/vmrs.bed` into `meth_matrix`, or require explicit BED path?
5. **GPL:** Confirm clean-room policy with maintainers before copying numba kernels verbatim.
6. **gexcb:** Same meth stages for gexcb path with merged barcode gather — any RNA-specific QC joins needed?

---

## References

- MethSCAn documentation: <https://anders-biostat.github.io/MethSCAn/>
- MethSCAn commands (local): `MethSCAn/docs/commands.md`
- seeksoul-matrix stage contracts: [`docs/developers/contracts.md`](developers/contracts.md)
- Engineering patterns: `dbit-matrix/` (`make_cmd.py`, thin stages, JSON workflows)
- ALLC producer: [`docs/developers/stage_notes/bam_to_allc.md`](developers/stage_notes/bam_to_allc.md)
