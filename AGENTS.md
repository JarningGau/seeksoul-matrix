# Repository Guidelines

Contributor guide for **seeksoul-matrix** — SeepSpace spatial methylation (DD-MET5, Slide-Tag–like chemistry). This repository is under active development.

Two **reference templates** sit alongside the production codebase — study them when building seeksoul-matrix; do not treat either as the production codebase or modify them unless explicitly updating the template:

- **`dbit-matrix/`** — engineering patterns (thin stage scripts, JSON workflows, `make_cmd.py` drivers, stage contracts, docs layout).
- **`SeekSoulMethyl/`** — SeekSoul's official single-cell methylation + RNA data-processing pipeline; domain reference for DD-MET3/DD-MET5 chemistry, barcode rules, stage logic, and analysis steps.

## Project Structure & Module Organization

```
seeksoul-matrix/
├── README.md              # DD-MET5 chemistry, barcode rules, whitelist path
├── pixi.toml              # Root workspace (linux-64)
├── whitelist/             # Cell-barcode whitelists (e.g. DD-MET5/U3CB_methylation.txt.gz)
├── scripts/               # Stage scripts and workflow driver (make_cmd.py)
├── workflow/              # JSON workflow configs
├── examples/              # Example local and HPC run commands (see run_*_example.sh)
├── docs/                  # Project documentation and stage contracts
├── dbit-matrix/           # Reference template — engineering patterns
└── SeekSoulMethyl/        # Reference template — SeekSoul official methylation pipeline
```

Barcodes for methylation data must not contain `C` (enzymatic C→T conversion). Fixed sequences (TSO, 17L, ME) retain `C` for conversion-rate QC — see `README.md`.

When adding seeksoul-matrix code, mirror **dbit-matrix** engineering conventions (thin stage scripts, JSON-driven workflows, `make_cmd.py` drivers, explicit stage contracts) and align stage behavior with **SeekSoulMethyl** where chemistry and analysis steps overlap.

## Build, Test, and Development Commands

**Root environment** (primary):

```bash
pixi install
```

**Reference templates** (read-only guidance; not the seeksoul-matrix runtime):

```bash
cd dbit-matrix && pixi install        # engineering patterns
# SeekSoulMethyl uses conda — see SeekSoulMethyl/README.md
```

Useful references:

- **seeksoul-matrix:** `docs/developers/doc-system.md` (doc layout and update rules), `docs/developers/contracts.md` (stable stage I/O), `docs/developers/chunk_model.md`, `docs/developers/qc_metrics.md`, `docs/developers/stage_notes/` (implementation details), `docs/developers/status.md` (current reliability map), `docs/developers/logs.md` (validation history)
- **dbit-matrix:** `docs/developers/architecture.md`, `contracts.md`, `AGENTS.md`
- **SeekSoulMethyl:** `README.md` (chemistry, installation), `nf/` (Nextflow workflow), `nf/bin/` (stage scripts), `docs/` (BAM dedup, reference genome)

Dry-run helpers (workflow driver):

```bash
pixi run fastp-dry-run         # fastp_split script generation
pixi run demux-dry-run         # demux_extract_bc script generation
pixi run regroup-dry-run       # regroup_shards script generation
pixi run bismark-align-dry-run # bismark_align script generation
pixi run bam-sort-dry-run      # bam_sort script generation
pixi run count-mapped-reads-dry-run # count_mapped_reads script generation
pixi run estimated-cells-dry-run   # estimated_cells script generation
pixi run split-bams-dry-run        # split_bams script generation
pixi run merge-fr-bams-dry-run     # merge_fr_bams script generation
pixi run bam-to-allc-dry-run       # bam_to_allc script generation
pixi run saturation-dry-run        # saturation script generation
pixi run qc-summary-dry-run        # qc_summary script generation
pixi run meth-allc-to-matrix-dry-run # allc_to_matrix script generation
pixi run meth-smooth-dry-run        # meth_smooth script generation
pixi run meth-scan-dry-run         # meth_scan script generation
pixi run meth-matrix-dry-run       # meth_matrix script generation
pixi run meth-e2e-dry-run          # --stage all with run_meth_analysis + run_meth_matrix
pixi run e2e-dry-run           # --stage all (local run.sh) dry-run
pixi run e2e-slurm-dry-run     # --stage all (Slurm run.sbatch) dry-run
```

Bismark environment (seekgene fork; required for `bismark_align`):

```bash
pixi run setup-bismark      # install seekgene/Bismark into pixi env
pixi run check-bismark-env  # verify bismark, bowtie2, samtools, seekgene flags
```

ALLCools environment (seekgene fork; required for `bam_to_allc`):

```bash
pixi run setup-allcools      # install seekgene/ALLCools into pixi env
pixi run check-allcools-env  # verify allcools, samtools, and tabix
```

End-to-end script generation and run:

```bash
pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all
bash work/<sample>/commands/run.sh
# Slurm: same with --runner slurm → work/<sample>/commands/run.sbatch
# HPC submit (login node): bash work/<sample>/commands/run.sbatch
```

Workflow JSON configs (`workflow/`):

| Config | Runner default | Barcode mode | Intended use |
|--------|----------------|--------------|--------------|
| `dd_met5_test.json` | `local` | methylation-only (`expected_cell_num`; optional `force_cell_num`) | Local runs, pixi dry-run tasks, CI-style validation |
| `dd_met5_slurm.json` | `slurm` | methylation-only | Production HPC (`number_of_split_parts=8`, `split_fastq_prefix_bases=2`); see `examples/run_slurm_example.sh` |
| `dd_met5_gexcb_test.json` | `local` | `gexcb` (RNA barcodes) | gexcb path testing; skips `count_mapped_reads` / `estimated_cells` |

## Pipeline stages

Implemented stages (`scripts/make_cmd.py`; stable I/O in `docs/developers/contracts.md`, details in `docs/developers/doc-system.md`):

| Stage | Scripts |
|-------|---------|
| `fastp_split` | `scripts/fastp_split.py` |
| `demux_extract_bc` | `scripts/demux_extract_bc.py`, `scripts/aggregate_ct_qc.py` |
| `regroup_shards` | `scripts/regroup_shards.py` |
| `bismark_align` | `scripts/bismark_align.py` |
| `bam_sort` | `scripts/bam_sort.py` |
| `count_mapped_reads` | `scripts/count_mapped_reads.py` |
| `estimated_cells` | `scripts/estimated_cells.py` |
| `split_bams` | `scripts/split_bams.py` |
| `merge_fr_bams` | `scripts/merge_fr_bams.py` |
| `bam_to_allc` | `scripts/bam_to_allc.py` |
| `saturation` | `scripts/saturation.py` |
| `qc_summary` | `scripts/qc_summary.py` |
| `allc_to_matrix` | `scripts/allc_to_matrix.py` |
| `meth_smooth` | `scripts/meth_smooth.py` |
| `meth_scan` | `scripts/meth_scan.py` |
| `meth_matrix` | `scripts/meth_matrix.py` |

Per-stage and workflow validation detail: [`docs/developers/status.md`](docs/developers/status.md).

`--stage all` generates per-stage scripts under `work/<sample>/commands/` plus a driver: `run.sh` (local) or `run.sbatch` (Slurm DAG). Analysis chunks are keyed by barcode prefix (`split_fastq_prefix_bases`, default `1`); `number_of_split_parts` controls read-order demux parallelism only. Barcode selection is **mutually exclusive**: `expected_cell_num` (default 3000, methylation-only path: count → estimate → split → merge → allc → saturation → qc_summary) or `gexcb` (RNA barcodes, split → merge → allc → saturation → qc_summary). Optional `force_cell_num` in workflow JSON takes top N barcodes by `aligned_reads` and overrides `expected_cell_num` threshold filtering in `estimated_cells`. Slurm emits per-chunk sbatch files for parallel stages and aggregate jobs for `estimated_cells` / `aggregate_ct_qc`.

Twelve-stage driver (`fastp_split` → `qc_summary`) with barcode-prefix analysis chunks; optional meth analysis adds up to four stages through `meth_matrix` (`run_meth_analysis`, `run_meth_matrix`). MethSCAn `diff` / `profile` are out of scope. See [`docs/developers/status.md`](docs/developers/status.md) for methylation-only and gexcb validation posture.

## Coding Style & Naming Conventions

Follow **dbit-matrix** engineering patterns when implementing seeksoul-matrix; consult **SeekSoulMethyl** for methylation-specific stage logic and chemistry details:

- **Python 3.11**; 4-space indentation; `snake_case` for modules and functions.
- Thin, explicit single-stage scripts; workflow parameters in `workflow/*.json`.
- Every stage must support `--dry-run`.
- New docs in **English** under `docs/`.
- Update `pixi.lock` when root dependencies change.

## Testing Guidelines

No automated test suite yet. For current per-stage and workflow validation posture, see [`docs/developers/status.md`](docs/developers/status.md); historical evidence in [`docs/developers/logs.md`](docs/developers/logs.md). When adding tests, follow the template's regression style in `dbit-matrix/docs/maintenance/`. Before finishing workflow changes, run `scripts/make_cmd.py --version`, per-stage `--help` and `--dry-run`, and the relevant `pixi run *-dry-run` tasks (`fastp-dry-run`, `demux-dry-run`, `regroup-dry-run`, `bismark-align-dry-run`, `bam-sort-dry-run`, `count-mapped-reads-dry-run`, `estimated-cells-dry-run`, `split-bams-dry-run`, `merge-fr-bams-dry-run`, `bam-to-allc-dry-run`, `saturation-dry-run`, `qc-summary-dry-run`, `meth-allc-to-matrix-dry-run`, `meth-smooth-dry-run`, `meth-scan-dry-run`, `meth-matrix-dry-run`, `meth-e2e-dry-run`, `e2e-dry-run`, `e2e-slurm-dry-run`).

## Lightweight Development Loop

Minimal loop for pipeline work — not a heavyweight engineering process.

### Before editing

Briefly declare the task:

- what will be changed
- why the change is needed
- which files are likely to be touched
- how the result will be checked

### During implementation

- make the smallest change that solves the task
- avoid unrelated refactoring
- preserve existing input/output formats unless explicitly required
- preserve existing command-line behavior when possible
- keep code and documentation consistent with the current project structure

### After editing

1. Run the most relevant lightweight check available — e.g. script on a small input, one representative pipeline step, expected output files, existing tests, or generated reports/tables.
2. Append `docs/developers/logs.md` with a short entry:

   - **date**
   - **task name**
   - **files changed**
   - **summary of changes**
   - **check performed**
   - **status:** `done`, `needs_review`, or `blocked`
   - **notes:** remaining risks, if any

3. If the change affects validation posture (new check, regression, limitation, or contract-sensitive behavior), update [`docs/developers/status.md`](docs/developers/status.md) in the same change.

**Validation rule:** if no check was run, status must be `needs_review` (not `done`); record why. A `logs.md` entry alone does not update the reliability map — update `status.md` when posture changes.

## Commit & Pull Request Guidelines

Use [Conventional Commits](https://www.conventionalcommits.org/): `type: short description` in imperative mood, lowercase, under 72 characters. Common types: `feat`, `fix`, `refactor`, `docs`, `chore`, `data`, `analysis`.

```
feat: add DD-MET5 whitelist validation
fix: reject barcodes containing C in demux
docs: document library structure in README
```

One commit per logical change. PRs should link issues and describe chemistry or workflow impact.

## Agent-Specific Instructions

- Build new functionality at the **repository root**, adapting engineering patterns from `dbit-matrix/` and domain logic from `SeekSoulMethyl/` — do not extend either template in place.
- Do not change stage input/output contracts silently; update the matching layer in `docs/developers/` per `doc-system.md` (contracts, chunk model, QC metrics, or stage notes).
- Prefer the simplest fix that preserves required behavior.
- Do not overwrite user changes without explicit permission.
