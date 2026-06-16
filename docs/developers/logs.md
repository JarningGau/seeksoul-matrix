# Development log

## 2026-06-16 — seekgene Bismark environment setup

**Task:** Add bowtie2/samtools/perl/git to pixi; install seekgene/Bismark fork with `--add_barcode` / `--add_umi`.

**Files changed:**
- `pixi.toml`, `pixi.lock`
- `scripts/install_seekgene_bismark.sh`, `scripts/check_bismark_env.sh`
- `.gitignore`

**Summary:**
- Added bioconda deps: `bowtie2`, `samtools`, `perl`, `git`.
- `pixi run setup-bismark` clones [seekgene/Bismark](https://github.com/seekgene/Bismark) at pinned commit `363ea7a` into `$CONDA_PREFIX/bin` (same pattern as SeekSoulMethyl).
- `pixi run check-bismark-env` verifies `bismark`, `bowtie2`, `samtools` on PATH and seekgene-only flags.

**Checks performed:**
- `pixi install`
- `pixi run setup-bismark`
- `pixi run check-bismark-env`

**Status:** done

**Notes:** Re-run `pixi run setup-bismark` after `pixi install` recreates `.pixi/envs/default` (Bismark scripts are not a conda package). Mouse Bismark index: `/mnt/wd-4t/resource/mouse-reference-GRCm39/fasta/` (`--genome` parent of `Bisulfite_Genome/`).

## 2026-06-15 — end-to-end workflow driver (`--stage all`)

**Task:** Consolidate `fastp_split` + `demux_extract_bc` with `run.sh` / `run.sbatch` generation; validate local e2e run and Slurm script generation.

**Files changed:**
- `scripts/make_cmd.py`
- `pixi.toml`

**Summary:**
- Added `--stage all` to generate per-stage scripts plus `work/<sample>/commands/run.sh` (local) or `run.sbatch` (Slurm DAG driver).
- Slurm demux: per-chunk `02_demux_extract_bc_<chunk>.sbatch`, `02_aggregate_ct_qc.sbatch`, and standalone `02_demux_extract_bc_submit.sh`; `run.sbatch` chains fastp → parallel demux chunks → aggregate.
- Added `pixi run e2e-dry-run` and `pixi run e2e-slurm-dry-run` helpers.

**Checks performed:**
- `pixi run e2e-dry-run`
- `pixi run e2e-slurm-dry-run`
- `pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all` (generated `run.sh`)
- `pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all --runner slurm` (generated `run.sbatch` + chunk sbatch files)
- `bash work/dd-met5-example/commands/run.sh` on `data/test_R1.fastq.gz` / `data/test_R2.fastq.gz`
  - fastp: 992780 reads/chunk pass filter; 2 shards
  - demux chunk 0001: valid=339315, CtoT=0.997; chunk 0002: valid=339753, CtoT=0.997
  - `qc.CtoT.tsv`: 23610 barcodes

**Status:** done

**Notes:** Slurm `run.sbatch` was generated and inspected; not submitted (no Slurm cluster in dev environment). Use `bash work/<sample>/commands/run.sbatch` on a Slurm login node to submit the DAG.

## 2026-06-15 — demux_extract_bc stage

**Task:** Implement DD-MET5 `demux_extract_bc` stage with linker.tsv / stats.json / qc.CtoT.tsv outputs.

**Files changed:**
- `scripts/demux_extract_bc.py`, `scripts/aggregate_ct_qc.py`
- `scripts/make_cmd.py`, `scripts/workflow_input_checks.py`
- `pixi.toml`, `pixi.lock`
- `workflow/dd_met5_test.json`
- `whitelist/DD-MET5/U3CB_methylation.txt.gz`
- `docs/developers/contracts.md`

**Summary:**
- Added per-chunk demux: B17U12 parse, HD=1 CB correction, UMI-deduped C→T QC, cutadapt trim, forward/reverse FASTQ.
- Per-chunk `linker.tsv` (CR, UB, C, T) and funnel-shaped `stats.json` (`funnel` + `ct` sections).
- Sample-level `qc.CtoT.tsv` via `aggregate_ct_qc.py` after all chunks.
- Extended `make_cmd.py` for local batch and Slurm per-chunk jobs plus aggregate submit helper.
- C→T QC restricted to TTT-insert (forward) reads per SeekSoulMethyl DD-MET5 spec.

**Checks performed:**
- `pixi install`
- `pixi run python scripts/demux_extract_bc.py --help`
- `pixi run python scripts/aggregate_ct_qc.py --help`
- `pixi run demux-dry-run`
- Real data: `work/dd-met5-example/shard_fastq/` (0001/0002 chunks, ~496k reads each)
  - chunk 0001: valid=339315, CtoT=0.997; chunk 0002: valid=339753, CtoT=0.997
  - `stats.json` grouped sections `reads`, `barcode`, `ct`; chunk 0001 CtoT=0.997, ct_umi_dedup=78381

**Status:** done

**Notes:** Shard discovery supports both `0001.R1.fq.gz` and fastp-native `R1.fq.gz` / `R1_001.fq.gz` naming. C→T QC includes only forward (TTT insert) reads per SeekSoulMethyl DD-MET5 spec.

## 2026-06-15 — fastp_split stage scaffold

**Task:** Implement first pipeline stage `fastp_split` with minimal workflow driver.

**Files changed:**
- `pixi.toml`, `pixi.lock`
- `scripts/fastp_split.py`, `scripts/make_cmd.py`, `scripts/_version.py`, `scripts/workflow_input_checks.py`, `scripts/__init__.py`
- `workflow/dd_met5_test.json`
- `docs/developers/contracts.md`

**Summary:**
- Added `fastp_split` stage script (QC + chunked output under `shard_fastq/`).
- Added minimal `make_cmd.py` driver for local and Slurm script generation.
- Documented stage I/O contract; example workflow JSON for DD-MET5.

**Checks performed:**
- `pixi install`
- `pixi run python scripts/fastp_split.py --help`
- `pixi run python scripts/make_cmd.py --version`
- `pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --skip-workdir-input-checks --dry-run`
- `pixi run python scripts/fastp_split.py --r1 /path/to/R1.fq.gz --r2 /path/to/R2.fq.gz --work-path /tmp/fastp_test --number-of-split-parts 2 --dry-run`
- `pixi run fastp-dry-run`
- Script generation without `--dry-run` (local `01_fastp_split.sh`)

**Status:** done

**Notes:** Workflow JSON uses placeholder FASTQ paths; real runs need existing `r1`/`r2` or omit `--skip-workdir-input-checks` only when inputs are present.
