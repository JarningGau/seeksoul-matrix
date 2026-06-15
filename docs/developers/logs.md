# Development log

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
