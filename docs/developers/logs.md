# Development log

## 2026-06-23 — HPC nine-stage end-to-end (Slurm `--stage all`)

**Task:** Validate full methylation-only pipeline on HPC via Slurm DAG submit after PATH fixes for Bismark/bowtie2 and ALLCools/tabix.

**Files changed:**
- `docs/developers/logs.md`, `AGENTS.md`, `README.md`, `docs/developers/contracts.md`

**Summary:**
- HPC Slurm `--stage all --submit` completed all nine stages (`fastp_split` → `bam_to_allc`) on cluster compute nodes.
- Supersedes prior "Slurm generation only / cluster submit not tested" notes for the nine-stage driver.

**Checks performed:**
- HPC manual validation (user-reported; full Slurm DAG through `09_bam_to_allc`)

**Status:** done

**Notes:** gexcb path not exercised. Local e2e remains `work/dd-met5-example` via `run.sh` (see 2026-06-18 entry).

## 2026-06-23 — bam_to_allc PATH for tabix on Slurm compute nodes

**Task:** Fix HPC `bam_to_allc` failure: ALLCools found but internal `tabix` not on PATH on Slurm compute nodes (same class of issue as bowtie2 for Bismark).

**Files changed:**
- `pixi.toml`, `pixi.lock`
- `scripts/bam_to_allc.py`
- `scripts/check_allcools_env.sh`
- `docs/developers/logs.md`

**Summary:**
- Added explicit `htslib` conda dependency (provides `tabix` alongside `samtools`).
- Prepend pixi/conda `bin` to `PATH` in `bam_to_allc` subprocess env so ALLCools can invoke `tabix` on compute nodes without pixi activation.
- `check-allcools-env` now verifies `tabix` on PATH.

**Checks performed:**
- `pixi install`
- `pixi run check-allcools-env`
- `pixi run python scripts/bam_to_allc.py --help`
- `pixi run bam-to-allc-dry-run`
- HPC manual validation: `09_bam_to_allc` sbatch on Slurm compute nodes (user-reported; part of nine-stage HPC e2e)

**Status:** done

**Notes:** Subprocess `PATH` prepend mirrors `bismark_align` bowtie2 fix; Slurm sbatch scripts do not need a separate `export PATH=...` line.

## 2026-06-23 — bismark_align PATH for Slurm compute nodes

**Task:** Fix HPC `bismark_align` failure: Bismark found but internal `bowtie2 --version` returned 32512 (command not found) on Slurm compute nodes without pixi PATH.

**Files changed:**
- `scripts/bismark_align.py`
- `docs/developers/logs.md`

**Summary:**
- Prepend `bismark_bin` directory to `PATH` in the subprocess env when invoking Bismark so `bowtie2` (and `perl` via shebang) resolve on compute nodes that do not activate pixi.

**Checks performed:**
- `pixi run python scripts/bismark_align.py --help`
- `pixi run bismark-align-dry-run`
- HPC manual validation on `work/test` (user-reported; Slurm `03_bismark_align_0001.sbatch`)

**Status:** done

**Notes:** Slurm sbatch scripts do not need a separate `export PATH=...` line when this fix is present; optional sbatch PATH export on HPC was redundant with the subprocess env change.

## 2026-06-23 — Slurm upfront chunk planning for `--stage all`

**Task:** Fix `ValueError: no demux align inputs found for bismark script generation` when running fresh HPC `--stage all --submit` before prior-stage outputs exist.

**Files changed:**
- `scripts/workflow_input_checks.py`, `scripts/make_cmd.py`
- `docs/developers/logs.md`

**Summary:**
- Added `plan_*` helpers for demux FASTQs, Bismark BAMs, split/merge dirs (mirrors existing `plan_fastp_shards` / demux Slurm pattern).
- Slurm per-chunk script generation now discovers existing outputs when present, otherwise falls back to `--number-of-split-parts` for upfront DAG generation (`bismark_align` through `bam_to_allc`).

**Checks performed:**
- `pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all --runner slurm --dry-run --number-of-split-parts 2 --sample-id test-fresh` (empty work dir; all nine stages)
- `pixi run e2e-slurm-dry-run` (existing `work/dd-met5-example` discovery path)

**Status:** done

## 2026-06-23 — pin seekgene/ALLCools install commit

**Task:** Pin seekgene/ALLCools install commit for reproducible `setup-allcools` (mirror Bismark install script).

**Files changed:**
- `scripts/install_seekgene_allcools.sh`
- `docs/developers/logs.md`

**Summary:**
- Replaced floating `SEEKGENE_ALLCOOLS_REF=master` with pinned `SEEKGENE_ALLCOOLS_COMMIT=b84c180752c7bcc090994efc8534852d497f59d2` (allcools 1.2.0).

**Checks performed:**
- `pixi run setup-allcools`, `pixi run check-allcools-env`

**Status:** done

**Notes:** seekgene/ALLCools has no release tags; commit matches prior validation (allcools 1.2.0).

## 2026-06-18 — nine-stage end-to-end (`--stage all` / `run.sh`)

**Task:** Validate full methylation-only pipeline (`fastp_split` → `bam_to_allc`) via local `run.sh` driver.

**Files changed:**
- `docs/developers/logs.md`, `AGENTS.md`

**Summary:**
- Ran `bash work/dd-met5-example/commands/run.sh` on `data/test_R1.fastq.gz` / `data/test_R2.fastq.gz` with `workflow/dd_met5_test.json` (`force_cell_num=10`, `filter_ch=2`, 2 shards).
- All nine stages completed in order: fastp → demux → bismark → bam_sort → count_mapped_reads → estimated_cells → split_bams → merge_fr_bams → bam_to_allc.

**Checks performed:**
- `bash work/dd-met5-example/commands/run.sh` (local, all 9 stages)
- Outputs under `work/dd-met5-example/`: chunks `0001` and `0002` each with ~3,805 merged per-cell BAMs (`split_bams/merged/`) and matching ALLCools outputs (`allcools/<chunk>_merged_fr_bam_allcools/*_allc.gz` + `.count.csv` + `.tbi`)

**Status:** done

**Notes:** Supersedes the seven-stage e2e entry (stages 08–09 added). Slurm DAG validated on HPC (2026-06-23); see that entry. gexcb path not exercised.

## 2026-06-18 — estimated_cells force_cell_num

**Task:** Add `--force-cell-num` to `estimated_cells` for top-N barcode selection by aligned_reads.

**Files changed:**
- `scripts/estimated_cells.py`, `scripts/make_cmd.py`
- `workflow/dd_met5_test.json`
- `docs/developers/contracts.md`, `docs/developers/logs.md`

**Summary:**
- New `force_filter_barcodes()`: nonzero reads only, sorted by reads desc (barcode tie-break), take top N.
- When `force_cell_num` is set (CLI or workflow JSON), `make_cmd.py` emits `--force-cell-num` instead of `--expected-cell-num`.
- `expected_cell_num` threshold path unchanged when `force_cell_num` is absent.

**Checks performed:**
- `estimated_cells.py --help`
- `pixi run estimated-cells-dry-run` (generated command includes `--force-cell-num 10`)
- Real run on `work/dd-met5-example` with `--force-cell-num 10` → `filtered_barcodes=10`

**Status:** done

## 2026-06-18 — bam_to_allc stage

**Task:** Implement `bam_to_allc` stage (ALLCools bam-to-allc on merged per-cell BAMs) and wire into `make_cmd.py` / `--stage all`.

**Files changed:**
- `scripts/bam_to_allc.py`, `scripts/install_seekgene_allcools.sh`, `scripts/check_allcools_env.sh`
- `scripts/workflow_input_checks.py`, `scripts/make_cmd.py`
- `workflow/dd_met5_test.json`, `pixi.toml`, `pixi.lock`
- `docs/developers/contracts.md`, `AGENTS.md`, `docs/developers/logs.md`

**Summary:**
- Per-chunk ALLCools conversion: `split_bams/merged/<chunk>_merged_fr_bam/*.bam` → `allcools/<chunk>_merged_fr_bam_allcools/<barcode>_allc.gz` with `--convert_bam_strandness --tag UR --save_count_df`.
- Stage 09 (methy-only) / 07 (gexcb); requires `genome_fa` + `chrom_size_path` in workflow JSON; seekgene ALLCools via `pixi run setup-allcools`.
- Added `pip` to pixi dependencies for ALLCools install.

**Checks performed:**
- `pixi run setup-allcools`, `check-allcools-env` (allcools 1.2.0)
- `pixi run bam-to-allc-dry-run`, `bam_to_allc.py --help` / `--dry-run` (chunk 0001: 3,805 barcodes)
- Smoke real run: 1 barcode (`AAAGAAGAGAGGAATAA`) → `*_allc.gz` + `.count.csv` + `.tbi` in ~14s
- `pixi run e2e-dry-run` / `e2e-slurm-dry-run` — `09_bam_to_allc` in driver

**Status:** done

**Notes:** Full-chunk real run validated in nine-stage e2e (`work/dd-met5-example`, both chunks). Slurm cluster submit not tested.

## 2026-06-18 — merge_fr_bams stage

**Task:** Implement `merge_fr_bams` stage (forward/reverse per-cell BAM merge) and wire into `make_cmd.py` / `--stage all`.

**Files changed:**
- `scripts/merge_fr_bams.py`, `scripts/workflow_input_checks.py`, `scripts/make_cmd.py`
- `workflow/dd_met5_test.json`, `pixi.toml`
- `docs/developers/contracts.md`, `AGENTS.md`, `docs/developers/logs.md`

**Summary:**
- Per-chunk merge of `split_bams/<chunk>.{forward,reverse}_1/<barcode>.bam` → `split_bams/merged/<chunk>_merged_fr_bam/<barcode>.bam` via `samtools merge -n` (single-strand copy fallback).
- Aggregates `*_merge_filtered_barcode` and `*_merge_filtered_barcode_reads_counts.csv` across strands.
- Stage 08 (methy-only) / 06 (gexcb); local single script or Slurm per-chunk sbatch.

**Checks performed:**
- `pixi run merge-fr-bams-dry-run`, `merge_fr_bams.py --help` / `--dry-run`
- Real run chunk `0001`: 3,805 merged BAMs in ~56s; `samtools quickcheck` pass; sample barcode read counts forward+reverse=merged (10+10=20)
- `pixi run e2e-dry-run` / `e2e-slurm-dry-run` — `08_merge_fr_bams` in driver

**Status:** done

**Notes:** allcools not in scope. Full-sample merge (both chunks) not run in this check; chunk 0001 only.

## 2026-06-18 — seven-stage end-to-end (`--stage all` / `run.sh`)

**Task:** Validate full methylation-only pipeline (`fastp_split` → `split_bams`) via local `run.sh` driver.

**Files changed:**
- `docs/developers/logs.md`, `AGENTS.md`

**Summary:**
- Ran `bash work/dd-met5-example/commands/run.sh` on `data/test_R1.fastq.gz` / `data/test_R2.fastq.gz` with `workflow/dd_met5_test.json` (`expected_cell_num=10`, `filter_ch=2`, 2 shards).
- All seven stages completed in order: fastp → demux → bismark → bam_sort → count_mapped_reads → estimated_cells → split_bams.

**Checks performed:**
- `bash work/dd-met5-example/commands/run.sh` (local, all 7 stages)
- Outputs under `work/dd-met5-example/`: 2 demux chunks (~224k valid reads/chunk, C→T 0.997); `qc.CtoT.tsv` 23,610 barcodes; `cells/filtered_barcode` 3,804 cells; per-chunk split BAMs (~3,804 forward/reverse cells each in `split_bams/`)

**Status:** done

**Notes:** Slurm `run.sbatch` generation only (cluster submit not tested). gexcb path not exercised in this e2e run.

## 2026-06-18 — step3 count / estimate / split (methylation-only + gexcb)

**Task:** Implement `count_mapped_reads`, `estimated_cells`, and `split_bams` stages with mutually exclusive `expected_cell_num` (default 3000) vs `gexcb` barcode modes.

**Files changed:**
- `scripts/count_mapped_reads.py`, `scripts/estimated_cells.py`, `scripts/split_bams.py`
- `scripts/make_cmd.py`, `scripts/workflow_input_checks.py`
- `workflow/dd_met5_test.json`, `pixi.toml`, `pixi.lock`
- `docs/developers/contracts.md`, `AGENTS.md`

**Summary:**
- Methylation-only path: unsorted BAM CB-tag counts → 99th-percentile cell filter → name-sorted BAM split by barcode (`split_bams/<chunk>.{forward,reverse}_1/<CB>.bam`).
- RNA path: `gexcb` skips count/estimate; `split_bams` uses RNA barcodes directly.
- `make_cmd` dynamic stage sequence (5–7 stages); script prefixes `05`/`06`/`07` or `05` for gexcb-only split.
- Added `pysam`; pixi dry-run helpers for new stages.

**Checks performed:**
- `pixi install`
- `--help` and `--dry-run` on three new scripts and `e2e-dry-run`
- Mutual exclusion: `--gexcb` + workflow `expected_cell_num` raises error
- Real run on `work/dd-met5-example` chunk `0001`: 26,128 barcodes counted; 9,225 filtered; 9,216 forward cells with reads; `samtools view` confirms `CB:Z:` in split BAM

**Status:** done

**Notes:** `merge_fr_bams` / allcools not in scope. gexcb real-data run not tested (no RNA barcodes in dev). Reverse-strand split ~9 min for chunk 0001 with 9k barcodes.

## 2026-06-17 — bam_sort stage

**Task:** Implement `bam_sort` stage (samtools name-sort of Bismark PE BAMs); wire into `make_cmd.py` and `--stage all`.

**Files changed:**
- `scripts/bam_sort.py`, `scripts/make_cmd.py`, `scripts/workflow_input_checks.py`
- `workflow/dd_met5_test.json`, `pixi.toml`
- `docs/developers/contracts.md`, `AGENTS.md`

**Summary:**
- Added per-chunk name-sort: `samtools sort -n` on forward/reverse BAMs in `align/` → `*_sortbyname.bam`.
- `discover_bismark_pe_bams()` pairs `<chunk>.forward_1_bismark_bt2_pe.bam` / `.reverse_1_bismark_bt2_pe.bam`.
- `make_cmd` emits `04_bam_sort.sh` (local) or per-chunk `04_bam_sort_<chunk>.sbatch` (Slurm); `STAGE_SEQUENCE` now includes `bam_sort`.
- Skip re-sort when sortbyname output is newer than input.

**Checks performed:**
- `pixi run python scripts/bam_sort.py --help`
- `pixi run bam-sort-dry-run` (generated command includes `samtools sort -n`)
- `pixi run e2e-dry-run` / `e2e-slurm-dry-run` (`04_bam_sort` in stage list)
- Real run on `work/dd-met5-example/align/` chunk `0001`: forward 25M + reverse 28M sortbyname BAMs; `samtools view` shows adjacent read pairs share QNAME; re-run logs `skipped=1`

**Status:** done

**Notes:** Per-cell `split_bams` (SeekSoulMethyl `step3_split_bams.py`) is the natural follow-on stage.

## 2026-06-17 — demux filter_ch (CH chimeric filtering)

**Task:** Add SeekSoulMethyl-aligned `filter_ch` to `demux_extract_bc` (default 2).

**Files changed:**
- `scripts/demux_extract_bc.py`, `scripts/make_cmd.py`, `workflow/dd_met5_test.json`, `docs/developers/contracts.md`

**Summary:**
- Ported `should_filter_read_ch_pattern()` and strand-specific CH regexes; filter runs post-adapter-trim, pre-FASTQ output.
- CLI `--filter-ch` (default 2); workflow key `filter_ch`; `stats.json` records `filter_ch` and `funnel.barcode_passed.chimeric_filtered`.
- `make_cmd` passes `--filter-ch` in local/Slurm demux commands; default applies for `demux_extract_bc` and `--stage all`.

**Checks performed:**
- `pixi run python scripts/demux_extract_bc.py --help` (`--filter-ch` present)
- `pixi run demux-dry-run` (generated command includes `--filter-ch 2`)
- Chunk `0001` on `work/dd-met5-example/shard_fastq/`: `filter_ch=2` → `valid=224161`, `chimeric_filtered=115154` (fwd 45315, rev 69839); `filter_ch=0` → `valid=339315` (matches pre-change); `CtoT=0.997` unchanged

**Status:** done

**Notes:** `linker.tsv` / C→T QC unaffected (raw R1, pre-trim). Validation outputs under `work/dd-met5-example/demux_validate/`. Bismark on `fc2_0001` → `align_fc2/`: pooled CpG 77.3%, CHG 1.62%, CHH 1.91%, mapped 88.0%, confident 79.5% (vs pre-filter CHG 7.3%/CHH 8.1%/mapped 77.4%; SeekSoulMethyl ref ~77/1.6/1.9/89.9/81.3%).

## 2026-06-17 — demux_extract_bc tqdm progress bar

**Task:** Add read-pair progress bar to `demux_extract_bc` using tqdm; total from `fastp.json`.

**Files changed:**
- `scripts/demux_extract_bc.py`, `scripts/make_cmd.py`, `pixi.toml`, `pixi.lock`

**Summary:**
- Added `--fastp-json` and `--total-reads` CLI options; auto-detect `shard_fastq/fastp.json` when omitted.
- Per-chunk total = `read1_after_filtering.total_reads` ÷ shard count (fastp `--split` does not emit per-chunk stats).
- Wrapped main read-pair loop with tqdm on stderr; `make_cmd` passes `--fastp-json` in local and Slurm demux commands.

**Checks performed:**
- `pixi install` (tqdm added)
- `pixi run python scripts/demux_extract_bc.py --help`
- `pixi run demux-dry-run` (generated script includes `--fastp-json`)
- Single-chunk run on `work/dd-met5-example/shard_fastq/0001`: tqdm reached 496k pairs; `funnel.total=496450` (estimate 496390 from fastp)

**Status:** done

**Notes:** Uneven fastp shard split may leave progress at ~100.01%; non-fastp inputs fall back to indeterminate tqdm.

## 2026-06-17 — demux shard discovery fix (dbit-matrix pattern)

**Task:** Fix `FileNotFoundError` when `--stage all --submit` runs demux before fastp shard paths are known at script-generation time.

**Files changed:**
- `scripts/make_cmd.py`, `scripts/workflow_input_checks.py`

**Summary:**
- Local demux batch script now globs `shard_fastq/*.R1.fq.gz` at runtime (mirrors `dbit-matrix/scripts/make_cmd.py`).
- Slurm demux uses config-based `{chunk}.R1.fq.gz` paths via fixed `plan_fastp_shards` (fastp 0.24.3 `--split` output naming).
- Removed stale `R1.fq.gz` / `R1_001.fq.gz` fallback in `plan_fastp_shards`.

**Checks performed:**
- `pixi run demux-dry-run`
- `pixi run e2e-dry-run`
- `pixi run e2e-slurm-dry-run` (Slurm sbatch shows `0001.R1.fq.gz` paths)

**Status:** done

## 2026-06-16 — bismark_align stage

**Task:** Implement `bismark_align` stage (seekgene Bismark forward/reverse per demux chunk); wire into `make_cmd.py` and `--stage all`.

**Files changed:**
- `scripts/bismark_align.py`, `scripts/make_cmd.py`, `scripts/workflow_input_checks.py`
- `workflow/dd_met5_test.json`, `pixi.toml`
- `docs/developers/contracts.md`

**Summary:**
- Added per-chunk Bismark alignment: forward (`--add_barcode --add_umi`) and reverse (`--pbat`) into `work/<sample>/align/`.
- Extended `make_cmd.py` with `03_bismark_align` local batch and per-chunk Slurm sbatch; `STAGE_SEQUENCE` now includes `bismark_align`.
- Workflow JSON: `bismark_ref`, `bismark_parallel`, `bismark_max_insert`, `bismark_bin`; `pixi run bismark-align-dry-run` helper.

**Checks performed:**
- `pixi run check-bismark-env`
- `pixi run bismark-align-dry-run`
- `pixi run python scripts/bismark_align.py --help`
- `pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage bismark_align` (generated `03_bismark_align.sh`)
- `pixi run e2e-dry-run`
- Real run: chunk `0001` on `work/dd-met5-example/demux/` (~2 min, `--bismark-parallel 4`)
  - outputs: `0001.forward_1_bismark_bt2_pe.bam`, `0001.reverse_1_bismark_bt2_pe.bam` + PE reports
  - `samtools view` confirms `CB:Z:` and `UR:Z:` tags in BAM

**Status:** done

**Notes:** Sort/dedup/split BAM not in scope. Re-run `pixi run setup-bismark` after pixi env rebuild. Slurm `03_bismark_align_<chunk>.sbatch` generation follows demux pattern; cluster submit not tested.

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
