# Development status

Living snapshot of what is reliable today. History and per-task checks: [logs.md](logs.md).

## Validated

Per-stage confidence (what was actually exercised):

- fastp_split: local real run on `work/dd-met5-example`; dry-run + `e2e-dry-run`
- demux_extract_bc: local real run + HPC via methylation-only Slurm DAG (`dd_met5_test.json`)
- regroup_shards: local e2e + HPC e2e (prefix chunks `A`/`G`/`T`)
- bismark_align: local real run; HPC compute-node submit (bowtie2 PATH fix)
- bam_sort: local twelve-stage e2e (prefix chunks)
- count_mapped_reads / estimated_cells: local twelve-stage e2e
- split_bams / merge_fr_bams: local twelve-stage e2e
- bam_to_allc: local e2e + HPC e2e; smoke on single barcode
- saturation: local eleven-stage e2e + HPC eleven-stage e2e; algorithm review fixes applied
- qc_summary: local real run (`work/dd-met5-example`, 50 cells); dry-run in twelve-stage driver — **needs biological sanity check**; not HPC-submitted
- allc_to_matrix: local run on `work/dd-met5-example` (50 cells, `filtered_barcode`) + `make_cmd` dry-run; MethSCAn `prepare` parity **passed** (all-context: exact match on `dd-met5-example`). Default `meth_context=CG` is intentional (CpG-only). Also exercised on `work/C283_Brain_DNAme_S1` (300 cells); Slurm meth path not cluster-tested
- meth_smooth: local run on `work/dd-met5-example` (50 cells, CG matrix) + `meth-smooth-dry-run`; MethSCAn `smooth` parity **passed** (max abs diff `0.0` on main chroms vs v1.1.0 reference)
- meth_scan: local run on `work/dd-met5-example` (2 VMRs, default params) + `meth-scan-dry-run`; MethSCAn `scan` parity **passed** (exact BED match vs v1.1.0 reference)
- meth_matrix: local sparse run on `work/dd-met5-example` (2 VMRs, `vmrs.bed` fallback) + `meth-matrix-dry-run`; MethSCAn `matrix --sparse` parity **passed** (19 entries; exact `(row, col, mfrac)` vs v1.1.0). Dense mode smoke-tested (`--dense`)

Workflow drivers:

- methylation-only `run.sh` (`workflow/dd_met5_test.json`): twelve stages through `qc_summary` (stage script generation + local runs)
- methylation-only with `run_meth_analysis: true`: fifteen-stage script generation validated (`meth-e2e-dry-run`: `allc_to_matrix` → `meth_smooth` → `meth_scan`)
- methylation-only with `run_meth_analysis` + `run_meth_matrix`: sixteen-stage script generation validated (`meth-e2e-dry-run`); full sixteen-stage local `run.sh` execution not run end-to-end in one command
- methylation-only `run.sbatch` (`dd_met5_test.json`): HPC submit validated through **saturation** (eleven stages at time of HPC run); `qc_summary` Slurm path not cluster-tested; meth analysis Slurm paths dry-run only
- `workflow/dd_met5_slurm.json`: Slurm command generation dry-run only; production-scale cluster submit not validated

## Partially validated / not exercised

- **Meth analysis path** (`run_meth_analysis: true`): Phase 1–2 complete through `meth_matrix`; `meth_diff` / `meth_profile` not implemented (Phase 3)
- **gexcb path** (`workflow/dd_met5_gexcb_test.json`): local `split_bams` → `merge_fr_bams` → `bam_to_allc` on 22 barcodes (reuses prior align BAMs); full `--stage all` from raw FASTQ not run; Slurm gexcb not tested
- **Automated tests:** none yet (manual validation only)

## Known limitations (out of scope)

- Spike-in output ([`stage_notes/demux_extract_bc.md`](stage_notes/demux_extract_bc.md))
- Multi-barcode rescue (ambiguous HD=1 matches discarded)
- `generate-dataset`, merged ALLC matrix, per-cell JSON/HTML reports ([`stage_notes/bam_to_allc.md`](stage_notes/bam_to_allc.md), [`stage_notes/qc_summary.md`](stage_notes/qc_summary.md))
- Sample-wide barcode union across analysis chunks in per-chunk stages ([`chunk_model.md`](chunk_model.md))
- `meth_matrix_filter` (MethSCAn `filter`) — skipped; cell QC in `allc_to_matrix` via `filtered_barcode`

## Do not change silently

Contract-sensitive surfaces — update the canonical doc layer per [`doc-system.md`](doc-system.md) and bump this file if validation posture changes:

| Surface | Canonical reference |
|---------|---------------------|
| Barcode correction (HD=1, no multi-rescue) | [`stage_notes/demux_extract_bc.md`](stage_notes/demux_extract_bc.md) |
| `chunk_id` / analysis-chunk semantics | [`chunk_model.md`](chunk_model.md) |
| BAM read-name format (`CB`/`UR` from demux names) | [`stage_notes/bismark_align.md`](stage_notes/bismark_align.md), [`contracts.md`](contracts.md#split_bams) |
| `cells_summary.tsv` columns, `*_mc_rate` pooling, `mito_CG_mc_rate` | [`qc_metrics.md`](qc_metrics.md#cells_summarytsv) |
