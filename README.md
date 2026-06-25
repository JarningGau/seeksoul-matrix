# SeepSpace Spatial Methylation (Slide-Tag like)

## DD-MET5 Chemistry

### Library structure
`CB`-`UB`-`TSO`-`17L`-`ME`-`9bp`-`Insert`-`9bp`

- `CB`: 17 bp cell barcode (do not contain C)
- `UB`: 12 bp UMI sequence (do not contain C)
- `TSO`: 13 bp TSO sequence TTTCTTATATGGG
- `17L`: 17 bp fixed sequence CgtCCgtCgttgCtCgt
- `ME`: 19 bp fixed sequence AGATGTGTATAAGAGACAG
- `9bp`: extension sequence from the Tn5 insertion fragment

Since the enzymatic treatment converts unmethylated cytosines (C) to thymines (T), the barcodes used for methylation data do not contain any C bases. In contrast, the C bases in TSO, 17L, and ME are not methylated and will be converted to T during the enzymatic process; we use these fixed sequences to calculate the C-to-T conversion rate.

### Cell Barcode Whitelist Validation

**Whitelist**: `whitelist/DD-MET5/U3CB_methylation.txt.gz` (829,440 × 17 bp cell barcodes)

----

## Installation

Requires [pixi](https://pixi.sh/) on **linux-64**.

```bash
git clone https://github.com/JarningGau/seeksoul-matrix.git
cd seeksoul-matrix
pixi install
```

The pixi environment provides core tools (`fastp`, `cutadapt`, `bowtie2`, `samtools`, `htslib`/`tabix`, Python 3.11). Two seekgene forks are installed separately into the active environment (not conda packages):

```bash
# bismark_align — seekgene/Bismark with --add_barcode / --add_umi
pixi run setup-bismark
pixi run check-bismark-env

# bam_to_allc — seekgene/ALLCools with UR-tag UMI dedup (requires tabix from htslib)
pixi run setup-allcools
pixi run check-allcools-env
```

Re-run `setup-bismark` and `setup-allcools` after `pixi install` recreates `.pixi/envs/default`.

----

## Running the pipeline

Twelve stages from `fastp_split` through `qc_summary`; stable I/O contracts and validation posture live under [`docs/developers/`](docs/developers/).

| Resource | Purpose |
|----------|---------|
| [`docs/developers/contracts.md`](docs/developers/contracts.md) | Stage input/output paths |
| [`docs/developers/status.md`](docs/developers/status.md) | What is validated today |
| [`docs/developers/chunk_model.md`](docs/developers/chunk_model.md) | Read-order vs analysis chunks, barcode modes |
| [`workflow/dd_met5_test.json`](workflow/dd_met5_test.json) | Local / CI-style test config (methylation-only) |
| [`workflow/dd_met5_slurm.json`](workflow/dd_met5_slurm.json) | Production Slurm config |
| [`workflow/dd_met5_gexcb_test.json`](workflow/dd_met5_gexcb_test.json) | RNA-barcode (`gexcb`) path |
| [`examples/`](examples/) | Example local and HPC invocations |

Generate runnable scripts with the workflow driver:

```bash
pixi run python scripts/make_cmd.py --workflow-config workflow/dd_met5_test.json --stage all
bash work/<sample>/commands/run.sh
```

Per-stage and end-to-end dry-runs (no cluster submit): `pixi run fastp-dry-run`, `pixi run e2e-dry-run`, `pixi run e2e-slurm-dry-run`. Full task list: [`AGENTS.md`](AGENTS.md).
