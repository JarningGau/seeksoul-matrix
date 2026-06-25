# bismark_align

Implementation notes for [`bismark_align`](../contracts.md#bismark_align). Not normative for I/O.

## Behavior

- Forward: `bismark --genome <bismark_ref> --parallel <N> -1 <fwd_r1> -2 <fwd_r2> -o <align_dir> -X <max_insert> --add_barcode --add_umi`.
- Reverse: same with `--pbat`.
- BAM read names retain demux format; Bismark writes `CB` and `UR` tags from read names.

## Defaults and workflow keys

| Key | Default | Role |
|-----|---------|------|
| `bismark_parallel` | `8` | Bismark `--parallel` |
| `bismark_max_insert` | `1000` | Bismark `-X` max insert size |
| `bismark_ref` | (required) | parent directory of `Bisulfite_Genome/` |

## Skip / incremental re-run

Per analysis chunk and strand; re-run replaces BAMs when invoked.

## Toolchain and environment

- Uses seekgene [Bismark](https://github.com/seekgene/Bismark) fork (stock bioconda Bismark lacks `--add_barcode` / `--add_umi`); requires `bowtie2` on PATH (`pixi run check-bismark-env`).
- `bismark_align.py` prepends the Bismark executable directory to subprocess `PATH` so Bismark can invoke `bowtie2` on Slurm compute nodes without pixi activation.

## SeekSoulMethyl / dbit-matrix alignment

Defaults aligned with SeekSoulMethyl `step2.nf`.

## Out of scope

`samtools sort`, UMI dedup, and per-cell BAM splitting.
