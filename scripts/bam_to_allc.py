#!/usr/bin/env python3
"""Convert merged per-cell BAMs to ALLC format via ALLCools bam-to-allc."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import workflow_input_checks as wic

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")


@dataclass(frozen=True)
class AllcChunk:
    chunk_id: str
    input_dir: Path
    filtered_barcode: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run ALLCools bam-to-allc on merged per-cell BAMs "
            "(DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help="Sample work directory; scans split_bams/merged/*_merged_fr_bam.",
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument(
        "--input-dir",
        help="Merged per-cell BAM directory for single-chunk mode.",
    )
    parser.add_argument(
        "--filtered-barcode",
        help="Filtered barcode list for single-chunk mode.",
    )
    parser.add_argument(
        "--output-dir",
        help="ALLCools output root for single-chunk mode (default: work/allcools).",
    )
    parser.add_argument(
        "--genome-fa",
        required=True,
        help="Reference genome FASTA for allcools --reference_fasta.",
    )
    parser.add_argument(
        "--chrom-size-path",
        required=True,
        help="Chromosome sizes BED (workflow parity; not passed to bam-to-allc).",
    )
    parser.add_argument(
        "--align-method",
        default="bismark",
        choices=["bismark"],
        help="Alignment method; bismark adds --convert_bam_strandness. Default: bismark.",
    )
    parser.add_argument(
        "--allcools-tag",
        default="UR",
        help="BAM tag for UMI correction. Default: UR.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=8,
        help="CPU cores for parallel per-barcode conversion. Default: 8.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or name. Default: samtools.",
    )
    parser.add_argument(
        "--allcools-bin",
        default="allcools",
        help="allcools executable path or name. Default: allcools.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths and actions without running conversion.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def resolve_env_executable(name: str) -> str:
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.is_file():
        return str(candidate)
    return name


def load_filtered_barcodes(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    barcodes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            barcode = line.strip()
            if barcode:
                barcodes.add(barcode)
    return barcodes


def allc_output_path(outdir: Path, barcode: str) -> Path:
    return outdir / f"{barcode}_allc.gz"


def allc_output_exists(outdir: Path, barcode: str) -> bool:
    path = allc_output_path(outdir, barcode)
    return path.is_file() and path.stat().st_size > 0


def collect_bams(input_dir: Path, filtered: set[str] | None) -> list[Path]:
    bams = sorted(input_dir.glob("*.bam"))
    if filtered is None:
        return bams
    return [bam for bam in bams if bam.stem in filtered]


def run_allcools_for_barcode(
    bam_path: Path,
    *,
    genome_fa: str,
    barcode: str,
    outdir: Path,
    align_method: str,
    tag: str | None,
    samtools_bin: str,
    allcools_bin: str,
) -> str:
    if allc_output_exists(outdir, barcode):
        return f"{barcode} skipped"

    sorted_bam = outdir / f"{barcode}_sort.bam"
    sorted_index = outdir / f"{barcode}_sort.bam.bai"
    allc_prefix = outdir / f"{barcode}_allc"

    subprocess.run(
        [samtools_bin, "sort", "-o", str(sorted_bam), str(bam_path)],
        check=True,
    )
    subprocess.run(
        [samtools_bin, "index", str(sorted_bam)],
        check=True,
    )

    cmd = [
        allcools_bin,
        "bam-to-allc",
        "--bam_path",
        str(sorted_bam),
        "--reference_fasta",
        genome_fa,
        "--output_path",
        str(allc_prefix),
        "--save_count_df",
    ]
    if align_method == "bismark":
        cmd.append("--convert_bam_strandness")
    if tag:
        cmd.extend(["--tag", tag])

    subprocess.run(cmd, check=True)

    if sorted_bam.exists():
        sorted_bam.unlink()
    if sorted_index.exists():
        sorted_index.unlink()

    return f"{barcode} done"


def allcools_outdir(output_root: Path, input_dir: Path) -> Path:
    return output_root / f"{input_dir.name}_allcools"


def process_chunk(
    chunk: AllcChunk,
    *,
    output_root: Path,
    genome_fa: str,
    chrom_size_path: str,
    align_method: str,
    tag: str | None,
    cores: int,
    samtools_bin: str,
    allcools_bin: str,
    dry_run: bool,
) -> None:
    del chrom_size_path  # required for workflow parity; not used by bam-to-allc CLI

    outdir = allcools_outdir(output_root, chunk.input_dir)
    filtered = load_filtered_barcodes(chunk.filtered_barcode)
    bams = collect_bams(chunk.input_dir, filtered)

    print(f"[bam_to_allc] chunk={chunk.chunk_id}")
    print(f"[bam_to_allc] input_dir={chunk.input_dir}")
    print(f"[bam_to_allc] filtered_barcode={chunk.filtered_barcode}")
    print(f"[bam_to_allc] output_dir={outdir}")
    print(f"[bam_to_allc] genome_fa={genome_fa}")
    print(f"[bam_to_allc] barcode_count={len(bams)} cores={cores}")

    if dry_run:
        if bams:
            sample = bams[0].stem
            print(
                "[bam_to_allc] sample_command="
                + quoted(
                    [
                        allcools_bin,
                        "bam-to-allc",
                        "--bam_path",
                        f"{outdir}/{sample}_sort.bam",
                        "--reference_fasta",
                        genome_fa,
                        "--output_path",
                        f"{outdir}/{sample}_allc",
                        "--save_count_df",
                        "--convert_bam_strandness",
                        "--tag",
                        tag or "",
                    ]
                )
            )
        return

    outdir.mkdir(parents=True, exist_ok=True)
    workers = max(1, min(cores, len(bams) or 1))

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                run_allcools_for_barcode,
                bam,
                genome_fa=genome_fa,
                barcode=bam.stem,
                outdir=outdir,
                align_method=align_method,
                tag=tag,
                samtools_bin=samtools_bin,
                allcools_bin=allcools_bin,
            )
            for bam in bams
        ]
        for future in as_completed(futures):
            future.result()

    kept = sum(1 for path in outdir.glob("*_allc.gz") if path.stat().st_size > 0)
    print(f"[bam_to_allc] allc_files={kept}")


def resolve_jobs(args: argparse.Namespace) -> tuple[Path, list[AllcChunk]]:
    has_work_mode = args.work_path is not None
    has_single_mode = args.input_dir is not None
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or --input-dir mode")
    if not has_work_mode and not has_single_mode:
        raise ValueError("missing input mode: provide --work-path or --input-dir")

    if has_work_mode:
        work_path = Path(args.work_path)
        merged_root = work_path / "split_bams" / "merged"
        output_root = work_path / "allcools"
        chunks = [
            AllcChunk(chunk_id, bam_dir, filtered_barcode)
            for chunk_id, bam_dir, filtered_barcode in wic.discover_merged_fr_bam_chunks(
                merged_root
            )
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(f"no merged FR BAM chunks found under {merged_root}")
        return output_root, chunks

    if args.filtered_barcode is None:
        raise ValueError("single-chunk mode requires --filtered-barcode")
    input_dir = Path(args.input_dir)
    filtered_barcode = Path(args.filtered_barcode)
    if args.output_dir:
        output_root = Path(args.output_dir)
    else:
        output_root = input_dir.parent.parent / "allcools"
    chunk_id = args.chunk_id or input_dir.name.replace("_merged_fr_bam", "")
    return output_root, [
        AllcChunk(chunk_id, input_dir, filtered_barcode),
    ]


def main() -> int:
    args = parse_args()
    if args.cores <= 0:
        raise ValueError("cores must be > 0")

    samtools_bin = args.samtools_bin
    if samtools_bin == "samtools":
        samtools_bin = resolve_env_executable("samtools")
    allcools_bin = args.allcools_bin
    if allcools_bin == "allcools":
        allcools_bin = resolve_env_executable("allcools")

    wic.require_file("genome_fa", wic.resolve_config_path(args.genome_fa))
    wic.require_file("chrom_size_path", wic.resolve_config_path(args.chrom_size_path))
    wic.require_optional_executable_path("samtools_bin", samtools_bin)
    wic.require_optional_executable_path("allcools_bin", allcools_bin)

    output_root, chunks = resolve_jobs(args)
    print(f"[bam_to_allc] output_root={output_root}")
    print(f"[bam_to_allc] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[bam_to_allc] progress={index}/{len(chunks)}")
        process_chunk(
            chunk,
            output_root=output_root,
            genome_fa=str(wic.resolve_config_path(args.genome_fa)),
            chrom_size_path=str(wic.resolve_config_path(args.chrom_size_path)),
            align_method=args.align_method,
            tag=args.allcools_tag or None,
            cores=args.cores,
            samtools_bin=samtools_bin,
            allcools_bin=allcools_bin,
            dry_run=args.dry_run,
        )

    print("[bam_to_allc] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
