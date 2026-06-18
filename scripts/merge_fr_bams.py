#!/usr/bin/env python3
"""Merge forward and reverse per-cell split BAMs into one BAM per barcode."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import workflow_input_checks as wic


@dataclass(frozen=True)
class MergeChunk:
    chunk_id: str
    forward_dir: Path
    reverse_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge forward and reverse per-cell split BAMs "
            "(DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help="Sample work directory; scans split_bams/*.forward_1 dirs.",
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument(
        "--forward-dir",
        help="Forward split BAM directory for single-chunk mode.",
    )
    parser.add_argument(
        "--reverse-dir",
        help="Reverse split BAM directory for single-chunk mode.",
    )
    parser.add_argument(
        "--output-dir",
        help="Merge output root for single-chunk mode (default: split_bams/merged).",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=8,
        help="CPU cores for parallel per-barcode merges. Default: 8.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or name. Default: samtools.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths and actions without merging BAMs.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def resolve_samtools_bin(raw: str) -> str:
    if raw != "samtools":
        return raw
    candidate = Path(sys.executable).resolve().parent / "samtools"
    if candidate.is_file():
        return str(candidate)
    return raw


def collect_bam_map(split_dir: Path | None) -> dict[str, Path]:
    if split_dir is None or not split_dir.is_dir():
        return {}
    return {path.stem: path for path in sorted(split_dir.glob("*.bam"))}


def bam_is_valid(bam_path: Path, samtools_bin: str) -> bool:
    if not bam_path.is_file():
        return False
    try:
        subprocess.run(
            [samtools_bin, "quickcheck", str(bam_path)],
            check=True,
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def safe_merge_bam(
    f_bam: str | None,
    r_bam: str | None,
    out_bam: str,
    samtools_bin: str,
) -> None:
    output_path = Path(out_bam)
    if bam_is_valid(output_path, samtools_bin):
        return

    if output_path.exists():
        output_path.unlink()

    temp_bam = output_path.with_suffix(".tmp.bam")
    if temp_bam.exists():
        temp_bam.unlink()

    try:
        if f_bam and r_bam:
            cmd = [
                samtools_bin,
                "merge",
                "-n",
                "-f",
                "-@",
                "1",
                str(temp_bam),
                f_bam,
                r_bam,
            ]
            subprocess.run(cmd, check=True)
        elif f_bam:
            shutil.copy2(f_bam, temp_bam)
        elif r_bam:
            shutil.copy2(r_bam, temp_bam)
        else:
            raise ValueError(f"no input BAMs for {out_bam}")

        temp_bam.rename(output_path)
    except Exception:
        if temp_bam.exists():
            temp_bam.unlink()
        raise


def collect_filtered_barcode_files(split_dir: Path | None) -> list[Path]:
    if split_dir is None or not split_dir.is_dir():
        return []
    return sorted(split_dir.glob("*_filtered_barcode"))


def collect_filtered_counts_files(split_dir: Path | None) -> list[Path]:
    if split_dir is None or not split_dir.is_dir():
        return []
    return sorted(split_dir.glob("*_filtered_barcode_reads_counts.csv"))


def merge_filtered_barcodes(
    forward_dir: Path,
    reverse_dir: Path,
    output_path: Path,
) -> None:
    barcodes: set[str] = set()
    for barcode_file in collect_filtered_barcode_files(
        forward_dir
    ) + collect_filtered_barcode_files(reverse_dir):
        with barcode_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                barcode = line.strip()
                if barcode:
                    barcodes.add(barcode)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for barcode in sorted(barcodes):
            handle.write(f"{barcode}\n")


def merge_filtered_counts(
    forward_dir: Path,
    reverse_dir: Path,
    output_path: Path,
) -> None:
    barcode_counts: dict[str, int] = {}
    for counts_file in collect_filtered_counts_files(
        forward_dir
    ) + collect_filtered_counts_files(reverse_dir):
        with counts_file.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            for row in reader:
                if len(row) < 2:
                    continue
                reads = int(row[0])
                barcode = row[1]
                barcode_counts[barcode] = barcode_counts.get(barcode, 0) + reads

    rows = sorted(barcode_counts.items(), key=lambda item: item[1], reverse=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["reads_counts", "barcode"])
        for barcode, reads in rows:
            writer.writerow([reads, barcode])


def merge_chunk_pair(
    chunk: MergeChunk,
    *,
    merged_root: Path,
    cores: int,
    samtools_bin: str,
    dry_run: bool,
) -> None:
    chunk_id = chunk.chunk_id
    merged_bam_dir = merged_root / f"{chunk_id}_merged_fr_bam"
    merge_barcode_path = merged_root / f"{chunk_id}_merge_filtered_barcode"
    merge_counts_path = (
        merged_root / f"{chunk_id}_merge_filtered_barcode_reads_counts.csv"
    )

    fmap = collect_bam_map(chunk.forward_dir)
    rmap = collect_bam_map(chunk.reverse_dir)
    all_barcodes = sorted(set(fmap) | set(rmap))

    print(f"[merge_fr_bams] chunk={chunk.chunk_id}")
    print(f"[merge_fr_bams] forward_dir={chunk.forward_dir}")
    print(f"[merge_fr_bams] reverse_dir={chunk.reverse_dir}")
    print(f"[merge_fr_bams] merged_bam_dir={merged_bam_dir}")
    print(f"[merge_fr_bams] barcode_count={len(all_barcodes)} cores={cores}")

    if dry_run:
        return

    merged_root.mkdir(parents=True, exist_ok=True)
    merged_bam_dir.mkdir(parents=True, exist_ok=True)

    workers = max(1, min(cores, len(all_barcodes) or 1))
    merge_tasks = [
        (
            str(fmap[barcode]) if barcode in fmap else None,
            str(rmap[barcode]) if barcode in rmap else None,
            str(merged_bam_dir / f"{barcode}.bam"),
            samtools_bin,
        )
        for barcode in all_barcodes
    ]

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(safe_merge_bam, f_bam, r_bam, out_bam, samtools)
            for f_bam, r_bam, out_bam, samtools in merge_tasks
        ]
        for future in as_completed(futures):
            future.result()

    merge_filtered_barcodes(chunk.forward_dir, chunk.reverse_dir, merge_barcode_path)
    merge_filtered_counts(chunk.forward_dir, chunk.reverse_dir, merge_counts_path)

    kept = sum(1 for path in merged_bam_dir.glob("*.bam") if path.is_file())
    print(f"[merge_fr_bams] merged_bams={kept}")


def resolve_jobs(args: argparse.Namespace) -> tuple[Path, list[MergeChunk]]:
    has_work_mode = args.work_path is not None
    has_single_mode = any(
        value is not None for value in (args.forward_dir, args.reverse_dir)
    )
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or --forward-dir/--reverse-dir")
    if not has_work_mode and not has_single_mode:
        raise ValueError(
            "missing input mode: provide --work-path or "
            "--forward-dir/--reverse-dir"
        )

    if has_work_mode:
        work_path = Path(args.work_path)
        split_root = work_path / "split_bams"
        merged_root = split_root / "merged"
        chunks = [
            MergeChunk(chunk_id, forward_dir, reverse_dir)
            for chunk_id, forward_dir, reverse_dir in wic.discover_split_bam_chunk_pairs(
                split_root
            )
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(f"no split BAM chunk pairs found under {split_root}")
        return merged_root, chunks

    if not all(value is not None for value in (args.forward_dir, args.reverse_dir)):
        raise ValueError(
            "single-chunk mode requires --forward-dir and --reverse-dir"
        )
    forward_dir = Path(args.forward_dir)
    reverse_dir = Path(args.reverse_dir)
    if args.output_dir:
        merged_root = Path(args.output_dir)
    else:
        merged_root = forward_dir.parent / "merged"
    chunk_id = args.chunk_id or forward_dir.name.split(".", 1)[0]
    return merged_root, [MergeChunk(chunk_id, forward_dir, reverse_dir)]


def main() -> int:
    args = parse_args()
    if args.cores <= 0:
        raise ValueError("cores must be > 0")

    samtools_bin = resolve_samtools_bin(args.samtools_bin)
    wic.require_optional_executable_path("samtools_bin", samtools_bin)

    merged_root, chunks = resolve_jobs(args)
    print(f"[merge_fr_bams] merged_root={merged_root}")
    print(f"[merge_fr_bams] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[merge_fr_bams] progress={index}/{len(chunks)}")
        merge_chunk_pair(
            chunk,
            merged_root=merged_root,
            cores=args.cores,
            samtools_bin=samtools_bin,
            dry_run=args.dry_run,
        )

    print("[merge_fr_bams] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
