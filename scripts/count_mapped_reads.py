#!/usr/bin/env python3
"""Count aligned reads per cell barcode from unsorted Bismark BAMs (CB tag)."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam

import workflow_input_checks as wic


@dataclass(frozen=True)
class CountChunk:
    chunk_id: str
    forward_bam: Path
    reverse_bam: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count aligned reads per barcode from unsorted Bismark PE BAMs "
            "using the CB tag (DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help=(
            "Sample work directory. Scans align/*.forward_1_bismark_bt2_pe.bam "
            "and writes counts CSVs alongside inputs."
        ),
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument("--bam", help="Single BAM path for one-strand mode.")
    parser.add_argument(
        "--output-dir",
        help="Output directory for single-BAM mode (default: align/ under work-path).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths and actions without reading BAMs or writing CSVs.",
    )
    return parser.parse_args()


def should_skip_count(input_bam: Path, output_csv: Path) -> bool:
    if not output_csv.is_file():
        return False
    if not input_bam.is_file():
        return False
    return output_csv.stat().st_mtime >= input_bam.stat().st_mtime


def count_bam_barcode_reads(bam_path: Path) -> dict[str, int]:
    barcode_counts: dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.has_tag("CB"):
                barcode_counts[read.get_tag("CB")] += 1
    return dict(barcode_counts)


def write_counts_csv(output_csv: Path, barcode_counts: dict[str, int]) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(
        barcode_counts.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "aligned_reads"])
        for barcode, aligned_reads in rows:
            writer.writerow([barcode, aligned_reads])


def count_one_bam(
    input_bam: Path,
    *,
    dry_run: bool,
    chunk_id: str,
    strand: str,
) -> None:
    output_csv = wic.counts_output_path(input_bam)
    print(f"[count_mapped_reads] chunk={chunk_id} strand={strand}")
    print(f"[count_mapped_reads] input_bam={input_bam}")
    print(f"[count_mapped_reads] output_csv={output_csv}")
    if dry_run:
        return
    if should_skip_count(input_bam, output_csv):
        print("[count_mapped_reads] skipped=1")
        return
    barcode_counts = count_bam_barcode_reads(input_bam)
    write_counts_csv(output_csv, barcode_counts)
    print(f"[count_mapped_reads] barcodes={len(barcode_counts)}")


def resolve_jobs(args: argparse.Namespace) -> list[CountChunk]:
    has_work_mode = args.work_path is not None
    has_single_mode = args.bam is not None
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or --bam")
    if not has_work_mode and not has_single_mode:
        raise ValueError("missing input mode: provide --work-path or --bam")

    if has_work_mode:
        work_path = Path(args.work_path)
        chunks = [
            CountChunk(chunk_id, fwd_bam, rev_bam)
            for chunk_id, fwd_bam, rev_bam in wic.discover_bismark_pe_bams(work_path / "align")
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(f"no unsorted Bismark PE BAMs found under {work_path / 'align'}")
        return chunks

    bam_path = Path(args.bam)
    chunk_id = args.chunk_id or bam_path.name.split(".", 1)[0]
    return [CountChunk(chunk_id, bam_path, bam_path)]


def count_chunk(chunk: CountChunk, *, dry_run: bool) -> None:
    strands = (
        ("forward", chunk.forward_bam),
        ("reverse", chunk.reverse_bam),
    )
    seen: set[Path] = set()
    for strand, input_bam in strands:
        if input_bam in seen:
            continue
        seen.add(input_bam)
        count_one_bam(
            input_bam,
            dry_run=dry_run,
            chunk_id=chunk.chunk_id,
            strand=strand,
        )


def main() -> int:
    args = parse_args()
    chunks = resolve_jobs(args)
    print(f"[count_mapped_reads] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[count_mapped_reads] progress={index}/{len(chunks)}")
        count_chunk(chunk, dry_run=args.dry_run)

    print("[count_mapped_reads] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
