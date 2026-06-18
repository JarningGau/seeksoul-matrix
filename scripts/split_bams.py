#!/usr/bin/env python3
"""Split name-sorted Bismark BAMs into per-cell BAM files."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path

import pysam

import workflow_input_checks as wic


@dataclass(frozen=True)
class SplitChunk:
    chunk_id: str
    forward_bam: Path
    reverse_bam: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split name-sorted Bismark PE BAMs by cell barcode "
            "(DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help="Sample work directory; scans align/*sortbyname.bam.",
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument("--bam", help="Single sortbyname BAM for one-strand mode.")
    parser.add_argument(
        "--output-dir",
        help="Split output root for single-BAM mode.",
    )
    barcode_group = parser.add_mutually_exclusive_group()
    barcode_group.add_argument(
        "--filtered-barcode",
        help="Path to filtered barcode list (one barcode per line).",
    )
    barcode_group.add_argument(
        "--gexcb",
        help="Path to RNA filtered barcodes (one barcode per line or .tsv.gz).",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=8,
        help="CPU cores for parallel barcode batch splitting. Default: 8.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths and actions without splitting BAMs.",
    )
    return parser.parse_args()


def load_barcodes(path: Path) -> list[str]:
    barcodes: list[str] = []
    if path.suffix == ".gz":
        handle_ctx = gzip.open(path, "rt", encoding="utf-8")
    else:
        handle_ctx = path.open("r", encoding="utf-8")
    with handle_ctx as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            barcodes.append(line.split("\t")[0].split(",")[0])
    if not barcodes:
        raise ValueError(f"no barcodes found in {path}")
    return barcodes


def split_bams_dirname(sortbyname_bam: Path) -> str:
    stem = sortbyname_bam.stem
    return re.sub(r"_bismark_.*", "", stem)


def split_one_batch(
    bam_path: str,
    outdir: str,
    keep_barcodes: list[str],
    batch_id: int,
) -> dict[str, int]:
    barcode_set = set(keep_barcodes)
    processed_barcodes: set[str] = set()
    barcode_read_counts = {barcode: 0 for barcode in keep_barcodes}

    with pysam.AlignmentFile(bam_path, "rb") as input_bam:
        for barcode, reads_group in groupby(
            input_bam, key=lambda read: read.query_name.split("_", 1)[0]
        ):
            if len(processed_barcodes) >= len(barcode_set):
                break
            if barcode not in barcode_set or barcode in processed_barcodes:
                continue
            output_path = Path(outdir) / f"{barcode}.bam"
            read_count = 0
            with pysam.AlignmentFile(
                str(output_path), "wb", template=input_bam
            ) as outfh:
                for read in reads_group:
                    outfh.write(read)
                    read_count += 1
            barcode_read_counts[barcode] = read_count
            processed_barcodes.add(barcode)

    return barcode_read_counts


def write_filtered_outputs(
    output_dir: Path,
    prefix: str,
    barcode_read_counts: dict[str, int],
) -> None:
    rows = [
        (reads, barcode)
        for barcode, reads in barcode_read_counts.items()
        if reads > 0
    ]
    rows.sort(key=lambda item: item[0], reverse=True)

    counts_path = output_dir / f"{prefix}_filtered_barcode_reads_counts.csv"
    with counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["reads_counts", "barcode"])
        for reads, barcode in rows:
            writer.writerow([reads, barcode])

    barcode_path = output_dir / f"{prefix}_filtered_barcode"
    with barcode_path.open("w", encoding="utf-8", newline="") as handle:
        for _reads, barcode in rows:
            handle.write(f"{barcode}\n")


def split_one_bam(
    sortbyname_bam: Path,
    *,
    split_root: Path,
    keep_barcodes: list[str],
    cores: int,
    dry_run: bool,
    chunk_id: str,
    strand: str,
) -> None:
    prefix = split_bams_dirname(sortbyname_bam)
    output_dir = split_root / prefix
    print(f"[split_bams] chunk={chunk_id} strand={strand}")
    print(f"[split_bams] input_bam={sortbyname_bam}")
    print(f"[split_bams] output_dir={output_dir}")
    print(f"[split_bams] keep_barcodes={len(keep_barcodes)} cores={cores}")
    if dry_run:
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    workers = max(1, min(cores, len(keep_barcodes)))
    sorted_barcodes = sorted(keep_barcodes)
    batch_size = max(1, (len(sorted_barcodes) + workers - 1) // workers)
    batches = [
        sorted_barcodes[index : index + batch_size]
        for index in range(0, len(sorted_barcodes), batch_size)
    ]

    all_barcode_counts: dict[str, int] = {}
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                split_one_batch,
                str(sortbyname_bam),
                str(output_dir),
                batch,
                batch_id,
            )
            for batch_id, batch in enumerate(batches)
            if batch
        ]
        for future in as_completed(futures):
            all_barcode_counts.update(future.result())

    write_filtered_outputs(output_dir, prefix, all_barcode_counts)
    kept = sum(1 for reads in all_barcode_counts.values() if reads > 0)
    print(f"[split_bams] cells_with_reads={kept}")


def resolve_barcode_list(args: argparse.Namespace, work_path: Path | None) -> Path:
    if args.filtered_barcode:
        return Path(args.filtered_barcode)
    if args.gexcb:
        return Path(args.gexcb)
    if work_path is not None:
        return work_path / "cells" / "filtered_barcode"
    raise ValueError(
        "provide --filtered-barcode or --gexcb (or run estimated_cells first)"
    )


def resolve_jobs(args: argparse.Namespace) -> tuple[Path, list[SplitChunk]]:
    has_work_mode = args.work_path is not None
    has_single_mode = args.bam is not None
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or --bam")
    if not has_work_mode and not has_single_mode:
        raise ValueError("missing input mode: provide --work-path or --bam")

    if has_work_mode:
        work_path = Path(args.work_path)
        split_root = work_path / "split_bams"
        chunks = [
            SplitChunk(chunk_id, fwd_bam, rev_bam)
            for chunk_id, fwd_bam, rev_bam in wic.discover_bismark_sortbyname_bams(
                work_path / "align"
            )
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(
                f"no sortbyname Bismark PE BAMs found under {work_path / 'align'}"
            )
        return split_root, chunks

    bam_path = Path(args.bam)
    split_root = Path(args.output_dir) if args.output_dir else bam_path.parent.parent / "split_bams"
    chunk_id = args.chunk_id or bam_path.name.split(".", 1)[0]
    return split_root, [SplitChunk(chunk_id, bam_path, bam_path)]


def split_chunk(
    chunk: SplitChunk,
    *,
    split_root: Path,
    keep_barcodes: list[str],
    cores: int,
    dry_run: bool,
) -> None:
    strands = (
        ("forward", chunk.forward_bam),
        ("reverse", chunk.reverse_bam),
    )
    seen: set[Path] = set()
    for strand, sortbyname_bam in strands:
        if sortbyname_bam in seen:
            continue
        seen.add(sortbyname_bam)
        split_one_bam(
            sortbyname_bam,
            split_root=split_root,
            keep_barcodes=keep_barcodes,
            cores=cores,
            dry_run=dry_run,
            chunk_id=chunk.chunk_id,
            strand=strand,
        )


def main() -> int:
    args = parse_args()
    if args.cores <= 0:
        raise ValueError("cores must be > 0")

    work_path = Path(args.work_path) if args.work_path else None
    barcode_path = resolve_barcode_list(args, work_path)
    print(f"[split_bams] barcode_list={barcode_path}")

    keep_barcodes = load_barcodes(barcode_path) if not args.dry_run else []
    if args.dry_run:
        keep_barcodes = ["DRYRUN"]

    split_root, chunks = resolve_jobs(args)
    print(f"[split_bams] split_root={split_root}")
    print(f"[split_bams] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[split_bams] progress={index}/{len(chunks)}")
        split_chunk(
            chunk,
            split_root=split_root,
            keep_barcodes=keep_barcodes,
            cores=args.cores,
            dry_run=args.dry_run,
        )

    print("[split_bams] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
