#!/usr/bin/env python3
"""Name-sort Bismark PE BAMs per chunk (samtools sort -n)."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import workflow_input_checks as wic


@dataclass(frozen=True)
class SortChunk:
    chunk_id: str
    forward_bam: Path
    reverse_bam: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Name-sort Bismark paired-end BAMs with samtools sort -n "
            "(DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help=(
            "Sample work directory. Scans align/*.forward_1_bismark_bt2_pe.bam "
            "and writes sortbyname BAMs alongside inputs."
        ),
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument(
        "--forward-bam",
        help="Forward BAM for single-chunk mode.",
    )
    parser.add_argument(
        "--reverse-bam",
        help="Reverse BAM for single-chunk mode.",
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for single-chunk mode (default: align/ under work-path).",
    )
    parser.add_argument(
        "--sort-threads",
        type=int,
        default=6,
        help="Thread count for samtools sort -@. Default: 6.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or name. Default: samtools.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and expected outputs without running samtools.",
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


def sortbyname_output_path(bam_path: Path) -> Path:
    return bam_path.parent / f"{bam_path.stem}_sortbyname.bam"


def build_sort_cmd(
    *,
    samtools_bin: str,
    sort_threads: int,
    input_bam: Path,
    output_bam: Path,
) -> list[str]:
    return [
        samtools_bin,
        "sort",
        "-n",
        "-@",
        str(sort_threads),
        "-o",
        str(output_bam),
        str(input_bam),
    ]


def should_skip_sort(input_bam: Path, output_bam: Path) -> bool:
    if not output_bam.is_file():
        return False
    if not input_bam.is_file():
        return False
    return output_bam.stat().st_mtime >= input_bam.stat().st_mtime


def resolve_jobs(args: argparse.Namespace) -> tuple[Path, list[SortChunk]]:
    has_work_mode = args.work_path is not None
    has_single_mode = any(
        value is not None
        for value in (args.forward_bam, args.reverse_bam, args.output_dir)
    )
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or explicit BAM paths")
    if not has_work_mode and not has_single_mode:
        raise ValueError(
            "missing input mode: provide --work-path or "
            "--forward-bam/--reverse-bam/--output-dir"
        )

    if has_work_mode:
        work_path = Path(args.work_path)
        align_dir = work_path / "align"
        chunks = [
            SortChunk(chunk_id, fwd_bam, rev_bam)
            for chunk_id, fwd_bam, rev_bam in wic.discover_bismark_pe_bams(align_dir)
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(f"no Bismark PE BAMs found under {align_dir}")
        return align_dir, chunks

    if not all(
        value is not None
        for value in (args.forward_bam, args.reverse_bam, args.output_dir)
    ):
        raise ValueError(
            "single-chunk mode requires --forward-bam, --reverse-bam, and --output-dir"
        )
    chunk_id = args.chunk_id or Path(args.forward_bam).name.split(".", 1)[0]
    chunk = SortChunk(
        chunk_id,
        Path(args.forward_bam),
        Path(args.reverse_bam),
    )
    return Path(args.output_dir), [chunk]


def sort_bam(
    *,
    input_bam: Path,
    samtools_bin: str,
    sort_threads: int,
    dry_run: bool,
    chunk_id: str,
    strand: str,
) -> None:
    output_bam = sortbyname_output_path(input_bam)
    cmd = build_sort_cmd(
        samtools_bin=samtools_bin,
        sort_threads=sort_threads,
        input_bam=input_bam,
        output_bam=output_bam,
    )
    print(f"[bam_sort] chunk={chunk_id} strand={strand}")
    print(f"[bam_sort] input_bam={input_bam}")
    print(f"[bam_sort] output_bam={output_bam}")
    print(f"[bam_sort] command={quoted(cmd)}")
    if dry_run:
        return
    if should_skip_sort(input_bam, output_bam):
        print("[bam_sort] skipped=1")
        return
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


def sort_chunk(
    chunk: SortChunk,
    *,
    samtools_bin: str,
    sort_threads: int,
    dry_run: bool,
) -> None:
    strands = (
        ("forward", chunk.forward_bam),
        ("reverse", chunk.reverse_bam),
    )
    for strand, input_bam in strands:
        sort_bam(
            input_bam=input_bam,
            samtools_bin=samtools_bin,
            sort_threads=sort_threads,
            dry_run=dry_run,
            chunk_id=chunk.chunk_id,
            strand=strand,
        )


def main() -> int:
    args = parse_args()
    if args.sort_threads <= 0:
        raise ValueError("sort_threads must be > 0")

    samtools_bin = resolve_samtools_bin(args.samtools_bin)
    wic.require_optional_executable_path("samtools_bin", samtools_bin)

    output_dir, chunks = resolve_jobs(args)
    print(f"[bam_sort] output_dir={output_dir}")
    print(f"[bam_sort] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[bam_sort] progress={index}/{len(chunks)}")
        sort_chunk(
            chunk,
            samtools_bin=samtools_bin,
            sort_threads=args.sort_threads,
            dry_run=args.dry_run,
        )

    print("[bam_sort] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
