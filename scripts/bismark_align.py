#!/usr/bin/env python3
"""Run seekgene Bismark on demux forward/reverse paired FASTQ per chunk."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import workflow_input_checks as wic


@dataclass(frozen=True)
class AlignChunk:
    chunk_id: str
    forward_r1: Path
    forward_r2: Path
    reverse_r1: Path
    reverse_r2: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align demux forward/reverse FASTQ pairs with seekgene Bismark "
            "(DD-MET5 workflow)."
        )
    )
    parser.add_argument(
        "--work-path",
        help=(
            "Sample work directory. Scans demux/*.forward_1.fq.gz and writes "
            "outputs to <work-path>/align/."
        ),
    )
    parser.add_argument("--chunk-id", help="Optional chunk filter in work-path mode.")
    parser.add_argument("--forward-r1", help="Forward R1 FASTQ for single-chunk mode.")
    parser.add_argument("--forward-r2", help="Forward R2 FASTQ for single-chunk mode.")
    parser.add_argument("--reverse-r1", help="Reverse R1 FASTQ for single-chunk mode.")
    parser.add_argument("--reverse-r2", help="Reverse R2 FASTQ for single-chunk mode.")
    parser.add_argument(
        "--output-dir",
        help="Output directory for single-chunk mode (default: align/ under work-path).",
    )
    parser.add_argument(
        "--bismark-ref",
        required=True,
        help="Bismark --genome path (parent of Bisulfite_Genome/).",
    )
    parser.add_argument(
        "--bismark-parallel",
        type=int,
        default=8,
        help="Bismark --parallel value. Default: 8.",
    )
    parser.add_argument(
        "--bismark-max-insert",
        type=int,
        default=1000,
        help="Bismark -X max insert size. Default: 1000.",
    )
    parser.add_argument(
        "--bismark-bin",
        default="bismark",
        help="bismark executable path or name. Default: bismark.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and expected outputs without running Bismark.",
    )
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def resolve_bismark_bin(raw: str) -> str:
    if raw != "bismark":
        return raw
    candidate = Path(sys.executable).resolve().parent / "bismark"
    if candidate.is_file():
        return str(candidate)
    return raw


def ensure_seekgene_bismark(bismark_bin: str) -> None:
    try:
        completed = subprocess.run(
            [bismark_bin, "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise ValueError(
            f"bismark executable not found: {bismark_bin} "
            "(run 'pixi run setup-bismark')"
        ) from exc
    help_text = (completed.stdout or "") + (completed.stderr or "")
    if "--add_barcode" not in help_text or "--add_umi" not in help_text:
        raise ValueError(
            "bismark on PATH is not the seekgene fork (missing --add_barcode/--add_umi); "
            "run 'pixi run setup-bismark'"
        )


def bismark_output_bam(r1_path: Path) -> str:
    stem = r1_path.name
    for suffix in (".fq.gz", ".fastq.gz"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return f"{stem}_bismark_bt2_pe.bam"


def bismark_output_report(r1_path: Path) -> str:
    stem = r1_path.name
    for suffix in (".fq.gz", ".fastq.gz"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return f"{stem}_bismark_bt2_PE_report.txt"


def build_bismark_cmd(
    *,
    bismark_bin: str,
    bismark_ref: str,
    bismark_parallel: int,
    bismark_max_insert: int,
    r1_path: Path,
    r2_path: Path,
    output_dir: Path,
    temp_dir: Path,
    reverse: bool,
) -> list[str]:
    cmd = [
        bismark_bin,
        "--genome",
        bismark_ref,
        "--parallel",
        str(bismark_parallel),
        "-1",
        str(r1_path),
        "-2",
        str(r2_path),
        "-o",
        str(output_dir),
        "-X",
        str(bismark_max_insert),
        "--add_barcode",
        "--add_umi",
        "--temp_dir",
        str(temp_dir),
    ]
    if reverse:
        cmd.insert(cmd.index("-X"), "--pbat")
    return cmd


def resolve_jobs(args: argparse.Namespace) -> tuple[Path, list[AlignChunk]]:
    has_work_mode = args.work_path is not None
    has_single_mode = any(
        value is not None
        for value in (
            args.forward_r1,
            args.forward_r2,
            args.reverse_r1,
            args.reverse_r2,
            args.output_dir,
        )
    )
    if has_work_mode and has_single_mode:
        raise ValueError("use either --work-path mode or explicit FASTQ paths")
    if not has_work_mode and not has_single_mode:
        raise ValueError(
            "missing input mode: provide --work-path or "
            "--forward-r1/--forward-r2/--reverse-r1/--reverse-r2/--output-dir"
        )

    if has_work_mode:
        work_path = Path(args.work_path)
        demux_dir = work_path / "demux"
        output_dir = work_path / "align"
        chunks = [
            AlignChunk(chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2)
            for chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2 in wic.discover_demux_align_chunks(
                demux_dir
            )
        ]
        if args.chunk_id:
            chunks = [chunk for chunk in chunks if chunk.chunk_id == args.chunk_id]
        if not chunks:
            raise ValueError(f"no demux align chunks found under {demux_dir}")
        return output_dir, chunks

    if not all(
        value is not None
        for value in (
            args.forward_r1,
            args.forward_r2,
            args.reverse_r1,
            args.reverse_r2,
            args.output_dir,
        )
    ):
        raise ValueError(
            "single-chunk mode requires --forward-r1, --forward-r2, "
            "--reverse-r1, --reverse-r2, and --output-dir"
        )
    chunk_id = args.chunk_id or Path(args.forward_r1).name.split(".", 1)[0]
    chunk = AlignChunk(
        chunk_id,
        Path(args.forward_r1),
        Path(args.forward_r2),
        Path(args.reverse_r1),
        Path(args.reverse_r2),
    )
    return Path(args.output_dir), [chunk]


def align_chunk(
    chunk: AlignChunk,
    *,
    output_dir: Path,
    bismark_bin: str,
    bismark_ref: str,
    bismark_parallel: int,
    bismark_max_insert: int,
    dry_run: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    strands = (
        ("forward", chunk.forward_r1, chunk.forward_r2, False),
        ("reverse", chunk.reverse_r1, chunk.reverse_r2, True),
    )
    for strand, r1_path, r2_path, is_reverse in strands:
        temp_dir = output_dir / f"tmp_{chunk.chunk_id}_{strand}"
        cmd = build_bismark_cmd(
            bismark_bin=bismark_bin,
            bismark_ref=bismark_ref,
            bismark_parallel=bismark_parallel,
            bismark_max_insert=bismark_max_insert,
            r1_path=r1_path,
            r2_path=r2_path,
            output_dir=output_dir,
            temp_dir=temp_dir,
            reverse=is_reverse,
        )
        bam_name = bismark_output_bam(r1_path)
        report_name = bismark_output_report(r1_path)
        print(f"[bismark_align] chunk={chunk.chunk_id} strand={strand}")
        print(f"[bismark_align] input_r1={r1_path}")
        print(f"[bismark_align] input_r2={r2_path}")
        print(f"[bismark_align] output_bam={output_dir / bam_name}")
        print(f"[bismark_align] output_report={output_dir / report_name}")
        print(f"[bismark_align] command={quoted(cmd)}")
        if dry_run:
            continue
        temp_dir.mkdir(parents=True, exist_ok=True)
        log_path = output_dir / f"bismark_{chunk.chunk_id}_{strand}.log"
        with log_path.open("w", encoding="utf-8") as log_handle:
            subprocess.run(cmd, check=True, stdout=log_handle, stderr=subprocess.STDOUT)


def main() -> int:
    args = parse_args()
    bismark_bin = resolve_bismark_bin(args.bismark_bin)
    ensure_seekgene_bismark(bismark_bin)
    wic.require_bismark_ref(wic.resolve_config_path(args.bismark_ref))

    output_dir, chunks = resolve_jobs(args)
    print(f"[bismark_align] output_dir={output_dir}")
    print(f"[bismark_align] chunk_count={len(chunks)}")

    for index, chunk in enumerate(chunks, start=1):
        if len(chunks) > 1:
            print(f"[bismark_align] progress={index}/{len(chunks)}")
        align_chunk(
            chunk,
            output_dir=output_dir,
            bismark_bin=bismark_bin,
            bismark_ref=str(wic.resolve_config_path(args.bismark_ref)),
            bismark_parallel=args.bismark_parallel,
            bismark_max_insert=args.bismark_max_insert,
            dry_run=args.dry_run,
        )

    print("[bismark_align] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
