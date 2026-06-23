#!/usr/bin/env python3
"""Concatenate demux prefix sub-shards into one shard per barcode prefix."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

STREAMS = ("forward_1", "forward_2", "reverse_1", "reverse_2")
SUBSHARD_RE = re.compile(
    r"^(?P<readchunk>[^_]+)__(?P<prefix>[A-Z]+)\."
    r"(?P<stream>forward_1|forward_2|reverse_1|reverse_2)\.fq\.gz$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regroup demux/shards/*__<prefix>.* sub-shards into "
            "demux/<prefix>.* analysis shards."
        )
    )
    parser.add_argument(
        "--work-path",
        help="Sample work directory containing demux/shards/.",
    )
    parser.add_argument(
        "--prefix",
        help="Regroup a single barcode prefix only (for Slurm per-prefix jobs).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned outputs without reading or writing files.",
    )
    return parser.parse_args()


def discover_subshards(shard_dir: Path) -> dict[str, dict[str, list[Path]]]:
    """Return {prefix: {stream: [sorted sub-shard paths]}}."""
    if not shard_dir.is_dir():
        return {}

    grouped: dict[str, dict[str, list[Path]]] = defaultdict(lambda: defaultdict(list))
    for path in sorted(shard_dir.glob("*__*.fq.gz")):
        match = SUBSHARD_RE.match(path.name)
        if not match:
            continue
        prefix = match.group("prefix")
        stream = match.group("stream")
        if stream not in STREAMS:
            continue
        grouped[prefix][stream].append(path)

    return {prefix: dict(streams) for prefix, streams in sorted(grouped.items())}


def output_path(demux_dir: Path, prefix: str, stream: str) -> Path:
    return demux_dir / f"{prefix}.{stream}.fq.gz"


def concat_gzip_members(sources: list[Path], destination: Path) -> None:
    """Append gzip members (binary concat is valid multi-member gzip)."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as out_handle:
        for source in sources:
            with source.open("rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)


def write_chunks_manifest(
    demux_dir: Path,
    grouped: dict[str, dict[str, list[Path]]],
    *,
    dry_run: bool,
) -> Path:
    manifest = demux_dir / "chunks.tsv"
    lines = ["chunk_id\tprefix\tstream\tsubshard_count\tsubshard_paths"]
    for prefix, streams in sorted(grouped.items()):
        for stream in STREAMS:
            sources = streams.get(stream, [])
            if not sources:
                continue
            paths = ";".join(str(p) for p in sources)
            lines.append(f"{prefix}\t{prefix}\t{stream}\t{len(sources)}\t{paths}")
    content = "\n".join(lines) + "\n"
    print(f"[regroup_shards] manifest={manifest}")
    if not dry_run:
        manifest.write_text(content, encoding="utf-8")
    return manifest


def regroup_prefix(
    demux_dir: Path,
    prefix: str,
    streams: dict[str, list[Path]],
    *,
    dry_run: bool,
) -> list[Path]:
    outputs: list[Path] = []
    for stream in STREAMS:
        sources = streams.get(stream, [])
        if not sources:
            continue
        dest = output_path(demux_dir, prefix, stream)
        print(
            f"[regroup_shards] prefix={prefix} stream={stream} "
            f"sources={len(sources)} output={dest}"
        )
        if not dry_run:
            concat_gzip_members(sources, dest)
        outputs.append(dest)
    return outputs


def main() -> int:
    args = parse_args()
    if not args.work_path:
        raise ValueError("--work-path is required")

    demux_dir = Path(args.work_path) / "demux"
    shard_dir = demux_dir / "shards"
    grouped = discover_subshards(shard_dir)

    print(f"[regroup_shards] demux_dir={demux_dir}")
    print(f"[regroup_shards] shard_dir={shard_dir}")
    print(f"[regroup_shards] prefix_count={len(grouped)}")

    if not grouped:
        raise ValueError(f"no sub-shards found under {shard_dir}")

    if args.prefix is not None:
        if args.prefix not in grouped:
            raise ValueError(f"prefix not found in sub-shards: {args.prefix}")
        grouped = {args.prefix: grouped[args.prefix]}

    total_outputs = 0
    for prefix, streams in sorted(grouped.items()):
        outputs = regroup_prefix(demux_dir, prefix, streams, dry_run=args.dry_run)
        total_outputs += len(outputs)

    if args.prefix is None:
        write_chunks_manifest(demux_dir, grouped, dry_run=args.dry_run)

    print(f"[regroup_shards] outputs={total_outputs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
