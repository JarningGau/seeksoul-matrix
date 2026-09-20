#!/usr/bin/env python3
"""Sum fastp reports from multiple lanes into one shard_fastq/fastp.json.

qc_summary and saturation read a single work/<sample>/shard_fastq/fastp.json.
This glue script does not change those stages: it writes a minimal report with
summed total_reads / total_bases and a read-weighted duplication rate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge lane fastp JSON reports by summing read and base counts."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more fastp JSON files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Merged fastp JSON path.",
    )
    return parser.parse_args()


def load_counts(path: Path) -> tuple[int, int, int, int, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    summary = payload.get("summary") or {}
    after = summary.get("after_filtering") or {}
    before = summary.get("before_filtering") or {}
    after_reads = int(after.get("total_reads") or before.get("total_reads") or 0)
    after_bases = int(after.get("total_bases") or before.get("total_bases") or 0)
    before_reads = int(before.get("total_reads") or after.get("total_reads") or 0)
    before_bases = int(before.get("total_bases") or after.get("total_bases") or 0)
    duplication = payload.get("duplication") or {}
    dup_rate = float(duplication.get("rate") or 0.0)
    return after_reads, after_bases, before_reads, before_bases, dup_rate


def main() -> int:
    args = parse_args()
    after_reads = after_bases = before_reads = before_bases = 0
    weighted_dup = 0.0
    sources: list[str] = []

    for path in args.inputs:
        if not path.is_file():
            raise FileNotFoundError(f"fastp JSON not found: {path}")
        a_reads, a_bases, b_reads, b_bases, dup_rate = load_counts(path)
        after_reads += a_reads
        after_bases += a_bases
        before_reads += b_reads
        before_bases += b_bases
        weighted_dup += dup_rate * a_reads
        sources.append(str(path))

    dup_rate = (weighted_dup / after_reads) if after_reads else 0.0
    merged = {
        "summary": {
            "before_filtering": {
                "total_reads": before_reads,
                "total_bases": before_bases,
            },
            "after_filtering": {
                "total_reads": after_reads,
                "total_bases": after_bases,
            },
        },
        "duplication": {"rate": dup_rate},
        "seeksoul_matrix_merge": {"sources": sources},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(merged, indent=2) + "\n", encoding="utf-8")
    print(f"[merge_fastp_json] output={args.output}")
    print(f"[merge_fastp_json] lanes={len(sources)}")
    print(f"[merge_fastp_json] total_reads={after_reads}")
    print(f"[merge_fastp_json] total_bases={after_bases}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
