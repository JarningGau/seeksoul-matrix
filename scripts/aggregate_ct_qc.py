#!/usr/bin/env python3
"""Aggregate per-chunk linker.tsv files into sample-level qc.CtoT.tsv."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge chunk linker.tsv files and summarize C→T by corrected CB."
    )
    parser.add_argument(
        "--demux-dir",
        required=True,
        help="Directory containing <chunk>.linker.tsv files.",
    )
    parser.add_argument(
        "--output",
        help="Output TSV path. Default: <demux-dir>/qc.CtoT.tsv.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print output path without reading or writing files.",
    )
    return parser.parse_args()


def round_ctot(c_sum: int, t_sum: int) -> float:
    return round(t_sum / (c_sum + t_sum), 3)


def aggregate_linker_files(linker_files: list[Path]) -> dict[str, tuple[int, int]]:
    totals: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for path in linker_files:
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header.startswith("CR\t"):
                raise ValueError(f"unexpected linker.tsv header in {path}")
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != 4:
                    raise ValueError(f"invalid linker.tsv row in {path}: {line}")
                cr, _ub, c_str, t_str = parts
                totals[cr][0] += int(c_str)
                totals[cr][1] += int(t_str)
    return {cr: (vals[0], vals[1]) for cr, vals in totals.items()}


def write_qc_ctot(path: Path, totals: dict[str, tuple[int, int]]) -> int:
    rows_written = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write("CR\tC\tT\tCtoT\n")
        for cr in sorted(totals):
            c_sum, t_sum = totals[cr]
            if c_sum + t_sum == 0:
                continue
            handle.write(f"{cr}\t{c_sum}\t{t_sum}\t{round_ctot(c_sum, t_sum):.3f}\n")
            rows_written += 1
    return rows_written


def main() -> int:
    args = parse_args()
    demux_dir = Path(args.demux_dir)
    output_path = Path(args.output) if args.output else demux_dir / "qc.CtoT.tsv"
    linker_files = sorted(demux_dir.glob("*.linker.tsv"))

    print(f"[aggregate_ct_qc] demux_dir={demux_dir}")
    print(f"[aggregate_ct_qc] linker_count={len(linker_files)}")
    print(f"[aggregate_ct_qc] output={output_path}")

    if args.dry_run:
        return 0

    if not linker_files:
        raise ValueError(f"no *.linker.tsv files found under {demux_dir}")

    totals = aggregate_linker_files(linker_files)
    rows = write_qc_ctot(output_path, totals)
    print(f"[aggregate_ct_qc] barcodes={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
