#!/usr/bin/env python3
"""Merge per-chunk barcode counts and filter called cells (methylation-only path)."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate *_cb_aligned_reads_counts.csv files and filter barcodes "
            "by expected cell count (SeekSoulMethyl step3_estimated_cells)."
        )
    )
    parser.add_argument(
        "--align-dir",
        help="Directory containing *_cb_aligned_reads_counts.csv files.",
    )
    parser.add_argument(
        "--work-path",
        help="Sample work directory; reads align/ and writes cells/ outputs.",
    )
    parser.add_argument(
        "--expected-cell-num",
        type=int,
        default=3000,
        help="Expected number of cells for 99th-percentile threshold. Default: 3000.",
    )
    parser.add_argument(
        "--force-cell-num",
        type=int,
        default=None,
        help=(
            "Take top N barcodes by aligned_reads (nonzero only). "
            "When set, overrides --expected-cell-num threshold filtering."
        ),
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory (default: work-path/cells or align-dir parent/cells).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths without reading or writing files.",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    if args.work_path:
        work_path = Path(args.work_path)
        align_dir = work_path / "align"
        output_dir = Path(args.output_dir) if args.output_dir else work_path / "cells"
        return align_dir, output_dir
    if args.align_dir:
        align_dir = Path(args.align_dir)
        output_dir = (
            Path(args.output_dir)
            if args.output_dir
            else align_dir.parent / "cells"
        )
        return align_dir, output_dir
    raise ValueError("provide --work-path or --align-dir")


def load_barcode_totals(align_dir: Path) -> dict[str, int]:
    csv_files = sorted(align_dir.glob("*_cb_aligned_reads_counts.csv"))
    if not csv_files:
        raise ValueError(
            f"no *_cb_aligned_reads_counts.csv files found under {align_dir}"
        )

    barcode_totals: dict[str, int] = defaultdict(int)
    for csv_path in csv_files:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "barcode" not in reader.fieldnames:
                raise ValueError(f"missing barcode column in {csv_path}")
            if "aligned_reads" not in reader.fieldnames:
                raise ValueError(f"missing aligned_reads column in {csv_path}")
            for row in reader:
                barcode = row["barcode"]
                barcode_totals[barcode] += int(row["aligned_reads"])
    return dict(barcode_totals)


def filter_barcodes(
    barcode_totals: dict[str, int],
    expected_cell_num: int,
) -> list[tuple[str, int]]:
    if expected_cell_num <= 0:
        raise ValueError("expected_cell_num must be > 0")
    if not barcode_totals:
        return []

    merged_rows = sorted(
        barcode_totals.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    percentile = 99
    threshold_index = int(expected_cell_num * (1 - percentile / 100))
    if threshold_index >= len(merged_rows):
        threshold_index = len(merged_rows) - 1
    readscut = merged_rows[threshold_index][1] * 0.1
    return [(barcode, reads) for barcode, reads in merged_rows if reads > readscut]


def force_filter_barcodes(
    barcode_totals: dict[str, int],
    force_cell_num: int,
) -> list[tuple[str, int]]:
    if force_cell_num <= 0:
        raise ValueError("force_cell_num must be > 0")

    nonzero_rows = [
        (barcode, reads)
        for barcode, reads in barcode_totals.items()
        if reads > 0
    ]
    if not nonzero_rows:
        return []

    ranked_rows = sorted(
        nonzero_rows,
        key=lambda item: (-item[1], item[0]),
    )
    return ranked_rows[:force_cell_num]


def write_outputs(
    output_dir: Path,
    barcode_totals: dict[str, int],
    filtered_rows: list[tuple[str, int]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_path = output_dir / "merged_barcode_counts.csv"
    merged_rows = sorted(
        barcode_totals.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    with merged_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "aligned_reads"])
        for barcode, aligned_reads in merged_rows:
            writer.writerow([barcode, aligned_reads])

    filtered_counts_path = output_dir / "filtered_barcode_read_counts.csv"
    with filtered_counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["aligned_reads", "barcode"])
        for barcode, aligned_reads in filtered_rows:
            writer.writerow([aligned_reads, barcode])

    filtered_barcode_path = output_dir / "filtered_barcode"
    with filtered_barcode_path.open("w", encoding="utf-8", newline="") as handle:
        for barcode, _aligned_reads in filtered_rows:
            handle.write(f"{barcode}\n")


def main() -> int:
    args = parse_args()
    align_dir, output_dir = resolve_paths(args)
    print(f"[estimated_cells] align_dir={align_dir}")
    print(f"[estimated_cells] output_dir={output_dir}")
    if args.force_cell_num is not None:
        print(f"[estimated_cells] force_cell_num={args.force_cell_num}")
    else:
        print(f"[estimated_cells] expected_cell_num={args.expected_cell_num}")

    if args.dry_run:
        csv_files = list(align_dir.glob("*_cb_aligned_reads_counts.csv"))
        print(f"[estimated_cells] input_csv_count={len(csv_files)}")
        print("[estimated_cells] done")
        return 0

    barcode_totals = load_barcode_totals(align_dir)
    if args.force_cell_num is not None:
        filtered_rows = force_filter_barcodes(barcode_totals, args.force_cell_num)
    else:
        filtered_rows = filter_barcodes(barcode_totals, args.expected_cell_num)
    write_outputs(output_dir, barcode_totals, filtered_rows)
    print(f"[estimated_cells] total_barcodes={len(barcode_totals)}")
    print(f"[estimated_cells] filtered_barcodes={len(filtered_rows)}")
    print("[estimated_cells] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
