#!/usr/bin/env python3
"""Per-region methylation matrices from CSR store and BED intervals."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.meth_matrix import build_region_matrices, resolve_regions_label


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build per-region methylation matrices from a MethSCAn-compatible CSR "
            "store and BED regions. Default output is sparse (matrix.mtx.gz)."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--work-path",
        help="Sample work directory; matrix at <work_path>/meth/matrix/.",
    )
    source.add_argument(
        "--data-dir",
        help="Explicit CSR matrix directory (test override).",
    )
    parser.add_argument(
        "--regions-bed",
        default=None,
        help="BED file of regions (chrom, start, end). Default: meth/vmr/vmrs.bed.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/meth/regions/<label>/.",
    )
    parser.add_argument(
        "--regions-label",
        default=None,
        help="Output subdirectory label under meth/regions/. Default: BED basename.",
    )
    parser.add_argument(
        "--dense",
        action="store_true",
        help="Write four dense cell x region .csv.gz tables instead of sparse output.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=-1,
        help="CPU threads for region aggregation. Default: all available.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and exit without writing files.",
    )
    return parser.parse_args()


def resolve_regions_bed(work_path: Path | None, regions_bed: str | None) -> Path:
    if regions_bed:
        bed_path = Path(regions_bed)
        if not bed_path.is_file():
            raise FileNotFoundError(f"regions BED not found: {bed_path}")
        return bed_path
    if work_path is None:
        raise ValueError("--regions-bed is required when using --data-dir")
    fallback = work_path / "meth" / "vmr" / "vmrs.bed"
    if fallback.is_file():
        return fallback
    raise FileNotFoundError(
        "no regions BED provided and default vmrs.bed not found at "
        f"{fallback}; pass --regions-bed explicitly"
    )


def resolve_output_dir(
    work_path: Path | None,
    regions_bed: Path,
    output_dir: str | None,
    regions_label: str | None,
) -> Path:
    if output_dir:
        return Path(output_dir)
    if work_path is None:
        raise ValueError("--output-dir is required when using --data-dir")
    label = resolve_regions_label(regions_bed, regions_label)
    return work_path / "meth" / "regions" / label


def resolve_matrix_dir(args: argparse.Namespace) -> tuple[Path, str]:
    if args.work_path is not None:
        work_path = Path(args.work_path)
        return work_path / "meth" / "matrix", "work_path"
    assert args.data_dir is not None
    return Path(args.data_dir), "data_dir"


def main() -> int:
    args = parse_args()
    if args.threads == 0 or args.threads < -1:
        raise ValueError("--threads must be -1 or >= 1")

    work_path = Path(args.work_path) if args.work_path else None
    matrix_dir, gather_mode = resolve_matrix_dir(args)
    regions_bed = resolve_regions_bed(work_path, args.regions_bed)
    output_dir = resolve_output_dir(
        work_path, regions_bed, args.output_dir, args.regions_label
    )
    label = resolve_regions_label(regions_bed, args.regions_label)
    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)

    print(f"[meth_matrix] gather_mode={gather_mode}")
    if work_path is not None:
        print(f"[meth_matrix] work_path={work_path}")
    if args.data_dir is not None:
        print(f"[meth_matrix] data_dir={args.data_dir}")
    print(f"[meth_matrix] matrix_dir={matrix_dir}")
    print(f"[meth_matrix] regions_bed={regions_bed}")
    print(f"[meth_matrix] output_dir={output_dir}")
    print(f"[meth_matrix] regions_label={label}")
    print(f"[meth_matrix] dense={int(args.dense)}")
    print(f"[meth_matrix] threads={threads}")

    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")
    smoothed_dir = matrix_dir / "smoothed"
    if not smoothed_dir.is_dir():
        raise FileNotFoundError(f"smoothed directory not found: {smoothed_dir}")

    if args.dry_run:
        print("[meth_matrix] dry_run=1")
        return 0

    outputs = build_region_matrices(
        matrix_dir,
        regions_bed,
        output_dir,
        dense=args.dense,
        threads=args.threads,
        regions_label=label,
        run_info_extra={"gather_mode": gather_mode},
    )
    print(f"[meth_matrix] output_dir={outputs['output_dir']}")
    print(f"[meth_matrix] run_info={outputs['run_info']}")
    print("[meth_matrix] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
