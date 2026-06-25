#!/usr/bin/env python3
"""Sliding-window VMR scan over smoothed methylation matrices."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.meth_matrix import scan_vmrs

DEFAULT_BANDWIDTH = 2000
DEFAULT_STEPSIZE = 100
DEFAULT_VAR_THRESHOLD = 0.02
DEFAULT_MIN_CELLS = 6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan the genome for variably methylated regions (VMRs) using "
            "smoothed CSR matrices under <work_path>/meth/matrix/."
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
        "--output",
        default=None,
        help="Output BED path. Default: <work_path>/meth/vmr/vmrs.bed.",
    )
    parser.add_argument(
        "--bandwidth",
        type=int,
        default=DEFAULT_BANDWIDTH,
        help="Sliding-window bandwidth in bp. Default: 2000.",
    )
    parser.add_argument(
        "--stepsize",
        type=int,
        default=DEFAULT_STEPSIZE,
        help="Sliding-window step size in bp. Default: 100.",
    )
    parser.add_argument(
        "--var-threshold",
        type=float,
        default=DEFAULT_VAR_THRESHOLD,
        help="Top fraction of variable windows to merge as VMRs. Default: 0.02.",
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=DEFAULT_MIN_CELLS,
        help="Minimum cells with coverage to report a VMR. Default: 6.",
    )
    parser.add_argument(
        "--bridge-gaps",
        type=int,
        default=0,
        help="Merge VMRs within this gap distance (bp). Default: 0 (off).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=-1,
        help="CPU threads for window scan. Default: all available.",
    )
    parser.add_argument(
        "--write-header",
        action="store_true",
        help="Write a header row to the output BED.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and exit without writing files.",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, str]:
    if args.work_path is not None:
        work_path = Path(args.work_path)
        matrix_dir = work_path / "meth" / "matrix"
        output_bed = (
            Path(args.output) if args.output else work_path / "meth" / "vmr" / "vmrs.bed"
        )
        return matrix_dir, output_bed, "work_path"
    assert args.data_dir is not None
    matrix_dir = Path(args.data_dir)
    if args.output is None:
        raise ValueError("--output is required when using --data-dir")
    return matrix_dir, Path(args.output), "data_dir"


def main() -> int:
    args = parse_args()
    if args.bandwidth < 1:
        raise ValueError("--bandwidth must be >= 1")
    if args.stepsize < 1:
        raise ValueError("--stepsize must be >= 1")
    if not 0 <= args.var_threshold <= 1:
        raise ValueError("--var-threshold must be between 0 and 1")
    if args.min_cells < 1:
        raise ValueError("--min-cells must be >= 1")
    if args.bridge_gaps < 0:
        raise ValueError("--bridge-gaps must be >= 0")

    matrix_dir, output_bed, gather_mode = resolve_paths(args)
    smoothed_dir = matrix_dir / "smoothed"
    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)

    print(f"[meth_scan] gather_mode={gather_mode}")
    if args.work_path is not None:
        print(f"[meth_scan] work_path={args.work_path}")
    if args.data_dir is not None:
        print(f"[meth_scan] data_dir={args.data_dir}")
    print(f"[meth_scan] matrix_dir={matrix_dir}")
    print(f"[meth_scan] smoothed_dir={smoothed_dir}")
    print(f"[meth_scan] output_bed={output_bed}")
    print(f"[meth_scan] bandwidth={args.bandwidth}")
    print(f"[meth_scan] stepsize={args.stepsize}")
    print(f"[meth_scan] var_threshold={args.var_threshold}")
    print(f"[meth_scan] min_cells={args.min_cells}")
    print(f"[meth_scan] bridge_gaps={args.bridge_gaps}")
    print(f"[meth_scan] threads={threads}")

    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")
    if not smoothed_dir.is_dir():
        raise FileNotFoundError(f"smoothed directory not found: {smoothed_dir}")
    smoothed_files = list(smoothed_dir.glob("*.csv.gz")) + list(smoothed_dir.glob("*.csv"))
    if not smoothed_files:
        raise FileNotFoundError(f"no smoothed chromosome files under {smoothed_dir}")

    if args.dry_run:
        print("[meth_scan] dry_run=1")
        return 0

    outputs = scan_vmrs(
        matrix_dir,
        output_bed,
        bandwidth=args.bandwidth,
        stepsize=args.stepsize,
        var_threshold=args.var_threshold,
        min_cells=args.min_cells,
        bridge_gaps=args.bridge_gaps,
        threads=threads,
        write_header=args.write_header,
        run_info_extra={"gather_mode": gather_mode},
    )
    print(f"[meth_scan] output_bed={outputs['output_bed']}")
    print(f"[meth_scan] run_info={outputs['run_info']}")
    print("[meth_scan] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
