#!/usr/bin/env python3
"""Tricube pseudobulk smoothing over a MethSCAn-compatible CSR matrix store."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.meth_matrix import smooth_matrix_store

DEFAULT_BANDWIDTH = 1000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Smooth per-position methylation fractions over a CSR matrix store "
            "under <work_path>/meth/matrix/smoothed/."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--work-path",
        help="Sample work directory; reads <work_path>/meth/matrix/.",
    )
    source.add_argument(
        "--data-dir",
        help="Explicit CSR matrix directory (test override).",
    )
    parser.add_argument(
        "--bandwidth",
        type=int,
        default=DEFAULT_BANDWIDTH,
        help="Tricube smoothing bandwidth in bp. Default: 1000.",
    )
    parser.add_argument(
        "--use-weights",
        action="store_true",
        help="Weight sites by log1p(coverage). Default: off.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and exit without writing files.",
    )
    return parser.parse_args()


def resolve_matrix_dir(args: argparse.Namespace) -> tuple[Path, str]:
    if args.work_path is not None:
        work_path = Path(args.work_path)
        return work_path / "meth" / "matrix", "work_path"
    assert args.data_dir is not None
    return Path(args.data_dir), "data_dir"


def main() -> int:
    args = parse_args()
    if args.bandwidth < 1:
        raise ValueError("--bandwidth must be >= 1")

    matrix_dir, gather_mode = resolve_matrix_dir(args)
    smoothed_dir = matrix_dir / "smoothed"

    print(f"[meth_smooth] gather_mode={gather_mode}")
    if args.work_path is not None:
        print(f"[meth_smooth] work_path={args.work_path}")
    if args.data_dir is not None:
        print(f"[meth_smooth] data_dir={args.data_dir}")
    print(f"[meth_smooth] matrix_dir={matrix_dir}")
    print(f"[meth_smooth] smoothed_dir={smoothed_dir}")
    print(f"[meth_smooth] bandwidth={args.bandwidth}")
    print(f"[meth_smooth] use_weights={int(args.use_weights)}")

    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")

    npz_files = sorted(matrix_dir.glob("*.npz"))
    print(f"[meth_smooth] chrom_count={len(npz_files)}")
    for npz_path in npz_files[:5]:
        print(f"[meth_smooth] sample_chrom={npz_path.name}")

    if args.dry_run:
        print("[meth_smooth] dry_run=1")
        return 0

    outputs = smooth_matrix_store(
        matrix_dir,
        bandwidth=args.bandwidth,
        use_weights=args.use_weights,
        run_info_extra={"gather_mode": gather_mode},
    )
    print(f"[meth_smooth] smoothed_dir={outputs['smoothed_dir']}")
    print(f"[meth_smooth] run_info={outputs['run_info']}")
    print("[meth_smooth] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
