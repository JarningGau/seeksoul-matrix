#!/usr/bin/env python3
"""Per-cell methylation profiles around user-supplied BED features."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.meth_matrix import build_profile, resolve_regions_label
from lib.meth_matrix.profile import DEFAULT_WIDTH


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute the average per-cell methylation profile around BED "
            "features from a CSR store under <work_path>/meth/matrix/."
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
        required=True,
        help="BED file of features (chrom, start, end). Required.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/meth/profile/<label>/.",
    )
    parser.add_argument(
        "--regions-label",
        default=None,
        help="Output subdirectory label under meth/profile/. Default: BED basename.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=DEFAULT_WIDTH,
        help=(
            "Profile width in bp. Each feature is centered and extended or "
            f"trimmed to this width. Default: {DEFAULT_WIDTH}."
        ),
    )
    parser.add_argument(
        "--strand-column",
        type=int,
        default=None,
        help=(
            "1-indexed BED column with strand (+, -, or .). "
            "Minus-strand features are flipped so upstream is negative."
        ),
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Optional constant column added to the output profile.csv.gz.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=-1,
        help="CPU threads for profile accumulation. Default: all available.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and exit without writing files.",
    )
    return parser.parse_args()


def resolve_matrix_dir(args: argparse.Namespace) -> tuple[Path, str]:
    if args.work_path is not None:
        return Path(args.work_path) / "meth" / "matrix", "work_path"
    assert args.data_dir is not None
    return Path(args.data_dir), "data_dir"


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
    return work_path / "meth" / "profile" / label


def main() -> int:
    args = parse_args()
    if args.width < 1:
        raise ValueError("--width must be >= 1")
    if args.strand_column is not None and args.strand_column < 1:
        raise ValueError("--strand-column must be >= 1")
    if args.threads == 0 or args.threads < -1:
        raise ValueError("--threads must be -1 or >= 1")

    work_path = Path(args.work_path) if args.work_path else None
    matrix_dir, gather_mode = resolve_matrix_dir(args)
    regions_bed = Path(args.regions_bed)
    if not regions_bed.is_file():
        raise FileNotFoundError(f"regions BED not found: {regions_bed}")
    output_dir = resolve_output_dir(
        work_path, regions_bed, args.output_dir, args.regions_label
    )
    label = resolve_regions_label(regions_bed, args.regions_label)
    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)
    column_label = args.label.strip() if args.label else ""

    print(f"[meth_profile] gather_mode={gather_mode}")
    if work_path is not None:
        print(f"[meth_profile] work_path={work_path}")
    if args.data_dir is not None:
        print(f"[meth_profile] data_dir={args.data_dir}")
    print(f"[meth_profile] matrix_dir={matrix_dir}")
    print(f"[meth_profile] regions_bed={regions_bed}")
    print(f"[meth_profile] output_dir={output_dir}")
    print(f"[meth_profile] regions_label={label}")
    print(f"[meth_profile] profile_csv={output_dir / 'profile.csv.gz'}")
    print(f"[meth_profile] width={args.width}")
    print(f"[meth_profile] strand_column={args.strand_column}")
    print(f"[meth_profile] label={column_label}")
    print(f"[meth_profile] threads={threads}")

    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")
    if not (matrix_dir / "column_header.txt").is_file():
        raise FileNotFoundError(
            f"column_header.txt not found under {matrix_dir}"
        )
    npz_files = list(matrix_dir.glob("*.npz"))
    if not npz_files:
        raise FileNotFoundError(f"no CSR matrix files found under {matrix_dir}")

    if args.dry_run:
        print("[meth_profile] dry_run=1")
        return 0

    outputs = build_profile(
        matrix_dir,
        regions_bed,
        output_dir,
        width=args.width,
        strand_column=args.strand_column,
        label=column_label or None,
        threads=args.threads,
        regions_label=label,
        run_info_extra={"gather_mode": gather_mode},
    )
    print(f"[meth_profile] profile_csv={outputs['profile_csv']}")
    print(f"[meth_profile] run_info={outputs['run_info']}")
    print("[meth_profile] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
