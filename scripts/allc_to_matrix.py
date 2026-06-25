#!/usr/bin/env python3
"""Convert per-cell ALLC files into a MethSCAn-compatible CSR matrix store."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.meth_matrix import barcode_from_allc_path, build_matrix_store

DEFAULT_METH_CONTEXT = "CG"
DEFAULT_CHUNKSIZE = 10_000_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Gather per-cell ALLC from seeksoul-matrix contract paths and write "
            "a MethSCAn-compatible CSR matrix store under <work_path>/meth/matrix."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--work-path",
        help="Sample work directory containing allcools/*_merged_fr_bam_allcools/.",
    )
    source.add_argument(
        "--allc-dir",
        help="Directory of flat <barcode>_allc.gz files (test override).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/meth/matrix.",
    )
    parser.add_argument(
        "--cell-names",
        default=None,
        help="Optional explicit barcode list file (one barcode per line).",
    )
    parser.add_argument(
        "--barcode-mode",
        choices=["methylation_only", "gexcb"],
        default="methylation_only",
        help="Barcode list source when --cell-names is not set. Default: methylation_only.",
    )
    parser.add_argument(
        "--meth-context",
        default=DEFAULT_METH_CONTEXT,
        help="Context filter prefix (CG, CHG, CHH, CH, all). Default: CG.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=DEFAULT_CHUNKSIZE,
        help="Chromosome chunk size in bp for COO temp files. Default: 10000000.",
    )
    parser.add_argument(
        "--round-sites",
        action="store_true",
        help="Round ambiguous sites (0 < mc < cov) to majority vote.",
    )
    parser.add_argument(
        "--exclude-contigs",
        default="",
        help="Comma-separated contig names to skip.",
    )
    parser.add_argument(
        "--main-chroms-only",
        action="store_true",
        help="Keep only main chromosomes (chr1-19, chrX, chrY, chrM).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


def read_barcode_list(path: Path) -> list[str]:
    barcodes: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            barcode = raw_line.strip()
            if barcode and not barcode.startswith("#"):
                barcodes.append(barcode)
    return barcodes


def discover_work_allc_map(work_path: Path) -> dict[str, Path]:
    allc_map: dict[str, Path] = {}
    duplicates: dict[str, list[Path]] = defaultdict(list)
    for allc_path in sorted(work_path.glob("allcools/*_merged_fr_bam_allcools/*_allc.gz")):
        barcode = barcode_from_allc_path(allc_path)
        duplicates[barcode].append(allc_path)
        if barcode in allc_map:
            continue
        allc_map[barcode] = allc_path
    dup_only = {barcode: paths for barcode, paths in duplicates.items() if len(paths) > 1}
    if dup_only:
        details = "; ".join(
            f"{barcode} ({', '.join(str(path) for path in paths)})"
            for barcode, paths in sorted(dup_only.items())
        )
        raise ValueError(f"duplicate barcodes across ALLC chunks: {details}")
    return allc_map


def discover_allc_dir_map(allc_dir: Path) -> dict[str, Path]:
    allc_map: dict[str, Path] = {}
    for allc_path in sorted(allc_dir.glob("*_allc.gz")):
        barcode = barcode_from_allc_path(allc_path)
        if barcode in allc_map:
            raise ValueError(f"duplicate barcode in --allc-dir: {barcode}")
        allc_map[barcode] = allc_path
    return allc_map


def resolve_barcode_list_file(work_path: Path, barcode_mode: str) -> Path | None:
    if barcode_mode == "gexcb":
        merged_root = work_path / "split_bams" / "merged"
        paths = sorted(merged_root.glob("*_merge_filtered_barcode"))
        if not paths:
            return None
        return paths[0] if len(paths) == 1 else _union_barcode_files(paths)
    filtered = work_path / "cells" / "filtered_barcode"
    return filtered if filtered.is_file() else None


def _union_barcode_files(paths: list[Path]) -> Path:
    barcodes: list[str] = []
    seen: set[str] = set()
    for path in paths:
        for barcode in read_barcode_list(path):
            if barcode not in seen:
                seen.add(barcode)
                barcodes.append(barcode)
    union_path = paths[0].parent / "_allc_to_matrix_union_filtered_barcode"
    union_path.write_text("".join(f"{barcode}\n" for barcode in barcodes), encoding="utf-8")
    return union_path


def select_cells(
    allc_map: dict[str, Path],
    *,
    cell_names_path: Path | None,
    barcode_list_path: Path | None,
) -> tuple[list[str], list[Path], str]:
    if cell_names_path is not None:
        selected = read_barcode_list(cell_names_path)
        source = f"cell_names:{cell_names_path}"
    elif barcode_list_path is not None:
        selected = read_barcode_list(barcode_list_path)
        source = f"barcode_list:{barcode_list_path}"
    else:
        selected = sorted(allc_map.keys())
        source = "all_discovered_allc"

    missing = [barcode for barcode in selected if barcode not in allc_map]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "..." if len(missing) > 5 else ""
        raise FileNotFoundError(
            f"{len(missing)} selected barcodes have no ALLC file (e.g. {preview}{suffix})"
        )

    ordered = [barcode for barcode in selected if barcode in allc_map]
    paths = [allc_map[barcode] for barcode in ordered]
    return ordered, paths, source


def parse_exclude_contigs(value: str) -> set[str]:
    if not value.strip():
        return set()
    return {part.strip() for part in value.split(",") if part.strip()}


def main() -> int:
    args = parse_args()
    if args.chunksize <= 0:
        raise ValueError("--chunksize must be > 0")

    work_path = Path(args.work_path) if args.work_path else None
    allc_dir = Path(args.allc_dir) if args.allc_dir else None

    if work_path is not None:
        output_dir = (
            Path(args.output_dir) if args.output_dir else work_path / "meth" / "matrix"
        )
        allc_map = discover_work_allc_map(work_path)
        cell_names_path = Path(args.cell_names) if args.cell_names else None
        barcode_list_path = None if cell_names_path else resolve_barcode_list_file(
            work_path, args.barcode_mode
        )
        gather_mode = "work_path"
    else:
        assert allc_dir is not None
        if args.output_dir is None:
            raise ValueError("--output-dir is required when using --allc-dir")
        output_dir = Path(args.output_dir)
        allc_map = discover_allc_dir_map(allc_dir)
        cell_names_path = Path(args.cell_names) if args.cell_names else None
        barcode_list_path = None
        gather_mode = "allc_dir"

    if not allc_map:
        raise FileNotFoundError("no *_allc.gz files discovered")

    cell_names, allc_paths, barcode_source = select_cells(
        allc_map,
        cell_names_path=cell_names_path,
        barcode_list_path=barcode_list_path,
    )
    exclude_contigs = parse_exclude_contigs(args.exclude_contigs)

    print(f"[allc_to_matrix] gather_mode={gather_mode}")
    if work_path is not None:
        print(f"[allc_to_matrix] work_path={work_path}")
    if allc_dir is not None:
        print(f"[allc_to_matrix] allc_dir={allc_dir}")
    print(f"[allc_to_matrix] output_dir={output_dir}")
    print(f"[allc_to_matrix] discovered_allc_count={len(allc_map)}")
    print(f"[allc_to_matrix] selected_cell_count={len(cell_names)}")
    print(f"[allc_to_matrix] barcode_source={barcode_source}")
    print(f"[allc_to_matrix] meth_context={args.meth_context}")
    print(f"[allc_to_matrix] chunksize={args.chunksize}")
    print(f"[allc_to_matrix] round_sites={int(args.round_sites)}")
    print(f"[allc_to_matrix] main_chroms_only={int(args.main_chroms_only)}")
    print(f"[allc_to_matrix] exclude_contigs={sorted(exclude_contigs)}")

    if barcode_source == "all_discovered_allc":
        print("[allc_to_matrix] warning=no_barcode_list_using_all")

    if args.dry_run:
        for barcode, allc_path in list(zip(cell_names, allc_paths))[:5]:
            print(f"[allc_to_matrix] sample_cell={barcode} allc={allc_path}")
        print("[allc_to_matrix] dry_run=1")
        return 0

    outputs = build_matrix_store(
        allc_paths,
        cell_names,
        output_dir,
        meth_context=args.meth_context,
        chunksize=args.chunksize,
        round_sites=args.round_sites,
        exclude_contigs=exclude_contigs,
        main_chroms_only=args.main_chroms_only,
        run_info_extra={
            "gather_mode": gather_mode,
            "barcode_source": barcode_source,
            "barcode_mode": args.barcode_mode,
        },
    )
    print(f"[allc_to_matrix] matrix_dir={outputs['matrix_dir']}")
    print(f"[allc_to_matrix] column_header={outputs['column_header']}")
    print(f"[allc_to_matrix] cell_stats={outputs['cell_stats']}")
    print(f"[allc_to_matrix] run_info={outputs['run_info']}")
    print("[allc_to_matrix] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
