"""Per-cell mean methylation profiles around BED features."""

from __future__ import annotations

import csv
import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numba
import numpy as np
from numba import njit, prange

from .region_matrix import read_column_header, resolve_regions_label
from .smooth import load_chrom_csr

DEFAULT_WIDTH = 4000


def parse_profile_regions(
    bed_path: Path,
    strand_column: int | None = None,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Parse BED into per-chromosome start, end, and minus-strand arrays.

    ``strand_column`` is 1-indexed. When it is unset, every region is treated
    as ``+``. ``.`` is also treated as ``+``.
    """
    starts: dict[str, list[int]] = {}
    ends: dict[str, list[int]] = {}
    minus: dict[str, list[int]] = {}
    is_empty = True

    with bed_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            is_empty = False
            values = line.rstrip("\n").split("\t")
            if len(values) < 3:
                raise ValueError(
                    f"{bed_path}:{line_number}: need chrom, start, end"
                )
            chrom = values[0]
            start = int(values[1])
            end = int(values[2])
            if start < 0 or end < start:
                raise ValueError(
                    f"{bed_path}:{line_number}: invalid interval {start}-{end}"
                )
            strand_minus = 0
            if strand_column is not None:
                if len(values) < strand_column:
                    raise ValueError(
                        f"{bed_path}:{line_number}: missing strand column "
                        f"{strand_column}"
                    )
                strand = values[strand_column - 1]
                if strand == "-":
                    strand_minus = 1
                elif strand not in {"+", "."}:
                    raise ValueError(
                        f"{bed_path}:{line_number}: strand must be +, -, or . "
                        f"(got {strand!r})"
                    )
            starts.setdefault(chrom, []).append(start)
            ends.setdefault(chrom, []).append(end)
            minus.setdefault(chrom, []).append(strand_minus)

    if is_empty:
        raise ValueError(f"BED file is empty: {bed_path}")

    packed: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom in starts:
        packed[chrom] = (
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(minus[chrom], dtype=np.uint8),
        )
    return packed


@njit(parallel=True)
def _accumulate_profile(
    data,
    indices,
    indptr,
    starts,
    ends,
    strand_minus,
    chrom_len,
    half,
    width,
    partial_meth,
    partial_total,
):
    """Add per-cell counts into thread-local (cell, relative position) buffers.

    The window is ``width`` bp starting at ``center - half``, with
    ``center = (start + end) // 2`` and ``half = width // 2``. Relative
    position 0 is the center. Minus-strand regions negate that offset so
    upstream stays negative. The stored axis is ``[-half, width - half)``.
    A minus-strand site that lands on ``+half`` (even widths only) is outside
    that axis and is dropped.
    """
    n_regions = starts.shape[0]
    n_threads = partial_meth.shape[0]
    for thread_i in prange(n_threads):
        region_lo = n_regions * thread_i // n_threads
        region_hi = n_regions * (thread_i + 1) // n_threads
        for region_i in range(region_lo, region_hi):
            center = (starts[region_i] + ends[region_i]) // 2
            origin = center - half
            win_lo = origin
            if win_lo < 0:
                win_lo = 0
            win_hi = origin + width
            if win_hi > chrom_len + 1:
                win_hi = chrom_len + 1
            if win_lo >= win_hi:
                continue
            row_start = indptr[win_lo]
            row_end = indptr[win_hi]
            if row_start == row_end:
                continue
            region_data = data[row_start:row_end]
            region_indices = indices[row_start:row_end]
            local_indptr = indptr[win_lo : win_hi + 1] - row_start
            indptr_diff = np.diff(local_indptr)
            cpg_idx = 0
            nobs = indptr_diff[0]
            minus = strand_minus[region_i] != 0
            for i in range(region_data.shape[0]):
                while nobs == 0:
                    cpg_idx += 1
                    nobs = indptr_diff[cpg_idx]
                nobs -= 1
                rel = (win_lo + cpg_idx) - center
                if minus:
                    rel = -rel
                col = rel + half
                if col < 0 or col >= width:
                    continue
                cell_i = region_indices[i]
                partial_total[thread_i, cell_i, col] += 1
                if region_data[i] > 0:
                    partial_meth[thread_i, cell_i, col] += 1


def _accumulate_chromosome(
    mat,
    starts: np.ndarray,
    ends: np.ndarray,
    strand_minus: np.ndarray,
    *,
    width: int,
    n_threads: int,
) -> tuple[np.ndarray, np.ndarray]:
    half = width // 2
    chrom_len = mat.shape[0] - 1
    n_cells = mat.shape[1]
    partial_meth = np.zeros((n_threads, n_cells, width), dtype=np.int32)
    partial_total = np.zeros((n_threads, n_cells, width), dtype=np.int32)
    _accumulate_profile(
        mat.data,
        mat.indices,
        mat.indptr,
        starts,
        ends,
        strand_minus,
        chrom_len,
        half,
        width,
        partial_meth,
        partial_total,
    )
    n_meth = partial_meth.sum(axis=0, dtype=np.int64)
    n_total = partial_total.sum(axis=0, dtype=np.int64)
    return n_meth, n_total


def _write_profile_csv(
    path: Path,
    *,
    n_meth: np.ndarray,
    n_total: np.ndarray,
    cell_names: list[str],
    width: int,
    label: str | None,
) -> int:
    half = width // 2
    cell_idx, pos_idx = np.nonzero(n_total)
    order = np.lexsort((cell_idx, pos_idx))
    cell_idx = cell_idx[order]
    pos_idx = pos_idx[order]
    header = ["position", "cell_name", "meth_frac", "n_meth", "n_total"]
    if label:
        header.append("label")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for cell_i, pos_i in zip(cell_idx.tolist(), pos_idx.tolist()):
            meth = int(n_meth[cell_i, pos_i])
            total = int(n_total[cell_i, pos_i])
            row = [
                int(pos_i) - half,
                cell_names[cell_i],
                f"{meth / total:.6g}",
                meth,
                total,
            ]
            if label:
                row.append(label)
            writer.writerow(row)
    return int(cell_idx.size)


def _write_run_info(
    path: Path,
    *,
    begin_time: datetime,
    regions_bed: Path,
    regions_label: str,
    width: int,
    strand_column: int | None,
    label: str | None,
    n_cells: int,
    n_regions: int,
    n_regions_used: int,
    n_rows: int,
    threads: int,
    extra: dict | None = None,
) -> Path:
    end_time = datetime.now(timezone.utc)
    payload = {
        "stage": "meth_profile",
        "regions_bed": str(regions_bed),
        "regions_label": regions_label,
        "width": width,
        "strand_column": strand_column,
        "label": label or "",
        "n_cells": n_cells,
        "n_regions": n_regions,
        "n_regions_used": n_regions_used,
        "n_rows": n_rows,
        "threads": threads,
        "begin_time_utc": begin_time.isoformat(),
        "end_time_utc": end_time.isoformat(),
        "runtime_seconds": (end_time - begin_time).total_seconds(),
        **(extra or {}),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def build_profile(
    matrix_dir: Path,
    regions_bed: Path,
    output_dir: Path,
    *,
    width: int = DEFAULT_WIDTH,
    strand_column: int | None = None,
    label: str | None = None,
    threads: int = -1,
    regions_label: str | None = None,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Write a long per-cell methylation profile CSV for BED features."""
    matrix_dir = Path(matrix_dir)
    regions_bed = Path(regions_bed)
    output_dir = Path(output_dir)
    if width < 1:
        raise ValueError("width must be >= 1")
    if strand_column is not None and strand_column < 1:
        raise ValueError("strand_column must be >= 1")
    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")
    if not regions_bed.is_file():
        raise FileNotFoundError(f"regions BED not found: {regions_bed}")

    if threads > 0:
        numba.set_num_threads(threads)
    n_threads = numba.get_num_threads() if threads > 0 else (os.cpu_count() or 1)
    if n_threads < 1:
        n_threads = 1

    begin_time = datetime.now(timezone.utc)
    cell_names = read_column_header(matrix_dir)
    regions = parse_profile_regions(regions_bed, strand_column)
    n_regions = sum(starts.shape[0] for starts, _ends, _minus in regions.values())
    n_cells = len(cell_names)
    partial_bytes = n_threads * n_cells * width * 8
    print(
        f"[meth_profile] n_cells={n_cells} n_regions={n_regions} "
        f"width={width} threads={n_threads} partial_bytes={partial_bytes}"
    )

    n_meth = np.zeros((n_cells, width), dtype=np.int64)
    n_total = np.zeros((n_cells, width), dtype=np.int64)
    n_regions_used = 0
    for chrom in sorted(regions):
        starts, ends, strand_minus = regions[chrom]
        mat_path = matrix_dir / f"{chrom}.npz"
        if not mat_path.is_file():
            print(f"[meth_profile] warning=missing_chrom_npz chrom={chrom}")
            continue
        print(f"[meth_profile] chrom={chrom} regions={starts.size}")
        mat = load_chrom_csr(matrix_dir, chrom)
        if mat.shape[1] != n_cells:
            raise ValueError(
                f"{chrom}.npz has {mat.shape[1]} cells, column_header has {n_cells}"
            )
        meth_chrom, total_chrom = _accumulate_chromosome(
            mat,
            starts,
            ends,
            strand_minus,
            width=width,
            n_threads=n_threads,
        )
        n_meth += meth_chrom
        n_total += total_chrom
        n_regions_used += int(starts.size)

    if n_regions_used == 0:
        raise ValueError(
            "no BED regions matched a chromosome matrix under "
            f"{matrix_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / "profile.csv.gz"
    print(f"[meth_profile] writing {profile_path}")
    n_rows = _write_profile_csv(
        profile_path,
        n_meth=n_meth,
        n_total=n_total,
        cell_names=cell_names,
        width=width,
        label=label or None,
    )
    resolved_label = resolve_regions_label(regions_bed, regions_label)
    info_path = output_dir / "run_info.json"
    _write_run_info(
        info_path,
        begin_time=begin_time,
        regions_bed=regions_bed,
        regions_label=resolved_label,
        width=width,
        strand_column=strand_column,
        label=label,
        n_cells=n_cells,
        n_regions=n_regions,
        n_regions_used=n_regions_used,
        n_rows=n_rows,
        threads=n_threads,
        extra=run_info_extra,
    )
    return {
        "output_dir": output_dir,
        "profile_csv": profile_path,
        "run_info": info_path,
    }
