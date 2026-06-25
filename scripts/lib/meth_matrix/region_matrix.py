"""Per-region methylation matrices from CSR store and BED regions."""

from __future__ import annotations

import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numba
import numpy as np
import pandas as pd
from numba import njit, prange

from .smooth import list_chrom_npz_paths, load_chrom_csr, load_smoothed_chrom

SPARSE_REGION_CHUNK_SIZE = 10_000


@njit(parallel=True)
def _calc_mean_mfracs(
    data_chrom,
    indices_chrom,
    indptr_chrom,
    starts,
    ends,
    chrom_len,
    n_cells,
    smoothed_vals,
    chunk_size=500,
):
    """Per-region methylated counts, totals, and mean shrunken residuals."""
    n_regions = starts.shape[0]
    ends = ends + 1

    n_meth = np.zeros((n_cells, n_regions), dtype=np.int64)
    n_total = np.zeros((n_cells, n_regions), dtype=np.int64)
    smooth_sums = np.full((n_cells, n_regions), np.nan, dtype=np.float32)

    chunks = np.arange(0, n_regions, chunk_size)
    for chunk_i in prange(chunks.shape[0]):
        chunk_start = chunks[chunk_i]
        chunk_end = chunk_start + chunk_size
        if chunk_end > n_regions:
            chunk_end = n_regions

        for region_i in range(chunk_start, chunk_end):
            start = starts[region_i]
            if start > chrom_len:
                continue
            end = ends[region_i]
            if end > chrom_len:
                end = chrom_len
            data = data_chrom[indptr_chrom[start] : indptr_chrom[end]]
            if data.size == 0:
                continue
            indices = indices_chrom[indptr_chrom[start] : indptr_chrom[end]]
            indptr = indptr_chrom[start : end + 1] - indptr_chrom[start]
            indptr_diff = np.diff(indptr)
            cpg_idx = 0
            nobs_cpg = indptr_diff[cpg_idx]
            for i in range(data.shape[0]):
                while nobs_cpg == 0:
                    cpg_idx += 1
                    nobs_cpg = indptr_diff[cpg_idx]
                nobs_cpg -= 1
                cell_i = indices[i]
                meth_value = data[i]
                n_total[cell_i, region_i] += 1
                if np.isnan(smooth_sums[cell_i, region_i]):
                    smooth_sums[cell_i, region_i] = 0.0
                if meth_value == -1:
                    smooth_sums[cell_i, region_i] -= smoothed_vals[start + cpg_idx]
                    continue
                smooth_sums[cell_i, region_i] += 1 - smoothed_vals[start + cpg_idx]
                n_meth[cell_i, region_i] += meth_value
    mean_shrunk_res = smooth_sums / (n_total + 1)
    return n_meth, n_total, mean_shrunk_res


@njit
def _dense_to_sparse(n_meth, n_total, mean_shrunk_res):
    row_i, col_i = n_total.nonzero()
    residuals = np.empty(row_i.size, dtype=np.float32)
    mfracs = np.empty(row_i.size, dtype=np.float32)
    coverage = np.empty(row_i.size, dtype=np.int64)
    for i, (r, c) in enumerate(zip(row_i, col_i)):
        residuals[i] = mean_shrunk_res[r, c]
        mfracs[i] = n_meth[r, c] / n_total[r, c]
        coverage[i] = n_total[r, c]
    return row_i, col_i, residuals, mfracs, coverage


def read_column_header(matrix_dir: Path) -> list[str]:
    header_path = matrix_dir / "column_header.txt"
    if not header_path.is_file():
        raise FileNotFoundError(f"column_header.txt not found: {header_path}")
    cell_names: list[str] = []
    with header_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            name = line.strip()
            if name:
                cell_names.append(name)
    if not cell_names:
        raise ValueError(f"no cell names in {header_path}")
    return cell_names


def parse_bed_regions(bed_path: Path) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    """Parse BED; return per-chromosome start/end coordinate lists."""
    start_dict: dict[str, list[int]] = {}
    end_dict: dict[str, list[int]] = {}
    is_empty = True

    with bed_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            is_empty = False
            values = line.strip().split("\t")
            if len(values) < 3:
                raise ValueError(f"invalid BED line (need chrom,start,end): {line.strip()}")
            chrom = values[0]
            start = int(values[1])
            end = int(values[2])
            start_dict.setdefault(chrom, []).append(start)
            end_dict.setdefault(chrom, []).append(end)

    if is_empty:
        raise ValueError(f"BED file is empty: {bed_path}")
    return start_dict, end_dict


def resolve_regions_label(bed_path: Path, override: str | None = None) -> str:
    if override and override.strip():
        return override.strip()
    return bed_path.stem or "regions"


def _write_dense_mtx(mtx_list, row_names, col_names, out_dir: Path, fname: str) -> Path:
    out_path = out_dir / fname
    print(f"[meth_matrix] writing {out_path}")
    frame = pd.DataFrame(
        data=np.hstack(mtx_list),
        index=row_names,
        columns=col_names,
    )
    frame.to_csv(out_path)
    return out_path


def _write_sparse_mtx_chunk(
    out_path: Path,
    row_i,
    col_i,
    residuals,
    mfracs,
    coverage,
    *,
    region_n_offset: int = 0,
) -> None:
    frame = pd.DataFrame(
        {
            "row_i": row_i + 1,
            "col_i": col_i + 1 + region_n_offset,
            "residuals": residuals,
            "mfracs": mfracs,
            "coverage": coverage,
        }
    )
    frame.to_csv(
        out_path,
        mode="a",
        header=False,
        index=False,
        compression="gzip",
        sep=" ",
        float_format="%.4g",
    )


def _finalize_sparse_mtx(
    out_dir: Path, region_names: list[str], cell_names: list[str]
) -> dict[str, Path]:
    features_path = out_dir / "features.tsv.gz"
    barcodes_path = out_dir / "barcodes.tsv.gz"
    with gzip.open(features_path, "wt", encoding="utf-8") as region_out:
        region_out.write("\n".join(region_names))
    with gzip.open(barcodes_path, "wt", encoding="utf-8") as bc_out:
        bc_out.write("\n".join(cell_names))
    return {"features": features_path, "barcodes": barcodes_path}


def _prepare_sparse_out_dir(out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    mtx_path = out_dir / "matrix.mtx.gz"
    for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
        old = out_dir / name
        if old.is_file():
            print(f"[meth_matrix] removing previous {old}")
            old.unlink()
    return mtx_path


def _write_run_info(
    path: Path,
    *,
    begin_time: datetime,
    output_format: str,
    regions_bed: Path,
    regions_label: str,
    n_cells: int,
    n_regions: int,
    dense: bool,
    threads: int,
    extra: dict | None = None,
) -> Path:
    end_time = datetime.now(timezone.utc)
    payload = {
        "stage": "meth_matrix",
        "output_format": output_format,
        "dense": dense,
        "regions_bed": str(regions_bed),
        "regions_label": regions_label,
        "n_cells": n_cells,
        "n_regions": n_regions,
        "threads": threads,
        "begin_time_utc": begin_time.isoformat(),
        "end_time_utc": end_time.isoformat(),
        "runtime_seconds": (end_time - begin_time).total_seconds(),
        **(extra or {}),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def build_region_matrices(
    matrix_dir: Path,
    regions_bed: Path,
    output_dir: Path,
    *,
    dense: bool = False,
    threads: int = -1,
    regions_label: str | None = None,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Build per-region methylation matrices (sparse default, dense optional)."""
    matrix_dir = Path(matrix_dir)
    regions_bed = Path(regions_bed)
    output_dir = Path(output_dir)
    smoothed_dir = matrix_dir / "smoothed"

    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")
    if not smoothed_dir.is_dir():
        raise FileNotFoundError(f"smoothed directory not found: {smoothed_dir}")
    if not regions_bed.is_file():
        raise FileNotFoundError(f"regions BED not found: {regions_bed}")

    npz_paths = list_chrom_npz_paths(matrix_dir)
    if not npz_paths:
        raise FileNotFoundError(f"no *.npz matrices found under {matrix_dir}")

    if threads != -1:
        numba.set_num_threads(threads)
    n_threads = numba.get_num_threads() if threads != -1 else (os.cpu_count() or 1)

    begin_time = datetime.now(timezone.utc)
    cell_names = read_column_header(matrix_dir)
    start_dict, end_dict = parse_bed_regions(regions_bed)
    region_names: list[str] = []
    label = resolve_regions_label(regions_bed, regions_label)

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {"output_dir": output_dir}

    if dense:
        meth_parts: list[np.ndarray] = []
        total_parts: list[np.ndarray] = []
        msr_parts: list[np.ndarray] = []

        for chrom in sorted(start_dict.keys()):
            mat_path = matrix_dir / f"{chrom}.npz"
            if not mat_path.is_file():
                print(f"[meth_matrix] warning=missing_chrom_npz chrom={chrom}")
                continue
            starts = np.asarray(start_dict[chrom], dtype=np.int64)
            ends = np.asarray(end_dict[chrom], dtype=np.int64)
            print(f"[meth_matrix] chrom={chrom} regions={starts.size}")
            mat = load_chrom_csr(matrix_dir, chrom)
            chrom_len, n_cells = mat.shape[0] - 1, mat.shape[1]
            smoothed_vals = load_smoothed_chrom(smoothed_dir, chrom)
            meth_chrom, total_chrom, msr_chrom = _calc_mean_mfracs(
                mat.data,
                mat.indices,
                mat.indptr,
                starts,
                ends,
                chrom_len,
                n_cells,
                smoothed_vals,
                chunk_size=500,
            )
            meth_parts.append(meth_chrom)
            total_parts.append(total_chrom)
            msr_parts.append(msr_chrom)
            for start, end in zip(starts, ends):
                region_names.append(f"{chrom}:{start}-{end}")

        print(f"[meth_matrix] writing dense tables to {output_dir}")
        mfrac_parts = [np.divide(m, t, out=np.zeros_like(m, dtype=float), where=t > 0)
                       for m, t in zip(meth_parts, total_parts)]
        outputs["methylation_fractions"] = _write_dense_mtx(
            mfrac_parts, cell_names, region_names, output_dir, "methylation_fractions.csv.gz"
        )
        outputs["mean_shrunken_residuals"] = _write_dense_mtx(
            msr_parts, cell_names, region_names, output_dir, "mean_shrunken_residuals.csv.gz"
        )
        outputs["total_sites"] = _write_dense_mtx(
            total_parts, cell_names, region_names, output_dir, "total_sites.csv.gz"
        )
        outputs["methylated_sites"] = _write_dense_mtx(
            meth_parts, cell_names, region_names, output_dir, "methylated_sites.csv.gz"
        )
        output_format = "dense"
    else:
        mtx_path = _prepare_sparse_out_dir(output_dir)
        outputs["matrix_mtx"] = mtx_path
        n_processed_regions = 0

        for chrom in sorted(start_dict.keys()):
            mat_path = matrix_dir / f"{chrom}.npz"
            if not mat_path.is_file():
                print(f"[meth_matrix] warning=missing_chrom_npz chrom={chrom}")
                continue

            starts = np.asarray(start_dict[chrom], dtype=np.int64)
            ends = np.asarray(end_dict[chrom], dtype=np.int64)
            n_regions = starts.size
            print(f"[meth_matrix] chrom={chrom} regions={n_regions}")
            mat = load_chrom_csr(matrix_dir, chrom)
            chrom_len, n_cells = mat.shape[0] - 1, mat.shape[1]
            smoothed_vals = load_smoothed_chrom(smoothed_dir, chrom)

            chunk_size = SPARSE_REGION_CHUNK_SIZE
            chunks = np.arange(0, n_regions, chunk_size)
            for chunk_i in range(chunks.shape[0]):
                chunk_start = int(chunks[chunk_i])
                chunk_end = min(chunk_start + chunk_size, n_regions)
                n_meth, n_total, mean_shrunk_res = _calc_mean_mfracs(
                    mat.data,
                    mat.indices,
                    mat.indptr,
                    starts[chunk_start:chunk_end],
                    ends[chunk_start:chunk_end],
                    chrom_len,
                    n_cells,
                    smoothed_vals,
                    chunk_size=(chunk_size // n_threads) + 1,
                )
                row_i, col_i, residuals, mfracs, coverage = _dense_to_sparse(
                    n_meth, n_total, mean_shrunk_res
                )
                _write_sparse_mtx_chunk(
                    mtx_path,
                    row_i,
                    col_i,
                    residuals,
                    mfracs,
                    coverage,
                    region_n_offset=n_processed_regions,
                )
                n_processed_regions += chunk_end - chunk_start

            for start, end in zip(starts, ends):
                region_names.append(f"{chrom}:{start}-{end}")

        sidecars = _finalize_sparse_mtx(output_dir, region_names, cell_names)
        outputs.update(sidecars)
        output_format = "sparse"

    if not region_names:
        raise ValueError(
            "no regions produced methylation matrix output "
            "(missing chromosome matrices for all BED regions?)"
        )

    info_path = output_dir / "run_info.json"
    outputs["run_info"] = _write_run_info(
        info_path,
        begin_time=begin_time,
        output_format=output_format,
        regions_bed=regions_bed,
        regions_label=label,
        n_cells=len(cell_names),
        n_regions=len(region_names),
        dense=dense,
        threads=n_threads,
        extra=run_info_extra,
    )
    return outputs
