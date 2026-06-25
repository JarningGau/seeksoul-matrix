"""Sliding-window VMR scan over smoothed methylation matrices."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numba
import numpy as np
from numba import njit, prange

from .numerics import calc_mean_shrunken_residuals, count_n_cells, count_n_cpg
from .smooth import list_chrom_npz_paths, load_chrom_csr, load_smoothed_chrom

np.seterr(divide="ignore", invalid="ignore")


@njit
def find_peaks(smoothed_vars, swindow_centers, var_cutoff, half_bw, bridge_gaps=0):
    """Merge overlapping high-variance windows into VMR intervals."""
    peak_starts = []
    peak_ends = []
    prev_pos = 0
    in_peak = False
    for var, pos in zip(smoothed_vars, swindow_centers):
        if var > var_cutoff:
            if not in_peak:
                in_peak = True
                if peak_ends and pos - half_bw - bridge_gaps <= max(peak_ends):
                    peak_ends.pop()
                else:
                    peak_starts.append(pos - half_bw)
        else:
            if in_peak:
                in_peak = False
                peak_ends.append(prev_pos + half_bw)
        prev_pos = pos
    if in_peak:
        peak_ends.append(pos + half_bw)
    return peak_starts, peak_ends


@njit(parallel=True)
def move_windows(
    start,
    end,
    stepsize,
    half_bw,
    data_chrom,
    indices_chrom,
    indptr_chrom,
    smoothed_vals,
    n_cells,
    chrom_len,
):
    """Slide windows and compute variance of shrunken residuals per window."""
    windows = np.arange(start, end, stepsize)
    smoothed_var = np.empty(windows.shape, dtype=np.float64)
    for i in prange(windows.shape[0]):
        pos = windows[i]
        mean_shrunk_resid = calc_mean_shrunken_residuals(
            data_chrom,
            indices_chrom,
            indptr_chrom,
            pos - half_bw,
            pos + half_bw,
            smoothed_vals,
            n_cells,
            chrom_len,
        )
        smoothed_var[i] = np.nanvar(mean_shrunk_resid)
    return windows, smoothed_var


def scan_vmrs(
    matrix_dir: Path,
    output_bed: Path,
    *,
    bandwidth: int = 2000,
    stepsize: int = 100,
    var_threshold: float = 0.02,
    min_cells: int = 6,
    bridge_gaps: int = 0,
    threads: int = -1,
    write_header: bool = False,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Scan genome for variably methylated regions; write BED output."""
    matrix_dir = Path(matrix_dir)
    output_bed = Path(output_bed)
    smoothed_dir = matrix_dir / "smoothed"
    if not smoothed_dir.is_dir():
        raise FileNotFoundError(f"smoothed directory not found: {smoothed_dir}")

    mat_paths = sorted(
        list_chrom_npz_paths(matrix_dir),
        key=lambda path: path.stat().st_size,
        reverse=True,
    )
    if not mat_paths:
        raise FileNotFoundError(f"no *.npz matrices found under {matrix_dir}")

    if threads != -1:
        numba.set_num_threads(threads)
    n_threads = numba.get_num_threads()
    half_bw = bandwidth // 2

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    begin_time = datetime.now(timezone.utc)
    var_threshold_value = None
    n_vmrs_total = 0
    n_vmrs_hq = 0
    chrom_stats: list[dict] = []

    with output_bed.open("w", encoding="utf-8") as bed_out:
        if write_header:
            bed_out.write(
                "chromosome\tVMR_start\tVMR_end\tvariance\tn_sites\tn_cells\n"
            )

        for mat_path in mat_paths:
            chrom = mat_path.stem
            mat = load_chrom_csr(matrix_dir, chrom)
            smoothed_cpg_vals = load_smoothed_chrom(smoothed_dir, chrom)
            chrom_len, n_cells = mat.shape
            cpg_pos_chrom = np.nonzero(mat.getnnz(axis=1))[0]
            if cpg_pos_chrom.size == 0:
                print(f"[meth_scan] chrom={chrom} skipped=no_cpg_sites")
                continue

            if n_threads > 1:
                print(f"[meth_scan] chrom={chrom} threads={n_threads}")
            else:
                print(f"[meth_scan] chrom={chrom}")

            start = int(cpg_pos_chrom[0] + half_bw + 1)
            end = int(cpg_pos_chrom[-1] - half_bw - 1)
            if start >= end:
                print(f"[meth_scan] chrom={chrom} skipped=window_out_of_range")
                continue

            genomic_positions, window_variances = move_windows(
                start,
                end,
                stepsize,
                half_bw,
                mat.data,
                mat.indices,
                mat.indptr,
                smoothed_cpg_vals,
                n_cells,
                chrom_len,
            )

            if var_threshold_value is None:
                var_threshold_value = float(
                    np.nanquantile(window_variances, 1 - var_threshold)
                )
                print(f"[meth_scan] var_threshold_value={var_threshold_value}")

            peak_starts, peak_ends = find_peaks(
                window_variances,
                genomic_positions,
                var_threshold_value,
                half_bw,
                bridge_gaps,
            )

            n_vmrs_chrom_hq = 0
            for ps, pe in zip(peak_starts, peak_ends):
                n_vmrs_total += 1
                region_indices = mat.indices[mat.indptr[ps] : mat.indptr[pe + 1]]
                n_obs_cells = int(count_n_cells(region_indices))
                if n_obs_cells < min_cells:
                    continue
                region_indptr = mat.indptr[ps : pe + 2] - mat.indptr[ps]
                n_cpg = int(count_n_cpg(region_indptr))
                peak_var = float(
                    np.nanvar(
                        calc_mean_shrunken_residuals(
                            mat.data,
                            mat.indices,
                            mat.indptr,
                            ps,
                            pe,
                            smoothed_cpg_vals,
                            n_cells,
                            chrom_len,
                        )
                    )
                )
                bed_out.write(
                    f"{chrom}\t{ps}\t{pe}\t{peak_var}\t{n_cpg}\t{n_obs_cells}\n"
                )
                n_vmrs_hq += 1
                n_vmrs_chrom_hq += 1

            chrom_stats.append(
                {
                    "chrom": chrom,
                    "vmrs_reported": n_vmrs_chrom_hq,
                }
            )
            print(f"[meth_scan] chrom={chrom} vmrs_reported={n_vmrs_chrom_hq}")

    if not n_vmrs_hq:
        raise RuntimeError(
            f"found {n_vmrs_total} potential VMRs but none with coverage in "
            f"at least {min_cells} cells; try lowering --min-cells"
        )

    info_path = output_bed.parent / "run_info.json"
    end_time = datetime.now(timezone.utc)
    payload = {
        "stage": "meth_scan",
        "matrix_dir": str(matrix_dir),
        "output_bed": str(output_bed),
        "bandwidth": bandwidth,
        "stepsize": stepsize,
        "var_threshold": var_threshold,
        "var_threshold_value": var_threshold_value,
        "min_cells": min_cells,
        "bridge_gaps": bridge_gaps,
        "threads": n_threads,
        "n_vmrs_total": n_vmrs_total,
        "n_vmrs_reported": n_vmrs_hq,
        "chromosomes": chrom_stats,
        "begin_time_utc": begin_time.isoformat(),
        "end_time_utc": end_time.isoformat(),
        "runtime_seconds": (end_time - begin_time).total_seconds(),
        **(run_info_extra or {}),
    }
    info_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"[meth_scan] vmrs_reported={n_vmrs_hq} output={output_bed}")
    return {
        "output_bed": output_bed,
        "run_info": info_path,
    }
