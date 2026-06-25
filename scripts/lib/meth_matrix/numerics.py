"""Numba kernels for shrunken residuals and region statistics."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(nogil=True)
def calc_mean_shrunken_residuals(
    data_chrom,
    indices_chrom,
    indptr_chrom,
    start,
    end,
    smoothed_vals,
    n_cells,
    chrom_len,
    shrinkage_factor=1,
):
    """Per-cell mean shrunken residuals for a genomic window."""
    shrunken_resid = np.full(n_cells, np.nan)
    if start > chrom_len:
        return shrunken_resid
    end += 1
    if end > chrom_len:
        end = chrom_len
    data = data_chrom[indptr_chrom[start] : indptr_chrom[end]]
    if data.size == 0:
        return shrunken_resid
    indices = indices_chrom[indptr_chrom[start] : indptr_chrom[end]]
    indptr = indptr_chrom[start : end + 1] - indptr_chrom[start]
    indptr_diff = np.diff(indptr)

    n_obs = np.zeros(n_cells, dtype=np.int64)
    n_obs_start = np.bincount(indices)
    n_obs[0 : n_obs_start.shape[0]] = n_obs_start

    meth_sums = np.zeros(n_cells, dtype=np.int64)
    smooth_sums = np.zeros(n_cells, dtype=np.float64)
    cpg_idx = 0
    nobs_cpg = indptr_diff[cpg_idx]
    for i in range(data.shape[0]):
        while nobs_cpg == 0:
            cpg_idx += 1
            nobs_cpg = indptr_diff[cpg_idx]
        nobs_cpg -= 1
        cell_idx = indices[i]
        smooth_sums[cell_idx] += smoothed_vals[start + cpg_idx]
        meth_value = data[i]
        if meth_value == -1:
            continue
        meth_sums[cell_idx] += 1

    for i in range(n_cells):
        if n_obs[i] > 0:
            shrunken_resid[i] = (meth_sums[i] - smooth_sums[i]) / (
                n_obs[i] + shrinkage_factor
            )
    return shrunken_resid


@njit
def count_n_cpg(region_indptr):
    """Count CpG sites in a region from a sliced CSR indptr."""
    prev_val = 0
    n_cpg = 0
    for val in region_indptr:
        if val != prev_val:
            n_cpg += 1
        prev_val = val
    return n_cpg


@njit
def count_n_cells(region_indices):
    """Count distinct cells with coverage in a region."""
    return np.unique(region_indices).size
