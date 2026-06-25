"""Tricube pseudobulk smoothing over CSR methylation matrices."""

from __future__ import annotations

import gzip
import json
from datetime import datetime, timezone
from glob import glob
from pathlib import Path

import numpy as np
import scipy.sparse as sp_sparse

np.seterr(divide="ignore", invalid="ignore")


def list_chrom_npz_paths(matrix_dir: Path) -> list[Path]:
    """Return chromosome CSR paths sorted by name."""
    return sorted(
        path
        for path in matrix_dir.glob("*.npz")
        if path.is_file()
    )


def load_chrom_csr(matrix_dir: Path, chrom: str) -> sp_sparse.csr_matrix:
    """Load one chromosome CSR matrix."""
    mat_path = matrix_dir / f"{chrom}.npz"
    if not mat_path.is_file():
        raise FileNotFoundError(f"CSR matrix not found: {mat_path}")
    return sp_sparse.load_npz(mat_path)


def smooth_chromosome(
    mat: sp_sparse.csr_matrix,
    *,
    bandwidth: int,
    use_weights: bool = False,
) -> dict[int, float]:
    """Smooth per-position methylation fractions with a tricube kernel."""
    if bandwidth < 1:
        raise ValueError("bandwidth must be >= 1")
    hbw = bandwidth // 2
    rel_dist = np.abs((np.arange(bandwidth) - hbw) / hbw)
    kernel = (1 - (rel_dist**3)) ** 3

    n_obs = mat.getnnz(axis=1)
    n_meth = np.ravel(np.sum(mat > 0, axis=1))
    mfracs = np.divide(n_meth, n_obs, out=np.full(n_meth.shape, np.nan), where=n_obs > 0)
    cpg_pos = (~np.isnan(mfracs)).nonzero()[0]
    weights = np.log1p(n_obs) if use_weights else None

    smoothed: dict[int, float] = {}
    for pos in cpg_pos:
        window = mfracs[pos - hbw : pos + hbw]
        nz = ~np.isnan(window)
        try:
            k = kernel[nz]
            if use_weights:
                w = weights[pos - hbw : pos + hbw][nz]
                smooth_val = np.divide(np.sum(window[nz] * k * w), np.sum(k * w))
            else:
                smooth_val = np.divide(np.sum(window[nz] * k), np.sum(k))
            smoothed[int(pos)] = float(smooth_val)
        except IndexError:
            smoothed[int(pos)] = float("nan")
    return smoothed


def write_smoothed_chrom(
    smoothed: dict[int, float],
    output_path: Path,
) -> None:
    """Write position,smoothed_mfrac rows as gzip CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        for pos in sorted(smoothed):
            handle.write(f"{pos},{smoothed[pos]}\n")


def load_smoothed_chrom(smoothed_dir: Path, chrom: str):
    """Load smoothed values for one chromosome into a numba typed dict."""
    import numba

    csv_gz = smoothed_dir / f"{chrom}.csv.gz"
    csv_plain = smoothed_dir / f"{chrom}.csv"
    if csv_gz.is_file():
        import pandas as pd

        frame = pd.read_csv(csv_gz, delimiter=",", header=None, dtype="float64")
    elif csv_plain.is_file():
        import pandas as pd

        frame = pd.read_csv(csv_plain, delimiter=",", header=None, dtype="float64")
    else:
        raise FileNotFoundError(
            f"smoothed methylation not found for {chrom} under {smoothed_dir}"
        )
    typed_dict = numba.typed.Dict.empty(
        key_type=numba.types.int64,
        value_type=numba.types.float64,
    )
    values = frame.values
    for row_idx in range(values.shape[0]):
        typed_dict[int(values[row_idx, 0])] = values[row_idx, 1]
    return typed_dict


def smooth_matrix_store(
    matrix_dir: Path,
    *,
    bandwidth: int = 1000,
    use_weights: bool = False,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Smooth all chromosome CSR matrices; write smoothed/*.csv.gz."""
    matrix_dir = Path(matrix_dir)
    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"matrix directory not found: {matrix_dir}")

    mat_paths = list_chrom_npz_paths(matrix_dir)
    if not mat_paths:
        raise FileNotFoundError(f"no *.npz matrices found under {matrix_dir}")

    begin_time = datetime.now(timezone.utc)
    smoothed_dir = matrix_dir / "smoothed"
    smoothed_dir.mkdir(parents=True, exist_ok=True)
    chromosomes: list[str] = []

    for mat_path in mat_paths:
        chrom = mat_path.stem
        chromosomes.append(chrom)
        print(f"[meth_smooth] chrom={chrom} reading={mat_path}")
        mat = sp_sparse.load_npz(mat_path)
        print(f"[meth_smooth] chrom={chrom} smoothing rows={mat.shape[0]} cols={mat.shape[1]}")
        smoothed = smooth_chromosome(mat, bandwidth=bandwidth, use_weights=use_weights)
        out_path = smoothed_dir / f"{chrom}.csv.gz"
        write_smoothed_chrom(smoothed, out_path)
        print(f"[meth_smooth] chrom={chrom} sites={len(smoothed)} output={out_path}")

    info_path = matrix_dir / "smoothed" / "run_info.json"
    end_time = datetime.now(timezone.utc)
    payload = {
        "stage": "meth_smooth",
        "matrix_dir": str(matrix_dir),
        "smoothed_dir": str(smoothed_dir),
        "bandwidth": bandwidth,
        "use_weights": use_weights,
        "chromosomes": chromosomes,
        "begin_time_utc": begin_time.isoformat(),
        "end_time_utc": end_time.isoformat(),
        "runtime_seconds": (end_time - begin_time).total_seconds(),
        **(run_info_extra or {}),
    }
    info_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "matrix_dir": matrix_dir,
        "smoothed_dir": smoothed_dir,
        "run_info": info_path,
    }
