"""COO chunking and CSR matrix store for per-cell methylation data."""

from __future__ import annotations

import json
from collections.abc import Sequence
from datetime import datetime, timezone
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from numba import njit

from .allc import iter_allc_sites


@njit
def _process_chunk(positions, indptr, last_pos, indptr_counter, indptr_i):
    for pos in positions:
        if pos > last_pos:
            for _ in range(pos - last_pos):
                indptr[indptr_i] = indptr_counter
                indptr_i += 1
            last_pos = pos
        indptr_counter += 1
    return last_pos, indptr_counter, indptr_i


def build_matrix_store(
    allc_paths: Sequence[Path],
    cell_names: Sequence[str],
    output_dir: Path,
    *,
    meth_context: str = "CG",
    chunksize: int = 10_000_000,
    round_sites: bool = False,
    exclude_contigs: set[str] | None = None,
    main_chroms_only: bool = False,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Build CSR store from per-cell ALLC files; return output paths."""
    if len(allc_paths) != len(cell_names):
        raise ValueError("allc_paths and cell_names length mismatch")
    if not allc_paths:
        raise ValueError("no ALLC inputs provided")

    output_dir.mkdir(parents=True, exist_ok=True)
    begin_time = datetime.now(timezone.utc)
    n_cells = len(cell_names)
    cell_index = {name: index for index, name in enumerate(cell_names)}

    coo_handles: dict[tuple[str, int], object] = {}
    chrom_sizes: dict[str, int] = {}

    for cell_n, allc_path in enumerate(allc_paths):
        if cell_n % 50 == 0:
            print(
                f"[allc_to_matrix] progress={100 * cell_n / n_cells:.2f}% "
                f"cell_index={cell_n}/{n_cells}"
            )
        for chrom, genomic_pos, meth_value in iter_allc_sites(
            allc_path,
            meth_context=meth_context,
            round_sites=round_sites,
            exclude_contigs=exclude_contigs,
            main_chroms_only=main_chroms_only,
        ):
            chrom_chunk = int(genomic_pos // chunksize)
            coo_key = (chrom, chrom_chunk)
            if coo_key not in coo_handles:
                coo_path = output_dir / f"{chrom}_chunk{chrom_chunk:07}.coo"
                coo_handles[coo_key] = coo_path.open("w", encoding="utf-8")
                chrom_sizes.setdefault(chrom, 0)
            if genomic_pos > chrom_sizes[chrom]:
                chrom_sizes[chrom] = genomic_pos
            coo_handles[coo_key].write(f"{genomic_pos},{cell_n},{meth_value}\n")

    for handle in coo_handles.values():
        handle.close()

    print("[allc_to_matrix] progress=100.00% coo_dump_done=1")

    n_obs_cell = np.zeros(n_cells, dtype=np.int64)
    n_meth_cell = np.zeros(n_cells, dtype=np.int64)

    for chrom, chrom_size in chrom_sizes.items():
        print(
            f"[allc_to_matrix] csr_chrom={chrom} rows={chrom_size + 1} "
            f"cols={n_cells}"
        )
        mat = _load_csr_from_coo(output_dir, chrom, chrom_size, n_cells)
        n_obs_cell += np.asarray(mat.getnnz(axis=0)).ravel()
        n_meth_cell += np.ravel(np.sum(mat > 0, axis=0))
        mat_path = output_dir / f"{chrom}.npz"
        sp_sparse.save_npz(mat_path, mat)
        _delete_coo_chunks(output_dir, chrom)

    colname_path = _write_column_names(output_dir, cell_names)
    stats_path = _write_summary_stats(output_dir, cell_names, n_obs_cell, n_meth_cell)
    info_path = _write_run_info(
        output_dir / "run_info.json",
        begin_time=begin_time,
        meth_context=meth_context,
        chunksize=chunksize,
        round_sites=round_sites,
        exclude_contigs=sorted(exclude_contigs or []),
        main_chroms_only=main_chroms_only,
        cell_names=list(cell_names),
        allc_paths=[str(path) for path in allc_paths],
        chromosomes=sorted(chrom_sizes.keys()),
        extra=run_info_extra or {},
    )

    return {
        "matrix_dir": output_dir,
        "column_header": colname_path,
        "cell_stats": stats_path,
        "run_info": info_path,
    }


def _iter_chunks(data_dir: Path, chrom: str):
    chunk_paths = sorted(glob(str(data_dir / f"{chrom}_chunk*.coo")))
    for chunk_path in chunk_paths:
        print(f"[allc_to_matrix] coo_chunk={Path(chunk_path).name}")
        chunk = pd.read_csv(chunk_path, delimiter=",", header=None).values
        yield chunk


def _delete_coo_chunks(data_dir: Path, chrom: str) -> None:
    for chunk_path in glob(str(data_dir / f"{chrom}_chunk*.coo")):
        Path(chunk_path).unlink(missing_ok=True)


def _load_csr_from_coo(
    data_dir: Path, chrom: str, chrom_size: int, n_cells: int
) -> sp_sparse.csr_matrix:
    data_chunks: list[np.ndarray] = []
    indices_chunks: list[np.ndarray] = []
    indptr = np.empty(chrom_size + 2, dtype=np.int64)

    last_pos = -1
    indptr_counter = 0
    indptr_i = 0
    for chunk in _iter_chunks(data_dir, chrom):
        sorting_idx = np.lexsort((chunk[:, 1], chunk[:, 0]))
        last_pos, indptr_counter, indptr_i = _process_chunk(
            chunk[sorting_idx, 0], indptr, last_pos, indptr_counter, indptr_i
        )
        data_chunks.append(chunk[sorting_idx, 2].astype(np.int8))
        indices_chunks.append(chunk[sorting_idx, 1].astype(np.int64))
    indptr[indptr_i] = indptr_counter
    data = np.concatenate(data_chunks) if data_chunks else np.array([], dtype=np.int8)
    indices = (
        np.concatenate(indices_chunks) if indices_chunks else np.array([], dtype=np.int64)
    )
    return sp_sparse.csr_matrix(
        (data, indices, indptr), shape=(chrom_size + 1, n_cells)
    )


def _write_column_names(output_dir: Path, cell_names: Sequence[str]) -> Path:
    out_path = output_dir / "column_header.txt"
    out_path.write_text("".join(f"{name}\n" for name in cell_names), encoding="utf-8")
    return out_path


def _write_summary_stats(
    output_dir: Path,
    cell_names: Sequence[str],
    n_obs: np.ndarray,
    n_meth: np.ndarray,
) -> Path:
    stats_df = pd.DataFrame(
        {
            "cell_name": list(cell_names),
            "n_obs": n_obs,
            "n_meth": n_meth,
            "global_meth_frac": np.divide(
                n_meth,
                n_obs,
                out=np.zeros_like(n_meth, dtype=float),
                where=n_obs > 0,
            ),
        }
    )
    out_path = output_dir / "cell_stats.csv"
    out_path.write_text(stats_df.to_csv(index=False), encoding="utf-8")
    return out_path


def _write_run_info(path: Path, *, begin_time: datetime, extra: dict, **kwargs) -> Path:
    end_time = datetime.now(timezone.utc)
    payload = {
        "stage": "allc_to_matrix",
        "begin_time_utc": begin_time.isoformat(),
        "end_time_utc": end_time.isoformat(),
        "runtime_seconds": (end_time - begin_time).total_seconds(),
        **kwargs,
        **extra,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path
