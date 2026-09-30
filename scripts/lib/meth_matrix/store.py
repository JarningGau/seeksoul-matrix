"""Per-cell ALLC assembly into a CSR matrix store."""

from __future__ import annotations

import json
import multiprocessing as mp
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse

from .allc import read_cell_sites


def _read_cell_worker(
    task: tuple[str, str, bool, tuple[str, ...], bool],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    allc_path, meth_context, round_sites, exclude_contigs, main_chroms_only = task
    return read_cell_sites(
        Path(allc_path),
        meth_context=meth_context,
        round_sites=round_sites,
        exclude_contigs=set(exclude_contigs),
        main_chroms_only=main_chroms_only,
    )


def _iter_cell_sites(
    tasks: Sequence[tuple[str, str, bool, tuple[str, ...], bool]],
    *,
    threads: int,
):
    n_cells = len(tasks)
    if threads == 1 or n_cells <= 1:
        for cell_n, task in enumerate(tasks):
            if cell_n % 50 == 0:
                print(
                    f"[allc_to_matrix] progress={100 * cell_n / n_cells:.2f}% "
                    f"cell_index={cell_n}/{n_cells}"
                )
            yield _read_cell_worker(task)
        return

    ctx = mp.get_context("fork")
    with ctx.Pool(processes=min(threads, n_cells)) as pool:
        for cell_n, cell_sites in enumerate(
            pool.imap(_read_cell_worker, tasks, chunksize=1)
        ):
            if cell_n % 50 == 0:
                print(
                    f"[allc_to_matrix] progress={100 * cell_n / n_cells:.2f}% "
                    f"cell_index={cell_n}/{n_cells}"
                )
            yield cell_sites


def _csr_from_cells(
    cell_results: Sequence[dict[str, tuple[np.ndarray, np.ndarray]]],
    chrom: str,
    chrom_size: int,
    n_cells: int,
) -> sp_sparse.csr_matrix:
    """Build one chromosome CSR via CSC (columns = cells) then ``tocsr``."""
    counts = np.zeros(n_cells, dtype=np.int64)
    pos_chunks: list[np.ndarray] = []
    val_chunks: list[np.ndarray] = []
    for cell_i, cell_sites in enumerate(cell_results):
        packed = cell_sites.get(chrom)
        if packed is None:
            continue
        pos, val = packed
        if pos.size == 0:
            continue
        counts[cell_i] = pos.size
        pos_chunks.append(pos)
        val_chunks.append(val)
    indptr = np.empty(n_cells + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(counts, out=indptr[1:])
    if pos_chunks:
        indices = np.concatenate(pos_chunks)
        data = np.concatenate(val_chunks)
    else:
        indices = np.array([], dtype=np.int64)
        data = np.array([], dtype=np.int8)
    csc = sp_sparse.csc_matrix(
        (data, indices, indptr),
        shape=(chrom_size + 1, n_cells),
        dtype=np.int8,
    )
    csr = csc.tocsr()
    if csr.data.dtype != np.int8:
        csr = csr.astype(np.int8)
    return csr


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
    threads: int = 1,
    run_info_extra: dict | None = None,
) -> dict[str, Path]:
    """Build CSR store from per-cell ALLC files; return output paths.

    ``chunksize`` is accepted for CLI compatibility and written to
    ``run_info.json``. It does not affect the matrix.
    """
    if len(allc_paths) != len(cell_names):
        raise ValueError("allc_paths and cell_names length mismatch")
    if not allc_paths:
        raise ValueError("no ALLC inputs provided")
    if threads < 1:
        raise ValueError("threads must be >= 1")

    output_dir.mkdir(parents=True, exist_ok=True)
    begin_time = datetime.now(timezone.utc)
    n_cells = len(cell_names)
    exclude_key = tuple(sorted(exclude_contigs or []))
    tasks = [
        (str(path), meth_context, round_sites, exclude_key, main_chroms_only)
        for path in allc_paths
    ]
    cell_results = list(_iter_cell_sites(tasks, threads=threads))
    print("[allc_to_matrix] progress=100.00% cell_read_done=1")

    chrom_sizes: dict[str, int] = {}
    for cell_sites in cell_results:
        for chrom, (pos, _val) in cell_sites.items():
            if pos.size == 0:
                continue
            peak = int(pos.max())
            previous = chrom_sizes.get(chrom)
            if previous is None or peak > previous:
                chrom_sizes[chrom] = peak

    n_obs_cell = np.zeros(n_cells, dtype=np.int64)
    n_meth_cell = np.zeros(n_cells, dtype=np.int64)

    for chrom, chrom_size in chrom_sizes.items():
        print(
            f"[allc_to_matrix] csr_chrom={chrom} rows={chrom_size + 1} "
            f"cols={n_cells}"
        )
        mat = _csr_from_cells(cell_results, chrom, chrom_size, n_cells)
        n_obs_cell += np.asarray(mat.getnnz(axis=0)).ravel()
        n_meth_cell += np.ravel(np.sum(mat > 0, axis=0))
        mat_path = output_dir / f"{chrom}.npz"
        sp_sparse.save_npz(mat_path, mat)

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
        threads=threads,
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
