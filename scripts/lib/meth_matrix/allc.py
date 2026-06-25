"""ALLC reader and site encoding for seeksoul-matrix methylation matrices."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path


def barcode_from_allc_path(path: Path) -> str:
    """Derive cell barcode from ``<barcode>_allc.gz`` filename."""
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    if name.endswith("_allc"):
        return name[: -len("_allc")]
    return Path(name).stem


def context_matches(context: str, meth_context: str) -> bool:
    """Return whether an ALLC context trinucleotide passes the filter."""
    normalized = meth_context.strip().upper()
    if normalized in ("ALL", "*"):
        return True
    return context.upper().startswith(normalized)


def encode_site(mc: int, cov: int, *, round_sites: bool) -> int | None:
    """Encode one site as +1 (methylated), -1 (unmethylated), or None (skip)."""
    if cov <= 0:
        return None
    n_meth = mc
    n_unmeth = cov - mc
    if n_meth != 0 and n_unmeth != 0:
        if round_sites:
            if n_meth == n_unmeth:
                return None
            return 1 if n_meth > n_unmeth else -1
        return None
    return 1 if n_meth > 0 else -1


def iter_allc_sites(
    allc_path: Path,
    *,
    meth_context: str = "CG",
    round_sites: bool = False,
    exclude_contigs: set[str] | None = None,
    main_chroms_only: bool = False,
) -> Iterator[tuple[str, int, int]]:
    """Yield ``(chrom, pos, meth_value)`` tuples from one gzipped ALLC file."""
    exclude = exclude_contigs or set()
    opener = gzip.open if str(allc_path).lower().endswith(".gz") else open
    with opener(allc_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            chrom = fields[0]
            if main_chroms_only and not _is_main_chrom(chrom):
                continue
            if chrom in exclude:
                continue
            try:
                pos = int(fields[1])
                mc = int(fields[4])
                cov = int(fields[5])
            except ValueError:
                continue
            context = fields[3]
            if not context_matches(context, meth_context):
                continue
            meth_value = encode_site(mc, cov, round_sites=round_sites)
            if meth_value is None:
                continue
            yield chrom, pos, meth_value


def _is_main_chrom(chrom: str) -> bool:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        return suffix.isdigit() or suffix in {"X", "Y", "M"}
    return False
