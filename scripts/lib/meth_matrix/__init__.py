"""MethSCAn-compatible sparse methylation matrix store (clean-room port)."""

from .allc import barcode_from_allc_path, iter_allc_sites
from .store import build_matrix_store

__all__ = [
    "barcode_from_allc_path",
    "build_matrix_store",
    "iter_allc_sites",
]
