"""MethSCAn-compatible sparse methylation matrix store (clean-room port)."""

from .allc import barcode_from_allc_path, iter_allc_sites
from .region_matrix import build_region_matrices, parse_bed_regions, resolve_regions_label
from .scan import scan_vmrs
from .smooth import smooth_matrix_store
from .store import build_matrix_store

__all__ = [
    "barcode_from_allc_path",
    "build_matrix_store",
    "build_region_matrices",
    "iter_allc_sites",
    "parse_bed_regions",
    "resolve_regions_label",
    "scan_vmrs",
    "smooth_matrix_store",
]
