#!/usr/bin/env python3
"""Generate and optionally submit seeksoul-matrix workflow commands."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

from _version import __version__
import workflow_input_checks as wic

BASE_STAGE_SEQUENCE = [
    "fastp_split",
    "demux_extract_bc",
    "regroup_shards",
    "bismark_align",
    "bam_sort",
]
POST_BAM_SORT_METHY_ONLY = [
    "count_mapped_reads",
    "estimated_cells",
    "split_bams",
    "merge_fr_bams",
    "bam_to_allc",
    "saturation",
    "qc_summary",
]
POST_BAM_SORT_GEXCB = [
    "split_bams",
    "merge_fr_bams",
    "bam_to_allc",
    "saturation",
    "qc_summary",
]
METH_STAGE_SEQUENCE = [
    "allc_to_matrix",
    "meth_smooth",
    "meth_scan",
]
ALL_STAGE_NAMES = [
    *BASE_STAGE_SEQUENCE,
    *POST_BAM_SORT_METHY_ONLY,
    *METH_STAGE_SEQUENCE,
]
STAGE_CHOICES = [*ALL_STAGE_NAMES, "all"]
SLURM_NEST_STAGE_KEYS = frozenset(ALL_STAGE_NAMES)
STAGE_REQUIRED_FIELDS = {
    "fastp_split": ["r1", "r2", "number_of_split_parts"],
    "demux_extract_bc": ["barcode_whitelist"],
    "regroup_shards": [],
    "bismark_align": ["bismark_ref"],
    "bam_sort": [],
    "count_mapped_reads": [],
    "estimated_cells": [],
    "split_bams": [],
    "merge_fr_bams": [],
    "bam_to_allc": ["genome_fa", "chrom_size_path"],
    "saturation": ["chrom_size_path"],
    "qc_summary": [],
    "allc_to_matrix": [],
    "meth_smooth": [],
    "meth_scan": [],
}
DEFAULT_BARCODE_WHITELIST = "whitelist/DD-MET5/U3CB_methylation.txt.gz"
DEFAULT_EXPECTED_CELL_NUM = 3000
DEFAULT_SPLIT_FASTQ_PREFIX_BASES = 1
DEFAULT_METH_CONTEXT = "CG"
DEFAULT_METH_CHUNKSIZE = 10_000_000
DEFAULT_METH_SMOOTH_BANDWIDTH = 1000
DEFAULT_METH_SCAN_BANDWIDTH = 2000
DEFAULT_METH_SCAN_STEPSIZE = 100
DEFAULT_METH_SCAN_VAR_THRESHOLD = 0.02
DEFAULT_METH_SCAN_MIN_CELLS = 6
DEFAULT_METH_SCAN_BRIDGE_GAPS = 0
DEFAULT_METH_MATRIX_CORES = 8


def build_stage_sequence(settings: dict) -> list[str]:
    mode = settings.get("_barcode_mode") or wic.resolve_barcode_mode(settings)
    if mode == "gexcb":
        sequence = [*BASE_STAGE_SEQUENCE, *POST_BAM_SORT_GEXCB]
    else:
        sequence = [*BASE_STAGE_SEQUENCE, *POST_BAM_SORT_METHY_ONLY]
    if settings.get("run_meth_analysis"):
        sequence = [*sequence, *METH_STAGE_SEQUENCE]
    return sequence


def stage_prefix_map(stage_sequence: list[str]) -> dict[str, str]:
    return {name: f"{index + 1:02d}" for index, name in enumerate(stage_sequence)}


def stage_script_name(
    settings: dict, stage_name: str, *, suffix: str = "sh", chunk_id: str = ""
) -> str:
    sequence = build_stage_sequence(settings)
    if stage_name not in sequence and stage_name in METH_STAGE_SEQUENCE:
        sequence = [*sequence, stage_name]
    prefix = stage_prefix_map(sequence)[stage_name]
    chunk_part = f"_{chunk_id}" if chunk_id else ""
    return f"{prefix}_{stage_name}{chunk_part}.{suffix}"


def resolve_env_executable(name: str) -> str:
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.is_file():
        return str(candidate)
    return name


def normalize_executable_setting(value: str | None, default_name: str) -> str:
    if not value or value == default_name:
        return resolve_env_executable(default_name)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate executable command scripts for seeksoul-matrix workflow."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--workflow-config",
        help="JSON config path for workflow/sample settings.",
    )
    parser.add_argument(
        "--runner",
        choices=["local", "slurm"],
        help="Command target: local shell or slurm sbatch.",
    )
    parser.add_argument(
        "--stage",
        choices=STAGE_CHOICES,
        help="Workflow stage to generate command script for. Default: fastp_split.",
    )
    parser.add_argument("--sample-id", help="Sample identifier.")
    parser.add_argument("--r1", help="Input R1 FASTQ(.gz).")
    parser.add_argument("--r2", help="Input R2 FASTQ(.gz).")
    parser.add_argument("--work-root", help="Work root directory. Default: work.")
    parser.add_argument(
        "--fastp-threads",
        type=int,
        help="Thread count for fastp split. Default: 8.",
    )
    parser.add_argument(
        "--number-of-split-parts",
        type=int,
        help="Value passed to fastp --split.",
    )
    parser.add_argument(
        "--fastp-bin",
        help=(
            "fastp executable path or command name. "
            "Default: fastp from current Python env if available, else fastp."
        ),
    )
    parser.add_argument(
        "--barcode-whitelist",
        help=f"Cell barcode whitelist for demux. Default: {DEFAULT_BARCODE_WHITELIST}.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        help="Hamming distance for demux barcode correction. Default: 1.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        help="gzip level for demux output FASTQ. Default: 6.",
    )
    parser.add_argument(
        "--filter-ch",
        type=int,
        help=(
            "CH chimeric filter threshold for demux (0=disabled). Default: 2."
        ),
    )
    parser.add_argument(
        "--split-fastq-prefix-bases",
        type=int,
        help=(
            "Barcode prefix length for demux sharding (maps to SeekSoul split_fastq). "
            f"Default: {DEFAULT_SPLIT_FASTQ_PREFIX_BASES}."
        ),
    )
    parser.add_argument(
        "--bismark-ref",
        help="Bismark --genome path (parent of Bisulfite_Genome/).",
    )
    parser.add_argument(
        "--bismark-parallel",
        type=int,
        help="Bismark --parallel value. Default: 8.",
    )
    parser.add_argument(
        "--bismark-max-insert",
        type=int,
        help="Bismark -X max insert size. Default: 1000.",
    )
    parser.add_argument(
        "--bismark-bin",
        help=(
            "bismark executable path or command name. "
            "Default: bismark from current Python env if available, else bismark."
        ),
    )
    parser.add_argument(
        "--sort-threads",
        type=int,
        help="Thread count for samtools sort -@ in bam_sort. Default: 6.",
    )
    parser.add_argument(
        "--samtools-bin",
        help=(
            "samtools executable path or command name. "
            "Default: samtools from current Python env if available, else samtools."
        ),
    )
    parser.add_argument(
        "--gexcb",
        help=(
            "RNA filtered barcodes for split_bams (mutually exclusive with "
            "--expected-cell-num)."
        ),
    )
    parser.add_argument(
        "--expected-cell-num",
        type=int,
        help=(
            "Expected cell count for estimated_cells threshold. Default: 3000. "
            "Mutually exclusive with --gexcb."
        ),
    )
    parser.add_argument(
        "--force-cell-num",
        type=int,
        help=(
            "Top N barcodes by aligned_reads for estimated_cells. "
            "When set, overrides --expected-cell-num threshold filtering."
        ),
    )
    parser.add_argument(
        "--split-bams-cores",
        type=int,
        help="CPU cores for split_bams parallel batches. Default: 8.",
    )
    parser.add_argument(
        "--merge-fr-bams-cores",
        type=int,
        help="CPU cores for merge_fr_bams parallel per-barcode merges. Default: 8.",
    )
    parser.add_argument(
        "--genome-fa",
        help="Reference genome FASTA for bam_to_allc.",
    )
    parser.add_argument(
        "--chrom-size-path",
        help="Chromosome sizes BED for bam_to_allc workflow parity.",
    )
    parser.add_argument(
        "--bam-to-allc-cores",
        type=int,
        help="CPU cores for bam_to_allc parallel per-barcode conversion. Default: 8.",
    )
    parser.add_argument(
        "--allcools-tag",
        help="BAM tag for ALLCools UMI correction. Default: UR.",
    )
    parser.add_argument(
        "--allcools-bin",
        help=(
            "allcools executable path or command name. "
            "Default: allcools from current Python env if available, else allcools."
        ),
    )
    parser.add_argument(
        "--saturation-script",
        help="Path to saturation script. Default: scripts/saturation.py.",
    )
    parser.add_argument(
        "--saturation-reads-threshold",
        type=float,
        help="HQ cell reads threshold for saturation stage. Default: 100.",
    )
    parser.add_argument(
        "--saturation-max-cells",
        type=int,
        help="Maximum HQ cells for saturation estimation. Default: 100.",
    )
    parser.add_argument(
        "--saturation-sample-seed",
        type=int,
        help="Random seed for saturation HQ cell sampling. Default: 42.",
    )
    parser.add_argument(
        "--saturation-linear-r2-threshold",
        type=float,
        help=(
            "Linear-fit R^2 above which saturation uses linear extrapolation "
            "instead of the saturation curve. Default: 0.99."
        ),
    )
    parser.add_argument(
        "--qc-summary-script",
        help="Path to qc_summary script. Default: scripts/qc_summary.py.",
    )
    parser.add_argument(
        "--cbcsv",
        help=(
            "Optional methylation ↔ GEX barcode map for qc_summary gex_cb column. "
            "Overrides workflow JSON cbcsv when set."
        ),
    )
    parser.add_argument(
        "--mito-chromosomes",
        help=(
            "Comma-separated mitochondrial contigs for qc_summary mito_CG_mc_rate. "
            "Default: chrM. Overrides workflow JSON mito_chromosomes when set."
        ),
    )
    parser.add_argument(
        "--run-meth-analysis",
        action="store_true",
        default=None,
        help=(
            "Append optional meth analysis stages "
            "(allc_to_matrix, meth_smooth, meth_scan) after qc_summary."
        ),
    )
    parser.add_argument(
        "--allc-to-matrix-script",
        help="Path to allc_to_matrix script. Default: scripts/allc_to_matrix.py.",
    )
    parser.add_argument(
        "--meth-context",
        help=f"ALLC context filter for allc_to_matrix. Default: {DEFAULT_METH_CONTEXT}.",
    )
    parser.add_argument(
        "--meth-chunksize",
        type=int,
        help=f"COO chunk size for allc_to_matrix. Default: {DEFAULT_METH_CHUNKSIZE}.",
    )
    parser.add_argument(
        "--meth-round-sites",
        action="store_true",
        default=None,
        help="Round ambiguous ALLC sites in allc_to_matrix.",
    )
    parser.add_argument(
        "--meth-main-chroms-only",
        action="store_true",
        default=None,
        help="Restrict allc_to_matrix to main chromosomes only.",
    )
    parser.add_argument(
        "--meth-exclude-contigs",
        help="Comma-separated contigs to exclude in allc_to_matrix.",
    )
    parser.add_argument(
        "--meth-smooth-script",
        help="Path to meth_smooth script. Default: scripts/meth_smooth.py.",
    )
    parser.add_argument(
        "--meth-smooth-bandwidth",
        type=int,
        help=f"Smoothing bandwidth for meth_smooth. Default: {DEFAULT_METH_SMOOTH_BANDWIDTH}.",
    )
    parser.add_argument(
        "--meth-smooth-use-weights",
        action="store_true",
        default=None,
        help="Weight methylation sites by log1p(coverage) in meth_smooth.",
    )
    parser.add_argument(
        "--meth-scan-script",
        help="Path to meth_scan script. Default: scripts/meth_scan.py.",
    )
    parser.add_argument(
        "--meth-scan-bandwidth",
        type=int,
        help=f"Sliding-window bandwidth for meth_scan. Default: {DEFAULT_METH_SCAN_BANDWIDTH}.",
    )
    parser.add_argument(
        "--meth-scan-stepsize",
        type=int,
        help=f"Sliding-window step size for meth_scan. Default: {DEFAULT_METH_SCAN_STEPSIZE}.",
    )
    parser.add_argument(
        "--meth-scan-var-threshold",
        type=float,
        help=(
            "Top variable-window fraction for meth_scan VMR merge. "
            f"Default: {DEFAULT_METH_SCAN_VAR_THRESHOLD}."
        ),
    )
    parser.add_argument(
        "--meth-scan-min-cells",
        type=int,
        help=f"Minimum cells per VMR for meth_scan. Default: {DEFAULT_METH_SCAN_MIN_CELLS}.",
    )
    parser.add_argument(
        "--meth-scan-bridge-gaps",
        type=int,
        help=(
            "Merge neighboring VMRs within this gap (bp) in meth_scan. "
            f"Default: {DEFAULT_METH_SCAN_BRIDGE_GAPS} (off)."
        ),
    )
    parser.add_argument(
        "--meth-matrix-cores",
        type=int,
        help=f"CPU threads for meth_scan. Default: {DEFAULT_METH_MATRIX_CORES}.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Submit immediately after generating command file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command and output path without writing files.",
    )
    parser.add_argument(
        "--skip-workdir-input-checks",
        action="store_true",
        help=(
            "Do not require prior-stage outputs under the sample work directory. "
            "For fastp_split, still skips r1/r2 existence checks when set. "
            "When generating --stage all, this is passed to each per-stage subprocess automatically."
        ),
    )
    parser.add_argument("--slurm-partition")
    parser.add_argument("--slurm-mem")
    parser.add_argument("--slurm-cpus-per-task", type=int)
    parser.add_argument("--slurm-output")
    parser.add_argument("--slurm-error")
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_fastp_split_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/fastp_split.py")
    command = [
        sys.executable,
        str(script_path),
        "--r1",
        args.r1,
        "--r2",
        args.r2,
        "--work-path",
        str(sample_work),
        "--fastp-threads",
        str(args.fastp_threads),
        "--number-of-split-parts",
        str(args.number_of_split_parts),
        "--fastp-bin",
        args.fastp_bin,
    ]
    return quoted(command)


def build_demux_chunk_command(
    args: argparse.Namespace, r1_path: Path, r2_path: Path, out_prefix: Path
) -> str:
    fastp_json = r1_path.parent / "fastp.json"
    command = [
        sys.executable,
        "scripts/demux_extract_bc.py",
        str(r1_path),
        str(r2_path),
        "--barcode-whitelist",
        args.barcode_whitelist,
        "--output-prefix",
        str(out_prefix),
        "--barcode-hamming-distance",
        str(args.barcode_hamming_distance),
        "--gzip-level",
        str(args.gzip_level),
        "--filter-ch",
        str(args.filter_ch),
        "--fastp-json",
        str(fastp_json),
        "--split-fastq-prefix-bases",
        str(args.split_fastq_prefix_bases),
    ]
    return quoted(command)


def build_regroup_work_command(sample_work: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/regroup_shards.py",
            "--work-path",
            str(sample_work),
        ]
    )


def build_regroup_prefix_command(sample_work: Path, prefix: str) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/regroup_shards.py",
            "--work-path",
            str(sample_work),
            "--prefix",
            prefix,
        ]
    )


def build_aggregate_ct_command(demux_dir: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/aggregate_ct_qc.py",
            "--demux-dir",
            str(demux_dir),
        ]
    )


def build_bismark_align_work_command(args: argparse.Namespace, sample_work: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/bismark_align.py",
            "--work-path",
            str(sample_work),
            "--bismark-ref",
            args.bismark_ref,
            "--bismark-parallel",
            str(args.bismark_parallel),
            "--bismark-max-insert",
            str(args.bismark_max_insert),
            "--bismark-bin",
            args.bismark_bin,
        ]
    )


def build_bismark_align_chunk_command(
    args: argparse.Namespace,
    sample_work: Path,
    chunk_id: str,
    fwd_r1: Path,
    fwd_r2: Path,
    rev_r1: Path,
    rev_r2: Path,
) -> str:
    align_dir = sample_work / "align"
    return quoted(
        [
            sys.executable,
            "scripts/bismark_align.py",
            "--chunk-id",
            chunk_id,
            "--forward-r1",
            str(fwd_r1),
            "--forward-r2",
            str(fwd_r2),
            "--reverse-r1",
            str(rev_r1),
            "--reverse-r2",
            str(rev_r2),
            "--output-dir",
            str(align_dir),
            "--bismark-ref",
            args.bismark_ref,
            "--bismark-parallel",
            str(args.bismark_parallel),
            "--bismark-max-insert",
            str(args.bismark_max_insert),
            "--bismark-bin",
            args.bismark_bin,
        ]
    )


def discover_bismark_align_chunks(sample_work: Path) -> list[tuple[str, Path, Path, Path, Path]]:
    return wic.discover_demux_align_chunks(sample_work / "demux")


def discover_bam_sort_chunks(sample_work: Path) -> list[tuple[str, Path, Path]]:
    return wic.discover_bismark_pe_bams(sample_work / "align")


def resolve_slurm_chunks(
    *,
    discover: callable,
    plan: callable,
    number_of_split_parts: int | None,
    label: str,
) -> list:
    chunks = discover()
    if chunks:
        return chunks
    if number_of_split_parts is None:
        raise ValueError(
            f"no {label} found; pass --number-of-split-parts for upfront Slurm script generation"
        )
    return plan(number_of_split_parts)


def resolve_prefix_chunks(
    *,
    discover: callable,
    plan_by_prefix: callable,
    base_dir: Path,
    settings: dict,
    label: str,
) -> list:
    chunks = discover()
    if chunks:
        return chunks
    prefix_bases = settings.get("split_fastq_prefix_bases")
    whitelist = settings.get("barcode_whitelist")
    if prefix_bases is None or int(prefix_bases) <= 0:
        raise ValueError(
            f"no {label} found; set split_fastq_prefix_bases and barcode_whitelist "
            "for upfront Slurm generation"
        )
    if not whitelist:
        raise ValueError(
            f"no {label} found; barcode_whitelist is required for prefix chunk planning"
        )
    prefixes = wic.plan_prefix_chunks(
        wic.resolve_config_path(whitelist), int(prefix_bases)
    )
    return plan_by_prefix(base_dir, prefixes)


def build_bam_sort_work_command(args: argparse.Namespace, sample_work: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/bam_sort.py",
            "--work-path",
            str(sample_work),
            "--sort-threads",
            str(args.sort_threads),
            "--samtools-bin",
            args.samtools_bin,
        ]
    )


def build_bam_sort_chunk_command(
    args: argparse.Namespace,
    sample_work: Path,
    chunk_id: str,
    forward_bam: Path,
    reverse_bam: Path,
) -> str:
    align_dir = sample_work / "align"
    return quoted(
        [
            sys.executable,
            "scripts/bam_sort.py",
            "--chunk-id",
            chunk_id,
            "--forward-bam",
            str(forward_bam),
            "--reverse-bam",
            str(reverse_bam),
            "--output-dir",
            str(align_dir),
            "--sort-threads",
            str(args.sort_threads),
            "--samtools-bin",
            args.samtools_bin,
        ]
    )


def build_count_mapped_reads_work_command(sample_work: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/count_mapped_reads.py",
            "--work-path",
            str(sample_work),
        ]
    )


def build_count_mapped_reads_chunk_command(sample_work: Path, chunk_id: str) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/count_mapped_reads.py",
            "--work-path",
            str(sample_work),
            "--chunk-id",
            chunk_id,
        ]
    )


def build_saturation_command(args: argparse.Namespace, sample_work: Path) -> str:
    fastp_json = sample_work / "shard_fastq" / "fastp.json"
    command = [
        sys.executable,
        str(args.saturation_script),
        "--work-path",
        str(sample_work),
        "--fastp-json",
        str(fastp_json),
        "--chrom-size-path",
        str(args.chrom_size_path),
        "--reads-threshold",
        str(args.saturation_reads_threshold),
        "--max-cells",
        str(args.saturation_max_cells),
        "--sample-seed",
        str(args.saturation_sample_seed),
        "--linear-r2-threshold",
        str(args.saturation_linear_r2_threshold),
    ]
    return quoted(command)


def build_qc_summary_command(
    args: argparse.Namespace,
    sample_work: Path,
    *,
    sample_id: str,
    barcode_mode: str,
) -> str:
    command = [
        sys.executable,
        str(args.qc_summary_script),
        "--work-path",
        str(sample_work),
        "--sample-id",
        sample_id,
        "--barcode-mode",
        "gexcb" if barcode_mode == "gexcb" else "methylation_only",
    ]
    if getattr(args, "cbcsv", None):
        command.extend(["--cbcsv", str(args.cbcsv)])
    mito_chromosomes = getattr(args, "mito_chromosomes", None) or "chrM"
    command.extend(["--mito-chromosomes", str(mito_chromosomes)])
    return quoted(command)


def build_allc_to_matrix_command(
    args: argparse.Namespace,
    sample_work: Path,
    *,
    barcode_mode: str,
) -> str:
    command = [
        sys.executable,
        str(args.allc_to_matrix_script),
        "--work-path",
        str(sample_work),
        "--barcode-mode",
        "gexcb" if barcode_mode == "gexcb" else "methylation_only",
        "--meth-context",
        str(args.meth_context),
        "--chunksize",
        str(args.meth_chunksize),
    ]
    if args.meth_round_sites:
        command.append("--round-sites")
    if args.meth_main_chroms_only:
        command.append("--main-chroms-only")
    if getattr(args, "meth_exclude_contigs", None):
        command.extend(["--exclude-contigs", str(args.meth_exclude_contigs)])
    return quoted(command)


def build_meth_smooth_command(args: argparse.Namespace, sample_work: Path) -> str:
    command = [
        sys.executable,
        str(args.meth_smooth_script),
        "--work-path",
        str(sample_work),
        "--bandwidth",
        str(args.meth_smooth_bandwidth),
    ]
    if args.meth_smooth_use_weights:
        command.append("--use-weights")
    return quoted(command)


def build_meth_scan_command(args: argparse.Namespace, sample_work: Path) -> str:
    command = [
        sys.executable,
        str(args.meth_scan_script),
        "--work-path",
        str(sample_work),
        "--bandwidth",
        str(args.meth_scan_bandwidth),
        "--stepsize",
        str(args.meth_scan_stepsize),
        "--var-threshold",
        str(args.meth_scan_var_threshold),
        "--min-cells",
        str(args.meth_scan_min_cells),
        "--bridge-gaps",
        str(args.meth_scan_bridge_gaps),
        "--threads",
        str(args.meth_matrix_cores),
    ]
    return quoted(command)


def build_estimated_cells_command(
    sample_work: Path,
    expected_cell_num: int,
    force_cell_num: int | None = None,
) -> str:
    command = [
        sys.executable,
        "scripts/estimated_cells.py",
        "--work-path",
        str(sample_work),
    ]
    if force_cell_num is not None:
        command.extend(["--force-cell-num", str(force_cell_num)])
    else:
        command.extend(["--expected-cell-num", str(expected_cell_num)])
    return quoted(command)


def build_split_bams_work_command(args: argparse.Namespace, sample_work: Path) -> str:
    command = [
        sys.executable,
        "scripts/split_bams.py",
        "--work-path",
        str(sample_work),
        "--cores",
        str(args.split_bams_cores),
    ]
    if args.gexcb:
        command.extend(["--gexcb", args.gexcb])
    else:
        command.extend(
            [
                "--filtered-barcode",
                str(sample_work / "cells" / "filtered_barcode"),
            ]
        )
    return quoted(command)


def build_split_bams_chunk_command(
    args: argparse.Namespace,
    sample_work: Path,
    chunk_id: str,
) -> str:
    command = [
        sys.executable,
        "scripts/split_bams.py",
        "--work-path",
        str(sample_work),
        "--chunk-id",
        chunk_id,
        "--cores",
        str(args.split_bams_cores),
    ]
    if args.gexcb:
        command.extend(["--gexcb", args.gexcb])
    else:
        command.extend(
            [
                "--filtered-barcode",
                str(sample_work / "cells" / "filtered_barcode"),
            ]
        )
    return quoted(command)


def build_merge_fr_bams_work_command(
    args: argparse.Namespace,
    sample_work: Path,
) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/merge_fr_bams.py",
            "--work-path",
            str(sample_work),
            "--cores",
            str(args.merge_fr_bams_cores),
            "--samtools-bin",
            args.samtools_bin,
        ]
    )


def build_merge_fr_bams_chunk_command(
    args: argparse.Namespace,
    sample_work: Path,
    chunk_id: str,
) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/merge_fr_bams.py",
            "--work-path",
            str(sample_work),
            "--chunk-id",
            chunk_id,
            "--cores",
            str(args.merge_fr_bams_cores),
            "--samtools-bin",
            args.samtools_bin,
        ]
    )


def build_bam_to_allc_work_command(
    args: argparse.Namespace,
    sample_work: Path,
) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/bam_to_allc.py",
            "--work-path",
            str(sample_work),
            "--genome-fa",
            args.genome_fa,
            "--chrom-size-path",
            args.chrom_size_path,
            "--cores",
            str(args.bam_to_allc_cores),
            "--allcools-tag",
            args.allcools_tag,
            "--samtools-bin",
            args.samtools_bin,
            "--allcools-bin",
            args.allcools_bin,
        ]
    )


def build_bam_to_allc_chunk_command(
    args: argparse.Namespace,
    sample_work: Path,
    chunk_id: str,
) -> str:
    return quoted(
        [
            sys.executable,
            "scripts/bam_to_allc.py",
            "--work-path",
            str(sample_work),
            "--chunk-id",
            chunk_id,
            "--genome-fa",
            args.genome_fa,
            "--chrom-size-path",
            args.chrom_size_path,
            "--cores",
            str(args.bam_to_allc_cores),
            "--allcools-tag",
            args.allcools_tag,
            "--samtools-bin",
            args.samtools_bin,
            "--allcools-bin",
            args.allcools_bin,
        ]
    )


def build_demux_chunks_from_config(
    sample_work: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path, Path]]:
    shard_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    return [
        (chunk_id, r1_path, r2_path, demux_dir / chunk_id)
        for chunk_id, r1_path, r2_path in wic.plan_fastp_shards(
            shard_dir, number_of_split_parts
        )
    ]


def build_demux_local_batch_command(args: argparse.Namespace, sample_work: Path) -> str:
    chunk_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    fastp_json = chunk_dir / "fastp.json"
    chunk_dir_q = shlex.quote(str(chunk_dir))
    demux_dir_q = shlex.quote(str(demux_dir))
    fastp_json_q = shlex.quote(str(fastp_json))
    py = quoted([sys.executable, "scripts/demux_extract_bc.py"])
    aggregate_cmd = build_aggregate_ct_command(demux_dir)
    return (
        f"chunk_dir={chunk_dir_q}\n"
        f"demux_dir={demux_dir_q}\n"
        f"fastp_json={fastp_json_q}\n"
        "\n"
        "mkdir -p \"$demux_dir\"\n"
        "shopt -s nullglob\n"
        "r1_files=(\"$chunk_dir\"/*.R1.fq.gz)\n"
        "total=${#r1_files[@]}\n"
        "if [ \"$total\" -eq 0 ]; then\n"
        '  echo "[demux] no chunk found under shard_fastq"\n'
        "  exit 1\n"
        "fi\n"
        "\n"
        "idx=0\n"
        "for r1 in \"${r1_files[@]}\"; do\n"
        "  idx=$((idx + 1))\n"
        "  percent=$((idx * 100 / total))\n"
        '  chunk="$(basename "$r1" .R1.fq.gz)"\n'
        "  r2=\"$chunk_dir/${chunk}.R2.fq.gz\"\n"
        "  printf '[demux] %3d%% (%d/%d) %s\\n' \"$percent\" \"$idx\" \"$total\" \"$chunk\"\n"
        '  [ -f "$r2" ] || { echo "[demux] missing pair for $r1: $r2"; exit 1; }\n'
        f"  {py} "
        '"$r1" "$r2" '
        f"--barcode-whitelist {shlex.quote(args.barcode_whitelist)} "
        '--output-prefix "$demux_dir/$chunk" '
        f"--barcode-hamming-distance {int(args.barcode_hamming_distance)} "
        f"--gzip-level {int(args.gzip_level)} "
        f"--filter-ch {int(args.filter_ch)} "
        f"--split-fastq-prefix-bases {int(args.split_fastq_prefix_bases)} "
        '--fastp-json "$fastp_json"\n'
        "done\n"
        "\n"
        f"{aggregate_cmd}\n"
        "\n"
        'echo "[demux] done"'
    )


def build_demux_slurm_submit_command(sample_work: Path) -> str:
    command_dir = sample_work / "commands"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f'SCRIPT_DIR={shlex.quote(str(command_dir))}',
        'job_ids=""',
        'for script in "$SCRIPT_DIR"/02_demux_extract_bc_*.sbatch; do',
        '  [ -f "$script" ] || continue',
        '  jid=$(sbatch --parsable "$script")',
        '  job_ids="${job_ids}:${jid}"',
        "done",
        'if [ -z "$job_ids" ]; then',
        '  echo "[demux] no chunk sbatch scripts found"',
        "  exit 1",
        "fi",
        'sbatch --dependency=afterok"${job_ids}" "$SCRIPT_DIR/02_aggregate_ct_qc.sbatch"',
        'echo "[demux] submitted aggregate job with dependency afterok${job_ids}"',
    ]
    return "\n".join(lines)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def load_workflow_config(path: str) -> dict:
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("workflow config must be a JSON object")
    return data


def pick(cli_value, cfg_value):
    return cli_value if cli_value is not None else cfg_value


def validate_required_for_stage(stage: str, settings: dict) -> None:
    required = ["runner", "sample_id", *STAGE_REQUIRED_FIELDS[stage]]
    missing = [key for key in required if settings.get(key) in (None, "")]
    if missing:
        raise ValueError(f"missing required settings: {', '.join(missing)}")


def validate_inputs_for_stage(
    stage: str,
    settings: dict,
    sample_work: Path,
    *,
    skip_workdir_inputs: bool = False,
) -> None:
    if skip_workdir_inputs:
        return
    if stage == "fastp_split":
        wic.require_file("r1", wic.resolve_config_path(settings["r1"]))
        wic.require_file("r2", wic.resolve_config_path(settings["r2"]))
        wic.require_optional_executable_path("fastp_bin", settings["fastp_bin"])
    elif stage == "demux_extract_bc":
        wic.require_file(
            "barcode_whitelist",
            wic.resolve_config_path(settings["barcode_whitelist"]),
        )
        if not skip_workdir_inputs:
            shard_dir = sample_work / "shard_fastq"
            shards = wic.discover_fastp_shards(shard_dir)
            if not shards:
                raise ValueError(f"no fastp shards found under {shard_dir}")
            for _chunk_id, r1_path, r2_path in shards:
                wic.require_file(f"shard_fastq/{r1_path.name}", r1_path)
                wic.require_file(f"shard_fastq/{r2_path.name}", r2_path)
    elif stage == "regroup_shards":
        if int(settings.get("split_fastq_prefix_bases") or 0) <= 0:
            return
        demux_dir = sample_work / "demux"
        subshards = wic.discover_demux_subshards(demux_dir)
        if not subshards:
            raise ValueError(f"no demux sub-shards found under {demux_dir}/shards")
    elif stage == "bismark_align":
        wic.require_bismark_ref(wic.resolve_config_path(settings["bismark_ref"]))
        wic.require_optional_executable_path("bismark_bin", settings["bismark_bin"])
        demux_dir = sample_work / "demux"
        chunks = wic.discover_demux_align_chunks(demux_dir)
        if not chunks:
            raise ValueError(f"no demux align inputs found under {demux_dir}")
        for _chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2 in chunks:
            wic.require_file(f"demux/{fwd_r1.name}", fwd_r1)
            wic.require_file(f"demux/{fwd_r2.name}", fwd_r2)
            wic.require_file(f"demux/{rev_r1.name}", rev_r1)
            wic.require_file(f"demux/{rev_r2.name}", rev_r2)
    elif stage == "bam_sort":
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
        align_dir = sample_work / "align"
        chunks = wic.discover_bismark_pe_bams(align_dir)
        if not chunks:
            raise ValueError(f"no Bismark PE BAMs found under {align_dir}")
        for _chunk_id, fwd_bam, rev_bam in chunks:
            wic.require_file(f"align/{fwd_bam.name}", fwd_bam)
            wic.require_file(f"align/{rev_bam.name}", rev_bam)
    elif stage == "count_mapped_reads":
        align_dir = sample_work / "align"
        chunks = wic.discover_bismark_pe_bams(align_dir)
        if not chunks:
            raise ValueError(f"no unsorted Bismark PE BAMs found under {align_dir}")
        for _chunk_id, fwd_bam, rev_bam in chunks:
            wic.require_file(f"align/{fwd_bam.name}", fwd_bam)
            wic.require_file(f"align/{rev_bam.name}", rev_bam)
    elif stage == "estimated_cells":
        align_dir = sample_work / "align"
        count_files = list(align_dir.glob("*_cb_aligned_reads_counts.csv"))
        if not count_files:
            raise ValueError(
                f"no *_cb_aligned_reads_counts.csv files found under {align_dir}"
            )
    elif stage == "split_bams":
        align_dir = sample_work / "align"
        chunks = wic.discover_bismark_sortbyname_bams(align_dir)
        if not chunks:
            raise ValueError(f"no sortbyname Bismark PE BAMs found under {align_dir}")
        for _chunk_id, fwd_bam, rev_bam in chunks:
            wic.require_file(f"align/{fwd_bam.name}", fwd_bam)
            wic.require_file(f"align/{rev_bam.name}", rev_bam)
        mode = settings.get("_barcode_mode") or wic.resolve_barcode_mode(settings)
        if mode == "gexcb":
            wic.require_file("gexcb", wic.resolve_config_path(settings["gexcb"]))
        else:
            filtered_barcode = sample_work / "cells" / "filtered_barcode"
            wic.require_file("cells/filtered_barcode", filtered_barcode)
    elif stage == "merge_fr_bams":
        split_root = sample_work / "split_bams"
        pairs = wic.discover_split_bam_chunk_pairs(split_root)
        if not pairs:
            raise ValueError(f"no split BAM chunk pairs found under {split_root}")
        for _chunk_id, forward_dir, reverse_dir in pairs:
            wic.require_dir(f"split_bams/{forward_dir.name}", forward_dir)
            wic.require_dir(f"split_bams/{reverse_dir.name}", reverse_dir)
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "bam_to_allc":
        merged_root = sample_work / "split_bams" / "merged"
        chunks = wic.discover_merged_fr_bam_chunks(merged_root)
        if not chunks:
            raise ValueError(f"no merged FR BAM chunks found under {merged_root}")
        for _chunk_id, bam_dir, filtered_barcode in chunks:
            wic.require_dir(f"split_bams/merged/{bam_dir.name}", bam_dir)
            wic.require_file(
                f"split_bams/merged/{filtered_barcode.name}",
                filtered_barcode,
            )
        wic.require_file(
            "genome_fa",
            wic.resolve_config_path(settings["genome_fa"]),
        )
        wic.require_file(
            "chrom_size_path",
            wic.resolve_config_path(settings["chrom_size_path"]),
        )
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
        wic.require_optional_executable_path("allcools_bin", settings["allcools_bin"])
    elif stage == "saturation":
        script_path = Path(settings["saturation_script"])
        if not script_path.is_file():
            raise FileNotFoundError(f"saturation_script not found: {script_path}")
        wic.require_file(
            "chrom_size_path",
            wic.resolve_config_path(settings["chrom_size_path"]),
        )
        merged_root = sample_work / "split_bams" / "merged"
        bam_files = list(merged_root.glob("*_merged_fr_bam/*.bam"))
        if not bam_files:
            raise FileNotFoundError(
                f"no per-cell BAM files found under {merged_root}/*_merged_fr_bam/"
            )
    elif stage == "qc_summary":
        script_path = Path(settings["qc_summary_script"])
        if not script_path.is_file():
            raise FileNotFoundError(f"qc_summary_script not found: {script_path}")
        saturation_summary = sample_work / "qc" / "saturation" / "saturation_summary.tsv"
        if not saturation_summary.is_file():
            raise FileNotFoundError(
                f"saturation summary not found: {saturation_summary}"
            )
        allcools_dir = sample_work / "allcools"
        if not allcools_dir.is_dir():
            raise FileNotFoundError(f"allcools directory not found: {allcools_dir}")
        count_files = list(allcools_dir.glob("*/*_allc.gz.count.csv"))
        if not count_files:
            raise FileNotFoundError(
                f"no per-cell ALLC count files found under {allcools_dir}"
            )
        mode = settings.get("_barcode_mode") or wic.resolve_barcode_mode(settings)
        if mode == "gexcb":
            gexcb_reads = list(
                sample_work.glob(
                    "split_bams/merged/*_merge_filtered_barcode_reads_counts.csv"
                )
            )
            if not gexcb_reads:
                raise FileNotFoundError(
                    "no gexcb merge read-count tables found under "
                    f"{sample_work}/split_bams/merged/"
                )
        else:
            reads_path = sample_work / "cells" / "filtered_barcode_read_counts.csv"
            if not reads_path.is_file():
                raise FileNotFoundError(f"cell reads table not found: {reads_path}")
    elif stage == "allc_to_matrix":
        script_path = Path(settings["allc_to_matrix_script"])
        if not script_path.is_file():
            raise FileNotFoundError(f"allc_to_matrix_script not found: {script_path}")
        allcools_dir = sample_work / "allcools"
        allc_files = list(allcools_dir.glob("*_merged_fr_bam_allcools/*_allc.gz"))
        if not allc_files:
            raise FileNotFoundError(
                f"no per-cell ALLC files found under {allcools_dir}"
            )
    elif stage == "meth_smooth":
        script_path = Path(settings["meth_smooth_script"])
        if not script_path.is_file():
            raise FileNotFoundError(f"meth_smooth_script not found: {script_path}")
        matrix_dir = sample_work / "meth" / "matrix"
        npz_files = list(matrix_dir.glob("*.npz"))
        if not npz_files:
            raise FileNotFoundError(
                f"no CSR matrix files found under {matrix_dir}"
            )
    elif stage == "meth_scan":
        script_path = Path(settings["meth_scan_script"])
        if not script_path.is_file():
            raise FileNotFoundError(f"meth_scan_script not found: {script_path}")
        matrix_dir = sample_work / "meth" / "matrix"
        npz_files = list(matrix_dir.glob("*.npz"))
        if not npz_files:
            raise FileNotFoundError(
                f"no CSR matrix files found under {matrix_dir}"
            )
        smoothed_dir = matrix_dir / "smoothed"
        smoothed_files = list(smoothed_dir.glob("*.csv.gz")) + list(
            smoothed_dir.glob("*.csv")
        )
        if not smoothed_files:
            raise FileNotFoundError(
                f"no smoothed chromosome files found under {smoothed_dir}"
            )
    else:
        raise ValueError(f"unsupported stage for input validation: {stage}")


def select_stage_slurm_cfg(slurm_cfg_raw: dict, stage: str) -> dict:
    if any(key in slurm_cfg_raw for key in SLURM_NEST_STAGE_KEYS):
        stage_slurm_cfg = slurm_cfg_raw.get(stage, {})
    else:
        stage_slurm_cfg = slurm_cfg_raw
    if stage_slurm_cfg is None:
        stage_slurm_cfg = {}
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    return stage_slurm_cfg


def resolve_settings(args: argparse.Namespace) -> dict:
    cfg: dict = {}
    if args.workflow_config:
        cfg = load_workflow_config(args.workflow_config)

    slurm_cfg_raw = cfg.get("slurm", {})
    if slurm_cfg_raw is None:
        slurm_cfg_raw = {}
    if not isinstance(slurm_cfg_raw, dict):
        raise ValueError("workflow config key 'slurm' must be an object")

    stage = pick(args.stage, cfg.get("stage")) or "fastp_split"
    if stage not in STAGE_CHOICES:
        raise ValueError(f"unsupported stage: {stage}")

    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, stage)
    settings = {
        "runner": pick(args.runner, cfg.get("runner")),
        "stage": stage,
        "sample_id": pick(args.sample_id, cfg.get("sample_id")),
        "r1": pick(args.r1, cfg.get("r1")),
        "r2": pick(args.r2, cfg.get("r2")),
        "work_root": pick(args.work_root, cfg.get("work_root")),
        "fastp_threads": pick(args.fastp_threads, cfg.get("fastp_threads")),
        "number_of_split_parts": pick(
            args.number_of_split_parts, cfg.get("number_of_split_parts")
        ),
        "fastp_bin": pick(args.fastp_bin, cfg.get("fastp_bin")),
        "barcode_whitelist": pick(args.barcode_whitelist, cfg.get("barcode_whitelist")),
        "barcode_hamming_distance": pick(
            args.barcode_hamming_distance, cfg.get("barcode_hamming_distance")
        ),
        "gzip_level": pick(args.gzip_level, cfg.get("gzip_level")),
        "filter_ch": pick(args.filter_ch, cfg.get("filter_ch")),
        "split_fastq_prefix_bases": pick(
            args.split_fastq_prefix_bases, cfg.get("split_fastq_prefix_bases")
        ),
        "bismark_ref": pick(args.bismark_ref, cfg.get("bismark_ref")),
        "bismark_parallel": pick(args.bismark_parallel, cfg.get("bismark_parallel")),
        "bismark_max_insert": pick(
            args.bismark_max_insert, cfg.get("bismark_max_insert")
        ),
        "bismark_bin": pick(args.bismark_bin, cfg.get("bismark_bin")),
        "sort_threads": pick(args.sort_threads, cfg.get("sort_threads")),
        "samtools_bin": pick(args.samtools_bin, cfg.get("samtools_bin")),
        "gexcb": pick(args.gexcb, cfg.get("gexcb")),
        "expected_cell_num": pick(args.expected_cell_num, cfg.get("expected_cell_num")),
        "force_cell_num": pick(args.force_cell_num, cfg.get("force_cell_num")),
        "split_bams_cores": pick(args.split_bams_cores, cfg.get("split_bams_cores")),
        "merge_fr_bams_cores": pick(
            args.merge_fr_bams_cores, cfg.get("merge_fr_bams_cores")
        ),
        "genome_fa": pick(args.genome_fa, cfg.get("genome_fa")),
        "chrom_size_path": pick(args.chrom_size_path, cfg.get("chrom_size_path")),
        "bam_to_allc_cores": pick(args.bam_to_allc_cores, cfg.get("bam_to_allc_cores")),
        "allcools_tag": pick(args.allcools_tag, cfg.get("allcools_tag")),
        "allcools_bin": pick(args.allcools_bin, cfg.get("allcools_bin")),
        "saturation_script": pick(args.saturation_script, cfg.get("saturation_script")),
        "saturation_reads_threshold": pick(
            args.saturation_reads_threshold,
            cfg.get("saturation_reads_threshold"),
        ),
        "saturation_max_cells": pick(
            args.saturation_max_cells,
            cfg.get("saturation_max_cells"),
        ),
        "saturation_sample_seed": pick(
            args.saturation_sample_seed,
            cfg.get("saturation_sample_seed"),
        ),
        "saturation_linear_r2_threshold": pick(
            args.saturation_linear_r2_threshold,
            cfg.get("saturation_linear_r2_threshold"),
        ),
        "qc_summary_script": pick(args.qc_summary_script, cfg.get("qc_summary_script")),
        "cbcsv": pick(args.cbcsv, cfg.get("cbcsv")),
        "mito_chromosomes": pick(args.mito_chromosomes, cfg.get("mito_chromosomes")),
        "run_meth_analysis": (
            args.run_meth_analysis
            if args.run_meth_analysis is not None
            else bool(cfg.get("run_meth_analysis"))
        ),
        "allc_to_matrix_script": pick(
            args.allc_to_matrix_script, cfg.get("allc_to_matrix_script")
        ),
        "meth_context": pick(args.meth_context, cfg.get("meth_context")),
        "meth_chunksize": pick(args.meth_chunksize, cfg.get("meth_chunksize")),
        "meth_round_sites": (
            args.meth_round_sites
            if args.meth_round_sites is not None
            else bool(cfg.get("meth_round_sites"))
        ),
        "meth_main_chroms_only": (
            args.meth_main_chroms_only
            if args.meth_main_chroms_only is not None
            else bool(cfg.get("meth_main_chroms_only"))
        ),
        "meth_exclude_contigs": pick(
            args.meth_exclude_contigs, cfg.get("meth_exclude_contigs")
        ),
        "meth_smooth_script": pick(args.meth_smooth_script, cfg.get("meth_smooth_script")),
        "meth_smooth_bandwidth": pick(
            args.meth_smooth_bandwidth, cfg.get("meth_smooth_bandwidth")
        ),
        "meth_smooth_use_weights": (
            args.meth_smooth_use_weights
            if args.meth_smooth_use_weights is not None
            else bool(cfg.get("meth_smooth_use_weights"))
        ),
        "meth_scan_script": pick(args.meth_scan_script, cfg.get("meth_scan_script")),
        "meth_scan_bandwidth": pick(
            args.meth_scan_bandwidth, cfg.get("meth_scan_bandwidth")
        ),
        "meth_scan_stepsize": pick(args.meth_scan_stepsize, cfg.get("meth_scan_stepsize")),
        "meth_scan_var_threshold": pick(
            args.meth_scan_var_threshold, cfg.get("meth_scan_var_threshold")
        ),
        "meth_scan_min_cells": pick(
            args.meth_scan_min_cells, cfg.get("meth_scan_min_cells")
        ),
        "meth_scan_bridge_gaps": pick(
            args.meth_scan_bridge_gaps, cfg.get("meth_scan_bridge_gaps")
        ),
        "meth_matrix_cores": pick(args.meth_matrix_cores, cfg.get("meth_matrix_cores")),
        "slurm_partition": pick(args.slurm_partition, stage_slurm_cfg.get("partition")),
        "slurm_mem": pick(args.slurm_mem, stage_slurm_cfg.get("mem")),
        "slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, stage_slurm_cfg.get("cpus_per_task")
        ),
        "slurm_output": pick(args.slurm_output, stage_slurm_cfg.get("output")),
        "slurm_error": pick(args.slurm_error, stage_slurm_cfg.get("error")),
        "submit": args.submit,
        "dry_run": args.dry_run,
        "_slurm_cfg_raw": slurm_cfg_raw,
    }

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = settings["fastp_threads"] or 8
    if settings["number_of_split_parts"] is not None:
        settings["number_of_split_parts"] = int(settings["number_of_split_parts"])
        if settings["number_of_split_parts"] <= 0:
            raise ValueError("number_of_split_parts must be > 0")
    settings["fastp_bin"] = normalize_executable_setting(
        settings["fastp_bin"], "fastp"
    )
    if stage in ("demux_extract_bc", "regroup_shards", "all"):
        settings["barcode_whitelist"] = (
            settings["barcode_whitelist"] or DEFAULT_BARCODE_WHITELIST
        )
    if stage in ("demux_extract_bc", "regroup_shards", "bismark_align", "all"):
        settings["split_fastq_prefix_bases"] = int(
            settings["split_fastq_prefix_bases"]
            if settings["split_fastq_prefix_bases"] is not None
            else DEFAULT_SPLIT_FASTQ_PREFIX_BASES
        )
        if settings["split_fastq_prefix_bases"] < 0:
            raise ValueError("split_fastq_prefix_bases must be >= 0")
    if stage in ("demux_extract_bc", "all"):
        settings["barcode_hamming_distance"] = int(
            settings["barcode_hamming_distance"] or 1
        )
        settings["gzip_level"] = int(settings["gzip_level"] or 6)
        settings["filter_ch"] = int(settings["filter_ch"] if settings["filter_ch"] is not None else 2)
    if stage == "bismark_align" or stage == "all":
        settings["bismark_parallel"] = int(settings["bismark_parallel"] or 8)
        settings["bismark_max_insert"] = int(settings["bismark_max_insert"] or 1000)
        settings["bismark_bin"] = normalize_executable_setting(
            settings["bismark_bin"], "bismark"
        )
    if stage == "bam_sort" or stage == "all":
        settings["sort_threads"] = int(settings["sort_threads"] or 6)
        settings["samtools_bin"] = normalize_executable_setting(
            settings["samtools_bin"], "samtools"
        )
    settings["_barcode_mode"] = wic.resolve_barcode_mode(settings)
    if settings["_barcode_mode"] == "expected_cell_num":
        if settings["expected_cell_num"] is None:
            settings["expected_cell_num"] = DEFAULT_EXPECTED_CELL_NUM
        settings["expected_cell_num"] = int(settings["expected_cell_num"])
    if settings.get("force_cell_num") is not None:
        settings["force_cell_num"] = int(settings["force_cell_num"])
        if settings["force_cell_num"] <= 0:
            raise ValueError("force_cell_num must be > 0")
    if stage in ("split_bams", "all") or settings["_barcode_mode"] == "gexcb":
        settings["split_bams_cores"] = int(settings["split_bams_cores"] or 8)
    if stage in ("merge_fr_bams", "all"):
        settings["merge_fr_bams_cores"] = int(settings["merge_fr_bams_cores"] or 8)
        settings["samtools_bin"] = normalize_executable_setting(
            settings["samtools_bin"], "samtools"
        )
    if stage in ("bam_to_allc", "all"):
        settings["bam_to_allc_cores"] = int(settings["bam_to_allc_cores"] or 8)
        settings["allcools_tag"] = settings["allcools_tag"] or "UR"
        settings["samtools_bin"] = normalize_executable_setting(
            settings["samtools_bin"], "samtools"
        )
        settings["allcools_bin"] = normalize_executable_setting(
            settings["allcools_bin"], "allcools"
        )
    if stage in ("saturation", "all"):
        settings["saturation_script"] = (
            settings["saturation_script"] or "scripts/saturation.py"
        )
        settings["saturation_reads_threshold"] = (
            float(settings["saturation_reads_threshold"])
            if settings["saturation_reads_threshold"] is not None
            else 100.0
        )
        settings["saturation_max_cells"] = int(
            settings["saturation_max_cells"]
            if settings["saturation_max_cells"] is not None
            else 100
        )
        settings["saturation_sample_seed"] = int(
            settings["saturation_sample_seed"]
            if settings["saturation_sample_seed"] is not None
            else 42
        )
        settings["saturation_linear_r2_threshold"] = (
            float(settings["saturation_linear_r2_threshold"])
            if settings["saturation_linear_r2_threshold"] is not None
            else 0.99
        )
        if settings["saturation_reads_threshold"] <= 0:
            raise ValueError("saturation_reads_threshold must be > 0")
        if settings["saturation_max_cells"] <= 0:
            raise ValueError("saturation_max_cells must be > 0")
    if stage in ("qc_summary", "all"):
        settings["qc_summary_script"] = (
            settings["qc_summary_script"] or "scripts/qc_summary.py"
        )
        settings["mito_chromosomes"] = settings["mito_chromosomes"] or "chrM"
    if (
        stage in ("allc_to_matrix", "meth_smooth", "meth_scan", "all")
        or settings.get("run_meth_analysis")
    ):
        settings["allc_to_matrix_script"] = (
            settings["allc_to_matrix_script"] or "scripts/allc_to_matrix.py"
        )
        settings["meth_context"] = settings["meth_context"] or DEFAULT_METH_CONTEXT
        settings["meth_chunksize"] = int(
            settings["meth_chunksize"]
            if settings["meth_chunksize"] is not None
            else DEFAULT_METH_CHUNKSIZE
        )
        if settings["meth_chunksize"] <= 0:
            raise ValueError("meth_chunksize must be > 0")
        settings["meth_round_sites"] = bool(settings["meth_round_sites"])
        settings["meth_main_chroms_only"] = bool(settings["meth_main_chroms_only"])
        settings["meth_exclude_contigs"] = settings["meth_exclude_contigs"] or ""
        settings["meth_smooth_script"] = (
            settings["meth_smooth_script"] or "scripts/meth_smooth.py"
        )
        settings["meth_smooth_bandwidth"] = int(
            settings["meth_smooth_bandwidth"]
            if settings["meth_smooth_bandwidth"] is not None
            else DEFAULT_METH_SMOOTH_BANDWIDTH
        )
        if settings["meth_smooth_bandwidth"] < 1:
            raise ValueError("meth_smooth_bandwidth must be >= 1")
        settings["meth_smooth_use_weights"] = bool(settings["meth_smooth_use_weights"])
        settings["meth_scan_script"] = settings["meth_scan_script"] or "scripts/meth_scan.py"
        settings["meth_scan_bandwidth"] = int(
            settings["meth_scan_bandwidth"]
            if settings["meth_scan_bandwidth"] is not None
            else DEFAULT_METH_SCAN_BANDWIDTH
        )
        settings["meth_scan_stepsize"] = int(
            settings["meth_scan_stepsize"]
            if settings["meth_scan_stepsize"] is not None
            else DEFAULT_METH_SCAN_STEPSIZE
        )
        settings["meth_scan_var_threshold"] = float(
            settings["meth_scan_var_threshold"]
            if settings["meth_scan_var_threshold"] is not None
            else DEFAULT_METH_SCAN_VAR_THRESHOLD
        )
        settings["meth_scan_min_cells"] = int(
            settings["meth_scan_min_cells"]
            if settings["meth_scan_min_cells"] is not None
            else DEFAULT_METH_SCAN_MIN_CELLS
        )
        settings["meth_scan_bridge_gaps"] = int(
            settings["meth_scan_bridge_gaps"]
            if settings["meth_scan_bridge_gaps"] is not None
            else DEFAULT_METH_SCAN_BRIDGE_GAPS
        )
        settings["meth_matrix_cores"] = int(
            settings["meth_matrix_cores"]
            if settings["meth_matrix_cores"] is not None
            else DEFAULT_METH_MATRIX_CORES
        )
        if settings["meth_scan_bandwidth"] < 1:
            raise ValueError("meth_scan_bandwidth must be >= 1")
        if settings["meth_scan_stepsize"] < 1:
            raise ValueError("meth_scan_stepsize must be >= 1")
        if not 0 <= settings["meth_scan_var_threshold"] <= 1:
            raise ValueError("meth_scan_var_threshold must be between 0 and 1")
        if settings["meth_scan_min_cells"] < 1:
            raise ValueError("meth_scan_min_cells must be >= 1")
        if settings["meth_scan_bridge_gaps"] < 0:
            raise ValueError("meth_scan_bridge_gaps must be >= 0")
        if settings["meth_matrix_cores"] < 1:
            raise ValueError("meth_matrix_cores must be >= 1")
    settings["_stage_sequence"] = build_stage_sequence(settings)
    if stage in ("count_mapped_reads", "estimated_cells") and settings["_barcode_mode"] == "gexcb":
        raise ValueError(
            f"stage {stage} is not used when gexcb is set; use split_bams with --gexcb"
        )
    if (
        settings.get("runner") == "slurm"
        and any(key in slurm_cfg_raw for key in SLURM_NEST_STAGE_KEYS)
        and stage != "all"
        and not settings["slurm_partition"]
    ):
        raise ValueError(
            f"workflow slurm config missing partition for stage '{stage}'; "
            "add slurm.<stage>.partition or pass --slurm-partition"
        )
    settings["slurm_partition"] = settings["slurm_partition"] or "cpu"
    settings["slurm_mem"] = settings["slurm_mem"] or "16G"
    settings["slurm_cpus_per_task"] = settings["slurm_cpus_per_task"] or 8
    settings["slurm_output"] = settings["slurm_output"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.out"
    )
    settings["slurm_error"] = settings["slurm_error"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.err"
    )

    if stage == "all":
        for stage_name in settings["_stage_sequence"]:
            validate_required_for_stage(stage_name, settings)
    else:
        validate_required_for_stage(stage, settings)
    return settings


def generate_local_script(command: str, output_path: Path) -> None:
    content = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"{command}\n"
    )
    write_text(output_path, content)
    output_path.chmod(0o755)


def generate_slurm_script(
    command: str, output_path: Path, log_dir: Path, args: argparse.Namespace
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={args.job_name}",
        f"#SBATCH --partition={args.slurm_partition}",
        f"#SBATCH --cpus-per-task={args.slurm_cpus_per_task}",
        f"#SBATCH --mem={args.slurm_mem}",
        f"#SBATCH --output={args.slurm_output}",
        f"#SBATCH --error={args.slurm_error}",
        "",
        "set -euo pipefail",
        "",
        command,
        "",
    ]
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def submit_script(path: Path, runner: str) -> None:
    if runner == "local":
        subprocess.run(["bash", str(path)], check=True)
    else:
        subprocess.run(["sbatch", str(path)], check=True)


def parse_generated_paths(command_output: str) -> list[Path]:
    generated: list[Path] = []
    for line in command_output.splitlines():
        prefix = "[make_cmd] generated="
        if line.startswith(prefix):
            generated.append(Path(line[len(prefix) :].strip()))
    return generated


def build_stage_passthrough_args(argv: list[str]) -> list[str]:
    passthrough: list[str] = []
    flags_without_value = {
        "--submit",
        "--dry-run",
        "--skip-workdir-input-checks",
    }
    index = 0
    while index < len(argv):
        token = argv[index]
        if token in {"--stage", "--runner"}:
            index += 2
            continue
        if token.startswith("--stage=") or token.startswith("--runner="):
            index += 1
            continue
        if token in {"--submit", "--dry-run", "--skip-workdir-input-checks"}:
            index += 1
            continue
        if token in flags_without_value:
            passthrough.append(token)
            index += 1
            continue
        if token.startswith("--"):
            passthrough.append(token)
            if index + 1 < len(argv):
                passthrough.append(argv[index + 1])
            index += 2
            continue
        passthrough.append(token)
        index += 1
    return passthrough


def driver_scripts_for_stage(
    stage_name: str,
    scripts: list[Path],
    *,
    runner: str,
    stage_sequence: list[str],
    settings: dict | None = None,
) -> list[Path]:
    prefix = stage_prefix_map(stage_sequence)[stage_name]
    if stage_name == "demux_extract_bc":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_demux_extract_bc.sh"
            ]
        return [script for script in scripts if script.suffix == ".sbatch"]
    if stage_name == "regroup_shards":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_regroup_shards.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_regroup_shards_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "bismark_align":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_bismark_align.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_bismark_align_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "bam_sort":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_bam_sort.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_bam_sort_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "count_mapped_reads":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_count_mapped_reads.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_count_mapped_reads_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "estimated_cells":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_estimated_cells")
        ]
    if stage_name == "split_bams":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_split_bams.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_split_bams_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "merge_fr_bams":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_merge_fr_bams.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_merge_fr_bams_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "bam_to_allc":
        if runner == "local":
            return [
                script
                for script in scripts
                if script.name == f"{prefix}_bam_to_allc.sh"
            ]
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_bam_to_allc_")
            and script.suffix == ".sbatch"
        ]
    if stage_name == "saturation":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_saturation")
        ]
    if stage_name == "qc_summary":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_qc_summary")
        ]
    if stage_name == "allc_to_matrix":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_allc_to_matrix")
        ]
    if stage_name == "meth_smooth":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_meth_smooth")
        ]
    if stage_name == "meth_scan":
        return [
            script
            for script in scripts
            if script.name.startswith(f"{prefix}_meth_scan")
        ]
    return scripts


def generate_local_driver_script(
    stage_scripts: list[tuple[str, list[Path]]],
    output_path: Path,
    stage_sequence: list[str],
) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        "",
    ]
    for stage_name, scripts in stage_scripts:
        runnable = driver_scripts_for_stage(
            stage_name, scripts, runner="local", stage_sequence=stage_sequence
        )
        if not runnable:
            continue
        for script_path in runnable:
            lines.append(f'bash "$SCRIPT_DIR/{script_path.name}"')
    lines.append("")
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def generate_slurm_driver_script(
    stage_scripts: list[tuple[str, list[Path]]],
    output_path: Path,
    log_dir: Path,
    settings: dict,
) -> None:
    stage_sequence = settings["_stage_sequence"]
    lines = [
        "submit_with_dep() {",
        '  local script_path="$1"',
        '  local dep_chain="$2"',
        "  local out",
        '  if [[ -n "$dep_chain" ]]; then',
        '    out="$(sbatch --dependency=afterok:${dep_chain} "$script_path")"',
        "  else",
        '    out="$(sbatch "$script_path")"',
        "  fi",
        '  echo "$out" >&2',
        '  echo "${out##* }"',
        "}",
        "",
        "join_deps() {",
        "  local joined=''",
        '  for item in "$@"; do',
        '    if [[ -z "$item" ]]; then',
        "      continue",
        "    fi",
        '    if [[ -z "$joined" ]]; then',
        '      joined="$item"',
        "    else",
        '      joined="${joined}:$item"',
        "    fi",
        "  done",
        '  echo "$joined"',
        "}",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'prev_stage_deps=""',
        "",
    ]
    scripts_by_stage = dict(stage_scripts)
    for stage_name, scripts in stage_scripts:
        runnable = driver_scripts_for_stage(
            stage_name, scripts, runner="slurm", stage_sequence=stage_sequence
        )
        if not runnable:
            continue
        lines.append(f'echo "[run.sbatch] stage={stage_name}"')
        if stage_name == "demux_extract_bc":
            demux_prefix = stage_prefix_map(stage_sequence)["demux_extract_bc"]
            chunk_scripts = [
                script
                for script in runnable
                if script.name.startswith(f"{demux_prefix}_demux_extract_bc_")
            ]
            aggregate_scripts = [
                script
                for script in runnable
                if script.name.startswith(f"{demux_prefix}_aggregate_")
            ]
            chunk_job_vars: list[str] = []
            for index, script_path in enumerate(chunk_scripts):
                var_name = f"jid_demux_chunk_{index}"
                lines.append(
                    f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$prev_stage_deps")"'
                )
                chunk_job_vars.append(var_name)
            if chunk_job_vars:
                deps_join = " ".join(f"${var_name}" for var_name in chunk_job_vars)
                lines.append(f'chunk_deps="$(join_deps {deps_join})"')
                post_demux_job_vars: list[str] = []
                for script_path in aggregate_scripts:
                    var_name = f"jid_demux_agg_{len(post_demux_job_vars)}"
                    lines.append(
                        f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$chunk_deps")"'
                    )
                    post_demux_job_vars.append(var_name)
                regroup_scripts = driver_scripts_for_stage(
                    "regroup_shards",
                    scripts_by_stage.get("regroup_shards", []),
                    runner="slurm",
                    stage_sequence=stage_sequence,
                )
                for index, script_path in enumerate(regroup_scripts):
                    var_name = f"jid_regroup_{index}"
                    lines.append(
                        f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$chunk_deps")"'
                    )
                    post_demux_job_vars.append(var_name)
                deps_join = " ".join(f"${var_name}" for var_name in post_demux_job_vars)
                lines.append(f'prev_stage_deps="$(join_deps {deps_join})"')
            continue
        if stage_name == "regroup_shards":
            continue
        job_vars: list[str] = []
        for index, script_path in enumerate(runnable):
            var_name = f"jid_{stage_name}_{index}".replace("-", "_")
            lines.append(
                f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$prev_stage_deps")"'
            )
            job_vars.append(var_name)
        deps_join = " ".join(f"${var_name}" for var_name in job_vars)
        lines.append(f'prev_stage_deps="$(join_deps {deps_join})"')
        lines.append("")
    lines.append('echo "[run.sbatch] done final_dep=${prev_stage_deps}"')
    driver_command = "\n".join(lines)
    slurm_args = argparse.Namespace(
        job_name=f"seeksoul_all_driver_{settings['sample_id']}",
        slurm_partition=settings["slurm_partition"],
        slurm_mem=settings["slurm_mem"],
        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
        slurm_output=settings["slurm_output"].replace(
            "%x", f"seeksoul_all_driver_{settings['sample_id']}"
        ),
        slurm_error=settings["slurm_error"].replace(
            "%x", f"seeksoul_all_driver_{settings['sample_id']}"
        ),
    )
    generate_slurm_script(driver_command, output_path, log_dir, slurm_args)


def apply_stage_slurm_settings(settings: dict, stage: str) -> dict:
    slurm_cfg_raw = settings.get("_slurm_cfg_raw", {})
    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, stage)
    updated = dict(settings)
    updated["slurm_partition"] = stage_slurm_cfg.get("partition") or settings["slurm_partition"]
    updated["slurm_mem"] = stage_slurm_cfg.get("mem") or settings["slurm_mem"]
    updated["slurm_cpus_per_task"] = (
        stage_slurm_cfg.get("cpus_per_task") or settings["slurm_cpus_per_task"]
    )
    updated["slurm_output"] = stage_slurm_cfg.get("output") or settings["slurm_output"]
    updated["slurm_error"] = stage_slurm_cfg.get("error") or settings["slurm_error"]
    return updated


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"

    if settings["stage"] != "all":
        validate_inputs_for_stage(
            settings["stage"],
            settings,
            sample_work,
            skip_workdir_inputs=bool(args.skip_workdir_input_checks),
        )

    if settings["stage"] == "all":
        stage_sequence = settings["_stage_sequence"]
        for stage_name in stage_sequence:
            validate_inputs_for_stage(
                stage_name,
                settings,
                sample_work,
                skip_workdir_inputs=True,
            )
        passthrough_args = build_stage_passthrough_args(sys.argv[1:])
        stage_scripts: list[tuple[str, list[Path]]] = []
        for stage_name in stage_sequence:
            stage_argv = [
                sys.executable,
                __file__,
                *passthrough_args,
                "--runner",
                settings["runner"],
                "--stage",
                stage_name,
            ]
            if settings["dry_run"]:
                stage_argv.append("--dry-run")
            stage_argv.append("--skip-workdir-input-checks")
            completed = subprocess.run(
                stage_argv, check=False, capture_output=True, text=True
            )
            if completed.stdout:
                print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
            if completed.returncode != 0:
                return completed.returncode
            stage_scripts.append((stage_name, parse_generated_paths(completed.stdout)))

        driver_path: Path
        if settings["runner"] == "local":
            driver_path = command_dir / "run.sh"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                generate_local_driver_script(stage_scripts, driver_path, stage_sequence)
        else:
            driver_path = command_dir / "run.sbatch"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                driver_settings = apply_stage_slurm_settings(
                    settings, stage_sequence[-1]
                )
                generate_slurm_driver_script(
                    stage_scripts, driver_path, log_dir, driver_settings
                )

        if not settings["dry_run"] and driver_path.exists():
            print(f"[make_cmd] generated={driver_path}")
        if settings["submit"] and not settings["dry_run"]:
            if settings["runner"] == "slurm":
                subprocess.run(["bash", str(driver_path)], check=True)
                print("[make_cmd] submitted_driver=1")
                print("[make_cmd] submit_mode=client_side_sbatch_dag")
            else:
                submit_script(driver_path, settings["runner"])
                print("[make_cmd] submitted_driver=1")

        print("[make_cmd] stage=all helper generation complete")
        return 0

    generated_scripts: list[Path] = []
    if settings["stage"] == "fastp_split":
        base_name = f"{stage_prefix_map(settings['_stage_sequence'])['fastp_split']}_fastp_split"
        command_args = argparse.Namespace(
            r1=settings["r1"],
            r2=settings["r2"],
            fastp_threads=settings["fastp_threads"],
            number_of_split_parts=settings["number_of_split_parts"],
            fastp_bin=settings["fastp_bin"],
        )
        command = build_fastp_split_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / f"{base_name}.sh"
        else:
            script_path = command_dir / f"{base_name}.sbatch"

        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")

        if settings["dry_run"]:
            return 0

        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_fastp_split_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"],
                slurm_error=settings["slurm_error"],
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "demux_extract_bc":
        command_args = argparse.Namespace(
            barcode_whitelist=settings["barcode_whitelist"],
            barcode_hamming_distance=settings["barcode_hamming_distance"],
            gzip_level=settings["gzip_level"],
            filter_ch=settings["filter_ch"],
            split_fastq_prefix_bases=settings["split_fastq_prefix_bases"],
        )
        demux_dir = sample_work / "demux"
        demux_prefix = stage_prefix_map(settings["_stage_sequence"])["demux_extract_bc"]
        if settings["runner"] == "local":
            script_path = command_dir / f"{demux_prefix}_demux_extract_bc.sh"
            command = build_demux_local_batch_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            if settings.get("number_of_split_parts") is None:
                raise ValueError(
                    "number_of_split_parts is required for slurm demux script generation"
                )
            chunks = build_demux_chunks_from_config(
                sample_work, settings["number_of_split_parts"]
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, r1_path, r2_path, out_prefix in chunks:
                base_name = f"{demux_prefix}_demux_extract_bc_{chunk_id}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_demux_chunk_command(
                    command_args, r1_path, r2_path, out_prefix
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"seeksoul_demux_{settings['sample_id']}_{chunk_id}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"seeksoul_demux_{settings['sample_id']}_{chunk_id}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_demux_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)

            aggregate_script = command_dir / f"{demux_prefix}_aggregate_ct_qc.sbatch"
            aggregate_command = build_aggregate_ct_command(demux_dir)
            aggregate_output = settings["slurm_output"].replace(
                "%x", f"seeksoul_aggregate_ct_{settings['sample_id']}"
            )
            aggregate_error = settings["slurm_error"].replace(
                "%x", f"seeksoul_aggregate_ct_{settings['sample_id']}"
            )
            print(f"[make_cmd] script={aggregate_script}")
            print(f"[make_cmd] command={aggregate_command}")
            if not settings["dry_run"]:
                slurm_args = argparse.Namespace(
                    job_name=f"seeksoul_aggregate_ct_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=aggregate_output,
                    slurm_error=aggregate_error,
                )
                generate_slurm_script(
                    aggregate_command, aggregate_script, log_dir, slurm_args
                )
                submit_script_path = command_dir / f"{demux_prefix}_demux_extract_bc_submit.sh"
                generate_local_script(
                    build_demux_slurm_submit_command(sample_work),
                    submit_script_path,
                )
                generated_scripts.append(aggregate_script)
                generated_scripts.append(submit_script_path)
            else:
                generated_scripts.append(aggregate_script)
    elif settings["stage"] == "regroup_shards":
        regroup_prefix = stage_prefix_map(settings["_stage_sequence"])["regroup_shards"]
        if settings["runner"] == "local":
            script_path = command_dir / f"{regroup_prefix}_regroup_shards.sh"
            command = build_regroup_work_command(sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            prefixes = resolve_prefix_chunks(
                discover=lambda: list(wic.discover_demux_subshards(sample_work / "demux").keys()),
                plan_by_prefix=lambda _demux_dir, prefix_list: prefix_list,
                base_dir=sample_work / "demux",
                settings=settings,
                label="regroup prefixes",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] prefix_count={len(prefixes)}")
            for prefix in prefixes:
                base_name = f"{regroup_prefix}_regroup_shards_{prefix}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_regroup_prefix_command(sample_work, prefix)
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"seeksoul_regroup_{settings['sample_id']}_{prefix}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"seeksoul_regroup_{settings['sample_id']}_{prefix}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_regroup_{settings['sample_id']}_{prefix}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "bismark_align":
        command_args = argparse.Namespace(
            bismark_ref=settings["bismark_ref"],
            bismark_parallel=settings["bismark_parallel"],
            bismark_max_insert=settings["bismark_max_insert"],
            bismark_bin=settings["bismark_bin"],
        )
        align_prefix = stage_prefix_map(settings["_stage_sequence"])["bismark_align"]
        if settings["runner"] == "local":
            script_path = command_dir / f"{align_prefix}_bismark_align.sh"
            command = build_bismark_align_work_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = resolve_prefix_chunks(
                discover=lambda: discover_bismark_align_chunks(sample_work),
                plan_by_prefix=lambda demux_dir, prefixes: wic.plan_demux_align_chunks_by_prefix(
                    demux_dir, prefixes
                ),
                base_dir=sample_work / "demux",
                settings=settings,
                label="demux align inputs",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2 in chunks:
                base_name = f"{align_prefix}_bismark_align_{chunk_id}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_bismark_align_chunk_command(
                    command_args,
                    sample_work,
                    chunk_id,
                    fwd_r1,
                    fwd_r2,
                    rev_r1,
                    rev_r2,
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"seeksoul_bismark_{settings['sample_id']}_{chunk_id}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"seeksoul_bismark_{settings['sample_id']}_{chunk_id}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_bismark_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "bam_sort":
        command_args = argparse.Namespace(
            sort_threads=settings["sort_threads"],
            samtools_bin=settings["samtools_bin"],
        )
        sort_prefix = stage_prefix_map(settings["_stage_sequence"])["bam_sort"]
        if settings["runner"] == "local":
            script_path = command_dir / f"{sort_prefix}_bam_sort.sh"
            command = build_bam_sort_work_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = resolve_prefix_chunks(
                discover=lambda: discover_bam_sort_chunks(sample_work),
                plan_by_prefix=lambda align_dir, prefixes: wic.plan_bismark_pe_bams_by_prefix(
                    align_dir, prefixes
                ),
                base_dir=sample_work / "align",
                settings=settings,
                label="Bismark PE BAMs",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, forward_bam, reverse_bam in chunks:
                base_name = f"{sort_prefix}_bam_sort_{chunk_id}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_bam_sort_chunk_command(
                    command_args,
                    sample_work,
                    chunk_id,
                    forward_bam,
                    reverse_bam,
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"seeksoul_bamsort_{settings['sample_id']}_{chunk_id}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"seeksoul_bamsort_{settings['sample_id']}_{chunk_id}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_bamsort_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "count_mapped_reads":
        script_name = stage_script_name(settings, "count_mapped_reads")
        if settings["runner"] == "local":
            script_path = command_dir / script_name
            command = build_count_mapped_reads_work_command(sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = resolve_prefix_chunks(
                discover=lambda: discover_bam_sort_chunks(sample_work),
                plan_by_prefix=lambda align_dir, prefixes: wic.plan_bismark_pe_bams_by_prefix(
                    align_dir, prefixes
                ),
                base_dir=sample_work / "align",
                settings=settings,
                label="unsorted Bismark PE BAMs",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, _forward_bam, _reverse_bam in chunks:
                base_name = stage_script_name(
                    settings, "count_mapped_reads", suffix="sbatch", chunk_id=chunk_id
                ).removesuffix(".sbatch")
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_count_mapped_reads_chunk_command(sample_work, chunk_id)
                chunk_output = settings["slurm_output"].replace(
                    "%x",
                    f"seeksoul_count_{settings['sample_id']}_{chunk_id}",
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x",
                    f"seeksoul_count_{settings['sample_id']}_{chunk_id}",
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_count_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "estimated_cells":
        script_name = stage_script_name(settings, "estimated_cells")
        command = build_estimated_cells_command(
            sample_work,
            settings["expected_cell_num"],
            settings.get("force_cell_num"),
        )
        if settings["runner"] == "local":
            script_path = command_dir / script_name
        else:
            script_path = command_dir / script_name.replace(".sh", ".sbatch")
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_estcells_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_estcells_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_estcells_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "split_bams":
        command_args = argparse.Namespace(
            gexcb=settings.get("gexcb"),
            split_bams_cores=settings["split_bams_cores"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "split_bams")
            command = build_split_bams_work_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = resolve_prefix_chunks(
                discover=lambda: wic.discover_bismark_sortbyname_bams(sample_work / "align"),
                plan_by_prefix=lambda align_dir, prefixes: wic.plan_bismark_sortbyname_bams_by_prefix(
                    align_dir, prefixes
                ),
                base_dir=sample_work / "align",
                settings=settings,
                label="sortbyname Bismark PE BAMs",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, _forward_bam, _reverse_bam in chunks:
                base_name = stage_script_name(
                    settings, "split_bams", suffix="sbatch", chunk_id=chunk_id
                ).removesuffix(".sbatch")
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_split_bams_chunk_command(
                    command_args, sample_work, chunk_id
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x",
                    f"seeksoul_split_{settings['sample_id']}_{chunk_id}",
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x",
                    f"seeksoul_split_{settings['sample_id']}_{chunk_id}",
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_split_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "merge_fr_bams":
        command_args = argparse.Namespace(
            merge_fr_bams_cores=settings["merge_fr_bams_cores"],
            samtools_bin=settings["samtools_bin"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "merge_fr_bams")
            command = build_merge_fr_bams_work_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            pairs = resolve_prefix_chunks(
                discover=lambda: wic.discover_split_bam_chunk_pairs(sample_work / "split_bams"),
                plan_by_prefix=lambda split_root, prefixes: wic.plan_split_bam_chunk_pairs_by_prefix(
                    split_root, prefixes
                ),
                base_dir=sample_work / "split_bams",
                settings=settings,
                label="split BAM chunk pairs",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(pairs)}")
            for chunk_id, _forward_dir, _reverse_dir in pairs:
                base_name = stage_script_name(
                    settings, "merge_fr_bams", suffix="sbatch", chunk_id=chunk_id
                ).removesuffix(".sbatch")
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_merge_fr_bams_chunk_command(
                    command_args, sample_work, chunk_id
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x",
                    f"seeksoul_merge_{settings['sample_id']}_{chunk_id}",
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x",
                    f"seeksoul_merge_{settings['sample_id']}_{chunk_id}",
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_merge_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "bam_to_allc":
        command_args = argparse.Namespace(
            genome_fa=settings["genome_fa"],
            chrom_size_path=settings["chrom_size_path"],
            bam_to_allc_cores=settings["bam_to_allc_cores"],
            allcools_tag=settings["allcools_tag"],
            samtools_bin=settings["samtools_bin"],
            allcools_bin=settings["allcools_bin"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "bam_to_allc")
            command = build_bam_to_allc_work_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = resolve_prefix_chunks(
                discover=lambda: wic.discover_merged_fr_bam_chunks(
                    sample_work / "split_bams" / "merged"
                ),
                plan_by_prefix=lambda merged_root, prefixes: wic.plan_merged_fr_bam_chunks_by_prefix(
                    merged_root, prefixes
                ),
                base_dir=sample_work / "split_bams" / "merged",
                settings=settings,
                label="merged FR BAM chunks",
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk_id, _bam_dir, _filtered_barcode in chunks:
                base_name = stage_script_name(
                    settings, "bam_to_allc", suffix="sbatch", chunk_id=chunk_id
                ).removesuffix(".sbatch")
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_bam_to_allc_chunk_command(
                    command_args, sample_work, chunk_id
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x",
                    f"seeksoul_allc_{settings['sample_id']}_{chunk_id}",
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x",
                    f"seeksoul_allc_{settings['sample_id']}_{chunk_id}",
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"seeksoul_allc_{settings['sample_id']}_{chunk_id}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "saturation":
        command_args = argparse.Namespace(
            saturation_script=settings["saturation_script"],
            chrom_size_path=str(wic.resolve_config_path(settings["chrom_size_path"])),
            saturation_reads_threshold=settings["saturation_reads_threshold"],
            saturation_max_cells=settings["saturation_max_cells"],
            saturation_sample_seed=settings["saturation_sample_seed"],
            saturation_linear_r2_threshold=settings[
                "saturation_linear_r2_threshold"
            ],
        )
        command = build_saturation_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "saturation")
        else:
            script_path = command_dir / stage_script_name(
                settings, "saturation", suffix="sbatch"
            )
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_saturation_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_saturation_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_saturation_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "qc_summary":
        command_args = argparse.Namespace(
            qc_summary_script=settings["qc_summary_script"],
            cbcsv=settings.get("cbcsv"),
            mito_chromosomes=settings.get("mito_chromosomes"),
        )
        command = build_qc_summary_command(
            command_args,
            sample_work,
            sample_id=settings["sample_id"],
            barcode_mode=settings["_barcode_mode"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "qc_summary")
        else:
            script_path = command_dir / stage_script_name(
                settings, "qc_summary", suffix="sbatch"
            )
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_qc_summary_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_qc_summary_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_qc_summary_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "allc_to_matrix":
        command_args = argparse.Namespace(
            allc_to_matrix_script=settings["allc_to_matrix_script"],
            meth_context=settings["meth_context"],
            meth_chunksize=settings["meth_chunksize"],
            meth_round_sites=settings["meth_round_sites"],
            meth_main_chroms_only=settings["meth_main_chroms_only"],
            meth_exclude_contigs=settings["meth_exclude_contigs"],
        )
        command = build_allc_to_matrix_command(
            command_args,
            sample_work,
            barcode_mode=settings["_barcode_mode"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "allc_to_matrix")
        else:
            script_path = command_dir / stage_script_name(
                settings, "allc_to_matrix", suffix="sbatch"
            )
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_allc_to_matrix_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_allc_to_matrix_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_allc_to_matrix_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "meth_smooth":
        command_args = argparse.Namespace(
            meth_smooth_script=settings["meth_smooth_script"],
            meth_smooth_bandwidth=settings["meth_smooth_bandwidth"],
            meth_smooth_use_weights=settings["meth_smooth_use_weights"],
        )
        command = build_meth_smooth_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "meth_smooth")
        else:
            script_path = command_dir / stage_script_name(
                settings, "meth_smooth", suffix="sbatch"
            )
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_meth_smooth_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_meth_smooth_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_meth_smooth_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "meth_scan":
        command_args = argparse.Namespace(
            meth_scan_script=settings["meth_scan_script"],
            meth_scan_bandwidth=settings["meth_scan_bandwidth"],
            meth_scan_stepsize=settings["meth_scan_stepsize"],
            meth_scan_var_threshold=settings["meth_scan_var_threshold"],
            meth_scan_min_cells=settings["meth_scan_min_cells"],
            meth_scan_bridge_gaps=settings["meth_scan_bridge_gaps"],
            meth_matrix_cores=settings["meth_matrix_cores"],
        )
        command = build_meth_scan_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / stage_script_name(settings, "meth_scan")
        else:
            script_path = command_dir / stage_script_name(
                settings, "meth_scan", suffix="sbatch"
            )
        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")
        if settings["dry_run"]:
            return 0
        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"seeksoul_meth_scan_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"].replace(
                    "%x", f"seeksoul_meth_scan_{settings['sample_id']}"
                ),
                slurm_error=settings["slurm_error"].replace(
                    "%x", f"seeksoul_meth_scan_{settings['sample_id']}"
                ),
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    else:
        raise ValueError(f"unsupported stage: {settings['stage']}")

    for script_path in generated_scripts:
        print(f"[make_cmd] generated={script_path}")

    if settings["submit"]:
        for script_path in generated_scripts:
            submit_script(script_path, settings["runner"])
        print(f"[make_cmd] submitted_count={len(generated_scripts)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
