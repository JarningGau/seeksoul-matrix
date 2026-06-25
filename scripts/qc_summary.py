#!/usr/bin/env python3
"""Gather per-cell and sample-level QC metrics into summary tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from _version import __version__

CPG_CONTEXTS = ("CGT", "CGG", "CGC", "CGA")
NA = "NA"

AGGREGATED_MC_RATE_COLUMNS = (
    "CG_mc_rate",
    "CH_mc_rate",
    "CHG_mc_rate",
    "CHH_mc_rate",
    "CA_mc_rate",
    "CT_mc_rate",
    "CC_mc_rate",
)

CORE_CELL_COLUMNS = (
    "cell_barcode",
    "aligned_reads",
    "CtoT",
    "total_cpg_number",
    "genome_cov",
    "genome_cov_raw_umi",
    "genome_cov_new_umi",
    "cell_saturation",
    *AGGREGATED_MC_RATE_COLUMNS,
)

SAMPLE_SUMMARY_COLUMNS = (
    "sample_id",
    "raw_reads",
    "total_bases",
    "duplication_rate",
    "valid_barcode_rate",
    "valid_demux_rate",
    "barcode_corrected_fraction",
    "dropped_too_short",
    "dropped_chimeric",
    "forward_reads",
    "reverse_reads",
    "rate_17lme",
    "CtoT",
    "rate_7f",
    "rate_7f17lme",
    "cc_mean",
    "mapped_to_genome",
    "confidently_mapped",
    "cpg_methylation_rate",
    "chg_methylation_rate",
    "chh_methylation_rate",
    "estimated_cells",
    "reads_in_cells",
    "fraction_reads_in_cells",
    "observed_median_genome_fraction",
    "theoretical_max_median_genome_fraction",
    "saturation_rate",
    "extrapolation_model",
    "sampled_cell_count",
    "sample_seed",
    "median_genome_cov",
    "median_total_cpg_number",
    "median_aligned_reads",
    "median_cell_saturation",
)

WGS_HEADER = (
    "Samplename",
    "Estimated_Number_of_Cells",
    "Number_of_Reads",
    "Valid_Barcode_Ratio",
    "Dropped_Too_Short",
    "Dropped_Chimeric",
    "Valid_7F_Reads_Rate",
    "Valid_17LME_Reads_Rate",
    "Valid_7F17LME_Reads_Rate",
    "C-T_Conversion",
    "C-C_Ratio",
    "Reads_Mapped_to_Genome",
    "Reads_Mapped_Confidently_to_Genome",
    "CpG_Methylation_Rate",
    "CHG_Methylation_Rate",
    "CHH_Methylation_Rate",
    "Unknown_Methylation_Rate",
    "CpG_Coverage_Rate",
    "Total_CPGs_Detected",
    "Genome_Coverage_Rate_of_Max_Cell",
    "CPGs_of_Max_Cell",
    "Reads_of_Max_Cell",
    "Saturation_of_Max_Cell",
    "Genome_Coverage_Rate_of_Median_Cell",
    "CPGs_of_Median_Cell",
    "Reads_of_Median_Cell",
    "Saturation_of_Median_Cell",
    "Fraction_Reads_in_Cells",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Gather per-cell and sample-level QC metrics under "
            "<work_path>/summary/."
        )
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing prior stage outputs.",
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Sample identifier written to summary tables.",
    )
    parser.add_argument(
        "--output-dir",
        default="summary",
        help="Output directory relative to work-path. Default: summary.",
    )
    parser.add_argument(
        "--barcode-mode",
        choices=("methylation_only", "gexcb"),
        default="methylation_only",
        help="Cell read-count source. Default: methylation_only.",
    )
    parser.add_argument(
        "--cbcsv",
        default=None,
        help="Optional gexcb map CSV with columns m_cb,gex_cb.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


def nested_get(payload: dict, *keys, default=None):
    current = payload
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def format_rate(value: float | None) -> str:
    if value is None:
        return NA
    return f"{value:.6f}"


def format_pct(value: float | None) -> str:
    if value is None:
        return NA
    return f"{value:.2%}"


def format_pct_one_decimal(value: float | None) -> str:
    if value is None:
        return NA
    return f"{value:.1f}%"


def parse_fastp_metrics(fastp_path: Path) -> dict[str, str | int | float]:
    with fastp_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    summary = payload.get("summary") or {}
    after = summary.get("after_filtering") or summary.get("before_filtering") or {}
    before = summary.get("before_filtering") or after
    duplication = payload.get("duplication") or {}
    return {
        "raw_reads": int(after.get("total_reads") or before.get("total_reads") or 0),
        "total_bases": int(after.get("total_bases") or before.get("total_bases") or 0),
        "duplication_rate": float(duplication.get("rate") or 0.0),
    }


def aggregate_demux_stats(demux_dir: Path) -> dict[str, str | int | float]:
    stats_files = sorted(demux_dir.glob("*.stats.json"))
    totals = {
        "total": 0,
        "barcode_passed": 0,
        "corrected": 0,
        "too_short": 0,
        "chimeric": 0,
        "valid": 0,
        "forward": 0,
        "reverse": 0,
        "num_17lme": 0,
        "ct_convertible": 0,
        "ct_converted": 0,
    }
    for stats_path in stats_files:
        with stats_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        funnel = payload.get("funnel") or {}
        passed = funnel.get("barcode_passed") or {}
        valid = passed.get("valid") or {}
        chimeric = passed.get("chimeric_filtered") or {}
        ct = payload.get("ct") or {}
        totals["total"] += int(funnel.get("total") or 0)
        totals["barcode_passed"] += int(passed.get("total") or 0)
        totals["corrected"] += int(passed.get("corrected") or 0)
        totals["too_short"] += int(passed.get("too_short") or 0)
        totals["chimeric"] += int(chimeric.get("total") or 0)
        totals["valid"] += int(valid.get("total") or 0)
        totals["forward"] += int(valid.get("forward") or 0)
        totals["reverse"] += int(valid.get("reverse") or 0)
        totals["num_17lme"] += int(ct.get("num_17lme") or 0)
        totals["ct_convertible"] += int(ct.get("ct_convertible_bases") or 0)
        totals["ct_converted"] += int(ct.get("ct_converted_bases") or 0)

    valid_barcode_rate = (
        totals["barcode_passed"] / totals["total"] if totals["total"] else None
    )
    valid_demux_rate = totals["valid"] / totals["total"] if totals["total"] else None
    corrected_fraction = (
        totals["corrected"] / totals["barcode_passed"]
        if totals["barcode_passed"]
        else None
    )
    dropped_too_short = (
        totals["too_short"] / totals["total"] if totals["total"] else None
    )
    denom_chimeric = totals["barcode_passed"] - totals["too_short"]
    dropped_chimeric = (
        totals["chimeric"] / denom_chimeric if denom_chimeric > 0 else None
    )
    rate_17lme = totals["num_17lme"] / totals["total"] if totals["total"] else None
    ctot = (
        totals["ct_converted"] / totals["ct_convertible"]
        if totals["ct_convertible"]
        else None
    )
    return {
        "valid_barcode_rate": format_rate(valid_barcode_rate),
        "valid_demux_rate": format_rate(valid_demux_rate),
        "barcode_corrected_fraction": format_rate(corrected_fraction),
        "dropped_too_short": format_rate(dropped_too_short),
        "dropped_chimeric": format_rate(dropped_chimeric),
        "forward_reads": totals["forward"],
        "reverse_reads": totals["reverse"],
        "rate_17lme": format_rate(rate_17lme),
        "CtoT": format_rate(ctot),
        "rate_7f": NA,
        "rate_7f17lme": NA,
        "cc_mean": NA,
    }


def load_ctot_map(ctot_path: Path) -> dict[str, str]:
    if not ctot_path.is_file():
        return {}
    mapping: dict[str, str] = {}
    with ctot_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            barcode = (row.get("CR") or "").strip()
            if not barcode:
                continue
            mapping[barcode] = (row.get("CtoT") or "").strip()
    return mapping


def extract_percentage(text: str, pattern: str) -> float:
    match = re.search(pattern, text)
    return float(match.group(1)) if match else 0.0


def extract_number(text: str, pattern: str) -> int:
    match = re.search(pattern, text)
    return int(match.group(1)) if match else 0


def parse_bismark_report(report_path: Path) -> dict[str, int | float] | None:
    try:
        content = report_path.read_text(encoding="utf-8")
    except OSError:
        return None
    total_pairs = extract_number(content, r"Sequence pairs analysed in total:\s*(\d+)")
    unique_pairs = extract_number(
        content,
        r"Number of paired-end alignments with a unique best hit:\s*(\d+)",
    )
    not_unique_pairs = extract_number(
        content,
        r"Sequence pairs did not map uniquely:\s*(\d+)",
    )
    total_reads = total_pairs * 2
    unique_reads = unique_pairs * 2
    aligned_reads = unique_reads + not_unique_pairs * 2
    return {
        "totalreads": total_reads,
        "alignedreads": aligned_reads,
        "uniquereads": unique_reads,
        "methylated_cpg": extract_number(
            content, r"Total methylated C's in CpG context:\s*(\d+)"
        ),
        "unmethylated_cpg": extract_number(
            content, r"Total unmethylated C's in CpG context:\s*(\d+)"
        ),
        "methylated_chg": extract_number(
            content, r"Total methylated C's in CHG context:\s*(\d+)"
        ),
        "unmethylated_chg": extract_number(
            content, r"Total unmethylated C's in CHG context:\s*(\d+)"
        ),
        "methylated_chh": extract_number(
            content, r"Total methylated C's in CHH context:\s*(\d+)"
        ),
        "unmethylated_chh": extract_number(
            content, r"Total unmethylated C's in CHH context:\s*(\d+)"
        ),
        "cpg_methylation_rate": extract_percentage(
            content, r"C methylated in CpG context:\s*([\d.]+)%"
        ),
        "chg_methylation_rate": extract_percentage(
            content, r"C methylated in CHG context:\s*([\d.]+)%"
        ),
        "chh_methylation_rate": extract_percentage(
            content, r"C methylated in CHH context:\s*([\d.]+)%"
        ),
    }


def context_rate(methylated: int, unmethylated: int) -> float:
    total = methylated + unmethylated
    if total <= 0:
        return 0.0
    return methylated / total * 100.0


def aggregate_bismark_reports(align_dir: Path) -> dict[str, str | float | int]:
    report_paths = sorted(align_dir.glob("*_bismark_bt2_PE_report.txt"))
    pooled = {
        "totalreads": 0,
        "alignedreads": 0,
        "uniquereads": 0,
        "methylated_cpg": 0,
        "unmethylated_cpg": 0,
        "methylated_chg": 0,
        "unmethylated_chg": 0,
        "methylated_chh": 0,
        "unmethylated_chh": 0,
    }
    parsed_any = False
    for report_path in report_paths:
        metrics = parse_bismark_report(report_path)
        if metrics is None:
            continue
        parsed_any = True
        for key in pooled:
            pooled[key] += int(metrics[key])

    if not parsed_any:
        return {
            "mapped_to_genome": NA,
            "confidently_mapped": NA,
            "cpg_methylation_rate": NA,
            "chg_methylation_rate": NA,
            "chh_methylation_rate": NA,
            "uniquereads": 0,
        }

    totalreads = pooled["totalreads"]
    alignedreads = pooled["alignedreads"]
    uniquereads = pooled["uniquereads"]
    mapped = alignedreads / totalreads if totalreads else None
    confident = uniquereads / totalreads if totalreads else None
    return {
        "mapped_to_genome": format_rate(mapped),
        "confidently_mapped": format_rate(confident),
        "cpg_methylation_rate": format_rate(
            context_rate(pooled["methylated_cpg"], pooled["unmethylated_cpg"]) / 100.0
            if pooled["methylated_cpg"] + pooled["unmethylated_cpg"]
            else None
        ),
        "chg_methylation_rate": format_rate(
            context_rate(pooled["methylated_chg"], pooled["unmethylated_chg"]) / 100.0
            if pooled["methylated_chg"] + pooled["unmethylated_chg"]
            else None
        ),
        "chh_methylation_rate": format_rate(
            context_rate(pooled["methylated_chh"], pooled["unmethylated_chh"]) / 100.0
            if pooled["methylated_chh"] + pooled["unmethylated_chh"]
            else None
        ),
        "uniquereads": uniquereads,
        "cpg_methylation_pct": context_rate(
            pooled["methylated_cpg"], pooled["unmethylated_cpg"]
        ),
        "chg_methylation_pct": context_rate(
            pooled["methylated_chg"], pooled["unmethylated_chg"]
        ),
        "chh_methylation_pct": context_rate(
            pooled["methylated_chh"], pooled["unmethylated_chh"]
        ),
    }


def load_cell_reads_methylation_only(work_path: Path) -> dict[str, int]:
    reads_path = work_path / "cells" / "filtered_barcode_read_counts.csv"
    mapping: dict[str, int] = {}
    with reads_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            barcode = (row.get("barcode") or "").strip()
            if not barcode:
                continue
            mapping[barcode] = int(float(row.get("aligned_reads") or 0))
    return mapping


def load_cell_reads_gexcb(work_path: Path) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for reads_path in sorted(
        work_path.glob("split_bams/merged/*_merge_filtered_barcode_reads_counts.csv")
    ):
        with reads_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                barcode = (row.get("barcode") or "").strip()
                if not barcode:
                    continue
                reads = int(float(row.get("reads_counts") or 0))
                mapping[barcode] = mapping.get(barcode, 0) + reads
    return mapping


def load_gex_cb_map(cbcsv: str | None) -> dict[str, str] | None:
    if not cbcsv:
        return None
    mapping: dict[str, str] = {}
    with Path(cbcsv).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            m_cb = (row.get("m_cb") or "").strip()
            gex_cb = (row.get("gex_cb") or "").strip()
            if m_cb:
                mapping[m_cb] = gex_cb
    return mapping


def extract_cell_barcode(count_path: Path) -> str:
    name = count_path.name
    if name.endswith("_allc.gz.count.csv"):
        return name[: -len("_allc.gz.count.csv")]
    if name.endswith(".count.csv"):
        return name[: -len(".count.csv")]
    return count_path.stem


def read_count_csv(count_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with count_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header:
            return rows
        for line in reader:
            if not line:
                continue
            context = line[0].strip()
            if not context:
                continue
            row = {"context": context}
            for index, column in enumerate(header[1:], start=1):
                if index < len(line):
                    row[column] = line[index]
            rows.append(row)
    return rows


def to_float(value: str | None) -> float:
    if value is None or value == "":
        return 0.0
    return float(value)


def to_int(value: str | None) -> int:
    if value is None or value == "":
        return 0
    return int(float(value))


def aggregated_mc_rate(mc: float, cov: float) -> float:
    return mc / cov if cov > 0 else 0.0


def compute_aggregated_mc_rates(context_rows: list[dict[str, str]]) -> dict[str, float]:
    totals = {column: {"mc": 0.0, "cov": 0.0} for column in AGGREGATED_MC_RATE_COLUMNS}
    for row in context_rows:
        context = row["context"]
        if len(context) < 2:
            continue
        mc = to_float(row.get("mc"))
        cov = to_float(row.get("cov"))
        second = context[1]
        third = context[2] if len(context) >= 3 else ""
        if second == "G":
            totals["CG_mc_rate"]["mc"] += mc
            totals["CG_mc_rate"]["cov"] += cov
        else:
            totals["CH_mc_rate"]["mc"] += mc
            totals["CH_mc_rate"]["cov"] += cov
            if third == "G":
                totals["CHG_mc_rate"]["mc"] += mc
                totals["CHG_mc_rate"]["cov"] += cov
            else:
                totals["CHH_mc_rate"]["mc"] += mc
                totals["CHH_mc_rate"]["cov"] += cov
            if second == "A":
                totals["CA_mc_rate"]["mc"] += mc
                totals["CA_mc_rate"]["cov"] += cov
            elif second == "T":
                totals["CT_mc_rate"]["mc"] += mc
                totals["CT_mc_rate"]["cov"] += cov
            elif second == "C":
                totals["CC_mc_rate"]["mc"] += mc
                totals["CC_mc_rate"]["cov"] += cov
    return {
        column: aggregated_mc_rate(totals[column]["mc"], totals[column]["cov"])
        for column in AGGREGATED_MC_RATE_COLUMNS
    }


def build_cells_summary(
    allcools_dir: Path,
    reads_map: dict[str, int],
    ctot_map: dict[str, str],
    gex_cb_map: dict[str, str] | None,
) -> list[dict[str, str]]:
    count_paths = sorted(allcools_dir.glob("*/*_allc.gz.count.csv"))
    cell_rows: list[dict[str, str]] = []

    for count_path in count_paths:
        barcode = extract_cell_barcode(count_path)
        context_rows = read_count_csv(count_path)
        if not context_rows:
            continue

        total_cpg_number = 0
        genome_cov = ""
        genome_cov_raw_umi = ""
        genome_cov_new_umi = ""
        cell_saturation = ""

        for row in context_rows:
            context = row["context"]
            number = to_int(row.get("number"))
            if context in CPG_CONTEXTS:
                total_cpg_number += number
            if not genome_cov:
                genome_cov = row.get("genome_cov") or ""
                genome_cov_raw_umi = row.get("genome_cov_raw_umi") or ""
                genome_cov_new_umi = row.get("genome_cov_new_umi") or ""
                cell_saturation = row.get("cell_saturation") or ""

        aggregated_rates = compute_aggregated_mc_rates(context_rows)
        row_out: dict[str, str] = {
            "cell_barcode": barcode,
            "aligned_reads": str(reads_map.get(barcode, "")),
            "CtoT": ctot_map.get(barcode, ""),
            "total_cpg_number": str(total_cpg_number),
            "genome_cov": genome_cov,
            "genome_cov_raw_umi": genome_cov_raw_umi,
            "genome_cov_new_umi": genome_cov_new_umi,
            "cell_saturation": cell_saturation,
            **{
                column: format_rate(aggregated_rates[column])
                for column in AGGREGATED_MC_RATE_COLUMNS
            },
        }
        if gex_cb_map is not None:
            row_out["gex_cb"] = gex_cb_map.get(barcode, "")
        cell_rows.append(row_out)

    return cell_rows


def median_cell_stats(
    cell_rows: list[dict[str, str]],
    reads_map: dict[str, int],
) -> dict[str, str]:
    eligible = [
        row
        for row in cell_rows
        if row["cell_barcode"] in reads_map and row.get("genome_cov") not in ("", NA)
    ]
    if not eligible:
        return {
            "median_genome_cov": NA,
            "median_total_cpg_number": NA,
            "median_aligned_reads": NA,
            "median_cell_saturation": NA,
        }
    eligible.sort(key=lambda row: float(row["genome_cov"]), reverse=True)
    median_row = eligible[len(eligible) // 2]
    barcode = median_row["cell_barcode"]
    return {
        "median_genome_cov": median_row.get("genome_cov", NA),
        "median_total_cpg_number": median_row.get("total_cpg_number", NA),
        "median_aligned_reads": str(reads_map.get(barcode, NA)),
        "median_cell_saturation": median_row.get("cell_saturation", NA),
    }


def load_saturation_summary(saturation_path: Path, sample_id: str) -> dict[str, str]:
    if not saturation_path.is_file():
        return {
            "observed_median_genome_fraction": NA,
            "theoretical_max_median_genome_fraction": NA,
            "saturation_rate": NA,
            "extrapolation_model": NA,
            "sampled_cell_count": NA,
            "sample_seed": NA,
        }
    with saturation_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if (row.get("sample_id") or "").strip() != sample_id:
                continue
            return {
                "observed_median_genome_fraction": (
                    row.get("observed_median_genome_fraction") or NA
                ),
                "theoretical_max_median_genome_fraction": (
                    row.get("theoretical_max_median_genome_fraction") or NA
                ),
                "saturation_rate": row.get("saturation_rate") or NA,
                "extrapolation_model": row.get("extrapolation_model") or NA,
                "sampled_cell_count": row.get("sampled_cell_count") or NA,
                "sample_seed": row.get("sample_seed") or NA,
            }
    return {
        "observed_median_genome_fraction": NA,
        "theoretical_max_median_genome_fraction": NA,
        "saturation_rate": NA,
        "extrapolation_model": NA,
        "sampled_cell_count": NA,
        "sample_seed": NA,
    }


def build_sample_summary_row(
    sample_id: str,
    fastp: dict,
    demux: dict,
    bismark: dict,
    reads_map: dict[str, int],
    saturation: dict[str, str],
    median_stats: dict[str, str],
) -> dict[str, str]:
    reads_in_cells = sum(reads_map.values())
    uniquereads = int(bismark.get("uniquereads") or 0)
    fraction = reads_in_cells / uniquereads if uniquereads > 0 else None
    row = {
        "sample_id": sample_id,
        **{key: str(fastp[key]) for key in ("raw_reads", "total_bases", "duplication_rate")},
        **{key: str(demux[key]) for key in demux},
        "mapped_to_genome": str(bismark.get("mapped_to_genome", NA)),
        "confidently_mapped": str(bismark.get("confidently_mapped", NA)),
        "cpg_methylation_rate": str(bismark.get("cpg_methylation_rate", NA)),
        "chg_methylation_rate": str(bismark.get("chg_methylation_rate", NA)),
        "chh_methylation_rate": str(bismark.get("chh_methylation_rate", NA)),
        "estimated_cells": str(len(reads_map)),
        "reads_in_cells": str(reads_in_cells),
        "fraction_reads_in_cells": format_rate(fraction),
        **saturation,
        **median_stats,
    }
    return row


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_rate_text(value: str) -> float | None:
    if value in ("", NA):
        return None
    return float(value)


def write_wgs_summary_csv(path: Path, sample_row: dict[str, str], bismark: dict) -> None:
    valid_ratio = parse_rate_text(sample_row.get("valid_barcode_rate", NA))
    dropped_too_short = parse_rate_text(sample_row.get("dropped_too_short", NA))
    dropped_chimeric = parse_rate_text(sample_row.get("dropped_chimeric", NA))
    rate_17lme = parse_rate_text(sample_row.get("rate_17lme", NA))
    ctot = parse_rate_text(sample_row.get("CtoT", NA))
    mapped = parse_rate_text(sample_row.get("mapped_to_genome", NA))
    confident = parse_rate_text(sample_row.get("confidently_mapped", NA))
    fraction = parse_rate_text(sample_row.get("fraction_reads_in_cells", NA))
    median_genome_cov = parse_rate_text(sample_row.get("median_genome_cov", NA))
    median_saturation = parse_rate_text(sample_row.get("median_cell_saturation", NA))

    wgs_row = [
        sample_row["sample_id"],
        sample_row.get("estimated_cells", NA),
        sample_row.get("raw_reads", NA),
        format_pct(valid_ratio),
        format_pct(dropped_too_short),
        format_pct(dropped_chimeric),
        NA,
        format_pct(rate_17lme),
        NA,
        format_pct(ctot),
        NA,
        format_pct(mapped),
        format_pct(confident),
        format_pct_one_decimal(bismark.get("cpg_methylation_pct")),
        format_pct_one_decimal(bismark.get("chg_methylation_pct")),
        format_pct_one_decimal(bismark.get("chh_methylation_pct")),
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        format_pct(median_genome_cov),
        sample_row.get("median_total_cpg_number", NA),
        sample_row.get("median_aligned_reads", NA),
        format_pct(median_saturation),
        format_pct(fraction),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(",".join(WGS_HEADER) + "\n")
        handle.write(",".join(str(value).replace(",", "") for value in wgs_row) + "\n")


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    output_dir = work_path / args.output_dir
    fastp_json = work_path / "shard_fastq" / "fastp.json"
    demux_dir = work_path / "demux"
    align_dir = work_path / "align"
    allcools_dir = work_path / "allcools"
    saturation_summary = work_path / "qc" / "saturation" / "saturation_summary.tsv"
    ctot_path = demux_dir / "qc.CtoT.tsv"

    cells_out = output_dir / "cells_summary.tsv"
    sample_out = output_dir / "sample_summary.tsv"
    wgs_out = output_dir / "wgs_summary.csv"

    print(f"[qc_summary] work_path={work_path}")
    print(f"[qc_summary] sample_id={args.sample_id}")
    print(f"[qc_summary] barcode_mode={args.barcode_mode}")
    print(f"[qc_summary] fastp_json={fastp_json}")
    print(f"[qc_summary] demux_dir={demux_dir}")
    print(f"[qc_summary] align_dir={align_dir}")
    print(f"[qc_summary] allcools_dir={allcools_dir}")
    print(f"[qc_summary] saturation_summary={saturation_summary}")
    print(f"[qc_summary] output_cells={cells_out}")
    print(f"[qc_summary] output_sample={sample_out}")
    print(f"[qc_summary] output_wgs={wgs_out}")

    if args.dry_run:
        print("[qc_summary] dry_run=1")
        return 0

    if not fastp_json.is_file():
        raise FileNotFoundError(f"fastp JSON not found: {fastp_json}")
    if not allcools_dir.is_dir():
        raise FileNotFoundError(f"allcools directory not found: {allcools_dir}")
    if not saturation_summary.is_file():
        raise FileNotFoundError(f"saturation summary not found: {saturation_summary}")

    if args.barcode_mode == "gexcb":
        reads_map = load_cell_reads_gexcb(work_path)
    else:
        reads_path = work_path / "cells" / "filtered_barcode_read_counts.csv"
        if not reads_path.is_file():
            raise FileNotFoundError(f"cell reads table not found: {reads_path}")
        reads_map = load_cell_reads_methylation_only(work_path)

    fastp = parse_fastp_metrics(fastp_json)
    demux = aggregate_demux_stats(demux_dir)
    bismark = aggregate_bismark_reports(align_dir)
    ctot_map = load_ctot_map(ctot_path)
    gex_cb_map = load_gex_cb_map(args.cbcsv)

    cell_rows = build_cells_summary(
        allcools_dir,
        reads_map,
        ctot_map,
        gex_cb_map,
    )
    median_stats = median_cell_stats(cell_rows, reads_map)
    saturation = load_saturation_summary(saturation_summary, args.sample_id)
    sample_row = build_sample_summary_row(
        args.sample_id,
        fastp,
        demux,
        bismark,
        reads_map,
        saturation,
        median_stats,
    )

    cell_columns = list(CORE_CELL_COLUMNS)
    if gex_cb_map is not None:
        cell_columns = [*cell_columns, "gex_cb"]

    write_tsv(cells_out, cell_columns, cell_rows)
    write_tsv(sample_out, list(SAMPLE_SUMMARY_COLUMNS), [sample_row])
    write_wgs_summary_csv(wgs_out, sample_row, bismark)

    print(f"[qc_summary] cells_rows={len(cell_rows)}")
    print("[qc_summary] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
