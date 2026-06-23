#!/usr/bin/env python3
"""Estimate genome-fraction saturation curve from pre-dedup per-cell BAMs."""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


DEFAULT_FRACTIONS = [
    0.01,
    0.02,
    0.05,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    1.00,
]
DEFAULT_MAX_CELLS = 100
DEFAULT_SAMPLE_SEED = 42
DEFAULT_LINEAR_R2_THRESHOLD = 0.99


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate per-cell genome-fraction saturation curve from pre-dedup "
            "merged BAMs and write plot/summary under <work_path>/qc/saturation."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help=(
            "Sample work directory containing split_bams/merged/*_merged_fr_bam "
            "and cell read-count tables."
        ),
    )
    parser.add_argument(
        "--chrom-size-path",
        required=True,
        help="Chromosome sizes BED used to compute total genome size.",
    )
    parser.add_argument(
        "--reads-threshold",
        type=float,
        default=100.0,
        help="HQ cell threshold for aligned reads. Default: 100.",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=DEFAULT_MAX_CELLS,
        help="Maximum HQ cells used for estimation. Default: 100.",
    )
    parser.add_argument(
        "--sample-seed",
        type=int,
        default=DEFAULT_SAMPLE_SEED,
        help="Random seed when sampling HQ cells above --max-cells. Default: 42.",
    )
    parser.add_argument(
        "--pred-fraction",
        type=float,
        default=2.0,
        help="Coverage fraction used for forward prediction point. Default: 2.0.",
    )
    parser.add_argument(
        "--linear-r2-threshold",
        type=float,
        default=DEFAULT_LINEAR_R2_THRESHOLD,
        help=(
            "If the linear (through-origin) fit reaches this R^2, use linear "
            "extrapolation; otherwise fall back to the saturation curve. "
            "Default: 0.99."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <work_path>/qc/saturation.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved inputs/outputs and exit without writing files.",
    )
    return parser.parse_args()


def format_optional_float(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.6f}"


def parse_reads_csv(path: Path) -> dict[str, float]:
    rows: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return rows
        for row in reader:
            barcode = (row.get("barcode") or "").strip()
            if not barcode:
                continue
            reads_text = (
                (row.get("aligned_reads") or row.get("reads_counts") or "").strip()
            )
            if not reads_text:
                continue
            try:
                reads = float(reads_text)
            except ValueError:
                continue
            if reads > 0:
                rows[barcode] = rows.get(barcode, 0.0) + reads
    return rows


def parse_cell_reads(work_path: Path) -> dict[str, float]:
    filtered_counts = work_path / "cells" / "filtered_barcode_read_counts.csv"
    if filtered_counts.is_file():
        return parse_reads_csv(filtered_counts)

    merged_root = work_path / "split_bams" / "merged"
    totals: dict[str, float] = {}
    for path in sorted(merged_root.glob("*_merge_filtered_barcode_reads_counts.csv")):
        for barcode, reads in parse_reads_csv(path).items():
            totals[barcode] = totals.get(barcode, 0.0) + reads
    return totals


def discover_cell_bam_map(work_path: Path) -> dict[str, list[Path]]:
    merged_root = work_path / "split_bams" / "merged"
    cell_bams: dict[str, list[Path]] = defaultdict(list)
    for bam_path in sorted(merged_root.glob("*_merged_fr_bam/*.bam")):
        barcode = bam_path.stem
        cell_bams[barcode].append(bam_path)
    return dict(cell_bams)


def parse_genome_size(chrom_size_path: Path) -> int:
    total = 0
    with chrom_size_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                if len(fields) >= 3:
                    start = int(fields[1])
                    end = int(fields[2])
                    total += end - start
                else:
                    total += int(fields[1])
            except ValueError:
                continue
    if total <= 0:
        raise ValueError(f"invalid or empty genome size from {chrom_size_path}")
    return total


def build_genome_depth_histogram(bam_paths: list[Path]) -> dict[int, int]:
    # Depth is counted per sequencing FRAGMENT (read pair), not per mate.
    # In PE data the two mates of a fragment share a query name and are
    # sampled together, so overlapping mates must contribute at most one to a
    # position's depth; counting each mate separately inflates depth and makes
    # the subsampling curve look more saturated than it is.
    fragment_positions: dict[str, set[tuple[str, int]]] = defaultdict(set)
    for bam_path in bam_paths:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary:
                    continue
                if read.is_unmapped or read.is_qcfail:
                    continue
                ref_name = read.reference_name
                if ref_name is None:
                    continue
                fragment = fragment_positions[read.query_name]
                for pos in read.get_reference_positions(full_length=False):
                    if pos is None:
                        continue
                    fragment.add((ref_name, pos))

    position_depth: dict[tuple[str, int], int] = defaultdict(int)
    for positions in fragment_positions.values():
        for key in positions:
            position_depth[key] += 1

    hist: dict[int, int] = defaultdict(int)
    for depth in position_depth.values():
        if depth > 0:
            hist[depth] += 1
    return dict(hist)


def expected_unique(hist: dict[int, int], fraction: float) -> float:
    total = 0.0
    for depth, count in hist.items():
        total += count * (1.0 - (1.0 - fraction) ** depth)
    return total


def genome_fraction(hist: dict[int, int], fraction: float, genome_size: int) -> float:
    return expected_unique(hist, fraction) / float(genome_size)


def select_sample_cells(
    hq_cells: list[str],
    max_cells: int,
    seed: int,
) -> list[str]:
    if len(hq_cells) <= max_cells:
        return hq_cells
    return sorted(random.Random(seed).sample(hq_cells, max_cells))


def median_and_iqr(values: list[float]) -> tuple[float, float, float]:
    if not values:
        raise ValueError("empty values for median_and_iqr")
    median = float(statistics.median(values))
    if len(values) == 1:
        return median, 0.0, 0.0
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return median, max(0.0, median - q1), max(0.0, q3 - median)


def sat_func(fraction: float, a: float, b: float) -> float:
    return a * (1.0 - math.exp(-b * fraction))


def fit_linear_through_origin(fractions: list[float], y_values: list[float]) -> float:
    denominator = sum(f * f for f in fractions)
    if denominator <= 0:
        return 0.0
    return sum(f * y for f, y in zip(fractions, y_values)) / denominator


def r_squared(y_values: list[float], predictions: list[float]) -> float:
    if not y_values:
        return 0.0
    mean_y = sum(y_values) / len(y_values)
    sst = sum((y - mean_y) ** 2 for y in y_values)
    sse = sum((y - p) ** 2 for y, p in zip(y_values, predictions))
    if sst <= 0:
        return 1.0 if sse <= 0 else 0.0
    return 1.0 - sse / sst


def fit_saturation_curve(fractions: list[float], y_values: list[float]) -> tuple[float, float]:
    if len(fractions) != len(y_values):
        raise ValueError("fractions and y_values length mismatch")
    if not fractions:
        raise ValueError("empty inputs for fitting")

    best_a = 0.0
    best_b = 1.0
    best_error = math.inf

    def search_b(log10_low: float, log10_high: float, count: int) -> None:
        nonlocal best_a, best_b, best_error
        if count <= 1:
            return
        step = (log10_high - log10_low) / float(count - 1)
        for index in range(count):
            b = 10 ** (log10_low + step * index)
            transformed = [1.0 - math.exp(-b * f) for f in fractions]
            denominator = sum(value * value for value in transformed)
            if denominator <= 0:
                continue
            a = sum(y * value for y, value in zip(y_values, transformed)) / denominator
            a = max(0.0, a)
            error = sum(
                (y - a * value) * (y - a * value)
                for y, value in zip(y_values, transformed)
            )
            if error < best_error:
                best_error = error
                best_a = a
                best_b = b

    search_b(log10_low=-4.0, log10_high=2.0, count=800)
    for _ in range(3):
        center = math.log10(best_b if best_b > 0 else 1.0)
        search_b(log10_low=center - 0.7, log10_high=center + 0.7, count=240)

    return best_a, best_b


def write_summary_tsv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "observed_median_genome_fraction",
        "theoretical_max_median_genome_fraction",
        "predicted_median_genome_fraction_at_2x",
        "saturation_rate",
        "extrapolation_model",
        "hq_cell_count",
        "sampled_cell_count",
        "sample_seed",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_empty_plot(path: Path, sample_id: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: NA")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def write_plot(
    path: Path,
    sample_id: str,
    fractions: list[float],
    medians: list[float],
    err_low: list[float],
    err_high: list[float],
    model: str,
    a: float,
    b: float,
    slope: float,
    pred_fraction: float,
    predicted: float | None,
    theoretical: float | None,
    saturation_rate: float | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scale = 100.0

    max_x = max(max(fractions), pred_fraction, 1.0)
    fit_x = [max_x * index / 250.0 for index in range(251)]
    if model == "linear":
        fit_y = [slope * value * scale for value in fit_x]
        fit_label = "Fitted line (linear)"
    else:
        fit_y = [sat_func(value, a, b) * scale for value in fit_x]
        fit_label = "Fitted curve (saturation)"
    pred_y = (predicted if predicted is not None else 0.0) * scale
    medians_pct = [value * scale for value in medians]
    yerr = [
        [value * scale for value in err_low],
        [value * scale for value in err_high],
    ]

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(
        fractions,
        medians_pct,
        yerr=yerr,
        fmt="o-",
        color="blue",
        linewidth=2,
        markersize=5,
        capsize=4,
        label="Observed median (IQR)",
    )
    ax.plot(fit_x, fit_y, "r--", linewidth=2, label=fit_label)
    ax.scatter(
        [pred_fraction],
        [pred_y],
        color="green",
        s=55,
        zorder=4,
        label=f"Prediction at {pred_fraction:g}x",
    )
    if model != "linear" and theoretical is not None:
        ax.axhline(
            theoretical * scale,
            color="purple",
            linestyle="--",
            linewidth=1.8,
            label="Max genome fraction",
        )
    if saturation_rate is not None:
        sat_text = f"{saturation_rate:.2f}%"
    elif model == "linear":
        sat_text = "NA (linear, unsaturated)"
    else:
        sat_text = "NA"
    ax.set_title(f"Saturation analysis ({sample_id})\nSaturation rate: {sat_text}")
    ax.set_xlabel("Coverage Fraction")
    ax.set_ylabel("Median Genome Fraction (%)")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def empty_summary_row(
    sample_id: str,
    hq_cell_count: int,
    sample_seed: int,
    sampled_cell_count: int = 0,
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "observed_median_genome_fraction": "NA",
        "theoretical_max_median_genome_fraction": "NA",
        "predicted_median_genome_fraction_at_2x": "NA",
        "saturation_rate": "NA",
        "extrapolation_model": "NA",
        "hq_cell_count": str(hq_cell_count),
        "sampled_cell_count": str(sampled_cell_count),
        "sample_seed": str(sample_seed),
    }


def main() -> int:
    args = parse_args()
    work_path = Path(args.work_path)
    chrom_size_path = Path(args.chrom_size_path)
    sample_id = work_path.name
    output_dir = Path(args.output_dir) if args.output_dir else work_path / "qc" / "saturation"
    plot_path = output_dir / "saturation_curve.png"
    summary_path = output_dir / "saturation_summary.tsv"

    if args.max_cells <= 0:
        raise ValueError("--max-cells must be > 0")

    genome_size = parse_genome_size(chrom_size_path)
    read_counts = parse_cell_reads(work_path)
    cell_bam_map = discover_cell_bam_map(work_path)

    print(f"[saturation] work_path={work_path}")
    print(f"[saturation] chrom_size_path={chrom_size_path}")
    print(f"[saturation] genome_size={genome_size}")
    print(f"[saturation] cell_reads_count={len(read_counts)}")
    print(f"[saturation] cell_bam_count={len(cell_bam_map)}")
    print(f"[saturation] reads_threshold={args.reads_threshold}")
    print(f"[saturation] max_cells={args.max_cells}")
    print(f"[saturation] sample_seed={args.sample_seed}")
    print(f"[saturation] output_plot={plot_path}")
    print(f"[saturation] output_summary={summary_path}")

    hq_cells = sorted(
        barcode
        for barcode, reads in read_counts.items()
        if reads > args.reads_threshold and barcode in cell_bam_map
    )
    sampled_cells = select_sample_cells(hq_cells, args.max_cells, args.sample_seed)

    print(f"[saturation] hq_cell_count={len(hq_cells)}")
    print(f"[saturation] sampled_cell_count={len(sampled_cells)}")

    if args.dry_run:
        for barcode in sampled_cells[:5]:
            chunk_count = len(cell_bam_map[barcode])
            print(f"[saturation] sampled_cell={barcode} chunk_bams={chunk_count}")
        print("[saturation] dry_run=1")
        return 0

    if not hq_cells:
        print("[saturation] warning=no_hq_cells_after_filter")
        write_summary_tsv(summary_path, empty_summary_row(sample_id, 0, args.sample_seed))
        write_empty_plot(plot_path, sample_id, "No HQ cells after reads filter")
        print("[saturation] done")
        return 0

    per_fraction_values: dict[float, list[float]] = {
        fraction: [] for fraction in DEFAULT_FRACTIONS
    }
    for barcode in sampled_cells:
        hist = build_genome_depth_histogram(cell_bam_map[barcode])
        if not hist:
            continue
        for fraction in DEFAULT_FRACTIONS:
            per_fraction_values[fraction].append(
                genome_fraction(hist, fraction, genome_size)
            )

    usable_fractions: list[float] = []
    medians: list[float] = []
    err_low: list[float] = []
    err_high: list[float] = []
    for fraction in DEFAULT_FRACTIONS:
        values = per_fraction_values[fraction]
        if not values:
            continue
        median, low, high = median_and_iqr(values)
        usable_fractions.append(fraction)
        medians.append(median)
        err_low.append(low)
        err_high.append(high)

    if not usable_fractions:
        print("[saturation] warning=no_genome_fraction_values_for_sampled_cells")
        write_summary_tsv(
            summary_path,
            empty_summary_row(
                sample_id,
                len(hq_cells),
                args.sample_seed,
                len(sampled_cells),
            ),
        )
        write_empty_plot(plot_path, sample_id, "No valid genome depth histograms")
        print("[saturation] done")
        return 0

    observed = medians[-1]
    a, b = fit_saturation_curve(usable_fractions, medians)
    slope = fit_linear_through_origin(usable_fractions, medians)
    linear_r2 = r_squared(medians, [slope * f for f in usable_fractions])
    exp_r2 = r_squared(medians, [sat_func(f, a, b) for f in usable_fractions])

    # Prefer the linear (through-origin) extrapolation when the observed curve
    # is essentially linear (undersaturated): in that regime the exponential
    # asymptote/saturation rate is unidentifiable and unreliable. Only when the
    # linear fit is poor (clear curvature) do we trust the saturation curve.
    use_linear = linear_r2 >= args.linear_r2_threshold
    if use_linear:
        model = "linear"
        theoretical = None
        predicted = slope * args.pred_fraction
        saturation_rate = None
    else:
        model = "saturation"
        theoretical = a if a > 0 else None
        predicted = sat_func(args.pred_fraction, a, b) if a > 0 else None
        saturation_rate = None
        if theoretical is not None and theoretical > 0:
            saturation_rate = observed / theoretical * 100.0

    print(f"[saturation] linear_r2={linear_r2:.6f} exp_r2={exp_r2:.6f}")
    print(f"[saturation] extrapolation_model={model}")
    print(f"[saturation] saturation_rate={format_optional_float(saturation_rate)}")

    write_plot(
        path=plot_path,
        sample_id=sample_id,
        fractions=usable_fractions,
        medians=medians,
        err_low=err_low,
        err_high=err_high,
        model=model,
        a=a,
        b=b,
        slope=slope,
        pred_fraction=args.pred_fraction,
        predicted=predicted,
        theoretical=theoretical,
        saturation_rate=saturation_rate,
    )
    row = {
        "sample_id": sample_id,
        "observed_median_genome_fraction": format_optional_float(observed),
        "theoretical_max_median_genome_fraction": format_optional_float(theoretical),
        "predicted_median_genome_fraction_at_2x": format_optional_float(predicted),
        "saturation_rate": format_optional_float(saturation_rate),
        "extrapolation_model": model,
        "hq_cell_count": str(len(hq_cells)),
        "sampled_cell_count": str(len(sampled_cells)),
        "sample_seed": str(args.sample_seed),
    }
    write_summary_tsv(summary_path, row)
    print("[saturation] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
