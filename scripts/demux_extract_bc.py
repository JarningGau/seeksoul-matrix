#!/usr/bin/env python3
"""DD-MET5 demux: extract CB/UB, C→T QC, adapter trim, forward/reverse FASTQ output."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Iterator, TextIO

from cutadapt.adapters import BackAdapter, FrontAdapter
from tqdm import tqdm

ME5_POSITIONS = [42, 45, 46, 49, 54, 56, 75]
ME5_NON_INSERT_LEN = 70
ME5_L17ME_START = 57
ME5_L17ME_LEN = 18
ME5_L17ME_SEQ = "GTAGATGTGTATAAGAGA"
ME5_MAX_READ1_LEN = 63
R1_MINLEN = 20
R2_MINLEN = 60
STRUCTURE = "B17U12"

DD_MET5_ADAPTER1 = [
    ["AGATGTGTATAAGAGAYAG", "5", 0.1, 9],
    ["CTGTCTCTTATACACATCT", "3", 0.1, 9],
]
DD_MET5_ADAPTER2 = [
    ["AGATGTGTATAAGAGACAG", "5", 0.1, 9],
    ["CTRTCTCTTATACACATCT", "3", 0.1, 9],
]

# CH chimeric filtering (SeekSoulMethyl filter_ch)
R1_FORWARD_CH_PATTERN = re.compile(r"C[ATC]")
R2_FORWARD_CH_PATTERN = re.compile(r"[ATG]G")
R1_REVERSE_CH_PATTERN = re.compile(r"[ATG]G")
R2_REVERSE_CH_PATTERN = re.compile(r"C[ATC]")


@dataclass
class FastqRead:
    name: str
    sequence: str
    qualities: str

    def __getitem__(self, slc: slice) -> FastqRead:
        return FastqRead(self.name, self.sequence[slc], self.qualities[slc])


class AdapterFilter:
    """Trim ME adapters and crop non-insert regions (ported from SeekSoulMethyl)."""

    def __init__(
        self,
        adapter1: list[list],
        adapter2: list[list],
        non_insert_len: int = ME5_NON_INSERT_LEN,
    ) -> None:
        self.adapter1 = self._build_adapters(adapter1)
        self.adapter2 = self._build_adapters(adapter2)
        self.non_insert_len = non_insert_len

    @staticmethod
    def _build_adapters(specs: list[list]) -> list:
        adapters = []
        for p in specs:
            if len(p) >= 4:
                if p[1] == "3":
                    adapters.append(
                        BackAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3])
                    )
                elif p[1] == "5":
                    adapters.append(
                        FrontAdapter(sequence=p[0], max_errors=p[2], min_overlap=p[3])
                    )
            elif p[1] == "3":
                adapters.append(BackAdapter(sequence=p[0], min_overlap=10))
            elif p[1] == "5":
                adapters.append(FrontAdapter(sequence=p[0], min_overlap=7))
        return adapters

    def filter(self, r1: FastqRead, r2: FastqRead) -> tuple[bool, FastqRead, FastqRead]:
        flag = False
        r1_me_left = False
        r1_me_right = False
        if self.adapter1:
            for adapter in self.adapter1:
                match = adapter.match_to(r1.sequence)
                if match:
                    if adapter.sequence == "AGATGTGTATAAGAGAYAG":
                        r1_me_left = True
                    if adapter.sequence == "CTGTCTCTTATACACATCT":
                        r1_me_right = True
                    flag = True
                    trimmed = match.trimmed(r1)
                    r1 = FastqRead(r1.name, trimmed.sequence, trimmed.qualities)
            r1_start = 0
            r1_end = len(r1.sequence)
            if len(r1.sequence) > 18:
                if r1_me_left:
                    r1_start = 9
                else:
                    r1_start = self.non_insert_len
                if r1_me_right:
                    r1_end = len(r1.sequence) - 9
                if r1_end - r1_start > ME5_MAX_READ1_LEN:
                    r1_end = r1_start + ME5_MAX_READ1_LEN
                r1 = FastqRead(
                    r1.name, r1.sequence[r1_start:r1_end], r1.qualities[r1_start:r1_end]
                )

        r2_me_right = False
        if self.adapter2:
            for adapter in self.adapter2:
                match = adapter.match_to(r2.sequence)
                if match:
                    if adapter.sequence == "CTRTCTCTTATACACATCT":
                        r2_me_right = True
                    flag = True
                    trimmed = match.trimmed(r2)
                    r2 = FastqRead(r2.name, trimmed.sequence, trimmed.qualities)
            r2_start = 9
            r2_end = len(r2.sequence)
            if len(r2.sequence) > 18 and r2_me_right:
                r2_end = len(r2.sequence) - 9
            r2 = FastqRead(
                r2.name, r2.sequence[r2_start:r2_end], r2.qualities[r2_start:r2_end]
            )
        return flag, r1, r2


def parse_structure(string: str) -> tuple[tuple[str, int], ...]:
    regex = re.compile(r"([BLUXT])(\d+)")
    return tuple((code, int(length)) for code, length in regex.findall(string))


def open_text(path: str | Path, mode: str) -> TextIO:
    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, mode + "t", encoding="utf-8")
    return open(path_str, mode, encoding="utf-8")


def read_whitelist(path: str | Path) -> set[str]:
    barcodes: set[str] = set()
    with open_text(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip().upper()
            if line and not line.startswith("#"):
                barcodes.add(line)
    if not barcodes:
        raise ValueError(f"empty barcode whitelist: {path}")
    return barcodes


def get_new_bc(bc: str, white_list: set[str], distance: int = 1) -> set[str]:
    if distance == 1:
        base_list = ["T", "C", "G", "A"]
        mm_dict: dict[str, str] = {}
        for i, base in enumerate(bc):
            if base == "N":
                mm_dict = {bc[:i] + alt + bc[i + 1 :]: f"{i}{alt}" for alt in base_list}
                break
            for alt in base_list:
                if alt != base:
                    mm_dict[bc[:i] + alt + bc[i + 1 :]] = f"{i}{alt}"
        return set(mm_dict.keys()).intersection(white_list)

    bc_dict: dict[int, set[str]] = defaultdict(set)
    for bc_true in white_list:
        hmm = sum(ch1 != ch2 for ch1, ch2 in zip(bc_true, bc))
        if hmm <= distance:
            bc_dict[hmm].add(bc_true)
    if not bc_dict:
        return set()
    best_dist = min(bc_dict)
    return bc_dict[best_dist]


def match_barcode(
    observed: str, whitelist: set[str], hamming_distance: int
) -> tuple[str | None, str, str]:
    observed = observed.upper()
    if observed in whitelist:
        return observed, "M", "exact"
    if hamming_distance <= 0:
        return None, "", "not_found"
    candidates = get_new_bc(observed, whitelist, hamming_distance)
    if not candidates:
        return None, "", "not_found"
    if len(candidates) == 1:
        corrected = next(iter(candidates))
        alt = "".join(
            f"{i}{old}"
            for i, (old, new) in enumerate(zip(observed, corrected))
            if old != new
        )
        return corrected, alt, "corrected"
    return None, "", "ambiguous"


def extract_ct_counts(sequence: str, qualities: str) -> tuple[int, int] | None:
    """Return (C, T) counts for C→T QC; only forward (TTT insert) reads qualify."""
    l17me = sequence[ME5_L17ME_START : ME5_L17ME_START + ME5_L17ME_LEN]
    if l17me != ME5_L17ME_SEQ:
        return None
    max_position = max(ME5_POSITIONS)
    if len(sequence) <= max_position:
        return None

    bases = [sequence[pos] for pos in ME5_POSITIONS]
    quals = [qualities[pos] for pos in ME5_POSITIONS]
    allseq = "".join(bases)
    countseq = allseq[1:5]
    insert = bases[0] + bases[-2] + bases[-1]
    avg_q = sum(ord(char) - 33 for char in quals) / len(quals)
    if avg_q < 30:
        return None
    # DD-MET5: only TTT insert reads are used for C→T conversion QC (SeekSoulMethyl).
    if insert != "TTT":
        return None
    c_count = countseq.count("C")
    t_count = countseq.count("T")
    if c_count + t_count != len(countseq):
        return None
    return c_count, t_count


def determine_chain_direction(sequence: str) -> str:
    max_position = max(ME5_POSITIONS)
    if len(sequence) <= max_position:
        return "unknown"
    insert_positions = [ME5_POSITIONS[0], ME5_POSITIONS[5], ME5_POSITIONS[6]]
    insert = "".join(sequence[pos] for pos in insert_positions)
    if insert.count("C") == 3:
        return "reverse"
    return "forward"


def should_filter_read_ch_pattern(
    r1_sequence: str,
    r2_sequence: str,
    chain_direction: str,
    threshold: int,
) -> bool:
    """Return True when read pair exceeds filter_ch CH pattern count (SeekSoulMethyl)."""
    if threshold <= 0:
        return False
    if chain_direction == "forward":
        r1_count = len(R1_FORWARD_CH_PATTERN.findall(r1_sequence))
        r2_count = len(R2_FORWARD_CH_PATTERN.findall(r2_sequence))
    elif chain_direction == "reverse":
        r1_count = len(R1_REVERSE_CH_PATTERN.findall(r1_sequence))
        r2_count = len(R2_REVERSE_CH_PATTERN.findall(r2_sequence))
    else:
        return False
    return r1_count > threshold or r2_count > threshold


def fastq_iter(handle: TextIO) -> Iterator[tuple[str, str, str, str]]:
    while True:
        head = handle.readline()
        if not head:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not (seq and plus and qual):
            raise ValueError("incomplete FASTQ record detected")
        yield head.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def format_read_name(
    header: str, cb: str, umi: str, chain: str, alt: str
) -> str:
    if not header.startswith("@"):
        raise ValueError(f"invalid FASTQ header: {header}")
    parts = header.split(" ", 1)
    orig = parts[0][1:]
    suffix = f" {parts[1]}" if len(parts) > 1 else ""
    new_core = f"{cb}_{umi}_{chain}_{alt}_{orig}:{cb}"
    return f"@{new_core}{suffix}"


def has_17lme(sequence: str) -> bool:
    if len(sequence) < ME5_L17ME_START + ME5_L17ME_LEN:
        return False
    return (
        sequence[ME5_L17ME_START : ME5_L17ME_START + ME5_L17ME_LEN] == ME5_L17ME_SEQ
    )


def safe_fraction(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return round(numerator / denominator, 6)


def build_stats_payload(
    *,
    chunk_id: str,
    input_r1: str,
    input_r2: str,
    stats: dict[str, int],
    ct_c_total: int,
    ct_t_total: int,
    filter_ch: int,
) -> dict:
    total = stats["total"]
    valid = stats["valid"]
    b_exact = stats["B_exact"]
    b_corrected = stats["B_corrected"]
    b_whitelist = b_exact + b_corrected
    ct_convertible = ct_c_total + ct_t_total
    ctot = round_ctot(ct_c_total, ct_t_total)

    chimeric = stats["chimeric_filtered"]
    return {
        "chunk_id": chunk_id,
        "chemistry": "DD-MET5",
        "filter_ch": filter_ch,
        "input_r1": input_r1,
        "input_r2": input_r2,
        "funnel": {
            "total": total,
            "barcode_rejected": stats["B_rejected"],
            "barcode_ambiguous": stats["B_ambiguous"],
            "barcode_passed": {
                "total": b_whitelist,
                "exact": b_exact,
                "corrected": b_corrected,
                "corrected_fraction": safe_fraction(b_corrected, b_whitelist),
                "unknown_chain": stats["unknown_chain"],
                "too_short": stats["too_short"],
                "chimeric_filtered": {
                    "total": chimeric,
                    "forward": stats["forward_chimeric"],
                    "reverse": stats["reverse_chimeric"],
                },
                "valid": {
                    "total": valid,
                    "fraction_of_input": safe_fraction(valid, total),
                    "forward": stats["forward"],
                    "reverse": stats["reverse"],
                },
            },
        },
        "ct": {
            "num_17lme": stats["num_17lme"],
            "rate_17lme": safe_fraction(stats["num_17lme"], total),
            "ct_reads": stats["ct_reads"],
            "ct_umi_dedup": stats["ct_umi_dedup"],
            "ct_convertible_bases": ct_convertible,
            "ct_converted_bases": ct_t_total,
            "CtoT": ctot if ctot is not None else 0.0,
        },
    }


def round_ctot(c_sum: int, t_sum: int) -> float | None:
    total = c_sum + t_sum
    if total == 0:
        return None
    return round(t_sum / total, 3)


def output_paths(prefix: str | Path) -> dict[str, Path]:
    prefix = Path(prefix)
    return {
        "forward_1": prefix.with_name(f"{prefix.name}.forward_1.fq.gz"),
        "forward_2": prefix.with_name(f"{prefix.name}.forward_2.fq.gz"),
        "reverse_1": prefix.with_name(f"{prefix.name}.reverse_1.fq.gz"),
        "reverse_2": prefix.with_name(f"{prefix.name}.reverse_2.fq.gz"),
        "linker": prefix.with_name(f"{prefix.name}.linker.tsv"),
        "stats": prefix.with_name(f"{prefix.name}.stats.json"),
    }


def load_fastp_pair_total(fastp_json: Path) -> int:
    data = json.loads(fastp_json.read_text(encoding="utf-8"))
    read1 = data.get("read1_after_filtering", {}).get("total_reads")
    if read1 is not None:
        return int(read1)
    after = data.get("summary", {}).get("after_filtering", {}).get("total_reads")
    if after is not None:
        return int(after) // 2
    raise ValueError(f"no read count found in fastp report: {fastp_json}")


def count_shard_r1_files(shard_dir: Path) -> int:
    numbered = [
        p
        for p in shard_dir.glob("*.R1.fq.gz")
        if re.fullmatch(r"\d+\.R1\.fq\.gz", p.name)
    ]
    if numbered:
        return len(numbered)
    legacy = list(shard_dir.glob("R1*.fq.gz"))
    if legacy:
        return len(legacy)
    raise ValueError(f"no fastp shard R1 files found under {shard_dir}")


def estimate_chunk_pair_total(fastp_json: Path) -> int:
    shard_dir = fastp_json.parent
    shard_count = count_shard_r1_files(shard_dir)
    pair_total = load_fastp_pair_total(fastp_json)
    return pair_total // shard_count


def discover_fastp_json(r1_path: str | Path) -> Path | None:
    shard_dir = Path(r1_path).resolve().parent
    if shard_dir.name != "shard_fastq":
        return None
    candidate = shard_dir / "fastp.json"
    return candidate if candidate.is_file() else None


def resolve_progress_total(args: argparse.Namespace) -> int | None:
    if args.total_reads is not None:
        return int(args.total_reads)
    fastp_json = Path(args.fastp_json) if args.fastp_json else discover_fastp_json(args.r1)
    if fastp_json is None:
        return None
    try:
        return estimate_chunk_pair_total(fastp_json)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(
            f"[demux_extract_bc] warning: could not estimate progress total "
            f"from {fastp_json}: {exc}",
            file=sys.stderr,
        )
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="DD-MET5 demux: extract barcodes, C→T QC, and paired FASTQ output."
    )
    parser.add_argument("r1", help="Input chunk R1 FASTQ(.gz).")
    parser.add_argument("r2", help="Input chunk R2 FASTQ(.gz).")
    parser.add_argument(
        "--barcode-whitelist",
        required=True,
        help="Cell barcode whitelist (one barcode per line, optionally gzipped).",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix, e.g. work/<sample>/demux/<chunk>.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        default=1,
        help="Max Hamming distance for whitelist correction. Default: 1.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        help="gzip compression level for output FASTQ (0-9). Default: 6.",
    )
    parser.add_argument(
        "--fastp-json",
        help="fastp JSON report for progress-bar total (default: auto-detect in shard_fastq/).",
    )
    parser.add_argument(
        "--total-reads",
        type=int,
        help="Override expected read-pair count for the progress bar.",
    )
    parser.add_argument(
        "--filter-ch",
        type=int,
        default=2,
        help=(
            "CH chimeric filter threshold (0=disabled; >0 drops pairs when R1 or R2 "
            "CH pattern count exceeds this value). Default: 2."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print output paths without reading FASTQ.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paths = output_paths(args.output_prefix)

    print(f"[demux_extract_bc] input_r1={args.r1}")
    print(f"[demux_extract_bc] input_r2={args.r2}")
    print(f"[demux_extract_bc] filter_ch={args.filter_ch}")
    for label, path in paths.items():
        print(f"[demux_extract_bc] output_{label}={path}")

    if args.dry_run:
        return 0

    whitelist = read_whitelist(args.barcode_whitelist)
    structure = parse_structure(STRUCTURE)
    adapter_filter = AdapterFilter(DD_MET5_ADAPTER1, DD_MET5_ADAPTER2)

    paths["forward_1"].parent.mkdir(parents=True, exist_ok=True)

    stats = defaultdict(int)
    ct_c_total = 0
    ct_t_total = 0
    seen_umi: set[tuple[str, str]] = set()
    chunk_id = Path(args.output_prefix).name
    total_pairs = resolve_progress_total(args)

    t0 = time.monotonic()
    with (
        open_text(args.r1, "r") as r1_in,
        open_text(args.r2, "r") as r2_in,
        gzip.open(paths["forward_1"], "wt", encoding="utf-8", compresslevel=args.gzip_level) as fwd1_out,
        gzip.open(paths["forward_2"], "wt", encoding="utf-8", compresslevel=args.gzip_level) as fwd2_out,
        gzip.open(paths["reverse_1"], "wt", encoding="utf-8", compresslevel=args.gzip_level) as rev1_out,
        gzip.open(paths["reverse_2"], "wt", encoding="utf-8", compresslevel=args.gzip_level) as rev2_out,
        open(paths["linker"], "w", encoding="utf-8") as linker_out,
        tqdm(
            total=total_pairs,
            desc=f"demux {chunk_id}",
            unit="pair",
            unit_scale=True,
            file=sys.stderr,
        ) as progress,
    ):
        linker_out.write("CR\tUB\tC\tT\n")
        for rec in zip_longest(fastq_iter(r1_in), fastq_iter(r2_in)):
            if rec[0] is None or rec[1] is None:
                raise ValueError("R1/R2 read count mismatch")
            h1, s1, p1, q1 = rec[0]
            h2, s2, p2, q2 = rec[1]
            stats["total"] += 1
            progress.update(1)
            if has_17lme(s1):
                stats["num_17lme"] += 1

            start_pos = 0
            observed_cb = ""
            umi = ""
            corrected_cb: str | None = None
            cb_alt = ""
            cb_status = "not_found"

            for code, length in structure:
                end_pos = start_pos + length
                seq = s1[start_pos:end_pos]
                if code == "B":
                    observed_cb = seq.upper()
                    corrected_cb, cb_alt, cb_status = match_barcode(
                        observed_cb, whitelist, args.barcode_hamming_distance
                    )
                    if corrected_cb is None:
                        break
                elif code == "U":
                    umi = seq.upper()
                start_pos = end_pos

            if corrected_cb is None:
                if cb_status == "ambiguous":
                    stats["B_ambiguous"] += 1
                else:
                    stats["B_rejected"] += 1
                continue
            if cb_status == "corrected":
                stats["B_corrected"] += 1
            elif cb_status == "exact":
                stats["B_exact"] += 1

            ct_counts = extract_ct_counts(s1, q1)
            if ct_counts is not None:
                stats["ct_reads"] += 1
                c_count, t_count = ct_counts
                umi_key = (corrected_cb, umi)
                if umi_key not in seen_umi:
                    seen_umi.add(umi_key)
                    linker_out.write(
                        f"{corrected_cb}\t{umi}\t{c_count}\t{t_count}\n"
                    )
                    stats["ct_umi_dedup"] += 1
                    ct_c_total += c_count
                    ct_t_total += t_count

            chain = determine_chain_direction(s1)
            if chain == "unknown":
                stats["unknown_chain"] += 1
                continue

            r1 = FastqRead(h1[1:].split()[0], s1, q1)
            r2 = FastqRead(h2[1:].split()[0], s2, q2)
            _, r1, r2 = adapter_filter.filter(r1, r2)

            if len(r1.sequence) < R1_MINLEN or len(r2.sequence) < R2_MINLEN:
                stats["too_short"] += 1
                continue

            if should_filter_read_ch_pattern(
                r1.sequence, r2.sequence, chain, args.filter_ch
            ):
                stats["chimeric_filtered"] += 1
                if chain == "forward":
                    stats["forward_chimeric"] += 1
                else:
                    stats["reverse_chimeric"] += 1
                continue

            stats["valid"] += 1
            if chain == "forward":
                stats["forward"] += 1
                out1, out2 = fwd1_out, fwd2_out
            else:
                stats["reverse"] += 1
                out1, out2 = rev1_out, rev2_out

            new_h1 = format_read_name(h1, corrected_cb, umi, chain, cb_alt)
            new_h2 = format_read_name(h2, corrected_cb, umi, chain, cb_alt)
            out1.write(f"{new_h1}\n{r1.sequence}\n{p1}\n{r1.qualities}\n")
            out2.write(f"{new_h2}\n{r2.sequence}\n{p2}\n{r2.qualities}\n")

    stats_payload = build_stats_payload(
        chunk_id=chunk_id,
        input_r1=args.r1,
        input_r2=args.r2,
        stats=stats,
        ct_c_total=ct_c_total,
        ct_t_total=ct_t_total,
        filter_ch=args.filter_ch,
    )
    paths["stats"].write_text(
        json.dumps(stats_payload, indent=2) + "\n", encoding="utf-8"
    )

    elapsed = max(time.monotonic() - t0, 1e-9)
    print(
        f"[demux_extract_bc] total={stats['total']} valid={stats['valid']} "
        f"forward={stats['forward']} reverse={stats['reverse']} "
        f"B_corrected={stats['B_corrected']} B_rejected={stats['B_rejected']} "
        f"B_ambiguous={stats['B_ambiguous']} too_short={stats['too_short']} "
        f"chimeric_filtered={stats['chimeric_filtered']} "
        f"forward_chimeric={stats['forward_chimeric']} "
        f"reverse_chimeric={stats['reverse_chimeric']} "
        f"ct_umi_dedup={stats['ct_umi_dedup']} "
        f"CtoT={stats_payload['ct']['CtoT']:.3f} "
        f"speed={stats['total'] / elapsed:.1f} reads/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
