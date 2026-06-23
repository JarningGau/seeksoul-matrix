"""Shared path existence checks for workflow drivers (make_cmd entrypoints)."""

from __future__ import annotations

import gzip
import re
from pathlib import Path


def resolve_config_path(s: str) -> Path:
    p = Path(s).expanduser()
    if p.is_absolute():
        return p
    return (Path.cwd() / p).resolve()


def require_file(label: str, path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"missing required file ({label}): {path}")


def require_dir(label: str, path: Path) -> None:
    if not path.is_dir():
        raise ValueError(f"missing required directory ({label}): {path}")


def looks_like_filesystem_path(s: str) -> bool:
    return "/" in s or s.startswith(".") or s.startswith("~")


def require_optional_executable_path(label: str, raw: str) -> None:
    if looks_like_filesystem_path(raw):
        require_file(label, resolve_config_path(raw))


def r1_name_to_chunk_id(name: str) -> str:
    """Map fastp shard R1 filename to zero-padded chunk id."""
    stem = name
    for suffix in (".fq.gz", ".fastq.gz"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    if stem == "R1":
        return "0001"
    match = re.fullmatch(r"R1_(\d+)", stem)
    if match:
        return f"{int(match.group(1)) + 1:04d}"
    raise ValueError(f"unsupported fastp shard R1 name: {name}")


def discover_fastp_shards(shard_dir: Path) -> list[tuple[str, Path, Path]]:
    """Return sorted (chunk_id, r1_path, r2_path) for fastp split outputs."""
    shard_dir = Path(shard_dir)
    if not shard_dir.is_dir():
        return []

    numbered = sorted(
        (
            p
            for p in shard_dir.glob("*.R1.fq.gz")
            if re.fullmatch(r"\d+\.R1\.fq\.gz", p.name)
        ),
        key=lambda p: p.name,
    )
    if numbered:
        shards: list[tuple[str, Path, Path]] = []
        for r1_path in numbered:
            chunk_id = r1_path.name.split(".", 1)[0]
            r2_path = shard_dir / r1_path.name.replace(".R1.", ".R2.", 1)
            shards.append((chunk_id, r1_path, r2_path))
        return shards

    shards = []
    for r1_path in sorted(shard_dir.glob("R1*.fq.gz")):
        chunk_id = r1_name_to_chunk_id(r1_path.name)
        r2_name = r1_path.name.replace("R1", "R2", 1)
        r2_path = shard_dir / r2_name
        shards.append((chunk_id, r1_path, r2_path))
    return sorted(shards, key=lambda item: item[0])


def build_chunk_names(number_of_split_parts: int) -> list[str]:
    if number_of_split_parts <= 0:
        raise ValueError("number_of_split_parts must be > 0")
    width = max(4, len(str(number_of_split_parts)))
    return [f"{index:0{width}d}" for index in range(1, number_of_split_parts + 1)]


def read_whitelist_barcodes(path: Path) -> set[str]:
    """Load barcode whitelist (plain or gzipped)."""
    path = Path(path)
    barcodes: set[str] = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt"
    with opener(path, mode, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip().upper()
            if line and not line.startswith("#"):
                barcodes.add(line)
    if not barcodes:
        raise ValueError(f"empty barcode whitelist: {path}")
    return barcodes


def plan_prefix_chunks(whitelist_path: Path, prefix_bases: int) -> list[str]:
    """Return sorted barcode prefixes for analysis chunk planning."""
    if prefix_bases <= 0:
        raise ValueError("split_fastq_prefix_bases must be > 0 for prefix chunk planning")
    barcodes = read_whitelist_barcodes(whitelist_path)
    prefixes = sorted(
        {bc[:prefix_bases] for bc in barcodes if len(bc) >= prefix_bases}
    )
    if not prefixes:
        raise ValueError(
            f"no barcode prefixes from whitelist {whitelist_path} "
            f"with split_fastq_prefix_bases={prefix_bases}"
        )
    return prefixes


def discover_demux_subshards(
    demux_dir: Path,
) -> dict[str, dict[str, list[Path]]]:
    """Return {prefix: {stream: [sub-shard paths]}} under demux/shards/."""
    shard_dir = Path(demux_dir) / "shards"
    if not shard_dir.is_dir():
        return {}

    subshard_re = re.compile(
        r"^(?P<readchunk>[^_]+)__(?P<prefix>[A-Z]+)\."
        r"(?P<stream>forward_1|forward_2|reverse_1|reverse_2)\.fq\.gz$"
    )
    streams = ("forward_1", "forward_2", "reverse_1", "reverse_2")
    grouped: dict[str, dict[str, list[Path]]] = {}
    for path in sorted(shard_dir.glob("*__*.fq.gz")):
        match = subshard_re.match(path.name)
        if not match:
            continue
        prefix = match.group("prefix")
        stream = match.group("stream")
        if stream not in streams:
            continue
        grouped.setdefault(prefix, {}).setdefault(stream, []).append(path)
    return grouped


def plan_demux_align_chunks_by_prefix(
    demux_dir: Path, prefixes: list[str]
) -> list[tuple[str, Path, Path, Path, Path]]:
    """Return expected regrouped demux FASTQ paths for prefix analysis chunks."""
    demux_dir = Path(demux_dir)
    return [
        (
            prefix,
            demux_dir / f"{prefix}.forward_1.fq.gz",
            demux_dir / f"{prefix}.forward_2.fq.gz",
            demux_dir / f"{prefix}.reverse_1.fq.gz",
            demux_dir / f"{prefix}.reverse_2.fq.gz",
        )
        for prefix in prefixes
    ]


def plan_fastp_shards(shard_dir: Path, number_of_split_parts: int) -> list[tuple[str, Path, Path]]:
    """Return expected shard paths when fastp outputs are not present yet."""
    shard_dir = Path(shard_dir)
    return [
        (
            chunk_id,
            shard_dir / f"{chunk_id}.R1.fq.gz",
            shard_dir / f"{chunk_id}.R2.fq.gz",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def plan_demux_align_chunks(
    demux_dir: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path, Path, Path]]:
    """Return expected demux FASTQ paths when demux outputs are not present yet."""
    demux_dir = Path(demux_dir)
    return [
        (
            chunk_id,
            demux_dir / f"{chunk_id}.forward_1.fq.gz",
            demux_dir / f"{chunk_id}.forward_2.fq.gz",
            demux_dir / f"{chunk_id}.reverse_1.fq.gz",
            demux_dir / f"{chunk_id}.reverse_2.fq.gz",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def plan_bismark_pe_bams_by_prefix(
    align_dir: Path, prefixes: list[str]
) -> list[tuple[str, Path, Path]]:
    """Return expected unsorted Bismark PE BAM paths for prefix chunks."""
    align_dir = Path(align_dir)
    return [
        (
            prefix,
            align_dir / f"{prefix}.forward_1_bismark_bt2_pe.bam",
            align_dir / f"{prefix}.reverse_1_bismark_bt2_pe.bam",
        )
        for prefix in prefixes
    ]


def plan_bismark_pe_bams(
    align_dir: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path]]:
    """Return expected unsorted Bismark PE BAM paths before alignment completes."""
    align_dir = Path(align_dir)
    return [
        (
            chunk_id,
            align_dir / f"{chunk_id}.forward_1_bismark_bt2_pe.bam",
            align_dir / f"{chunk_id}.reverse_1_bismark_bt2_pe.bam",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def plan_bismark_sortbyname_bams_by_prefix(
    align_dir: Path, prefixes: list[str]
) -> list[tuple[str, Path, Path]]:
    """Return expected sortbyname Bismark PE BAM paths for prefix chunks."""
    align_dir = Path(align_dir)
    return [
        (
            prefix,
            align_dir / f"{prefix}.forward_1_bismark_bt2_pe_sortbyname.bam",
            align_dir / f"{prefix}.reverse_1_bismark_bt2_pe_sortbyname.bam",
        )
        for prefix in prefixes
    ]


def plan_bismark_sortbyname_bams(
    align_dir: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path]]:
    """Return expected sortbyname Bismark PE BAM paths before bam_sort completes."""
    align_dir = Path(align_dir)
    return [
        (
            chunk_id,
            align_dir / f"{chunk_id}.forward_1_bismark_bt2_pe_sortbyname.bam",
            align_dir / f"{chunk_id}.reverse_1_bismark_bt2_pe_sortbyname.bam",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def plan_split_bam_chunk_pairs_by_prefix(
    split_root: Path, prefixes: list[str]
) -> list[tuple[str, Path, Path]]:
    """Return expected split BAM strand dirs for prefix chunks."""
    split_root = Path(split_root)
    return [
        (
            prefix,
            split_root / f"{prefix}.forward_1",
            split_root / f"{prefix}.reverse_1",
        )
        for prefix in prefixes
    ]


def plan_split_bam_chunk_pairs(
    split_root: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path]]:
    """Return expected split BAM strand dirs before split_bams completes."""
    split_root = Path(split_root)
    return [
        (
            chunk_id,
            split_root / f"{chunk_id}.forward_1",
            split_root / f"{chunk_id}.reverse_1",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def plan_merged_fr_bam_chunks_by_prefix(
    merged_root: Path, prefixes: list[str]
) -> list[tuple[str, Path, Path]]:
    """Return expected merged FR BAM dirs for prefix chunks."""
    merged_root = Path(merged_root)
    return [
        (
            prefix,
            merged_root / f"{prefix}_merged_fr_bam",
            merged_root / f"{prefix}_merge_filtered_barcode",
        )
        for prefix in prefixes
    ]


def plan_merged_fr_bam_chunks(
    merged_root: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path]]:
    """Return expected merged FR BAM dirs before merge_fr_bams completes."""
    merged_root = Path(merged_root)
    return [
        (
            chunk_id,
            merged_root / f"{chunk_id}_merged_fr_bam",
            merged_root / f"{chunk_id}_merge_filtered_barcode",
        )
        for chunk_id in build_chunk_names(number_of_split_parts)
    ]


def discover_demux_align_chunks(
    demux_dir: Path,
) -> list[tuple[str, Path, Path, Path, Path]]:
    """Return sorted (chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2) from demux outputs."""
    demux_dir = Path(demux_dir)
    if not demux_dir.is_dir():
        return []

    chunks: list[tuple[str, Path, Path, Path, Path]] = []
    for fwd_r1 in sorted(demux_dir.glob("*.forward_1.fq.gz")):
        chunk_id = fwd_r1.name[: -len(".forward_1.fq.gz")]
        fwd_r2 = demux_dir / f"{chunk_id}.forward_2.fq.gz"
        rev_r1 = demux_dir / f"{chunk_id}.reverse_1.fq.gz"
        rev_r2 = demux_dir / f"{chunk_id}.reverse_2.fq.gz"
        chunks.append((chunk_id, fwd_r1, fwd_r2, rev_r1, rev_r2))
    return chunks


def discover_bismark_pe_bams(
    align_dir: Path,
) -> list[tuple[str, Path, Path]]:
    """Return sorted (chunk_id, forward_bam, reverse_bam) from Bismark align outputs."""
    align_dir = Path(align_dir)
    if not align_dir.is_dir():
        return []

    chunks: list[tuple[str, Path, Path]] = []
    for fwd_bam in sorted(align_dir.glob("*.forward_1_bismark_bt2_pe.bam")):
        if "_sortbyname" in fwd_bam.name:
            continue
        chunk_id = fwd_bam.name[: -len(".forward_1_bismark_bt2_pe.bam")]
        rev_bam = align_dir / f"{chunk_id}.reverse_1_bismark_bt2_pe.bam"
        chunks.append((chunk_id, fwd_bam, rev_bam))
    return chunks


def discover_bismark_sortbyname_bams(
    align_dir: Path,
) -> list[tuple[str, Path, Path]]:
    """Return sorted (chunk_id, forward_sortbyname_bam, reverse_sortbyname_bam)."""
    align_dir = Path(align_dir)
    if not align_dir.is_dir():
        return []

    chunks: list[tuple[str, Path, Path]] = []
    for fwd_bam in sorted(align_dir.glob("*.forward_1_bismark_bt2_pe_sortbyname.bam")):
        chunk_id = fwd_bam.name[: -len(".forward_1_bismark_bt2_pe_sortbyname.bam")]
        rev_bam = align_dir / f"{chunk_id}.reverse_1_bismark_bt2_pe_sortbyname.bam"
        chunks.append((chunk_id, fwd_bam, rev_bam))
    return chunks


def counts_output_path(bam_path: Path) -> Path:
    """Map unsorted BAM to per-barcode counts CSV path."""
    return bam_path.parent / f"{bam_path.stem}_cb_aligned_reads_counts.csv"


def discover_split_bam_chunk_pairs(
    split_root: Path,
) -> list[tuple[str, Path, Path]]:
    """Return sorted (chunk_id, forward_dir, reverse_dir) from split_bams outputs."""
    split_root = Path(split_root)
    if not split_root.is_dir():
        return []

    pairs: list[tuple[str, Path, Path]] = []
    for forward_dir in sorted(split_root.glob("*.forward_1")):
        if not forward_dir.is_dir():
            continue
        chunk_id = forward_dir.name[: -len(".forward_1")]
        reverse_dir = split_root / f"{chunk_id}.reverse_1"
        pairs.append((chunk_id, forward_dir, reverse_dir))
    return pairs


def split_bams_strand_dir(split_root: Path, sortbyname_bam: Path) -> Path:
    """Directory name for one strand's split BAM outputs."""
    stem = sortbyname_bam.name[: -len(".bam")] if sortbyname_bam.name.endswith(".bam") else sortbyname_bam.name
    prefix = re.sub(r"_bismark_.*", "", stem)
    return split_root / prefix


def discover_merged_fr_bam_chunks(
    merged_root: Path,
) -> list[tuple[str, Path, Path]]:
    """Return sorted (chunk_id, bam_dir, filtered_barcode_path) from merge_fr_bams."""
    merged_root = Path(merged_root)
    if not merged_root.is_dir():
        return []

    chunks: list[tuple[str, Path, Path]] = []
    for bam_dir in sorted(merged_root.glob("*_merged_fr_bam")):
        if not bam_dir.is_dir():
            continue
        chunk_id = bam_dir.name[: -len("_merged_fr_bam")]
        filtered_barcode = merged_root / f"{chunk_id}_merge_filtered_barcode"
        chunks.append((chunk_id, bam_dir, filtered_barcode))
    return chunks


def resolve_barcode_mode(settings: dict) -> str:
    """Return 'gexcb' or 'expected_cell_num' (mutually exclusive)."""
    gexcb = settings.get("gexcb")
    expected_cell_num = settings.get("expected_cell_num")
    has_gexcb = gexcb not in (None, "")
    has_expected = expected_cell_num is not None
    if has_gexcb and has_expected:
        raise ValueError("gexcb and expected_cell_num are mutually exclusive")
    if has_gexcb:
        return "gexcb"
    return "expected_cell_num"


def require_bismark_ref(path: Path) -> None:
    ref = Path(path)
    if not ref.is_dir():
        raise ValueError(f"missing bismark_ref directory: {ref}")
    bisulfite = ref / "Bisulfite_Genome"
    for subdir in ("CT_conversion", "GA_conversion"):
        candidate = bisulfite / subdir
        if not candidate.is_dir():
            raise ValueError(
                f"invalid bismark_ref (missing {candidate}); "
                "pass the parent of Bisulfite_Genome/"
            )
