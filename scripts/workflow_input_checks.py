"""Shared path existence checks for workflow drivers (make_cmd entrypoints)."""

from __future__ import annotations

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


def plan_fastp_shards(shard_dir: Path, number_of_split_parts: int) -> list[tuple[str, Path, Path]]:
    """Return expected shard paths when fastp outputs are not present yet."""
    if number_of_split_parts <= 0:
        raise ValueError("number_of_split_parts must be > 0")
    shard_dir = Path(shard_dir)
    shards: list[tuple[str, Path, Path]] = []
    for index in range(number_of_split_parts):
        chunk_id = f"{index + 1:04d}"
        if index == 0:
            r1_name, r2_name = "R1.fq.gz", "R2.fq.gz"
        else:
            r1_name = f"R1_{index:03d}.fq.gz"
            r2_name = f"R2_{index:03d}.fq.gz"
        shards.append((chunk_id, shard_dir / r1_name, shard_dir / r2_name))
    return shards


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
