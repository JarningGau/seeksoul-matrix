"""Shared path existence checks for workflow drivers (make_cmd entrypoints)."""

from __future__ import annotations

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
