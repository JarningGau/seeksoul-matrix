#!/usr/bin/env bash
set -euo pipefail

require_cmd() {
    local name="$1"
    if ! command -v "${name}" >/dev/null 2>&1; then
        echo "error: ${name} not found on PATH" >&2
        exit 1
    fi
    echo "[check-allcools-env] ${name}: $(command -v "${name}")"
}

require_cmd allcools
require_cmd samtools

version="$(allcools --version 2>&1 || true)"
if [[ -z "${version}" ]]; then
    echo "error: allcools --version produced no output" >&2
    echo "hint: run 'pixi run setup-allcools'" >&2
    exit 1
fi

echo "[check-allcools-env] allcools version: ${version}"
echo "[check-allcools-env] ok"
