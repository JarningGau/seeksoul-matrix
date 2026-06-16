#!/usr/bin/env bash
set -euo pipefail

require_cmd() {
    local name="$1"
    if ! command -v "${name}" >/dev/null 2>&1; then
        echo "error: ${name} not found on PATH" >&2
        exit 1
    fi
    echo "[check-bismark-env] ${name}: $(command -v "${name}")"
}

require_cmd bismark
require_cmd bowtie2
require_cmd samtools

if ! bismark --help 2>&1 | grep -q -- '--add_barcode'; then
    echo "error: bismark on PATH is not seekgene fork (missing --add_barcode)" >&2
    echo "hint: run 'pixi run setup-bismark'" >&2
    exit 1
fi

if ! bismark --help 2>&1 | grep -q -- '--add_umi'; then
    echo "error: bismark on PATH is not seekgene fork (missing --add_umi)" >&2
    exit 1
fi

echo "[check-bismark-env] seekgene Bismark flags present (--add_barcode, --add_umi)"
echo "[check-bismark-env] ok"
