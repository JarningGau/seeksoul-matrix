#!/usr/bin/env bash
# Install seekgene/Bismark Perl scripts into the active pixi/conda environment.
# Mirrors SeekSoulMethyl: https://github.com/seekgene/SeekSoulMethyl README install steps.
set -euo pipefail

SEEKGENE_BISMARK_COMMIT="${SEEKGENE_BISMARK_COMMIT:-363ea7ab5b3aa329e363b6de1778f52e9004a9e0}"
SEEKGENE_BISMARK_REPO="${SEEKGENE_BISMARK_REPO:-https://github.com/seekgene/Bismark.git}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "error: CONDA_PREFIX is not set; run via 'pixi run setup-bismark'" >&2
    exit 1
fi

dest="${CONDA_PREFIX}/bin"
tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

echo "[setup-bismark] cloning ${SEEKGENE_BISMARK_REPO} @ ${SEEKGENE_BISMARK_COMMIT}"
git clone --quiet "${SEEKGENE_BISMARK_REPO}" "${tmpdir}/Bismark"
git -C "${tmpdir}/Bismark" checkout --quiet "${SEEKGENE_BISMARK_COMMIT}"

echo "[setup-bismark] installing scripts into ${dest}"
cp -r "${tmpdir}/Bismark"/* "${dest}/"
chmod +x "${dest}"/bismark* "${dest}"/deduplicate_bismark

if ! grep -q -- '--add_barcode' "${dest}/bismark"; then
    echo "error: installed bismark does not expose --add_barcode (wrong fork?)" >&2
    exit 1
fi

echo "[setup-bismark] ok: $(command -v bismark)"
