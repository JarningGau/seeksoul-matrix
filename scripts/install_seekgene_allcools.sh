#!/usr/bin/env bash
# Install seekgene/ALLCools into the active pixi/conda environment.
set -euo pipefail

SEEKGENE_ALLCOOLS_COMMIT="${SEEKGENE_ALLCOOLS_COMMIT:-b84c180752c7bcc090994efc8534852d497f59d2}"
SEEKGENE_ALLCOOLS_REPO="${SEEKGENE_ALLCOOLS_REPO:-https://github.com/seekgene/ALLCools.git}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "error: CONDA_PREFIX is not set; run via 'pixi run setup-allcools'" >&2
    exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

echo "[setup-allcools] cloning ${SEEKGENE_ALLCOOLS_REPO} @ ${SEEKGENE_ALLCOOLS_COMMIT}"
git clone --quiet "${SEEKGENE_ALLCOOLS_REPO}" "${tmpdir}/ALLCools"
git -C "${tmpdir}/ALLCools" checkout --quiet "${SEEKGENE_ALLCOOLS_COMMIT}"

echo "[setup-allcools] pip install into ${CONDA_PREFIX}"
if [[ -x "${CONDA_PREFIX}/bin/pip" ]]; then
    "${CONDA_PREFIX}/bin/pip" install --quiet "${tmpdir}/ALLCools"
else
    "${CONDA_PREFIX}/bin/python" -m pip install --quiet "${tmpdir}/ALLCools"
fi

if ! command -v allcools >/dev/null 2>&1; then
    echo "error: allcools not found on PATH after install" >&2
    exit 1
fi

echo "[setup-allcools] ok: $(command -v allcools)"
allcools --version
