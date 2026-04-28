#!/usr/bin/env bash
# Compile paper.md to paper.pdf.
#
# Two paths are supported, in order of preference:
#   1. JOSS Inara container (Docker or Podman): produces the exact
#      JOSS-styled PDF that reviewers will see.
#   2. Pandoc fallback: produces a generic PDF for local preview.
#      Requires pandoc + a LaTeX engine (xelatex preferred).
#
# Run from the source/ directory (the directory containing paper.md).

set -euo pipefail

cd "$(dirname "$0")"

if [ ! -f paper.md ] || [ ! -f paper.bib ]; then
    echo "error: paper.md and paper.bib must be present in $(pwd)" >&2
    exit 1
fi

run_inara() {
    local engine="$1"
    echo "[build] using ${engine} + openjournals/inara (JOSS-styled output)"
    local extra=()
    if [ "${engine}" = "podman" ]; then
        # Rootless podman idmaps --user through subuid, which prevents
        # writing to the host-owned bind mount. keep-id maps the
        # container UID 1:1 to the host UID, so output files land in
        # the source/ folder owned by the invoking user.
        extra+=(--userns=keep-id)
    fi
    "${engine}" run --rm \
        --volume "${PWD}:/data" \
        --user "$(id -u):$(id -g)" \
        "${extra[@]}" \
        --env JOURNAL=joss \
        docker.io/openjournals/inara paper.md
}

run_pandoc() {
    echo "[build] Docker/Podman not found; falling back to pandoc (generic PDF)"
    local engine
    if command -v xelatex >/dev/null 2>&1; then
        engine=xelatex
    elif command -v pdflatex >/dev/null 2>&1; then
        engine=pdflatex
    else
        echo "error: need xelatex or pdflatex installed for pandoc fallback" >&2
        exit 1
    fi
    pandoc paper.md \
        --output paper.pdf \
        --citeproc \
        --bibliography paper.bib \
        --pdf-engine="${engine}" \
        --metadata link-citations=true
}

if command -v docker >/dev/null 2>&1; then
    run_inara docker
elif command -v podman >/dev/null 2>&1; then
    run_inara podman
else
    run_pandoc
fi

echo "[build] wrote $(pwd)/paper.pdf"
