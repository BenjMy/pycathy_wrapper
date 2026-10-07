#!/usr/bin/env bash
# Build the pyCATHY Sphinx docs locally and publish them to the gh-pages branch.
#
# Usage (from the repository root):
#   ./publish_docs.sh            # build, check, publish
#   ./publish_docs.sh --no-push  # build and check only, do not publish
#
# The script activates the micromamba environment ENV_NAME and saves everything
# printed to the terminal in logs/publish_docs_<timestamp>.log

set -euo pipefail

ENV_NAME="pycathy_doc"
DOC_DIR="doc"
BUILD_DIR="${DOC_DIR}/_build"
REMOTE="origin"
BRANCH="gh-pages"
LOG_DIR="logs"
PUSH=true

if [[ "${1:-}" == "--no-push" ]]; then
    PUSH=false
fi

# Must be run from the repository root
if [[ ! -d "${DOC_DIR}" || ! -d ".git" ]]; then
    echo "Error: run this script from the repository root (where '${DOC_DIR}/' and '.git/' are)." >&2
    exit 1
fi

# ---------------------------------------------------------------- logging
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/publish_docs_$(date +%Y%m%d_%H%M%S).log"
# Send stdout and stderr to the terminal AND to the log file
exec > >(tee -a "${LOG_FILE}") 2>&1
trap 'status=$?; echo; if [[ ${status} -ne 0 ]]; then echo "==> FAILED (exit code ${status}). See log: ${LOG_FILE}"; else echo "==> Log saved to: ${LOG_FILE}"; fi' EXIT

echo "==> $(date)"
echo "==> Log file: ${LOG_FILE}"

# ------------------------------------------------- micromamba environment
if ! command -v micromamba >/dev/null 2>&1; then
    echo "Error: micromamba not found in PATH." >&2
    exit 1
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "${ENV_NAME}" ]]; then
    echo "==> Activating micromamba environment '${ENV_NAME}'"
    # Activation scripts can reference unset variables, so relax 'set -u' here
    set +u
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "${ENV_NAME}"
    set -u
else
    echo "==> Environment '${ENV_NAME}' already active"
fi

echo "==> Python: $(which python) ($(python --version 2>&1))"

# ------------------------------------------------------------ tools check
if ! python -c "import sphinx" 2>/dev/null; then
    echo "Error: sphinx is not installed in environment '${ENV_NAME}'." >&2
    exit 1
fi

if ! python -c "import ghp_import" 2>/dev/null; then
    echo "==> Installing ghp-import..."
    python -m pip install ghp-import
fi

# Warn about uncommitted changes (the docs may not match what is on GitHub)
if [[ -n "$(git status --porcelain)" ]]; then
    echo "Warning: you have uncommitted changes. The docs will be built from your working tree."
    printf "Continue? [y/N] "
    read -r answer
    [[ "${answer}" =~ ^[Yy]$ ]] || exit 1
fi

# ------------------------------------------------------------------ build
echo "==> Cleaning ${BUILD_DIR}"
rm -rf "${BUILD_DIR}"

echo "==> Building docs with Sphinx"
python -m sphinx -b html "${DOC_DIR}" "${BUILD_DIR}"

# Check the output before publishing anything
if [[ ! -f "${BUILD_DIR}/index.html" ]]; then
    echo "Error: ${BUILD_DIR}/index.html was not created. Not publishing." >&2
    exit 1
fi
echo "==> Build OK (${BUILD_DIR}/index.html found)"

if [[ "${PUSH}" == false ]]; then
    echo "==> --no-push given: skipping publish. Open ${BUILD_DIR}/index.html to preview."
    exit 0
fi

# ---------------------------------------------------------------- publish
# -n : add .nojekyll so the _static folder is served
# -p : push to the remote
# -f : force push, so the branch contains exactly this build
echo "==> Publishing to ${REMOTE}/${BRANCH}"
python -m ghp_import -n -p -f -r "${REMOTE}" -b "${BRANCH}" "${BUILD_DIR}"

echo "==> Done. The site should update in a minute or two at:"
echo "    https://benjmy.github.io/pycathy_wrapper/"
