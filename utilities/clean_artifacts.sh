#!/usr/bin/env bash
# Remove every run directory under artifacts/, keeping the directory itself.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
artifacts="${CGLBM_ARTIFACTS_DIR:-${root}/artifacts}"

read -rp "Delete every run under ${artifacts}? [y/N] " answer
if [[ "${answer}" =~ ^[Yy]$ ]]; then
    find "${artifacts}" -mindepth 1 -not -name '.gitkeep' -delete
    echo "Cleaned ${artifacts}"
else
    echo "Nothing removed."
fi
