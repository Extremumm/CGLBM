#!/usr/bin/env bash
# Run one CGLBM program inside a fresh directory under artifacts/.
#
#   utilities/run_case.sh laplace [variant] [rundir]
#
# The programs write their CSV output to the current working directory with
# fixed names, so every run needs a directory of its own.
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <program> [opt|dbg] [rundir]" >&2
    exit 1
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
program="$1"
variant="${2:-opt}"
rundir="${3:-${root}/artifacts/${program}}"
bindir="${CGLBM_BIN_DIR:-${root}/bin}"

executable="$(find "${bindir}" -type f -name "${program}_${variant}" -print -quit)"
if [[ -z "${executable}" ]]; then
    echo "No executable ${program}_${variant} under ${bindir}." >&2
    echo "Build it first: cmake --build --preset gnu --target ${program}_${variant}" >&2
    exit 1
fi

mkdir -p "${rundir}"
echo "Running ${executable}"
echo "     in ${rundir}"
cd "${rundir}"
exec "${executable}" 2>&1 | tee run.log
