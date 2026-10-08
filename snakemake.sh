#!/usr/bin/env bash
set -euo pipefail
WORKFLOW_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MODE="${1:-dry}"
if (($#)); then shift; fi
case "$MODE" in
  dry|dag|run|rerun|unlock|pilot|controller) ;;
  *) printf 'Unknown mode: %s\n' "$MODE" >&2; exit 2 ;;
esac
# Additive infrastructure modules only; biological tools remain the user's PATH.
if ! command -v snakemake >/dev/null 2>&1 || [[ "$(snakemake --version)" != 9.4.0 ]]; then
  if type module >/dev/null 2>&1; then
    module load bioinfo-ifb
    module load snakemake/9.4.0
  fi
fi
if ! command -v snakemake >/dev/null 2>&1; then
  printf 'Required Snakemake 9.4.0 is unavailable; preload it or expose the infrastructure modules.\n' >&2
  exit 1
fi
exec python3 "$WORKFLOW_ROOT/profile/workflow.py" "$MODE" "$@"
