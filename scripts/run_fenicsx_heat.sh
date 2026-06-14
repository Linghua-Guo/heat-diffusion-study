#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FENICSX_PYTHON="${FENICSX_PYTHON:-/opt/anaconda3/envs/fenicsx/bin/python}"

export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/private/tmp/pde_fenics_cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/private/tmp/pde_mpl}"
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"

unset CC CXX CFLAGS LDFLAGS
export PATH="$(dirname "$FENICSX_PYTHON"):/opt/anaconda3/bin:/usr/bin:/bin:/usr/sbin:/sbin"

exec "$FENICSX_PYTHON" "$PROJECT_DIR/scripts/solve_fenicsx_heat.py" "$@"
