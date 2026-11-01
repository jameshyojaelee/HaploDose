#!/usr/bin/env bash
# Build the Rust extension into the haplodose env and run Rust + Python tests.
# Run on a compute node. Usage: scripts/build_and_test.sh [pytest args]
set -euo pipefail
cd "$(dirname "$0")/.."
ENV=/gpfs/commons/home/jameslee/micromamba/envs/haplodose
PY=$ENV/bin/python
export PYTHONNOUSERSITE=1 JAX_PLATFORMS=${JAX_PLATFORMS:-cpu} XLA_PYTHON_CLIENT_PREALLOCATE=false
export PATH="$HOME/.cargo/bin:$ENV/bin:$PATH" PYO3_PYTHON=$PY
export CARGO_NET_OFFLINE=true
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
nproc
rustc --version
(cd crates/core && cargo test --release --offline)
"$PY" -s -m pip install --no-user --no-build-isolation --no-deps -e . 2>&1 | tail -2
"$PY" -s -c "import haplodose._native as n; print('native', n.native_available())"
"$PY" -s -m pytest -q "$@"
"$PY" -s -m pip freeze > "${FREEZE_OUT:-/dev/null}"
