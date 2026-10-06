#!/bin/bash
# Evaluation-only package reuse with the exact frozen CPython patch release.
set -euo pipefail
TREDNET_REPO=/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow
TREDNET_SITE_PACKAGES="$TREDNET_REPO/models/TREDNET_v2/.venv/lib/python3.13/site-packages"
TREDNET_PYTHON=/home/maheshwarany2/.local/share/uv/python/cpython-3.13.0-linux-x86_64-gnu/bin/python3.13
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="$TREDNET_SITE_PACKAGES"
export TF_DETERMINISTIC_OPS=1
export TF_CUDNN_DETERMINISTIC=1
export TF_XLA_FLAGS=--tf_xla_enable_xla_devices=false
export OMP_NUM_THREADS=8
export TF_NUM_INTRAOP_THREADS=8
export TF_NUM_INTEROP_THREADS=2
for TREDNET_LIBRARY_DIR in "$TREDNET_SITE_PACKAGES"/nvidia/*/lib; do
    if [[ -d "$TREDNET_LIBRARY_DIR" ]]; then
        export LD_LIBRARY_PATH="$TREDNET_LIBRARY_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
done
exec "$TREDNET_PYTHON" "$@"
