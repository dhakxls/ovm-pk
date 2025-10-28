#!/usr/bin/env bash
# Generate a lock file for fully pinned, reproducible installs
set -euxo pipefail
if ! command -v conda-lock >/dev/null 2>&1; then
  mamba install -y -n base conda-lock -c conda-forge || conda install -y -n base conda-lock -c conda-forge
fi
conda-lock -f environment.yml -p linux-64
echo "Created lock: conda-lock-linux-64.yml"
