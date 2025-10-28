#!/usr/bin/env bash
# One-shot: system deps → conda env → hooks
set -euxo pipefail
./scripts/install_system_deps.sh
./scripts/create_or_update_env.sh ovmpk

# activate to install hooks
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook)"; mamba activate ovmpk
else
  eval "$(conda shell.bash hook)"; conda activate ovmpk
fi
./scripts/install_hooks.sh

echo "Bootstrap complete. Next: 'mamba activate ovmpk' and run 'pytest -q'"
