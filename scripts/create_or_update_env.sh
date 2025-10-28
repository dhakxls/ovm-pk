#!/usr/bin/env bash
# Usage: scripts/create_or_update_env.sh [env_name]   (default: ovmpk)
set -euxo pipefail
ENV_NAME="${1:-ovmpk}"

# prefer mamba if available
if command -v mamba >/dev/null 2>&1; then
  SHELL_HOOK='eval "$(mamba shell hook)"'
  PM=mamba
else
  SHELL_HOOK='eval "$(conda shell.bash hook)"'
  PM=conda
fi
bash -lc "$SHELL_HOOK; $PM --version"

if ls conda-lock-*.yml >/dev/null 2>&1; then
  LOCKFILE="$(ls conda-lock-*.yml | head -n1)"
  echo "Using lock file: $LOCKFILE"
  if command -v micromamba >/dev/null 2>&1; then
    micromamba create -y -n "$ENV_NAME" -f "$LOCKFILE"
  else
    $PM create -y -n "$ENV_NAME" -f "$LOCKFILE"
  fi
elif [ -f environment.yml ]; then
  # create if missing, otherwise update
  if $PM env list | grep -qE "^$ENV_NAME\s"; then
    $PM env update -n "$ENV_NAME" -f environment.yml
  else
    $PM env create -n "$ENV_NAME" -f environment.yml
  fi
else
  echo "No environment.yml or lock file found." >&2
  exit 1
fi

# make repo importable
bash -lc "$SHELL_HOOK; $PM activate $ENV_NAME; pip install -e ."
