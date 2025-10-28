#!/usr/bin/env bash
# Snapshot the current env spec (human-readable)
set -euxo pipefail
ENV_NAME="${1:-ovmpk}"
mkdir -p runlogs
mamba env export -n "$ENV_NAME" --no-builds > runlogs/env.export.yml || conda env export -n "$ENV_NAME" --no-builds > runlogs/env.export.yml
echo "Wrote runlogs/env.export.yml"
