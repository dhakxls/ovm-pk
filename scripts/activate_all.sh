# shellcheck shell=bash
set -e
MINIFORGE_DIR="$HOME/miniforge3"
if ! command -v conda >/dev/null 2>&1; then
  eval "$("$MINIFORGE_DIR/bin/conda" shell.bash hook)"
else
  eval "$(conda shell.bash hook)"
fi
conda activate ovmpk
# add GNINA bindings
source "$(dirname "${BASH_SOURCE[0]}")/activate_gnina.sh"
echo "[ovmpk] environment ready. Use: pytest -q, or run scripts."
