#!/usr/bin/env bash
# Run AFTER: mamba/conda activate ovmpk
set -euo pipefail
: "${CONDA_PREFIX:?Activate your env first (e.g., 'mamba activate ovmpk')}"

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"

cat > "$CONDA_PREFIX/etc/conda/activate.d/ovmpk.sh" <<'EOS'
# GNINA paths (adjust GNINA_HOME if you installed elsewhere)
export GNINA_HOME="${GNINA_HOME:-$HOME/src/gnina/build}"
export PATH="$GNINA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GNINA_HOME/_deps/libtorch-src/lib:$GNINA_HOME/external/lib:${LD_LIBRARY_PATH:-}"
# CNN models
if [ -d "$GNINA_HOME/share/gnina/models" ]; then
  export GNINA_MODEL_DIR="$GNINA_HOME/share/gnina/models"
else
  export GNINA_MODEL_DIR="$HOME/src/gnina/gninasrc/lib/models"
fi
# Default to CPU unless you set this yourself
: "${CUDA_VISIBLE_DEVICES:=}"; export CUDA_VISIBLE_DEVICES
EOS

cat > "$CONDA_PREFIX/etc/conda/deactivate.d/ovmpk.sh" <<'EOS'
unset GNINA_HOME GNINA_MODEL_DIR CUDA_VISIBLE_DEVICES
EOS

echo "Activation hooks installed under: $CONDA_PREFIX/etc/conda/{activate.d,deactivate.d}/"
