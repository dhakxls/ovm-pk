# source this after conda/mamba activation
export GNINA_BUILD="${HOME}/src/gnina/build"
# Prefer GNINA’s vendored libs (libtorch + extras)
if [ -d "${GNINA_BUILD}" ]; then
  export LD_LIBRARY_PATH="${GNINA_BUILD}/_deps/libtorch-src/lib:${GNINA_BUILD}/external/lib:${LD_LIBRARY_PATH:-}"
  # Models (fallback to repo location if share dir not present)
  if [ -d "${GNINA_BUILD}/share/gnina/models" ]; then
    export GNINA_MODEL_DIR="${GNINA_BUILD}/share/gnina/models"
  fi
fi
