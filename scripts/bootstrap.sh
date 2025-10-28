#!/usr/bin/env bash
set -euo pipefail

# -------- helpers --------
need() { command -v "$1" >/dev/null 2>&1; }
log()  { printf "\n\033[1;34m[ovmpk]\033[0m %s\n" "$*"; }
err()  { printf "\n\033[1;31m[ovmpk ERROR]\033[0m %s\n" "$*" >&2; }

OS="$(uname -s)"
ARCH="$(uname -m)"
IS_WSL="0"
if [ "$OS" = "Linux" ] && grep -qi microsoft /proc/version 2>/dev/null; then IS_WSL="1"; fi

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
SCRIPTS_DIR="$REPO_ROOT/scripts"
MINIFORGE_DIR="$HOME/miniforge3"

mkdir -p "$SCRIPTS_DIR"

# -------- install system deps --------
install_linux_deps() {
  log "Installing Linux/WSL deps (sudo apt)…"
  sudo apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential git curl wget cmake ninja-build pkg-config \
    mesa-utils libglu1-mesa-dev \
    openbabel pymol || true   # PyMOL may not exist on every Ubuntu; ignore failure
}

install_macos_deps() {
  log "Installing macOS deps (Homebrew)…"
  if ! need brew; then
    log "Homebrew not found — installing Homebrew (will prompt if needed)…"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    if [ -x /opt/homebrew/bin/brew ]; then
      eval "$(/opt/homebrew/bin/brew shellenv)"
    else
      eval "$(/usr/local/bin/brew shellenv)"
    fi
  fi
  brew update
  brew install git curl wget cmake ninja pkg-config open-babel || true
  # Optional GUI PyMOL:
  # brew install --cask pymol || true
}

# -------- install miniforge/mamba --------
install_miniforge() {
  if need mamba; then
    log "mamba found — skipping Miniforge install."
    return 0
  fi
  log "Installing Miniforge (mamba)…"
  case "$OS-$ARCH" in
    Linux-x86_64)  MF=Miniforge3-Linux-x86_64.sh ;;
    Linux-aarch64) MF=Miniforge3-Linux-aarch64.sh ;;
    Darwin-arm64)  MF=Miniforge3-MacOSX-arm64.sh ;;
    Darwin-x86_64) MF=Miniforge3-MacOSX-x86_64.sh ;;
    *) err "Unsupported platform: $OS $ARCH"; exit 1 ;;
  esac
  curl -fsSL -o /tmp/$MF "https://github.com/conda-forge/miniforge/releases/latest/download/$MF"
  bash /tmp/$MF -b -p "$MINIFORGE_DIR"
  eval "$("$MINIFORGE_DIR/bin/conda" shell.bash hook)"
  "$MINIFORGE_DIR/bin/conda" config --set auto_activate_base false
}

activate_conda() {
  eval "$("$MINIFORGE_DIR/bin/conda" shell.bash hook)"
  if ! need mamba; then
    log "Installing mamba…"
    conda install -y -n base -c conda-forge mamba
  fi
}

# -------- conda env --------
create_update_env() {
  log "Creating/updating conda env ovmpk from environment.yml…"
  if conda env list | grep -q "^ovmpk\s"; then
    mamba env update -n ovmpk -f "$REPO_ROOT/environment.yml"
  else
    mamba env create -f "$REPO_ROOT/environment.yml"
  fi
}

# -------- GNINA build --------
build_gnina() {
  local SRC="${HOME}/src/gnina"
  local USE_CUDA_FLAG="-DUSE_CUDA=OFF"
  if [ "$OS" = "Linux" ] && need nvidia-smi; then
    log "NVIDIA GPU detected — building GNINA with CUDA"
    USE_CUDA_FLAG="-DUSE_CUDA=ON"
  else
    log "CUDA not detected — building GNINA CPU-only"
  fi

  if [ ! -d "$SRC" ]; then
    log "Cloning GNINA…"
    git clone --recursive https://github.com/gnina/gnina.git "$SRC"
  else
    log "Updating GNINA…"
    (cd "$SRC" && git fetch -q && git pull -q && git submodule update --init --recursive)
  fi

  log "Configuring & building GNINA…"
  cmake -S "$SRC" -B "$SRC/build" -DCMAKE_BUILD_TYPE=Release $USE_CUDA_FLAG
  cmake --build "$SRC/build" -j
}

# -------- write activation stubs --------
write_activation() {
  log "Writing activation helper scripts…"

  # GNINA env snippet
  cat > "$SCRIPTS_DIR/activate_gnina.sh" << 'EOF'
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
EOF

  # one-shot activator for everything
  cat > "$SCRIPTS_DIR/activate_all.sh" << 'EOF'
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
EOF

  # small PyMOL view macro (optional)
  cat > "$REPO_ROOT/view.pml" << 'EOF'
# Quick view: receptor surface + top docked pose (if present)
load data/work/docking/5VCC_fixed_ph7.4_cleaned_fixed.pdbqt, rec
# Prefer selected pose if available
load data/output/pose_selection_test/selected_pose/ligand_best_pose_1_from_5VCC_fixed_ph7.4_cleaned_fixed_ketoconazole_prepared_ph7.4_s5573673_poses.sdf, lig
hide everything
as surface, rec
set transparency, 0.25, rec
split_states lig
disable lig
as sticks, lig_0001
util.cbaw lig_0001
orient lig_0001
select pocket, byres (rec within 4 of lig_0001)
show sticks, pocket
color tv_yellow, pocket
EOF

  chmod +x "$SCRIPTS_DIR/activate_all.sh" || true
}

# -------- main --------
main() {
  case "$OS" in
    Linux) install_linux_deps ;;
    Darwin) install_macos_deps ;;
    *) err "Unsupported OS: $OS"; exit 1 ;;
  esac

  install_miniforge
  activate_conda
  create_update_env
  build_gnina
  write_activation

  log "Done. Next steps:"
  echo "  source scripts/activate_all.sh"
  echo "  pytest -q                    # run tests"
  echo "  pymol -cq view.pml           # quick visualization (if PyMOL available)"
}

main "$@"
