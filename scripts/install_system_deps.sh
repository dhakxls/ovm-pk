#!/usr/bin/env bash
set -euxo pipefail
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential git cmake ninja-build pkg-config \
  wget curl unzip ca-certificates \
  libgl1 libegl1 libx11-6 libxi6 libxext6 libxrender1 libxtst6 \
  mesa-utils
# Keep CUDA drivers/toolkit in apt (if you use them) and OUT of conda.
