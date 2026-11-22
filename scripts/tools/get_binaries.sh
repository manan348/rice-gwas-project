#!/usr/bin/env bash
# Download the pinned PLINK 1.9 and GEMMA binaries used by the pipeline into ./bin
# (Linux x86_64 only — run inside Linux, WSL, or Google Colab).
set -euo pipefail

BIN_DIR="$(cd "$(dirname "$0")/../.." && pwd)/bin"
mkdir -p "$BIN_DIR"
cd "$BIN_DIR"

PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"   # PLINK v1.90b7.2
GEMMA_URL="https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz"

if [ ! -x plink ]; then
  wget -q "$PLINK_URL" -O plink.zip
  unzip -o -q plink.zip plink
  rm -f plink.zip
  chmod +x plink
fi

if [ ! -x gemma ]; then
  wget -q "$GEMMA_URL" -O gemma.gz
  gunzip -f gemma.gz
  chmod +x gemma
fi

./plink --version
./gemma -h 2>&1 | head -n 4 | sed -n '/GEMMA/p'
