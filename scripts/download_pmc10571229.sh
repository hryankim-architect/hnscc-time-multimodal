#!/usr/bin/env bash
# Download the PMC10571229 (DeepLIIF / Ghahremani 2023 MICCAI) HNSCC IHC dataset.
# Run in the background; logs to /tmp/p4_data_pmc.log.
#
# Usage:
#   nohup bash scripts/download_pmc10571229.sh >/tmp/p4_data_pmc.log 2>&1 &
#
# This script is intentionally idempotent: if the destination already exists
# and looks complete, it skips. Re-running on top of a partial download is safe.

set -uo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DEST="$REPO_ROOT/data/pmc10571229"
LOG_PREFIX="[p4_pmc_download $(date -u +%Y-%m-%dT%H:%M:%SZ)]"

mkdir -p "$DEST"
echo "$LOG_PREFIX starting; dest=$DEST"

# ---- Step 1: shallow-clone the source repo (small, ~50 MB) -----------------
if [[ ! -d "$DEST/source-repo/.git" ]]; then
  echo "$LOG_PREFIX cloning DeepLIIF source repo"
  git clone --depth 1 https://github.com/nadeemlab/DeepLIIF "$DEST/source-repo" \
    || { echo "$LOG_PREFIX clone failed"; exit 2; }
else
  echo "$LOG_PREFIX source-repo already present, skipping clone"
fi

# ---- Step 2: fetch the ROI image set (RESOLVED 2026-06-07: TCIA via pathdb) -
# The canonical home of the PMC10571229 HNSCC mIF/mIHC dataset is The Cancer
# Imaging Archive, NOT Zenodo (the earlier Zenodo record 8367318 was wrong).
#
#   Collection : HNSCC-mIF-mIHC-comparison (TCIA, Version 2, 2023-08-31)
#   DOI        : 10.7937/TCIA.2020.T90F-WB82  |  CC BY 4.0  |  8 patients
#   Contents   : 3212 PNG ROIs (~1.01 GB); naming Case[id]_[T/M/S][1-3]_[idx]_[marker]
#
# TCIA's web bundle is faspex/Aspera, but every ROI is *also* directly HTTPS-
# fetchable from the TCIA pathology host (pathdb). data/pmc10571229/rois_manifest.tsv
# (committed) pins all 3212: rel_path <TAB> sha256 <TAB> url. We download each by
# URL and verify its sha256, so no Aspera client is needed.
LEDGER="$REPO_ROOT/data/pmc10571229/rois_manifest.tsv"
if [[ ! -f "$LEDGER" ]]; then
  echo "$LOG_PREFIX ERROR: $LEDGER missing (it is committed; check out the repo)."; exit 2
fi
mkdir -p "$DEST/rois"
n=0; ok=0; bad=0
while IFS=$'\t' read -r REL SHA URL; do
  [[ "$REL" == \#* || -z "$REL" ]] && continue
  n=$((n+1))
  OUT="$DEST/rois/$REL"; mkdir -p "$(dirname "$OUT")"
  if [[ -f "$OUT" ]] && [[ "$(shasum -a 256 "$OUT" | cut -d' ' -f1)" == "$SHA" ]]; then
    ok=$((ok+1)); continue
  fi
  curl -fsSL --retry 3 -o "$OUT" "${URL/http:\/\//https://}" || { echo "$LOG_PREFIX WARN download failed: $REL"; bad=$((bad+1)); continue; }
  got="$(shasum -a 256 "$OUT" | cut -d' ' -f1)"
  if [[ "$got" == "$SHA" ]]; then ok=$((ok+1)); else echo "$LOG_PREFIX SHA MISMATCH: $REL"; bad=$((bad+1)); fi
  [[ $((n % 200)) -eq 0 ]] && echo "$LOG_PREFIX ...$n verified=$ok"
done < "$LEDGER"
echo "$LOG_PREFIX ROI fetch done: $ok/$n verified, $bad failed"

# ---- Step 3: integrity check of the committed ledger itself -----------------
# (The ihc.checksums_file_sha256 in data/manifest.yaml should match this.)
echo "$LOG_PREFIX ledger sha256: $(shasum -a 256 "$LEDGER" | cut -d' ' -f1)"

echo "$LOG_PREFIX done"
