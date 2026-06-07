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

# ---- Step 2: the ROI image set (RESOLVED 2026-06-07: TCIA, not Zenodo) ------
# The canonical home of the PMC10571229 HNSCC mIF/mIHC dataset is The Cancer
# Imaging Archive, NOT Zenodo (the earlier Zenodo record 8367318 was wrong --
# it yielded a 7-KB stray .py):
#
#   Collection : HNSCC-mIF-mIHC-comparison (TCIA, Version 2, 2023-08-31)
#   DOI        : 10.7937/TCIA.2020.T90F-WB82
#   Landing    : https://www.cancerimagingarchive.net/collection/hnscc-mif-mihc-comparison/
#   Contents   : 8 patients, 3,216 PNG images (~1.01 GB), CC BY 4.0
#   Naming      : Case[id]_[T/M/S][1-3]_[ROI_index]_[marker]
#
# TCIA distributes this package via faspex/Aspera (IBM Aspera Connect), not a
# plain HTTP(S) URL, so it cannot be curl'd unattended here. Download it once
# with Aspera Connect from the landing page above and unpack into:
#
#   data/pmc10571229/rois/
#
# Then run Step 3 below to record sha256 sums. (Per-file checksums are therefore
# NOT pinned in data/manifest.yaml; the ihc: block there cites the DOI as the
# authoritative source. Wiring an automated TCIA fetch is tracked in ROADMAP.)
echo "$LOG_PREFIX IHC ROI set lives on TCIA (DOI 10.7937/TCIA.2020.T90F-WB82),"
echo "$LOG_PREFIX downloaded via Aspera Connect into data/pmc10571229/rois/ -- see header."
echo "$LOG_PREFIX done (source-repo cloned; ROI images are a manual TCIA/Aspera step)"

# ---- Step 3: record sha256 sums for whatever ROI images are present ---------
mkdir -p "$DEST/rois"
if find "$DEST/rois" -type f \( -name '*.png' -o -name '*.tif' -o -name '*.tiff' \) | head -1 | grep -q .; then
  echo "$LOG_PREFIX computing sha256 sums for ROI images"
  cd "$DEST"
  find . -type f \( -name '*.tif' -o -name '*.tiff' -o -name '*.png' -o -name '*.zip' -o -name '*.tar.gz' \) \
    -exec shasum -a 256 {} \; > "$DEST/sha256sums.txt"
  wc -l "$DEST/sha256sums.txt"
else
  echo "$LOG_PREFIX no ROI images present yet -- skipping sha256 (download from TCIA first)"
fi

echo "$LOG_PREFIX done"
