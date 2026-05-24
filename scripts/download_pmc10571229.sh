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

# ---- Step 2: pull the Zenodo deposit (the actual ROI images, ~5 GB) --------
# The deposit is at https://zenodo.org/records/8367318 (or the canonical
# record referenced by the PMC paper). Use curl with --continue-at - so a
# resumed download survives Wi-Fi blips on chi-mac-p overnight.
#
# The Zenodo record-id is parameterised so v0.1 can swap to a mirror if
# needed.
# v0.0 note: record 8367318 was wrong (yielded a 7-KB stray .py). The correct
# PMC10571229 HNSCC ROI archive URL needs to be confirmed in v0.2 prep.
# Until then, refuse to run unless the caller overrides via env var.
ZENODO_RECORD="${PMC_ZENODO_RECORD:-}"
if [[ -z "$ZENODO_RECORD" ]]; then
  echo "$LOG_PREFIX SKIP Zenodo step: PMC10571229 record id unresolved (v0.0)."
  echo "$LOG_PREFIX Set PMC_ZENODO_RECORD=<id> to enable. Cloning source repo only."
  echo "$LOG_PREFIX done (source-repo only)"
  exit 0
fi
ZENODO_API="https://zenodo.org/api/records/${ZENODO_RECORD}"

mkdir -p "$DEST/rois"
echo "$LOG_PREFIX fetching Zenodo record metadata for record=$ZENODO_RECORD"
META=$(curl -fsSL "$ZENODO_API" || true)
if [[ -z "$META" ]]; then
  echo "$LOG_PREFIX WARNING: Zenodo API unreachable. v0.0 will ship without ROI data;"
  echo "$LOG_PREFIX v0.2 (Arm 1, IHC) will retry. Source repo above is still usable."
  echo "$LOG_PREFIX done (partial)"
  exit 0
fi

# Extract download URLs from the JSON metadata. Use python for robustness
# (jq is not assumed installed).
python3 - <<PY >/tmp/p4_pmc_files.txt
import json, sys
meta = json.loads("""$META""") if False else None
PY

# Simpler: use python3 with stdin
python3 - "$META" >/tmp/p4_pmc_files.txt <<'PY'
import json, sys
meta = json.loads(sys.argv[1])
for f in meta.get("files", []):
    print(f["links"]["self"], f["key"], f.get("size", 0))
PY

if [[ ! -s /tmp/p4_pmc_files.txt ]]; then
  echo "$LOG_PREFIX WARNING: no files in Zenodo record. Aborting Step 2."
  exit 0
fi

echo "$LOG_PREFIX Zenodo file list:"
cat /tmp/p4_pmc_files.txt | tee -a "$DEST/rois/_zenodo_manifest.txt"

while IFS=' ' read -r URL NAME SIZE; do
  OUT="$DEST/rois/$NAME"
  if [[ -f "$OUT" ]] && [[ "$(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT")" == "$SIZE" ]]; then
    echo "$LOG_PREFIX skip $NAME (size match)"
    continue
  fi
  echo "$LOG_PREFIX downloading $NAME ($SIZE bytes)"
  curl -fL --continue-at - -o "$OUT" "$URL" || {
    echo "$LOG_PREFIX WARNING: $NAME download failed; continuing"
  }
done < /tmp/p4_pmc_files.txt

# ---- Step 3: record sha256 sums for everything that landed -----------------
echo "$LOG_PREFIX computing sha256 sums"
cd "$DEST"
find . -type f \( -name '*.tif' -o -name '*.tiff' -o -name '*.png' -o -name '*.zip' -o -name '*.tar.gz' \) \
  -exec shasum -a 256 {} \; > "$DEST/sha256sums.txt"
wc -l "$DEST/sha256sums.txt"

echo "$LOG_PREFIX done"
