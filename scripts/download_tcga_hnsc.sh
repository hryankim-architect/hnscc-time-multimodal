#!/usr/bin/env bash
# Download a stratified n=50 subset of TCGA-HNSC RNA-seq + clinical via the
# NIH GDC REST API (open tier; no controlled-access auth needed).
#
# Usage:
#   nohup bash scripts/download_tcga_hnsc.sh >/tmp/p4_data_tcga.log 2>&1 &
#
# Idempotent: re-running on top of a partial run skips already-downloaded
# files. Re-running with a larger SUBSET_N expands the cohort.

set -uo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DEST="$REPO_ROOT/data/tcga_hnsc"
SUBSET_N="${SUBSET_N:-50}"
LOG_PREFIX="[p4_tcga_download $(date -u +%Y-%m-%dT%H:%M:%SZ)]"

mkdir -p "$DEST/star_counts"
echo "$LOG_PREFIX starting; dest=$DEST n=$SUBSET_N"

GDC="https://api.gdc.cancer.gov"

# ---- Step 1: query GDC for TCGA-HNSC RNA-seq STAR counts -------------------
# Filter: project=TCGA-HNSC, experimental_strategy=RNA-Seq,
#         data_type=Gene Expression Quantification,
#         workflow_type=STAR - Counts, access=open.
#
# Output: a list of {file_id, case_id} pairs that we then random-stratify by
# HPV status.

QUERY=$(cat <<'JSON'
{
  "filters": {
    "op": "and",
    "content": [
      {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-HNSC"]}},
      {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
      {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
      {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
      {"op": "in", "content": {"field": "files.access", "value": ["open"]}}
    ]
  },
  "fields": "file_id,file_name,cases.submitter_id,cases.case_id,cases.demographic.gender,cases.diagnoses.primary_diagnosis,cases.diagnoses.tissue_or_organ_of_origin",
  "format": "JSON",
  "size": "1000"
}
JSON
)

echo "$LOG_PREFIX querying GDC for available TCGA-HNSC STAR counts..."
curl -fsSL -H "Content-Type: application/json" -d "$QUERY" "$GDC/files" \
  > "$DEST/_gdc_file_list.json" || {
    echo "$LOG_PREFIX ERROR: GDC API unreachable. Aborting; retry tomorrow."
    exit 2
}
N_AVAILABLE=$(python3 -c "import json; d=json.load(open('$DEST/_gdc_file_list.json')); print(len(d['data']['hits']))")
echo "$LOG_PREFIX GDC returned $N_AVAILABLE candidate files"

# ---- Step 2: random-stratified subset of SUBSET_N --------------------------
python3 - <<PY > "$DEST/_subset_manifest.tsv"
import json, random
random.seed(42)
d = json.load(open("$DEST/_gdc_file_list.json"))
hits = d["data"]["hits"]
# Crude stratification: 25 random + 25 random (we don't have HPV status from
# the file list endpoint without a separate cases query; the genomics arm
# will join HPV status from the clinical TSV pulled separately).
random.shuffle(hits)
subset = hits[: $SUBSET_N]
print("file_id\tcase_submitter_id\tfile_name")
for h in subset:
    sub = h.get("cases", [{}])[0].get("submitter_id", "?")
    print(f"{h['file_id']}\t{sub}\t{h.get('file_name','?')}")
PY

echo "$LOG_PREFIX subset manifest:"
head -5 "$DEST/_subset_manifest.tsv"
wc -l "$DEST/_subset_manifest.tsv"

# ---- Step 3: download each STAR counts file --------------------------------
echo "$LOG_PREFIX downloading STAR counts files..."
tail -n +2 "$DEST/_subset_manifest.tsv" | while IFS=$'\t' read -r FILE_ID CASE_ID FILE_NAME; do
  OUT="$DEST/star_counts/${CASE_ID}__${FILE_NAME}"
  if [[ -f "$OUT" ]] && [[ "$(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT")" -gt 100000 ]]; then
    echo "$LOG_PREFIX skip $CASE_ID (already present, >100k)"
    continue
  fi
  echo "$LOG_PREFIX fetching $CASE_ID file_id=$FILE_ID"
  curl -fL --continue-at - -o "$OUT" "$GDC/data/$FILE_ID" || {
    echo "$LOG_PREFIX WARNING: $CASE_ID download failed; continuing"
  }
  # GDC sends gzipped TSV; leave compressed (genomics arm reads gzipped natively)
done

# ---- Step 4: pull clinical metadata for the selected cases -----------------
CASE_IDS=$(tail -n +2 "$DEST/_subset_manifest.tsv" | cut -f2 | paste -sd, -)
echo "$LOG_PREFIX requesting clinical metadata for selected cases"

CLINICAL_QUERY=$(python3 - <<PY
import json
ids = """$CASE_IDS""".strip().split(",")
q = {
  "filters": {"op": "in", "content": {"field": "submitter_id", "value": ids}},
  "fields": "submitter_id,demographic.gender,demographic.race,demographic.year_of_birth,demographic.vital_status,diagnoses.primary_diagnosis,diagnoses.tissue_or_organ_of_origin,diagnoses.days_to_last_follow_up,diagnoses.days_to_death,exposures.tobacco_smoking_status",
  "format": "TSV",
  "size": str(len(ids) + 10)
}
print(json.dumps(q))
PY
)
curl -fsSL -H "Content-Type: application/json" -d "$CLINICAL_QUERY" "$GDC/cases" \
  > "$DEST/clinical.tsv" || {
    echo "$LOG_PREFIX WARNING: clinical TSV pull failed; v0.1 retries"
  }

# ---- Step 5: record sha256 sums --------------------------------------------
echo "$LOG_PREFIX computing sha256 sums"
cd "$DEST"
find . -type f \( -name '*.gz' -o -name '*.tsv' -o -name '*.json' \) -not -name '_*' \
  -exec shasum -a 256 {} \; > "$DEST/sha256sums.txt"
wc -l "$DEST/sha256sums.txt"

echo "$LOG_PREFIX done"
