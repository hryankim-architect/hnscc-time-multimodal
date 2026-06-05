# Architecture, `hnscc-time-multimodal`

This repo builds a three-arm multimodal pipeline for HNSCC tumor immune
microenvironment (TIME) profiling. The three arms are released sequentially
across v0.1, v0.2, and v0.3; the substrate wiring and the per-arm audit
contract are locked at v0.0 so each arm can land without touching the others.

---

## Three-arm control flow (target state, v0.3)

The operator runs `scripts/run_lab.sh`, which calls `make run`, which calls
`pipeline.run_pipeline(...)` in `src/hnscc_time/pipeline.py`.

At startup the pipeline fires `audit.emit('pipeline_start')` and opens an
MLflow run via `tracking.run().start_run()`. Then the three arms execute in
sequence:

**Arm 2 — Genomics deconvolution (v0.1)**
TCGA-HNSC bulk RNA-seq counts (~50 patients) -> normalize to TPM ->
deconvolution (xCell-equivalent or EPIC-equivalent gene-set scoring) ->
per-patient TIME profile JSON -> `audit.emit('genomics.time_profile.computed')`

**Arm 1 — IHC cell quantification (v0.2)**
PMC10571229 ROIs (72 ROIs, 8 patients) -> Cellpose nuclei segmentation ->
per-cell marker classification (CD3 / CD8 / FoxP3 / PanCK) from mIF channels ->
per-ROI and per-patient TIME profile JSON ->
`audit.emit('ihc.time_profile.computed')`

**Arm 3 — Cross-cohort calibration (v0.3)**
For each PMC patient: nearest-neighbor lookup in TCGA-HNSC on subsite + age +
HPV -> fit a per-cell-type calibration mapping (genomics_estimate ->
ihc_equivalent) -> hold-out validation (leave-one-PMC-patient-out) ->
`audit.emit('calibration.trained')` -> `audit.emit('multimodal.prediction.served')`

After all arms finish, `tracking.log_metrics(...)` and `log_artifact` record
the run in MLflow, then `audit.emit('pipeline_end', status, metrics)` closes
the ledger. The final audit record lands in `audit/local-demo.ndjson`.

At v0.0 the pipeline body emits only `pipeline_start` and `pipeline_end`; the
arms are stubs. Each subsequent release replaces one stub with a real
implementation while preserving the audit contract above.

---

## Common-schema TIME profile

Both arms produce the same JSON schema so Arm 3 can compare them without
per-arm adapters:

```json
{
  "patient_id": "TCGA-CV-xxxx | Case3",
  "cohort": "tcga_hnsc | pmc10571229",
  "modality": "rna_seq | mIF | mIHC",
  "regions": {
    "tumor_core":     {"CD3_density": 0.0, "CD8_density": 0.0, "FoxP3_density": 0.0, "PanCK_density": 0.0},
    "tumor_margin":   {"CD3_density": 0.0, "CD8_density": 0.0, "FoxP3_density": 0.0, "PanCK_density": 0.0},
    "adjacent_stroma":{"CD3_density": 0.0, "CD8_density": 0.0, "FoxP3_density": 0.0, "PanCK_density": 0.0}
  },
  "TIL_score": 0.0,
  "immune_phenotype": "inflamed | excluded | desert | unknown",
  "provenance": {
    "method": "xCell-equivalent | Cellpose+marker_classifier | calibrated_predictor",
    "version": "v0.1.0",
    "ledger_id": "<sha256 of the audit entry that produced this profile>"
  }
}
```

TCGA-HNSC bulk RNA-seq has no per-region resolution, so Arm 2 populates
`"tumor_core"` only and leaves margin and stroma as null. Arm 3 expands the
per-region prediction by transferring the spatial pattern from the nearest
IHC neighbor.

---

## Substrate integration points

The pipeline uses four channels, all of which degrade to no-ops when the
relevant endpoint is absent:

**Audit** (`hnscc_time.audit`, env var `AUDIT_HOST`): posts events to
`http://${AUDIT_HOST}/events`. When `AUDIT_HOST` is unset, events are written
only to the local NDJSON file. The local file is always the authoritative
record; remote POST is secondary.

**MLflow** (`hnscc_time.tracking`, env var `MLFLOW_TRACKING_URI`): logs
metrics and artifacts. When the tracking URI points nowhere, all `log_*` calls
are no-ops. No run data is lost because the audit ledger captures the same
metrics independently.

**Canary** (`hnscc_time.canary`, env var `HEALTHOMICS_LAB_CANARY_FIXTURE`):
a daily probe invoked by `lab_semantic_check.py`. Checks that the pipeline
can complete a minimal synthetic run. Fails silently if the fixture env var
is absent.

**Cohort and ROI registry** (`data/manifest.yaml` at v0.0, SQLite at v0.1+):
local only, no env var. Tracks which cohorts are present and where the raw
files live.

---

## Audit ledger mechanics

Each entry in `audit/local-demo.ndjson` contains a `prev_hash` field set to
the SHA-256 of the preceding entry's canonical JSON serialization. This makes
the file a hash-chain: truncating or editing any entry invalidates all
subsequent hashes. Appending a new entry takes roughly 6.19 µs on a modern
laptop. The local file is the source of truth for audit regardless of whether
the remote POST succeeded.

---

## Data layout (v0.0 -> v0.3 progression)

```
data/
├── manifest.yaml               # human-editable: which cohorts, where they live
├── .gitignore                  # excludes raw downloads
├── pmc10571229/                # IHC dataset (v0.0 D download)
│   ├── README.md               # provenance + checksum
│   ├── rois/                   # 72 ROI images (1356x1012 px)
│   └── patches/                # 268 patches at 512x512
├── tcga_hnsc/                  # RNA-seq subset (v0.0 D download)
│   ├── README.md               # provenance + GDC query parameters
│   ├── clinical.tsv            # n=50 patient metadata
│   └── star_counts/            # per-patient STAR gene counts
└── time_profiles/              # produced by Arms 1 & 2
    ├── pmc10571229/
    │   └── Case<n>.json        # 8 files (v0.2)
    └── tcga_hnsc/
        └── TCGA-CV-xxxx.json   # 50 files (v0.1)
```

Only `manifest.yaml`, the per-cohort `README.md` files, and (eventually) the
50 + 8 = 58 small `time_profile.json` files are tracked in git. Raw images
and raw RNA-seq counts stay local; only manifest, checksums, and derived
artifacts are committed.

---

## Per-arm release contract

Each arm release must:

1. Add its real implementation under `src/hnscc_time/` (e.g. `genomics.py`,
   `ihc.py`, `calibrate.py`).
2. Wire it into `pipeline.run_pipeline(...)` with the audit emissions listed
   in the control-flow section above.
3. Add unit tests under `tests/` using synthetic fixtures (no real data
   needed in CI).
4. Add a one-page section to README "Climax" describing what just landed and
   how it ranks against the previous arm.
5. Bump `pyproject.toml` version and tag the release with `v0.1 / v0.2 / v0.3`.

Cloning the repo at any tag gives a complete, self-consistent snapshot of
what existed at that point.

---

## CI checks

The CI suite runs: ruff (lint + format), pytest (unit tests with synthetic
fixtures), and the canary probe (synthetic pipeline run).

---

## What this architecture intentionally avoids

**No DAG engine.** No Nextflow / Airflow / Prefect / Dagster. Three arms run
in sequence inside one Python process. A separate orchestration repo handles
DAG-engine work.

**No GPU dependency.** PyTorch MPS is opportunistic for Cellpose on Apple
Silicon; CPU is sufficient for n=8 ROIs and n=50 patients.

**No paired multimodal training.** No public dataset has paired patient-level
IHC + RNA-seq for HNSCC at this scale. The repo uses cross-cohort calibration
(Approach B) explicitly rather than treating the data gap as invisible.

**No deep generative model of histology.** The PMC dataset (n=8 patients,
72 ROIs) is too small for foundation-model fine-tuning. The architecture stays
at segmentation + classification.

The implementation stays small through v0.1 -> v0.3; each release adds real
analysis, not new infrastructure.
