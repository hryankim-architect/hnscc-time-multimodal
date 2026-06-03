# Architecture, `hnscc-time-multimodal`

This repo stages a three-arm multimodal pipeline on top of the shared scaffold
substrate. The arms are released independently across v0.1 / v0.2 / v0.3, but
the substrate and the per-arm contract are fixed at v0.0 so each arm can land
without re-architecting.

---

## Three-arm control flow (target state, v0.3)

```
                                Operator
                                   |
                                   v
                       scripts/run_lab.sh   (substrate endpoints)
                                   |
                                   v
                              make run
                                   |
                                   v
                +---------------------------------+
                |  pipeline.run_pipeline(...)     |
                |  (src/hnscc_time/pipeline.py)   |
                +---------------------------------+
                                   |
                                   +-->  audit.emit('pipeline_start')
                                   |
                                   +-->  tracking.run().start_run()
                                   |
                                   v
        +---------------------------------------------------+
        |  Arm 2 (v0.1), Genomics deconvolution           |
        |    TCGA-HNSC bulk RNA-seq counts (~50 patients)   |
        |    -> normalize (TPM)                             |
        |    -> deconvolution (xCell-equiv / EPIC-equiv)    |
        |    -> per-patient TIME profile JSON               |
        |    audit.emit('genomics.time_profile.computed')   |
        +---------------------------------------------------+
                                   |
                                   v
        +---------------------------------------------------+
        |  Arm 1 (v0.2), IHC cell quantification          |
        |    PMC10571229 ROIs (72 ROIs, 8 patients)         |
        |    -> Cellpose nuclei segmentation                |
        |    -> per-cell marker classification (CD3/CD8/    |
        |       FoxP3/PanCK) from mIF channels              |
        |    -> per-ROI then per-patient TIME profile JSON  |
        |    audit.emit('ihc.time_profile.computed')        |
        +---------------------------------------------------+
                                   |
                                   v
        +---------------------------------------------------+
        |  Arm 3 (v0.3), Cross-cohort calibration         |
        |    For each PMC patient: nearest-neighbor lookup  |
        |      in TCGA-HNSC on subsite + age + HPV          |
        |    Fit a per-cell-type calibration mapping        |
        |      genomics_estimate -> ihc_equivalent          |
        |    Hold-out validation: leave-one-PMC-patient-out |
        |    audit.emit('calibration.trained')              |
        |    audit.emit('multimodal.prediction.served')     |
        +---------------------------------------------------+
                                   |
                                   +-->  tracking.log_metrics(...) + log_artifact
                                   |
                                   +-->  audit.emit('pipeline_end', status, metrics)
                                   |
                                   v
                         audit/local-demo.ndjson
                         (hash-chained, end-to-end verifiable)
```

At v0.0 the pipeline body only emits `pipeline_start` and `pipeline_end`;
the substrate is wired but the arms are stubs. Each subsequent release
replaces one stub with a real implementation while preserving the audit
contract above.

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

TCGA-HNSC bulk RNA-seq has no per-region resolution, so the genomics arm
populates the `"tumor_core"` slot only and leaves margin / stroma as null.
The calibration arm then expands the per-region prediction by transferring
the spatial pattern from the nearest IHC neighbour.

---

## Substrate integration points

Same four-channel substrate as the other portfolio repos:

| Channel | Module | Env var | Substrate endpoint |
|---|---|---|---|
| Audit (immutable record) | `hnscc_time.audit` | `AUDIT_HOST` | `http://${AUDIT_HOST}/events` |
| MLflow (experiment tracking) | `hnscc_time.tracking` | `MLFLOW_TRACKING_URI` | configurable |
| Canary (daily probe) | `hnscc_time.canary` | `HEALTHOMICS_LAB_CANARY_FIXTURE` | invoked by `lab_semantic_check.py` |
| Cohort + ROI registry | `data/manifest.yaml` (v0.0) -> SQLite (v0.1+) | (none) | local only |

All channels degrade to no-ops when the substrate is absent. The deterministic
local NDJSON ledger remains the source of truth for audit even when remote
POST fails.

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
and raw RNA-seq counts stay local, only manifest + checksums + derived
artifacts.

---

## Per-arm release contract

Each arm release MUST:

1. Add its real implementation under `src/hnscc_time/` (e.g. `genomics.py`,
   `ihc.py`, `calibrate.py`).
2. Wire it into `pipeline.run_pipeline(...)` with the audit emissions
   listed in the control-flow diagram above.
3. Add unit tests under `tests/` using synthetic fixtures (no real data
   needed in CI).
4. Add a one-page section to README "Climax" describing what just
   landed and how it ranks against the previous arm.
5. Bump `pyproject.toml` version and tag the release with `v0.1 / v0.2 /
   v0.3`.

The release contract guarantees that anyone cloning the repo at any tag
sees a complete, self-consistent snapshot of what existed at that point.

---

## What this architecture intentionally avoids

- **No DAG engine.** No Nextflow / Airflow / Prefect / Dagster. Three
  arms run in sequence inside one Python process. P1
  (`healthomics-lab-orchestrator`) handles DAG-engine orchestration.
- **No GPU dependency.** PyTorch MPS is opportunistic for Cellpose on
  Apple Silicon; CPU is sufficient for n=8 ROIs and n=50 patients.
- **No paired multimodal training.** No public dataset has paired
  patient-level IHC + RNA-seq for HNSCC at this scale. The repo openly
  uses cross-cohort calibration (Approach B) instead of pretending.
- **No deep generative model of histology.** The PMC dataset is too small
  (n=8 patients, 72 ROIs) for foundation-model fine-tuning to do anything
  except memorise. The honest architecture stays at segmentation +
  classification.

The contract is small and the implementation is small; the v0.1 -> v0.3
releases add real analysis, not architectural complexity.
