# Roadmap, `hnscc-time-multimodal`

This file is the public face of the sprint plan that turns v0.0 (scaffold
only) into v0.3 (three-arm multimodal pipeline). It is intentionally short
and dated so anyone clicking the repo can see exactly where the work
stands.

The full design rationale lives in the private design doc
`~/Downloads/AI/P4-IHC-Genomics-TIME-Plan.md` (Saturday 2026-05-23, v0.1).
This file mirrors only the public-facing engineering plan.

---

## v0.0, scaffold + data downloads (today, Sun 2026-05-24)

**Goal**: ship a substrate-clean repo on GitHub with the two real datasets
downloading in the background so Tuesday's work can start instantly.

- [x] Inherit scaffold-template; rename `bioscaffold` -> `hnscc_time`
- [x] P4-specific README, `docs/architecture.md`, `docs/what-is-out-of-scope.md`
- [x] `ROADMAP.md` (this file)
- [ ] Push to GitHub (`hryankim-architect/hnscc-time-multimodal`)
- [ ] Scaffold-level CI green (`ci` + `english-only`)
- [ ] Kick off background downloads:
  - PMC10571229 IHC dataset (~5 GB, via `scripts/download_pmc10571229.sh`)
  - TCGA-HNSC RNA-seq subset, n=50 (~2 GB, via `scripts/download_tcga_hnsc.sh`)
- [ ] v0.0 tag + release

---

## v0.1, Arm 2, Genomics deconvolution (Tue-Wed 2026-05-26/27)

**Goal**: ingest the TCGA-HNSC subset, produce per-patient TIME profile
JSON files in the common schema, wire into pipeline + audit.

- [ ] `src/hnscc_time/cohort.py`, TCGA-HNSC manifest loader + HPV /
      subsite stratification
- [ ] `src/hnscc_time/genomics.py`, TPM normalisation + gene-set
      deconvolution (xCell-equivalent or EPIC port; Python-resident
      so `uv sync` stays small)
- [ ] `src/hnscc_time/time_schema.py`, Pydantic model for the common
      TIME profile JSON (matches `docs/architecture.md` §Common-schema)
- [ ] Wire `pipeline.run_pipeline()` to call cohort -> genomics ->
      `time_profile.write()` for each patient
- [ ] Audit emissions: `cohort.tcga_hnsc.assembled`,
      `genomics.time_profile.computed.<patient_id>`
- [ ] Tests: synthetic-fixture deconvolution test + per-patient JSON
      schema validation
- [ ] README climax: per-cohort summary table (TIL score distribution,
      immune-phenotype call counts)
- [ ] v0.1 tag + release

---

## v0.2, Arm 1, IHC cell quantification (Thu 2026-05-28)

**Goal**: process the 72 PMC ROIs through Cellpose + per-marker
classification, produce per-patient TIME profile JSON files in the same
common schema.

- [ ] `src/hnscc_time/ihc.py`, Cellpose nuclei segmentation on DAPI
      channel + per-cell marker classification from CD3 / CD8 / FoxP3 /
      PanCK mIF channels
- [ ] Per-ROI aggregation -> per-region aggregation -> per-patient TIME
      profile JSON
- [ ] Audit emissions: `image.segmented.<patient_id>.<roi_id>`,
      `ihc.time_profile.computed.<patient_id>`
- [ ] Tests: synthetic 32x32 single-cell image fixture; assert
      segmentation + classification round-trip
- [ ] README climax: PMC8 vs paper Fig. 2 concordance (qualitative;
      formal concordance lands in v0.3 validation)
- [ ] v0.2 tag + release

---

## v0.3, Arm 3, Cross-cohort calibration (Fri 2026-05-29)

**Goal**: fit a calibration mapping from genomics-only TIME profiles to
IHC-equivalent TIME profiles, validate on held-out PMC patients, surface
a `predict_time_from_genomics()` callable.

- [ ] `src/hnscc_time/calibrate.py`, nearest-neighbor TCGA->PMC matcher
      on subsite + age + HPV; per-cell-type linear calibration mapping
- [ ] `predict_time_from_genomics(rna_seq_counts) -> TIMEProfile`
- [ ] Held-out validation: leave-one-PMC-patient-out; report
      per-cell-type MAE and immune-phenotype agreement
- [ ] Audit emissions: `calibration.trained.<version_id>`,
      `multimodal.prediction.served.<request_id>`
- [ ] Tests: end-to-end pipeline on synthetic fixture; tamper-detection
      test on the audit chain
- [ ] README climax: cross-cohort validation table (genomics-only vs
      calibrated vs IHC ground truth)
- [ ] v0.3 tag + release + final v0.3 release notes summarising all
      three arms

---

## After v0.3

Out of scope for this sprint (see `docs/what-is-out-of-scope.md`):

- HPV+ vs HPV- subgroup survival analysis (v0.4 candidate)
- Full ~530-patient TCGA-HNSC run (v0.4 candidate)
- Additional deconvolution tools via R bridge (v0.4 candidate)
- Production-style interactive inference UI (separate deployment project)
- Foundation-model fine-tuning on histology patches (out of scope
  permanently, n=8 patients is too small)

---

## Why this sequence

Arm 2 (Genomics) lands first because:
1. TCGA-HNSC data is more uniform than the PMC ROIs (no per-image
   pre-processing risk).
2. Python-resident deconvolution is install-light (no Cellpose / GPU
   dependency yet).
3. It establishes the common-schema TIME profile JSON, which Arms 1 and 3
   then conform to.

Arm 1 (IHC) lands second because Cellpose introduces the heaviest install
footprint (PyTorch + segmentation models), and it benefits from having the
schema already locked.

Arm 3 (Calibration) lands last because it depends on outputs from both
arms.

Each arm is a clean release; the repo at any tag is internally consistent.
