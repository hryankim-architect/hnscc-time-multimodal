# Roadmap, `hnscc-time-multimodal`

This file is the public face of the sprint plan that turns v0.0 (scaffold
only) into v0.3 (three-arm multimodal pipeline). It is intentionally short
and dated so anyone clicking the repo can see exactly where the work
stands.

The full design rationale lives in `docs/architecture.md`.
This file mirrors only the public-facing engineering plan.

> **Data-layer follow-up (2026-06-07).** `data/manifest.yaml` is now a real
> checksum ledger: the genomics arm pins a deterministic n=50 TCGA-HNSC
> STAR-Counts subset + canonical clinical TSV by sha256, and the IHC arm is fully
> resolved + checksummed — all 3212 TCIA ROIs (DOI `10.7937/TCIA.2020.T90F-WB82`)
> are pinned by sha256 in `data/pmc10571229/rois_manifest.tsv` and fetched by
> direct HTTPS (no Aspera). **`pipeline fetch` is now wired end-to-end** — it
> builds `tcga_hnsc/_subset_manifest.tsv` from the STAR inputs and fetches the
> canonical `clinical:` block, so `make data && make run` is reproducible from
> the manifest alone (the `scripts/download_*.sh` scripts remain as an equivalent
> path).

> **Status: shipped.** Milestones v0.0–v0.3 are complete and tagged (`v0.0`,
> `v0.3`); the three arms run end to end and the as-built result tables are in
> the README. The checklist below is kept as the original sprint plan. One
> deviation from it: the IHC arm (v0.2) shipped on 5 DeepLIIF
> Sample_Large_Tissues ROIs rather than the originally scoped PMC10571229 set;
> the README carries the as-built numbers (6,725 nuclei segmented; cross-cohort
> LOO MAE 0.210 vs 0.466 intercept-only).

---

## v0.0, scaffold + data downloads (today, Sun 2026-05-24)

**Goal**: ship a substrate-clean repo on GitHub with the two real datasets
downloading in the background so Tuesday's work can start instantly.

- [x] Inherit scaffold-template; rename `bioscaffold` -> `hnscc_time`
- [x] P4-specific README, `docs/architecture.md`, `docs/what-is-out-of-scope.md`
- [x] `ROADMAP.md` (this file)
- [x] Push to GitHub (`hryankim-architect/hnscc-time-multimodal`)
- [x] Scaffold-level CI green (`ci`)
- [x] Kick off background downloads:
  - PMC10571229 IHC dataset (~5 GB, via `scripts/download_pmc10571229.sh`)
  - TCGA-HNSC RNA-seq subset, n=50 (~2 GB, via `scripts/download_tcga_hnsc.sh`)
- [x] v0.0 tag + release

---

## v0.1, Arm 2, Genomics deconvolution (Tue-Wed 2026-05-26/27)

**Goal**: ingest the TCGA-HNSC subset, produce per-patient TIME profile
JSON files in the common schema, wire into pipeline + audit.

- [x] `src/hnscc_time/cohort.py`, TCGA-HNSC manifest loader + HPV /
      subsite stratification
- [x] `src/hnscc_time/genomics.py`, TPM normalisation + gene-set
      deconvolution (xCell-equivalent or EPIC port; Python-resident
      so `uv sync` stays small)
- [x] `src/hnscc_time/time_schema.py`, Pydantic model for the common
      TIME profile JSON (matches `docs/architecture.md` §Common-schema)
- [x] Wire `pipeline.run_pipeline()` to call cohort -> genomics ->
      `time_profile.write()` for each patient
- [x] Audit emissions: `cohort.tcga_hnsc.assembled`,
      `genomics.time_profile.computed.<patient_id>`
- [x] Tests: synthetic-fixture deconvolution test + per-patient JSON
      schema validation
- [x] README climax: per-cohort summary table (TIL score distribution,
      immune-phenotype call counts)
- [x] v0.1 tag + release

---

## v0.2, Arm 1, IHC cell quantification (Thu 2026-05-28)

**Goal**: process the 72 PMC ROIs through Cellpose + per-marker
classification, produce per-patient TIME profile JSON files in the same
common schema.

- [x] `src/hnscc_time/ihc.py`, Cellpose nuclei segmentation on DAPI
      channel + per-cell marker classification from CD3 / CD8 / FoxP3 /
      PanCK mIF channels
- [x] Per-ROI aggregation -> per-region aggregation -> per-patient TIME
      profile JSON
- [x] Audit emissions: `image.segmented.<patient_id>.<roi_id>`,
      `ihc.time_profile.computed.<patient_id>`
- [x] Tests: synthetic 32x32 single-cell image fixture; assert
      segmentation + classification round-trip
- [x] README climax: PMC8 vs paper Fig. 2 concordance (qualitative;
      formal concordance lands in v0.3 validation)
- [x] v0.2 tag + release

---

## v0.3, Arm 3, Cross-cohort calibration (Fri 2026-05-29)

**Goal**: fit a calibration mapping from genomics-only TIME profiles to
IHC-equivalent TIME profiles, validate on held-out PMC patients, surface
a `predict_time_from_genomics()` callable.

- [x] `src/hnscc_time/calibrate.py`, nearest-neighbor TCGA->PMC matcher
      on subsite + age + HPV; per-cell-type linear calibration mapping
- [x] `predict_time_from_genomics(rna_seq_counts) -> TIMEProfile`
- [x] Held-out validation: leave-one-PMC-patient-out; report
      per-cell-type MAE and immune-phenotype agreement
- [x] Audit emissions: `calibration.trained.<version_id>`,
      `multimodal.prediction.served.<request_id>`
- [x] Tests: end-to-end pipeline on synthetic fixture; tamper-detection
      test on the audit chain
- [x] README climax: cross-cohort validation table (genomics-only vs
      calibrated vs IHC ground truth)
- [x] v0.3 tag + release + final v0.3 release notes summarising all
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
