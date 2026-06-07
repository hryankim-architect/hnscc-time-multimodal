# `data/`

This directory is **not** for committed data. `manifest.yaml` is the committed
**checksum ledger** for the demo's public inputs; raw data is git-ignored.

Two arms (schema documented inline in `manifest.yaml`):

- **genomics — TCGA-HNSC (NIH GDC open tier):** a deterministic n=50 subset of
  RNA-Seq STAR-Counts (all open STAR-Counts file_ids sorted, then `seed=42`
  shuffle, then first 50) plus a clinical `/cases` query. Every STAR file is
  pinned by sha256 (`inputs:`); the clinical TSV is canonicalized (header +
  sorted rows) and pinned (`clinical:`). STAR counts are **plain TSV** — GDC
  does not gzip them. This is a plain random subset, **not** HPV-stratified.
- **IHC — TCIA (resolved):** the `HNSCC-mIF-mIHC-comparison` collection
  (DOI `10.7937/TCIA.2020.T90F-WB82`, 8 patients, ~1.01 GB, CC BY 4.0). TCIA
  distributes via Aspera, so per-file sha256 are **not** pinned here; the `ihc:`
  block records the DOI as the authoritative source.

**How the data is actually prepared:** run `scripts/download_tcga_hnsc.sh`, which
writes `tcga_hnsc/_subset_manifest.tsv` + `star_counts/<case>__<file>.tsv` +
`clinical.tsv` — the layout `cohort.load_cohort` reads. `make data`
(`pipeline fetch`) downloads the STAR `inputs` by URL to the same paths and
verifies sha256, but does not yet build the subset index or fetch clinical;
wiring `pipeline fetch` end-to-end is a tracked follow-up (see `ROADMAP.md`).
The IHC ROI images are a manual TCIA/Aspera step (see `scripts/download_pmc10571229.sh`).

Tiny fixtures for tests go under `tests/fixtures/`, not here. Keep
`manifest.yaml` lean: **every input must be small and necessary**. Adding one is
a PR-sized decision, not a casual edit.
