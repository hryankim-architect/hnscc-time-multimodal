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
- **IHC — TCIA (resolved + checksummed):** the `HNSCC-mIF-mIHC-comparison`
  collection (DOI `10.7937/TCIA.2020.T90F-WB82`, 8 patients, ~1.01 GB, CC BY 4.0).
  TCIA's web bundle is Aspera, but every ROI is also directly HTTPS-fetchable from
  the TCIA pathology host, so all **3212 ROI PNGs are pinned by sha256** in
  `pmc10571229/rois_manifest.tsv` (rel_path + sha256 + direct URL; that ledger's
  own sha256 is pinned in the `ihc:` block). `scripts/download_pmc10571229.sh`
  fetches each by URL and verifies — no Aspera client needed.

**How the data is prepared:** `make data` (`pipeline fetch`) is now end-to-end —
it downloads the STAR `inputs` (sha256-verified), **derives
`tcga_hnsc/_subset_manifest.tsv`** from them, and **fetches + canonicalizes the
`clinical` block** (byte-stable, sha256-verified). That is exactly the layout
`cohort.load_cohort` reads, so `make data && make run` is reproducible from this
manifest alone. (`scripts/download_tcga_hnsc.sh` remains as an equivalent shell
path.) The IHC ROIs are fetched + sha256-verified by
`scripts/download_pmc10571229.sh` straight from `pmc10571229/rois_manifest.tsv`
(direct HTTPS, no Aspera).

Tiny fixtures for tests go under `tests/fixtures/`, not here. Keep
`manifest.yaml` lean: **every input must be small and necessary**. Adding one is
a PR-sized decision, not a casual edit.
