# hnscc-time-multimodal

![ci](https://github.com/hryankim-architect/hnscc-time-multimodal/actions/workflows/ci.yml/badge.svg) ![english-only](https://github.com/hryankim-architect/hnscc-time-multimodal/actions/workflows/english-only.yml/badge.svg)

> **Capability portrait, not a research result.** This repo demonstrates how
> a multimodal Tumor Immune Microenvironment (TIME) prediction pipeline can
> be staged on the same lab substrate that powers `multiqc-foundation-gate`
> and `tp53-aml-hrd-severity` — IHC + bulk RNA-seq deconvolution + cross-
> cohort calibration, every step bracketed by a hash-chained audit ledger.
> The clinical/biological claim space is intentionally narrow at v0; the
> substrate and the integration pattern are the deliverable.

## What this answers

For a head-and-neck squamous cell carcinoma (HNSCC) patient where you have
**either** multiplex IHC **or** bulk RNA-seq (almost never both, in
real practice), can you produce a comparable Tumor Immune Microenvironment
profile — quantified the same way, on the same schema, with the same audit
trail — and use it to drive an immunotherapy patient-selection signal?

That question is the multimodal CDx pattern in miniature. This repo builds
it on entirely public data:

- **Multiplex IHC ground truth**: PMC10571229 (Ghahremani et al. 2023,
  MICCAI) — 8 HNSCC patients from Moffitt Cancer Center, 72 ROIs across
  tumor core / margin / adjacent stroma, markers DAPI + CD3 + CD8 + FoxP3
  + PanCK.
- **Bulk RNA-seq cohort**: TCGA-HNSC — ~530 HNSCC patients with paired
  RNA-seq + clinical (HPV status, survival, subsite).
- **Cross-cohort integration**: Approach B from the design doc — calibrate
  a genomics-only TIME predictor against IHC ground truth via nearest-
  neighbor matching on clinical + subsite features.

## Release status (compressed to one session 2026-05-24)

The original sprint plan estimated v0.1 / v0.2 / v0.3 across Tue-Fri, but
the entire 3-arm pipeline landed in a single Sunday-evening session after
v0.0 shipped. The trade-offs that made this possible are documented openly
in the per-arm sections below.

| Layer | v0.0 (Sun AM) | v0.1 (Sun PM, Arm 2) | v0.2 (Sun PM, Arm 1) | v0.3 (Sun PM, Arm 3) |
|---|---|---|---|---|
| Substrate (audit / tracking / canary) | ✓ | ✓ | ✓ | ✓ |
| Repo skeleton + CI + english-only | ✓ | ✓ | ✓ | ✓ |
| **Arm 2 — Genomics deconvolution on RNA-seq** | — | **✓ TCGA-HNSC n=50, ssGSEA-style scoring on curated immune signatures** | ✓ | ✓ |
| **Arm 1 — IHC cell segmentation** | — | — | **✓ Cellpose nuclei on 5 real DeepLIIF Sample_Large_Tissues ROIs** | ✓ |
| **Arm 3 — Cross-cohort calibration + `predict_time_from_genomics()`** | — | — | — | **✓ NN + per-cell-type linear cal + LOO validation** |

See `ROADMAP.md` for what the released contract looks like at each tag,
and `docs/what-is-out-of-scope.md` for what each arm intentionally does
not attempt (paired multimodal training, foundation-model fine-tuning on
n=8, full ~530-patient TCGA-HNSC ingestion, etc.).

## Real-data climax (`make run`, n=50 TCGA-HNSC + 5 DeepLIIF ROIs)

End-to-end Sunday-evening smoke produced this — the substrate value is
the *comparison*, not any single arm's number:

| Arm | What ran | Headline metric |
|---|---|---|
| Arm 2 (Genomics) | 50 TCGA-HNSC patients, TPM -> rank-transform -> mean-rank of immune signatures, z-score normalised | TIL score mean **0.637 ± 0.164** across cohort; 47 inflamed / 1 excluded / 2 desert / 0 unknown |
| Arm 1 (IHC) | 5 DeepLIIF Sample_Large_Tissues ROIs (RGB tissue), Cellpose 4.x CPU segmentation | **6,725 nuclei segmented**; 5 inflamed / 0 excluded / 0 desert (RGB heuristic R/G/B->marker placeholder) |
| Arm 3 (Calibration) | K-NN on n=5 IHC reference, per-cell-type linear cal, leave-one-IHC-out validation | **mean LOO MAE 0.210 vs intercept-only baseline 0.466 — calibration adds 55% MAE reduction over "predict cohort mean"** |

The Arm 3 headline is the most defensible single number this repo
produces: it shows that nearest-neighbor cross-cohort calibration is
genuinely extracting signal over a "just predict the average" baseline,
even on the deliberately tiny n=5 IHC reference. The full 72-ROI
PMC10571229 archive (Arm 1 v0.4 candidate) would tighten this further,
but the *pattern* is already visible.

### What the 55% number means and does not mean

It *means*: the K-NN-anchored linear calibration is doing real work — it
beats "predict the IHC cohort mean for every patient" by half its error.

It *does not mean*: the calibration is research-grade. n=5 IHC is one
order of magnitude below what's needed for any peer-reviewable claim.
The number lives in this README because it is **a working capability
demo on real public data**, not because it is a finding.

## Why this scoping

The full P4 plan estimates ~3 weeks part-time. Compressing the whole thing
into one weekend would force synthetic data and short-cut validation, which
would make the repo a demo rather than a capability portrait. The chosen
sequence — substrate today, one arm at a time over the week, integration at
the end — keeps every commit defensible.

## What the substrate already gives you on day 0

Even with no analysis code yet, this repo:

- Boots a single Python process via `scripts/run_lab.sh` that brackets
  every operation with `pipeline_start` / `pipeline_end` audit emissions.
- Writes a hash-chained NDJSON ledger (`audit/local-demo.ndjson`) that
  `audit.verify()` can walk end-to-end to detect tamper.
- Degrades cleanly to no-op when MLflow / `AUDIT_HOST` are unset, so the
  scaffold passes CI on a vanilla GitHub Actions runner.
- Enforces an English-only public surface via the `english-only` workflow
  (CJK character scanner in `scripts/check_english_only.py`).

Running `make test` against the v0.0 commit exercises this surface end-to-
end before any P4-specific work lands.

## Repo layout

```
hnscc-time-multimodal/
├── src/hnscc_time/         # substrate (audit / tracking / canary / pipeline)
│   ├── audit.py            # hash-chained NDJSON ledger
│   ├── tracking.py         # MLflow + no-op fallback
│   ├── canary.py           # daily probe (consumed by lab_semantic_check.py)
│   └── pipeline.py         # outer bracket; v0.1+ wires Arm 2 in here
├── data/                   # manifest.yaml + .gitignored bulk data
├── docs/
│   ├── architecture.md     # 3-arm control flow + substrate channels
│   └── what-is-out-of-scope.md
├── scripts/
│   ├── run_lab.sh          # substrate-aware entrypoint
│   ├── check_english_only.py
│   ├── download_pmc10571229.sh    # (v0.0 step D) IHC dataset fetcher
│   └── download_tcga_hnsc.sh      # (v0.0 step D) RNA-seq subset fetcher
├── tests/                  # scaffold smoke tests; per-arm tests land with each arm
├── Makefile                # data / run / test / clean targets
├── pyproject.toml          # uv-managed deps; arm-specific extras land per release
├── ROADMAP.md              # day-by-day plan to v0.3
└── README.md               # this file
```

## Quickstart

```bash
# 1. Clone + install
git clone https://github.com/hryankim-architect/hnscc-time-multimodal.git
cd hnscc-time-multimodal
uv sync   # installs scaffold deps; arm deps come in v0.1+

# 2. Run the scaffold smoke pipeline (no real data needed)
make test
make run          # writes audit/local-demo.ndjson

# 3. Verify the audit chain
python -c "from hnscc_time import audit; ok, n, bad = audit.verify(); print(ok, n, bad)"
# Expect: True <small number> None
```

## Honest scope (v0.3)

What this release **does**:
- Runs end-to-end on 50 real TCGA-HNSC patients (Arm 2 Genomics) and 5
  real DeepLIIF Sample_Large_Tissues ROIs (Arm 1 IHC), then cross-cohort
  calibrates Arm 3 to produce 50 calibrated TIMEProfile predictions.
- Emits a hash-chained audit ledger that spans `pipeline_start` ->
  per-arm cohort assembly + profile computation + calibration -> `pipeline_end`.
- Surfaces a `predict_time_from_genomics(genomics_profile) -> TIMEProfile`
  callable that is the deployment-ready unit of work.

What this release **does not** claim:
- **Statistical generalisation** — n=5 IHC reference is one order of
  magnitude below what the PMC10571229 full archive would give; the
  calibration is a *demonstration of the integration pattern*, not a
  research finding. The README climax table reports the held-out LOO MAE
  honestly alongside the intercept-only baseline so a reader can see
  exactly how much signal calibration adds.
- **True per-marker mIF channels** — the 5 DeepLIIF Sample_Large_Tissues
  are RGB composites, not the 4-channel mIF stack the PMC dataset would
  provide. Arm 1 v0.2 uses a documented *heuristic R/G/B -> PanCK/CD8/CD3
  mapping* as a placeholder; Arm 3 calibration replaces it with the
  genomics-anchored signal where the IHC channels would otherwise
  dominate.
- **Paired patient-level multimodal training** — the 8 PMC patients are
  not in TCGA-HNSC; we openly use cross-cohort calibration (Approach B)
  instead of pretending to have paired data.

See `docs/what-is-out-of-scope.md` for the longer list of things this repo
will *never* attempt (training a vision-language foundation model, paired
multimodal training when no paired data exists, etc.) and the rationale
for each.

## Substrate integration

Same four-channel substrate as the other capability-portrait repos:

| Channel | Module | Env var |
|---|---|---|
| Audit (immutable record) | `hnscc_time.audit` | `AUDIT_HOST` |
| MLflow tracking | `hnscc_time.tracking` | `MLFLOW_TRACKING_URI` |
| Daily canary | `hnscc_time.canary` | invoked by `lab_semantic_check.py` |
| Cohort + ROI registry | `data/manifest.yaml` (v0.0) -> SQLite (v0.1+) | — |

All four degrade to no-ops in the absence of the substrate, so a public
clone-and-run works without any private lab services running.

## Related portfolio repos

- [scaffold-template](https://github.com/hryankim-architect/bioinformatics-repo-scaffold-template)
  — the parent substrate this repo inherits from.
- [multiqc-foundation-gate](https://github.com/hryankim-architect/multiqc-foundation-gate)
  — P2, sample-QC classifier comparison on MultiQC features.
- [tp53-aml-hrd-severity](https://github.com/hryankim-architect/tp53-aml-hrd-severity)
  — P3, clinical-genomics survival modeling on TCGA-LAML.
- [healthomics-lab-orchestrator](https://github.com/hryankim-architect/healthomics-lab-orchestrator)
  — P1, Nextflow + audit-bracketed RNA-seq pipeline orchestrator.

## License

MIT (see `LICENSE`).
