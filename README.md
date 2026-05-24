# hnscc-time-multimodal

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

## Release status (v0.0 — scaffold only)

This release ships **the substrate, not the analysis**. It is intentionally
labelled v0.0 to make that obvious: the v0.1 / v0.2 / v0.3 sprints over the
following week build the three arms on top.

| Layer | v0.0 (today) | Lands in |
|---|---|---|
| Substrate (audit / tracking / canary) | ✓ inherited from scaffold-template | — |
| Repo skeleton (pyproject + Makefile + CI + english-only) | ✓ | — |
| Background data downloads (PMC10571229 + TCGA-HNSC subset) | ✓ kicked off | — |
| Arm 2 — Genomics deconvolution on RNA-seq | — | **v0.1 (Tue-Wed)** |
| Arm 1 — IHC cell segmentation + TIME profile (Cellpose on 72 ROIs) | — | **v0.2 (Thu)** |
| Arm 3 — Cross-cohort calibration (Approach B + held-out validation) | — | **v0.3 (Fri)** |

See `ROADMAP.md` for the day-by-day sprint plan and the cross-link to the
full design doc.

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

## Honest scope (v0.0)

This release does not yet:

- Touch a single multiplex-IHC image or RNA-seq count file.
- Run Cellpose, xCell, EPIC, quanTIseq, or any deconvolution method.
- Produce a TIME profile JSON.
- Predict anything clinical.

What it does is establish that the next three sprints (v0.1 / v0.2 / v0.3)
have somewhere to land that is already CI-green, audit-verifiable, and
plugged into the same substrate as P1 / P2 / P3.

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
