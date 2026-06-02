# What is out of scope (P4, `hnscc-time-multimodal`)

This file is the anti-scope-creep ledger for the P4 capability portrait.
The repo's value comes from being *small and complete*, every item below
is something a reviewer might reasonably ask for that the v0.0 -> v0.3
release plan deliberately does not attempt.

If a future PR proposes any of these, the contributor must answer one
question: **why is this still out of scope?** If the answer is good, edit
this file in the same PR. If not, the PR doesn't land.

---

## Paired patient-level multimodal training

The ideal multimodal model would train on the same patient's IHC + RNA-seq
+ outcome triplet. Such a dataset does not exist in public form at the
scale this repo could use.

**Why out of scope**: this is the central honesty point of the entire P4
plan. The PMC10571229 dataset has 8 HNSCC patients with multiplex IHC but
no genomics; TCGA-HNSC has ~530 HNSCC patients with RNA-seq but no IHC.
The repo openly uses Approach B (cross-cohort calibration via nearest-
neighbor matching) instead of pretending to have paired data. Any future
PR claiming paired multimodal training must produce the actual paired
dataset first.

---

## Foundation-model fine-tuning on the IHC images

UNI, Phikon, PathDino, and several other histology foundation models could
in principle be fine-tuned on the 268 PMC patches.

**Why out of scope**: 268 patches across 8 patients is two orders of
magnitude below the smallest defensible fine-tuning corpus for a histology
foundation model. The result would memorise patient identity, not learn
generalisable features. The honest architecture stays at off-the-shelf
Cellpose segmentation + per-cell marker classification.

---

## Vision-language model on IHC + clinical reports

VL models that take a histology patch + a clinical-report sentence and
produce structured output (TIME profile, immune phenotype call) are
plausible for HNSCC.

**Why out of scope**: the PMC dataset does not ship paired clinical
reports, and free-text TCGA pathology reports are not aligned to the
PMC patients. Building the report-image correspondence by hand for 8
patients would be a research project of its own, not a capability
portrait.

---

## Real-time interactive inference UI

A clinician-facing UI that takes a new RNA-seq panel and returns a TIME
profile would be the natural deployment of Arm 3.

**Why out of scope**: clinician-facing UIs require IRB framing, an actual
clinician collaborator, and a deployment substrate beyond the lab. The
substrate hooks (`audit.emit`, `tracking.log_*`) are the building blocks;
a real UI belongs to a downstream deployment project.

---

## Full TCGA-HNSC ingestion (~530 patients)

The first cut uses ~50 TCGA-HNSC patients (random subset, stratified by
HPV status + subsite). The full ~530 would give better deconvolution-tool
agreement statistics.

**Why out of scope (for v0.1)**: ~50 patients is enough to demonstrate the
deconvolution pipeline + calibration approach end-to-end on a workstation
in minutes. The ~530-patient run is a Phase 2 / v0.4 expansion that the
substrate is already sized for; it does not change the architecture.

---

## Additional deconvolution tools (CIBERSORTx absolute mode, MCPcounter)

Arm 2 v0.1 ships with one Python-resident deconvolution method
(xCell-equivalent gene-set scoring or EPIC port). The full plan envisions
2-3 methods for cross-tool agreement.

**Why out of scope (for v0.1)**: CIBERSORTx requires either a hosted
account or a non-trivial R toolchain install; MCPcounter is R-only. Adding
either pushes the v0.1 install footprint past "uv sync" friendliness.
v0.2 or v0.4 can add an R bridge if cross-tool agreement becomes a
substrate question.

---

## HPV+ vs HPV- subgroup survival modelling

TCGA-HNSC has HPV status and overall survival. A clinically interesting
extension is "does TIME profile predict survival differently for HPV+ vs
HPV-?" The data supports this analysis.

**Why out of scope (for v0.3)**: P3 (`tp53-aml-hrd-severity`) already
demonstrates the substrate's survival-modelling pattern. P4's contribution
is the multimodal integration, not a second survival analysis. The
survival follow-on belongs to a v0.4 / "applied" tag.

---

## ROI-level vs patch-level vs single-cell TIME profile aggregation

The plan reports per-ROI and per-patient TIME profiles. Per-patch or
per-single-cell aggregation surfaces would give finer-grained spatial
analysis.

**Why out of scope**: the common-schema TIME profile in
`docs/architecture.md` is a per-region (core / margin / stroma)
aggregation. Per-cell spatial statistics (cell-cell distances, infiltrate
gradients) are a substrate-extending PR. The current schema is the
minimum needed for cross-cohort calibration.

---

## Production hardening (HA, RBAC, multi-tenant, model registry)

The pipeline runs in a single Python process. There is no HA, no RBAC,
no per-tenant isolation, no streaming, no retry/backoff, no model
registry beyond MLflow artifacts.

**Why out of scope**: production hardening belongs to P1 (the
orchestration capability portrait) and to a separate deployment project,
not to the analytical-method portrait.

---

## Adding an item

Open a PR that:

1. Adds the item to the appropriate section above (or creates a new
   section if none fits).
2. Adds a one-sentence reason in italics for why it remains out of scope.
3. Links to the upstream feature request or issue if there is one.

That's it. The friction is intentional.
