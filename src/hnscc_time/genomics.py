"""Genomics arm — bulk RNA-seq -> per-patient TIME profile.

Pipeline:
    STAR-counts TSV (per patient, has both raw counts and TPM_unstranded)
    -> load TPM matrix (gene_name index, single-sample column)
    -> log1p(TPM) for variance stabilisation
    -> rank-transform genes within sample (ssGSEA-style)
    -> for each immune-cell signature, mean rank of signature genes
    -> normalise to z-scores across the cohort
    -> per-patient TIMEProfile (tumor_core only; margin + stroma None)

Design choices:
    - No external deconvolution dependency (no decoupler-py, xCell, EPIC).
      We hand-roll a minimal ssGSEA-style scorer + ship the immune-cell
      signature gene lists inline. This keeps `uv sync` slim, makes the
      whole arm Python-only, and is reproducible across runs.
    - Signature genes are well-established human immune-cell markers.
      Sources documented per-list. v0.1 deliberately uses canonical
      single-gene markers (the strongest signal) rather than long
      curated panels; this matches the 4-marker IHC schema (CD3 / CD8 /
      FoxP3 / PanCK) and keeps cross-arm comparison tight.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from hnscc_time.time_schema import (
    CELL_TYPES,
    Provenance,
    RegionProfile,
    TIMEProfile,
    derive_immune_phenotype,
    derive_til_score,
)

# ---------------------------------------------------------------------------
# Immune-cell signature gene lists (HGNC symbols, GENCODE v36 compatible)
# ---------------------------------------------------------------------------
# Each list expands the IHC-arm single-marker into a short transcript
# signature so the scoring is robust to sample-level noise. The first
# entry of each list is the canonical IHC marker.
#
# Sources:
#   - CD3 (pan T-cell):        CD3D / CD3E / CD3G — universally cited
#   - CD8 (cytotoxic T-cell):  CD8A / CD8B + GZMB / PRF1 (cytotoxic effector
#                              programme; Charoentong et al. 2017 immune
#                              landscape paper)
#   - FoxP3 (Treg):            FOXP3 + IL2RA / CTLA4 (canonical Treg panel)
#   - PanCK (epithelial/tumour): KRT5 / KRT6A / KRT14 / KRT17 (squamous
#                              keratins relevant for HNSCC; widely cited)
#
# Keep the lists short (3-5 genes) so a single missing/lowly-expressed
# transcript does not destroy the signal.
IMMUNE_SIGNATURES: dict[str, list[str]] = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B", "GZMB", "PRF1"],
    "FoxP3": ["FOXP3", "IL2RA", "CTLA4"],
    "PanCK": ["KRT5", "KRT6A", "KRT14", "KRT17"],
}

# Pipeline version recorded in TIMEProfile provenance.
GENOMICS_METHOD = "ssGSEA_immune_signatures"
GENOMICS_VERSION = "v0.1.0"


def _open_star_counts(path: Path):
    """Open a STAR counts file (gzipped or plain) for read."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_tpm(path: Path) -> pd.Series:
    """Load gene_name -> tpm_unstranded from one STAR counts TSV.

    The GDC STAR counts TSV format:
        # gene-model: GENCODE v36
        gene_id  gene_name  gene_type  unstranded  stranded_first  stranded_second  tpm_unstranded  fpkm_unstranded  fpkm_uq_unstranded
        N_unmapped     ...
        N_multimapping ...
        N_noFeature    ...
        N_ambiguous    ...
        ENSG00...      TSPAN6  protein_coding  1597  770  827  21.5889  ...
        ...

    The first 4 N_* rows have NaN in numeric columns; we drop those.
    Multiple Ensembl IDs can map to the same gene_name (e.g. PAR genes);
    we sum TPMs so each gene_name is a single value.
    """
    with _open_star_counts(path) as fh:
        # Skip the "# gene-model" comment line.
        df = pd.read_csv(fh, sep="\t", comment="#")
    # Drop the N_* QC rows (their gene_name is NaN).
    df = df.dropna(subset=["gene_name"])
    # Sum across Ensembl-ID-level duplicates per gene_name.
    s = df.groupby("gene_name")["tpm_unstranded"].sum()
    return s.astype(float)


def rank_transform(tpm: pd.Series) -> pd.Series:
    """Log1p + rank-transform genes within a single sample.

    Returns ranks in [1, n_genes]. Higher rank = higher expression.
    """
    log_tpm = np.log1p(tpm)
    return log_tpm.rank(method="average", ascending=True)


def signature_score(ranks: pd.Series, sig_genes: list[str]) -> float:
    """Mean rank of the signature genes within one sample.

    Missing genes are dropped from the mean (not zero-filled — a missing
    gene is "no information", not "no expression").
    """
    present = [g for g in sig_genes if g in ranks.index]
    if not present:
        return 0.0
    return float(ranks.loc[present].mean())


def score_sample(tpm: pd.Series) -> dict[str, float]:
    """Per-cell-type signature score for a single sample."""
    ranks = rank_transform(tpm)
    return {cell: signature_score(ranks, IMMUNE_SIGNATURES[cell]) for cell in CELL_TYPES}


def cohort_score_matrix(cohort_df: pd.DataFrame) -> pd.DataFrame:
    """Score every patient in the cohort. Returns DataFrame [submitter_id x cell_type]."""
    rows = {}
    for sid, row in cohort_df.iterrows():
        path = row.get("star_counts_path")
        if path is None or pd.isna(path):
            continue
        path = Path(path)
        if not path.exists():
            continue
        tpm = load_tpm(path)
        rows[sid] = score_sample(tpm)
    return pd.DataFrame.from_dict(rows, orient="index", columns=list(CELL_TYPES))


def cohort_zscore_normalize(scores: pd.DataFrame) -> pd.DataFrame:
    """Z-score across the cohort per cell type so 0 = cohort mean.

    Then shift+clip into [0, 5] so the values are non-negative
    (TIMEProfile.RegionProfile requires ge=0). The shift is +2.5 so a
    cohort-mean sample lands at 2.5 — well inside the valid range.
    """
    z = (scores - scores.mean()) / scores.std(ddof=0).replace(0, 1.0)
    shifted = (z + 2.5).clip(lower=0.0, upper=5.0)
    return shifted


def to_time_profile(submitter_id: str, normalised: pd.Series, ledger_id: str | None = None) -> TIMEProfile:
    """Convert one row of the normalised score matrix into a TIMEProfile."""
    region = RegionProfile(
        CD3_density=float(normalised["CD3"]),
        CD8_density=float(normalised["CD8"]),
        FoxP3_density=float(normalised["FoxP3"]),
        PanCK_density=float(normalised["PanCK"]),
    )
    til = derive_til_score(region)
    phen = derive_immune_phenotype(til, region.CD8_density, region.PanCK_density)
    return TIMEProfile(
        patient_id=submitter_id,
        cohort="tcga_hnsc",
        modality="rna_seq",
        regions={"tumor_core": region, "tumor_margin": None, "adjacent_stroma": None},
        TIL_score=til,
        immune_phenotype=phen,
        provenance=Provenance(
            method=GENOMICS_METHOD,
            version=GENOMICS_VERSION,
            ledger_id=ledger_id,
        ),
    )


def run_genomics_arm(cohort_df: pd.DataFrame, out_dir: Path) -> dict:
    """End-to-end Arm 2: cohort -> TIMEProfile JSONs in out_dir/tcga_hnsc/.

    Returns a summary dict suitable for inclusion in pipeline metrics:
        {n_patients_scored, n_inflamed, n_excluded, n_desert, n_unknown,
         til_mean, til_std}
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_dir = out_dir / "tcga_hnsc"
    cohort_dir.mkdir(parents=True, exist_ok=True)

    scores = cohort_score_matrix(cohort_df)
    if scores.empty:
        return {"n_patients_scored": 0}
    normalised = cohort_zscore_normalize(scores)

    profiles: list[TIMEProfile] = []
    for sid in normalised.index:
        prof = to_time_profile(str(sid), normalised.loc[sid])
        profiles.append(prof)
        (cohort_dir / f"{sid}.json").write_text(
            __import__("json").dumps(prof.to_dict(), indent=2)
        )

    phen_counts = pd.Series([p.immune_phenotype for p in profiles]).value_counts().to_dict()
    return {
        "n_patients_scored": len(profiles),
        "n_inflamed": int(phen_counts.get("inflamed", 0)),
        "n_excluded": int(phen_counts.get("excluded", 0)),
        "n_desert": int(phen_counts.get("desert", 0)),
        "n_unknown": int(phen_counts.get("unknown", 0)),
        "til_mean": float(np.mean([p.TIL_score for p in profiles])),
        "til_std": float(np.std([p.TIL_score for p in profiles])),
    }
