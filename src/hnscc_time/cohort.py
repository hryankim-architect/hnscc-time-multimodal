"""TCGA-HNSC cohort loader.

Reads:
    - `data/tcga_hnsc/_subset_manifest.tsv` — file_id / case_submitter_id / file_name
      produced by `scripts/download_tcga_hnsc.sh` (v0.0 download).
    - `data/tcga_hnsc/clinical.tsv` — GDC clinical metadata for the
      selected cases.

Joins the two into a single tidy DataFrame indexed by submitter_id with
columns the genomics arm needs for stratification + downstream analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class CohortRow:
    """Per-patient metadata + STAR-counts file location."""

    submitter_id: str
    star_counts_path: Path
    primary_diagnosis: str | None
    subsite: str | None
    gender: str | None
    vital_status: str | None
    tobacco_smoking_status: str | None


def _normalise_subsite(raw: str | None) -> str | None:
    """Collapse the long TCGA subsite labels into 5 broad HNSCC sites.

    The TCGA strings are messy: 'Larynx, NOS', 'Mouth, NOS',
    'Overlapping lesion of lip, oral cavity and pharynx', etc. The 5
    canonical HNSCC subsites are oral cavity, oropharynx, larynx,
    hypopharynx, nasopharynx — we map best-effort.
    """
    if raw is None or not isinstance(raw, str):
        return None
    r = raw.lower()
    if "larynx" in r:
        return "larynx"
    if "hypopharynx" in r:
        return "hypopharynx"
    if "nasopharynx" in r:
        return "nasopharynx"
    if "oropharynx" in r or "tonsil" in r or "base of tongue" in r:
        return "oropharynx"
    if "mouth" in r or "tongue" in r or "lip" in r or "oral cavity" in r or "palate" in r or "gum" in r:
        return "oral_cavity"
    return "other"


def load_cohort(data_dir: Path) -> pd.DataFrame:
    """Load the TCGA-HNSC subset cohort.

    Returns a DataFrame indexed by `submitter_id` with columns:
        - star_counts_path (Path)
        - primary_diagnosis (str)
        - subsite (5-class normalised string)
        - gender, vital_status, tobacco_smoking_status (raw)
    """
    manifest_path = data_dir / "tcga_hnsc" / "_subset_manifest.tsv"
    clinical_path = data_dir / "tcga_hnsc" / "clinical.tsv"
    star_dir = data_dir / "tcga_hnsc" / "star_counts"

    if not manifest_path.exists():
        raise FileNotFoundError(
            f"{manifest_path} missing. Run scripts/download_tcga_hnsc.sh first."
        )

    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest = manifest.rename(columns={"case_submitter_id": "submitter_id"})

    # Resolve per-patient STAR counts file path. The download script writes
    # `${CASE_ID}__${ORIGINAL_FILE_NAME}` so we glob.
    paths: list[Path | None] = []
    for sid in manifest["submitter_id"]:
        candidates = sorted(star_dir.glob(f"{sid}__*.tsv"))
        paths.append(candidates[0] if candidates else None)
    manifest["star_counts_path"] = paths

    # Load and align clinical
    if clinical_path.exists():
        clinical = pd.read_csv(clinical_path, sep="\t")
        # GDC TSV uses dotted-field names; we keep only what the genomics
        # arm + Arm 3 actually need.
        rename = {
            "demographic.gender": "gender",
            "demographic.vital_status": "vital_status",
            "diagnoses.0.primary_diagnosis": "primary_diagnosis",
            "diagnoses.0.tissue_or_organ_of_origin": "subsite_raw",
            "exposures.0.tobacco_smoking_status": "tobacco_smoking_status",
            "submitter_id": "submitter_id",
        }
        present = {k: v for k, v in rename.items() if k in clinical.columns}
        clinical = clinical[list(present.keys())].rename(columns=present)
        clinical["subsite"] = clinical["subsite_raw"].map(_normalise_subsite) if "subsite_raw" in clinical else None
        if "subsite_raw" in clinical.columns:
            clinical = clinical.drop(columns=["subsite_raw"])
    else:
        clinical = pd.DataFrame(columns=["submitter_id", "primary_diagnosis", "subsite", "gender", "vital_status", "tobacco_smoking_status"])

    df = manifest.merge(clinical, on="submitter_id", how="left")
    return df.set_index("submitter_id")


def cohort_rows(df: pd.DataFrame) -> list[CohortRow]:
    """Iterate the cohort as typed rows."""
    out: list[CohortRow] = []
    for sid, row in df.iterrows():
        out.append(
            CohortRow(
                submitter_id=str(sid),
                star_counts_path=Path(row["star_counts_path"]) if pd.notna(row["star_counts_path"]) else Path(""),
                primary_diagnosis=row.get("primary_diagnosis") if pd.notna(row.get("primary_diagnosis")) else None,
                subsite=row.get("subsite") if pd.notna(row.get("subsite")) else None,
                gender=row.get("gender") if pd.notna(row.get("gender")) else None,
                vital_status=row.get("vital_status") if pd.notna(row.get("vital_status")) else None,
                tobacco_smoking_status=row.get("tobacco_smoking_status") if pd.notna(row.get("tobacco_smoking_status")) else None,
            )
        )
    return out
