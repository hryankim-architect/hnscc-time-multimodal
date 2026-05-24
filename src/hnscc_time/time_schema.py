"""Common-schema Tumor Immune Microenvironment (TIME) profile.

Both the Genomics arm (Arm 2) and the IHC arm (Arm 1) emit per-patient
TIMEProfile objects in this exact schema, so the Cross-cohort Calibration
arm (Arm 3) can compare them without per-arm adapters.

See `docs/architecture.md` (Common-schema TIME profile section) for the
field-by-field rationale.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, field_validator

# The four cell-type slots tracked across both modalities.
# Names chosen to match the markers in the PMC10571229 multiplex IHC panel
# so a future paired-cohort drop-in does not need a schema migration.
CELL_TYPES: tuple[str, ...] = ("CD3", "CD8", "FoxP3", "PanCK")

# The three regions the IHC arm produces. The Genomics arm has only
# whole-tumor resolution, so it populates `tumor_core` and leaves
# `tumor_margin` + `adjacent_stroma` as None; Arm 3 expands by transferring
# the spatial pattern from the nearest IHC neighbour.
REGION_NAMES: tuple[str, ...] = ("tumor_core", "tumor_margin", "adjacent_stroma")

# Three-class immune-phenotype call used in immunotherapy patient
# selection. Defined per Galon et al. 2016 "The immune contexture" review.
ImmunePhenotype = Literal["inflamed", "excluded", "desert", "unknown"]

# Per-arm modality tag.
Modality = Literal["rna_seq", "mIF", "mIHC", "calibrated_prediction"]


class RegionProfile(BaseModel):
    """Per-cell-type density in one tissue region.

    Units are arbitrary in v0.1 (Arm 2 uses ssGSEA-style enrichment scores;
    Arm 1 uses cells-per-mm^2). Arm 3 calibrates between them.
    """

    CD3_density: float = Field(ge=0.0)
    CD8_density: float = Field(ge=0.0)
    FoxP3_density: float = Field(ge=0.0)
    PanCK_density: float = Field(ge=0.0)


class Provenance(BaseModel):
    """How a TIMEProfile was produced — for the audit trail."""

    method: str  # e.g. "ssGSEA_immune_signatures_v0.1", "cellpose_nuclei_v3"
    version: str  # semver of the producing module
    ledger_id: str | None = None  # sha256 of the audit entry that produced this


class TIMEProfile(BaseModel):
    """Tumor Immune Microenvironment profile in the cross-arm common schema."""

    patient_id: str
    cohort: Literal["tcga_hnsc", "pmc10571229", "deepliif_sample"]
    modality: Modality

    # Region-resolved profile. Genomics arm only fills tumor_core; IHC arm
    # fills all three; calibrated predictions fill all three (transferred).
    regions: dict[str, RegionProfile | None]

    # Tumor-immune lymphocyte score (0-1), derived from cell-type densities.
    TIL_score: float = Field(ge=0.0, le=1.0)
    immune_phenotype: ImmunePhenotype

    # Provenance
    provenance: Provenance

    @field_validator("regions")
    @classmethod
    def _regions_schema(cls, v: dict[str, RegionProfile | None]) -> dict[str, RegionProfile | None]:
        unknown = set(v.keys()) - set(REGION_NAMES)
        if unknown:
            raise ValueError(f"unexpected region keys: {sorted(unknown)}")
        # Ensure all expected keys are present (None is acceptable for
        # genomics-only profiles that lack spatial resolution).
        out = {r: v.get(r) for r in REGION_NAMES}
        return out

    def to_dict(self) -> dict:
        """Stable JSON-serialisable dict matching the architecture spec."""
        return self.model_dump(mode="json")


def derive_til_score(region: RegionProfile) -> float:
    """Heuristic TIL score in [0, 1] from a RegionProfile.

    Computed as the fraction of CD3+ + CD8+ density relative to
    (CD3+ + CD8+ + PanCK+) — i.e. immune density vs tumour-cell density.
    FoxP3 (regulatory T-cells) is intentionally excluded from the
    numerator since high Tregs suppress the inflamed phenotype.

    Edge cases:
        - if total is 0 -> 0.0 (no signal)
        - clipped to [0, 1]
    """
    total = region.CD3_density + region.CD8_density + region.PanCK_density
    if total <= 0:
        return 0.0
    raw = (region.CD3_density + region.CD8_density) / total
    return max(0.0, min(1.0, raw))


def derive_immune_phenotype(
    til_score: float,
    cd8_density: float,
    panck_density: float,
) -> ImmunePhenotype:
    """3-class call from numeric profile, following Galon et al. 2016.

    Decision rule (v0.1, deliberately simple):
        - TIL >= 0.4 AND CD8 within tumor zone -> "inflamed"
        - TIL >= 0.2 but CD8 mostly excluded     -> "excluded"
        - TIL < 0.2                              -> "desert"
        - else                                   -> "unknown"

    The exact thresholds are placeholders; calibration in Arm 3 may
    re-fit them against the IHC reference set.
    """
    if til_score >= 0.4 and cd8_density > 0:
        return "inflamed"
    if til_score >= 0.2:
        return "excluded"
    if til_score < 0.2 and panck_density > 0:
        return "desert"
    return "unknown"
