"""Unit tests for the common-schema TIMEProfile."""

from __future__ import annotations

import pytest

from hnscc_time.time_schema import (
    CELL_TYPES,
    REGION_NAMES,
    Provenance,
    RegionProfile,
    TIMEProfile,
    derive_immune_phenotype,
    derive_til_score,
)


def test_cell_types_match_ihc_panel():
    assert CELL_TYPES == ("CD3", "CD8", "FoxP3", "PanCK")


def test_region_names_canonical_three():
    assert REGION_NAMES == ("tumor_core", "tumor_margin", "adjacent_stroma")


def test_region_profile_rejects_negative():
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        RegionProfile(CD3_density=-1.0, CD8_density=0.0, FoxP3_density=0.0, PanCK_density=0.0)


def test_til_score_zero_when_no_signal():
    r = RegionProfile(CD3_density=0.0, CD8_density=0.0, FoxP3_density=0.0, PanCK_density=0.0)
    assert derive_til_score(r) == 0.0


def test_til_score_one_when_only_immune():
    r = RegionProfile(CD3_density=1.0, CD8_density=1.0, FoxP3_density=0.0, PanCK_density=0.0)
    assert derive_til_score(r) == 1.0


def test_til_score_excludes_foxp3_from_numerator():
    # FoxP3 should not push TIL up — only CD3 + CD8 are in the numerator.
    no_foxp3 = RegionProfile(CD3_density=1.0, CD8_density=1.0, FoxP3_density=0.0, PanCK_density=1.0)
    with_foxp3 = RegionProfile(CD3_density=1.0, CD8_density=1.0, FoxP3_density=10.0, PanCK_density=1.0)
    assert derive_til_score(no_foxp3) == derive_til_score(with_foxp3)


def test_immune_phenotype_inflamed():
    assert derive_immune_phenotype(til_score=0.7, cd8_density=2.0, panck_density=1.0) == "inflamed"


def test_immune_phenotype_excluded():
    assert derive_immune_phenotype(til_score=0.25, cd8_density=0.0, panck_density=3.0) == "excluded"


def test_immune_phenotype_desert():
    assert derive_immune_phenotype(til_score=0.05, cd8_density=0.0, panck_density=2.0) == "desert"


def test_time_profile_roundtrip():
    region = RegionProfile(CD3_density=1.0, CD8_density=0.5, FoxP3_density=0.1, PanCK_density=2.0)
    prof = TIMEProfile(
        patient_id="TEST-001",
        cohort="tcga_hnsc",
        modality="rna_seq",
        regions={"tumor_core": region, "tumor_margin": None, "adjacent_stroma": None},
        TIL_score=0.43,
        immune_phenotype="excluded",
        provenance=Provenance(method="test_method", version="v0.0.0"),
    )
    d = prof.to_dict()
    assert d["patient_id"] == "TEST-001"
    assert d["regions"]["tumor_core"]["CD3_density"] == 1.0
    # Round-trip back
    prof2 = TIMEProfile.model_validate(d)
    assert prof2.patient_id == prof.patient_id
    assert prof2.regions["tumor_core"].CD3_density == 1.0


def test_time_profile_rejects_unknown_region():
    from pydantic import ValidationError

    region = RegionProfile(CD3_density=0.0, CD8_density=0.0, FoxP3_density=0.0, PanCK_density=0.0)
    with pytest.raises(ValidationError):
        TIMEProfile(
            patient_id="X",
            cohort="tcga_hnsc",
            modality="rna_seq",
            regions={"some_unknown_region": region},
            TIL_score=0.0,
            immune_phenotype="unknown",
            provenance=Provenance(method="m", version="v"),
        )
