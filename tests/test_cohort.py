"""Unit tests for src/hnscc_time/cohort.py."""

from __future__ import annotations

import pytest

from hnscc_time import cohort


def test_normalise_subsite_canonical():
    assert cohort._normalise_subsite("Larynx, NOS") == "larynx"
    assert cohort._normalise_subsite("Mouth, NOS") == "oral_cavity"
    assert cohort._normalise_subsite("Tongue, NOS") == "oral_cavity"
    assert cohort._normalise_subsite("Base of tongue, NOS") == "oropharynx"
    assert cohort._normalise_subsite("Tonsil, NOS") == "oropharynx"
    assert cohort._normalise_subsite("Hypopharynx, NOS") == "hypopharynx"
    assert cohort._normalise_subsite("Nasopharynx, NOS") == "nasopharynx"


def test_normalise_subsite_unknown_returns_other():
    assert cohort._normalise_subsite("Lymph nodes of pelvis") == "other"


def test_normalise_subsite_none_passthrough():
    assert cohort._normalise_subsite(None) is None
    assert cohort._normalise_subsite("") == "other"


def test_load_cohort_missing_manifest_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        cohort.load_cohort(tmp_path)


def test_load_cohort_minimal_no_clinical(tmp_path):
    data = tmp_path / "tcga_hnsc"
    (data / "star_counts").mkdir(parents=True)
    # Manifest with 2 patients, both have a STAR file
    (data / "_subset_manifest.tsv").write_text(
        "file_id\tcase_submitter_id\tfile_name\n"
        "uuid-a\tTCGA-AA-0001\tA.rna_seq.augmented_star_gene_counts.tsv\n"
        "uuid-b\tTCGA-BB-0002\tB.rna_seq.augmented_star_gene_counts.tsv\n"
    )
    (data / "star_counts" / "TCGA-AA-0001__A.rna_seq.augmented_star_gene_counts.tsv").write_text("")
    (data / "star_counts" / "TCGA-BB-0002__B.rna_seq.augmented_star_gene_counts.tsv").write_text("")

    df = cohort.load_cohort(tmp_path)
    assert len(df) == 2
    assert "TCGA-AA-0001" in df.index
    assert df.loc["TCGA-AA-0001", "star_counts_path"].name.endswith(".tsv")
