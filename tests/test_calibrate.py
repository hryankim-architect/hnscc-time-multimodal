"""Unit tests for src/hnscc_time/calibrate.py."""

from __future__ import annotations

import json

import pytest

from hnscc_time import calibrate
from hnscc_time.time_schema import (
    CELL_TYPES,
    Provenance,
    RegionProfile,
    TIMEProfile,
)


def _make_profile(pid, cohort, modality, c3=1.0, c8=1.0, fp=0.5, pck=2.0):
    return TIMEProfile(
        patient_id=pid,
        cohort=cohort,
        modality=modality,
        regions={
            "tumor_core": RegionProfile(
                CD3_density=c3, CD8_density=c8, FoxP3_density=fp, PanCK_density=pck,
            ),
            "tumor_margin": None,
            "adjacent_stroma": None,
        },
        TIL_score=0.5,
        immune_phenotype="excluded",
        provenance=Provenance(method="test", version="v0"),
    )


def test_load_profiles_groups_by_cohort(tmp_path):
    d = tmp_path / "profiles"
    (d / "tcga_hnsc").mkdir(parents=True)
    (d / "deepliif_sample").mkdir(parents=True)
    p_g = _make_profile("TCGA-A-1", "tcga_hnsc", "rna_seq")
    p_i = _make_profile("deepliif_ROI_1", "deepliif_sample", "mIHC")
    (d / "tcga_hnsc" / "TCGA-A-1.json").write_text(json.dumps(p_g.to_dict()))
    (d / "deepliif_sample" / "deepliif_ROI_1.json").write_text(json.dumps(p_i.to_dict()))

    grouped = calibrate.load_profiles(d)
    assert "tcga_hnsc" in grouped
    assert "deepliif_sample" in grouped
    assert len(grouped["tcga_hnsc"]) == 1
    assert len(grouped["deepliif_sample"]) == 1


def test_fit_calibration_requires_both_cohorts():
    p_g = [_make_profile("g1", "tcga_hnsc", "rna_seq")]
    with pytest.raises(ValueError):
        calibrate.fit_calibration([], p_g)
    with pytest.raises(ValueError):
        calibrate.fit_calibration(p_g, [])


def test_fit_calibration_returns_all_cell_types():
    ihc = [
        _make_profile("i1", "deepliif_sample", "mIHC", c3=2.0, c8=1.5, pck=3.0),
        _make_profile("i2", "deepliif_sample", "mIHC", c3=1.0, c8=0.5, pck=4.0),
    ]
    gen = [
        _make_profile(f"g{i}", "tcga_hnsc", "rna_seq", c3=1.0 + i * 0.1, c8=1.0, pck=2.0)
        for i in range(5)
    ]
    coefs = calibrate.fit_calibration(ihc, gen, k=3)
    assert set(coefs.keys()) == set(CELL_TYPES)
    for cell in CELL_TYPES:
        assert "a" in coefs[cell]
        assert "b" in coefs[cell]


def test_loo_validation_falls_back_for_small_n():
    ihc_one = [_make_profile("i1", "deepliif_sample", "mIHC")]
    gen = [_make_profile(f"g{i}", "tcga_hnsc", "rna_seq") for i in range(3)]
    loo = calibrate.loo_validate(ihc_one, gen)
    # With n=1 IHC we cannot LOO; should return NaN per cell
    for cell in CELL_TYPES:
        assert loo[cell]["loo_mae"] != loo[cell]["loo_mae"]  # NaN != NaN


def test_predict_time_from_genomics_yields_calibrated_profile():
    coefs = {c: {"a": 1.0, "b": 0.0, "n_pairs": 5} for c in CELL_TYPES}
    g = _make_profile("g1", "tcga_hnsc", "rna_seq", c3=2.0, c8=1.5)
    pred = calibrate.predict_time_from_genomics(g, coefs)
    assert pred.modality == "calibrated_prediction"
    assert pred.regions["tumor_core"].CD3_density == 2.0
    assert pred.provenance.method == calibrate.CALIBRATION_METHOD


def test_run_calibration_arm_end_to_end(tmp_path):
    d = tmp_path / "profiles"
    (d / "tcga_hnsc").mkdir(parents=True)
    (d / "deepliif_sample").mkdir(parents=True)
    # 5 genomics + 3 IHC for a real LOO
    for i in range(5):
        p = _make_profile(f"g{i}", "tcga_hnsc", "rna_seq", c3=1.0 + i * 0.2, c8=1.0)
        (d / "tcga_hnsc" / f"g{i}.json").write_text(json.dumps(p.to_dict()))
    for i in range(3):
        p = _make_profile(f"i{i}", "deepliif_sample", "mIHC", c3=2.0 - i * 0.3, c8=1.5)
        (d / "deepliif_sample" / f"i{i}.json").write_text(json.dumps(p.to_dict()))

    out_dir = tmp_path / "calibration"
    summary = calibrate.run_calibration_arm(
        profiles_dir=d, cohort_df_or_dir=tmp_path, out_dir=out_dir,
    )
    assert summary["n_ihc_reference"] == 3
    assert summary["n_genomics_input"] == 5
    assert summary["n_calibrated_predictions"] == 5
    # Calibrated predictions written
    assert len(list((out_dir / "calibrated_predictions").glob("*.json"))) == 5
    # Calibration report exists
    assert (out_dir / "calibration_report.json").exists()
