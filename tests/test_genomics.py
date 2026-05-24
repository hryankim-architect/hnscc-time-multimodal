"""Unit tests for src/hnscc_time/genomics.py."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from hnscc_time import genomics
from hnscc_time.time_schema import CELL_TYPES


@pytest.fixture
def tiny_star_counts(tmp_path):
    """Synthetic STAR-counts TSV with 4 immune-signature genes + 30 background."""
    rng = np.random.default_rng(42)
    genes = (
        # Strong signal genes for each immune type
        ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMB", "PRF1", "FOXP3",
         "IL2RA", "CTLA4", "KRT5", "KRT6A", "KRT14", "KRT17"]
        + [f"BG{i:03d}" for i in range(30)]
    )
    tpm = list(rng.gamma(1.5, 5.0, len(genes)))

    rows = ["# gene-model: GENCODE v36",
            "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded",
            "N_unmapped\t\t\t1000\t1000\t1000\t\t\t",
            "N_multimapping\t\t\t500\t500\t500\t\t\t",
            "N_noFeature\t\t\t100\t100\t100\t\t\t",
            "N_ambiguous\t\t\t50\t50\t50\t\t\t",
            ]
    for i, (g, t) in enumerate(zip(genes, tpm, strict=True)):
        rows.append(f"ENSG{i:011d}.1\t{g}\tprotein_coding\t10\t5\t5\t{t:.4f}\t0\t0")
    path = tmp_path / "tiny.tsv"
    path.write_text("\n".join(rows) + "\n")
    return path


def test_load_tpm_drops_n_rows(tiny_star_counts):
    s = genomics.load_tpm(tiny_star_counts)
    # 14 immune genes + 30 background = 44 unique gene_names
    assert len(s) == 44
    assert "CD8A" in s.index
    # N_* rows must not appear
    assert "N_unmapped" not in s.index
    assert (s >= 0).all()


def test_load_tpm_handles_gzipped(tmp_path, tiny_star_counts):
    # Re-write as .gz and confirm same content
    gz_path = tmp_path / "tiny.tsv.gz"
    text = tiny_star_counts.read_text()
    with gzip.open(gz_path, "wt") as fh:
        fh.write(text)
    s_gz = genomics.load_tpm(gz_path)
    s_plain = genomics.load_tpm(tiny_star_counts)
    pd.testing.assert_series_equal(s_gz.sort_index(), s_plain.sort_index())


def test_rank_transform_returns_ranks_in_range(tiny_star_counts):
    s = genomics.load_tpm(tiny_star_counts)
    r = genomics.rank_transform(s)
    assert r.min() >= 1.0
    assert r.max() <= len(s)


def test_signature_score_ignores_missing_genes(tiny_star_counts):
    s = genomics.load_tpm(tiny_star_counts)
    r = genomics.rank_transform(s)
    # Sig with one real gene + two missing -> uses the real one only
    score = genomics.signature_score(r, ["CD8A", "DOES_NOT_EXIST_1", "DOES_NOT_EXIST_2"])
    assert score == r.loc["CD8A"]


def test_score_sample_returns_all_cell_types(tiny_star_counts):
    s = genomics.load_tpm(tiny_star_counts)
    scores = genomics.score_sample(s)
    assert set(scores.keys()) == set(CELL_TYPES)
    for v in scores.values():
        assert v >= 0


def test_cohort_zscore_normalize_within_bounds():
    raw = pd.DataFrame(
        {"CD3": [10, 20, 30], "CD8": [5, 5, 5], "FoxP3": [1, 2, 3], "PanCK": [100, 200, 300]},
        index=["pA", "pB", "pC"],
    )
    z = genomics.cohort_zscore_normalize(raw)
    assert (z >= 0).all().all()
    assert (z <= 5).all().all()
    # Constant column (CD8 = 5 everywhere) should land at the shift midpoint
    assert (z["CD8"] == 2.5).all()


def test_run_genomics_arm_writes_per_patient_jsons(tiny_star_counts, tmp_path):
    # Build a 3-patient mini cohort all pointing at the same fixture
    cohort_df = pd.DataFrame(
        {
            "star_counts_path": [tiny_star_counts, tiny_star_counts, tiny_star_counts],
            "subsite": ["larynx", "oral_cavity", "oropharynx"],
            "primary_diagnosis": ["SCC"] * 3,
            "gender": ["m"] * 3,
            "vital_status": ["alive"] * 3,
            "tobacco_smoking_status": [None] * 3,
        },
        index=["TCGA-AA-0001", "TCGA-BB-0002", "TCGA-CC-0003"],
    )
    out_dir = tmp_path / "profiles"
    summary = genomics.run_genomics_arm(cohort_df, out_dir)
    assert summary["n_patients_scored"] == 3
    # 3 JSONs written
    written = sorted((out_dir / "tcga_hnsc").glob("*.json"))
    assert len(written) == 3
    # Spot-check one
    import json as _json
    data = _json.loads(written[0].read_text())
    assert data["cohort"] == "tcga_hnsc"
    assert data["modality"] == "rna_seq"
    assert "tumor_core" in data["regions"]
    assert data["regions"]["tumor_margin"] is None
