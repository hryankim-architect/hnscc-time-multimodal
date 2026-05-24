"""End-to-end smoke tests for the P4 three-arm pipeline.

The pipeline is built to degrade gracefully: each arm is skipped (with an
audit-emit recording the skip) if its input data is missing or its optional
extra is uninstalled. The smoke test exercises the "all arms skipped"
no-data path to keep CI fast and dependency-light, then verifies the audit
chain spans `pipeline_start` -> per-arm skip -> `pipeline_end` and is
tamper-detectable.
"""

from __future__ import annotations

import json
from pathlib import Path

from hnscc_time import audit, pipeline


def test_pipeline_runs_and_produces_artifact(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    monkeypatch.delenv("MLFLOW_TRACKING_URI", raising=False)

    out_dir = tmp_path / "artifacts"
    # No data/ dir set up -> all three arms gracefully skip.
    result = pipeline.run_pipeline("smoke", out_dir, data_dir=tmp_path / "nonexistent")

    assert "job_id" in result
    assert "arms" in result
    # All arms recorded as skipped
    assert "skipped" in result["arms"]["arm2_genomics"]
    assert "skipped" in result["arms"]["arm1_ihc"]
    assert "skipped" in result["arms"]["arm3_calibration"]
    # Pipeline elapsed metric is present
    assert result["metrics"]["pipeline_elapsed_ms"] >= 0.0

    artifact_path = Path(result["artifact_path"])
    assert artifact_path.exists()
    payload = json.loads(artifact_path.read_text())
    assert payload["job_id"] == result["job_id"]


def test_audit_chain_is_valid_after_pipeline(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)

    pipeline.run_pipeline("smoke", tmp_path / "artifacts", data_dir=tmp_path / "nonexistent")

    ok, n_entries, first_bad = audit.verify()
    assert ok, f"audit chain invalid at {first_bad}"
    # pipeline_start + 3 skip emits + pipeline_end = 5
    assert n_entries >= 5


def test_audit_chain_detects_tamper(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)

    pipeline.run_pipeline("smoke", tmp_path / "artifacts", data_dir=tmp_path / "nonexistent")
    ledger = audit.DEFAULT_LEDGER

    lines = ledger.read_text().splitlines()
    assert len(lines) >= 2
    tampered = json.loads(lines[0])
    tampered["fields"]["out_dir"] = "/etc/evil"
    lines[0] = json.dumps(tampered, sort_keys=True, separators=(",", ":"))
    ledger.write_text("\n".join(lines) + "\n")

    ok, _, first_bad = audit.verify()
    assert not ok
    assert first_bad is not None


def test_pipeline_arm2_runs_with_real_data(tmp_path: Path, monkeypatch) -> None:
    """When data/tcga_hnsc/ is present, Arm 2 actually scores."""
    import numpy as np

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)

    data_dir = tmp_path / "data"
    star_dir = data_dir / "tcga_hnsc" / "star_counts"
    star_dir.mkdir(parents=True)
    # Manifest
    (data_dir / "tcga_hnsc" / "_subset_manifest.tsv").write_text(
        "file_id\tcase_submitter_id\tfile_name\n"
        "uuid-a\tTCGA-AA-0001\tA.tsv\n"
        "uuid-b\tTCGA-BB-0002\tB.tsv\n"
        "uuid-c\tTCGA-CC-0003\tC.tsv\n"
    )
    # 3 STAR counts files with the same synthetic content
    rng = np.random.default_rng(42)
    genes = (
        ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMB", "PRF1", "FOXP3",
         "IL2RA", "CTLA4", "KRT5", "KRT6A", "KRT14", "KRT17"]
        + [f"BG{i:03d}" for i in range(30)]
    )
    for pid in ["TCGA-AA-0001", "TCGA-BB-0002", "TCGA-CC-0003"]:
        tpm = rng.gamma(1.5, 5.0, len(genes))
        rows = ["# gene-model: GENCODE v36",
                "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded"]
        rows += [f"N_{x}\t\t\t1\t1\t1\t\t\t" for x in ("unmapped", "multi", "no", "amb")]
        for i, (g, t) in enumerate(zip(genes, tpm, strict=True)):
            rows.append(f"ENSG{i:011d}.1\t{g}\tprotein_coding\t1\t1\t1\t{t:.4f}\t0\t0")
        path = star_dir / f"{pid}__test.tsv"
        path.write_text("\n".join(rows) + "\n")
        del t

    result = pipeline.run_pipeline("smoke", tmp_path / "artifacts", data_dir=data_dir)
    arm2 = result["arms"]["arm2_genomics"]
    assert arm2.get("n_patients_scored") == 3
    # 3 per-patient JSONs should exist
    profiles = list((tmp_path / "artifacts" / "time_profiles" / "tcga_hnsc").glob("*.json"))
    assert len(profiles) == 3
