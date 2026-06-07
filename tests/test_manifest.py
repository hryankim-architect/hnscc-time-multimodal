"""Static guards on data/manifest.yaml (the committed checksum ledger).

Runs offline. Asserts the genomics arm is real (no v0.0 ``TBD`` placeholders,
no leftover wrong-Zenodo guess), the n=50 STAR inputs are GDC-open and
sha256-pinned, the clinical block is pinned, and the IHC arm cites its resolved
TCIA source. The network fetch path is exercised by ``make data`` /
``scripts/download_*.sh``, not unit-tested here.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "data" / "manifest.yaml"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
TCIA_DOI = "10.7937/TCIA.2020.T90F-WB82"


def _text() -> str:
    return MANIFEST.read_text(encoding="utf-8")


def _load() -> dict:
    return yaml.safe_load(_text())


class TestNoPlaceholders:
    def test_no_tbd_or_wrong_source(self):
        t = _text().lower()
        assert "tbd" not in t, "manifest still has v0.0 TBD placeholders"
        # The wrong Zenodo *guess* must not be used as a source URL/record.
        # (The word may still appear in a comment explaining the correction.)
        assert "zenodo.org" not in t, "manifest still points at a Zenodo URL"
        assert "8367318" not in t, "manifest still references the wrong Zenodo record id"

    def test_top_level_shape(self):
        m = _load()
        assert set(m) >= {"clinical", "ihc", "inputs"}


class TestClinicalBlock:
    def test_pinned(self):
        c = _load()["clinical"]
        assert c["endpoint"] == "https://api.gdc.cancer.gov/cases"
        assert c["path"] == "tcga_hnsc/clinical.tsv"
        assert HEX64.match(c["sha256"]), "clinical sha256 must be a real 64-hex digest"


class TestIhcBlock:
    def test_resolved_to_tcia(self):
        ihc = _load()["ihc"]
        assert ihc["doi"] == TCIA_DOI
        assert "cancerimagingarchive.net" in ihc["url"]
        assert str(ihc["license"]).upper().startswith("CC-BY")
        assert ihc["patients"] == 8

    def test_not_pending_and_pinned(self):
        import hashlib

        ihc = _load()["ihc"]
        # IHC checksums are resolved, not "pending".
        assert "pending" not in str(ihc).lower()
        assert ihc["rois"] == 3212
        assert HEX64.match(ihc["checksums_file_sha256"])
        # The committed ROI ledger exists and its sha256 matches the pinned value.
        ledger = REPO / "data" / ihc["checksums_file"]
        assert ledger.exists(), f"missing committed IHC ledger {ledger}"
        got = hashlib.sha256(ledger.read_bytes()).hexdigest()
        assert got == ihc["checksums_file_sha256"], "IHC ledger sha256 != pinned value"

    def test_ledger_rows_are_real(self):
        ihc = _load()["ihc"]
        ledger = REPO / "data" / ihc["checksums_file"]
        rows = [
            ln for ln in ledger.read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#")
        ]
        assert len(rows) == 3212, f"expected 3212 ROI rows, got {len(rows)}"
        seen: set[str] = set()
        for ln in rows:
            rel, sha, url = ln.split("\t")
            assert HEX64.match(sha), f"bad sha for {rel}"
            assert "cancerimagingarchive.net" in url
            assert rel not in seen, f"duplicate rel_path {rel}"
            seen.add(rel)


class TestGenomicsInputs:
    def test_fifty_real_open_tier_entries(self):
        inputs = _load()["inputs"]
        assert len(inputs) == 50, f"expected n=50 STAR inputs, got {len(inputs)}"
        seen: set[str] = set()
        for e in inputs:
            assert e["url"].startswith("https://api.gdc.cancer.gov/data/"), e["url"]
            assert e["path"].startswith("tcga_hnsc/star_counts/")
            assert "__" in e["path"], "path must use the <case>__<file> layout"
            assert e["path"].endswith(".tsv"), "STAR counts are plain TSV (not .gz)"
            assert HEX64.match(e["sha256"]), f"{e['path']} has a non-real sha256"
            assert e["path"] not in seen, f"duplicate path {e['path']}"
            seen.add(e["path"])

    def test_inputs_sha256_are_distinct(self):
        shas = [e["sha256"] for e in _load()["inputs"]]
        assert len(set(shas)) == len(shas), "duplicate sha256 across distinct STAR files"

    def test_file_id_in_provenance(self):
        # url ends with the GDC file_id; that id should be recorded in source too.
        for e in _load()["inputs"]:
            file_id = e["url"].rsplit("/", 1)[1]
            assert file_id in e["source"], f"{file_id} missing from its source line"


class TestFetchWiring:
    """make data (pipeline fetch) reproduces the layout make run reads — offline."""

    def test_star_meta_parses_path_and_url(self):
        from hnscc_time import pipeline

        url = "https://api.gdc.cancer.gov/data/abc-123"
        path = "tcga_hnsc/star_counts/TCGA-XX-1234__file.rna_seq.augmented_star_gene_counts.tsv"
        fid, case, fn = pipeline._star_meta(path, url)
        assert fid == "abc-123"
        assert case == "TCGA-XX-1234"
        assert fn == "file.rna_seq.augmented_star_gene_counts.tsv"

    def test_subset_manifest_matches_cohort_loader_columns(self, tmp_path):
        from hnscc_time import pipeline

        rows = [("fid2", "TCGA-B", "b.tsv"), ("fid1", "TCGA-A", "a.tsv")]
        sm = pipeline.write_subset_manifest(rows, tmp_path)
        lines = sm.read_text().splitlines()
        # Header is exactly what cohort.load_cohort expects (it renames
        # case_submitter_id -> submitter_id).
        assert lines[0] == "file_id\tcase_submitter_id\tfile_name"
        assert len(lines) - 1 == 2
        # Rows are sorted (deterministic output).
        assert lines[1].startswith("fid1\t") and lines[2].startswith("fid2\t")

    def test_subset_manifest_built_from_real_inputs(self, tmp_path):
        from hnscc_time import pipeline

        inputs = _load()["inputs"]
        rows = [pipeline._star_meta(e["path"], e["url"]) for e in inputs]
        sm = pipeline.write_subset_manifest(rows, tmp_path)
        data = sm.read_text().splitlines()[1:]
        assert len(data) == 50
        assert len({ln.split("\t")[1] for ln in data}) == 50  # 50 distinct cases
