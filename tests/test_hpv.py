"""Tests for Arm 4 — HPV± overall-survival stratification.

Survival stats run on small SYNTHETIC cohorts with a planted signal, so the test
is deterministic and offline. The real GDC fetch is exercised by ``make data``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("lifelines")

from hnscc_time import hpv  # noqa: E402


def _planted_cohort() -> pd.DataFrame:
    """HPV+ clearly protective, but NOT completely separated.

    Both groups have events across overlapping time ranges (HPV+ fewer, HPV-
    more), so the Cox model converges normally and the CI is finite — a strong
    but well-posed planted signal.
    """
    rng = np.random.default_rng(42)
    n = 40
    pos = pd.DataFrame({
        "submitter_id": [f"P{i:02d}" for i in range(n)],
        "hpv_status": "positive",
        "os_days": rng.uniform(300, 2200, n).round(1),
        "os_event": rng.binomial(1, 0.25, n),  # ~few events
    })
    neg = pd.DataFrame({
        "submitter_id": [f"N{i:02d}" for i in range(n)],
        "hpv_status": "negative",
        "os_days": rng.uniform(150, 1600, n).round(1),
        "os_event": rng.binomial(1, 0.75, n),  # ~many events
    })
    return pd.concat([pos, neg], ignore_index=True)


class TestHpvSurvivalSummary:
    def test_planted_signal_is_protective(self):
        s = hpv.hpv_survival_summary(_planted_cohort())
        assert s["n_hpv_positive"] == 40
        assert s["n_hpv_negative"] == 40
        assert s["cox_hr_hpv_pos_vs_neg"] < 1.0  # HPV+ protective
        assert s["direction_protective"] is True
        assert s["significant_p05"] is True  # planted signal is strong
        lo, hi = s["cox_hr_95ci"]
        # Well-posed cohort -> finite CI that brackets the point estimate.
        assert isinstance(lo, float) and isinstance(hi, float)
        assert lo <= s["cox_hr_hpv_pos_vs_neg"] <= hi

    def test_summary_is_json_serializable(self):
        import json

        json.dumps(hpv.hpv_survival_summary(_planted_cohort()))  # must not raise


class TestLoadHpvSurvival:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            hpv.load_hpv_survival(tmp_path / "nope.tsv")

    def test_filters_to_known_status(self, tmp_path):
        p = tmp_path / "hpv_status.tsv"
        p.write_text(
            "submitter_id\thpv_status\tos_days\tos_event\n"
            "A\tpositive\t100.0\t0\n"
            "B\tnegative\t200.0\t1\n"
            "C\tunknown\t300.0\t1\n",
            encoding="utf-8",
        )
        df = hpv.load_hpv_survival(p)
        assert set(df["hpv_status"]) == {"positive", "negative"}  # 'unknown' dropped
        assert len(df) == 2
