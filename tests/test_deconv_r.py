"""Tests for Arm 5 — Python<->R bridge deconvolution cross-check.

The pure scorer (`cross_method_agreement`) is tested offline. The R bridge itself
runs only where `Rscript` is installed; those tests skip otherwise (the arm is
designed to skip-if-no-R, like the IHC arm skips without the cellpose extra).
"""

from __future__ import annotations

import shutil

import pandas as pd
import pytest

from hnscc_time import deconv_r
from hnscc_time.time_schema import CELL_TYPES

HAS_R = shutil.which("Rscript") is not None


class TestCrossMethodAgreement:
    def _frame(self, scale: float) -> pd.DataFrame:
        # 5 samples, monotonically increasing scores per cell type.
        return pd.DataFrame(
            {c: [scale * v for v in (1, 2, 3, 4, 5)] for c in CELL_TYPES},
            index=[f"S{i}" for i in range(5)],
        )

    def test_monotone_methods_agree(self):
        py = self._frame(1.0)
        r = self._frame(10.0)  # same ordering, different magnitude
        a = deconv_r.cross_method_agreement(py, r)
        assert a["n_samples"] == 5
        assert a["mean_spearman"] == pytest.approx(1.0)
        assert set(a["cell_types"]) <= set(CELL_TYPES)

    def test_too_few_samples_is_none(self):
        idx = ["S0", "S1"]
        py = pd.DataFrame({c: [1.0, 2.0] for c in CELL_TYPES}, index=idx)
        a = deconv_r.cross_method_agreement(py, py)
        # <3 shared samples -> per-cell rho is None, mean is None.
        assert a["mean_spearman"] is None

    def test_zero_variance_is_none(self):
        idx = [f"S{i}" for i in range(5)]
        flat = pd.DataFrame({c: [1.0] * 5 for c in CELL_TYPES}, index=idx)
        a = deconv_r.cross_method_agreement(flat, flat)
        assert all(v is None for v in a["spearman_per_cell"].values())


class TestSkipWithoutR:
    @pytest.mark.skipif(HAS_R, reason="R present; skip path not exercised")
    def test_run_deconv_arm_raises_without_r(self, tmp_path):
        with pytest.raises(deconv_r.RDeconvUnavailable):
            deconv_r.run_deconv_arm(pd.DataFrame(), tmp_path)

    def test_r_available_is_bool(self):
        assert isinstance(deconv_r.r_available(), bool)


@pytest.mark.skipif(not HAS_R, reason="Rscript not installed")
class TestRBridgeIntegration:
    def test_r_script_runs_on_tiny_matrix(self):
        # genes x samples TPM; include some marker genes so scores are non-trivial.
        genes = ["CD3D", "CD3E", "CD8A", "FOXP3", "KRT5", "ACTB", "GAPDH"]
        tpm = pd.DataFrame(
            {f"S{j}": [float((i + 1) * (j + 1)) for i in range(len(genes))] for j in range(4)},
            index=genes,
        )
        scores = deconv_r.run_r_deconvolution(tpm)
        assert scores.shape[0] == 4  # 4 samples
        assert set(CELL_TYPES) <= set(scores.columns)
